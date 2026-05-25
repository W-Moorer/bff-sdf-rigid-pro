#include "atlas/BFFAdapter.h"
#include "atlas/ContactAtlas.h"
#include "atlas/PatchBuilder.h"
#include "contact/ChartProjector.h"
#include "contact/ContactDetector.h"
#include "core/MeshRepair.h"
#include "experiments/ProceduralMeshes.h"
#include "io/DebugExport.h"
#include "io/MeshIO.h"
#include "sdf/AnalyticSDF.h"
#include "sdf/GridSDF.h"
#include "sdf/MeshSDF.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace bff_sdf;

namespace {

enum class ChartMode {
    OfficialBff,
    PlanarFallback
};

struct Row {
    std::string scene;
    std::string sceneClass;
    std::string method;
    int meshFaces = 0;
    int sdfResolution = 0;
    double contactArea = 0.0;
    double referenceArea = 0.0;
    double areaError = 0.0;
    double maxPenetration = 0.0;
    double referencePenetration = 0.0;
    double penetrationError = 0.0;
    double normalErrorDeg = 0.0;
    double contourHausdorff = 0.0;
    double runtimeMs = 0.0;
    double integrationMs = 0.0;
    double projectionMs = 0.0;
    double otherMs = 0.0;
    int sdfQueries = 0;
    int projectionCount = 0;
    int projectionFailures = 0;
    double memoryMb = 0.0;
    std::string referenceType;
    std::string chartMethod;
};

struct AssetStat {
    std::string scene;
    int fullFaces = 0;
    int patchFaces = 0;
    int patchVertices = 0;
    int boundaryLoops = 0;
    int eulerCharacteristic = 0;
    bool diskLike = false;
    std::string chartMethod;
    bool remeshed = false;
    int repairRemovedDegenerate = 0;
    int repairRemovedDuplicate = 0;
    int repairFlippedFaces = 0;
    int projectionFailures = 0;
    std::string note;
};

struct ChartQualityRow {
    std::string scene;
    std::string chartMethod;
    int patchFaces = 0;
    int uvFlippedTriangles = 0;
    double meanAngleDistortionDeg = 0.0;
    double maxAngleDistortionDeg = 0.0;
    double areaLogStd = 0.0;
    double maxJacobianCondition = 0.0;
    bool diskLike = false;
    bool usedFallback = false;
    bool remeshed = false;
};

struct ReferenceResult {
    double area = 0.0;
    double maxPenetration = 0.0;
    double meanPenetration = 0.0;
};

struct PolyVertex {
    Eigen::Vector3d x = Eigen::Vector3d::Zero();
    double g = 0.0;
};

std::string boolText(bool v)
{
    return v ? "true" : "false";
}

std::string chartSuffix(const atlas::ContactChart& chart)
{
    return chart.parameterizationMethod == "official-bff" ? "bff" : "planar";
}

double triangleArea3D(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
{
    return 0.5 * (b - a).cross(c - a).norm();
}

double triangleArea2D(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c)
{
    return 0.5 * std::abs((b - a).x() * (c - a).y() - (b - a).y() * (c - a).x());
}

double angleBetween(const Eigen::Vector3d& a, const Eigen::Vector3d& b)
{
    const double denom = a.norm() * b.norm();
    if (denom <= 1e-15) return 0.0;
    const double c = std::max(-1.0, std::min(1.0, a.dot(b) / denom));
    return std::acos(c);
}

double angleBetween2D(const Eigen::Vector2d& a, const Eigen::Vector2d& b)
{
    const double denom = a.norm() * b.norm();
    if (denom <= 1e-15) return 0.0;
    const double c = std::max(-1.0, std::min(1.0, a.dot(b) / denom));
    return std::acos(c);
}

atlas::ContactChart buildChartFromPatch(const atlas::Patch& patch, ChartMode mode)
{
    atlas::BFFAdapter adapter;
    atlas::BFFAdapterOptions options;
    options.preferOfficialBff = mode == ChartMode::OfficialBff;
    const auto result = adapter.parameterize(patch, options);
    if (result.success) {
        return atlas::makeContactChart(0, patch, result.UV, result.Fuv,
                                       result.usedOfficialBff ? "official-bff" : "fallback-planar",
                                       result.usedFallback);
    }

    return atlas::makeContactChart(0, patch, experiments::radialUV(patch.localMesh), patch.localMesh.F,
                                   "fallback-planar", true);
}

atlas::ContactChart chartFromMesh(const core::Mesh& mesh, ChartMode mode)
{
    auto patch = atlas::PatchBuilder::buildWholeMeshPatch(mesh);
    return buildChartFromPatch(patch, patch.diskLike ? mode : ChartMode::PlanarFallback);
}

core::Mesh normalizeMesh(core::Mesh mesh, double targetExtent)
{
    const auto box = mesh.bounds();
    const Eigen::Vector3d center = box.center();
    const double scale = targetExtent / std::max(1e-12, box.extent().maxCoeff());
    for (int i = 0; i < mesh.V.rows(); ++i) {
        Eigen::Vector3d p = mesh.V.row(i);
        p = (p - center) * scale;
        mesh.V.row(i) = p.transpose();
    }
    mesh.computeNormals();
    return mesh;
}

core::Mesh makeSphericalCapPatch(double radius, double capAngleDeg, int rings, int segments)
{
    core::Mesh m;
    m.V.resize(1 + rings * segments, 3);
    m.V.row(0) = Eigen::Vector3d(0, 0, -radius);
    const double capAngle = capAngleDeg * experiments::pi / 180.0;
    int row = 1;
    for (int r = 1; r <= rings; ++r) {
        const double theta = experiments::pi - capAngle * static_cast<double>(r) / static_cast<double>(rings);
        const double z = radius * std::cos(theta);
        const double xy = radius * std::sin(theta);
        for (int s = 0; s < segments; ++s) {
            const double phi = 2.0 * experiments::pi * static_cast<double>(s) / static_cast<double>(segments);
            m.V.row(row++) = Eigen::Vector3d(xy * std::cos(phi), xy * std::sin(phi), z);
        }
    }
    m.F.resize(segments + (rings - 1) * 2 * segments, 3);
    row = 0;
    auto vid = [segments](int ringOneBased, int s) {
        return 1 + (ringOneBased - 1) * segments + ((s + segments) % segments);
    };
    for (int s = 0; s < segments; ++s) m.F.row(row++) = Eigen::Vector3i(0, vid(1, s + 1), vid(1, s));
    for (int r = 1; r < rings; ++r) {
        for (int s = 0; s < segments; ++s) {
            const int a = vid(r, s), b = vid(r, s + 1), c = vid(r + 1, s), d = vid(r + 1, s + 1);
            m.F.row(row++) = Eigen::Vector3i(a, b, c);
            m.F.row(row++) = Eigen::Vector3i(b, d, c);
        }
    }
    m.computeNormals();
    return m;
}

core::Mesh makeEllipsoidCapPatch(const Eigen::Vector3d& radii, double capAngleDeg, int rings, int segments)
{
    auto m = makeSphericalCapPatch(1.0, capAngleDeg, rings, segments);
    for (int i = 0; i < m.V.rows(); ++i) {
        m.V(i, 0) *= radii.x();
        m.V(i, 1) *= radii.y();
        m.V(i, 2) *= radii.z();
    }
    m.computeNormals();
    return m;
}

std::vector<std::vector<int>> selectedFaceComponents(const core::Mesh& mesh, const std::vector<int>& selectedFaces)
{
    const auto adjacency = mesh.faceAdjacency();
    std::set<int> selected(selectedFaces.begin(), selectedFaces.end());
    std::set<int> visited;
    std::vector<std::vector<int>> components;

    for (int seed : selectedFaces) {
        if (visited.count(seed)) continue;
        std::vector<int> component;
        std::queue<int> q;
        q.push(seed);
        visited.insert(seed);
        while (!q.empty()) {
            const int f = q.front();
            q.pop();
            component.push_back(f);
            for (int e = 0; e < 3; ++e) {
                const int nb = adjacency[f][e];
                if (selected.count(nb) && !visited.count(nb)) {
                    visited.insert(nb);
                    q.push(nb);
                }
            }
        }
        components.push_back(component);
    }

    std::sort(components.begin(), components.end(), [](const auto& a, const auto& b) {
        return a.size() > b.size();
    });
    return components;
}

Eigen::Vector3d faceCentroid(const core::Mesh& mesh, int faceId)
{
    const Eigen::Vector3i f = mesh.F.row(faceId);
    return (mesh.V.row(f[0]) + mesh.V.row(f[1]) + mesh.V.row(f[2])) / 3.0;
}

atlas::Patch remeshPatchBoundaryFan(const atlas::Patch& patch)
{
    if (patch.boundaryLoop.size() < 3) return patch;

    core::Mesh remeshed;
    remeshed.V.resize(static_cast<int>(patch.boundaryLoop.size()) + 1, 3);
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (int vid : patch.boundaryLoop) center += patch.localMesh.V.row(vid).transpose();
    center /= static_cast<double>(patch.boundaryLoop.size());
    remeshed.V.row(0) = center.transpose();
    for (int i = 0; i < static_cast<int>(patch.boundaryLoop.size()); ++i) {
        remeshed.V.row(i + 1) = patch.localMesh.V.row(patch.boundaryLoop[i]);
    }

    remeshed.F.resize(static_cast<int>(patch.boundaryLoop.size()), 3);
    for (int i = 0; i < static_cast<int>(patch.boundaryLoop.size()); ++i) {
        const int j = (i + 1) % static_cast<int>(patch.boundaryLoop.size());
        remeshed.F.row(i) = Eigen::Vector3i(0, i + 1, j + 1);
    }
    remeshed.computeNormals();

    Eigen::Vector3d originalNormal = Eigen::Vector3d::Zero();
    for (int i = 0; i < patch.localMesh.FN.rows(); ++i) originalNormal += patch.localMesh.FN.row(i).transpose();
    Eigen::Vector3d remeshNormal = Eigen::Vector3d::Zero();
    for (int i = 0; i < remeshed.FN.rows(); ++i) remeshNormal += remeshed.FN.row(i).transpose();
    if (originalNormal.norm() > 1e-12 && remeshNormal.norm() > 1e-12 && originalNormal.dot(remeshNormal) < 0.0) {
        for (int i = 0; i < remeshed.F.rows(); ++i) std::swap(remeshed.F(i, 1), remeshed.F(i, 2));
        remeshed.computeNormals();
    }

    return atlas::PatchBuilder::buildWholeMeshPatch(remeshed);
}

bool officialBffPasses(const atlas::Patch& patch)
{
    atlas::BFFAdapter adapter;
    atlas::BFFAdapterOptions options;
    options.preferOfficialBff = true;
    const auto result = adapter.parameterize(patch, options);
    return result.success && result.usedOfficialBff;
}

bool considerPatchCandidate(const atlas::Patch& candidate,
                            ChartMode chartMode,
                            atlas::Patch& bestPatch,
                            bool& bestRemeshed,
                            int& bestScore,
                            bool& bestBffPass)
{
    if (!candidate.diskLike) return false;

    const int score = candidate.localMesh.F.rows();
    if (chartMode != ChartMode::OfficialBff) {
        if (score > bestScore) {
            bestPatch = candidate;
            bestRemeshed = false;
            bestScore = score;
            bestBffPass = false;
        }
        return true;
    }

    const bool bffPass = officialBffPasses(candidate);
    if (bffPass && (!bestBffPass || score > bestScore)) {
        bestPatch = candidate;
        bestRemeshed = false;
        bestScore = score;
        bestBffPass = true;
        return true;
    }

    const auto remeshed = remeshPatchBoundaryFan(candidate);
    const bool remeshBffPass = remeshed.diskLike && officialBffPasses(remeshed);
    if (remeshBffPass && (!bestBffPass || remeshed.localMesh.F.rows() > bestScore)) {
        bestPatch = remeshed;
        bestRemeshed = true;
        bestScore = remeshed.localMesh.F.rows();
        bestBffPass = true;
        return true;
    }

    if (!bestBffPass && score > bestScore) {
        bestPatch = candidate;
        bestRemeshed = false;
        bestScore = score;
    }
    return true;
}

bool findContactDiskPatch(const core::Mesh& mesh,
                          ChartMode chartMode,
                          atlas::Patch& outPatch,
                          bool& usedRemesh,
                          std::string& note)
{
    const auto box = mesh.bounds();
    const Eigen::Vector3d contactNormal = Eigen::Vector3d::UnitZ();
    const Eigen::Vector3d contactSideNormal = -contactNormal;
    const double minPlane = box.min.z();
    const double zRange = std::max(1e-12, box.extent().z());
    const double xyExtent = std::max(box.extent().x(), box.extent().y());

    atlas::Patch best;
    int bestScore = -1;
    bool bestBffPass = false;
    bool bestRemesh = false;
    int planarCandidates = 0;

    for (double normalTol : {0.98, 0.95, 0.90, 0.80}) {
        for (double planeFrac : {1e-5, 1e-4, 1e-3, 0.005, 0.01, 0.02}) {
            std::vector<int> selected;
            const double planeTol = std::max(1e-8, planeFrac * zRange);
            for (int fi = 0; fi < mesh.F.rows(); ++fi) {
                const Eigen::Vector3d n = mesh.FN.row(fi).transpose();
                const Eigen::Vector3d c = faceCentroid(mesh, fi);
                if (n.dot(contactSideNormal) >= normalTol && c.z() <= minPlane + planeTol) {
                    selected.push_back(fi);
                }
            }
            planarCandidates += static_cast<int>(selected.size());
            for (const auto& component : selectedFaceComponents(mesh, selected)) {
                if (component.size() < 3) continue;
                auto patch = atlas::PatchBuilder::buildFromFaces(mesh, component);
                considerPatchCandidate(patch, chartMode, best, bestRemesh, bestScore, bestBffPass);
            }
        }
    }

    int seedFace = 0;
    double seedZ = std::numeric_limits<double>::infinity();
    Eigen::Vector3d seedCentroid = Eigen::Vector3d::Zero();
    for (int fi = 0; fi < mesh.F.rows(); ++fi) {
        const Eigen::Vector3d c = faceCentroid(mesh, fi);
        if (c.z() < seedZ) {
            seedZ = c.z();
            seedFace = fi;
            seedCentroid = c;
        }
    }

    for (double zFrac : {0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35}) {
        for (double radiusFrac : {0.12, 0.18, 0.25, 0.35, 0.5, 0.75, 1.0}) {
            std::vector<int> selected;
            const double zLimit = minPlane + zFrac * zRange;
            const double radius = radiusFrac * xyExtent;
            for (int fi = 0; fi < mesh.F.rows(); ++fi) {
                const Eigen::Vector3d c = faceCentroid(mesh, fi);
                if (c.z() <= zLimit && (c.head<2>() - seedCentroid.head<2>()).norm() <= radius) {
                    selected.push_back(fi);
                }
            }
            if (selected.empty()) selected.push_back(seedFace);

            for (const auto& component : selectedFaceComponents(mesh, selected)) {
                if (component.size() < 20) continue;
                auto patch = atlas::PatchBuilder::buildFromFaces(mesh, component);
                considerPatchCandidate(patch, chartMode, best, bestRemesh, bestScore, bestBffPass);
            }
        }
    }

    if (bestScore > 0) {
        outPatch = best;
        usedRemesh = bestRemesh;
        note = "plane-normal plus local disk search; bff_pass=" + boolText(bestBffPass) +
               "; remeshed=" + boolText(bestRemesh) +
               "; planar_candidate_face_evaluations=" + std::to_string(planarCandidates);
        return true;
    }

    std::vector<int> oneFace = {seedFace};
    outPatch = atlas::PatchBuilder::buildFromFaces(mesh, oneFace);
    usedRemesh = false;
    note = "single-face fallback";
    return outPatch.diskLike;
}

ChartQualityRow measureChartQuality(const std::string& scene,
                                    const atlas::ContactChart& chart,
                                    bool remeshed = false)
{
    ChartQualityRow row;
    row.scene = scene;
    row.chartMethod = chart.parameterizationMethod;
    row.patchFaces = chart.patch.localMesh.F.rows();
    row.uvFlippedTriangles = atlas::countFlippedTriangles(chart.UV, chart.Fuv);
    row.diskLike = chart.patch.diskLike;
    row.usedFallback = chart.usedFallbackParameterization;
    row.remeshed = remeshed;

    std::vector<double> angleDiffs;
    std::vector<double> logRatios;
    double maxCond = 0.0;

    for (int fi = 0; fi < chart.patch.localMesh.F.rows(); ++fi) {
        const Eigen::Vector3i f = chart.patch.localMesh.F.row(fi);
        const Eigen::Vector3d x0 = chart.patch.localMesh.V.row(f[0]);
        const Eigen::Vector3d x1 = chart.patch.localMesh.V.row(f[1]);
        const Eigen::Vector3d x2 = chart.patch.localMesh.V.row(f[2]);
        const Eigen::Vector2d u0 = chart.UV.row(f[0]);
        const Eigen::Vector2d u1 = chart.UV.row(f[1]);
        const Eigen::Vector2d u2 = chart.UV.row(f[2]);

        const double area3 = triangleArea3D(x0, x1, x2);
        const double area2 = triangleArea2D(u0, u1, u2);
        if (area3 <= 1e-15 || area2 <= 1e-15) continue;

        const std::array<double, 3> a3 = {
            angleBetween(x1 - x0, x2 - x0),
            angleBetween(x2 - x1, x0 - x1),
            angleBetween(x0 - x2, x1 - x2)
        };
        const std::array<double, 3> a2 = {
            angleBetween2D(u1 - u0, u2 - u0),
            angleBetween2D(u2 - u1, u0 - u1),
            angleBetween2D(u0 - u2, u1 - u2)
        };
        for (int k = 0; k < 3; ++k) {
            const double d = std::abs(a3[k] - a2[k]) * 180.0 / experiments::pi;
            angleDiffs.push_back(d);
            row.maxAngleDistortionDeg = std::max(row.maxAngleDistortionDeg, d);
        }
        logRatios.push_back(std::log(area2 / area3));

        Eigen::Matrix<double, 3, 2> Xedge;
        Xedge.col(0) = x1 - x0;
        Xedge.col(1) = x2 - x0;
        Eigen::Matrix2d Uedge;
        Uedge.col(0) = u1 - u0;
        Uedge.col(1) = u2 - u0;
        if (std::abs(Uedge.determinant()) > 1e-14) {
            const Eigen::Matrix<double, 3, 2> J = Xedge * Uedge.inverse();
            const Eigen::JacobiSVD<Eigen::Matrix<double, 3, 2>> svd(J);
            const auto sv = svd.singularValues();
            if (sv.size() == 2 && sv[1] > 1e-14) {
                maxCond = std::max(maxCond, sv[0] / sv[1]);
            }
        }
    }

    if (!angleDiffs.empty()) {
        double sum = 0.0;
        for (double d : angleDiffs) sum += d;
        row.meanAngleDistortionDeg = sum / static_cast<double>(angleDiffs.size());
    }
    if (!logRatios.empty()) {
        double mean = 0.0;
        for (double v : logRatios) mean += v;
        mean /= static_cast<double>(logRatios.size());
        double var = 0.0;
        for (double v : logRatios) var += (v - mean) * (v - mean);
        row.areaLogStd = std::sqrt(var / static_cast<double>(logRatios.size()));
    }
    row.maxJacobianCondition = maxCond;
    return row;
}

double clippedNegativeArea(const std::array<PolyVertex, 3>& tri)
{
    std::vector<PolyVertex> poly(tri.begin(), tri.end());
    std::vector<PolyVertex> out;
    for (int pass = 0; pass < 1; ++pass) {
        (void)pass;
        out.clear();
        for (int i = 0; i < static_cast<int>(poly.size()); ++i) {
            const auto& a = poly[i];
            const auto& b = poly[(i + 1) % poly.size()];
            const bool ain = a.g <= 0.0;
            const bool bin = b.g <= 0.0;
            if (ain && bin) {
                out.push_back(b);
            } else if (ain && !bin) {
                const double t = a.g / (a.g - b.g);
                out.push_back({(1.0 - t) * a.x + t * b.x, 0.0});
            } else if (!ain && bin) {
                const double t = a.g / (a.g - b.g);
                out.push_back({(1.0 - t) * a.x + t * b.x, 0.0});
                out.push_back(b);
            }
        }
        poly.swap(out);
    }

    if (poly.size() < 3) return 0.0;
    double area = 0.0;
    for (int i = 1; i + 1 < static_cast<int>(poly.size()); ++i) {
        area += triangleArea3D(poly[0].x, poly[i].x, poly[i + 1].x);
    }
    return area;
}

ReferenceResult subdividedReference(const core::Mesh& sourceMesh,
                                    const core::RigidPose& sourcePose,
                                    const sdf::ISDF& targetSdf,
                                    const core::RigidPose& targetPose,
                                    int subdivisions)
{
    ReferenceResult ref;
    double penetrationAreaSum = 0.0;
    double weightedPenetration = 0.0;
    const int n = std::max(1, subdivisions);

    auto bary = [n](int i, int j) {
        const double b = static_cast<double>(i) / static_cast<double>(n);
        const double c = static_cast<double>(j) / static_cast<double>(n);
        return Eigen::Vector3d(1.0 - b - c, b, c);
    };

    auto eval = [&](const Eigen::Vector3d& xLocal) {
        const Eigen::Vector3d xWorld = sourcePose.localToWorld(xLocal);
        const auto q = targetSdf.evalLocal(targetPose.worldToLocal(xWorld), false, false);
        return q.phi;
    };

    for (int fi = 0; fi < sourceMesh.F.rows(); ++fi) {
        const Eigen::Vector3i f = sourceMesh.F.row(fi);
        const Eigen::Vector3d v0 = sourceMesh.V.row(f[0]);
        const Eigen::Vector3d v1 = sourceMesh.V.row(f[1]);
        const Eigen::Vector3d v2 = sourceMesh.V.row(f[2]);
        auto point = [&](const Eigen::Vector3d& b) {
            return b[0] * v0 + b[1] * v1 + b[2] * v2;
        };

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n - i; ++j) {
                const Eigen::Vector3d b0 = bary(i, j);
                const Eigen::Vector3d b1 = bary(i + 1, j);
                const Eigen::Vector3d b2 = bary(i, j + 1);
                std::array<PolyVertex, 3> t0 = {{
                    {sourcePose.localToWorld(point(b0)), eval(point(b0))},
                    {sourcePose.localToWorld(point(b1)), eval(point(b1))},
                    {sourcePose.localToWorld(point(b2)), eval(point(b2))}
                }};
                const double a0 = clippedNegativeArea(t0);
                ref.area += a0;
                if (a0 > 0.0) {
                    const double p = std::max({0.0, -t0[0].g, -t0[1].g, -t0[2].g});
                    ref.maxPenetration = std::max(ref.maxPenetration, p);
                    penetrationAreaSum += a0;
                    weightedPenetration += a0 * p;
                }

                if (i + j < n - 1) {
                    const Eigen::Vector3d b3 = bary(i + 1, j + 1);
                    std::array<PolyVertex, 3> t1 = {{
                        {sourcePose.localToWorld(point(b1)), eval(point(b1))},
                        {sourcePose.localToWorld(point(b3)), eval(point(b3))},
                        {sourcePose.localToWorld(point(b2)), eval(point(b2))}
                    }};
                    const double a1 = clippedNegativeArea(t1);
                    ref.area += a1;
                    if (a1 > 0.0) {
                        const double p = std::max({0.0, -t1[0].g, -t1[1].g, -t1[2].g});
                        ref.maxPenetration = std::max(ref.maxPenetration, p);
                        penetrationAreaSum += a1;
                        weightedPenetration += a1 * p;
                    }
                }
            }
        }
    }

    if (penetrationAreaSum > 0.0) ref.meanPenetration = weightedPenetration / penetrationAreaSum;
    return ref;
}

Row makeRow(const std::string& scene,
            const std::string& sceneClass,
            const std::string& method,
            const core::Mesh& sourceMesh,
            const atlas::ContactChart& chart,
            const contact::ContactResult& result,
            double referenceArea,
            double referencePenetration,
            const std::string& referenceType,
            int sdfResolution = 0,
            double memoryMb = 0.0)
{
    Row row;
    row.scene = scene;
    row.sceneClass = sceneClass;
    row.method = method;
    row.meshFaces = sourceMesh.F.rows();
    row.sdfResolution = sdfResolution;
    row.contactArea = result.area;
    row.referenceArea = referenceArea;
    row.areaError = referenceArea > 1e-14 ? std::abs(result.area - referenceArea) / referenceArea : 0.0;
    row.maxPenetration = result.maxPenetration;
    row.referencePenetration = referencePenetration;
    row.penetrationError = referencePenetration > 1e-14 ? std::abs(result.maxPenetration - referencePenetration) / referencePenetration : 0.0;
    row.runtimeMs = result.timing.totalMs;
    row.integrationMs = result.timing.integrationMs;
    row.projectionMs = result.timing.projectionMs;
    row.otherMs = std::max(0.0, result.timing.totalMs - result.timing.integrationMs - result.timing.projectionMs);
    row.sdfQueries = result.sdfQueries;
    row.projectionCount = result.projectionCount;
    row.projectionFailures = result.projectionFailures;
    row.memoryMb = memoryMb;
    row.referenceType = referenceType;
    row.chartMethod = chart.parameterizationMethod;
    return row;
}

contact::ContactResult detectContact(const core::Mesh& sourceMesh,
                                     const atlas::ContactChart& chart,
                                     const core::RigidPose& sourcePose,
                                     const core::Mesh* targetMesh,
                                     const sdf::ISDF& targetSdf,
                                     const core::RigidPose& targetPose,
                                     bool projection,
                                     int sampleStride)
{
    atlas::ContactAtlas atlasData{{chart}};
    contact::RigidObject source{"source", &sourceMesh, &atlasData, nullptr, sourcePose};
    contact::RigidObject target{"target", targetMesh, nullptr, &targetSdf, targetPose};
    contact::ContactDetectorOptions options;
    options.refineProjection = projection;
    options.sampleStride = sampleStride;
    return contact::ContactDetector().detect(source, target, options);
}

Row runSpherePlaneScene(const std::string& scene,
                        double delta,
                        ChartMode chartMode,
                        bool projection,
                        int rings,
                        int segments)
{
    const double R = 1.0;
    auto sourceMesh = experiments::makeLowerHemispherePatch(R, rings, segments);
    auto chart = chartFromMesh(sourceMesh, chartMode);
    core::RigidPose sourcePose;
    sourcePose.t = Eigen::Vector3d(0, 0, -sourceMesh.bounds().min.z() - delta);
    auto targetMesh = experiments::makePlane(3.0);
    sdf::PlaneSDF plane(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
    const auto result = detectContact(sourceMesh, chart, sourcePose, &targetMesh, plane, core::RigidPose{}, projection, 14);
    const double analyticArea = 2.0 * experiments::pi * R * delta;
    return makeRow(scene, "analytic", projection ? "ours_sdf_projection_" + chartSuffix(chart)
                                                  : "ours_sdf_pullback_only_" + chartSuffix(chart),
                   sourceMesh, chart, result, analyticArea, delta, "closed_form_sphere_plane");
}

Row runTiltedSpherePlaneScene()
{
    const double R = 1.0;
    const double delta = 0.16;
    auto sourceMesh = experiments::makeLowerHemispherePatch(R, 42, 120);
    auto chart = chartFromMesh(sourceMesh, ChartMode::OfficialBff);
    core::RigidPose sourcePose;
    const Eigen::Vector3d n = Eigen::Vector3d(0.28, 0.12, 1.0).normalized();
    sourcePose.t = n * (R - delta);
    auto targetMesh = experiments::makePlane(3.0);
    sdf::PlaneSDF plane(Eigen::Vector3d::Zero(), n);
    const auto result = detectContact(sourceMesh, chart, sourcePose, &targetMesh, plane, core::RigidPose{}, true, 14);
    const double analyticArea = 2.0 * experiments::pi * R * delta;
    return makeRow("sphere_tilted_plane", "analytic", "ours_sdf_projection_" + chartSuffix(chart),
                   sourceMesh, chart, result, analyticArea, delta, "closed_form_tilted_sphere_plane");
}

Row runSphereSphereScene(const std::string& scene,
                         double sourceRadius,
                         double targetRadius,
                         double delta,
                         const sdf::ISDF& targetSdf,
                         const core::Mesh* targetMesh,
                         const std::string& methodPrefix,
                         ChartMode chartMode,
                         int sdfResolution = 0,
                         double memoryMb = 0.0)
{
    const double d = sourceRadius + targetRadius - delta;
    auto sourceMesh = experiments::makeLowerHemispherePatch(sourceRadius, 40, 96);
    auto chart = chartFromMesh(sourceMesh, chartMode);
    core::RigidPose targetPose;
    targetPose.t = Eigen::Vector3d(0, 0, -d);
    const auto result = detectContact(sourceMesh, chart, core::RigidPose{}, targetMesh, targetSdf, targetPose, targetMesh != nullptr, 16);
    const double z0 = (targetRadius * targetRadius - sourceRadius * sourceRadius - d * d) / (2.0 * d);
    const double capHeight = std::max(0.0, std::min(2.0 * sourceRadius, z0 + sourceRadius));
    const double analyticArea = 2.0 * experiments::pi * sourceRadius * capHeight;
    return makeRow(scene, "analytic", methodPrefix + "_" + chartSuffix(chart),
                   sourceMesh, chart, result, analyticArea, delta, "closed_form_sphere_sphere",
                   sdfResolution, memoryMb);
}

Row runEllipsoidPlaneReferenceScene()
{
    const double delta = 0.18;
    auto sourceMesh = experiments::makeEllipsoidLowerPatch(Eigen::Vector3d(1.2, 0.75, 1.0), 38, 108);
    auto chart = chartFromMesh(sourceMesh, ChartMode::OfficialBff);
    core::RigidPose sourcePose;
    sourcePose.t = Eigen::Vector3d(0, 0, -sourceMesh.bounds().min.z() - delta);
    auto targetMesh = experiments::makePlane(4.0);
    sdf::PlaneSDF plane(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
    const auto result = detectContact(sourceMesh, chart, sourcePose, &targetMesh, plane, core::RigidPose{}, true, 12);
    const auto ref = subdividedReference(sourceMesh, sourcePose, plane, core::RigidPose{}, 6);
    return makeRow("ellipsoid_plane_reference", "analytic", "ours_sdf_projection_" + chartSuffix(chart),
                   sourceMesh, chart, result, ref.area, delta, "subdivided_surface_reference");
}

Row runProceduralPlaneReferenceScene(const std::string& scene, const core::Mesh& sourceMesh, double delta, ChartMode chartMode)
{
    auto chart = chartFromMesh(sourceMesh, chartMode);
    core::RigidPose sourcePose;
    sourcePose.t = Eigen::Vector3d(0, 0, -sourceMesh.bounds().min.z() - delta);
    auto targetMesh = experiments::makePlane(4.0);
    sdf::PlaneSDF plane(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
    const auto result = detectContact(sourceMesh, chart, sourcePose, &targetMesh, plane, core::RigidPose{}, true, 12);
    const auto ref = subdividedReference(sourceMesh, sourcePose, plane, core::RigidPose{}, 5);
    return makeRow(scene, "stress", "ours_sdf_projection_" + chartSuffix(chart),
                   sourceMesh, chart, result, ref.area, delta, "subdivided_surface_reference");
}

bool loadBffOfficialObj(const std::string& name, core::Mesh& mesh)
{
    const std::filesystem::path path = std::filesystem::path("bff_official") / "input" / name;
    std::string error;
    return std::filesystem::exists(path) && io::readObj(path.string(), mesh, &error);
}

bool loadMechanicalGear(core::Mesh& mesh)
{
    const std::filesystem::path path = std::filesystem::path("data") / "meshes" / "involute_gear_teeth_16_angle_20_cc0.stl";
    std::string error;
    return std::filesystem::exists(path) && io::readStl(path.string(), mesh, &error);
}

bool prepareRealPatch(const std::string& scene,
                      core::Mesh fullMesh,
                      ChartMode chartMode,
                      core::Mesh& repairedFull,
                      atlas::Patch& patch,
                      atlas::ContactChart& chart,
                      AssetStat& stat)
{
    repairedFull = normalizeMesh(fullMesh, 1.4);
    core::MeshRepairReport repairReport;
    repairedFull = core::repairMesh(repairedFull, {}, &repairReport);
    bool usedRemesh = false;
    std::string patchNote;
    if (!findContactDiskPatch(repairedFull, chartMode, patch, usedRemesh, patchNote)) {
        stat.scene = scene;
        stat.note = "no disk-like contact patch found";
        return false;
    }
    chart = buildChartFromPatch(patch, chartMode);
    stat.scene = scene;
    stat.fullFaces = repairedFull.F.rows();
    stat.patchFaces = patch.localMesh.F.rows();
    stat.patchVertices = patch.localMesh.V.rows();
    stat.boundaryLoops = static_cast<int>(patch.localMesh.boundaryLoops().size());
    stat.eulerCharacteristic = patch.localMesh.eulerCharacteristic();
    stat.diskLike = patch.diskLike;
    stat.chartMethod = chart.parameterizationMethod;
    stat.remeshed = usedRemesh;
    stat.repairRemovedDegenerate = repairReport.removedDegenerateFaces;
    stat.repairRemovedDuplicate = repairReport.removedDuplicateFaces;
    stat.repairFlippedFaces = repairReport.flippedFaces;
    stat.note = patchNote;
    return true;
}

Row runRealPlaneScene(const std::string& scene,
                      const atlas::Patch& patch,
                      const atlas::ContactChart& chart,
                      double delta,
                      AssetStat& stat)
{
    core::RigidPose sourcePose;
    sourcePose.t = Eigen::Vector3d(0, 0, -patch.localMesh.bounds().min.z() - delta);
    auto targetMesh = experiments::makePlane(3.0);
    sdf::PlaneSDF plane(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
    const auto result = detectContact(patch.localMesh, chart, sourcePose, &targetMesh, plane, core::RigidPose{}, true, 8);
    const auto ref = subdividedReference(patch.localMesh, sourcePose, plane, core::RigidPose{}, 6);
    stat.projectionFailures += result.projectionFailures;
    return makeRow(scene, "real_mesh", "ours_sdf_projection_" + chartSuffix(chart),
                   patch.localMesh, chart, result, ref.area, delta, "subdivided_surface_reference");
}

Row runRealSphereScene(const std::string& scene,
                       const atlas::Patch& patch,
                       const atlas::ContactChart& chart,
                       double sphereRadius,
                       double delta,
                       AssetStat& stat)
{
    const auto box = patch.localMesh.bounds();
    core::RigidPose targetPose;
    targetPose.t = Eigen::Vector3d(0, 0, box.min.z() - sphereRadius + delta);
    auto targetMesh = experiments::makeSphere(sphereRadius, 36, 96);
    sdf::SphereSDF sphere(Eigen::Vector3d::Zero(), sphereRadius);
    const auto result = detectContact(patch.localMesh, chart, core::RigidPose{}, &targetMesh, sphere, targetPose, true, 1);
    const auto ref = subdividedReference(patch.localMesh, core::RigidPose{}, sphere, targetPose, 7);
    stat.projectionFailures += result.projectionFailures;
    return makeRow(scene, "real_mesh", "ours_sdf_projection_" + chartSuffix(chart),
                   patch.localMesh, chart, result, ref.area, ref.maxPenetration, "subdivided_surface_reference");
}

Row runRealTiltedPlaneScene(const std::string& scene,
                            const atlas::Patch& patch,
                            const atlas::ContactChart& chart,
                            double delta,
                            AssetStat& stat)
{
    const Eigen::Vector3d n = Eigen::Vector3d(0.22, -0.18, 1.0).normalized();
    core::RigidPose sourcePose;
    double minProjection = std::numeric_limits<double>::infinity();
    for (int i = 0; i < patch.localMesh.V.rows(); ++i) {
        minProjection = std::min(minProjection, n.dot(patch.localMesh.V.row(i).transpose()));
    }
    sourcePose.t = n * (-delta - minProjection);
    auto targetMesh = experiments::makePlane(3.0);
    sdf::PlaneSDF plane(Eigen::Vector3d::Zero(), n);
    const auto result = detectContact(patch.localMesh, chart, sourcePose, &targetMesh, plane, core::RigidPose{}, true, 8);
    const auto ref = subdividedReference(patch.localMesh, sourcePose, plane, core::RigidPose{}, 7);
    stat.projectionFailures += result.projectionFailures;
    return makeRow(scene, "real_mesh", "ours_sdf_projection_" + chartSuffix(chart),
                   patch.localMesh, chart, result, ref.area, ref.maxPenetration, "subdivided_surface_reference");
}

atlas::ContactAtlas splitHemisphereAtlas(const core::Mesh& mesh, ChartMode chartMode)
{
    std::vector<int> left;
    std::vector<int> right;
    for (int fi = 0; fi < mesh.F.rows(); ++fi) {
        const Eigen::Vector3d c = faceCentroid(mesh, fi);
        (c.x() <= 0.0 ? left : right).push_back(fi);
    }

    atlas::ContactAtlas atlasData;
    int id = 0;
    for (const auto& group : {left, right}) {
        if (group.empty()) continue;
        auto patch = atlas::PatchBuilder::buildFromFaces(mesh, group);
        if (!patch.diskLike) continue;
        auto chart = buildChartFromPatch(patch, chartMode);
        chart.id = id++;
        atlasData.charts.push_back(chart);
    }
    return atlasData;
}

Row runSeamSpherePlaneScene()
{
    const double R = 1.0;
    const double delta = 0.2;
    auto mesh = experiments::makeLowerHemispherePatch(R, 36, 108);
    auto atlasData = splitHemisphereAtlas(mesh, ChartMode::OfficialBff);
    core::RigidPose sourcePose;
    sourcePose.t = Eigen::Vector3d(0, 0, -mesh.bounds().min.z() - delta);
    auto targetMesh = experiments::makePlane(3.0);
    sdf::PlaneSDF plane(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
    contact::RigidObject source{"split_source", &mesh, &atlasData, nullptr, sourcePose};
    contact::RigidObject target{"plane", &targetMesh, nullptr, &plane, {}};
    contact::ContactDetectorOptions options;
    options.sampleStride = 12;
    const auto result = contact::ContactDetector().detect(source, target, options);
    atlas::ContactChart representative = atlasData.charts.empty() ? chartFromMesh(mesh, ChartMode::PlanarFallback) : atlasData.charts.front();
    const double analyticArea = 2.0 * experiments::pi * R * delta;
    return makeRow("seam_two_chart_sphere_plane", "stress", "ours_sdf_projection_multichart_bff",
                   mesh, representative, result, analyticArea, delta, "closed_form_sphere_plane_multichart");
}

Row runSeamSphereSphereScene()
{
    const double R = 1.0;
    const double delta = 0.16;
    const double d = 2.0 * R - delta;
    auto mesh = experiments::makeLowerHemispherePatch(R, 36, 108);
    auto atlasData = splitHemisphereAtlas(mesh, ChartMode::OfficialBff);
    auto targetMesh = experiments::makeSphere(R, 36, 96);
    sdf::SphereSDF sphere(Eigen::Vector3d::Zero(), R);
    core::RigidPose targetPose;
    targetPose.t = Eigen::Vector3d(0, 0, -d);
    contact::RigidObject source{"split_source", &mesh, &atlasData, nullptr, {}};
    contact::RigidObject target{"sphere", &targetMesh, nullptr, &sphere, targetPose};
    contact::ContactDetectorOptions options;
    options.sampleStride = 12;
    const auto result = contact::ContactDetector().detect(source, target, options);
    atlas::ContactChart representative = atlasData.charts.empty() ? chartFromMesh(mesh, ChartMode::PlanarFallback) : atlasData.charts.front();
    const double analyticArea = experiments::pi * R * delta;
    return makeRow("seam_two_chart_sphere_sphere", "stress", "ours_sdf_projection_multichart_bff",
                   mesh, representative, result, analyticArea, delta, "closed_form_sphere_sphere_multichart");
}

void writeRows(const std::string& path, const std::vector<Row>& rows)
{
    std::ofstream out(path);
    out << "scene,scene_class,method,mesh_faces,sdf_resolution,contact_area,reference_area,area_error,"
           "max_penetration,reference_penetration,penetration_error,normal_error_deg,contour_hausdorff,"
           "runtime_ms,integration_ms,projection_ms,other_ms,sdf_queries,projection_count,projection_failures,"
           "memory_mb,reference_type,chart_method\n";
    for (const auto& r : rows) {
        out << r.scene << "," << r.sceneClass << "," << r.method << "," << r.meshFaces << "," << r.sdfResolution << ","
            << r.contactArea << "," << r.referenceArea << "," << r.areaError << ","
            << r.maxPenetration << "," << r.referencePenetration << "," << r.penetrationError << ","
            << r.normalErrorDeg << "," << r.contourHausdorff << "," << r.runtimeMs << ","
            << r.integrationMs << "," << r.projectionMs << "," << r.otherMs << ","
            << r.sdfQueries << "," << r.projectionCount << "," << r.projectionFailures << ","
            << r.memoryMb << "," << r.referenceType << "," << r.chartMethod << "\n";
    }
}

void writeAssetStats(const std::string& path, const std::vector<AssetStat>& stats)
{
    std::ofstream out(path);
    out << "scene,full_faces,patch_faces,patch_vertices,boundary_loops,euler_characteristic,disk_like,"
           "chart_method,remeshed,repair_removed_degenerate,repair_removed_duplicate,repair_flipped_faces,"
           "projection_failures,note\n";
    for (const auto& s : stats) {
        out << s.scene << "," << s.fullFaces << "," << s.patchFaces << "," << s.patchVertices << ","
            << s.boundaryLoops << "," << s.eulerCharacteristic << "," << boolText(s.diskLike) << ","
            << s.chartMethod << "," << boolText(s.remeshed) << "," << s.repairRemovedDegenerate << ","
            << s.repairRemovedDuplicate << "," << s.repairFlippedFaces << "," << s.projectionFailures << ","
            << s.note << "\n";
    }
}

void writeChartQuality(const std::string& path, const std::vector<ChartQualityRow>& rows)
{
    std::ofstream out(path);
    out << "scene,chart_method,patch_faces,uv_flipped_triangles,mean_angle_distortion_deg,"
           "max_angle_distortion_deg,area_log_std,max_jacobian_condition,disk_like,used_fallback,remeshed\n";
    for (const auto& r : rows) {
        out << r.scene << "," << r.chartMethod << "," << r.patchFaces << "," << r.uvFlippedTriangles << ","
            << r.meanAngleDistortionDeg << "," << r.maxAngleDistortionDeg << "," << r.areaLogStd << ","
            << r.maxJacobianCondition << "," << boolText(r.diskLike) << "," << boolText(r.usedFallback) << ","
            << boolText(r.remeshed) << "\n";
    }
}

} // namespace

int main()
{
    const std::string dir = "results/benchmarks";
    io::ensureDirectory(dir);

    std::vector<Row> rows;
    std::vector<AssetStat> assetStats;
    std::vector<ChartQualityRow> qualityRows;

    rows.push_back(runSpherePlaneScene("sphere_plane_shallow", 0.05, ChartMode::OfficialBff, true, 44, 120));
    rows.push_back(runSpherePlaneScene("sphere_plane_medium", 0.20, ChartMode::OfficialBff, true, 44, 120));
    rows.push_back(runSpherePlaneScene("sphere_plane_medium", 0.20, ChartMode::OfficialBff, false, 44, 120));
    rows.push_back(runSpherePlaneScene("sphere_plane_medium", 0.20, ChartMode::PlanarFallback, true, 44, 120));
    rows.push_back(runTiltedSpherePlaneScene());

    auto targetSphereMesh = experiments::makeSphere(1.0, 40, 96);
    sdf::SphereSDF targetSphere(Eigen::Vector3d::Zero(), 1.0);
    rows.push_back(runSphereSphereScene("sphere_sphere_equal", 1.0, 1.0, 0.18, targetSphere, &targetSphereMesh,
                                        "ours_sdf_projection", ChartMode::OfficialBff));
    rows.push_back(runSphereSphereScene("sphere_sphere_equal", 1.0, 1.0, 0.18, sdf::MeshSDF(targetSphereMesh),
                                        &targetSphereMesh, "meshsdf_bvh_exact_triangle_closest_point", ChartMode::OfficialBff));
    rows.push_back(runSphereSphereScene("sphere_sphere_equal", 1.0, 1.0, 0.18, sdf::MeshSDF(targetSphereMesh),
                                        nullptr, "meshsdf_bvhsdf_pullback_only", ChartMode::OfficialBff));

    core::AABB gridBounds;
    gridBounds.expand(Eigen::Vector3d(-1.5, -1.5, -1.5));
    gridBounds.expand(Eigen::Vector3d(1.5, 1.5, 1.5));
    const auto gridSphere = sdf::GridSDF::sampleFrom(targetSphere, gridBounds, 32, 0.0);
    rows.push_back(runSphereSphereScene("sphere_sphere_equal", 1.0, 1.0, 0.18, gridSphere, nullptr,
                                        "grid_sdf_pullback_only", ChartMode::OfficialBff, 32,
                                        static_cast<double>(gridSphere.memoryBytes()) / (1024.0 * 1024.0)));

    auto targetLargeMesh = experiments::makeSphere(1.35, 40, 96);
    sdf::SphereSDF targetLargeSphere(Eigen::Vector3d::Zero(), 1.35);
    rows.push_back(runSphereSphereScene("sphere_sphere_unequal", 1.0, 1.35, 0.16, targetLargeSphere, &targetLargeMesh,
                                        "ours_sdf_projection", ChartMode::OfficialBff));
    rows.push_back(runEllipsoidPlaneReferenceScene());

    rows.push_back(runProceduralPlaneReferenceScene("torus_local_patch_plane_stress",
                                                    experiments::makeTorusPatch(0.75, 0.28, 96, 18), 0.32,
                                                    ChartMode::OfficialBff));
    rows.push_back(runSeamSpherePlaneScene());
    rows.push_back(runSeamSphereSphereScene());

    core::Mesh bunny;
    if (loadBffOfficialObj("bunny.obj", bunny)) {
        core::Mesh repaired;
        atlas::Patch patch;
        atlas::ContactChart chart;
        AssetStat stat;
        if (prepareRealPatch("bunny_real_contact_patch", bunny, ChartMode::OfficialBff, repaired, patch, chart, stat)) {
            rows.push_back(runRealPlaneScene("bunny_plane_real", patch, chart, 0.12, stat));
            rows.push_back(runRealSphereScene("bunny_sphere_real", patch, chart, 0.72, 0.10, stat));
            rows.push_back(runRealTiltedPlaneScene("bunny_tilted_plane_real", patch, chart, 0.10, stat));
            qualityRows.push_back(measureChartQuality("bunny_real_contact_patch", chart, stat.remeshed));
            assetStats.push_back(stat);
        } else {
            assetStats.push_back(stat);
        }
    } else {
        assetStats.push_back({"bunny_real_contact_patch", 0, 0, 0, 0, 0, false, "missing", false, 0, 0, 0, 0,
                              "missing bff_official/input/bunny.obj"});
    }

    core::Mesh gear;
    if (loadMechanicalGear(gear)) {
        core::Mesh repaired;
        atlas::Patch patch;
        atlas::ContactChart chart;
        AssetStat stat;
        if (prepareRealPatch("mechanical_gear_real_contact_patch", gear, ChartMode::OfficialBff, repaired, patch, chart, stat)) {
            rows.push_back(runRealPlaneScene("mechanical_gear_plane_real", patch, chart, 0.08, stat));
            rows.push_back(runRealSphereScene("mechanical_gear_sphere_real", patch, chart, 0.70, 0.07, stat));
            qualityRows.push_back(measureChartQuality("mechanical_gear_real_contact_patch", chart, stat.remeshed));
            assetStats.push_back(stat);
        } else {
            assetStats.push_back(stat);
        }
    } else {
        assetStats.push_back({"mechanical_gear_real_contact_patch", 0, 0, 0, 0, 0, false, "missing", false, 0, 0, 0, 0,
                              "missing tracked CC0 gear STL"});
    }

    const std::vector<std::pair<std::string, core::Mesh>> curvedPatches = {
        {"spherical_cap_30", makeSphericalCapPatch(1.0, 30.0, 24, 96)},
        {"spherical_cap_90", makeSphericalCapPatch(1.0, 90.0, 32, 112)},
        {"spherical_cap_135", makeSphericalCapPatch(1.0, 135.0, 36, 128)},
        {"ellipsoid_cap_curved", makeEllipsoidCapPatch(Eigen::Vector3d(1.3, 0.7, 1.0), 110.0, 34, 120)}
    };
    for (const auto& item : curvedPatches) {
        const auto patch = atlas::PatchBuilder::buildWholeMeshPatch(item.second);
        qualityRows.push_back(measureChartQuality(item.first, buildChartFromPatch(patch, ChartMode::OfficialBff)));
        qualityRows.push_back(measureChartQuality(item.first, buildChartFromPatch(patch, ChartMode::PlanarFallback)));
    }

    writeRows(dir + "/accuracy.csv", rows);
    writeRows(dir + "/runtime.csv", rows);
    writeAssetStats(dir + "/asset_stats.csv", assetStats);
    writeChartQuality(dir + "/chart_quality.csv", qualityRows);
    writeChartQuality(dir + "/bff_vs_planar_ablation.csv", qualityRows);

    sdf::SphereSDF sphere(Eigen::Vector3d::Zero(), 1.0);
    std::ofstream sdfStudy(dir + "/sdf_resolution_study.csv");
    sdfStudy << "scene,method,sdf_resolution,gap_error,memory_mb\n";
    const Eigen::Vector3d p(1.08, 0.31, 0.17);
    const double exact = sphere.evalLocal(p, false, false).phi;
    auto sphereMesh = experiments::makeSphere(1.0, 42, 96);
    for (int res : {10, 14, 20, 32, 48, 64}) {
        const auto grid = sdf::GridSDF::sampleFrom(sphere, gridBounds, res, 0.0);
        const auto q = grid.evalLocal(p, true, true);
        contact::ChartProjector projector;
        const auto pr = projector.project(p, sphereMesh, core::RigidPose{}, &q);
        sdfStudy << "sphere,grid_sdf," << res << "," << std::abs(q.phi - exact) << ","
                 << static_cast<double>(grid.memoryBytes()) / (1024.0 * 1024.0) << "\n";
        sdfStudy << "sphere,grid_sdf_projection," << res << "," << std::abs(pr.signedGap - exact) << ","
                 << static_cast<double>(grid.memoryBytes()) / (1024.0 * 1024.0) << "\n";
    }

    std::ofstream meshStudy(dir + "/mesh_resolution_study.csv");
    meshStudy << "scene,method,rings,segments,mesh_faces,area_error,runtime_ms\n";
    for (int rings : {12, 16, 28, 44, 56}) {
        const int segments = rings * 3;
        const auto row = runSpherePlaneScene("sphere_plane_mesh_sweep", 0.2, ChartMode::OfficialBff, true, rings, segments);
        meshStudy << "sphere_plane,ours_sdf_projection_bff," << rings << "," << segments << ","
                  << row.meshFaces << "," << row.areaError << "," << row.runtimeMs << "\n";
    }

    std::ofstream ablation(dir + "/ablation.csv");
    ablation << "question,variant,answer_metric,value,notes\n";
    const auto noProjection = rows[2];
    const auto withProjection = rows[1];
    ablation << "projection reduces low-resolution GridSDF error,no_projection,gap_error,see_sdf_resolution_study,pure grid row\n";
    ablation << "projection reduces low-resolution GridSDF error,with_projection,gap_error,see_sdf_resolution_study,projection row; high-resolution grid can beat the linear target-mesh projection floor\n";
    ablation << "SDF reduces expensive projections,no_sdf_projection_count,projection_count," << withProjection.meshFaces << ",brute-force would project every source triangle\n";
    ablation << "SDF reduces expensive projections,with_sdf_projection_count,projection_count," << withProjection.projectionCount << ",narrow-band projection samples\n";
    ablation << "BFF chart vs fallback planar,official_bff,quality_metrics,see_bff_vs_planar_ablation,curved patch chart quality rows\n";
    ablation << "BFF chart vs fallback planar,fallback_planar,quality_metrics,see_bff_vs_planar_ablation,curved patch chart quality rows\n";
    ablation << "BVH exact baseline,meshsdf_bvh_exact_triangle_closest_point,area_error," << rows[6].areaError << ",sphere_sphere closed target mesh baseline with target-mesh projection\n";
    ablation << "MeshSDF/BVHSDF baseline,meshsdf_bvhsdf_pullback_only,area_error," << rows[7].areaError << ",sphere_sphere BVH SDF pullback without projection\n";
    ablation << "GridSDF baseline,grid_sdf_pullback_only,area_error," << rows[8].areaError << ",sphere_sphere pure grid baseline\n";
    ablation << "mesh resolution affects area,coarse_to_fine,area_error,see_mesh_resolution_study,computed rows\n";

    std::ofstream report(dir + "/validation_report.md");
    const int analyticRows = static_cast<int>(std::count_if(rows.begin(), rows.end(), [](const Row& r) { return r.sceneClass == "analytic"; }));
    const int realRows = static_cast<int>(std::count_if(rows.begin(), rows.end(), [](const Row& r) { return r.sceneClass == "real_mesh"; }));
    const int stressRows = static_cast<int>(std::count_if(rows.begin(), rows.end(), [](const Row& r) { return r.sceneClass == "stress"; }));
    std::set<std::string> uniqueScenes;
    std::set<std::string> uniqueAnalyticScenes;
    std::set<std::string> uniqueRealScenes;
    std::set<std::string> uniqueStressScenes;
    for (const auto& row : rows) {
        uniqueScenes.insert(row.scene);
        if (row.sceneClass == "analytic") uniqueAnalyticScenes.insert(row.scene);
        if (row.sceneClass == "real_mesh") uniqueRealScenes.insert(row.scene);
        if (row.sceneClass == "stress") uniqueStressScenes.insert(row.scene);
    }
    report << "# Benchmark Validation Report\n\n"
           << "- Generated `accuracy.csv`, `runtime.csv`, `asset_stats.csv`, `chart_quality.csv`, `sdf_resolution_study.csv`, `mesh_resolution_study.csv`, `ablation.csv`, and `bff_vs_planar_ablation.csv`.\n"
           << "- Unique main scenes: " << uniqueScenes.size() << " (`analytic=" << uniqueAnalyticScenes.size()
           << "`, `real_mesh=" << uniqueRealScenes.size() << "`, `stress=" << uniqueStressScenes.size() << "`); benchmark rows including baselines: "
           << rows.size() << " (`analytic=" << analyticRows
           << "`, `real_mesh=" << realRows << "`, `stress=" << stressRows << "`).\n"
           << "- Baselines represented: BVH exact target-mesh projection, pure GridSDF, MeshSDF/BVHSDF pullback, planar fallback.\n"
           << "- SDF resolution study is interpreted as coarse-grid error reduction with a high-resolution linear projection floor.\n"
           << "- Complex mesh rows use `subdivided_surface_reference`; area errors are no longer default-filled zeros.\n"
           << "- `scripts/generate_paper_figures.py` writes pipeline, convergence, chart-quality, GridSDF-refinement, real-asset, and runtime SVG figures under `results/figures/`.\n"
           << "- Sphere-plane medium BFF area error: " << rows[1].areaError << "\n"
           << "- Sphere-sphere equal analytic area error: " << rows[5].areaError << "\n"
           << "- Sphere-sphere MeshSDF/BVH baseline area error: " << rows[6].areaError << "\n"
           << "- Sphere-sphere MeshSDF/BVHSDF pullback area error: " << rows[7].areaError << "\n"
           << "- Sphere-sphere pure GridSDF baseline area error: " << rows[8].areaError << "\n"
           << "- Curved BFF-vs-planar chart quality rows: " << qualityRows.size() << "\n";
    for (const auto& stat : assetStats) {
        report << "- " << stat.scene << ": full_faces=" << stat.fullFaces
               << ", patch_faces=" << stat.patchFaces
               << ", disk_like=" << boolText(stat.diskLike)
               << ", chart=" << stat.chartMethod
               << ", remeshed=" << boolText(stat.remeshed)
               << ", projection_failures=" << stat.projectionFailures
               << ", repair_removed_degenerate=" << stat.repairRemovedDegenerate
               << ", repair_flipped_faces=" << stat.repairFlippedFaces
               << ", note=" << stat.note << "\n";
    }
    report << "- Mechanical gear asset source: koppi/involute-gear-collection, CC0-1.0, `angle-20/teeth-16.stl`.\n";

    std::cout << "benchmarks written to " << dir << "\n";
    return 0;
}
