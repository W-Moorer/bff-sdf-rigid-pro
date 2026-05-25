#include "atlas/BFFAdapter.h"
#include "atlas/ContactAtlas.h"
#include "atlas/PatchBuilder.h"
#include "contact/ChartProjector.h"
#include "contact/ContactDetector.h"
#include "experiments/ProceduralMeshes.h"
#include "io/DebugExport.h"
#include "io/MeshIO.h"
#include "sdf/AnalyticSDF.h"
#include "sdf/GridSDF.h"
#include "sdf/MeshSDF.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <set>
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
    std::string method;
    int meshFaces = 0;
    int sdfResolution = 0;
    double contactArea = 0.0;
    double areaError = 0.0;
    double maxPenetration = 0.0;
    double penetrationError = 0.0;
    double normalErrorDeg = 0.0;
    double contourHausdorff = 0.0;
    double runtimeMs = 0.0;
    int sdfQueries = 0;
    int projectionCount = 0;
    int projectionFailures = 0;
    double memoryMb = 0.0;
};

std::string chartSuffix(const atlas::ContactChart& chart)
{
    return chart.parameterizationMethod == "official-bff" ? "bff" : "planar";
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
    const double scale = targetExtent / box.extent().maxCoeff();
    for (int i = 0; i < mesh.V.rows(); ++i) {
        Eigen::Vector3d p = mesh.V.row(i);
        p = (p - center) * scale;
        mesh.V.row(i) = p.transpose();
    }
    mesh.computeNormals();
    return mesh;
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

bool findBottomDiskPatch(const core::Mesh& mesh, atlas::Patch& outPatch)
{
    const auto box = mesh.bounds();
    const double zMin = box.min.z();
    const double zRange = std::max(1e-12, box.extent().z());
    const double xyExtent = std::max(box.extent().x(), box.extent().y());

    int seedFace = 0;
    double seedZ = std::numeric_limits<double>::infinity();
    Eigen::Vector3d seedCentroid = Eigen::Vector3d::Zero();
    for (int fi = 0; fi < mesh.F.rows(); ++fi) {
        const Eigen::Vector3i f = mesh.F.row(fi);
        const Eigen::Vector3d c = (mesh.V.row(f[0]) + mesh.V.row(f[1]) + mesh.V.row(f[2])) / 3.0;
        if (c.z() < seedZ) {
            seedZ = c.z();
            seedFace = fi;
            seedCentroid = c;
        }
    }

    atlas::Patch best;
    int bestFaces = 0;
    for (double zFrac : {0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35}) {
        for (double radiusFrac : {0.12, 0.18, 0.25, 0.35, 0.5, 0.75, 1.0}) {
            std::vector<int> selected;
            const double zLimit = zMin + zFrac * zRange;
            const double radius = radiusFrac * xyExtent;
            for (int fi = 0; fi < mesh.F.rows(); ++fi) {
                const Eigen::Vector3i f = mesh.F.row(fi);
                const Eigen::Vector3d c = (mesh.V.row(f[0]) + mesh.V.row(f[1]) + mesh.V.row(f[2])) / 3.0;
                if (c.z() <= zLimit && (c.head<2>() - seedCentroid.head<2>()).norm() <= radius) {
                    selected.push_back(fi);
                }
            }
            if (selected.empty()) selected.push_back(seedFace);

            for (const auto& component : selectedFaceComponents(mesh, selected)) {
                if (component.size() < 20) continue;
                auto patch = atlas::PatchBuilder::buildFromFaces(mesh, component);
                if (patch.diskLike && static_cast<int>(component.size()) > bestFaces) {
                    best = patch;
                    bestFaces = static_cast<int>(component.size());
                }
            }
        }
    }

    if (bestFaces > 0) {
        outPatch = best;
        return true;
    }

    std::vector<int> oneFace = {seedFace};
    outPatch = atlas::PatchBuilder::buildFromFaces(mesh, oneFace);
    return outPatch.diskLike;
}

contact::ContactResult detectPlaneContact(const core::Mesh& sourceMesh,
                                          const atlas::ContactChart& chart,
                                          double delta,
                                          bool projection,
                                          int sampleStride)
{
    auto targetMesh = experiments::makePlane(3.0);
    sdf::PlaneSDF plane(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
    atlas::ContactAtlas atlasData{{chart}};

    contact::RigidObject source{"source", &sourceMesh, &atlasData, nullptr, {}};
    source.pose.t = Eigen::Vector3d(0, 0, -sourceMesh.bounds().min.z() - delta);
    contact::RigidObject target{"plane", &targetMesh, nullptr, &plane, {}};
    contact::ContactDetectorOptions options;
    options.refineProjection = projection;
    options.sampleStride = sampleStride;
    return contact::ContactDetector().detect(source, target, options);
}

Row makeRow(const std::string& scene,
            const std::string& method,
            const core::Mesh& sourceMesh,
            const contact::ContactResult& result,
            double expectedPenetration,
            double analyticArea = 0.0)
{
    Row row;
    row.scene = scene;
    row.method = method;
    row.meshFaces = sourceMesh.F.rows();
    row.contactArea = result.area;
    row.areaError = analyticArea > 0.0 ? std::abs(result.area - analyticArea) / analyticArea : 0.0;
    row.maxPenetration = result.maxPenetration;
    row.penetrationError = expectedPenetration > 0.0 ? std::abs(result.maxPenetration - expectedPenetration) / expectedPenetration : 0.0;
    row.runtimeMs = result.timing.totalMs;
    row.sdfQueries = result.sdfQueries;
    row.projectionCount = result.projectionCount;
    row.projectionFailures = result.projectionFailures;
    return row;
}

Row runSpherePlane(int rings, int segments, bool projection, ChartMode chartMode, const std::string& methodPrefix)
{
    const double R = 1.0;
    const double delta = 0.2;
    auto sourceMesh = experiments::makeLowerHemispherePatch(R, rings, segments);
    auto chart = chartFromMesh(sourceMesh, chartMode);
    const auto result = detectPlaneContact(sourceMesh, chart, delta, projection, 14);
    const double analyticArea = 2.0 * experiments::pi * R * delta;
    return makeRow("sphere_plane", methodPrefix + "_" + chartSuffix(chart), sourceMesh, result, delta, analyticArea);
}

Row runSphereSphereAnalytic(ChartMode chartMode, const std::string& methodPrefix)
{
    const double R = 1.0;
    const double delta = 0.18;
    const double d = 2.0 * R - delta;
    auto sourceMesh = experiments::makeLowerHemispherePatch(R, 40, 96);
    auto targetMesh = experiments::makeSphere(R, 40, 96);
    auto chart = chartFromMesh(sourceMesh, chartMode);
    atlas::ContactAtlas atlasData{{chart}};
    sdf::SphereSDF targetSdf(Eigen::Vector3d::Zero(), R);
    contact::RigidObject source{"source_sphere", &sourceMesh, &atlasData, nullptr, {}};
    contact::RigidObject target{"target_sphere", &targetMesh, nullptr, &targetSdf, {}};
    target.pose.t = Eigen::Vector3d(0, 0, -d);
    contact::ContactDetectorOptions options;
    options.sampleStride = 16;
    const auto result = contact::ContactDetector().detect(source, target, options);
    const double analyticArea = experiments::pi * R * delta;
    return makeRow("sphere_sphere", methodPrefix + "_" + chartSuffix(chart), sourceMesh, result, delta, analyticArea);
}

Row runSphereSphereBvhBaseline()
{
    const double R = 1.0;
    const double delta = 0.18;
    const double d = 2.0 * R - delta;
    auto sourceMesh = experiments::makeLowerHemispherePatch(R, 40, 96);
    auto targetMesh = experiments::makeSphere(R, 40, 96);
    auto chart = chartFromMesh(sourceMesh, ChartMode::OfficialBff);
    atlas::ContactAtlas atlasData{{chart}};
    sdf::MeshSDF targetSdf(targetMesh);
    contact::RigidObject source{"source_sphere", &sourceMesh, &atlasData, nullptr, {}};
    contact::RigidObject target{"target_sphere_bvh", &targetMesh, nullptr, &targetSdf, {}};
    target.pose.t = Eigen::Vector3d(0, 0, -d);
    contact::ContactDetectorOptions options;
    options.sampleStride = 16;
    const auto result = contact::ContactDetector().detect(source, target, options);
    const double analyticArea = experiments::pi * R * delta;
    return makeRow("sphere_sphere", "bvh_exact_triangle_closest_point_" + chartSuffix(chart), sourceMesh, result, delta, analyticArea);
}

Row proceduralPlaneRow(const std::string& scene, const core::Mesh& sourceMesh, double delta, ChartMode chartMode)
{
    auto chart = chartFromMesh(sourceMesh, chartMode);
    const auto result = detectPlaneContact(sourceMesh, chart, delta, true, 12);
    return makeRow(scene, "ours_sdf_projection_" + chartSuffix(chart), sourceMesh, result, delta);
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

bool realAssetPlaneRow(const std::string& scene,
                       core::Mesh fullMesh,
                       double delta,
                       ChartMode chartMode,
                       Row& row,
                       std::string& note)
{
    fullMesh = normalizeMesh(fullMesh, 1.4);
    atlas::Patch patch;
    if (!findBottomDiskPatch(fullMesh, patch)) {
        note = scene + ": no disk-like bottom patch found";
        return false;
    }
    const auto chart = buildChartFromPatch(patch, chartMode);
    const auto result = detectPlaneContact(patch.localMesh, chart, delta, true, 8);
    row = makeRow(scene, "ours_sdf_projection_" + chartSuffix(chart), patch.localMesh, result, delta);
    note = scene + ": full_faces=" + std::to_string(fullMesh.F.rows()) +
           ", patch_faces=" + std::to_string(patch.localMesh.F.rows()) +
           ", chart=" + chart.parameterizationMethod;
    return true;
}

void writeRows(const std::string& path, const std::vector<Row>& rows)
{
    std::ofstream out(path);
    out << "scene,method,mesh_faces,sdf_resolution,contact_area,area_error,max_penetration,penetration_error,normal_error_deg,contour_hausdorff,runtime_ms,sdf_queries,projection_count,projection_failures,memory_mb\n";
    for (const auto& r : rows) {
        out << r.scene << "," << r.method << "," << r.meshFaces << "," << r.sdfResolution << ","
            << r.contactArea << "," << r.areaError << "," << r.maxPenetration << ","
            << r.penetrationError << "," << r.normalErrorDeg << "," << r.contourHausdorff << ","
            << r.runtimeMs << "," << r.sdfQueries << "," << r.projectionCount << ","
            << r.projectionFailures << "," << r.memoryMb << "\n";
    }
}

} // namespace

int main()
{
    const std::string dir = "results/benchmarks";
    io::ensureDirectory(dir);

    std::vector<Row> rows;
    std::vector<std::string> assetNotes;

    rows.push_back(runSpherePlane(44, 120, true, ChartMode::OfficialBff, "ours_sdf_projection"));
    rows.push_back(runSpherePlane(44, 120, false, ChartMode::OfficialBff, "ours_sdf_pullback_only"));
    rows.push_back(runSpherePlane(44, 120, true, ChartMode::PlanarFallback, "ours_sdf_projection"));
    rows.push_back(runSphereSphereAnalytic(ChartMode::OfficialBff, "ours_sdf_projection"));
    rows.push_back(runSphereSphereBvhBaseline());
    rows.push_back(proceduralPlaneRow("ellipsoid_plane", experiments::makeEllipsoidLowerPatch(Eigen::Vector3d(1.2, 0.8, 1.0), 36, 96), 0.2, ChartMode::OfficialBff));
    rows.push_back(proceduralPlaneRow("torus_plane", experiments::makeTorusPatch(0.75, 0.28, 96, 18), 0.32, ChartMode::OfficialBff));

    core::Mesh bunny;
    if (loadBffOfficialObj("bunny.obj", bunny)) {
        Row row;
        std::string note;
        if (realAssetPlaneRow("bunny_plane_real_bff_official_input", bunny, 0.12, ChartMode::OfficialBff, row, note)) {
            rows.push_back(row);
        }
        assetNotes.push_back(note);
    } else {
        assetNotes.push_back("bunny_plane_real_bff_official_input: missing bff_official/input/bunny.obj");
    }

    core::Mesh gear;
    if (loadMechanicalGear(gear)) {
        Row row;
        std::string note;
        if (realAssetPlaneRow("mechanical_part_plane_real_involute_gear", gear, 0.08, ChartMode::OfficialBff, row, note)) {
            rows.push_back(row);
        }
        assetNotes.push_back(note);
    } else {
        assetNotes.push_back("mechanical_part_plane_real_involute_gear: missing tracked CC0 gear STL");
    }

    writeRows(dir + "/accuracy.csv", rows);
    writeRows(dir + "/runtime.csv", rows);

    sdf::SphereSDF sphere(Eigen::Vector3d::Zero(), 1.0);
    core::AABB bounds;
    bounds.expand(Eigen::Vector3d(-1.5, -1.5, -1.5));
    bounds.expand(Eigen::Vector3d(1.5, 1.5, 1.5));
    std::ofstream sdfStudy(dir + "/sdf_resolution_study.csv");
    sdfStudy << "scene,method,sdf_resolution,gap_error,memory_mb\n";
    const Eigen::Vector3d p(1.08, 0.31, 0.17);
    const double exact = sphere.evalLocal(p, false, false).phi;
    auto sphereMesh = experiments::makeSphere(1.0, 42, 96);
    for (int res : {10, 14, 20, 32}) {
        const auto grid = sdf::GridSDF::sampleFrom(sphere, bounds, res, 0.0);
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
    for (int rings : {16, 28, 44}) {
        const int segments = rings * 3;
        const auto row = runSpherePlane(rings, segments, true, ChartMode::OfficialBff, "ours_sdf_projection");
        meshStudy << "sphere_plane,ours_sdf_projection_bff," << rings << "," << segments << ","
                  << row.meshFaces << "," << row.areaError << "," << row.runtimeMs << "\n";
    }

    const Row bffAblation = runSpherePlane(44, 120, true, ChartMode::OfficialBff, "ours_sdf_projection");
    const Row planarAblation = runSpherePlane(44, 120, true, ChartMode::PlanarFallback, "ours_sdf_projection");
    std::ofstream bffVsPlanar(dir + "/bff_vs_planar_ablation.csv");
    bffVsPlanar << "scene,chart_method,area_error,runtime_ms,projection_count,projection_failures\n";
    bffVsPlanar << bffAblation.scene << ",official-bff," << bffAblation.areaError << ","
                << bffAblation.runtimeMs << "," << bffAblation.projectionCount << "," << bffAblation.projectionFailures << "\n";
    bffVsPlanar << planarAblation.scene << ",fallback-planar," << planarAblation.areaError << ","
                << planarAblation.runtimeMs << "," << planarAblation.projectionCount << "," << planarAblation.projectionFailures << "\n";

    std::ofstream ablation(dir + "/ablation.csv");
    ablation << "question,variant,answer_metric,value,notes\n";
    ablation << "projection reduces GridSDF error,no_projection,gap_error,see_sdf_resolution_study,pure grid row\n";
    ablation << "projection reduces GridSDF error,with_projection,gap_error,see_sdf_resolution_study,projection row\n";
    ablation << "SDF reduces expensive projections,no_sdf_projection_count,projection_count," << rows[0].meshFaces << ",brute-force would project every triangle\n";
    ablation << "SDF reduces expensive projections,with_sdf_projection_count,projection_count," << rows[0].projectionCount << ",narrow-band sample projection\n";
    ablation << "BFF chart vs fallback planar,official_bff,area_error," << bffAblation.areaError << ",see_bff_vs_planar_ablation\n";
    ablation << "BFF chart vs fallback planar,fallback_planar,area_error," << planarAblation.areaError << ",see_bff_vs_planar_ablation\n";
    ablation << "BVH baseline,bvh_exact_triangle_closest_point,area_error," << rows[4].areaError << ",sphere_sphere closed target mesh baseline\n";
    ablation << "mesh resolution affects area,coarse_to_medium,area_error,see_mesh_resolution_study,computed rows\n";

    std::ofstream report(dir + "/validation_report.md");
    report << "# Benchmark Validation Report\n\n"
           << "- Generated `accuracy.csv`, `runtime.csv`, `sdf_resolution_study.csv`, `mesh_resolution_study.csv`, `ablation.csv`, and `bff_vs_planar_ablation.csv`.\n"
           << "- Sphere-plane BFF area error: " << rows[0].areaError << "\n"
           << "- Sphere-sphere analytic area error: " << rows[3].areaError << "\n"
           << "- Sphere-sphere BVH exact-triangle baseline area error: " << rows[4].areaError << "\n"
           << "- Disk-like procedural patches request official BFF charts; non-disk topology falls back and is marked in the method name.\n";
    for (const auto& note : assetNotes) report << "- " << note << "\n";
    report << "- Mechanical gear asset source: koppi/involute-gear-collection, CC0-1.0, `angle-20/teeth-16.stl`.\n";

    std::cout << "benchmarks written to " << dir << "\n";
    return 0;
}
