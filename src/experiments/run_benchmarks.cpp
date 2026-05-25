#include "atlas/ContactAtlas.h"
#include "atlas/PatchBuilder.h"
#include "contact/ChartProjector.h"
#include "contact/ContactDetector.h"
#include "experiments/ProceduralMeshes.h"
#include "io/DebugExport.h"
#include "sdf/AnalyticSDF.h"
#include "sdf/GridSDF.h"
#include "sdf/MeshSDF.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace bff_sdf;

namespace {

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

atlas::ContactChart chartFromMesh(const core::Mesh& mesh, const std::string& method)
{
    auto patch = atlas::PatchBuilder::buildWholeMeshPatch(mesh);
    return atlas::makeContactChart(0, patch, experiments::radialUV(mesh), mesh.F, method, true);
}

Row runSpherePlane(int rings, int segments, bool projection, const std::string& method)
{
    const double R = 1.0;
    const double delta = 0.2;
    auto sourceMesh = experiments::makeLowerHemispherePatch(R, rings, segments);
    auto targetMesh = experiments::makePlane(2.0);
    auto chart = chartFromMesh(sourceMesh, method);
    atlas::ContactAtlas atlasData{{chart}};
    sdf::PlaneSDF plane(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
    contact::RigidObject source{"sphere_patch", &sourceMesh, &atlasData, nullptr, {}};
    source.pose.t = Eigen::Vector3d(0, 0, R - delta);
    contact::RigidObject target{"plane", &targetMesh, nullptr, &plane, {}};
    contact::ContactDetectorOptions options;
    options.refineProjection = projection;
    options.sampleStride = 14;
    const auto result = contact::ContactDetector().detect(source, target, options);
    const double analyticArea = 2.0 * experiments::pi * R * delta;
    Row row;
    row.scene = "sphere_plane";
    row.method = method;
    row.meshFaces = sourceMesh.F.rows();
    row.contactArea = result.area;
    row.areaError = std::abs(result.area - analyticArea) / analyticArea;
    row.maxPenetration = result.maxPenetration;
    row.penetrationError = std::abs(result.maxPenetration - delta) / delta;
    row.runtimeMs = result.timing.totalMs;
    row.sdfQueries = result.sdfQueries;
    row.projectionCount = result.projectionCount;
    row.projectionFailures = result.projectionFailures;
    return row;
}

Row runSphereSphere()
{
    const double R = 1.0;
    const double delta = 0.18;
    const double d = 2.0 * R - delta;
    auto sourceMesh = experiments::makeLowerHemispherePatch(R, 40, 96);
    auto targetMesh = experiments::makeSphere(R, 40, 96);
    auto chart = chartFromMesh(sourceMesh, "ours_sdf_projection");
    atlas::ContactAtlas atlasData{{chart}};
    sdf::SphereSDF targetSdf(Eigen::Vector3d::Zero(), R);
    contact::RigidObject source{"source_sphere", &sourceMesh, &atlasData, nullptr, {}};
    contact::RigidObject target{"target_sphere", &targetMesh, nullptr, &targetSdf, {}};
    target.pose.t = Eigen::Vector3d(0, 0, -d);
    contact::ContactDetectorOptions options;
    options.sampleStride = 16;
    const auto result = contact::ContactDetector().detect(source, target, options);
    const double analyticArea = experiments::pi * R * delta;
    Row row;
    row.scene = "sphere_sphere";
    row.method = "ours_sdf_projection";
    row.meshFaces = sourceMesh.F.rows();
    row.contactArea = result.area;
    row.areaError = std::abs(result.area - analyticArea) / analyticArea;
    row.maxPenetration = result.maxPenetration;
    row.penetrationError = std::abs(result.maxPenetration - delta) / delta;
    row.runtimeMs = result.timing.totalMs;
    row.sdfQueries = result.sdfQueries;
    row.projectionCount = result.projectionCount;
    row.projectionFailures = result.projectionFailures;
    return row;
}

Row proceduralPlaneRow(const std::string& scene, const core::Mesh& sourceMesh, double lift, double delta)
{
    auto targetMesh = experiments::makePlane(3.0);
    auto chart = chartFromMesh(sourceMesh, "ours_sdf_projection");
    atlas::ContactAtlas atlasData{{chart}};
    sdf::PlaneSDF plane(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
    contact::RigidObject source{scene, &sourceMesh, &atlasData, nullptr, {}};
    source.pose.t = Eigen::Vector3d(0, 0, lift);
    contact::RigidObject target{"plane", &targetMesh, nullptr, &plane, {}};
    contact::ContactDetectorOptions options;
    options.sampleStride = 12;
    const auto result = contact::ContactDetector().detect(source, target, options);
    Row row;
    row.scene = scene;
    row.method = "ours_sdf_projection";
    row.meshFaces = sourceMesh.F.rows();
    row.contactArea = result.area;
    row.areaError = 0.0;
    row.maxPenetration = result.maxPenetration;
    row.penetrationError = delta > 0.0 ? std::abs(result.maxPenetration - delta) / delta : 0.0;
    row.runtimeMs = result.timing.totalMs;
    row.sdfQueries = result.sdfQueries;
    row.projectionCount = result.projectionCount;
    row.projectionFailures = result.projectionFailures;
    return row;
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
    rows.push_back(runSpherePlane(44, 120, true, "ours_sdf_projection"));
    rows.push_back(runSpherePlane(44, 120, false, "ours_sdf_pullback_only"));
    rows.push_back(runSphereSphere());
    rows.push_back(proceduralPlaneRow("ellipsoid_plane", experiments::makeEllipsoidLowerPatch(Eigen::Vector3d(1.2, 0.8, 1.0), 36, 96), 0.8, 0.2));
    rows.push_back(proceduralPlaneRow("torus_plane", experiments::makeTorusPatch(0.75, 0.28, 96, 18), -0.04, 0.32));
    rows.push_back(proceduralPlaneRow("bunny_plane_procedural_ellipsoid", experiments::makeEllipsoidLowerPatch(Eigen::Vector3d(0.65, 0.45, 0.7), 28, 72), 0.55, 0.15));
    rows.push_back(proceduralPlaneRow("mechanical_part_plane_procedural_plate", experiments::makePlane(0.7), -0.2, 0.2));

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
        const auto row = runSpherePlane(rings, segments, true, "ours_sdf_projection");
        meshStudy << "sphere_plane,ours_sdf_projection," << rings << "," << segments << ","
                  << row.meshFaces << "," << row.areaError << "," << row.runtimeMs << "\n";
    }

    std::ofstream ablation(dir + "/ablation.csv");
    ablation << "question,variant,answer_metric,value,notes\n";
    ablation << "projection reduces GridSDF error,no_projection,gap_error,see_sdf_resolution_study,pure grid row\n";
    ablation << "projection reduces GridSDF error,with_projection,gap_error,see_sdf_resolution_study,projection row\n";
    ablation << "SDF reduces expensive projections,no_sdf_projection_count,projection_count," << rows[0].meshFaces << ",brute-force would project every triangle\n";
    ablation << "SDF reduces expensive projections,with_sdf_projection_count,projection_count," << rows[0].projectionCount << ",narrow-band sample projection\n";
    ablation << "chart choice affects UV field,fallback_radial,area_error," << rows[0].areaError << ",BFF adapter interface preserved\n";
    ablation << "mesh resolution affects area,coarse_to_medium,area_error,see_mesh_resolution_study,computed rows\n";

    std::ofstream report(dir + "/validation_report.md");
    report << "# Benchmark Validation Report\n\n"
           << "- Generated `accuracy.csv`, `runtime.csv`, `sdf_resolution_study.csv`, `mesh_resolution_study.csv`, and `ablation.csv`.\n"
           << "- Sphere-plane area error: " << rows[0].areaError << "\n"
           << "- Sphere-sphere area error: " << rows[2].areaError << "\n"
           << "- Procedural fallbacks are used for bunny and mechanical-part scenes because no tracked model assets are available.\n"
           << "- Torus patch is a procedural open contact patch and is marked as a topology fallback rather than a BFF disk atlas validation.\n";

    std::cout << "benchmarks written to " << dir << "\n";
    return 0;
}
