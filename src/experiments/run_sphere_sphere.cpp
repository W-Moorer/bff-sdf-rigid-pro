#include "atlas/BFFAdapter.h"
#include "atlas/PatchBuilder.h"
#include "contact/ContactDetector.h"
#include "experiments/ProceduralMeshes.h"
#include "io/DebugExport.h"
#include "sdf/AnalyticSDF.h"

#include <cmath>
#include <fstream>
#include <iostream>

using namespace bff_sdf;

namespace {

atlas::ContactChart buildChart(const core::Mesh& mesh)
{
    auto patch = atlas::PatchBuilder::buildWholeMeshPatch(mesh);
    atlas::BFFAdapter adapter;
    atlas::BFFAdapterOptions options;
    options.preferOfficialBff = false;
    const auto result = adapter.parameterize(patch, options);
    return atlas::makeContactChart(0, patch, result.UV, result.Fuv,
                                   result.usedOfficialBff ? "official-bff" : "fallback-radial",
                                   result.usedFallback);
}

double meanContourRadius(const contact::ContactResult& result)
{
    double sum = 0.0;
    int count = 0;
    for (const auto& contour : result.contours) {
        for (const auto& p : contour.xWorldPoints) {
            sum += p.head<2>().norm();
            ++count;
        }
    }
    return count > 0 ? sum / count : 0.0;
}

} // namespace

int main()
{
    const std::string dir = "results/sphere_sphere";
    io::ensureDirectory(dir);

    const double R = 1.0;
    const double delta = 0.18;
    const double d = 2.0 * R - delta;

    auto sourceMesh = experiments::makeLowerHemispherePatch(R, 56, 144);
    auto targetMesh = experiments::makeSphere(R, 56, 144);
    auto chart = buildChart(sourceMesh);
    atlas::ContactAtlas atlasData;
    atlasData.charts.push_back(chart);

    sdf::SphereSDF targetSdf(Eigen::Vector3d::Zero(), R);
    contact::RigidObject source;
    source.name = "source_sphere";
    source.mesh = &sourceMesh;
    source.atlas = &atlasData;
    contact::RigidObject target;
    target.name = "target_sphere";
    target.mesh = &targetMesh;
    target.sdf = &targetSdf;
    target.pose.t = Eigen::Vector3d(0, 0, -d);

    contact::ContactDetector detector;
    contact::ContactDetectorOptions options;
    options.sampleStride = 12;
    const auto result = detector.detect(source, target, options);

    const double analyticArea = experiments::pi * R * delta;
    const double areaError = std::abs(result.area - analyticArea) / analyticArea;
    const double analyticRadius = std::sqrt(R * R - (d * 0.5) * (d * 0.5));
    const double radius = meanContourRadius(result);
    const double radiusError = std::abs(radius - analyticRadius) / analyticRadius;

    std::ofstream config(dir + "/config.json");
    config << "{\n  \"scene\": \"sphere_sphere\",\n  \"radius\": " << R
           << ",\n  \"penetration_depth\": " << delta << "\n}\n";

    std::ofstream metrics(dir + "/metrics.csv");
    metrics << "scene,method,mesh_faces,sdf_resolution,contact_area,area_error,max_penetration,penetration_error,normal_error_deg,contour_hausdorff,runtime_ms,sdf_queries,projection_count,projection_failures,memory_mb\n";
    metrics << "sphere_sphere,ours_sdf_projection," << sourceMesh.F.rows() << ",0,"
            << result.area << "," << areaError << ","
            << result.maxPenetration << "," << std::abs(result.maxPenetration - delta) / delta << ",0,"
            << radiusError << "," << result.timing.totalMs << ","
            << result.sdfQueries << "," << result.projectionCount << ","
            << result.projectionFailures << ",0\n";

    io::writeGapFieldCsv(dir + "/gap_field.csv", chart, source.pose, targetSdf, target.pose);
    io::writeContactSamplesObj(dir + "/contact_samples.obj", result);
    io::writeContactSamplesPly(dir + "/contact_samples.ply", result);
    io::writeContactContourObj(dir + "/contact_contour_3d.obj", result);
    io::writeContactContourUvCsv(dir + "/contact_contour_uv.csv", result);
    io::writeAtlasDebugObj(dir + "/atlas_debug.obj", chart, source.pose);

    std::ofstream report(dir + "/report.md");
    report << "# Sphere-Sphere Report\n\n"
           << "- Analytic contact area: " << analyticArea << "\n"
           << "- Integrated contact area: " << result.area << "\n"
           << "- Area error: " << areaError << "\n"
           << "- Radius error: " << radiusError << "\n";

    std::cout << "sphere_sphere area_error=" << areaError
              << " radius_error=" << radiusError << "\n";
    return (areaError < 0.04 && radiusError < 0.03) ? 0 : 2;
}

