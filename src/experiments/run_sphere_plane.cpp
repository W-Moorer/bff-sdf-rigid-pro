#include "atlas/BFFAdapter.h"
#include "atlas/PatchBuilder.h"
#include "contact/ContactDetector.h"
#include "experiments/ProceduralMeshes.h"
#include "io/DebugExport.h"
#include "sdf/AnalyticSDF.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

using namespace bff_sdf;

namespace {

struct Metrics {
    double analyticArea = 0.0;
    double areaError = 0.0;
    double analyticRadius = 0.0;
    double radiusError = 0.0;
    double penetrationError = 0.0;
};

atlas::ContactChart buildSourceChart(const core::Mesh& mesh, bool tryOfficialBff)
{
    auto patch = atlas::PatchBuilder::buildWholeMeshPatch(mesh);
    atlas::BFFAdapter adapter;
    atlas::BFFAdapterOptions options;
    options.preferOfficialBff = tryOfficialBff;
    auto bff = adapter.parameterize(patch, options);
    if (!bff.success) {
        bff.UV = experiments::radialUV(mesh);
        bff.Fuv = mesh.F;
        bff.usedFallback = true;
    }
    return atlas::makeContactChart(0, patch, bff.UV, bff.Fuv,
                                   bff.usedOfficialBff ? "official-bff" : "fallback-radial",
                                   bff.usedFallback);
}

Metrics computeMetrics(const contact::ContactResult& result, double R, double delta)
{
    Metrics m;
    m.analyticArea = 2.0 * experiments::pi * R * delta;
    m.areaError = std::abs(result.area - m.analyticArea) / m.analyticArea;
    m.analyticRadius = std::sqrt(2.0 * R * delta - delta * delta);
    double radiusSum = 0.0;
    int radiusCount = 0;
    for (const auto& contour : result.contours) {
        for (const auto& p : contour.xWorldPoints) {
            radiusSum += p.head<2>().norm();
            ++radiusCount;
        }
    }
    const double meanRadius = radiusCount > 0 ? radiusSum / radiusCount : 0.0;
    m.radiusError = radiusCount > 0 ? std::abs(meanRadius - m.analyticRadius) / m.analyticRadius : 1.0;
    m.penetrationError = std::abs(result.maxPenetration - delta) / delta;
    return m;
}

void writeMetricCsv(const std::string& path,
                    const contact::ContactResult& result,
                    const Metrics& metrics,
                    int meshFaces,
                    int sdfResolution)
{
    std::ofstream out(path);
    out << "scene,method,mesh_faces,sdf_resolution,contact_area,area_error,max_penetration,penetration_error,normal_error_deg,contour_hausdorff,runtime_ms,sdf_queries,projection_count,projection_failures,memory_mb\n";
    out << "sphere_plane,ours_sdf_projection," << meshFaces << "," << sdfResolution << ","
        << result.area << "," << metrics.areaError << ","
        << result.maxPenetration << "," << metrics.penetrationError << ","
        << 0.0 << "," << metrics.radiusError << ","
        << result.timing.totalMs << "," << result.sdfQueries << ","
        << result.projectionCount << "," << result.projectionFailures << ",0\n";
}

void writeTiming(const std::string& csvPath, const std::string& jsonPath, const contact::ContactResult& result)
{
    std::ofstream csv(csvPath);
    csv << "stage,runtime_ms\n";
    csv << "integration," << result.timing.integrationMs << "\n";
    csv << "projection," << result.timing.projectionMs << "\n";
    csv << "total," << result.timing.totalMs << "\n";

    std::ofstream json(jsonPath);
    json << "{\n"
         << "  \"integration_ms\": " << result.timing.integrationMs << ",\n"
         << "  \"projection_ms\": " << result.timing.projectionMs << ",\n"
         << "  \"total_ms\": " << result.timing.totalMs << "\n"
         << "}\n";
}

} // namespace

int main()
{
    const std::string dir = "results/sphere_plane";
    io::ensureDirectory(dir);

    const double R = 1.0;
    const double delta = 0.2;
    auto sourceMesh = experiments::makeLowerHemispherePatch(R, 56, 144);
    auto targetMesh = experiments::makePlane(2.0);
    auto chart = buildSourceChart(sourceMesh, true);
    atlas::ContactAtlas atlasData;
    atlasData.charts.push_back(chart);

    sdf::PlaneSDF plane(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
    contact::RigidObject source;
    source.name = "source_sphere_patch";
    source.mesh = &sourceMesh;
    source.atlas = &atlasData;
    source.pose.t = Eigen::Vector3d(0, 0, R - delta);
    contact::RigidObject target;
    target.name = "target_plane";
    target.mesh = &targetMesh;
    target.sdf = &plane;

    contact::ContactDetector detector;
    contact::ContactDetectorOptions options;
    options.sampleStride = 10;
    const auto result = detector.detect(source, target, options);
    const auto metrics = computeMetrics(result, R, delta);

    std::ofstream config(dir + "/config.json");
    config << "{\n  \"scene\": \"sphere_plane\",\n  \"radius\": " << R
           << ",\n  \"penetration_depth\": " << delta
           << ",\n  \"chart_method\": \"" << chart.parameterizationMethod << "\"\n}\n";

    writeMetricCsv(dir + "/metrics.csv", result, metrics, sourceMesh.F.rows(), 0);
    writeTiming(dir + "/timing.csv", dir + "/timing.json", result);
    io::writeGapFieldCsv(dir + "/gap_field.csv", chart, source.pose, plane, target.pose);
    io::writeGapFieldCsv(dir + "/uv_gap_field.csv", chart, source.pose, plane, target.pose);
    io::writeContactSamplesObj(dir + "/contact_samples.obj", result);
    io::writeContactSamplesPly(dir + "/contact_samples.ply", result);
    io::writeContactContourObj(dir + "/contact_contour.obj", result);
    io::writeContactContourObj(dir + "/contact_contour_3d.obj", result);
    io::writeContactContourUvCsv(dir + "/contact_contour_uv.csv", result);
    io::writeAtlasDebugObj(dir + "/atlas_debug.obj", chart, source.pose);

    std::ofstream report(dir + "/report.md");
    report << "# Sphere-Plane Report\n\n"
           << "- Chart method: " << chart.parameterizationMethod << "\n"
           << "- Contact area: " << result.area << "\n"
           << "- Analytic area: " << metrics.analyticArea << "\n"
           << "- Area error: " << metrics.areaError << "\n"
           << "- Contact radius error: " << metrics.radiusError << "\n"
           << "- Projection failures: " << result.projectionFailures << "\n";

    std::cout << "sphere_plane area_error=" << metrics.areaError
              << " radius_error=" << metrics.radiusError << "\n";
    return (metrics.areaError < 0.03 && metrics.radiusError < 0.02) ? 0 : 2;
}

