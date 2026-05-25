#include "test_main.h"

#include "atlas/ContactAtlas.h"
#include "atlas/PatchBuilder.h"
#include "contact/ContactIntegrator.h"
#include "procedural_test_meshes.h"
#include "sdf/AnalyticSDF.h"

#include <cmath>

using bff_sdf::atlas::PatchBuilder;
using bff_sdf::atlas::makeContactChart;
using bff_sdf::contact::ContactIntegrator;
using bff_sdf::core::RigidPose;
using bff_sdf::sdf::PlaneSDF;

int main()
{
    return test::run("test_contact_area", [] {
        constexpr double pi = 3.141592653589793238462643383279502884;
        const double R = 1.0;
        const double delta = 0.2;
        const int rings = 64;
        const int segments = 160;

        auto mesh = test_meshes::makeLowerHemispherePatch(R, rings, segments);
        auto patch = PatchBuilder::buildWholeMeshPatch(mesh);
        REQUIRE_MSG(patch.diskLike, patch.diagnostic);

        Eigen::MatrixXd UV(mesh.V.rows(), 2);
        for (int i = 0; i < mesh.V.rows(); ++i) {
            UV(i, 0) = mesh.V(i, 0);
            UV(i, 1) = mesh.V(i, 1);
        }
        const auto chart = makeContactChart(0, patch, UV, mesh.F, "radial-test", false);

        RigidPose sourcePose;
        sourcePose.t = Eigen::Vector3d(0, 0, R - delta);
        RigidPose targetPose;
        PlaneSDF plane(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());

        const auto result = ContactIntegrator::integrateChart(chart, sourcePose, plane, targetPose);
        const double analyticArea = 2.0 * pi * R * delta;
        const double areaError = std::abs(result.area - analyticArea) / analyticArea;
        REQUIRE_MSG(areaError < 0.03, "sphere-plane cap area error too high");

        const double analyticRadius = std::sqrt(2.0 * R * delta - delta * delta);
        double radiusSum = 0.0;
        int radiusCount = 0;
        for (const auto& contour : result.contours) {
            for (const auto& p : contour.xWorldPoints) {
                radiusSum += p.head<2>().norm();
                ++radiusCount;
                REQUIRE(std::abs(p.z()) < 5e-3);
            }
        }
        REQUIRE(radiusCount > 0);
        const double meanRadius = radiusSum / static_cast<double>(radiusCount);
        const double radiusError = std::abs(meanRadius - analyticRadius) / analyticRadius;
        REQUIRE_MSG(radiusError < 0.02, "sphere-plane contact radius error too high");
        REQUIRE(result.maxPenetration > 0.0);
        REQUIRE(result.meanPenetration > 0.0);
    });
}

