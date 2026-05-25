#include "test_main.h"

#include "contact/ChartProjector.h"
#include "procedural_test_meshes.h"
#include "sdf/AnalyticSDF.h"
#include "sdf/GridSDF.h"
#include "sdf/MeshSDF.h"

#include <cmath>

using bff_sdf::contact::ChartProjector;
using bff_sdf::contact::ProjectionOptions;
using bff_sdf::core::Mesh;
using bff_sdf::core::RigidPose;
using bff_sdf::sdf::GridSDF;
using bff_sdf::sdf::MeshSDF;
using bff_sdf::sdf::PlaneSDF;
using bff_sdf::sdf::SphereSDF;

namespace {

Mesh makePlaneMesh()
{
    Mesh m;
    m.V.resize(4, 3);
    m.V << -2, -2, 0,
            2, -2, 0,
            2,  2, 0,
           -2,  2, 0;
    m.F.resize(2, 3);
    m.F << 0, 1, 2,
           0, 2, 3;
    m.computeNormals();
    return m;
}

} // namespace

int main()
{
    return test::run("test_projection", [] {
        RigidPose pose;
        auto planeMesh = makePlaneMesh();
        PlaneSDF plane(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
        const Eigen::Vector3d p(0.25, -0.4, 0.75);
        const auto hintPlane = plane.evalLocal(p, true, true);
        ChartProjector projector;
        const auto projectedPlane = projector.project(p, planeMesh, pose, &hintPlane, ProjectionOptions{});
        REQUIRE(projectedPlane.converged);
        REQUIRE_VEC_NEAR(projectedPlane.closestPointWorld, Eigen::Vector3d(0.25, -0.4, 0.0), 1e-8);
        REQUIRE_NEAR(projectedPlane.signedGap, 0.75, 1e-8);

        const auto sphereMesh = test_meshes::makeSphere(1.0, 40, 80);
        MeshSDF meshSphere(sphereMesh);
        const Eigen::Vector3d query(1.15, 0.1, 0.05);
        const auto meshHint = meshSphere.evalLocal(query, true, true);
        const auto projectedSphere = projector.project(query, sphereMesh, pose, &meshHint, ProjectionOptions{});
        REQUIRE(projectedSphere.converged);
        REQUIRE_VEC_NEAR(projectedSphere.closestPointWorld, meshHint.closestPoint, 1e-12);

        bff_sdf::core::AABB bounds;
        bounds.expand(Eigen::Vector3d(-1.5, -1.5, -1.5));
        bounds.expand(Eigen::Vector3d(1.5, 1.5, 1.5));
        SphereSDF analyticSphere(Eigen::Vector3d::Zero(), 1.0);
        const auto grid = GridSDF::sampleFrom(analyticSphere, bounds, 14, 0.0);
        const Eigen::Vector3d gridQuery(1.08, 0.31, 0.17);
        const double exactPhi = analyticSphere.evalLocal(gridQuery, false, false).phi;
        const auto gridHint = grid.evalLocal(gridQuery, true, true);
        const auto refined = projector.project(gridQuery, sphereMesh, pose, &gridHint, ProjectionOptions{});
        REQUIRE(refined.converged);
        const double pureGridError = std::abs(gridHint.phi - exactPhi);
        const double refinedError = std::abs(refined.signedGap - exactPhi);
        REQUIRE_MSG(refinedError < pureGridError, "projection refinement should reduce GridSDF gap error");

        const auto counters = projector.counters();
        REQUIRE(counters.projection_attempts >= 3);
        REQUIRE(counters.projection_converged == counters.projection_attempts);
    });
}

