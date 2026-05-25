#include "test_main.h"

#include "core/RigidPose.h"
#include "procedural_test_meshes.h"
#include "sdf/AnalyticSDF.h"
#include "sdf/GridSDF.h"
#include "sdf/MeshSDF.h"

#include <vector>

using bff_sdf::core::RigidPose;
using bff_sdf::sdf::BoxSDF;
using bff_sdf::sdf::GridSDF;
using bff_sdf::sdf::MeshSDF;
using bff_sdf::sdf::PlaneSDF;
using bff_sdf::sdf::SphereSDF;

namespace {

Eigen::Vector3d finiteDifferenceGrad(const bff_sdf::sdf::ISDF& sdf, const Eigen::Vector3d& p)
{
    const double h = 1e-6;
    Eigen::Vector3d g;
    for (int a = 0; a < 3; ++a) {
        Eigen::Vector3d e = Eigen::Vector3d::Zero();
        e[a] = h;
        g[a] = (sdf.evalLocal(p + e, false, false).phi - sdf.evalLocal(p - e, false, false).phi) / (2.0 * h);
    }
    return g;
}

} // namespace

int main()
{
    return test::run("test_sdf", [] {
        SphereSDF sphere(Eigen::Vector3d::Zero(), 2.0);
        Eigen::Vector3d p(3.0, 4.0, 0.0);
        auto qs = sphere.evalLocal(p, true, true);
        REQUIRE_NEAR(qs.phi, 3.0, 1e-10);
        REQUIRE_VEC_NEAR(qs.grad, finiteDifferenceGrad(sphere, p), 1e-6);
        REQUIRE_NEAR(sphere.evalLocal(qs.closestPoint, false, false).phi, 0.0, 1e-10);

        PlaneSDF plane(Eigen::Vector3d(0.0, 0.0, 1.0), Eigen::Vector3d::UnitZ());
        auto qp = plane.evalLocal(Eigen::Vector3d(0.2, -0.3, 3.5), true, true);
        REQUIRE_NEAR(qp.phi, 2.5, 1e-10);
        REQUIRE_VEC_NEAR(qp.grad, finiteDifferenceGrad(plane, Eigen::Vector3d(0.2, -0.3, 3.5)), 1e-6);
        REQUIRE_NEAR(plane.evalLocal(qp.closestPoint, false, false).phi, 0.0, 1e-10);

        BoxSDF box(Eigen::Vector3d::Zero(), Eigen::Vector3d(1.0, 2.0, 3.0));
        auto qb = box.evalLocal(Eigen::Vector3d(2.5, 0.0, 0.0), true, true);
        REQUIRE_NEAR(qb.phi, 1.5, 1e-8);
        REQUIRE_VEC_NEAR(qb.grad, finiteDifferenceGrad(box, Eigen::Vector3d(2.5, 0.0, 0.0)), 1e-6);
        REQUIRE_NEAR(box.evalLocal(qb.closestPoint, false, false).phi, 0.0, 1e-8);

        RigidPose targetPose;
        targetPose.R = Eigen::AngleAxisd(0.25, Eigen::Vector3d::UnitY()).toRotationMatrix();
        targetPose.t = Eigen::Vector3d(0.4, 0.1, -0.2);
        const Eigen::Vector3d yLocal(0.0, 0.0, 2.75);
        const Eigen::Vector3d xWorld = targetPose.localToWorld(yLocal);
        const Eigen::Vector3d recovered = targetPose.worldToLocal(xWorld);
        REQUIRE_VEC_NEAR(recovered, yLocal, 1e-12);
        const auto qLocal = sphere.evalLocal(recovered, true, false);
        const Eigen::Vector3d normalWorld = targetPose.rotateLocalVectorToWorld(qLocal.grad);
        REQUIRE_NEAR(normalWorld.norm(), 1.0, 1e-12);

        const auto cubeMesh = test_meshes::makeCube();
        MeshSDF cubeSdf(cubeMesh);
        REQUIRE(cubeSdf.evalLocal(Eigen::Vector3d(0, 0, 0), true, true).phi < 0.0);
        REQUIRE(cubeSdf.evalLocal(Eigen::Vector3d(2, 0, 0), true, true).phi > 0.0);

        const auto sphereMesh = test_meshes::makeSphere(1.0, 24, 48);
        MeshSDF meshSphere(sphereMesh);
        const Eigen::Vector3d outside(1.4, 0.2, 0.1);
        const double analyticPhi = SphereSDF(Eigen::Vector3d::Zero(), 1.0).evalLocal(outside, false, false).phi;
        REQUIRE(std::abs(meshSphere.evalLocal(outside, true, true).phi - analyticPhi) < 0.02);

        bff_sdf::core::AABB bounds;
        bounds.expand(Eigen::Vector3d(-1.5, -1.5, -1.5));
        bounds.expand(Eigen::Vector3d(1.5, 1.5, 1.5));
        SphereSDF unitSphere(Eigen::Vector3d::Zero(), 1.0);
        const auto gridLow = GridSDF::sampleFrom(unitSphere, bounds, 12, 0.0);
        const auto gridHigh = GridSDF::sampleFrom(unitSphere, bounds, 36, 0.0);
        const std::vector<Eigen::Vector3d> samples = {
            {0.8, 0.3, 0.4},
            {1.2, 0.2, 0.1},
            {-0.4, 0.9, 0.25},
            {0.1, -1.15, 0.2}
        };
        double lowErr = 0.0;
        double highErr = 0.0;
        for (const auto& s : samples) {
            const double exact = unitSphere.evalLocal(s, false, false).phi;
            lowErr += std::abs(gridLow.evalLocal(s, false, false).phi - exact);
            highErr += std::abs(gridHigh.evalLocal(s, false, false).phi - exact);
        }
        REQUIRE(highErr < lowErr);
        const auto qGrid = gridHigh.evalLocal(Eigen::Vector3d(1.05, 0.15, 0.1), true, true);
        REQUIRE_NEAR(qGrid.grad.norm(), 1.0, 1e-8);
        REQUIRE(qGrid.hasClosest);
    });
}
