#include "test_main.h"

#include "atlas/ContactAtlas.h"
#include "atlas/PatchBuilder.h"
#include "contact/ContactGapField.h"
#include "sdf/AnalyticSDF.h"

using bff_sdf::atlas::PatchBuilder;
using bff_sdf::atlas::makeContactChart;
using bff_sdf::contact::ContactGapField;
using bff_sdf::core::Mesh;
using bff_sdf::core::RigidPose;
using bff_sdf::sdf::PlaneSDF;

int main()
{
    return test::run("test_gap_field", [] {
        Mesh mesh;
        mesh.V.resize(3, 3);
        mesh.V << 0, 0, 0.1,
                  1, 0, 1.1,
                  0, 1, 2.1;
        mesh.F.resize(1, 3);
        mesh.F << 0, 1, 2;
        mesh.computeNormals();

        auto patch = PatchBuilder::buildWholeMeshPatch(mesh);
        Eigen::MatrixXd UV(3, 2);
        UV << 0, 0,
              1, 0,
              0, 1;
        Eigen::MatrixXi Fuv = mesh.F;
        auto chart = makeContactChart(7, patch, UV, Fuv, "manual", false);

        PlaneSDF plane(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
        RigidPose sourcePose;
        RigidPose targetPose;
        const Eigen::Vector3d bary(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
        const auto sample = ContactGapField::sample(chart, 0, bary, sourcePose, plane, targetPose);

        REQUIRE(sample.sourceChartId == 7);
        REQUIRE_NEAR(sample.g, 1.1, 1e-12);
        REQUIRE_VEC_NEAR(sample.gradU, Eigen::Vector2d(1.0, 2.0), 1e-12);

        const double h = 1e-6;
        auto evalAtUv = [&](const Eigen::Vector2d& u) {
            int tri = -1;
            Eigen::Vector3d b;
            REQUIRE(chart.uvAccel.locate(u, tri, b));
            return ContactGapField::sample(chart, tri, b, sourcePose, plane, targetPose).g;
        };
        const Eigen::Vector2d u = sample.xWorld.head<2>();
        Eigen::Vector2d fd;
        fd.x() = (evalAtUv(u + Eigen::Vector2d(h, 0)) - evalAtUv(u - Eigen::Vector2d(h, 0))) / (2.0 * h);
        fd.y() = (evalAtUv(u + Eigen::Vector2d(0, h)) - evalAtUv(u - Eigen::Vector2d(0, h))) / (2.0 * h);
        REQUIRE_VEC_NEAR(sample.gradU, fd, 1e-6);
    });
}

