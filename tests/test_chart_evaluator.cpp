#include "test_main.h"

#include "atlas/ChartEvaluator.h"
#include "atlas/ContactAtlas.h"
#include "atlas/PatchBuilder.h"
#include "core/RigidPose.h"

using bff_sdf::atlas::ChartEvaluator;
using bff_sdf::atlas::PatchBuilder;
using bff_sdf::atlas::makeContactChart;
using bff_sdf::core::Mesh;
using bff_sdf::core::RigidPose;

namespace {

Mesh makeTriangle()
{
    Mesh m;
    m.V.resize(3, 3);
    m.V << 0, 0, 0,
           1, 0, 0,
           0, 1, 0;
    m.F.resize(1, 3);
    m.F << 0, 1, 2;
    m.computeNormals();
    return m;
}

} // namespace

int main()
{
    return test::run("test_chart_evaluator", [] {
        const Mesh mesh = makeTriangle();
        auto patch = PatchBuilder::buildWholeMeshPatch(mesh);
        REQUIRE(patch.diskLike);

        Eigen::MatrixXd UV(3, 2);
        UV << 0, 0,
              1, 0,
              0, 1;
        Eigen::MatrixXi Fuv = mesh.F;
        const auto chart = makeContactChart(0, patch, UV, Fuv, "manual", false);

        RigidPose pose;
        const auto v0 = ChartEvaluator::evaluate(chart, 0, Eigen::Vector3d(1, 0, 0), pose);
        REQUIRE(v0.valid);
        REQUIRE_VEC_NEAR(v0.X0, Eigen::Vector3d(0, 0, 0), 1e-12);

        const Eigen::Vector3d bary(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
        const auto center = ChartEvaluator::evaluate(chart, 0, bary, pose);
        REQUIRE(center.valid);
        REQUIRE_VEC_NEAR(center.X0, Eigen::Vector3d(1.0 / 3.0, 1.0 / 3.0, 0.0), 1e-12);
        REQUIRE_VEC_NEAR(center.u, Eigen::Vector2d(1.0 / 3.0, 1.0 / 3.0), 1e-12);
        REQUIRE_VEC_NEAR(center.normalWorld, Eigen::Vector3d::UnitZ(), 1e-12);

        Eigen::Matrix<double, 3, 2> expected;
        expected << 1, 0,
                    0, 1,
                    0, 0;
        REQUIRE_VEC_NEAR(Eigen::Map<const Eigen::VectorXd>(center.J0.data(), 6),
                         Eigen::Map<const Eigen::VectorXd>(expected.data(), 6),
                         1e-12);

        const double h = 1e-6;
        const Eigen::Vector2d u = center.u;
        int triId = -1;
        Eigen::Vector3d bPlusU, bMinusU, bPlusV, bMinusV;
        REQUIRE(chart.uvAccel.locate(u + Eigen::Vector2d(h, 0), triId, bPlusU));
        const auto xpU = ChartEvaluator::evaluate(chart, triId, bPlusU, pose).X0;
        REQUIRE(chart.uvAccel.locate(u - Eigen::Vector2d(h, 0), triId, bMinusU));
        const auto xmU = ChartEvaluator::evaluate(chart, triId, bMinusU, pose).X0;
        REQUIRE(chart.uvAccel.locate(u + Eigen::Vector2d(0, h), triId, bPlusV));
        const auto xpV = ChartEvaluator::evaluate(chart, triId, bPlusV, pose).X0;
        REQUIRE(chart.uvAccel.locate(u - Eigen::Vector2d(0, h), triId, bMinusV));
        const auto xmV = ChartEvaluator::evaluate(chart, triId, bMinusV, pose).X0;

        Eigen::Matrix<double, 3, 2> fd;
        fd.col(0) = (xpU - xmU) / (2.0 * h);
        fd.col(1) = (xpV - xmV) / (2.0 * h);
        REQUIRE_VEC_NEAR(Eigen::Map<const Eigen::VectorXd>(center.J0.data(), 6),
                         Eigen::Map<const Eigen::VectorXd>(fd.data(), 6),
                         1e-8);
    });
}

