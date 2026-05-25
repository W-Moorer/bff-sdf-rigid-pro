#include "atlas/ChartEvaluator.h"

#include <cmath>

namespace bff_sdf::atlas {

Eigen::Matrix<double, 3, 2> ChartEvaluator::restJacobian(const ContactChart& chart,
                                                         int triId,
                                                         bool& ok)
{
    ok = false;
    Eigen::Matrix<double, 3, 2> J = Eigen::Matrix<double, 3, 2>::Zero();
    if (triId < 0 || triId >= chart.patch.localMesh.F.rows()) return J;

    const Eigen::Vector3i f = chart.patch.localMesh.F.row(triId);
    const Eigen::Vector3i fu = chart.Fuv.row(triId);
    const Eigen::Vector3d X0 = chart.patch.localMesh.V.row(f[0]);
    const Eigen::Vector3d X1 = chart.patch.localMesh.V.row(f[1]);
    const Eigen::Vector3d X2 = chart.patch.localMesh.V.row(f[2]);
    const Eigen::Vector2d u0 = chart.UV.row(fu[0]);
    const Eigen::Vector2d u1 = chart.UV.row(fu[1]);
    const Eigen::Vector2d u2 = chart.UV.row(fu[2]);

    Eigen::Matrix2d dU;
    dU.col(0) = u1 - u0;
    dU.col(1) = u2 - u0;
    if (std::abs(dU.determinant()) < 1e-14) return J;

    Eigen::Matrix<double, 3, 2> dX;
    dX.col(0) = X1 - X0;
    dX.col(1) = X2 - X0;
    J = dX * dU.inverse();
    ok = J.array().isFinite().all();
    return J;
}

ChartEvaluation ChartEvaluator::evaluate(const ContactChart& chart,
                                          int triId,
                                          const Eigen::Vector3d& bary,
                                          const bff_sdf::core::RigidPose& pose)
{
    ChartEvaluation ev;
    if (triId < 0 || triId >= chart.patch.localMesh.F.rows()) {
        ev.message = "triangle id out of range";
        return ev;
    }
    if (std::abs(bary.sum() - 1.0) > 1e-8) {
        ev.message = "barycentric coordinates do not sum to 1";
        return ev;
    }

    const Eigen::Vector3i f = chart.patch.localMesh.F.row(triId);
    const Eigen::Vector3i fu = chart.Fuv.row(triId);
    for (int k = 0; k < 3; ++k) {
        ev.X0 += bary[k] * chart.patch.localMesh.V.row(f[k]).transpose();
        ev.u += bary[k] * chart.UV.row(fu[k]).transpose();
    }

    bool ok = false;
    ev.J0 = restJacobian(chart, triId, ok);
    if (!ok) {
        ev.message = "degenerate UV triangle";
        return ev;
    }

    ev.xWorld = pose.localToWorld(ev.X0);
    ev.JWorld = pose.R * ev.J0;

    Eigen::Vector3d n = ev.JWorld.col(0).cross(ev.JWorld.col(1));
    if (n.norm() < 1e-14) {
        const Eigen::Vector3d localN = chart.patch.localMesh.FN.row(triId);
        n = pose.rotateLocalVectorToWorld(localN);
    }
    if (n.norm() < 1e-14) {
        ev.message = "degenerate chart normal";
        return ev;
    }
    ev.normalWorld = n.normalized();
    ev.valid = true;
    ev.message = "ok";
    return ev;
}

} // namespace bff_sdf::atlas

