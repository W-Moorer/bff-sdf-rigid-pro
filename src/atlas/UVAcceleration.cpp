#include "atlas/UVAcceleration.h"

#include <cmath>

namespace bff_sdf::atlas {

double signedUVArea(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c)
{
    const Eigen::Vector2d ab = b - a;
    const Eigen::Vector2d ac = c - a;
    return 0.5 * (ab.x() * ac.y() - ab.y() * ac.x());
}

bool barycentricUV(const Eigen::Vector2d& p,
                   const Eigen::Vector2d& a,
                   const Eigen::Vector2d& b,
                   const Eigen::Vector2d& c,
                   Eigen::Vector3d& bary,
                   double eps)
{
    const double det = (b.y() - c.y()) * (a.x() - c.x()) + (c.x() - b.x()) * (a.y() - c.y());
    if (std::abs(det) < eps) return false;
    bary.x() = ((b.y() - c.y()) * (p.x() - c.x()) + (c.x() - b.x()) * (p.y() - c.y())) / det;
    bary.y() = ((c.y() - a.y()) * (p.x() - c.x()) + (a.x() - c.x()) * (p.y() - c.y())) / det;
    bary.z() = 1.0 - bary.x() - bary.y();
    return true;
}

void UVAcceleration::build(const Eigen::MatrixXd& UV, const Eigen::MatrixXi& Fuv)
{
    UV_ = UV;
    Fuv_ = Fuv;
    boxes_.resize(Fuv.rows());
    for (int i = 0; i < Fuv.rows(); ++i) {
        UVTriangleAABB box;
        box.min = UV.row(Fuv(i, 0));
        box.max = box.min;
        for (int k = 1; k < 3; ++k) {
            const Eigen::Vector2d u = UV.row(Fuv(i, k));
            box.min = box.min.cwiseMin(u);
            box.max = box.max.cwiseMax(u);
        }
        boxes_[i] = box;
    }
}

bool UVAcceleration::locate(const Eigen::Vector2d& u, int& triId, Eigen::Vector3d& bary) const
{
    constexpr double eps = 1e-10;
    for (int i = 0; i < Fuv_.rows(); ++i) {
        const auto& box = boxes_[i];
        if ((u.array() < (box.min.array() - eps)).any() || (u.array() > (box.max.array() + eps)).any()) {
            continue;
        }
        const Eigen::Vector2d a = UV_.row(Fuv_(i, 0));
        const Eigen::Vector2d b = UV_.row(Fuv_(i, 1));
        const Eigen::Vector2d c = UV_.row(Fuv_(i, 2));
        if (barycentricUV(u, a, b, c, bary) &&
            bary.minCoeff() >= -eps && bary.maxCoeff() <= 1.0 + eps) {
            triId = i;
            return true;
        }
    }
    triId = -1;
    bary.setZero();
    return false;
}

} // namespace bff_sdf::atlas

