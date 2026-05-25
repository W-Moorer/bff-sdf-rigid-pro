#include "sdf/AnalyticSDF.h"

#include <algorithm>
#include <cmath>

namespace bff_sdf::sdf {

PlaneSDF::PlaneSDF(const Eigen::Vector3d& pointOnPlane, const Eigen::Vector3d& normal)
    : point_(pointOnPlane), normal_(normal.normalized())
{
}

SDFQuery PlaneSDF::evalLocal(const Eigen::Vector3d& y, bool needGrad, bool needClosest) const
{
    SDFQuery q;
    q.phi = normal_.dot(y - point_);
    if (needGrad) q.grad = normal_;
    if (needClosest) {
        q.closestPoint = y - q.phi * normal_;
        q.hasClosest = true;
    }
    return q;
}

SphereSDF::SphereSDF(const Eigen::Vector3d& center, double radius)
    : center_(center), radius_(radius)
{
}

SDFQuery SphereSDF::evalLocal(const Eigen::Vector3d& y, bool needGrad, bool needClosest) const
{
    SDFQuery q;
    const Eigen::Vector3d d = y - center_;
    const double len = d.norm();
    q.phi = len - radius_;
    Eigen::Vector3d n = Eigen::Vector3d::UnitX();
    if (len > 1e-14 && std::isfinite(len)) n = d / len;
    if (needGrad) q.grad = n;
    if (needClosest) {
        q.closestPoint = center_ + radius_ * n;
        q.hasClosest = true;
    }
    return q;
}

BoxSDF::BoxSDF(const Eigen::Vector3d& center, const Eigen::Vector3d& halfExtents)
    : center_(center), halfExtents_(halfExtents.cwiseAbs())
{
}

double BoxSDF::phiOnly(const Eigen::Vector3d& y) const
{
    const Eigen::Vector3d p = y - center_;
    const Eigen::Vector3d q = p.cwiseAbs() - halfExtents_;
    const Eigen::Vector3d outside = q.cwiseMax(Eigen::Vector3d::Zero());
    const double inside = std::min(std::max(q.x(), std::max(q.y(), q.z())), 0.0);
    return outside.norm() + inside;
}

SDFQuery BoxSDF::evalLocal(const Eigen::Vector3d& y, bool needGrad, bool needClosest) const
{
    SDFQuery q;
    q.phi = phiOnly(y);

    if (needGrad) {
        const double h = 1e-6;
        q.grad = Eigen::Vector3d(
            (phiOnly(y + h * Eigen::Vector3d::UnitX()) - phiOnly(y - h * Eigen::Vector3d::UnitX())) / (2.0 * h),
            (phiOnly(y + h * Eigen::Vector3d::UnitY()) - phiOnly(y - h * Eigen::Vector3d::UnitY())) / (2.0 * h),
            (phiOnly(y + h * Eigen::Vector3d::UnitZ()) - phiOnly(y - h * Eigen::Vector3d::UnitZ())) / (2.0 * h));
        const double len = q.grad.norm();
        if (len > 1e-14 && std::isfinite(len)) q.grad /= len;
        else q.grad = Eigen::Vector3d::UnitX();
    }

    if (needClosest) {
        const Eigen::Vector3d p = y - center_;
        Eigen::Vector3d clamped = p.cwiseMax(-halfExtents_).cwiseMin(halfExtents_);
        if (q.phi < 0.0) {
            Eigen::Vector3d distToFace = halfExtents_ - p.cwiseAbs();
            int axis = 0;
            if (distToFace.y() < distToFace.x()) axis = 1;
            if (distToFace.z() < distToFace[axis]) axis = 2;
            clamped = p;
            clamped[axis] = (p[axis] >= 0.0 ? halfExtents_[axis] : -halfExtents_[axis]);
        }
        q.closestPoint = center_ + clamped;
        q.hasClosest = true;
    }

    return q;
}

} // namespace bff_sdf::sdf

