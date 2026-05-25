#include "sdf/GridSDF.h"

#include <algorithm>
#include <cmath>

namespace bff_sdf::sdf {

GridSDF GridSDF::sampleFrom(const ISDF& sdf,
                            const bff_sdf::core::AABB& bounds,
                            int resolution,
                            double padding)
{
    GridSDF grid;
    grid.n_ = std::max(2, resolution);
    const Eigen::Vector3d extent = bounds.extent();
    const double maxExtent = extent.maxCoeff() + 2.0 * padding;
    grid.h_ = maxExtent / static_cast<double>(grid.n_ - 1);
    grid.origin_ = bounds.center() - 0.5 * Eigen::Vector3d::Constant(maxExtent);
    grid.phi_.resize(static_cast<std::size_t>(grid.n_) * grid.n_ * grid.n_);

    for (int k = 0; k < grid.n_; ++k) {
        for (int j = 0; j < grid.n_; ++j) {
            for (int i = 0; i < grid.n_; ++i) {
                const Eigen::Vector3d p = grid.origin_ + grid.h_ * Eigen::Vector3d(i, j, k);
                grid.phi_[static_cast<std::size_t>((k * grid.n_ + j) * grid.n_ + i)] =
                    sdf.evalLocal(p, false, false).phi;
            }
        }
    }
    return grid;
}

double GridSDF::phiAt(int i, int j, int k) const
{
    i = std::clamp(i, 0, n_ - 1);
    j = std::clamp(j, 0, n_ - 1);
    k = std::clamp(k, 0, n_ - 1);
    return phi_[static_cast<std::size_t>((k * n_ + j) * n_ + i)];
}

double GridSDF::trilinear(const Eigen::Vector3d& y) const
{
    Eigen::Vector3d p = (y - origin_) / h_;
    p = p.cwiseMax(Eigen::Vector3d::Zero()).cwiseMin(Eigen::Vector3d::Constant(static_cast<double>(n_ - 1)));
    const int i0 = std::clamp(static_cast<int>(std::floor(p.x())), 0, n_ - 2);
    const int j0 = std::clamp(static_cast<int>(std::floor(p.y())), 0, n_ - 2);
    const int k0 = std::clamp(static_cast<int>(std::floor(p.z())), 0, n_ - 2);
    const double tx = p.x() - i0;
    const double ty = p.y() - j0;
    const double tz = p.z() - k0;

    const double c000 = phiAt(i0, j0, k0);
    const double c100 = phiAt(i0 + 1, j0, k0);
    const double c010 = phiAt(i0, j0 + 1, k0);
    const double c110 = phiAt(i0 + 1, j0 + 1, k0);
    const double c001 = phiAt(i0, j0, k0 + 1);
    const double c101 = phiAt(i0 + 1, j0, k0 + 1);
    const double c011 = phiAt(i0, j0 + 1, k0 + 1);
    const double c111 = phiAt(i0 + 1, j0 + 1, k0 + 1);

    const double c00 = c000 * (1.0 - tx) + c100 * tx;
    const double c10 = c010 * (1.0 - tx) + c110 * tx;
    const double c01 = c001 * (1.0 - tx) + c101 * tx;
    const double c11 = c011 * (1.0 - tx) + c111 * tx;
    const double c0 = c00 * (1.0 - ty) + c10 * ty;
    const double c1 = c01 * (1.0 - ty) + c11 * ty;
    return c0 * (1.0 - tz) + c1 * tz;
}

SDFQuery GridSDF::evalLocal(const Eigen::Vector3d& y, bool needGrad, bool needClosest) const
{
    SDFQuery q;
    q.phi = trilinear(y);
    if (needGrad) {
        const double step = h_;
        q.grad = Eigen::Vector3d(
            (trilinear(y + step * Eigen::Vector3d::UnitX()) - trilinear(y - step * Eigen::Vector3d::UnitX())) / (2.0 * step),
            (trilinear(y + step * Eigen::Vector3d::UnitY()) - trilinear(y - step * Eigen::Vector3d::UnitY())) / (2.0 * step),
            (trilinear(y + step * Eigen::Vector3d::UnitZ()) - trilinear(y - step * Eigen::Vector3d::UnitZ())) / (2.0 * step));
        const double len = q.grad.norm();
        if (len > 1e-12 && std::isfinite(len)) q.grad /= len;
        else q.grad = Eigen::Vector3d::UnitX();
    }
    if (needClosest) {
        if (!needGrad) {
            q.grad = evalLocal(y, true, false).grad;
        }
        q.closestPoint = y - q.phi * q.grad;
        q.hasClosest = true;
    }
    return q;
}

} // namespace bff_sdf::sdf

