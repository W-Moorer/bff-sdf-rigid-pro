#pragma once

#include "core/AABB.h"
#include "sdf/ISDF.h"

#include <vector>

namespace bff_sdf::sdf {

class GridSDF final : public ISDF {
public:
    static GridSDF sampleFrom(const ISDF& sdf,
                              const bff_sdf::core::AABB& bounds,
                              int resolution,
                              double padding);

    SDFQuery evalLocal(const Eigen::Vector3d& y, bool needGrad, bool needClosest) const override;
    double voxelSize() const { return h_; }
    int resolution() const { return n_; }
    std::size_t memoryBytes() const { return phi_.size() * sizeof(double); }

private:
    double phiAt(int i, int j, int k) const;
    double trilinear(const Eigen::Vector3d& y) const;

    Eigen::Vector3d origin_ = Eigen::Vector3d::Zero();
    double h_ = 1.0;
    int n_ = 0;
    std::vector<double> phi_;
};

} // namespace bff_sdf::sdf

