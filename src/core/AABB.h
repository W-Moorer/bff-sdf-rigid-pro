#pragma once

#include <Eigen/Dense>
#include <limits>

namespace bff_sdf::core {

struct AABB {
    Eigen::Vector3d min = Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity());
    Eigen::Vector3d max = Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity());

    void expand(const Eigen::Vector3d& p)
    {
        min = min.cwiseMin(p);
        max = max.cwiseMax(p);
    }

    Eigen::Vector3d center() const { return 0.5 * (min + max); }
    Eigen::Vector3d extent() const { return max - min; }

    bool valid() const
    {
        return min.array().isFinite().all() && max.array().isFinite().all() &&
               (min.array() <= max.array()).all();
    }
};

} // namespace bff_sdf::core

