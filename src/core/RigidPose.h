#pragma once

#include <Eigen/Dense>

namespace bff_sdf::core {

struct RigidPose {
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    Eigen::Vector3d t = Eigen::Vector3d::Zero();

    static RigidPose identity() { return {}; }

    Eigen::Vector3d localToWorld(const Eigen::Vector3d& X0) const
    {
        return R * X0 + t;
    }

    Eigen::Vector3d worldToLocal(const Eigen::Vector3d& xWorld) const
    {
        return R.transpose() * (xWorld - t);
    }

    Eigen::Vector3d rotateLocalVectorToWorld(const Eigen::Vector3d& vLocal) const
    {
        return R * vLocal;
    }

    Eigen::Vector3d rotateWorldVectorToLocal(const Eigen::Vector3d& vWorld) const
    {
        return R.transpose() * vWorld;
    }
};

} // namespace bff_sdf::core

