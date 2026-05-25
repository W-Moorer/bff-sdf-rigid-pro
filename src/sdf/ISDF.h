#pragma once

#include <Eigen/Dense>

namespace bff_sdf::sdf {

struct SDFQuery {
    double phi = 0.0;
    Eigen::Vector3d grad = Eigen::Vector3d::Zero();
    Eigen::Vector3d closestPoint = Eigen::Vector3d::Zero();
    int closestFaceId = -1;
    bool hasClosest = false;
};

class ISDF {
public:
    virtual ~ISDF() = default;

    // Coordinates: y is in the target object's local/rest frame.
    virtual SDFQuery evalLocal(const Eigen::Vector3d& y,
                               bool needGrad,
                               bool needClosest) const = 0;
};

} // namespace bff_sdf::sdf

