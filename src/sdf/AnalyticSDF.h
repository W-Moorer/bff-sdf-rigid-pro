#pragma once

#include "sdf/ISDF.h"

namespace bff_sdf::sdf {

class PlaneSDF final : public ISDF {
public:
    PlaneSDF(const Eigen::Vector3d& pointOnPlane, const Eigen::Vector3d& normal);
    SDFQuery evalLocal(const Eigen::Vector3d& y, bool needGrad, bool needClosest) const override;

private:
    Eigen::Vector3d point_;
    Eigen::Vector3d normal_;
};

class SphereSDF final : public ISDF {
public:
    SphereSDF(const Eigen::Vector3d& center, double radius);
    SDFQuery evalLocal(const Eigen::Vector3d& y, bool needGrad, bool needClosest) const override;

private:
    Eigen::Vector3d center_;
    double radius_;
};

class BoxSDF final : public ISDF {
public:
    BoxSDF(const Eigen::Vector3d& center, const Eigen::Vector3d& halfExtents);
    SDFQuery evalLocal(const Eigen::Vector3d& y, bool needGrad, bool needClosest) const override;

private:
    double phiOnly(const Eigen::Vector3d& y) const;

    Eigen::Vector3d center_;
    Eigen::Vector3d halfExtents_;
};

} // namespace bff_sdf::sdf

