#pragma once

#include "atlas/ContactAtlas.h"
#include "core/RigidPose.h"

namespace bff_sdf::atlas {

struct ChartEvaluation {
    bool valid = false;
    std::string message;
    Eigen::Vector3d X0 = Eigen::Vector3d::Zero();
    Eigen::Vector2d u = Eigen::Vector2d::Zero();
    Eigen::Matrix<double, 3, 2> J0 = Eigen::Matrix<double, 3, 2>::Zero();
    Eigen::Vector3d xWorld = Eigen::Vector3d::Zero();
    Eigen::Matrix<double, 3, 2> JWorld = Eigen::Matrix<double, 3, 2>::Zero();
    Eigen::Vector3d normalWorld = Eigen::Vector3d::UnitZ();
};

class ChartEvaluator {
public:
    static ChartEvaluation evaluate(const ContactChart& chart,
                                    int triId,
                                    const Eigen::Vector3d& bary,
                                    const bff_sdf::core::RigidPose& pose);

    static Eigen::Matrix<double, 3, 2> restJacobian(const ContactChart& chart,
                                                    int triId,
                                                    bool& ok);
};

} // namespace bff_sdf::atlas

