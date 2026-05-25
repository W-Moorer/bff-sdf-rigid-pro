#pragma once

#include <Eigen/Dense>
#include <vector>

namespace bff_sdf::atlas {

struct UVTriangleAABB {
    Eigen::Vector2d min = Eigen::Vector2d::Zero();
    Eigen::Vector2d max = Eigen::Vector2d::Zero();
};

class UVAcceleration {
public:
    void build(const Eigen::MatrixXd& UV, const Eigen::MatrixXi& Fuv);
    bool locate(const Eigen::Vector2d& u, int& triId, Eigen::Vector3d& bary) const;

private:
    Eigen::MatrixXd UV_;
    Eigen::MatrixXi Fuv_;
    std::vector<UVTriangleAABB> boxes_;
};

double signedUVArea(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c);
bool barycentricUV(const Eigen::Vector2d& p,
                   const Eigen::Vector2d& a,
                   const Eigen::Vector2d& b,
                   const Eigen::Vector2d& c,
                   Eigen::Vector3d& bary,
                   double eps = 1e-12);

} // namespace bff_sdf::atlas

