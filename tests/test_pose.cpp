#include "test_main.h"

#include "core/RigidPose.h"

using bff_sdf::core::RigidPose;

int main()
{
    return test::run("test_pose", [] {
        RigidPose pose;
        const double theta = 0.37;
        pose.R = Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitZ()).toRotationMatrix();
        pose.t = Eigen::Vector3d(1.0, -2.0, 0.5);

        const Eigen::Vector3d X0(0.2, 0.4, -0.7);
        const Eigen::Vector3d xWorld = pose.localToWorld(X0);
        REQUIRE_VEC_NEAR(pose.worldToLocal(xWorld), X0, 1e-12);

        const Eigen::Vector3d gLocal(0.0, 1.0, 0.0);
        const Eigen::Vector3d gWorld = pose.rotateLocalVectorToWorld(gLocal);
        REQUIRE_NEAR(gWorld.norm(), 1.0, 1e-12);
        REQUIRE_VEC_NEAR(pose.rotateWorldVectorToLocal(gWorld), gLocal, 1e-12);
    });
}

