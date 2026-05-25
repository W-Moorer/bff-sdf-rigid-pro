#pragma once

#include "core/AABB.h"
#include "core/Mesh.h"
#include "sdf/ISDF.h"

#include <vector>

namespace bff_sdf::sdf {

struct TriangleClosestPoint {
    Eigen::Vector3d point = Eigen::Vector3d::Zero();
    Eigen::Vector3d bary = Eigen::Vector3d::Zero();
    double squaredDistance = 0.0;
};

TriangleClosestPoint closestPointOnTriangle(const Eigen::Vector3d& p,
                                            const Eigen::Vector3d& a,
                                            const Eigen::Vector3d& b,
                                            const Eigen::Vector3d& c);

class MeshSDF final : public ISDF {
public:
    enum class SignMode {
        Unsigned,
        RaycastClosedMesh
    };

    explicit MeshSDF(const bff_sdf::core::Mesh& mesh, SignMode signMode = SignMode::RaycastClosedMesh);

    SDFQuery evalLocal(const Eigen::Vector3d& y, bool needGrad, bool needClosest) const override;
    const bff_sdf::core::Mesh& mesh() const { return mesh_; }
    double signedness(const Eigen::Vector3d& y) const;

private:
    struct BVHNode {
        bff_sdf::core::AABB bounds;
        int left = -1;
        int right = -1;
        int begin = 0;
        int end = 0;
        bool leaf = false;
    };

    void buildBvh();
    int buildNode(int begin, int end);
    bff_sdf::core::AABB faceBounds(int faceId) const;
    double aabbDistanceSquared(const bff_sdf::core::AABB& box, const Eigen::Vector3d& p) const;
    void closestInNode(int nodeId,
                       const Eigen::Vector3d& y,
                       double& best,
                       TriangleClosestPoint& bestCp,
                       int& bestFace) const;

    bool rayIntersectsTriangle(const Eigen::Vector3d& origin,
                               const Eigen::Vector3d& a,
                               const Eigen::Vector3d& b,
                               const Eigen::Vector3d& c) const;

    bff_sdf::core::Mesh mesh_;
    SignMode signMode_ = SignMode::RaycastClosedMesh;
    std::vector<int> faceOrder_;
    std::vector<BVHNode> bvh_;
};

} // namespace bff_sdf::sdf
