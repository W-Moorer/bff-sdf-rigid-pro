#pragma once

#include "core/Mesh.h"
#include "sdf/ISDF.h"

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
    bool rayIntersectsTriangle(const Eigen::Vector3d& origin,
                               const Eigen::Vector3d& a,
                               const Eigen::Vector3d& b,
                               const Eigen::Vector3d& c) const;

    bff_sdf::core::Mesh mesh_;
    SignMode signMode_ = SignMode::RaycastClosedMesh;
};

} // namespace bff_sdf::sdf

