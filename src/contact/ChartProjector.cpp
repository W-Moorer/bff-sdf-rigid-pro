#include "contact/ChartProjector.h"

#include "sdf/MeshSDF.h"

#include <limits>

namespace bff_sdf::contact {

ProjectionResult ChartProjector::project(const Eigen::Vector3d& xWorld,
                                         const bff_sdf::core::Mesh& targetMesh,
                                         const bff_sdf::core::RigidPose& targetPose,
                                         const bff_sdf::sdf::SDFQuery* hint,
                                         const ProjectionOptions& options)
{
    counters_.projection_attempts++;
    ProjectionResult result;

    if (options.mode == ProjectionMode::AtlasHO) {
        counters_.fallback_to_bvh++;
        result.usedFallback = true;
        result.message = "Atlas-HO is not implemented; falling back to Atlas-Linear closest triangle";
    }

    const Eigen::Vector3d y = targetPose.worldToLocal(xWorld);
    double best = std::numeric_limits<double>::infinity();
    bff_sdf::sdf::TriangleClosestPoint bestCp;
    int bestTri = -1;

    auto checkFace = [&](int fi) {
        if (fi < 0 || fi >= targetMesh.F.rows()) return;
        const Eigen::Vector3i f = targetMesh.F.row(fi);
        const auto cp = bff_sdf::sdf::closestPointOnTriangle(y, targetMesh.V.row(f[0]), targetMesh.V.row(f[1]), targetMesh.V.row(f[2]));
        if (cp.squaredDistance < best) {
            best = cp.squaredDistance;
            bestCp = cp;
            bestTri = fi;
        }
    };

    if (options.useHintFace && hint && hint->closestFaceId >= 0) checkFace(hint->closestFaceId);
    if (bestTri < 0) {
        counters_.fallback_to_bvh++;
        for (int fi = 0; fi < targetMesh.F.rows(); ++fi) checkFace(fi);
    } else {
        for (int fi = 0; fi < targetMesh.F.rows(); ++fi) checkFace(fi);
    }

    if (bestTri < 0) {
        counters_.projection_failures++;
        result.message = "target mesh has no projectable triangles";
        return result;
    }

    const Eigen::Vector3i f = targetMesh.F.row(bestTri);
    const Eigen::Vector3d a = targetMesh.V.row(f[0]);
    const Eigen::Vector3d b = targetMesh.V.row(f[1]);
    const Eigen::Vector3d c = targetMesh.V.row(f[2]);
    Eigen::Vector3d nLocal = (b - a).cross(c - a);
    if (nLocal.norm() > 1e-12) nLocal.normalize();
    else nLocal = Eigen::Vector3d::UnitZ();

    result.closestPointWorld = targetPose.localToWorld(bestCp.point);
    result.normalWorld = targetPose.rotateLocalVectorToWorld(nLocal).normalized();
    const double distance = (xWorld - result.closestPointWorld).norm();
    const double sign = hint ? (hint->phi < 0.0 ? -1.0 : 1.0)
                             : ((xWorld - result.closestPointWorld).dot(result.normalWorld) < 0.0 ? -1.0 : 1.0);
    result.signedGap = sign * distance;
    result.triId = bestTri;
    result.bary = bestCp.bary;
    result.iterations = 1;
    result.converged = true;
    if (result.message.empty()) result.message = "Atlas-Linear closest triangle projection";
    counters_.projection_converged++;
    return result;
}

} // namespace bff_sdf::contact
