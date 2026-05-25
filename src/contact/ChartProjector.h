#pragma once

#include "core/Mesh.h"
#include "core/RigidPose.h"
#include "sdf/ISDF.h"

#include <string>

namespace bff_sdf::contact {

enum class ProjectionMode {
    AtlasLinear,
    AtlasHO
};

struct ProjectionOptions {
    ProjectionMode mode = ProjectionMode::AtlasLinear;
    bool useHintFace = true;
};

struct ProjectionCounters {
    int projection_attempts = 0;
    int projection_converged = 0;
    int projection_failures = 0;
    int fallback_to_bvh = 0;
    int seam_crossing_count = 0;
};

struct ProjectionResult {
    bool converged = false;
    bool usedFallback = false;
    Eigen::Vector3d closestPointWorld = Eigen::Vector3d::Zero();
    Eigen::Vector3d normalWorld = Eigen::Vector3d::UnitZ();
    double signedGap = 0.0;
    int chartId = -1;
    int triId = -1;
    Eigen::Vector3d bary = Eigen::Vector3d::Zero();
    int iterations = 0;
    std::string message;
};

class ChartProjector {
public:
    ProjectionResult project(const Eigen::Vector3d& xWorld,
                             const bff_sdf::core::Mesh& targetMesh,
                             const bff_sdf::core::RigidPose& targetPose,
                             const bff_sdf::sdf::SDFQuery* hint,
                             const ProjectionOptions& options = {});

    const ProjectionCounters& counters() const { return counters_; }

private:
    ProjectionCounters counters_;
};

} // namespace bff_sdf::contact
