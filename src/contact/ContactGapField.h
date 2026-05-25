#pragma once

#include "atlas/ChartEvaluator.h"
#include "contact/ContactTypes.h"
#include "sdf/ISDF.h"

namespace bff_sdf::contact {

class ContactGapField {
public:
    static GapSample sample(const bff_sdf::atlas::ContactChart& sourceChart,
                            int sourceTriId,
                            const Eigen::Vector3d& bary,
                            const bff_sdf::core::RigidPose& sourcePose,
                            const bff_sdf::sdf::ISDF& targetSdf,
                            const bff_sdf::core::RigidPose& targetPose);
};

} // namespace bff_sdf::contact

