#pragma once

#include "atlas/ContactAtlas.h"
#include "contact/ContactTypes.h"
#include "sdf/ISDF.h"

namespace bff_sdf::contact {

struct ContactIntegrationOptions {
    double zeroEps = 1e-12;
};

class ContactIntegrator {
public:
    static ContactResult integrateChart(const bff_sdf::atlas::ContactChart& sourceChart,
                                        const bff_sdf::core::RigidPose& sourcePose,
                                        const bff_sdf::sdf::ISDF& targetSdf,
                                        const bff_sdf::core::RigidPose& targetPose,
                                        const ContactIntegrationOptions& options = {});
};

} // namespace bff_sdf::contact

