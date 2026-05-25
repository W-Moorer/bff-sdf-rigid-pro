#pragma once

#include "atlas/ContactAtlas.h"
#include "contact/ContactTypes.h"
#include "core/Mesh.h"
#include "sdf/ISDF.h"

#include <string>

namespace bff_sdf::contact {

struct RigidObject {
    std::string name;
    const bff_sdf::core::Mesh* mesh = nullptr;
    const bff_sdf::atlas::ContactAtlas* atlas = nullptr;
    const bff_sdf::sdf::ISDF* sdf = nullptr;
    bff_sdf::core::RigidPose pose;
};

struct ContactDetectorOptions {
    double contactMargin = 0.02;
    bool refineProjection = true;
    int sampleStride = 8;
};

class ContactDetector {
public:
    ContactResult detect(const RigidObject& source,
                         const RigidObject& target,
                         const ContactDetectorOptions& options = {});
};

} // namespace bff_sdf::contact

