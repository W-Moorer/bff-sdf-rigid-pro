#pragma once

#include "atlas/ContactAtlas.h"
#include "contact/ContactTypes.h"
#include "core/RigidPose.h"
#include "sdf/ISDF.h"

#include <string>

namespace bff_sdf::io {

bool ensureDirectory(const std::string& dir);
bool writeGapFieldCsv(const std::string& path,
                      const bff_sdf::atlas::ContactChart& chart,
                      const bff_sdf::core::RigidPose& sourcePose,
                      const bff_sdf::sdf::ISDF& targetSdf,
                      const bff_sdf::core::RigidPose& targetPose);
bool writeContactContourObj(const std::string& path, const bff_sdf::contact::ContactResult& result);
bool writeContactContourUvCsv(const std::string& path, const bff_sdf::contact::ContactResult& result);
bool writeContactSamplesObj(const std::string& path, const bff_sdf::contact::ContactResult& result);
bool writeContactSamplesPly(const std::string& path, const bff_sdf::contact::ContactResult& result);
bool writeAtlasDebugObj(const std::string& path,
                        const bff_sdf::atlas::ContactChart& chart,
                        const bff_sdf::core::RigidPose& pose);

} // namespace bff_sdf::io

