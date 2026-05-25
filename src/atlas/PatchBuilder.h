#pragma once

#include "atlas/Patch.h"

namespace bff_sdf::atlas {

class PatchBuilder {
public:
    static Patch buildWholeMeshPatch(const bff_sdf::core::Mesh& mesh);
    static Patch buildFromFaces(const bff_sdf::core::Mesh& mesh, const std::vector<int>& faceIds);
    static bool validateDiskLike(Patch& patch);

private:
    static bool isConnectedByFaces(const bff_sdf::core::Mesh& mesh);
};

} // namespace bff_sdf::atlas

