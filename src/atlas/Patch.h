#pragma once

#include "core/Mesh.h"

#include <string>
#include <vector>

namespace bff_sdf::atlas {

struct Patch {
    std::vector<int> globalFaceIds;
    std::vector<int> globalVertexIds;
    bff_sdf::core::Mesh localMesh;
    std::vector<int> boundaryLoop;
    bool connected = false;
    bool diskLike = false;
    int eulerCharacteristic = 0;
    std::string diagnostic;
};

} // namespace bff_sdf::atlas

