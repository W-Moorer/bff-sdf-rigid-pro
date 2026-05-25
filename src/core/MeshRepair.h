#pragma once

#include "core/Mesh.h"

#include <string>

namespace bff_sdf::core {

struct MeshRepairOptions {
    double weldEpsilon = 1e-9;
    double degeneracyEpsRelative = 1e-14;
    bool removeDuplicateFaces = true;
    bool orientConnectedComponents = true;
    bool orientClosedMeshOutward = true;
};

struct MeshRepairReport {
    int inputVertices = 0;
    int inputFaces = 0;
    int outputVertices = 0;
    int outputFaces = 0;
    int removedDegenerateFaces = 0;
    int removedDuplicateFaces = 0;
    int flippedFaces = 0;
    int connectedComponents = 0;
    int nonManifoldEdges = 0;
    std::string message;
};

Mesh repairMesh(const Mesh& input,
                const MeshRepairOptions& options = {},
                MeshRepairReport* report = nullptr);

bool sharedEdgeOrientationConsistent(const Mesh& mesh);

} // namespace bff_sdf::core

