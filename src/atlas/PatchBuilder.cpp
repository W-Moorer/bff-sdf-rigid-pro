#include "atlas/PatchBuilder.h"

#include <map>
#include <queue>
#include <set>

namespace bff_sdf::atlas {

Patch PatchBuilder::buildWholeMeshPatch(const bff_sdf::core::Mesh& mesh)
{
    std::vector<int> faces(mesh.F.rows());
    for (int i = 0; i < mesh.F.rows(); ++i) faces[i] = i;
    return buildFromFaces(mesh, faces);
}

Patch PatchBuilder::buildFromFaces(const bff_sdf::core::Mesh& mesh, const std::vector<int>& faceIds)
{
    Patch patch;
    patch.globalFaceIds = faceIds;

    std::map<int, int> globalToLocal;
    for (int fid : faceIds) {
        for (int k = 0; k < 3; ++k) {
            const int gv = mesh.F(fid, k);
            if (!globalToLocal.count(gv)) {
                const int lid = static_cast<int>(globalToLocal.size());
                globalToLocal[gv] = lid;
                patch.globalVertexIds.push_back(gv);
            }
        }
    }

    patch.localMesh.V.resize((int)patch.globalVertexIds.size(), 3);
    for (int i = 0; i < (int)patch.globalVertexIds.size(); ++i) {
        patch.localMesh.V.row(i) = mesh.V.row(patch.globalVertexIds[i]);
    }

    patch.localMesh.F.resize((int)faceIds.size(), 3);
    for (int i = 0; i < (int)faceIds.size(); ++i) {
        for (int k = 0; k < 3; ++k) {
            patch.localMesh.F(i, k) = globalToLocal[mesh.F(faceIds[i], k)];
        }
    }
    patch.localMesh.computeNormals();
    validateDiskLike(patch);
    return patch;
}

bool PatchBuilder::isConnectedByFaces(const bff_sdf::core::Mesh& mesh)
{
    if (mesh.F.rows() == 0) return false;
    const auto adjacency = mesh.faceAdjacency();
    std::vector<int> visited(mesh.F.rows(), 0);
    std::queue<int> q;
    q.push(0);
    visited[0] = 1;
    int count = 0;
    while (!q.empty()) {
        const int f = q.front();
        q.pop();
        ++count;
        for (int e = 0; e < 3; ++e) {
            const int nb = adjacency[f][e];
            if (nb >= 0 && !visited[nb]) {
                visited[nb] = 1;
                q.push(nb);
            }
        }
    }
    return count == mesh.F.rows();
}

bool PatchBuilder::validateDiskLike(Patch& patch)
{
    patch.eulerCharacteristic = patch.localMesh.eulerCharacteristic();
    patch.connected = isConnectedByFaces(patch.localMesh);
    const auto loops = patch.localMesh.boundaryLoops();
    if (!loops.empty()) patch.boundaryLoop = loops.front();

    patch.diskLike = patch.connected && loops.size() == 1 && patch.eulerCharacteristic == 1;
    if (!patch.connected) patch.diagnostic = "patch is not face-connected";
    else if (loops.size() != 1) patch.diagnostic = "patch must have exactly one boundary loop";
    else if (patch.eulerCharacteristic != 1) patch.diagnostic = "patch Euler characteristic is not 1";
    else patch.diagnostic = "disk-like patch";

    return patch.diskLike;
}

} // namespace bff_sdf::atlas

