#include "core/MeshRepair.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <vector>

namespace bff_sdf::core {
namespace {

struct EdgeIncident {
    int face = -1;
    int edge = -1;
    int sign = 0;
};

std::array<int, 3> sortedFaceKey(const Eigen::Vector3i& f)
{
    std::array<int, 3> key = {f[0], f[1], f[2]};
    std::sort(key.begin(), key.end());
    return key;
}

int edgeSignInFace(const Eigen::Vector3i& f, int edge)
{
    const int a = f[edge];
    const int b = f[(edge + 1) % 3];
    const auto ce = canonicalEdge(a, b);
    return (a == ce.first && b == ce.second) ? 1 : -1;
}

double signedVolume(const Mesh& mesh)
{
    double volume = 0.0;
    for (int i = 0; i < mesh.F.rows(); ++i) {
        const Eigen::Vector3i f = mesh.F.row(i);
        const Eigen::Vector3d a = mesh.V.row(f[0]);
        const Eigen::Vector3d b = mesh.V.row(f[1]);
        const Eigen::Vector3d c = mesh.V.row(f[2]);
        volume += a.dot(b.cross(c)) / 6.0;
    }
    return volume;
}

} // namespace

Mesh repairMesh(const Mesh& input, const MeshRepairOptions& options, MeshRepairReport* report)
{
    MeshRepairReport local;
    local.inputVertices = static_cast<int>(input.V.rows());
    local.inputFaces = static_cast<int>(input.F.rows());

    Mesh welded;
    std::map<std::array<long long, 3>, int> vertexMap;
    std::vector<Eigen::Vector3d> vertices;
    const double invEps = options.weldEpsilon > 0.0 ? 1.0 / options.weldEpsilon : 1e12;
    auto mapVertex = [&](const Eigen::Vector3d& p) {
        const std::array<long long, 3> key = {
            static_cast<long long>(std::llround(p.x() * invEps)),
            static_cast<long long>(std::llround(p.y() * invEps)),
            static_cast<long long>(std::llround(p.z() * invEps))
        };
        auto it = vertexMap.find(key);
        if (it != vertexMap.end()) return it->second;
        const int id = static_cast<int>(vertices.size());
        vertexMap[key] = id;
        vertices.push_back(p);
        return id;
    };

    std::vector<Eigen::Vector3i> faces;
    std::set<std::array<int, 3>> seenFaces;
    const double extent = input.bounds().extent().maxCoeff();
    const double areaEps = options.degeneracyEpsRelative * std::max(1.0, extent * extent);

    for (int i = 0; i < input.F.rows(); ++i) {
        Eigen::Vector3i f;
        for (int k = 0; k < 3; ++k) f[k] = mapVertex(input.V.row(input.F(i, k)));
        if (f[0] == f[1] || f[1] == f[2] || f[2] == f[0]) {
            local.removedDegenerateFaces++;
            continue;
        }
        const Eigen::Vector3d a = vertices[f[0]];
        const Eigen::Vector3d b = vertices[f[1]];
        const Eigen::Vector3d c = vertices[f[2]];
        if (0.5 * (b - a).cross(c - a).norm() <= areaEps) {
            local.removedDegenerateFaces++;
            continue;
        }
        const auto key = sortedFaceKey(f);
        if (options.removeDuplicateFaces && seenFaces.count(key)) {
            local.removedDuplicateFaces++;
            continue;
        }
        seenFaces.insert(key);
        faces.push_back(f);
    }

    std::vector<int> used(vertices.size(), 0);
    for (const auto& f : faces) {
        used[f[0]] = 1;
        used[f[1]] = 1;
        used[f[2]] = 1;
    }

    std::vector<int> remap(vertices.size(), -1);
    std::vector<Eigen::Vector3d> compactVertices;
    for (int i = 0; i < static_cast<int>(vertices.size()); ++i) {
        if (!used[i]) continue;
        remap[i] = static_cast<int>(compactVertices.size());
        compactVertices.push_back(vertices[i]);
    }
    for (auto& f : faces) {
        f[0] = remap[f[0]];
        f[1] = remap[f[1]];
        f[2] = remap[f[2]];
    }

    welded.V.resize(static_cast<int>(compactVertices.size()), 3);
    for (int i = 0; i < static_cast<int>(compactVertices.size()); ++i) welded.V.row(i) = compactVertices[i].transpose();
    welded.F.resize(static_cast<int>(faces.size()), 3);
    for (int i = 0; i < static_cast<int>(faces.size()); ++i) welded.F.row(i) = faces[i].transpose();

    if (options.orientConnectedComponents && welded.F.rows() > 0) {
        std::map<std::pair<int, int>, std::vector<EdgeIncident>> incidents;
        for (int fi = 0; fi < welded.F.rows(); ++fi) {
            const Eigen::Vector3i f = welded.F.row(fi);
            for (int e = 0; e < 3; ++e) {
                incidents[canonicalEdge(f[e], f[(e + 1) % 3])].push_back({fi, e, edgeSignInFace(f, e)});
            }
        }

        std::vector<std::vector<std::pair<int, int>>> adjacency(welded.F.rows());
        for (const auto& kv : incidents) {
            const auto& list = kv.second;
            if (list.size() > 2) {
                local.nonManifoldEdges++;
                continue;
            }
            if (list.size() != 2) continue;
            adjacency[list[0].face].push_back({list[1].face, list[0].sign == list[1].sign ? 1 : 0});
            adjacency[list[1].face].push_back({list[0].face, list[0].sign == list[1].sign ? 1 : 0});
        }

        std::vector<int> visited(welded.F.rows(), 0);
        std::vector<int> flip(welded.F.rows(), 0);
        for (int seed = 0; seed < welded.F.rows(); ++seed) {
            if (visited[seed]) continue;
            local.connectedComponents++;
            std::queue<int> q;
            q.push(seed);
            visited[seed] = 1;
            while (!q.empty()) {
                const int f = q.front();
                q.pop();
                for (const auto& nb : adjacency[f]) {
                    const int g = nb.first;
                    const int needOppositeFlip = nb.second;
                    const int desiredFlip = flip[f] ^ needOppositeFlip;
                    if (!visited[g]) {
                        visited[g] = 1;
                        flip[g] = desiredFlip;
                        q.push(g);
                    }
                }
            }
        }

        for (int fi = 0; fi < welded.F.rows(); ++fi) {
            if (!flip[fi]) continue;
            std::swap(welded.F(fi, 1), welded.F(fi, 2));
            local.flippedFaces++;
        }
    }

    welded.computeNormals();
    if (options.orientClosedMeshOutward && welded.boundaryEdges().empty() && signedVolume(welded) < 0.0) {
        for (int fi = 0; fi < welded.F.rows(); ++fi) std::swap(welded.F(fi, 1), welded.F(fi, 2));
        local.flippedFaces += static_cast<int>(welded.F.rows());
        welded.computeNormals();
    }

    local.outputVertices = static_cast<int>(welded.V.rows());
    local.outputFaces = static_cast<int>(welded.F.rows());
    local.message = "mesh repair completed";
    if (report) *report = local;
    return welded;
}

bool sharedEdgeOrientationConsistent(const Mesh& mesh)
{
    std::map<std::pair<int, int>, std::vector<int>> signs;
    for (int fi = 0; fi < mesh.F.rows(); ++fi) {
        const Eigen::Vector3i f = mesh.F.row(fi);
        for (int e = 0; e < 3; ++e) {
            signs[canonicalEdge(f[e], f[(e + 1) % 3])].push_back(edgeSignInFace(f, e));
        }
    }
    for (const auto& kv : signs) {
        if (kv.second.size() != 2) continue;
        if (kv.second[0] == kv.second[1]) return false;
    }
    return true;
}

} // namespace bff_sdf::core
