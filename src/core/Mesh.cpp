#include "core/Mesh.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <stdexcept>

namespace bff_sdf::core {

std::pair<int, int> canonicalEdge(int a, int b)
{
    return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
}

void Mesh::computeFaceNormals()
{
    FN.resize(F.rows(), 3);
    for (int i = 0; i < F.rows(); ++i) {
        const Eigen::Vector3d a = V.row(F(i, 0));
        const Eigen::Vector3d b = V.row(F(i, 1));
        const Eigen::Vector3d c = V.row(F(i, 2));
        Eigen::Vector3d n = (b - a).cross(c - a);
        const double len = n.norm();
        if (len > 0.0 && std::isfinite(len)) n /= len;
        else n.setZero();
        FN.row(i) = n.transpose();
    }
}

void Mesh::computeVertexNormals()
{
    if (FN.rows() != F.rows()) computeFaceNormals();
    VN = Eigen::MatrixXd::Zero(V.rows(), 3);
    for (int i = 0; i < F.rows(); ++i) {
        const Eigen::Vector3d n = FN.row(i);
        for (int k = 0; k < 3; ++k) VN.row(F(i, k)) += n.transpose();
    }
    for (int i = 0; i < VN.rows(); ++i) {
        Eigen::Vector3d n = VN.row(i);
        const double len = n.norm();
        if (len > 0.0 && std::isfinite(len)) n /= len;
        else n = Eigen::Vector3d::UnitZ();
        VN.row(i) = n.transpose();
    }
}

void Mesh::computeNormals()
{
    computeFaceNormals();
    computeVertexNormals();
}

AABB Mesh::bounds() const
{
    AABB box;
    for (int i = 0; i < V.rows(); ++i) box.expand(V.row(i));
    return box;
}

std::map<std::pair<int, int>, int> Mesh::edgeUseCounts() const
{
    std::map<std::pair<int, int>, int> counts;
    for (int i = 0; i < F.rows(); ++i) {
        for (int e = 0; e < 3; ++e) {
            const int a = F(i, e);
            const int b = F(i, (e + 1) % 3);
            counts[canonicalEdge(a, b)]++;
        }
    }
    return counts;
}

std::vector<std::pair<int, int>> Mesh::boundaryEdges() const
{
    std::vector<std::pair<int, int>> out;
    const auto counts = edgeUseCounts();
    for (const auto& kv : counts) {
        if (kv.second == 1) out.push_back(kv.first);
    }
    return out;
}

std::vector<std::vector<int>> Mesh::boundaryLoops() const
{
    const auto edges = boundaryEdges();
    std::map<int, std::vector<int>> adj;
    for (const auto& e : edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }

    std::set<std::pair<int, int>> used;
    std::vector<std::vector<int>> loops;

    for (const auto& edge : edges) {
        if (used.count(edge)) continue;
        std::vector<int> loop;
        int start = edge.first;
        int prev = -1;
        int cur = start;
        while (true) {
            loop.push_back(cur);
            int next = -1;
            for (int candidate : adj[cur]) {
                const auto ce = canonicalEdge(cur, candidate);
                if (candidate != prev && !used.count(ce)) {
                    next = candidate;
                    break;
                }
            }
            if (next < 0) {
                for (int candidate : adj[cur]) {
                    const auto ce = canonicalEdge(cur, candidate);
                    if (!used.count(ce)) {
                        next = candidate;
                        break;
                    }
                }
            }
            if (next < 0) break;
            used.insert(canonicalEdge(cur, next));
            prev = cur;
            cur = next;
            if (cur == start) break;
            if (loop.size() > static_cast<size_t>(edges.size() + 1)) break;
        }
        if (!loop.empty()) loops.push_back(loop);
    }

    return loops;
}

std::vector<Eigen::Vector3i> Mesh::faceAdjacency() const
{
    std::vector<Eigen::Vector3i> adjacency(F.rows(), Eigen::Vector3i::Constant(-1));
    std::map<std::pair<int, int>, std::pair<int, int>> owner;
    for (int fi = 0; fi < F.rows(); ++fi) {
        for (int e = 0; e < 3; ++e) {
            const auto key = canonicalEdge(F(fi, e), F(fi, (e + 1) % 3));
            auto it = owner.find(key);
            if (it == owner.end()) {
                owner[key] = {fi, e};
            } else {
                adjacency[fi][e] = it->second.first;
                adjacency[it->second.first][it->second.second] = fi;
            }
        }
    }
    return adjacency;
}

int Mesh::eulerCharacteristic() const
{
    return static_cast<int>(V.rows()) - static_cast<int>(edgeUseCounts().size()) +
           static_cast<int>(F.rows());
}

bool Mesh::normalsFiniteAndUnit(double tol) const
{
    auto checkRows = [tol](const Eigen::MatrixXd& M) {
        for (int i = 0; i < M.rows(); ++i) {
            const Eigen::Vector3d n = M.row(i);
            if (!n.array().isFinite().all()) return false;
            if (std::abs(n.norm() - 1.0) > tol) return false;
        }
        return true;
    };
    return checkRows(FN) && checkRows(VN);
}

} // namespace bff_sdf::core

