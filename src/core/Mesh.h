#pragma once

#include "core/AABB.h"

#include <Eigen/Dense>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace bff_sdf::core {

struct Mesh {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd VN;
    Eigen::MatrixXd FN;

    bool empty() const { return V.rows() == 0 || F.rows() == 0; }
    int vertexCount() const { return static_cast<int>(V.rows()); }
    int faceCount() const { return static_cast<int>(F.rows()); }

    void computeFaceNormals();
    void computeVertexNormals();
    void computeNormals();

    AABB bounds() const;
    std::map<std::pair<int, int>, int> edgeUseCounts() const;
    std::vector<std::pair<int, int>> boundaryEdges() const;
    std::vector<std::vector<int>> boundaryLoops() const;
    std::vector<Eigen::Vector3i> faceAdjacency() const;
    int eulerCharacteristic() const;
    bool normalsFiniteAndUnit(double tol = 1e-10) const;
};

std::pair<int, int> canonicalEdge(int a, int b);

} // namespace bff_sdf::core

