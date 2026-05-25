#pragma once

#include <Eigen/Dense>
#include <vector>

namespace bff_sdf::contact {

struct MarchingTriangleResult {
    std::vector<Eigen::Vector2d> negativePolygon;
    std::vector<Eigen::Vector2d> zeroSegment;
    bool hasContact = false;
    bool hasZeroSegment = false;
};

class MarchingTriangles {
public:
    static MarchingTriangleResult extract(const Eigen::Vector2d& u0,
                                          const Eigen::Vector2d& u1,
                                          const Eigen::Vector2d& u2,
                                          double g0,
                                          double g1,
                                          double g2,
                                          double eps = 1e-12);
};

double polygonArea2D(const std::vector<Eigen::Vector2d>& poly);
void cleanDuplicatePoints(std::vector<Eigen::Vector2d>& points, double eps = 1e-10);

} // namespace bff_sdf::contact

