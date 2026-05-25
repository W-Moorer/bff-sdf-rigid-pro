#include "contact/MarchingTriangles.h"

#include <algorithm>
#include <cmath>

namespace bff_sdf::contact {
namespace {

struct UVValue {
    Eigen::Vector2d u;
    double g = 0.0;
};

UVValue interpolateZero(const UVValue& a, const UVValue& b)
{
    const double denom = a.g - b.g;
    const double lambda = std::abs(denom) < 1e-30 ? 0.5 : a.g / denom;
    return {(1.0 - lambda) * a.u + lambda * b.u, 0.0};
}

bool containsPoint(const std::vector<Eigen::Vector2d>& points, const Eigen::Vector2d& p, double eps)
{
    for (const auto& q : points) {
        if ((q - p).norm() <= eps) return true;
    }
    return false;
}

} // namespace

double polygonArea2D(const std::vector<Eigen::Vector2d>& poly)
{
    if (poly.size() < 3) return 0.0;
    double twiceArea = 0.0;
    for (size_t i = 0; i < poly.size(); ++i) {
        const auto& a = poly[i];
        const auto& b = poly[(i + 1) % poly.size()];
        twiceArea += a.x() * b.y() - a.y() * b.x();
    }
    return 0.5 * twiceArea;
}

void cleanDuplicatePoints(std::vector<Eigen::Vector2d>& points, double eps)
{
    std::vector<Eigen::Vector2d> out;
    for (const auto& p : points) {
        if (!containsPoint(out, p, eps)) out.push_back(p);
    }
    points.swap(out);
}

MarchingTriangleResult MarchingTriangles::extract(const Eigen::Vector2d& u0,
                                                  const Eigen::Vector2d& u1,
                                                  const Eigen::Vector2d& u2,
                                                  double g0,
                                                  double g1,
                                                  double g2,
                                                  double eps)
{
    MarchingTriangleResult result;
    const std::vector<UVValue> input = {{u0, g0}, {u1, g1}, {u2, g2}};
    std::vector<UVValue> poly;

    for (int i = 0; i < 3; ++i) {
        const UVValue a = input[i];
        const UVValue b = input[(i + 1) % 3];
        const bool aInside = a.g <= eps;
        const bool bInside = b.g <= eps;

        if (aInside && bInside) {
            poly.push_back(b);
        } else if (aInside && !bInside) {
            poly.push_back(interpolateZero(a, b));
        } else if (!aInside && bInside) {
            poly.push_back(interpolateZero(a, b));
            poly.push_back(b);
        }
    }

    result.negativePolygon.reserve(poly.size());
    for (const auto& p : poly) result.negativePolygon.push_back(p.u);
    cleanDuplicatePoints(result.negativePolygon);
    result.hasContact = result.negativePolygon.size() >= 3 && std::abs(polygonArea2D(result.negativePolygon)) > eps;

    std::vector<Eigen::Vector2d> zeros;
    for (int i = 0; i < 3; ++i) {
        const UVValue a = input[i];
        const UVValue b = input[(i + 1) % 3];
        if (std::abs(a.g) <= eps) zeros.push_back(a.u);
        if ((a.g < -eps && b.g > eps) || (a.g > eps && b.g < -eps)) {
            zeros.push_back(interpolateZero(a, b).u);
        }
    }
    cleanDuplicatePoints(zeros);
    if (zeros.size() >= 2) {
        result.zeroSegment = {zeros[0], zeros[1]};
        result.hasZeroSegment = (zeros[0] - zeros[1]).norm() > eps;
    }

    return result;
}

} // namespace bff_sdf::contact

