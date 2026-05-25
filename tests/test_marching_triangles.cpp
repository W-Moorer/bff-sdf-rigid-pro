#include "test_main.h"

#include "contact/MarchingTriangles.h"

using bff_sdf::contact::MarchingTriangles;
using bff_sdf::contact::polygonArea2D;

int main()
{
    return test::run("test_marching_triangles", [] {
        const Eigen::Vector2d a(0, 0), b(1, 0), c(0, 1);

        auto allPositive = MarchingTriangles::extract(a, b, c, 1.0, 2.0, 3.0);
        REQUIRE(!allPositive.hasContact);
        REQUIRE(allPositive.negativePolygon.empty());

        auto allNegative = MarchingTriangles::extract(a, b, c, -1.0, -2.0, -3.0);
        REQUIRE(allNegative.hasContact);
        REQUIRE(allNegative.negativePolygon.size() == 3);
        REQUIRE_NEAR(std::abs(polygonArea2D(allNegative.negativePolygon)), 0.5, 1e-12);

        auto oneNegative = MarchingTriangles::extract(a, b, c, -1.0, 1.0, 1.0);
        REQUIRE(oneNegative.hasContact);
        REQUIRE(oneNegative.negativePolygon.size() == 3);
        REQUIRE(oneNegative.hasZeroSegment);
        REQUIRE_VEC_NEAR(oneNegative.zeroSegment[0], Eigen::Vector2d(0.5, 0.0), 1e-12);
        REQUIRE_VEC_NEAR(oneNegative.zeroSegment[1], Eigen::Vector2d(0.0, 0.5), 1e-12);

        auto twoNegative = MarchingTriangles::extract(a, b, c, -1.0, -1.0, 1.0);
        REQUIRE(twoNegative.hasContact);
        REQUIRE(twoNegative.negativePolygon.size() == 4);
        REQUIRE(twoNegative.hasZeroSegment);
        REQUIRE(std::abs(polygonArea2D(twoNegative.negativePolygon)) > 0.0);

        auto vertexZero = MarchingTriangles::extract(a, b, c, 0.0, 1.0, -1.0);
        REQUIRE(vertexZero.hasContact);
        REQUIRE(vertexZero.hasZeroSegment);
        REQUIRE(vertexZero.zeroSegment.size() == 2);
        REQUIRE(vertexZero.negativePolygon.size() >= 3);
    });
}

