#include "contact/ContactIntegrator.h"

#include "atlas/ChartEvaluator.h"
#include "atlas/UVAcceleration.h"
#include "contact/ContactGapField.h"
#include "contact/MarchingTriangles.h"

#include <algorithm>
#include <cmath>

namespace bff_sdf::contact {
namespace {

Eigen::Vector3d mapUvToWorld(const bff_sdf::atlas::ContactChart& chart,
                             int triId,
                             const Eigen::Vector2d& u,
                             const bff_sdf::core::RigidPose& pose)
{
    const Eigen::Vector3i fu = chart.Fuv.row(triId);
    Eigen::Vector3d bary;
    const bool ok = bff_sdf::atlas::barycentricUV(u,
                                                  chart.UV.row(fu[0]),
                                                  chart.UV.row(fu[1]),
                                                  chart.UV.row(fu[2]),
                                                  bary);
    if (!ok) return Eigen::Vector3d::Zero();
    return bff_sdf::atlas::ChartEvaluator::evaluate(chart, triId, bary, pose).xWorld;
}

double triangleArea3D(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
{
    return 0.5 * (b - a).cross(c - a).norm();
}

} // namespace

ContactResult ContactIntegrator::integrateChart(const bff_sdf::atlas::ContactChart& sourceChart,
                                                const bff_sdf::core::RigidPose& sourcePose,
                                                const bff_sdf::sdf::ISDF& targetSdf,
                                                const bff_sdf::core::RigidPose& targetPose,
                                                const ContactIntegrationOptions& options)
{
    ContactResult result;
    Eigen::Vector3d weightedCentroid = Eigen::Vector3d::Zero();
    Eigen::Vector3d weightedNormal = Eigen::Vector3d::Zero();
    double penetrationIntegral = 0.0;

    for (int tri = 0; tri < sourceChart.patch.localMesh.F.rows(); ++tri) {
        const Eigen::Vector3i fu = sourceChart.Fuv.row(tri);
        const Eigen::Vector2d u0 = sourceChart.UV.row(fu[0]);
        const Eigen::Vector2d u1 = sourceChart.UV.row(fu[1]);
        const Eigen::Vector2d u2 = sourceChart.UV.row(fu[2]);

        const auto g0 = ContactGapField::sample(sourceChart, tri, Eigen::Vector3d(1, 0, 0), sourcePose, targetSdf, targetPose);
        const auto g1 = ContactGapField::sample(sourceChart, tri, Eigen::Vector3d(0, 1, 0), sourcePose, targetSdf, targetPose);
        const auto g2 = ContactGapField::sample(sourceChart, tri, Eigen::Vector3d(0, 0, 1), sourcePose, targetSdf, targetPose);
        result.sdfQueries += 3;

        const auto mt = MarchingTriangles::extract(u0, u1, u2, g0.g, g1.g, g2.g, options.zeroEps);
        if (mt.hasZeroSegment) {
            ContactContour contour;
            contour.uvPoints = mt.zeroSegment;
            for (const auto& u : mt.zeroSegment) contour.xWorldPoints.push_back(mapUvToWorld(sourceChart, tri, u, sourcePose));
            result.contours.push_back(contour);
        }
        if (!mt.hasContact) continue;

        std::vector<Eigen::Vector3d> poly3;
        poly3.reserve(mt.negativePolygon.size());
        for (const auto& u : mt.negativePolygon) poly3.push_back(mapUvToWorld(sourceChart, tri, u, sourcePose));

        const Eigen::Vector3d base = poly3.front();
        for (size_t i = 1; i + 1 < poly3.size(); ++i) {
            const double area = triangleArea3D(base, poly3[i], poly3[i + 1]);
            if (area <= 0.0) continue;
            const Eigen::Vector3d centroid = (base + poly3[i] + poly3[i + 1]) / 3.0;
            Eigen::Vector3d bary;
            const Eigen::Vector2d uvCentroid =
                (mt.negativePolygon[0] + mt.negativePolygon[i] + mt.negativePolygon[i + 1]) / 3.0;
            bff_sdf::atlas::barycentricUV(uvCentroid, u0, u1, u2, bary);
            const auto gs = ContactGapField::sample(sourceChart, tri, bary, sourcePose, targetSdf, targetPose);
            result.sdfQueries += 1;
            const double penetration = std::max(0.0, -gs.g);

            result.area += area;
            weightedCentroid += area * centroid;
            weightedNormal += area * gs.normalWorld;
            penetrationIntegral += area * penetration;
            result.maxPenetration = std::max(result.maxPenetration, penetration);
        }
    }

    if (result.area > 0.0) {
        result.centroid = weightedCentroid / result.area;
        result.meanPenetration = penetrationIntegral / result.area;
        if (weightedNormal.norm() > 1e-12) result.meanNormal = weightedNormal.normalized();
    }
    return result;
}

} // namespace bff_sdf::contact

