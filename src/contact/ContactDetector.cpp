#include "contact/ContactDetector.h"

#include "atlas/ChartEvaluator.h"
#include "contact/ChartProjector.h"
#include "contact/ContactGapField.h"
#include "contact/ContactIntegrator.h"

#include <chrono>
#include <numeric>

namespace bff_sdf::contact {

namespace {

double elapsedMs(std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b)
{
    return std::chrono::duration<double, std::milli>(b - a).count();
}

void appendResult(ContactResult& dst, const ContactResult& src)
{
    const double oldArea = dst.area;
    dst.area += src.area;
    dst.maxPenetration = std::max(dst.maxPenetration, src.maxPenetration);
    if (dst.area > 0.0) {
        dst.centroid = (oldArea * dst.centroid + src.area * src.centroid) / dst.area;
        dst.meanNormal = oldArea * dst.meanNormal + src.area * src.meanNormal;
        if (dst.meanNormal.norm() > 1e-12) dst.meanNormal.normalize();
        dst.meanPenetration = (oldArea * dst.meanPenetration + src.area * src.meanPenetration) / dst.area;
    }
    dst.sdfQueries += src.sdfQueries;
    dst.contours.insert(dst.contours.end(), src.contours.begin(), src.contours.end());
}

} // namespace

ContactResult ContactDetector::detect(const RigidObject& source,
                                      const RigidObject& target,
                                      const ContactDetectorOptions& options)
{
    ContactResult result;
    const auto t0 = std::chrono::steady_clock::now();
    if (!source.atlas || !target.sdf) return result;

    ChartProjector projector;
    for (const auto& chart : source.atlas->charts) {
        const auto ti0 = std::chrono::steady_clock::now();
        auto integrated = ContactIntegrator::integrateChart(chart, source.pose, *target.sdf, target.pose);
        const auto ti1 = std::chrono::steady_clock::now();
        result.timing.integrationMs += elapsedMs(ti0, ti1);
        appendResult(result, integrated);

        const int stride = std::max(1, options.sampleStride);
        for (int tri = 0; tri < chart.patch.localMesh.F.rows(); tri += stride) {
            const Eigen::Vector3d bary(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
            const auto gs = ContactGapField::sample(chart, tri, bary, source.pose, *target.sdf, target.pose);
            result.sdfQueries += 1;
            if (gs.g > options.contactMargin) continue;

            ContactSample sample;
            sample.xA = gs.xWorld;
            sample.gap = gs.g;
            sample.penetration = std::max(0.0, -gs.g);
            sample.normal = gs.normalWorld;
            sample.sourceChartId = chart.id;
            sample.sourceTriId = tri;
            sample.sourceBary = bary;

            const auto qp = target.sdf->evalLocal(target.pose.worldToLocal(gs.xWorld), true, true);
            result.sdfQueries += 1;
            if (options.refineProjection && target.mesh) {
                const auto tp0 = std::chrono::steady_clock::now();
                const auto pr = projector.project(gs.xWorld, *target.mesh, target.pose, &qp);
                const auto tp1 = std::chrono::steady_clock::now();
                result.timing.projectionMs += elapsedMs(tp0, tp1);
                result.projectionCount += 1;
                if (pr.converged) {
                    sample.xB = pr.closestPointWorld;
                    sample.normal = pr.normalWorld;
                    sample.gap = pr.signedGap;
                    sample.penetration = std::max(0.0, -pr.signedGap);
                    sample.targetTriId = pr.triId;
                    sample.targetBary = pr.bary;
                } else {
                    result.projectionFailures += 1;
                    sample.xB = target.pose.localToWorld(qp.closestPoint);
                }
            } else {
                sample.xB = qp.hasClosest ? target.pose.localToWorld(qp.closestPoint) : gs.xWorld - gs.g * gs.normalWorld;
            }
            result.samples.push_back(sample);
        }
    }
    const auto t1 = std::chrono::steady_clock::now();
    result.timing.totalMs = elapsedMs(t0, t1);
    return result;
}

} // namespace bff_sdf::contact

