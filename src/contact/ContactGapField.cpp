#include "contact/ContactGapField.h"

namespace bff_sdf::contact {

GapSample ContactGapField::sample(const bff_sdf::atlas::ContactChart& sourceChart,
                                  int sourceTriId,
                                  const Eigen::Vector3d& bary,
                                  const bff_sdf::core::RigidPose& sourcePose,
                                  const bff_sdf::sdf::ISDF& targetSdf,
                                  const bff_sdf::core::RigidPose& targetPose)
{
    GapSample out;
    out.sourceChartId = sourceChart.id;
    out.sourceTriId = sourceTriId;
    out.bary = bary;

    const auto ev = bff_sdf::atlas::ChartEvaluator::evaluate(sourceChart, sourceTriId, bary, sourcePose);
    out.xWorld = ev.xWorld;
    const Eigen::Vector3d yLocal = targetPose.worldToLocal(ev.xWorld);
    const auto q = targetSdf.evalLocal(yLocal, true, true);
    out.g = q.phi;
    Eigen::Vector3d gradLocal = q.grad;
    if (gradLocal.norm() < 1e-12 || !gradLocal.array().isFinite().all()) {
        gradLocal = targetPose.rotateWorldVectorToLocal(ev.normalWorld);
    }
    out.normalWorld = targetPose.rotateLocalVectorToWorld(gradLocal).normalized();
    out.gradU = ev.JWorld.transpose() * out.normalWorld;
    return out;
}

} // namespace bff_sdf::contact

