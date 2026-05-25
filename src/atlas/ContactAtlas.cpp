#include "atlas/ContactAtlas.h"

namespace bff_sdf::atlas {

ContactChart makeContactChart(int id,
                              const Patch& patch,
                              const Eigen::MatrixXd& UV,
                              const Eigen::MatrixXi& Fuv,
                              const std::string& method,
                              bool usedFallback)
{
    ContactChart chart;
    chart.id = id;
    chart.patch = patch;
    chart.UV = UV;
    chart.Fuv = Fuv;
    chart.localToGlobalVertex = patch.globalVertexIds;
    chart.localFaceToGlobalFace = patch.globalFaceIds;
    chart.parameterizationMethod = method;
    chart.usedFallbackParameterization = usedFallback;
    chart.uvAccel.build(chart.UV, chart.Fuv);
    return chart;
}

} // namespace bff_sdf::atlas

