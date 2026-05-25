#pragma once

#include "atlas/Patch.h"
#include "atlas/UVAcceleration.h"

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace bff_sdf::atlas {

struct ContactChart {
    int id = -1;
    Patch patch;
    Eigen::MatrixXd UV;
    Eigen::MatrixXi Fuv;
    UVAcceleration uvAccel;
    std::vector<int> localToGlobalVertex;
    std::vector<int> localFaceToGlobalFace;
    std::vector<int> adjacentCharts;
    std::string parameterizationMethod;
    bool usedFallbackParameterization = false;
};

struct ContactAtlas {
    std::vector<ContactChart> charts;
};

ContactChart makeContactChart(int id,
                              const Patch& patch,
                              const Eigen::MatrixXd& UV,
                              const Eigen::MatrixXi& Fuv,
                              const std::string& method,
                              bool usedFallback);

} // namespace bff_sdf::atlas

