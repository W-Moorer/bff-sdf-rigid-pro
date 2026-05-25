#pragma once

#include "core/RigidPose.h"

#include <Eigen/Dense>
#include <vector>

namespace bff_sdf::contact {

struct TimingBreakdown {
    double broadPhaseMs = 0.0;
    double gapFieldMs = 0.0;
    double marchingMs = 0.0;
    double integrationMs = 0.0;
    double projectionMs = 0.0;
    double totalMs = 0.0;
};

struct ContactContour {
    std::vector<Eigen::Vector2d> uvPoints;
    std::vector<Eigen::Vector3d> xWorldPoints;
};

struct ContactSample {
    Eigen::Vector3d xA = Eigen::Vector3d::Zero();
    Eigen::Vector3d xB = Eigen::Vector3d::Zero();
    Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();
    double gap = 0.0;
    double penetration = 0.0;
    int sourceChartId = -1;
    int sourceTriId = -1;
    Eigen::Vector3d sourceBary = Eigen::Vector3d::Zero();
    int targetChartId = -1;
    int targetTriId = -1;
    Eigen::Vector3d targetBary = Eigen::Vector3d::Zero();
};

struct ContactResult {
    std::vector<ContactSample> samples;
    std::vector<ContactContour> contours;
    double area = 0.0;
    double maxPenetration = 0.0;
    double meanPenetration = 0.0;
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    Eigen::Vector3d meanNormal = Eigen::Vector3d::Zero();
    TimingBreakdown timing;
    int sdfQueries = 0;
    int projectionCount = 0;
    int projectionFailures = 0;
    int fallbackToBvh = 0;
    int seamCrossingCount = 0;
    int badSdfGradientCount = 0;
};

struct GapSample {
    double g = 0.0;
    Eigen::Vector3d xWorld = Eigen::Vector3d::Zero();
    Eigen::Vector3d normalWorld = Eigen::Vector3d::UnitZ();
    Eigen::Vector2d gradU = Eigen::Vector2d::Zero();
    int sourceChartId = -1;
    int sourceTriId = -1;
    Eigen::Vector3d bary = Eigen::Vector3d::Zero();
};

} // namespace bff_sdf::contact

