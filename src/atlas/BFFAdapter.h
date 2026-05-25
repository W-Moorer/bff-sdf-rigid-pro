#pragma once

#include "atlas/ContactAtlas.h"

#include <string>

namespace bff_sdf::atlas {

struct BFFAdapterOptions {
    bool preferOfficialBff = true;
    bool flattenToDisk = true;
    bool normalizeUVs = false;
    std::string executablePath;
    std::string tempDirectory;
};

struct BFFAdapterResult {
    bool success = false;
    bool usedOfficialBff = false;
    bool usedFallback = false;
    int flippedTriangleCount = 0;
    std::string message;
    Eigen::MatrixXd UV;
    Eigen::MatrixXi Fuv;
};

class BFFAdapter {
public:
    BFFAdapterResult parameterize(const Patch& patch, const BFFAdapterOptions& options = {}) const;
    ContactChart buildChart(int chartId, const Patch& patch, const BFFAdapterOptions& options = {}) const;

private:
    BFFAdapterResult runOfficialBff(const Patch& patch, const BFFAdapterOptions& options) const;
    BFFAdapterResult fallbackPlanar(const Patch& patch, const std::string& reason) const;
};

int countFlippedTriangles(const Eigen::MatrixXd& UV, const Eigen::MatrixXi& Fuv);
void orientUvToMatchMajority(Eigen::MatrixXd& UV, const Eigen::MatrixXi& Fuv);

} // namespace bff_sdf::atlas

