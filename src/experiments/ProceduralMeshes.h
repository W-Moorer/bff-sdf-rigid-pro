#pragma once

#include "core/Mesh.h"

namespace bff_sdf::experiments {

constexpr double pi = 3.141592653589793238462643383279502884;

bff_sdf::core::Mesh makePlane(double halfSize);
bff_sdf::core::Mesh makeCube(double halfSize);
bff_sdf::core::Mesh makeSphere(double radius, int rings, int segments);
bff_sdf::core::Mesh makeLowerHemispherePatch(double radius, int rings, int segments);
bff_sdf::core::Mesh makeEllipsoidLowerPatch(const Eigen::Vector3d& radii, int rings, int segments);
bff_sdf::core::Mesh makeTorusPatch(double majorRadius, double minorRadius, int majorSegments, int minorSegments);

Eigen::MatrixXd radialUV(const bff_sdf::core::Mesh& mesh);

} // namespace bff_sdf::experiments

