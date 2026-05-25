#pragma once

#include <Eigen/Dense>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace bff_sdf::core {

using Index = int;
using Face = Eigen::Vector3i;
using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;

struct OperationStatus {
    bool ok = true;
    std::string message;

    static OperationStatus success() { return {}; }
    static OperationStatus failure(std::string msg) { return {false, std::move(msg)}; }
};

} // namespace bff_sdf::core
