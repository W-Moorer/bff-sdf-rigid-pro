#pragma once

#include "core/Mesh.h"

#include <Eigen/Dense>
#include <string>

namespace bff_sdf::io {

bool readObj(const std::string& path, bff_sdf::core::Mesh& mesh, std::string* error = nullptr);
bool writeObj(const std::string& path, const bff_sdf::core::Mesh& mesh, std::string* error = nullptr);
bool readStl(const std::string& path, bff_sdf::core::Mesh& mesh, std::string* error = nullptr);

bool readObjWithUV(const std::string& path,
                   bff_sdf::core::Mesh& mesh,
                   Eigen::MatrixXd& UV,
                   Eigen::MatrixXi& Fuv,
                   std::string* error = nullptr);

bool writeObjWithUV(const std::string& path,
                    const bff_sdf::core::Mesh& mesh,
                    const Eigen::MatrixXd& UV,
                    const Eigen::MatrixXi& Fuv,
                    std::string* error = nullptr);

} // namespace bff_sdf::io
