#include "atlas/BFFAdapter.h"

#include "io/MeshIO.h"

#include <Eigen/Eigenvalues>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace bff_sdf::atlas {
namespace fs = std::filesystem;

namespace {

std::string quotePath(const fs::path& p)
{
    return "\"" + p.string() + "\"";
}

std::string findBffExecutable(const BFFAdapterOptions& options)
{
    if (!options.executablePath.empty()) return options.executablePath;
    if (const char* env = std::getenv("BFF_COMMAND")) return std::string(env);
    const fs::path local = fs::current_path() / "bff_official" / "binaries" / "windows-v1.6" / "bff-command-line.exe";
    if (fs::exists(local)) return local.string();
    return {};
}

fs::path makeTempBase(const BFFAdapterOptions& options)
{
    fs::path dir = options.tempDirectory.empty() ? fs::temp_directory_path() : fs::path(options.tempDirectory);
    const auto ticks = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    return dir / ("bff_sdf_patch_" + std::to_string(ticks));
}

bool convertCornerUvsToVertexUvs(const bff_sdf::core::Mesh& mesh,
                                 const Eigen::MatrixXd& rawUV,
                                 const Eigen::MatrixXi& rawFuv,
                                 Eigen::MatrixXd& UV,
                                 Eigen::MatrixXi& Fuv)
{
    UV = Eigen::MatrixXd::Zero(mesh.V.rows(), 2);
    Eigen::VectorXi counts = Eigen::VectorXi::Zero(mesh.V.rows());
    if (rawUV.rows() == 0 || rawFuv.rows() != mesh.F.rows()) return false;
    for (int fi = 0; fi < mesh.F.rows(); ++fi) {
        for (int k = 0; k < 3; ++k) {
            const int vi = mesh.F(fi, k);
            const int ui = rawFuv(fi, k);
            if (ui < 0 || ui >= rawUV.rows()) return false;
            UV.row(vi) += rawUV.row(ui);
            counts[vi] += 1;
        }
    }
    for (int i = 0; i < counts.size(); ++i) {
        if (counts[i] == 0) return false;
        UV.row(i) /= static_cast<double>(counts[i]);
    }
    Fuv = mesh.F;
    return true;
}

} // namespace

int countFlippedTriangles(const Eigen::MatrixXd& UV, const Eigen::MatrixXi& Fuv)
{
    int flipped = 0;
    for (int i = 0; i < Fuv.rows(); ++i) {
        const Eigen::Vector2d a = UV.row(Fuv(i, 0));
        const Eigen::Vector2d b = UV.row(Fuv(i, 1));
        const Eigen::Vector2d c = UV.row(Fuv(i, 2));
        if (signedUVArea(a, b, c) <= 1e-14) ++flipped;
    }
    return flipped;
}

void orientUvToMatchMajority(Eigen::MatrixXd& UV, const Eigen::MatrixXi& Fuv)
{
    double sum = 0.0;
    for (int i = 0; i < Fuv.rows(); ++i) {
        sum += signedUVArea(UV.row(Fuv(i, 0)), UV.row(Fuv(i, 1)), UV.row(Fuv(i, 2)));
    }
    if (sum < 0.0) UV.col(0) *= -1.0;
}

BFFAdapterResult BFFAdapter::parameterize(const Patch& patch, const BFFAdapterOptions& options) const
{
    if (!patch.diskLike) {
        return fallbackPlanar(patch, "patch is not disk-like: " + patch.diagnostic);
    }
    if (options.preferOfficialBff) {
        BFFAdapterResult official = runOfficialBff(patch, options);
        if (official.success) return official;
        return fallbackPlanar(patch, "official BFF failed: " + official.message);
    }
    return fallbackPlanar(patch, "official BFF disabled by options");
}

ContactChart BFFAdapter::buildChart(int chartId, const Patch& patch, const BFFAdapterOptions& options) const
{
    const BFFAdapterResult result = parameterize(patch, options);
    return makeContactChart(chartId,
                            patch,
                            result.UV,
                            result.Fuv,
                            result.usedOfficialBff ? "official-bff" : "fallback-planar",
                            result.usedFallback);
}

BFFAdapterResult BFFAdapter::runOfficialBff(const Patch& patch, const BFFAdapterOptions& options) const
{
    BFFAdapterResult result;
    const std::string exe = findBffExecutable(options);
    if (exe.empty()) {
        result.message = "BFF executable not found";
        return result;
    }

    const fs::path base = makeTempBase(options);
    const fs::path inPath = base.string() + "_in.obj";
    const fs::path outPath = base.string() + "_out.obj";
    std::string error;
    if (!bff_sdf::io::writeObj(inPath.string(), patch.localMesh, &error)) {
        result.message = error;
        return result;
    }

    std::ostringstream cmd;
    cmd << quotePath(exe) << " " << quotePath(inPath) << " " << quotePath(outPath);
    if (options.flattenToDisk) cmd << " --flattenToDisk";
    if (options.normalizeUVs) cmd << " --normalizeUVs";

    const int code = std::system(cmd.str().c_str());
    if (code != 0 || !fs::exists(outPath)) {
        result.message = "BFF command failed with code " + std::to_string(code);
        fs::remove(inPath);
        fs::remove(outPath);
        return result;
    }

    bff_sdf::core::Mesh outMesh;
    Eigen::MatrixXd rawUV;
    Eigen::MatrixXi rawFuv;
    if (!bff_sdf::io::readObjWithUV(outPath.string(), outMesh, rawUV, rawFuv, &error)) {
        result.message = error;
        fs::remove(inPath);
        fs::remove(outPath);
        return result;
    }

    if (!convertCornerUvsToVertexUvs(patch.localMesh, rawUV, rawFuv, result.UV, result.Fuv)) {
        result.message = "unable to convert BFF corner UVs to per-vertex UVs";
        fs::remove(inPath);
        fs::remove(outPath);
        return result;
    }

    orientUvToMatchMajority(result.UV, result.Fuv);
    result.flippedTriangleCount = countFlippedTriangles(result.UV, result.Fuv);
    result.success = result.UV.array().isFinite().all() && result.flippedTriangleCount == 0;
    result.usedOfficialBff = result.success;
    result.message = result.success ? "official BFF parameterization" : "official BFF produced flipped UV triangles";

    fs::remove(inPath);
    fs::remove(outPath);
    return result;
}

BFFAdapterResult BFFAdapter::fallbackPlanar(const Patch& patch, const std::string& reason) const
{
    BFFAdapterResult result;
    result.usedFallback = true;
    result.message = reason;
    result.Fuv = patch.localMesh.F;
    result.UV.resize(patch.localMesh.V.rows(), 2);

    Eigen::Vector3d normal = Eigen::Vector3d::Zero();
    for (int i = 0; i < patch.localMesh.FN.rows(); ++i) normal += patch.localMesh.FN.row(i).transpose();
    if (normal.norm() < 1e-12) normal = Eigen::Vector3d::UnitZ();
    normal.normalize();

    Eigen::Vector3d axis0 = normal.cross(Eigen::Vector3d::UnitX());
    if (axis0.norm() < 1e-8) axis0 = normal.cross(Eigen::Vector3d::UnitY());
    axis0.normalize();
    Eigen::Vector3d axis1 = normal.cross(axis0).normalized();

    for (int i = 0; i < patch.localMesh.V.rows(); ++i) {
        const Eigen::Vector3d p = patch.localMesh.V.row(i);
        result.UV(i, 0) = p.dot(axis0);
        result.UV(i, 1) = p.dot(axis1);
    }

    orientUvToMatchMajority(result.UV, result.Fuv);
    result.flippedTriangleCount = countFlippedTriangles(result.UV, result.Fuv);
    result.success = result.UV.array().isFinite().all() && result.flippedTriangleCount == 0;
    return result;
}

} // namespace bff_sdf::atlas

