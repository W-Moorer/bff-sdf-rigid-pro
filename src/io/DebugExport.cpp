#include "io/DebugExport.h"

#include "atlas/ChartEvaluator.h"
#include "contact/ContactGapField.h"

#include <filesystem>
#include <fstream>

namespace bff_sdf::io {

bool ensureDirectory(const std::string& dir)
{
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    return !ec;
}

bool writeGapFieldCsv(const std::string& path,
                      const bff_sdf::atlas::ContactChart& chart,
                      const bff_sdf::core::RigidPose& sourcePose,
                      const bff_sdf::sdf::ISDF& targetSdf,
                      const bff_sdf::core::RigidPose& targetPose)
{
    std::ofstream out(path);
    if (!out) return false;
    out << "chart_id,tri_id,vertex_id,u,v,g,x,y,z,nx,ny,nz,grad_u,grad_v\n";
    for (int tri = 0; tri < chart.patch.localMesh.F.rows(); ++tri) {
        for (int k = 0; k < 3; ++k) {
            Eigen::Vector3d bary = Eigen::Vector3d::Zero();
            bary[k] = 1.0;
            const auto gs = bff_sdf::contact::ContactGapField::sample(chart, tri, bary, sourcePose, targetSdf, targetPose);
            const auto ev = bff_sdf::atlas::ChartEvaluator::evaluate(chart, tri, bary, sourcePose);
            out << chart.id << "," << tri << "," << k << ","
                << ev.u.x() << "," << ev.u.y() << ","
                << gs.g << ","
                << gs.xWorld.x() << "," << gs.xWorld.y() << "," << gs.xWorld.z() << ","
                << gs.normalWorld.x() << "," << gs.normalWorld.y() << "," << gs.normalWorld.z() << ","
                << gs.gradU.x() << "," << gs.gradU.y() << "\n";
        }
    }
    return true;
}

bool writeContactContourObj(const std::string& path, const bff_sdf::contact::ContactResult& result)
{
    std::ofstream out(path);
    if (!out) return false;
    int index = 1;
    for (const auto& contour : result.contours) {
        for (const auto& p : contour.xWorldPoints) out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
        if (contour.xWorldPoints.size() == 2) out << "l " << index << " " << index + 1 << "\n";
        index += static_cast<int>(contour.xWorldPoints.size());
    }
    return true;
}

bool writeContactContourUvCsv(const std::string& path, const bff_sdf::contact::ContactResult& result)
{
    std::ofstream out(path);
    if (!out) return false;
    out << "contour_id,point_id,u,v\n";
    for (int ci = 0; ci < static_cast<int>(result.contours.size()); ++ci) {
        const auto& contour = result.contours[ci];
        for (int pi = 0; pi < static_cast<int>(contour.uvPoints.size()); ++pi) {
            out << ci << "," << pi << "," << contour.uvPoints[pi].x() << "," << contour.uvPoints[pi].y() << "\n";
        }
    }
    return true;
}

bool writeContactSamplesObj(const std::string& path, const bff_sdf::contact::ContactResult& result)
{
    std::ofstream out(path);
    if (!out) return false;
    for (const auto& s : result.samples) out << "v " << s.xA.x() << " " << s.xA.y() << " " << s.xA.z() << "\n";
    return true;
}

bool writeContactSamplesPly(const std::string& path, const bff_sdf::contact::ContactResult& result)
{
    std::ofstream out(path);
    if (!out) return false;
    out << "ply\nformat ascii 1.0\nelement vertex " << result.samples.size()
        << "\nproperty float x\nproperty float y\nproperty float z\nproperty float gap\nend_header\n";
    for (const auto& s : result.samples) out << s.xA.x() << " " << s.xA.y() << " " << s.xA.z() << " " << s.gap << "\n";
    return true;
}

bool writeAtlasDebugObj(const std::string& path,
                        const bff_sdf::atlas::ContactChart& chart,
                        const bff_sdf::core::RigidPose& pose)
{
    std::ofstream out(path);
    if (!out) return false;
    for (int i = 0; i < chart.patch.localMesh.V.rows(); ++i) {
        const Eigen::Vector3d p = pose.localToWorld(chart.patch.localMesh.V.row(i));
        out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
    for (int i = 0; i < chart.patch.localMesh.F.rows(); ++i) {
        out << "f " << chart.patch.localMesh.F(i, 0) + 1 << " "
            << chart.patch.localMesh.F(i, 1) + 1 << " "
            << chart.patch.localMesh.F(i, 2) + 1 << "\n";
    }
    return true;
}

} // namespace bff_sdf::io

