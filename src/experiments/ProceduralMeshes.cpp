#include "experiments/ProceduralMeshes.h"

#include <cmath>

namespace bff_sdf::experiments {

bff_sdf::core::Mesh makePlane(double halfSize)
{
    bff_sdf::core::Mesh m;
    m.V.resize(4, 3);
    m.V << -halfSize, -halfSize, 0,
            halfSize, -halfSize, 0,
            halfSize,  halfSize, 0,
           -halfSize,  halfSize, 0;
    m.F.resize(2, 3);
    m.F << 0, 1, 2,
           0, 2, 3;
    m.computeNormals();
    return m;
}

bff_sdf::core::Mesh makeCube(double halfSize)
{
    bff_sdf::core::Mesh m;
    m.V.resize(8, 3);
    const double h = halfSize;
    m.V << -h, -h, -h, h, -h, -h, h, h, -h, -h, h, -h,
           -h, -h, h, h, -h, h, h, h, h, -h, h, h;
    m.F.resize(12, 3);
    m.F << 0, 2, 1, 0, 3, 2,
           4, 5, 6, 4, 6, 7,
           0, 1, 5, 0, 5, 4,
           1, 2, 6, 1, 6, 5,
           2, 3, 7, 2, 7, 6,
           3, 0, 4, 3, 4, 7;
    m.computeNormals();
    return m;
}

bff_sdf::core::Mesh makeSphere(double radius, int rings, int segments)
{
    bff_sdf::core::Mesh m;
    const int ringCount = rings - 1;
    const int vertices = 2 + ringCount * segments;
    m.V.resize(vertices, 3);
    m.V.row(0) = Eigen::Vector3d(0, 0, radius);
    int row = 1;
    for (int r = 1; r < rings; ++r) {
        const double theta = pi * static_cast<double>(r) / static_cast<double>(rings);
        const double z = radius * std::cos(theta);
        const double xy = radius * std::sin(theta);
        for (int s = 0; s < segments; ++s) {
            const double phi = 2.0 * pi * static_cast<double>(s) / static_cast<double>(segments);
            m.V.row(row++) = Eigen::Vector3d(xy * std::cos(phi), xy * std::sin(phi), z);
        }
    }
    const int bottom = vertices - 1;
    m.V.row(bottom) = Eigen::Vector3d(0, 0, -radius);
    m.F.resize(2 * segments + (ringCount - 1) * 2 * segments, 3);
    row = 0;
    auto vid = [segments](int ringZeroBased, int s) {
        return 1 + ringZeroBased * segments + ((s + segments) % segments);
    };
    for (int s = 0; s < segments; ++s) m.F.row(row++) = Eigen::Vector3i(0, vid(0, s), vid(0, s + 1));
    for (int r = 0; r < ringCount - 1; ++r) {
        for (int s = 0; s < segments; ++s) {
            const int a = vid(r, s), b = vid(r, s + 1), c = vid(r + 1, s), d = vid(r + 1, s + 1);
            m.F.row(row++) = Eigen::Vector3i(a, c, b);
            m.F.row(row++) = Eigen::Vector3i(b, c, d);
        }
    }
    for (int s = 0; s < segments; ++s) m.F.row(row++) = Eigen::Vector3i(bottom, vid(ringCount - 1, s + 1), vid(ringCount - 1, s));
    m.computeNormals();
    return m;
}

bff_sdf::core::Mesh makeLowerHemispherePatch(double radius, int rings, int segments)
{
    bff_sdf::core::Mesh m;
    m.V.resize(1 + rings * segments, 3);
    m.V.row(0) = Eigen::Vector3d(0, 0, -radius);
    int row = 1;
    for (int r = 1; r <= rings; ++r) {
        const double theta = pi - 0.5 * pi * static_cast<double>(r) / static_cast<double>(rings);
        const double z = radius * std::cos(theta);
        const double xy = radius * std::sin(theta);
        for (int s = 0; s < segments; ++s) {
            const double phi = 2.0 * pi * static_cast<double>(s) / static_cast<double>(segments);
            m.V.row(row++) = Eigen::Vector3d(xy * std::cos(phi), xy * std::sin(phi), z);
        }
    }
    m.F.resize(segments + (rings - 1) * 2 * segments, 3);
    row = 0;
    auto vid = [segments](int ringOneBased, int s) {
        return 1 + (ringOneBased - 1) * segments + ((s + segments) % segments);
    };
    for (int s = 0; s < segments; ++s) m.F.row(row++) = Eigen::Vector3i(0, vid(1, s + 1), vid(1, s));
    for (int r = 1; r < rings; ++r) {
        for (int s = 0; s < segments; ++s) {
            const int a = vid(r, s), b = vid(r, s + 1), c = vid(r + 1, s), d = vid(r + 1, s + 1);
            m.F.row(row++) = Eigen::Vector3i(a, b, c);
            m.F.row(row++) = Eigen::Vector3i(b, d, c);
        }
    }
    m.computeNormals();
    return m;
}

bff_sdf::core::Mesh makeEllipsoidLowerPatch(const Eigen::Vector3d& radii, int rings, int segments)
{
    auto m = makeLowerHemispherePatch(1.0, rings, segments);
    for (int i = 0; i < m.V.rows(); ++i) {
        m.V(i, 0) *= radii.x();
        m.V(i, 1) *= radii.y();
        m.V(i, 2) *= radii.z();
    }
    m.computeNormals();
    return m;
}

bff_sdf::core::Mesh makeTorusPatch(double majorRadius, double minorRadius, int majorSegments, int minorSegments)
{
    bff_sdf::core::Mesh m;
    m.V.resize((minorSegments + 1) * majorSegments, 3);
    int row = 0;
    for (int i = 0; i <= minorSegments; ++i) {
        const double v = pi + pi * static_cast<double>(i) / static_cast<double>(minorSegments);
        for (int j = 0; j < majorSegments; ++j) {
            const double u = 2.0 * pi * static_cast<double>(j) / static_cast<double>(majorSegments);
            const double r = majorRadius + minorRadius * std::cos(v);
            m.V.row(row++) = Eigen::Vector3d(r * std::cos(u), r * std::sin(u), minorRadius * std::sin(v));
        }
    }
    m.F.resize(minorSegments * 2 * majorSegments, 3);
    row = 0;
    auto vid = [majorSegments](int i, int j) { return i * majorSegments + ((j + majorSegments) % majorSegments); };
    for (int i = 0; i < minorSegments; ++i) {
        for (int j = 0; j < majorSegments; ++j) {
            const int a = vid(i, j), b = vid(i, j + 1), c = vid(i + 1, j), d = vid(i + 1, j + 1);
            m.F.row(row++) = Eigen::Vector3i(a, b, c);
            m.F.row(row++) = Eigen::Vector3i(b, d, c);
        }
    }
    m.computeNormals();
    return m;
}

Eigen::MatrixXd radialUV(const bff_sdf::core::Mesh& mesh)
{
    Eigen::MatrixXd UV(mesh.V.rows(), 2);
    for (int i = 0; i < mesh.V.rows(); ++i) {
        UV(i, 0) = mesh.V(i, 0);
        UV(i, 1) = mesh.V(i, 1);
    }
    return UV;
}

} // namespace bff_sdf::experiments

