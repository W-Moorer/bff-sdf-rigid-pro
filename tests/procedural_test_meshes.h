#pragma once

#include "core/Mesh.h"

#include <cmath>

namespace test_meshes {

constexpr double pi = 3.141592653589793238462643383279502884;

inline bff_sdf::core::Mesh makeCube()
{
    bff_sdf::core::Mesh m;
    m.V.resize(8, 3);
    m.V << -1, -1, -1,
            1, -1, -1,
            1,  1, -1,
           -1,  1, -1,
           -1, -1,  1,
            1, -1,  1,
            1,  1,  1,
           -1,  1,  1;
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

inline bff_sdf::core::Mesh makeSphere(double radius, int rings, int segments)
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

    const int faces = 2 * segments + (ringCount - 1) * 2 * segments;
    m.F.resize(faces, 3);
    row = 0;
    auto vid = [segments](int ringZeroBased, int s) {
        return 1 + ringZeroBased * segments + ((s + segments) % segments);
    };
    for (int s = 0; s < segments; ++s) {
        m.F.row(row++) = Eigen::Vector3i(0, vid(0, s), vid(0, s + 1));
    }
    for (int r = 0; r < ringCount - 1; ++r) {
        for (int s = 0; s < segments; ++s) {
            const int a = vid(r, s);
            const int b = vid(r, s + 1);
            const int c = vid(r + 1, s);
            const int d = vid(r + 1, s + 1);
            m.F.row(row++) = Eigen::Vector3i(a, c, b);
            m.F.row(row++) = Eigen::Vector3i(b, c, d);
        }
    }
    for (int s = 0; s < segments; ++s) {
        m.F.row(row++) = Eigen::Vector3i(bottom, vid(ringCount - 1, s + 1), vid(ringCount - 1, s));
    }
    m.computeNormals();
    return m;
}

inline bff_sdf::core::Mesh makeLowerHemispherePatch(double radius, int rings, int segments)
{
    bff_sdf::core::Mesh m;
    const int vertices = 1 + rings * segments;
    m.V.resize(vertices, 3);
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
    for (int s = 0; s < segments; ++s) {
        m.F.row(row++) = Eigen::Vector3i(0, vid(1, s + 1), vid(1, s));
    }
    for (int r = 1; r < rings; ++r) {
        for (int s = 0; s < segments; ++s) {
            const int a = vid(r, s);
            const int b = vid(r, s + 1);
            const int c = vid(r + 1, s);
            const int d = vid(r + 1, s + 1);
            m.F.row(row++) = Eigen::Vector3i(a, b, c);
            m.F.row(row++) = Eigen::Vector3i(b, d, c);
        }
    }
    m.computeNormals();
    return m;
}

} // namespace test_meshes
