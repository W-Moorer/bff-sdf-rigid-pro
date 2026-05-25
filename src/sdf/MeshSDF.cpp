#include "sdf/MeshSDF.h"

#include <cmath>
#include <limits>

namespace bff_sdf::sdf {

TriangleClosestPoint closestPointOnTriangle(const Eigen::Vector3d& p,
                                            const Eigen::Vector3d& a,
                                            const Eigen::Vector3d& b,
                                            const Eigen::Vector3d& c)
{
    const Eigen::Vector3d ab = b - a;
    const Eigen::Vector3d ac = c - a;
    const Eigen::Vector3d ap = p - a;
    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0) return {a, Eigen::Vector3d(1, 0, 0), (p - a).squaredNorm()};

    const Eigen::Vector3d bp = p - b;
    const double d3 = ab.dot(bp);
    const double d4 = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3) return {b, Eigen::Vector3d(0, 1, 0), (p - b).squaredNorm()};

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        const Eigen::Vector3d q = a + v * ab;
        return {q, Eigen::Vector3d(1.0 - v, v, 0.0), (p - q).squaredNorm()};
    }

    const Eigen::Vector3d cp = p - c;
    const double d5 = ab.dot(cp);
    const double d6 = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6) return {c, Eigen::Vector3d(0, 0, 1), (p - c).squaredNorm()};

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        const Eigen::Vector3d q = a + w * ac;
        return {q, Eigen::Vector3d(1.0 - w, 0.0, w), (p - q).squaredNorm()};
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        const Eigen::Vector3d q = b + w * (c - b);
        return {q, Eigen::Vector3d(0.0, 1.0 - w, w), (p - q).squaredNorm()};
    }

    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom;
    const double w = vc * denom;
    const double u = 1.0 - v - w;
    const Eigen::Vector3d q = u * a + v * b + w * c;
    return {q, Eigen::Vector3d(u, v, w), (p - q).squaredNorm()};
}

MeshSDF::MeshSDF(const bff_sdf::core::Mesh& mesh, SignMode signMode)
    : mesh_(mesh), signMode_(signMode)
{
    mesh_.computeNormals();
}

bool MeshSDF::rayIntersectsTriangle(const Eigen::Vector3d& origin,
                                    const Eigen::Vector3d& a,
                                    const Eigen::Vector3d& b,
                                    const Eigen::Vector3d& c) const
{
    constexpr double eps = 1e-12;
    const Eigen::Vector3d dir = Eigen::Vector3d::UnitX();
    const Eigen::Vector3d e1 = b - a;
    const Eigen::Vector3d e2 = c - a;
    const Eigen::Vector3d h = dir.cross(e2);
    const double det = e1.dot(h);
    if (std::abs(det) < eps) return false;
    const double invDet = 1.0 / det;
    const Eigen::Vector3d s = origin - a;
    const double u = invDet * s.dot(h);
    if (u < -eps || u > 1.0 + eps) return false;
    const Eigen::Vector3d q = s.cross(e1);
    const double v = invDet * dir.dot(q);
    if (v < -eps || u + v > 1.0 + eps) return false;
    const double t = invDet * e2.dot(q);
    return t > eps;
}

double MeshSDF::signedness(const Eigen::Vector3d& y) const
{
    if (signMode_ == SignMode::Unsigned) return 1.0;
    int hits = 0;
    const Eigen::Vector3d origin = y + Eigen::Vector3d(0.0, 1e-9, 2e-9);
    for (int i = 0; i < mesh_.F.rows(); ++i) {
        const Eigen::Vector3i f = mesh_.F.row(i);
        if (rayIntersectsTriangle(origin, mesh_.V.row(f[0]), mesh_.V.row(f[1]), mesh_.V.row(f[2]))) {
            ++hits;
        }
    }
    return (hits % 2) == 1 ? -1.0 : 1.0;
}

SDFQuery MeshSDF::evalLocal(const Eigen::Vector3d& y, bool needGrad, bool needClosest) const
{
    SDFQuery q;
    double best = std::numeric_limits<double>::infinity();
    TriangleClosestPoint bestCp;
    int bestFace = -1;
    for (int i = 0; i < mesh_.F.rows(); ++i) {
        const Eigen::Vector3i f = mesh_.F.row(i);
        const auto cp = closestPointOnTriangle(y, mesh_.V.row(f[0]), mesh_.V.row(f[1]), mesh_.V.row(f[2]));
        if (cp.squaredDistance < best) {
            best = cp.squaredDistance;
            bestCp = cp;
            bestFace = i;
        }
    }

    const double sign = signedness(y);
    const double dist = std::sqrt(std::max(0.0, best));
    q.phi = sign * dist;
    q.closestFaceId = bestFace;
    if (needClosest) {
        q.closestPoint = bestCp.point;
        q.hasClosest = true;
    }
    if (needGrad) {
        Eigen::Vector3d dir = y - bestCp.point;
        if (dir.norm() > 1e-12) {
            dir.normalize();
            q.grad = sign > 0.0 ? dir : -dir;
        } else if (bestFace >= 0) {
            q.grad = mesh_.FN.row(bestFace).transpose().normalized();
        } else {
            q.grad = Eigen::Vector3d::UnitX();
        }
    }
    return q;
}

} // namespace bff_sdf::sdf

