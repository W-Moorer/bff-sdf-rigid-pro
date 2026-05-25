#include "sdf/MeshSDF.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

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
    buildBvh();
}

bff_sdf::core::AABB MeshSDF::faceBounds(int faceId) const
{
    bff_sdf::core::AABB box;
    const Eigen::Vector3i f = mesh_.F.row(faceId);
    box.expand(mesh_.V.row(f[0]));
    box.expand(mesh_.V.row(f[1]));
    box.expand(mesh_.V.row(f[2]));
    return box;
}

void MeshSDF::buildBvh()
{
    faceOrder_.resize(mesh_.F.rows());
    std::iota(faceOrder_.begin(), faceOrder_.end(), 0);
    bvh_.clear();
    if (!faceOrder_.empty()) buildNode(0, static_cast<int>(faceOrder_.size()));
}

int MeshSDF::buildNode(int begin, int end)
{
    const int nodeId = static_cast<int>(bvh_.size());
    bvh_.push_back({});
    BVHNode node;
    node.begin = begin;
    node.end = end;
    for (int i = begin; i < end; ++i) {
        const auto fb = faceBounds(faceOrder_[i]);
        node.bounds.expand(fb.min);
        node.bounds.expand(fb.max);
    }

    const int count = end - begin;
    if (count <= 8) {
        node.leaf = true;
        bvh_[nodeId] = node;
        return nodeId;
    }

    const Eigen::Vector3d extent = node.bounds.extent();
    int axis = 0;
    if (extent.y() > extent.x()) axis = 1;
    if (extent.z() > extent[axis]) axis = 2;
    const int mid = begin + count / 2;
    std::nth_element(faceOrder_.begin() + begin, faceOrder_.begin() + mid, faceOrder_.begin() + end,
                     [&](int lhs, int rhs) {
                         const auto lb = faceBounds(lhs);
                         const auto rb = faceBounds(rhs);
                         return lb.center()[axis] < rb.center()[axis];
                     });

    node.left = buildNode(begin, mid);
    node.right = buildNode(mid, end);
    bvh_[nodeId] = node;
    return nodeId;
}

double MeshSDF::aabbDistanceSquared(const bff_sdf::core::AABB& box, const Eigen::Vector3d& p) const
{
    double d2 = 0.0;
    for (int a = 0; a < 3; ++a) {
        if (p[a] < box.min[a]) {
            const double d = box.min[a] - p[a];
            d2 += d * d;
        } else if (p[a] > box.max[a]) {
            const double d = p[a] - box.max[a];
            d2 += d * d;
        }
    }
    return d2;
}

void MeshSDF::closestInNode(int nodeId,
                            const Eigen::Vector3d& y,
                            double& best,
                            TriangleClosestPoint& bestCp,
                            int& bestFace) const
{
    if (nodeId < 0) return;
    const BVHNode& node = bvh_[nodeId];
    if (aabbDistanceSquared(node.bounds, y) > best) return;

    if (node.leaf) {
        for (int i = node.begin; i < node.end; ++i) {
            const int faceId = faceOrder_[i];
            const Eigen::Vector3i f = mesh_.F.row(faceId);
            const auto cp = closestPointOnTriangle(y, mesh_.V.row(f[0]), mesh_.V.row(f[1]), mesh_.V.row(f[2]));
            if (cp.squaredDistance < best) {
                best = cp.squaredDistance;
                bestCp = cp;
                bestFace = faceId;
            }
        }
        return;
    }

    const double dl = node.left >= 0 ? aabbDistanceSquared(bvh_[node.left].bounds, y) : std::numeric_limits<double>::infinity();
    const double dr = node.right >= 0 ? aabbDistanceSquared(bvh_[node.right].bounds, y) : std::numeric_limits<double>::infinity();
    if (dl < dr) {
        closestInNode(node.left, y, best, bestCp, bestFace);
        closestInNode(node.right, y, best, bestCp, bestFace);
    } else {
        closestInNode(node.right, y, best, bestCp, bestFace);
        closestInNode(node.left, y, best, bestCp, bestFace);
    }
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
    closestInNode(bvh_.empty() ? -1 : 0, y, best, bestCp, bestFace);

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
