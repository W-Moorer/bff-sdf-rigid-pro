#include "test_main.h"

#include "atlas/BFFAdapter.h"
#include "atlas/PatchBuilder.h"

using bff_sdf::atlas::BFFAdapter;
using bff_sdf::atlas::BFFAdapterOptions;
using bff_sdf::atlas::PatchBuilder;
using bff_sdf::core::Mesh;

namespace {

Mesh makeDisk()
{
    Mesh m;
    m.V.resize(4, 3);
    m.V << 0, 0, 0,
           1, 0, 0,
           1, 1, 0,
           0, 1, 0;
    m.F.resize(2, 3);
    m.F << 0, 1, 2,
           0, 2, 3;
    m.computeNormals();
    return m;
}

} // namespace

int main()
{
    return test::run("test_bff_adapter", [] {
        const Mesh mesh = makeDisk();
        auto patch = PatchBuilder::buildWholeMeshPatch(mesh);
        REQUIRE(patch.diskLike);
        REQUIRE(patch.connected);
        REQUIRE(patch.eulerCharacteristic == 1);
        REQUIRE(patch.boundaryLoop.size() == 4);

        BFFAdapterOptions options;
        options.preferOfficialBff = false;
        const BFFAdapter adapter;
        const auto result = adapter.parameterize(patch, options);

        REQUIRE_MSG(result.success, result.message);
        REQUIRE(result.usedFallback);
        REQUIRE(result.UV.rows() == mesh.V.rows());
        REQUIRE(result.Fuv.rows() == mesh.F.rows());
        REQUIRE(result.UV.array().isFinite().all());
        REQUIRE(result.flippedTriangleCount == 0);

        const auto chart = adapter.buildChart(3, patch, options);
        REQUIRE(chart.id == 3);
        REQUIRE(chart.localToGlobalVertex.size() == 4);
        REQUIRE(chart.localFaceToGlobalFace.size() == 2);

        int triId = -1;
        Eigen::Vector3d bary;
        const Eigen::Vector2d center = (chart.UV.row(0) + chart.UV.row(1) + chart.UV.row(2)) / 3.0;
        REQUIRE(chart.uvAccel.locate(center, triId, bary));
        REQUIRE(triId == 0);
        REQUIRE_NEAR(bary.sum(), 1.0, 1e-12);
    });
}

