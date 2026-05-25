#include "test_main.h"

#include "core/MeshRepair.h"

using bff_sdf::core::Mesh;
using bff_sdf::core::MeshRepairReport;
using bff_sdf::core::repairMesh;
using bff_sdf::core::sharedEdgeOrientationConsistent;

int main()
{
    return test::run("test_mesh_repair", [] {
        Mesh dirty;
        dirty.V.resize(5, 3);
        dirty.V << 0, 0, 0,
                   1, 0, 0,
                   1, 1, 0,
                   0, 1, 0,
                   0, 0, 0;
        dirty.F.resize(4, 3);
        dirty.F << 0, 1, 2,
                   0, 3, 2,
                   0, 0, 1,
                   4, 1, 2;
        dirty.computeNormals();

        MeshRepairReport report;
        const Mesh repaired = repairMesh(dirty, {}, &report);

        REQUIRE(repaired.F.rows() == 2);
        REQUIRE(repaired.V.rows() == 4);
        REQUIRE(report.removedDegenerateFaces == 1);
        REQUIRE(report.removedDuplicateFaces == 1);
        REQUIRE(sharedEdgeOrientationConsistent(repaired));
        REQUIRE(repaired.boundaryLoops().size() == 1);
        REQUIRE(repaired.eulerCharacteristic() == 1);
        REQUIRE(repaired.normalsFiniteAndUnit());
    });
}

