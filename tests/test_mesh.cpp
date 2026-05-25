#include "test_main.h"

#include "core/Mesh.h"

using bff_sdf::core::Mesh;

namespace {

Mesh makeOpenDisk()
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

Mesh makeClosedCube()
{
    Mesh m;
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

} // namespace

int main()
{
    return test::run("test_mesh", [] {
        const Mesh disk = makeOpenDisk();
        REQUIRE(disk.boundaryEdges().size() == 4);
        REQUIRE(disk.boundaryLoops().size() == 1);
        REQUIRE(disk.eulerCharacteristic() == 1);
        REQUIRE(disk.normalsFiniteAndUnit());
        REQUIRE(disk.bounds().valid());

        const Mesh cube = makeClosedCube();
        REQUIRE(cube.boundaryEdges().empty());
        REQUIRE(cube.boundaryLoops().empty());
        REQUIRE(cube.eulerCharacteristic() == 2);
        REQUIRE(cube.normalsFiniteAndUnit());
    });
}

