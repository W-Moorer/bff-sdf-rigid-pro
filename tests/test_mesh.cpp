#include "test_main.h"

#include "core/Mesh.h"
#include "io/MeshIO.h"

#include <filesystem>

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

std::filesystem::path findRepoRoot()
{
    auto p = std::filesystem::current_path();
    for (int i = 0; i < 6; ++i) {
        if (std::filesystem::exists(p / "vcpkg.json")) return p;
        if (!p.has_parent_path() || p == p.parent_path()) break;
        p = p.parent_path();
    }
    return std::filesystem::current_path();
}

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

        Mesh gear;
        std::string error;
        const auto stlPath = findRepoRoot() / "data" / "meshes" / "involute_gear_teeth_16_angle_20_cc0.stl";
        REQUIRE_MSG(bff_sdf::io::readStl(stlPath.string(), gear, &error), error);
        REQUIRE(gear.V.rows() > 0);
        REQUIRE(gear.F.rows() > 0);
        REQUIRE(gear.normalsFiniteAndUnit());
    });
}
