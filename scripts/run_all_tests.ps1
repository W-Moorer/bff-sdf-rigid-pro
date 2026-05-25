param(
    [string]$Config = "Release"
)

$ErrorActionPreference = "Stop"

vcpkg install --triplet x64-windows
cmake --preset windows-vcpkg
cmake --build --preset "windows-vcpkg-$($Config.ToLower())"
ctest --preset windows-vcpkg-release

