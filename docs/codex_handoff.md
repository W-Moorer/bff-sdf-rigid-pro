# Codex Handoff

## Completed

- vcpkg manifest setup with `eigen3`.
- C++17 CMake build and test presets.
- Mesh topology utilities: normals, AABB, face adjacency, boundary loops, Euler characteristic.
- OBJ read/write and debug export helpers.
- RigidPose local/world transforms.
- Analytic Plane/Sphere/Box SDFs.
- MeshSDF BVH-accelerated closest triangle and closed-mesh ray sign.
- GridSDF uniform sampling, trilinear interpolation, central-difference gradient.
- PatchBuilder disk-like validation.
- BFFAdapter official command-line integration plus deterministic fallback.
- ContactChart, UVAcceleration, ChartEvaluator.
- Pulled-back SDF gap field and UV gradient.
- MarchingTriangles contact polygon and zero segment extraction.
- ContactIntegrator area/penetration/centroid/normal accumulation.
- Atlas-Linear ChartProjector with convergence counters and Atlas-HO fallback message.
- ContactDetector pipeline.
- Sphere-plane and sphere-sphere demos.
- Benchmark, resolution study, mesh study, and ablation CSV generation.
- BVH-accelerated exact triangle closest point inside `MeshSDF`.
- Real bunny benchmark patch from `bff_official/input/bunny.obj`.
- Tracked CC0 involute gear STL for mechanical-part benchmark.
- BFF-vs-planar same-scene ablation table.
- Documentation: theory, implementation plan, validation report, experiments, paper outline.

## BFF Status

Official BFF core code under `bff_official/` was not modified. The adapter invokes `bff_official/binaries/windows-v1.6/bff-command-line.exe` when present or a `BFF_COMMAND` environment variable is supplied. The Windows `cmd.exe` quote handling issue was fixed by wrapping the command line. The current sphere-plane run used `official-bff` and reports that in `results/sphere_plane/config.json`.

If the official binary is absent in another checkout, `BFFAdapter` falls back to deterministic radial/planar UVs and records the method as fallback. That fallback is an engineering continuity path, not a replacement for BFF in paper claims.

## Known Limitations

- Atlas-HO is not implemented.
- MeshSDF uses a deterministic in-memory AABB BVH, but it has not been benchmarked against industrial BVH libraries.
- Vertex-face/edge-edge discrete contact is represented by closest-triangle mesh projection, not a full pair-contact solver.
- Torus is non-disk-like and uses planar fallback in the benchmark.
- Mechanical gear is a real tracked CC0 STL asset, but the extracted STL patch currently falls back to planar because official BFF UV validation does not pass.
- Seam-crossing logic is counted in the interface but not fully implemented.
- No deformable dynamics, friction, self-contact, IPC, NCP/LCP, or FEM coupling.

## Reproduce

```powershell
vcpkg install --triplet x64-windows
cmake --preset windows-vcpkg
cmake --build --preset windows-vcpkg-release
ctest --preset windows-vcpkg-release
python scripts/run_rigid_benchmarks.py --config Release
```

## TODO

- Improve official BFF robustness for STL-derived mechanical patch charts.
- Implement seam-crossing continuity checks for multi-chart atlases.
- Add Atlas-HO projection through PN triangles, MLS, or subdivision surfaces.
- Add richer visualization scripts for paper figures.
- Expand BFF integration tests to parse official UVs on more patch topologies.
