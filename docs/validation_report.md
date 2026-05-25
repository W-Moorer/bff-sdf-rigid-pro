# Validation Report

This report is updated after each implementation milestone. Commands are run from the repository root.

## M0-M3: Build, Mesh, Pose, Analytic SDF, Patch/BFFAdapter, ChartEvaluator

Date: 2026-05-26

Commands:

```powershell
vcpkg install --triplet x64-windows
cmake --preset windows-vcpkg
cmake --build --preset windows-vcpkg-release
ctest --preset windows-vcpkg-release
```

Status:

- Build: passed with MSVC 19.44 and vcpkg `eigen3:x64-windows`.
- Tests: 5/5 passed.
- Passing tests:
  - `test_mesh`
  - `test_pose`
  - `test_sdf`
  - `test_bff_adapter`
  - `test_chart_evaluator`

Validated behavior:

- Closed cube has zero boundary edges and Euler characteristic 2.
- Open disk patch has one boundary loop and Euler characteristic 1.
- Face and vertex normals are finite and unit length.
- Sphere, plane, and box analytic SDF queries match formulas and finite-difference gradients.
- RigidPose local/world transform and gradient rotation round-trip correctly.
- PatchBuilder validates a disk-like patch.
- BFFAdapter exposes official-BFF integration and deterministic fallback parameterization; fallback UV for the disk patch is finite and unflipped.
- ChartEvaluator returns exact barycentric positions, UV coordinates, unit normals, and Jacobian finite-difference agreement.

Current limitations:

- The M0-M3 unit test intentionally disables official BFF execution to keep the test deterministic and independent of the local binary. Official command-line BFF integration is implemented in `BFFAdapter` but will be exercised in later integration/benchmark runs.
- High-order projection is not implemented yet.
- MeshSDF, GridSDF, contact integration, projection refinement, detector, and benchmark outputs are scheduled for M4-M12.

## M4-M8: Gap Field, Marching Triangles, Contact Area, Mesh/Grid SDF, Projection

Date: 2026-05-26

Commands:

```powershell
cmake --build --preset windows-vcpkg-release
ctest --preset windows-vcpkg-release
```

Status:

- Build: passed.
- Tests: 9/9 passed.
- New passing tests:
  - `test_gap_field`
  - `test_marching_triangles`
  - `test_contact_area`
  - `test_projection`

Validated behavior:

- Pulled-back gap field computes `g = Phi_B(R_B^T(x_A - t_B))`.
- `gradU = J_A^T nWorld` matches finite differences on a tilted chart.
- Marching triangles handles all-positive, all-negative, one-negative, two-negative, and exact-zero vertex cases without duplicate zero points.
- Sphere-plane cap integration meets analytic thresholds: contact area error below 3% and contour radius error below 2%.
- MeshSDF brute-force closest triangle and closed-mesh sign checks pass for procedural cube/sphere cases.
- GridSDF trilinear interpolation shows lower aggregate error at higher resolution and returns normalized gradients near the surface.
- Atlas-Linear projection agrees with exact triangle closest point and reduces low-resolution GridSDF gap error in the sphere test.

Current limitations:

- MeshSDF currently uses deterministic brute-force closest triangle search. It is API-compatible with a future BVH acceleration layer but not yet asymptotically optimized.
- GridSDF stores a dense cubic grid and is intended for controlled resolution studies, not production memory efficiency.
- Atlas-HO is an explicit fallback path to Atlas-Linear; no high-order PN/MLS/subdivision evaluator is implemented yet.

## M9-M12: Detector, Experiments, Ablations, Documentation

Date: 2026-05-26

Commands:

```powershell
cmake --build --preset windows-vcpkg-release
ctest --preset windows-vcpkg-release
build\Release\run_sphere_plane.exe
build\Release\run_sphere_sphere.exe
build\Release\run_benchmarks.exe
powershell -ExecutionPolicy Bypass -File scripts\run_all_tests.ps1
python scripts\run_rigid_benchmarks.py --config Release
```

Status:

- Build: passed.
- Tests: 9/9 passed.
- Sphere-plane demo: passed.
- Sphere-sphere demo: passed.
- Benchmark CSV generation: passed.
- One-command test script: passed.
- One-command benchmark script: passed.

Key generated outputs:

- `results/sphere_plane/config.json`
- `results/sphere_plane/metrics.csv`
- `results/sphere_plane/timing.csv`
- `results/sphere_plane/timing.json`
- `results/sphere_plane/contact_samples.ply`
- `results/sphere_plane/contact_samples.obj`
- `results/sphere_plane/contact_contour.obj`
- `results/sphere_plane/contact_contour_3d.obj`
- `results/sphere_plane/contact_contour_uv.csv`
- `results/sphere_plane/gap_field.csv`
- `results/sphere_plane/uv_gap_field.csv`
- `results/sphere_plane/atlas_debug.obj`
- `results/sphere_plane/report.md`
- `results/benchmarks/accuracy.csv`
- `results/benchmarks/runtime.csv`
- `results/benchmarks/sdf_resolution_study.csv`
- `results/benchmarks/mesh_resolution_study.csv`
- `results/benchmarks/ablation.csv`
- `results/benchmarks/validation_report.md`

Representative results:

- Sphere-plane area error: 0.000460444 in the standalone demo.
- Sphere-plane contact radius error: 0.0000872858 in the standalone demo.
- Sphere-sphere area error: 0.000859866 in the standalone demo.
- GridSDF resolution study shows projection-refined error below pure grid error for all reported resolutions.

BFF status:

- `run_sphere_plane.exe` used the official BFF command-line binary from `bff_official/binaries/windows-v1.6/`.
- Official BFF core files were not modified.
- Fallback radial/planar parameterization remains available and is explicitly reported when used.

Current limitations:

- The old procedural bunny/mechanical fallback has been superseded by the follow-up real-asset benchmark rows below.
- Torus-plane is a procedural open patch fallback and is not presented as a disk-like BFF topology validation.
- Atlas-HO remains a TODO.

## Follow-Up: BFF Coverage, Real Assets, BVH Baseline

Date: 2026-05-26

Commands:

```powershell
cmake --build --preset windows-vcpkg-release
ctest --preset windows-vcpkg-release
build\Release\run_benchmarks.exe
powershell -ExecutionPolicy Bypass -File scripts\run_all_tests.ps1
python scripts\run_rigid_benchmarks.py --config Release
```

Status:

- Build: passed.
- Tests: 9/9 passed.
- One-command test script: passed.
- One-command benchmark script: passed.
- BVH-accelerated MeshSDF closest-triangle query: implemented and checked against brute-force exact closest distance in `test_sdf`.
- STL reader: implemented for binary and ASCII STL, checked on tracked CC0 gear STL in `test_mesh`.
- Official BFF benchmark coverage:
  - sphere-plane disk patch: `official-bff`
  - sphere-sphere source disk patch: `official-bff`
  - ellipsoid-plane disk patch: `official-bff`
  - bunny real bottom patch: `official-bff`
  - torus patch: `fallback-planar` because the patch is not disk-like
  - mechanical gear patch: `fallback-planar` because official BFF UV validation does not pass on the STL-derived patch
- New output: `results/benchmarks/bff_vs_planar_ablation.csv`.
- New tracked real mechanical asset: `data/meshes/involute_gear_teeth_16_angle_20_cc0.stl`.

Representative results:

- Sphere-plane BFF area error: 0.00061283.
- Sphere-sphere analytic area error: 0.00140366.
- Sphere-sphere BVH exact-triangle baseline area error: 0.00585994.
- Bunny real patch: full mesh 28576 faces, extracted contact patch 3529 faces, official BFF chart succeeded.
- Mechanical gear: full mesh 1336 faces, extracted contact patch 448 faces, planar fallback recorded.
