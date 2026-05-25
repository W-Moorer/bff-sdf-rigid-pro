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
