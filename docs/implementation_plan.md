# Implementation Plan

Project: BFF-SDF Contact Atlas for Rigid Surface Contact

This plan follows the repository AGENTS.md requirements and the requested milestone order. The pre-existing root seed files and `bff_official/` remain ignored and are not uploaded. New implementation files live in tracked project directories.

## Repository and BFF Inspection

- Existing tracked root content before implementation: `README.md`, `LICENSE`, `.gitignore`.
- Existing local seed content: `AGENTS.md`, `BFF_SDF_Rigid_Contact_*.md`, `bff_official/`; all ignored.
- No existing `CMakeLists.txt`, `src/`, `tests/`, or `docs/` implementation existed before this plan.
- BFF official directory contains `include/bff`, `src/project/Bff.cpp`, `apps/command-line/src/CommandLine.cpp`, and Windows binaries in `bff_official/binaries/windows-v1.6/`.
- BFF code-level API exists through `bff::BFF`, but building the official library requires SuiteSparse/OpenBLAS setup. The local binary `bff-command-line.exe` is available.
- External libraries are introduced through vcpkg manifest mode. The local machine has `VCPKG_ROOT=C:\vcpkg`; `vcpkg.json` currently declares `eigen3`.

## BFF Adapter Strategy

1. Keep official BFF code untouched.
2. Implement `atlas::BFFAdapter` as a preprocessing adapter:
   - write a temporary patch OBJ;
   - call official `bff-command-line` when `BFF_COMMAND` is set or the local Windows binary exists;
   - parse `vt` UVs and face texture indices from the output OBJ;
   - validate finiteness and flipped UV triangles.
3. If the BFF executable is absent or fails, use deterministic fallback parameterization:
   - planar PCA-style projection for disk-like patches;
   - normalized orientation check and global reflection if needed;
   - record `usedFallback=true` and a diagnostic message.
4. Documentation and reports must clearly state whether BFF or fallback was used. Fallback is not claimed as final BFF integration.

## Milestones and Files

### M0 Build and Shared Test Harness

Files:
- `CMakeLists.txt`
- `CMakePresets.json`
- `vcpkg.json`
- `tests/test_main.h`
- `scripts/run_all_tests.ps1`
- `scripts/run_rigid_benchmarks.py`

Acceptance commands:
- `vcpkg install --triplet x64-windows`
- `cmake --preset windows-vcpkg`
- `cmake --build --preset windows-vcpkg-release`
- `ctest --preset windows-vcpkg-release`

### M1 Mesh, Geometry, Pose, Analytic SDF

Files:
- `src/core/GeometryTypes.h`
- `src/core/AABB.h`
- `src/core/Mesh.h`, `src/core/Mesh.cpp`
- `src/core/RigidPose.h`
- `src/io/MeshIO.h`, `src/io/MeshIO.cpp`
- `src/sdf/ISDF.h`
- `src/sdf/AnalyticSDF.h`, `src/sdf/AnalyticSDF.cpp`
- `tests/test_mesh.cpp`
- `tests/test_pose.cpp`
- `tests/test_sdf.cpp`

Acceptance:
- closed cube boundary edge count is zero;
- open disk has one boundary loop and Euler characteristic one;
- normals are finite and unit length;
- analytic SDF formulas and finite-difference gradients pass tolerance;
- world/local transform is exact within floating-point tolerance.

### M2 Patch and BFF Adapter

Files:
- `src/atlas/Patch.h`
- `src/atlas/PatchBuilder.h`, `src/atlas/PatchBuilder.cpp`
- `src/atlas/UVAcceleration.h`, `src/atlas/UVAcceleration.cpp`
- `src/atlas/BFFAdapter.h`, `src/atlas/BFFAdapter.cpp`
- `src/atlas/ContactAtlas.h`, `src/atlas/ContactAtlas.cpp`
- `tests/test_bff_adapter.cpp`

Acceptance:
- disk-like topology checks are enforced;
- UV output is finite;
- no flipped UV triangles after validation;
- topology coverage is complete;
- fallback status is exposed if official BFF cannot be invoked.

### M3 Chart Evaluator

Files:
- `src/atlas/ChartEvaluator.h`, `src/atlas/ChartEvaluator.cpp`
- `tests/test_chart_evaluator.cpp`

Acceptance:
- vertex and triangle-center barycentric evaluation is exact;
- `dX/du` matches finite differences;
- normals are consistently oriented and normalized.

### M4-M8 Contact Field, Marching, Integration, Mesh/Grid SDF, Projection

Files:
- `src/contact/ContactTypes.h`
- `src/contact/ContactGapField.h`, `src/contact/ContactGapField.cpp`
- `src/contact/MarchingTriangles.h`, `src/contact/MarchingTriangles.cpp`
- `src/contact/ContactIntegrator.h`, `src/contact/ContactIntegrator.cpp`
- `src/sdf/MeshSDF.h`, `src/sdf/MeshSDF.cpp`
- `src/sdf/GridSDF.h`, `src/sdf/GridSDF.cpp`
- `src/contact/ChartProjector.h`, `src/contact/ChartProjector.cpp`
- `tests/test_gap_field.cpp`
- `tests/test_marching_triangles.cpp`
- `tests/test_contact_area.cpp`
- `tests/test_projection.cpp`

Acceptance:
- pulled-back SDF matches analytic targets and UV gradients match finite differences;
- all marching triangle sign cases pass;
- sphere-plane contact area and radius meet thresholds on medium mesh;
- MeshSDF and GridSDF have deterministic convergence checks;
- Atlas-Linear projection agrees with exact triangle closest point.

### M9-M12 Detector, Experiments, Ablations, Documentation

Files:
- `src/contact/ContactDetector.h`, `src/contact/ContactDetector.cpp`
- `src/io/DebugExport.h`, `src/io/DebugExport.cpp`
- `src/experiments/ProceduralMeshes.h`, `src/experiments/ProceduralMeshes.cpp`
- `src/experiments/run_sphere_plane.cpp`
- `src/experiments/run_sphere_sphere.cpp`
- `src/experiments/run_benchmarks.cpp`
- `docs/theory.md`
- `docs/validation_report.md`
- `docs/paper_outline.md`
- `docs/codex_handoff.md`

Required generated outputs:
- `results/sphere_plane/config.json`
- `results/sphere_plane/metrics.csv`
- `results/sphere_plane/timing.csv`
- `results/sphere_plane/timing.json`
- `results/sphere_plane/contact_samples.ply`
- `results/sphere_plane/contact_samples.obj`
- `results/sphere_plane/contact_contour_3d.obj`
- `results/sphere_plane/contact_contour.obj`
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

Acceptance commands:
- `powershell -ExecutionPolicy Bypass -File scripts/run_all_tests.ps1`
- `python scripts/run_rigid_benchmarks.py --config Release`
- `build/Release/run_sphere_plane.exe` on Visual Studio generators, or the equivalent executable path reported by CMake.

## Validation Report Policy

After each milestone:
- build the affected targets;
- run the relevant tests;
- update `docs/validation_report.md` with command, status, and limitations;
- commit the milestone with a scoped message.

## Known Initial Risks

- Official BFF static library integration may be blocked by SuiteSparse availability. The command-line executable path is the first integration route.
- If the BFF binary is not present in a fresh clone, `BFFAdapter` will use fallback UV and report that state.
- Full high-order projection is out of current implementation scope; `Atlas-HO` extension points will be present and documented.
