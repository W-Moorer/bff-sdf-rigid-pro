# Experiments

## Commands

```powershell
vcpkg install --triplet x64-windows
cmake --preset windows-vcpkg
cmake --build --preset windows-vcpkg-release
ctest --preset windows-vcpkg-release
build\Release\run_sphere_plane.exe
build\Release\run_sphere_sphere.exe
build\Release\run_benchmarks.exe
python scripts/generate_paper_figures.py
```

or:

```powershell
python scripts/run_rigid_benchmarks.py --config Release
```

## Outputs

Sphere-plane writes all required files under `results/sphere_plane/`, including `metrics.csv`, `timing.csv`, `timing.json`, `contact_samples.ply`, `contact_contour_3d.obj`, `contact_contour_uv.csv`, `gap_field.csv`, `uv_gap_field.csv`, `atlas_debug.obj`, and `report.md`.

Benchmarks write summary CSV files under `results/benchmarks/`:

- `accuracy.csv`
- `runtime.csv`
- `asset_stats.csv`
- `chart_quality.csv`
- `sdf_resolution_study.csv`
- `mesh_resolution_study.csv`
- `ablation.csv`
- `bff_vs_planar_ablation.csv`
- `validation_report.md`

Paper SVG figures are generated under `results/figures/`:

- `pipeline.svg`
- `analytic_convergence.svg`
- `bff_vs_planar_quality.svg`
- `grid_sdf_refinement.svg`
- `real_asset_contact_visualization.svg`
- `runtime_breakdown.svg`

Real-asset benchmark notes are recorded in `docs/asset_sources.md` and `results/benchmarks/validation_report.md`.

## Current Benchmark Scope

The benchmark executable currently reports 14 unique main scenes: 6 analytic/reference scenes, 5 real-mesh scenes, and 3 stress/seam scenes. Complex real-mesh rows use a deterministic subdivided-surface reference instead of a default zero-error placeholder.
`accuracy.csv` and `runtime.csv` include `integration_ms`, `projection_ms`, and `other_ms` timing fields used by `results/figures/runtime_breakdown.svg`.

Implemented baseline families:

- BVH exact target-mesh projection.
- Pure GridSDF pullback.
- MeshSDF/BVHSDF target query.
- Planar fallback chart.

Implemented ablations:

- no projection versus projection;
- no SDF brute-force projection count versus SDF-guided narrow-band projection count;
- official BFF versus fallback planar chart quality on curved patches;
- SDF resolution;
- mesh resolution.
