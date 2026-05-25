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
- `sdf_resolution_study.csv`
- `mesh_resolution_study.csv`
- `ablation.csv`
- `bff_vs_planar_ablation.csv`
- `validation_report.md`

Real-asset benchmark notes are recorded in `docs/asset_sources.md` and `results/benchmarks/validation_report.md`.
