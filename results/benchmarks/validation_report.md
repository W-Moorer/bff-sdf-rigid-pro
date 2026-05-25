# Benchmark Validation Report

- Generated `accuracy.csv`, `runtime.csv`, `asset_stats.csv`, `chart_quality.csv`, `sdf_resolution_study.csv`, `mesh_resolution_study.csv`, `ablation.csv`, and `bff_vs_planar_ablation.csv`.
- Unique main scenes: 14 (`analytic=6`, `real_mesh=5`, `stress=3`); benchmark rows including baselines: 19 (`analytic=11`, `real_mesh=5`, `stress=3`).
- Baselines represented: BVH exact target-mesh projection, pure GridSDF, MeshSDF/BVHSDF pullback, planar fallback.
- SDF resolution study is interpreted as coarse-grid error reduction with a high-resolution linear projection floor.
- Complex mesh rows use `subdivided_surface_reference`; area errors are no longer default-filled zeros.
- `scripts/generate_paper_figures.py` writes pipeline, convergence, chart-quality, GridSDF-refinement, real-asset, and runtime SVG figures under `results/figures/`.
- Sphere-plane medium BFF area error: 0.00061283
- Sphere-sphere equal analytic area error: 0.00140366
- Sphere-sphere MeshSDF/BVH baseline area error: 0.00585994
- Sphere-sphere MeshSDF/BVHSDF pullback area error: 0.00585994
- Sphere-sphere pure GridSDF baseline area error: 0.0112312
- Curved BFF-vs-planar chart quality rows: 10
- bunny_real_contact_patch: full_faces=28576, patch_faces=3529, disk_like=true, chart=official-bff, remeshed=false, projection_failures=0, repair_removed_degenerate=0, repair_flipped_faces=0, note=plane-normal plus local disk search; bff_pass=true; remeshed=false; planar_candidate_face_evaluations=95
- mechanical_gear_real_contact_patch: full_faces=1336, patch_faces=430, disk_like=true, chart=official-bff, remeshed=true, projection_failures=0, repair_removed_degenerate=0, repair_flipped_faces=0, note=plane-normal plus local disk search; bff_pass=true; remeshed=true; planar_candidate_face_evaluations=8016
- Mechanical gear asset source: koppi/involute-gear-collection, CC0-1.0, `angle-20/teeth-16.stl`.
