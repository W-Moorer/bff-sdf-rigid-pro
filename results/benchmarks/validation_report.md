# Benchmark Validation Report

- Generated `accuracy.csv`, `runtime.csv`, `sdf_resolution_study.csv`, `mesh_resolution_study.csv`, `ablation.csv`, and `bff_vs_planar_ablation.csv`.
- Sphere-plane BFF area error: 0.00061283
- Sphere-sphere analytic area error: 0.00140366
- Sphere-sphere BVH exact-triangle baseline area error: 0.00585994
- Disk-like procedural patches request official BFF charts; non-disk topology falls back and is marked in the method name.
- bunny_plane_real_bff_official_input: full_faces=28576, patch_faces=3529, chart=official-bff
- mechanical_part_plane_real_involute_gear: full_faces=1336, patch_faces=448, chart=fallback-planar
- Mechanical gear asset source: koppi/involute-gear-collection, CC0-1.0, `angle-20/teeth-16.stl`.
