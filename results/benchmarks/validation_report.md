# Benchmark Validation Report

- Generated `accuracy.csv`, `runtime.csv`, `sdf_resolution_study.csv`, `mesh_resolution_study.csv`, `ablation.csv`, and `bff_vs_planar_ablation.csv`.
- Sphere-plane BFF area error: 0.00061283
- Sphere-sphere analytic area error: 0.00140366
- Sphere-sphere BVH exact-triangle baseline area error: 0.00585994
- Disk-like procedural patches request official BFF charts; non-disk topology falls back and is marked in the method name.
- bunny_plane_real_bff_official_input: full_faces=28576, patch_faces=3529, chart=official-bff, remeshed=false, repair_removed_degenerate=0, repair_flipped_faces=0, patch_selection=plane-normal plus local disk search; bff_pass=true; remeshed=false; planar_candidate_face_evaluations=95
- mechanical_part_plane_real_involute_gear: full_faces=1336, patch_faces=430, chart=official-bff, remeshed=true, repair_removed_degenerate=0, repair_flipped_faces=0, patch_selection=plane-normal plus local disk search; bff_pass=true; remeshed=true; planar_candidate_face_evaluations=8016
- Mechanical gear asset source: koppi/involute-gear-collection, CC0-1.0, `angle-20/teeth-16.stl`.
