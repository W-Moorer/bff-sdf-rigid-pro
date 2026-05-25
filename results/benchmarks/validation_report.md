# Benchmark Validation Report

- Generated `accuracy.csv`, `runtime.csv`, `sdf_resolution_study.csv`, `mesh_resolution_study.csv`, and `ablation.csv`.
- Sphere-plane area error: 0.00061283
- Sphere-sphere area error: 0.00140366
- Procedural fallbacks are used for bunny and mechanical-part scenes because no tracked model assets are available.
- Torus patch is a procedural open contact patch and is marked as a topology fallback rather than a BFF disk atlas validation.
