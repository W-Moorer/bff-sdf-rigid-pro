# Paper Outline

## Title

BFF-SDF Contact Atlas for Continuous Rigid Surface Contact

## Abstract Draft

We present a rigid geometric contact prototype that represents source surface patches as contact atlases and pulls target signed distance fields back to chart space. The method extracts contact regions and zero-level contours in two-dimensional charts, integrates contact area on the source surface, and refines candidate closest points by projecting SDF-guided samples back to a target triangle mesh. The current prototype focuses on rigid surface contact and provides analytic validations, GridSDF resolution studies, curved-chart quality ablations, real-mesh reference rows, and reproducible benchmark figures.

## Contributions

- A contact-atlas formulation for rigid surface contact using BFF charts as preprocessing.
- A pulled-back signed-gap field \(g(u)=\Phi_B(x_A(u))\) over source charts.
- Marching-triangle extraction of chart-space contact regions and contact boundaries.
- Atlas-Linear SDF-guided projection refinement for closest point, normal, and gap correction.
- Reproducible C++/CMake benchmark outputs for analytic, real-mesh, and stress/seam scenes.

## Introduction Logic

Rigid contact detection is often expressed as discrete mesh proximity. This prototype instead creates continuous chart-space fields over source patches. BFF is used as a preprocessing chart generator, while SDFs provide fast runtime evaluation. Projection refinement separates fast candidate detection from final closest-point correction.

## Related Work Placeholders

- Surface parameterization and Boundary First Flattening.
- Signed distance fields and grid SDF contact.
- BVH closest-point contact for triangle meshes.
- Continuous contact region extraction and contouring.
- Contact atlases for later deformable or material-coordinate extensions.

## Method Outline

- Mesh topology and disk-like patch validation.
- BFFAdapter and fallback parameterization.
- Chart evaluator and UV acceleration.
- Pulled-back SDF gap field and UV gradient.
- Marching-triangle contact extraction.
- Contact area integration on linear source triangles.
- MeshSDF, GridSDF, and Atlas-Linear projection refinement.
- ContactDetector runtime pipeline and debug exports.

## Experiment Outline

- Sphere-plane, tilted sphere-plane, sphere-sphere, unequal sphere-sphere, and ellipsoid-plane reference validation.
- Torus-plane procedural open patch fallback and two split-chart seam stress scenes.
- Bunny-plane, bunny-sphere, and bunny-tilted-plane real BFF input asset with extracted disk-like contact patch.
- Mechanical gear-plane and gear-sphere tracked CC0 asset with repaired/remeshed contact patch.
- GridSDF resolution study.
- Mesh resolution study.
- Ablations: no projection, no SDF, BFF versus planar curved chart quality, SDF resolution, mesh resolution.

## Limitations

- Atlas-HO is an interface placeholder; high-order surface projection is not implemented.
- MeshSDF uses a simple in-memory AABB BVH, not a tuned industrial BVH.
- Mechanical gear contact patch currently uses boundary-fan remeshing; a constrained planar triangulation would be stronger.
- Complex real-mesh references are deterministic subdivided-surface references, not experimental ground truth.
- The prototype does not solve deformable dynamics, friction, self-contact, NCP/LCP, IPC, or FEM coupling.

## Figures and Tables

- `results/figures/pipeline.svg`
- `results/figures/analytic_convergence.svg`
- `results/figures/bff_vs_planar_quality.svg`
- `results/figures/grid_sdf_refinement.svg`
- `results/figures/real_asset_contact_visualization.svg`
- `results/figures/runtime_breakdown.svg`
