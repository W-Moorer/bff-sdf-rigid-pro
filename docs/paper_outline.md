# Paper Outline

## Title

BFF-SDF Contact Atlas for Continuous Rigid Surface Contact

## Abstract Draft

We present a rigid geometric contact prototype that represents source surface patches as contact atlases and pulls target signed distance fields back to chart space. The method extracts contact regions and zero-level contours in two-dimensional charts, integrates contact area on the source surface, and optionally refines candidate closest points by projecting SDF-guided samples back to a target chart or triangle mesh. The current prototype focuses on rigid surface contact and provides analytic validations, grid-SDF resolution studies, and reproducible benchmark outputs.

## Contributions

- A contact-atlas formulation for rigid surface contact using BFF charts as preprocessing.
- A pulled-back signed-gap field \(g(u)=\Phi_B(x_A(u))\) over source charts.
- Marching-triangle extraction of chart-space contact regions and contact boundaries.
- Atlas-Linear SDF-guided projection refinement for closest point, normal, and gap correction.
- Reproducible C++/CMake benchmark outputs for analytic and procedural scenes.

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

- Sphere-plane analytic area/radius validation.
- Sphere-sphere analytic circle validation.
- Ellipsoid-plane procedural geometry.
- Torus-plane procedural open patch fallback.
- Bunny-plane procedural fallback because no tracked model asset is available.
- Mechanical part-plane procedural plate fallback.
- GridSDF resolution study.
- Mesh resolution study.
- Ablations: no projection, no SDF, fallback chart, SDF resolution, mesh resolution.

## Limitations

- Atlas-HO is an interface placeholder; high-order surface projection is not implemented.
- MeshSDF currently uses brute-force closest triangle search rather than a BVH.
- Procedural fallbacks replace bunny and mechanical assets until tracked meshes are added.
- The prototype does not solve deformable dynamics, friction, self-contact, NCP/LCP, IPC, or FEM coupling.

## Figures and Tables

- Contact atlas UV gap heatmap.
- Sphere-plane contact contour in 3D and UV.
- GridSDF resolution error table.
- Projection refinement ablation table.
- Runtime breakdown table.
- Mesh resolution convergence table.

