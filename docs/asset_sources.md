# Asset Sources

## Bunny

- Path used by benchmark: `bff_official/input/bunny.obj`
- Status: local BFF official input asset, kept ignored with `bff_official/`.
- Use: benchmark extracts a contact-side disk patch using plane/normal and local disk criteria, then parameterizes it through `BFFAdapter`.
- Current benchmark status: official BFF chart succeeded for the extracted bunny patch.

## Mechanical Gear

- Tracked path: `data/meshes/involute_gear_teeth_16_angle_20_cc0.stl`
- Source repository: `https://github.com/koppi/involute-gear-collection`
- Source file: `angle-20/teeth-16.stl`
- License: CC0-1.0. A copy of the license is stored in `data/meshes/involute_gear_collection_LICENSE_CC0.txt`.
- Use: benchmark repairs the STL mesh, extracts a contact-side disk patch using plane/normal and local disk criteria, and remeshes the patch with a boundary fan if the original STL triangulation fails UV quality checks.
- Current benchmark status: official BFF chart succeeds after repair and boundary-fan remeshing of the extracted gear patch. The report records `remeshed=true`.
