# Asset Sources

## Bunny

- Path used by benchmark: `bff_official/input/bunny.obj`
- Status: local BFF official input asset, kept ignored with `bff_official/`.
- Use: benchmark extracts a disk-like bottom contact patch and parameterizes it through `BFFAdapter`.
- Current benchmark status: official BFF chart succeeded for the extracted bunny patch.

## Mechanical Gear

- Tracked path: `data/meshes/involute_gear_teeth_16_angle_20_cc0.stl`
- Source repository: `https://github.com/koppi/involute-gear-collection`
- Source file: `angle-20/teeth-16.stl`
- License: CC0-1.0. A copy of the license is stored in `data/meshes/involute_gear_collection_LICENSE_CC0.txt`.
- Use: benchmark extracts a bottom contact patch from the gear STL.
- Current benchmark status: disk patch extraction succeeds, but official BFF UV validation falls back to planar on this STL-derived patch. The report records this explicitly.

