# Bundled fonts

Fonts shipped with Stibium for the **Text (font)** node and `fab.shapes.text_font`.
They are referenced by the display names below (the filename minus `.ttf`).
The built-in **Antimony** font (Matt Keeter's geometric typeface) is not here —
it lives in the shape library and is selected as `'Antimony'` (or an empty name).

These were chosen for a *machinist-label* aesthetic and for **meshing**: bold,
simple, geometric glyphs with generous counters bake into clean signed-distance
grids. Full license texts are in `licenses/`.

| Display name | Weight | License | Reserved Font Name | Source |
|---|---|---|---|---|
| Overpass | variable (renders at default/Regular via stb_truetype) | OFL 1.1 **+ Apache-2.0** | `Overpass` (OFL) — taken under Apache, so no RFN restriction | [google/fonts](https://github.com/google/fonts/tree/main/ofl/overpass) |
| Archivo Black | Black | OFL 1.1 | `Archivo` | [Omnibus-Type/Archivo](https://github.com/Omnibus-Type/Archivo) |
| B612 Mono | Bold | OFL 1.1 | none | [polarsys/b612](https://github.com/polarsys/b612) |
| DejaVu Sans Mono | Bold | Bitstream Vera + public-domain additions | none (OFL-style) | [dejavu-fonts](https://github.com/dejavu-fonts/dejavu-fonts) |
| Saira Stencil One | Regular | OFL 1.1 | `Saira` | [Omnibus-Type/Saira](https://github.com/Omnibus-Type/Saira) |

## On the Reserved Font Names

Three faces carry an OFL Reserved Font Name. Stibium ships the **pristine,
unmodified** TTF and bakes each glyph into a signed-distance grid **on the
user's machine** (cached under `~/.stibium/fonts/`). A baked SDF grid is
rendered output — like a glyph atlas or a mesh — not a *Modified Version of the
Font Software*, so the RFN does not attach to it. The original files are
redistributed under their own licenses (see `licenses/`), which OFL, Apache-2.0
and the Bitstream Vera license all permit for bundling inside an application.

## Notes

- **Overpass** is a variable font; `stb_truetype` reads its default instance
  (Regular weight). A static Bold/Heavy instance would be crisper for the
  machinist look — a nice follow-up if a static TTF is sourced.
- To add a font: drop a `.ttf` here named for its display name and add its
  license to `licenses/`. The resolver (`fab.shapes._resolve_font`) finds it by
  name automatically.
