# UI modernization plan

Prioritized glow-up plan, grounded in the actual source (paths cited).
Compiled 2026-07-12. Framework verdict up front: **keep Qt, port to Qt6
as its own later milestone, reject QML/full-rewrite** (the canvas is
deeply QGraphicsScene-native; everything modern below is achievable in
QPainter, and QGraphicsScene is fully supported in Qt6).

## Ground truth

- Qt5 Widgets + QGraphicsScene canvas; the 3D view already uses
  QOpenGLWidget (good for Qt6). GL 2.1 / GLSL 120 shaders.
- Palette (`app/app/colors.cpp`) is verbatim base16-default-dark (2012).
  It's exported to Python by name (`Colors::PyColors()`), so values can
  be retinted freely as long as hue identity is preserved (teal stays
  teal, green stays green — node scripts reference by name).
- No app-wide QSS/QPalette/Fusion → menus, dialogs, scrollbars render
  platform-native (usually light) against the dark canvas. Biggest
  single "dated" signal.
- No HiDPI flags at all. No QIcon usage anywhere. Script editor
  hardcodes "Courier".
- Node bodies: flat fill, faked title band (rounded rect + fat line),
  selection "glow" is a 10px 50%-white stroke. Ports are bare 10px
  squares. Wires are flat 4px solid strokes. Dot grid hard-vanishes
  past 5000px. Zoom unclamped.
- Qt6 debt markers (small!): qSort (inspector/frame.cpp), QWheelEvent
  delta()/pos() (canvas_view.cpp, viewport/view.cpp), QRegExp
  (populate.cpp), QTextCodec (main.cpp), setTabStopWidth +
  QFontMetrics::width (script/editor.cpp).

## Proposed palette

Keep the base16 structure (everything flows through `Colors::`); move to
a blue-tinted dark neutral ramp + cleaner accents. Hue identity
preserved for all semantic colors; blue becomes the single UI accent
for selection/focus/hover (replacing white-at-50%-alpha everywhere).

| Name | Now | Proposed | Role |
|---|---|---|---|
| base00 | #151515 | #111318 | canvas/editor background |
| base01 | #202020 | #191c23 | node body fill |
| base02 | #303030 | #232833 | editor fields |
| base03 | #505050 | #3a4150 | grid, borders, disabled |
| base04 | #b0b0b0 | #aab3c4 | body text |
| base05 | #d0d0d0 | #ccd3e0 | selected outline, bright text |
| base06 | #e0e0e0 | #e4e8f0 | titles |
| base07 | #f5f5f5 | #f5f7fb | highlights |
| red | #ac4142 | #e5484d | errors |
| orange | #d28445 | #e8883f | int |
| yellow | #f4bf75 | #f2c14e | float |
| green | #90a959 | #8fbf6b | Shape (stays green) |
| teal | #75b5aa | #4fc1b0 | handles (stays teal) |
| blue | #6a9fb5 | #5aa2e8 | links + NEW UI accent |
| violet | #aa759f | #b183d6 | reserved |
| brown | #8f5536 | #c08552 | string |

Also fix `Colors::adjust()` — multiplicative scaling can't brighten dark
colors properly; use HSL lightness shift.

## Tier A — weekend polish (each S; ~a week total, zero architecture)

- **A1 app-wide dark theme**: Fusion style + QPalette from `Colors::` +
  ~100-line QSS (menus/scrollbars/tooltips). Touch: `app/app/app.cpp`,
  new `app/app/theme.cpp`; delete the local hack in `window/quad.cpp`.
- **A2 palette swap**: the 16 hex literals in `colors.cpp` + adjust()
  fix. Do first; everything retints automatically.
- **A3 HiDPI flags**: AA_EnableHighDpiScaling + AA_UseHighDpiPixmaps in
  `main.cpp`; QSurfaceFormat setSamples(4) for the GL view.
- **A4 wires**: gradient stroke source-color → target-color, round caps,
  soft shadow pass, accent glow on hover (replaces the 20px white
  underlay). One function: `connection/base.cpp paint()`.
- **A5 ports**: circles not squares, contrast outline, hover ring,
  bigger hit area. `datum_port.cpp`.
- **A6 node bodies**: manual 3-layer drop shadow (NOT
  QGraphicsDropShadowEffect — rasterization foot-gun), subtle vertical
  gradient fill, proper top-rounded title band tinted by output type,
  accent-blue selection outline + soft glow. `inspector/frame.cpp
  paint()`.
- **A7 grid + rubber-band + zoom clamp**: two-level dot pitch with
  zoom-faded alpha (kill the hard cutoff), accent rubber-band, clamp
  zoom to [0.08, 4]. `canvas_view.cpp`.
- **A8 script editor**: bundle JetBrains Mono (OFL) with FixedFont
  fallback, retint highlighter, current-line highlight, line-number
  gutter. `script/editor.cpp`, `script/syntax.cpp`, new qrc.
- **A9 viewport chrome**: gradient background instead of pure black,
  saturated axis colors + X/Y/Z tip labels, monospace coordinate
  readout. `viewport/view.cpp`.

## Tier B — the big glow-up (M items, 2–4 weeks cumulative)

- **B1 shader modernization** (GLSL 120-compatible, high impact):
  hemispheric ambient, fresnel rim light, cheap depth-tap AO (depth
  texture already bound; needs a pixel_size uniform), gamma-correct
  output, optional matcap mode behind a menu toggle. Touch:
  `app/gl/shaded.frag`, `viewport/image.cpp`, `viewport/gl.cpp`.
- **B2 zoom LOD**: below ~0.4 draw nodes as solid type-tinted name
  cards, below ~0.25 wires drop to 1px / ports skip. Uses
  `option->levelOfDetailFromTransform`. How Blender/ComfyUI stay
  legible zoomed out.
- **B3 icons**: Lucide (ISC) or Phosphor (MIT) SVGs via Qt5Svg — menus,
  the four hand-painted inspector buttons, window icons. ~20 icons.
- **B4 UI font**: bundle Inter (OFL) as app font; graphics items
  inherit automatically; retest layout padding constants.
- **B5 fuzzy-search add-node palette**: Blender/ComfyUI-style popup
  (QLineEdit + filtered list) replacing the nested QMenu; refactor
  `populate.cpp` to expose a flat list, keep QMenu for the menubar.
  Biggest daily-use feel upgrade.
- **B6 motion**: 80–120ms QVariantAnimation on hover/selection alphas;
  zoomTo easing to ~180ms OutCubic.

## Tier C — structural

- **C1 minimap**: second non-interactive QGraphicsView sharing the
  scene, forced LOD-card mode, click-to-jump. Careful: CanvasView moves
  sceneRect while panning; minimap needs its own fitted rect.
- **C2 Qt6 port** — do it as a standalone milestone after A/B; nothing
  above requires it. Inventory is unusually small (list in Ground
  truth); QGraphicsScene survives; watch GL 2.1 context defaults
  (may need explicit CompatibilityProfile). Estimate M, not L. Buys
  per-monitor HiDPI, modern font rendering, platform health (esp.
  macOS, where Qt5 is legacy in homebrew).
- **C3 QML/full rewrite: REJECTED** — no QGraphicsScene equivalent;
  reimplementing items/hit-testing/undo for zero functional gain.
  Study (don't adopt): paceholder/nodeeditor (pipe painting,
  style-JSON), jchanvfx/NodeGraphQt (LOD, search-popup UX).
- **C4 docked single-window mode** (optional, discuss first): today
  every canvas/viewport/script is a separate QMainWindow — very 2014.
  QDockWidget default layout is pure Widgets work but touches window
  management broadly. L, medium-high risk.

## Sequencing

1. A2 palette → A1 theme → A3 HiDPI (foundation)
2. A4 wires → A5 ports → A6 nodes → A7 grid (canvas facelift)
3. A8 editor + A9 viewport
4. B1 shader ∥ B4 font → B3 icons → B6 motion
5. B5 search palette → B2 LOD → C1 minimap
6. C2 Qt6 milestone; C4 discuss.

Verify before vendoring (stated from model knowledge): licenses of
Lucide/Phosphor/Inter/JetBrains Mono/nodeeditor/NodeGraphQt, and Qt6
porting-guide specifics.
