# Workbench screenshot harness

> Companion to `VALIDATION.md` and `QA.md`. Renders every workbench
> panel + visualization mode to a PNG on disk so a reviewer has
> concrete artifacts for the **live visual aesthetic** check without
> running the app on a machine.

## What this harness does

Valenx already validates panel logic and the GPU render path headlessly:

- The `headless_ui_tests` name-filtered test set (165 tests as of this
  pass) drives every panel through `egui::Context::run` and asserts
  the draw fn never panics across representative states.
- The PBR shader is GPU-validated by
  `headless_pbr_render_shades_a_lit_quad` in
  `crates/valenx-app/src/wgpu_renderer.rs` — it builds the real
  `PBR_FORWARD_WGSL` pipeline on an off-screen wgpu device, shades a
  point-lit matte quad, reads back the pixels and asserts a real
  irradiance gradient. Ran on a real RTX 4070 (Vulkan).

**What was missing**: nobody could *see* what the workbench panels
look like without launching the app on graphics hardware. This
harness closes that gap by writing a PNG per panel to
`screenshots/` on demand.

## What this harness does NOT do

- It does **not** perform an aesthetic judgement. That is a human
  decision against the generated PNGs.
- It does **not** exercise every interaction state, every theme,
  every locale. One representative populated frame per panel.
- It does **not** screenshot the central 3-D viewport's wgpu PBR
  output — that path is covered by the GPU render-path validation
  test above.

## Running

```bash
cargo test -p valenx-app --test headless_screenshots
```

This is **scoped-safe** under the standing test lockdown — it
compiles + runs only `tests/headless_screenshots.rs`, never the
unfiltered `valenx-app` test set that holds the `rfd::FileDialog`
tests. The harness skips cleanly (no failure, no PNGs) when no wgpu
adapter is available — a CI sandbox without a GPU and without a
software fallback.

PNGs land under `screenshots/<workbench>/<panel>.png`. The
`screenshots/` directory is `.gitignore`-excluded — the harness
regenerates them on demand and they don't belong in git.

A `screenshots/INDEX.md` is emitted alongside, mapping each PNG to
its workbench, panel label, and a one-line description.

## What the harness renders

| Workbench | Panels |
|---|---|
| Mesh / CAD Toolbox | 1 host SidePanel + 10 sub-panels (Part, Draft, TechDraw, Assembly, Surface, CAM, Arch / BIM, Spreadsheet, Dock, Sketcher) |
| Genetics Workbench | 1 host SidePanel + 14 panels (Sequence, Alignment, Phylogenetics, Population Genetics, RNA Structure, RNA Designer, Molecular Dynamics, Cheminformatics, Macromolecular Structure, Quantum Chemistry, Genomics, Systems Biology, Docking, Gene Editing) |
| Aerodynamics / Wind Tunnel | 1 host SidePanel + 8 workflow sections (Body, Wind, Ground & wheels, Tunnel & mesh, Solver, Run, Results, Visualization) |
| **Total** | **35 PNGs** |

Each PNG is 1280 × 800 (`Rgba8Unorm`), one frame per panel.

## How it works

For each panel spec the harness:

1. Builds a fresh `ValenxApp` and toggles the relevant workbench on
   via the public `enable_*_workbench` setters.
2. Forces every `egui::CollapsingHeader` / popup / menu open through
   `egui::Memory::set_everything_is_visible(true)` so collapsed-by-
   default sections (every aero section, the mesh-toolbox CAD
   sub-panels) actually render their body and the screenshot is
   useful.
3. Runs one egui frame against an `egui::Context` with the screen
   rect bounded to 1280 × 800 so SidePanels mount at a stable width.
4. Tessellates the shapes, hands the paint jobs to a real
   `egui_wgpu::Renderer` pointed at an off-screen
   `Rgba8Unorm` colour texture (no surface, no window).
5. Encodes the render pass, submits, copies the texture to a
   256-byte-row-aligned buffer, maps the buffer, strips the row
   padding, encodes the tight RGBA buffer as a PNG with the `png`
   crate.

## Per-PNG assertions

- Every expected PNG exists on disk.
- Every PNG is > 1 KB (smaller would mean the encoder wrote a
  near-empty image).
- Every PNG was rendered at 1280 × 800.
- Every render produced more than 100 non-clear pixels (proves the
  draw fn actually drew something rather than no-op'ing into a
  uniform clear-colour frame).

## Honest caveats

- A panel that depends on a 3-D viewport overlay (the Aero
  Visualization section pushes fields *to* the viewport — it doesn't
  render them itself) renders only its controls + the status text,
  not the live overlay.
- `egui::ComboBox` dropdowns and menus may render their popup body
  off-screen even with `everything_is_visible` set; the dropdown
  trigger button is always visible.
- The harness uses `ValenxApp::default()` plus the panel's
  `Default` populators (the same state every `headless_ui_tests`
  uses). It does not run the panel's Run/Compute action before
  screenshotting — that would conflict with the Aero workbench's
  background-thread design. Post-run states with real solver output
  are out of scope for this single-frame harness.

## When to regenerate

Re-run the harness when:

- A panel's `draw` function changes layout, widget set, or section
  headings.
- The egui / eframe version is bumped.
- A workbench gains or loses a panel.
- The brand / theme tokens (`valenx-design-tokens`) change.

The output is deliberately not version-controlled — the source of
truth is the test, which produces deterministic output for a given
input state + harness version.
