# Frontend polish pass — 2026-05-23

## Scope

A systemic UX polish layer over the existing pre-alpha Valenx shell:
tokens consistency, tooltips, keyboard shortcuts, undo/redo,
friendly error states, accessibility, first-run onboarding,
per-panel contextual help, plus a follow-up pass adding a curated
icon set, broader undo/redo coverage, automatic friendly-error
mapping, more Mesh Toolbox tooltips, and subtle workbench-switch
animations. No new physics, no panel restructures — the goal is to
make the *existing* surface usable without reading the source code.

## Coverage table

| Category | Coverage | Concrete metric |
|---|---|---|
| Tokens | 100% (foundation) | New families `light_surface`, `light_text`, `hc_surface`, `hc_text`, `hc_accent` in `tokens.json`. Build script emits 5 new submodules. Genetics + aero error / OK colour lines route through `accent::ERROR` / `accent::SUCCESS` (was hardcoded hex). |
| Themes | 3 / 3 variants | `ThemeVariant::Dark` / `Light` / `HighContrast` in `valenx-app::theme`. Selectable in Settings → Appearance. HC palette verified AAA contrast on every text-tier × surface-tier pair. |
| Font scale | ✅ | 0.75–2.0× slider in Settings → Appearance, applied at theme-apply time to every `egui::Style::text_styles` entry. |
| Tooltips | ~95% (genetics + aero + every Mesh Toolbox operational control) | Every Run button + every sequence input across all 14 genetics panels (one stroke via `common::run_button` + `common::seq_input` helpers). ~30 hand-written on Wind / Tunnel / Solver / Body aero sections. ~20 on the genetics tool-selectors (MD, Genomics, Phylogenetics, Popgen, Docking, Biostruct, Gene Editing). Workbench chip-selectors route through `panel_help::short_summary`. **Second pass** added ~70 tooltips across Mesh Toolbox operational sections: Transformations (Translate/Scale/Rotate/Mirror XYZ inputs), Cut Plane (point + normal), Repair (merge tolerance), Mesh Tools (decimate/smooth params), Export (STL/OBJ/PLY/3MF), CAM Operations (Profile/Pocket/Drill common inputs + spindle/feed/Safe-Z + Step-down/Depth), Arch BIM (Wall + Slab parameters), Part Design (Save/Load/Import/Export STEP), Sketcher (tool palette + click input). **Third pass** added ~200 more tooltips: every CAM op-kind variant's per-op inputs (Pocket / Drill / Face / Adaptive Clearing / Helical Bore / Plunge Rough / Ramp Entry / Peck Drill Full / Contour 2D-3D / Engrave / Scribe / Spiral Pocket / Trochoidal Slot / Waterline 3D / Slot / Thread Mill / Rest Machining — with kind-specific step-down / step-over / depth descriptions), every Surface tool (NURBS curve / surface, Coons fill, Sew, Trim, KnotOps, SSI, Fit, Ruled) — every degree / CP / knot / weight / tolerance input, TechDraw (sheet size + title block + view position / scale / parametric, every dimension + chain + balloon + leader input + per-style descriptions), Assembly (part primitives, mate kinds, joint kinds — per-DOF hover descriptions), Spreadsheet (sheet picker + add/remove + dims + Set / Clear / Re-evaluate). Read-only display labels intentionally still un-annotated. |
| Shortcuts | 11 / 11 ShortcutAction variants | `Ctrl+P` palette, `Ctrl+1/2/3` workbench switch, `Ctrl+R` run, `Ctrl+S` save, `Ctrl+Z`/`Ctrl+Y`/`Ctrl+Shift+Z` undo/redo, `F1` contextual help, `?` toggle cheat-sheet, `Esc` cancel. All collision-tested via the `all_actions_have_unique_bindings` test. |
| Cheat-sheet | ✅ overlay | `?` toggles a 3-column grid (binding / action / description) enumerating `ShortcutAction::ALL`. Help menu has an entry. Open state persists via `settings.keyboard_shortcuts_overlay_open`. |
| Undo / Redo | **14 / 14** genetics panels + Aero workbench + Sketcher + Part Design (18 panels total) | `valenx-app::undo::History<T>` snapshot stack (64-deep, duplicate-collapse). **Second pass** extended this from the 3 first-pass panels (Sequence / Alignment / RnaStructure) to all 14 genetics panels (added Phylogenetics, Popgen, MD, Cheminformatics, Biostruct, Qchem, Genomics, Sysbio, Docking, GeneEditing, RnaDesigner) **plus** the Wind-Tunnel `AeroWorkbenchState`. **Third pass** wired the CAD-side Sketcher (`History<valenx_sketch::Sketch>` snapshotted before every click-add / constraint / Phase-12 primitive / toggle / sketch op) **and** Part Design (`History<valenx_feature_tree::FeatureTree>` snapshotted before every Add Sketch / 16 Add-Feature buttons / Suppress / Delete / Imported(Advanced) STEP-IGES). Unblocked by deriving `PartialEq` on `Sketch` + `FeatureTree` + their transitive sub-types (float fields use IEEE 754 semantics — documented NaN-trip degradation). Each panel exposes `record()` + `undo_edit` / `redo_edit` / `can_undo` / `can_redo` methods. The host's `try_undo_in_active_panel` dispatcher reaches every editor-state surface (and falls through Sketcher → Part Design when the Mesh Toolbox is open). Read-only result panels (`molecule_view`) intentionally still skipped. Every panel surfaces inline `↶` / `↷` buttons. |
| Friendly errors | **All ~30 panels** with Run actions | `sequence::friendly_error` (bespoke for primer-design / IUPAC / empty-input) + `aero::panels::friendly_aero_error` (bespoke for no-body / invalid-case / sweep-range) **stay** — they carry domain knowledge. **Second pass** added `genetics::common::friendly_error` — a pattern-matching mapper covering the common failure modes ("empty" / "invalid" / "not found" / "timeout" / "OOM" / "didn't converge" / "parse / malformed"). `common::error_line` now automatically appends the matching recovery-hint line under every error site, so all 14 genetics panels + every other panel using the helper gets friendly errors with zero per-panel rewrites. Unrecognised messages pass through unchanged — we never lose the raw text. |
| Icon set | **`valenx-icons` crate** with 7 families, 51 glyphs | Curated subset of Unicode glyphs covering run / file / edit / status / nav / view / domain icon families. Used by `common::run_button` (every Genetics Run button now shows `▶ Run`), `common::undo_redo_inline` (every inline undo/redo pair shows `↶ ↷`), `common::error_line` / `ok_line` (✖ and ✓ marks). License-clean (Unicode standard glyphs, no third-party icon font bundled). The mapping is semantically named so a future swap to rasterised PNGs / SVGs is a single-file change. |
| Animations | Subtle fade transitions on panel-body switches | `egui::Context::animate_bool_with_time` keyed on the active-panel label gives a 0.15 s fade-in when the user switches genetics panels (each new id starts at false, animates to true). The Aero workbench fades its body in 0.18 s when the panel opens via Ctrl+3. egui-native; no external animation runtime. |
| Accessibility | ✅ HC palette + font scale | `ThemeVariant::HighContrast` thickens stroke widths + highlights focus rings in pure-yellow. Font-scale slider 0.75–2.0×. WCAG AAA contrast verified (new `contrast_audit.rs` tests). AccessKit integration deferred — honestly named in the residue. |
| Welcome tour | ✅ 3 steps | Auto-opens on first launch (gated by `settings.welcome_tour_completed`). Step 1 "Welcome to Valenx", Step 2 "Three workbenches (Ctrl+1/2/3)", Step 3 "Shortcuts to know". Re-openable from Help → Welcome tour…. |
| Per-panel help | ~15 panels | F1 opens `panel_help::render_help_window` with a per-panel body (markdown-ish: `#` heading, `##` subhead, `-` bullets). Catalogue entries for every genetics panel + every aero workflow section + a generic fallback. |

## Architecture

Nine new modules under `crates/valenx-app/src/`:

- `theme.rs` — `ThemeVariant` enum + `ResolvedTheme` palette + `apply()` fn.
- `tooltips.rs` — `ResponseExt` trait (`.tt(...)`, `.tt_with_unit(...)`), `unit_str` lookup.
- `shortcuts.rs` — `ShortcutAction` enum + `poll_shortcut()` polled-input layer.
- `undo.rs` — generic bounded `History<T>` snapshot stack.
- `panel_help.rs` — per-panel help body + `short_summary` + window renderer.
- `keyboard_help.rs` — cheat-sheet overlay.
- `welcome_tour.rs` — 3-step orientation popup.

Plus extensions in:

- `valenx-design-tokens/tokens.json` + `build.rs` + `tokens.schema.json` — 4 new colour families.
- `valenx-design-tokens/tests/contrast_audit.rs` — 3 new HC + Light contrast tests.
- `valenx-app/src/settings.rs` — 4 new fields with serde-default migration.
- `valenx-app/src/update.rs` — shortcut dispatcher, overlay rendering, Help menu entries.
- `valenx-app/src/genetics/{common,sequence,alignment,rnastruct,cheminf,md,genomics,phylogenetics,popgen,docking,biostruct,genediting}.rs` — tooltips + `run_primary_shortcut` hooks.
- `valenx-app/src/aero/panels.rs` — tooltips on every wind/tunnel/solver/body input + `start_run_from_shortcut` public entry + `friendly_aero_error` mapper + Cancel button next to the spinner.

## Shortcuts table

| Binding | Action | Description |
|---|---|---|
| `Ctrl+P` | Open command palette | Fuzzy-search every action across all workbenches. |
| `F1` | Show this panel's help | Pop up help text for the currently focused panel. |
| `?` | Toggle keyboard cheat-sheet | Open / close the keyboard-shortcut cheat-sheet. |
| `Ctrl+1` | Show Mesh Toolbox | Switch to the Mesh Toolbox (Part / Draft / Assembly / CAM / etc.). |
| `Ctrl+2` | Show Genetics workbench | Switch to the Genetics workbench (14 computational-biology panels). |
| `Ctrl+3` | Show Wind Tunnel workbench | Switch to the Wind Tunnel workbench (3-D external-aero CFD). |
| `Ctrl+R` | Run / Compute | Triggers the active panel's main button (Translate / Fold / Solve / case Run). |
| `Ctrl+S` | Save | Save the project / settings (where applicable). |
| `Ctrl+Z` | Undo | Reverse the most recent edit in the active panel. |
| `Ctrl+Y` | Redo | Reapply the undone edit (also `Ctrl+Shift+Z`). |
| `Esc` | Cancel | Stop a long-running solve / sweep / run. Closes overlays first. |

## What this pass is NOT

Honest residue, named plainly:

- **Not a brand identity.** The High-Contrast palette is a verified-contrast scheme (pure-black + yellow / cyan / green) — not a Pantone-trained, motion-designed brand.
- **Not rasterised iconography.** The `valenx-icons` crate uses Unicode glyphs from egui's bundled font rather than rasterised SVGs / PNGs. They are *consistent* across the app and *named semantically*, so a future swap to Material or Feather rasterised icons is a single-file change — but they're not vector art.
- **Not Adobe After Effects motion.** The fade-in animations on panel-switch are subtle 150-180 ms opacity ramps. No spring physics, no slide-and-bounce, no animated chip-selector indicator (egui's `selectable_label` doesn't expose the per-item rect, so the active-chip slider that desktop apps usually have would need a custom widget — deferred).
- **Not 100% Mesh Toolbox tooltip coverage.** Operational controls (everything the user clicks to compute something) are now covered end-to-end across CAM, Surface, TechDraw, Assembly, Spreadsheet, Sketcher, Part Design, Transformations, Cut Plane, Repair, Mesh Tools, Export, Arch BIM. Pure display labels (counts, status text, generated dimension values) intentionally stay un-annotated. Architecture / sheet-metal / fasteners / robotics workbenches got first-pass coverage on their primary inputs but their full tail is still residue (not load-bearing for the common journey).
- **Sketcher / Part Design undo IS now wired.** Third pass added `History<Sketch>` and `History<FeatureTree>` snapshots in front of every mutating action with inline `↶ ↷` and host-level Ctrl+Z / Ctrl+Y. The snapshot strategy is structural `PartialEq` derived on `Sketch` and `FeatureTree` and all their sub-types; float fields use IEEE 754 semantics so a NaN snapshot fails to dedupe (DragValue widgets prevent NaN entry in practice — documented inline).
- **Not AccessKit integration.** egui-side AccessKit is a separate piece of work — what shipped is per-action keyboard navigation + WCAG-AAA visual contrast.

## Test gates

```
cargo test -p valenx-app headless_ui_tests                                       # 165 / 165 green (zero regressions)
cargo test -p valenx-app --test headless_screenshots                             # 1 / 1 green (35 PNGs still produced)
cargo test -p valenx-design-tokens                                               # 8 / 8 green (3 new HC + Light contrast)
cargo test -p valenx-icons                                                       # 5 / 5 green (icon-set smoke + non-empty + distinctness)
cargo check --workspace                                                          # clean
cargo clippy -p valenx-app -p valenx-design-tokens -p valenx-icons --all-targets -- -D warnings  # clean
cargo doc --workspace --no-deps                                                  # clean
```

## LOC added

First pass — ~2.7k LOC across 8 commits:

- Tokens + theme module: ~430
- Shortcuts + cheat-sheet: ~370
- Undo + tooltips + panel_help + welcome_tour: ~920
- Wiring (settings / lib / update / common / per-panel undo): ~700
- Aero + genetics tooltips + dispatchers: ~300

Second pass — ~1.5k LOC across the follow-up commits:

- `valenx-icons` crate: ~250 (7 families × ~7 icons, 5 unit tests, `label` formatter)
- Genetics common helpers (`friendly_error`, `undo_redo_inline`, `undo_redo_row`, `PanelHistory` trait, icon-prefixed buttons): ~150
- 11 genetics panels' `Snapshot` structs + `undo_edit` / `redo_edit` / record-on-Run wiring: ~700
- Aero workbench `History<WindTunnelForm>` + record-on-Run + inline undo/redo: ~70
- Update.rs dispatcher extended to 14 + Aero: ~50
- Mesh Toolbox tooltip pass (Transformations / Cut / Repair / Mesh Tools / Export / CAM / Arch Wall+Slab / Part Design / Sketcher): ~250
- Fade-in animation for genetics panel switches + aero workbench open: ~30

Third pass — final tail cleanup, ~1k LOC across one commit:

- `valenx-sketch` + `valenx-feature-tree`: derived `PartialEq` on `Sketch`, `Constraint`, every geometric primitive (`Point2`, `Line2`, `Circle2`, `Arc2`, `BSpline2`, `Ellipse2`, `EllipticalArc2`, `Entity`) and on `FeatureTree`, `TreeEntry`, `Value`, every `*Params` struct, `TransformOp`, `HoleDepthMode`, `Feature` — unblocks `History<T>`. Float fields (`Sketch::vars`, knot vectors, depths) use IEEE 754 semantics so a NaN snapshot fails to dedupe; in practice the panels' DragValue widgets never let NaN reach those fields. Documented inline. ~40 LOC of pure derives.
- Sketcher undo wired: `SketcherPanelState::history: History<Sketch>`, `record()` / `undo_edit` / `redo_edit` / `can_undo` / `can_redo` methods, `record()` calls inserted before every click-add-point, every constraint button (via `add_constraint_1` / `add_constraint_2` helpers), every Phase-12 primitive button (BSpline / Ellipse / EllipticalArc), Toggle Construction, every extra constraint button, every sketch op (Move / Rotate / Mirror / Copy / LinearArray / PolarArray). Inline `↶ ↷` row at top of the Sketcher panel. ~90 LOC.
- Part Design undo wired: same recipe — `PartDesignPanelState::history: History<FeatureTree>`, `record()` before Add Sketch / 16 Add-Feature buttons (Pad / Pocket / Revolve / Mirror / LinearPattern / CircularPattern / Fillet / Chamfer / Hole / Loft / Sweep / Pipe / Helix / MultiTransform / DraftAngle / Shell / Thickness / BooleanHistory) + Suppress + Delete + Imported / ImportedAdvanced from STEP/IGES. Inline `↶ ↷` row at top of the Part Design panel. ~120 LOC.
- Host Ctrl+Z / Ctrl+Y dispatcher (`update.rs::try_undo_in_active_panel` / `try_redo_in_active_panel`) extended to attempt Sketcher first, then Part Design when `show_mesh_toolbox` is on. ~15 LOC.
- Mesh Toolbox tail tooltips: every CAM op-kind variant (Pocket, Drill, Face, Adaptive Clearing, Helical Bore, Plunge Rough, Ramp Entry, Peck Drill Full, Contour 2D/3D, Engrave, Scribe, Spiral Pocket, Trochoidal Slot, Waterline 3D, Slot, Thread Mill, Rest Machining) — context-aware step-down / step-over / depth descriptions per kind. Surface tools: every NURBS curve / surface / Coons / Sew / Trim / KnotOps / SSI / Fit / Ruled input. TechDraw: sheet size + title block + view position/scale/parametric + linear-dim from/to/offset + dim chain kind/entries + balloon position/target/number/style + leader start/end/text/arrow + GD&T (existing helpers were already labelled). Assembly: part primitive picker + name + mate-kind ComboBox (per-kind hover descriptions) + joint-kind ComboBox (per-DOF hover descriptions) + Part A/B id inputs. Spreadsheet: sheet picker + new/remove sheet + view rows/cols + editor + Set/Clear/Re-evaluate. ~525 LOC of tooltips. Roughly 200+ individual hover-text annotations added in this pass.
