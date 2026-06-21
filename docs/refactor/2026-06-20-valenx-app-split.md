# valenx-app crate split — staged execution plan (2026-06-20)

Goal: editing one workbench should recompile a small leaf crate, not all ~115k
lines of `valenx-app`. Driven because incremental rebuilds are slow.

## Confirmed architecture (from the read-only analysis)
- Workbenches are **independent**: zero cross-workbench `use`/calls (only doc-link
  comments). Clean domain boundaries.
- Shared surface a workbench body touches beyond its own state: `app.mesh`,
  `app.stl`, `app.frame_current_mesh()`, `app.frame_current_stl()`,
  `app.aero_field_overlay` (fem/cad only) — and its own `show_<x>_workbench` +
  `<x>` state field. That's it.
- **Cycle-break (the only sound shape):** `ValenxApp` + ALL `*WorkbenchState`
  structs stay in `valenx-app-core`. Only the **draw logic** (`draw_*` / `*_body`
  / compute / mesh-emit fns + their `#[cfg(test)]`) moves to leaves.
  `frame_current_*` are inherent `impl ValenxApp` → must stay in core (orphan rule).
- **Dispatch must live in the BINARY** (`valenx-app`), not core: `update.rs`'s
  `impl eframe::App::update` (the 132 `draw_X` calls + View menu) and
  `dock_layout::render_panel_body` + `drain_workbench_deferred` call leaves →
  if they were in core, core→leaf→core cycle. Binary depends on core + all leaves.
- `TabKind`/`project_tabs` only flips `show_*` flags → stays in core.
- Stays in CORE as infra: `workbench_chrome` (keystone `workbench_shell`),
  `menu_ui`, `types`, `viewport*`, `wgpu_renderer`, `pbr_forward_pass`,
  `mesh_loader`, `mesh_toolbox` (default panel + impl ValenxApp), overlays,
  `assistant_workbench` (per-unit chat plumbing), `dock_layout` (state +
  `draw_dock_layout`; the `render_panel_body` MATCH goes to binary), settings,
  theme, landing, file_browser, run/sweep/audit (impl ValenxApp), setup, headless,
  rocket_mesh, state_paths, etc.

## 7 leaf crates (draw bodies grouped by domain; struct defs go to core)
- **valenx-app-aerospace**: aero, astro, rocket, engine, gasdynamics, combustion, fixedwing, drone, car, marine, rail, windturbine, solarpv, antenna (+ aero/ astro/ subdirs)
- **valenx-app-simulation**: cfd, fem, reactdyn, fields, fft, vibration, acoustics, optics, diffusion, dimensional, fanlaws, queueing, popdynamics, radioactivity, collision, projectile, statics, fluidstatics, animate, render
- **valenx-app-cad-mesh**: cad, sheetmetal, reverse, draft2d, interior, geomatics, frames, plate (mesh_toolbox stays in core)
- **valenx-app-mechanical**: springs, gears, fasteners, bearing, beltdrive, brake, buckling, fatigue, geartooth, gearbox, chaindrive, clutch, coil, bolt, rivet, shaftdesign, screwthread, pulley, springdesign, springcombination, flywheel, fourbar, camdynamics, leadscrew, leverage, inclinedplane, mohr, torsion, fracture, thermalexpansion, pressurevessel, creep, strainrosette, straingauge
- **valenx-app-electrical**: dcmotor, inductionmotor, threephase, transformer, powerfactor, capacitor, mosfet, bjt, opamp, led, rectifier, filter, resistornetwork, transmissionline, batteryecm, batterypack, thermistor, thermocouple
- **valenx-app-civil-aec**: piping, hvac, reinforcement, rcbeam, columnsteel, soilbearing, retainingwall, beam, truss, weir, openchannel, pipeflow, pipenetwork, hydraulics, pneumatics, orifice, pump, heatexchanger, heattransfer, heatpump, refrigeration, insulation, psychrometrics, thermocycle
- **valenx-app-life-sciences**: genetics (+genetics/ subdir), neuro, variant_effect, bonemech, hemodynamics, thermoreg, pharmacokinetics, enzymekinetics, acidbase, bmr, osmosis, electrochem

## Stages (build GREEN after every step; gate `cargo build -p valenx-app` + `cargo test -p <crate> --lib headless_ui_tests`; NEVER `cargo test --workspace`)
- **Stage 0** (in flight): widen `ValenxApp` struct fields `pub(crate)`→`pub`. No-op, validates.
- **Stage A** (the bulk, sub-stage it): create `crates/valenx-app-core`; move the core file set there; binary keeps all `*_workbench.rs` + `main.rs`/`setup` glue + `update.rs`'s `update()`+dispatch + `render_panel_body`. Split `update.rs` (impl-ValenxApp helpers→core, `update()`+dispatch→binary) and `dock_layout` (state/`draw_dock_layout`→core, `render_panel_body` match→binary). Rewrite `crate::`→`valenx_app_core::` across files. Win here = 0 (workbenches still in binary); this only sets the boundary.
- **Stage B** (proof + MEASURE + CHECKPOINT): extract ONE small clean leaf first (life-sciences or a small cluster). Per workbench: struct→core (`core/src/<x>_state.rs`), draw body+helpers+tests→leaf. Binary dispatch calls `valenx_app_<domain>::draw_X`. MEASURE: `touch` a leaf file → `cargo build -p valenx-app` vs the old whole-crate `touch`. **Report the number to the user before rolling out the rest.**
- **Stage C**: roll out remaining leaves one at a time, gate green + commit each. Order: life-sciences, aerospace, simulation, civil-aec, cad-mesh, mechanical, electrical.

## File-split mechanic
Each `<x>_workbench.rs` today = `struct XWorkbenchState{}+Default` + `draw_X`/`X_body`/helpers/tests. Split: struct+Default(+struct tests)→core; the rest→leaf. Struct blocks are short (~30 lines); the bulk (draw body) moves wholesale.

## Risks (ranked)
1. Splitting `update.rs` (4392 lines) — relocating `eframe::App::update` dispatch to the binary while core exposes its helpers `pub`. Do Stage 0 first (done) to de-risk visibility.
2. The ~131 per-file struct/body splits — high count, low individual risk; gate per leaf.
3. `drain_deferred` (rocket/engine/fem) + `aero_field_overlay` — the only non-trivial shared mutable surface; handled by dispatch-in-binary.

Commit each green stage LOCALLY (email-safe, leak==0); GitHub deferred.
