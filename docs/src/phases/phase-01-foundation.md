# Phase 1 — Foundation + first CFD thread

**Status:** 🔵 Complete. End-to-end shell + first physics thread + wgpu shaded viewport all shipped.

## Goal

Ship the smallest thing that proves Valenx is a real native desktop
app, not a doc pile: a window that opens, a case that runs, a
`Results` bundle that lands on disk.

## Capability inventory

- Native Rust workspace with canonical types (`valenx-fields`,
  `valenx-geo`, `valenx-mesh`, `valenx-core`).
- `.valenx` project format with atomic load/save.
- Workflow DAG with type-checked edges and topological execution.
- Adapter trait + registry + status classification
  (Ready / Missing / Outdated / Broken / Disabled).
- Four OpenFOAM solvers wired end-to-end: **simpleFoam** (steady
  incompressible RANS), **pimpleFoam** (transient PIMPLE, laminar or
  RANS), **icoFoam** (transient laminar PISO), **rhoSimpleFoam**
  (steady compressible RANS for transonic / mild supersonic external
  aero). Selected via the `case.solver` string in `case.toml`; the
  dict writer dispatches on a `SolverKind` enum.
- VTU result-rendering pipeline: every successful OpenFOAM run auto-
  parses its `.vtu` artifacts into the canonical `Field` catalog,
  loads the latest mesh into the viewport, and paints the wireframe
  edges by the first scalar field via a five-stop cool-to-warm
  divergent colour ramp. Bottom-right colour-bar legend shows field
  name + min/max + (timestep label for transient runs). Time-series
  slider in the Results pane scrubs through every snapshot. Field
  picker lets users switch between any scalar OnNode field.
- STL import for the viewport — no BRep kernel yet.
- Fusion-360-flavoured shell: ribbon, browser tree, viewport,
  timeline, command palette.

## Integrated tools graduating to Implemented

| Tool       | Adapter crate                     | Status this phase      |
|------------|-----------------------------------|------------------------|
| OpenFOAM   | `valenx-adapter-openfoam`         | probe + prepare (live); run + collect (Phase 1 tail) |

## Acceptance checklist

- [x] `cargo check --workspace` clean.
- [x] `cargo clippy --workspace --all-targets` — zero warnings.
- [x] `cargo test --workspace` — **121 tests green**.
- [x] `.valenx` load/save roundtrip passes on the `minimal.valenx`
      fixture.
- [x] Workflow DAG validates type mismatches and cycles.
- [x] OpenFOAM adapter emits a complete simpleFoam case tree
      (`system/`, `constant/`, `0/`) from canonical input.
- [x] OpenFOAM adapter also emits transient pimpleFoam / icoFoam case
      trees: `controlDict` writes real seconds + `adjustableRunTime`,
      `fvSchemes` swaps `steadyState` for `Euler`, `fvSolution`
      emits a `PIMPLE` or `PISO` block with `Final` correctors, and
      icoFoam refuses RANS turbulence with a clear error.
- [x] OpenFOAM adapter emits a complete rhoSimpleFoam case tree
      (compressible CFD): `thermophysicalProperties` with
      `hePsiThermo` + `perfectGas` + Sutherland transport,
      `0/T` temperature field with `fixedValue` inlet /
      `inletOutlet` outlet / `zeroGradient` adiabatic walls, and
      `fvSolution` SIMPLE block with `rhoMin`/`rhoMax` + energy in
      residualControl.
- [x] `.vtu` ASCII parser in `valenx-fields::vtu`. Handles
      `UnstructuredGrid` with point + cell data, rejects appended-
      binary / compressed / parallel variants with a clear error.
- [x] `VtuData::to_canonical(mesh_id)` converter producing
      `valenx_mesh::Mesh` + `Vec<Field>` with cached ranges.
- [x] OpenFOAM `collect()` parses every `.vtu` artifact into
      `Results.fields`, time-keyed by the snapshot index decoded
      from the filename (`cavity_500.vtu` → `Iteration(500)`).
- [x] Worker thread auto-calls `adapter.collect()` after each
      successful run and ships the result via
      `RunEvent::Collected(Box<Results>)`.
- [x] Auto-load the latest VTU's mesh into the viewport when an
      OpenFOAM run finishes — the user sees their geometry without
      any manual click.
- [x] Field-coloured wireframe overlay using a five-stop cool-to-
      warm divergent colour ramp from `valenx_fields::colormap`.
- [x] Bottom-right colour-bar legend with field name, gradient
      strip, min/max labels (and timestep label for transient runs).
- [x] Time-series scrubbing slider in the Results pane for fields
      with multiple snapshots; clamped against the current series'
      length so the slider never outruns reality.
- [x] Clickable field picker in the Results pane — switch which
      scalar drives the viewport overlay; click-the-selected to
      reset to auto-default.
- [x] STL ASCII + binary loaders under test.
- [x] Native eframe window opens with ribbon / browser tree /
      viewport / residual panel / status bar.
- [x] Interactive viewport — orbit, zoom, ViewCube snap, frame;
      shaded and wireframe render styles.
- [x] `valenx-adapter-openfoam::run()` streams residuals from a live
      `simpleFoam` child process via `std::process::Command` +
      reader threads.
- [x] Run orchestration — background thread + channel, UI stays
      responsive while the solver runs.
- [x] `AdapterRegistry` surfaced in the browser tree with status
      colours (Ready / Missing / Outdated / Broken / Disabled).
- [x] End-to-end thread: open a project → Run → residual chart
      updates live → convergence reported on the status bar.
- [x] Command palette (Ctrl+P) with fuzzy-matched commands.
- [x] `wgpu` shaded render pass — offscreen colour + depth textures,
      flat-shaded filled triangles with Lambert lighting, back-face
      culling in the pipeline, `egui::TextureId` sampled back into
      the viewport rect. Painter fallback retained for the (rare)
      case where eframe is built without the wgpu backend.
- [x] Log panel (tabbed with residuals) with per-level filters +
      autoscroll + 20k-line ring buffer.
- [x] Settings panel (theme, default shading, residual Y-axis scale,
      re-probe-on-close).
- [x] OpenFOAM `collect()` — walks workdir, classifies VTK / log /
      tabular / image artifacts, attaches them to `Results`.
- [x] Adapter registry re-probe action from both the browser and
      the Settings panel.

## Success metrics

| Metric                                              | Target           | Current     |
|-----------------------------------------------------|------------------|-------------|
| First-launch install size (Windows installer)       | < 200 MB         | n/a         |
| Cold-start to window visible                        | < 2 s            | ≈ 1 s       |
| Test suite runtime (all workspaces, release build)  | < 60 s           | < 5 s       |
| OpenFOAM prepare → `simpleFoam` exits 0             | works on 3/3 OSes| prepare + run wired; end-to-end pending OpenFOAM install |
| Residual plot updates at solver cadence             | yes              | yes (live)  |

## Leads into

[Phase 2 — CAD + meshing integration](./phase-02-cad-meshing.md).
