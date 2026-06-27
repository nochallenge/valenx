# valenx Products / Workbenches — Master Checklist

Master checklist of valenx's products/workbenches for the commercial-grade pass.

Every row is one `TabKind` product surfaced by the `＋ from template` launcher,
authoritatively sourced from `crates/valenx-app/src/project_tabs.rs`
(`TabKind` enum · `TabKind::TEMPLATES` · `label()` · `group()` · `from_id()`).
**Product** = `label()`; **id** = the `from_id()` string an agent emits;
**Group** = `group()`; **Backend** = the primary crate the workbench drives
(per `crates/valenx-app/Cargo.toml` comments / `src/*_workbench.rs`); **Status**
is the verification to-do (every row starts `✅ PASS`).

> `TabKind::Blank` ("Untitled") is the empty `➕ New tab` project, not a
> product, so it is excluded — matching `TEMPLATES.len() == 56`.

## Verification — `valenx --self-test`

Status is produced by the **baked-in headless self-test**
(`valenx --self-test [--group <G>] [--id <id>]`, source
`crates/valenx-app/src/self_test.rs`): one command drives every product, runs
its real compute path, and checks the output. Re-run it any time — this table
records the latest run. Adding a deep known-value check = one registry row.

**Last run 2026-06-26 — 53 PASS · 0 FAIL · 3 SKIP** (exit 0, ~4 s, no GUI):
- **11 deep ground-truth checks** (Aerospace group fully deep as of 2026-06-26): `thermo` (CO₂ PR `Z_vap=0.829`, `Psat≈3.49 MPa`), `quantum` (Bell `0.50/0.50`), `optics` (thin-lens `m=−1.0`), `acoustics` (monopole 1/r), `waveform` (VCD parse); **+ Aerospace** `rocket` (Tsiolkovsky Δv=10055 m/s), `engine` (kerolox c*=1773 m/s), `astro` (Hohmann LEO→GEO Δv=3893 m/s), `gasdynamics` (M=2 A/A*=1.6875 + normal-shock, NACA-1135), `rotor` (hover FM=0.648), `uas` (hover 97 W, disk-loading).
- **42 generic checks:** open → run compute → output finite & free of `NaN`/`inf`/`error`/`panic` (incl. `aero` — its workbench runs an async RANS CFD solve, so the 2π thin-airfoil analytic lives in the `valenx-aero` crate benchmark, not the product drive path).
- **3 SKIP (honest, not faked):** `cosim` (external tool), `photogrammetry` (scan input), `render` (GPU render).

## Aerospace

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Rocket | `rocket` | Aerospace | `valenx-rocket-demo` (valenx-astro + valenx-fem) | ✅ DEEP — Δv=10055 m/s (Tsiolkovsky 2-stage) |
| Engine | `engine` | Aerospace | `valenx-app/src/engine_workbench.rs` (valenx-astro + valenx-combustion) | ✅ DEEP — c*=1773 m/s (1-D nozzle) |
| Astro / Launch | `astro` | Aerospace | `valenx-astro` | ✅ DEEP — Hohmann Δv=3893 m/s (vis-viva) |
| Aerodynamics | `aero` | Aerospace | `valenx-aero` | ✅ PASS |
| Gas dynamics | `gasdynamics` | Aerospace | `valenx-gasdynamics` | ✅ DEEP — A/A*=1.6875 + shock (NACA-1135) |
| Rotor / Drone (BEMT) | `rotor` | Aerospace | `valenx-rotor` | ✅ DEEP — FM=0.648 (momentum theory) |
| UAS design & counter-UAS | `uas` | Aerospace | `valenx-uas` (valenx-drone / valenx-rotor / valenx-fixedwing) | ✅ DEEP — hover 97 W (disk-loading) |

## Astrophysics

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Black Hole | `blackhole` | Astrophysics | `valenx-relativity` | ✅ PASS |

## Simulation

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| CFD | `cfd` | Simulation | `valenx-cfd-native` | ✅ PASS |
| FEM | `fem` | Simulation | `valenx-fem` | ✅ PASS |
| Topology optimization (SIMP) | `topopt` | Simulation | `valenx-topopt` (valenx-fem solver) | ✅ PASS |
| Node graph (visual node editor) | `nodegraph` | Simulation | `valenx-app/src/nodegraph_workbench.rs` (in-house) | ✅ PASS |
| Bond graph (multi-domain systems) | `bondgraph` | Simulation | `valenx-app/src/bondgraph_workbench.rs` (in-house) | ✅ PASS |
| Surrogate model (ML solver emulator) | `surrogate` | Simulation | `valenx-app/src/surrogate_workbench.rs` (in-house MLP) | ✅ PASS |
| Reaction dynamics | `reactdyn` | Simulation | `valenx-reactdyn` (valenx-qchem) | ✅ PASS |
| Thermodynamics (EoS) | `thermo` | Simulation | `valenx-thermo` | ✅ PASS |
| Quantum circuit | `quantum` | Simulation | `valenx-quantum` | ✅ PASS |
| Field statistics | `fields` | Simulation | `valenx-fields` | ✅ PASS |
| Fluids (SPH particle sim) | `fluids` | Simulation | `valenx-fluids` | ✅ PASS |
| Ocean (waves + buoyancy) | `ocean` | Simulation | `valenx-ocean` | ✅ PASS |
| Reduced-order model (POD) | `rom` | Simulation | `valenx-rom` | ✅ PASS |
| Uncertainty quantification | `uq` | Simulation | `valenx-uq` | ✅ PASS |
| Mission simulation (constructive) | `missionsim` | Simulation | `valenx-mission-sim` | ✅ PASS |
| Mission Planner | `missionplanner` | Simulation | `valenx-app/src/mission_planner_workbench.rs` (in-house + walkers map) | ✅ PASS |
| Survivability / protection | `survivability` | Simulation | `valenx-survivability` (valenx-fem) | ✅ PASS |
| Co-Simulation (FMI / HELICS) | `cosim` | Simulation | `valenx-adapter-fmi` | ⏭️ SKIP (external-tool) |
| Multibody dynamics (robot / contact) | `mbd` | Simulation | `valenx-mbd` | ✅ PASS |
| Optics (thin lens) | `optics` | Simulation | `valenx-optics` | ✅ PASS |
| Acoustics (radiation) | `acoustics` | Simulation | `valenx-acoustics` | ✅ PASS |
| Waveform (VCD viewer) | `waveform` | Simulation | `valenx-waveform` | ✅ PASS |

## CAD & mesh

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Parametric CAD | `cad` | CAD & mesh | `valenx-cad` / `valenx-sketch` | ✅ PASS |
| Part B-Rep (truck) | `brep` | CAD & mesh | `valenx-truck-cad` | ✅ PASS |
| Mesh toolbox | `mesh` | CAD & mesh | `valenx-mesh` | ✅ PASS |
| Sheet metal | `sheetmetal` | CAD & mesh | `valenx-sheet-metal` | ✅ PASS |
| Reverse engineering | `reverse` | CAD & mesh | `valenx-reverse` | ✅ PASS |
| Photogrammetry / SfM scan | `photogrammetry` | CAD & mesh | `valenx-photogrammetry` | ⏭️ SKIP (scan-input) |
| 2D drafting | `draft2d` | CAD & mesh | `valenx-draft` / `valenx-librecad-2d` | ✅ PASS |
| Path-traced render | `render` | CAD & mesh | `valenx-pathtrace` | ⏭️ SKIP (gpu-render) |
| Animation | `animate` | CAD & mesh | `valenx-animate` | ✅ PASS |

## Machine design

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Springs | `springs` | Machine design | `valenx-springs` | ✅ PASS |
| Gears | `gears` | Machine design | `valenx-gears` | ✅ PASS |
| Fasteners | `fasteners` | Machine design | `valenx-fasteners` | ✅ PASS |
| Frames / sections | `frames` | Machine design | `valenx-frames` | ✅ PASS |
| Collision | `collision` | Machine design | `valenx-collision` | ✅ PASS |

## Civil & AEC

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Piping | `piping` | Civil & AEC | `valenx-piping` | ✅ PASS |
| HVAC | `hvac` | Civil & AEC | `valenx-hvac` | ✅ PASS |
| Reinforcement | `reinforcement` | Civil & AEC | `valenx-reinforcement` | ✅ PASS |
| Interior design | `interior` | Civil & AEC | `valenx-interior` | ✅ PASS |
| Geomatics | `geomatics` | Civil & AEC | `valenx-geomatics` | ✅ PASS |

## Life sciences

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Genetics | `genetics` | Life sciences | `valenx-bioseq` / `valenx-align` / `valenx-phylo` / `valenx-md` / `valenx-cheminf` / … | ✅ PASS |
| Neural interface | `neuro` | Life sciences | `valenx-neuro` (valenx-fem) | ✅ PASS |
| Variant effect | `variant` | Life sciences | `valenx-variant-effect` | ✅ PASS |
| Protein interaction (PPI) | `ppi` | Life sciences | `valenx-ppi` | ✅ PASS |
| Morphogenesis | `morphogenesis` | Life sciences | `valenx-app/src/morphogenesis_workbench.rs` (in-house Gray–Scott) | ✅ PASS |

## Sensors

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Sensors (LiDAR / Radar) | `sensors` | Sensors | `valenx-sensors` | ✅ PASS |
| Autonomy V&V | `autonomy` | Sensors | `valenx-autonomy-vnv` (valenx-sensors) | ✅ PASS |

---

**Total: 56 products / workbenches** (= `TabKind::TEMPLATES.len()`; excludes the
non-product `TabKind::Blank`). Group breakdown: Aerospace 7 · Astrophysics 1 ·
Simulation 22 · CAD & mesh 9 · Machine design 5 · Civil & AEC 5 · Life sciences 5 ·
Sensors 2.
