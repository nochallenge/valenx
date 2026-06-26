# valenx Products / Workbenches — Master Checklist

Master checklist of valenx's products/workbenches for the commercial-grade pass.

Every row is one `TabKind` product surfaced by the `＋ from template` launcher,
authoritatively sourced from `crates/valenx-app/src/project_tabs.rs`
(`TabKind` enum · `TabKind::TEMPLATES` · `label()` · `group()` · `from_id()`).
**Product** = `label()`; **id** = the `from_id()` string an agent emits;
**Group** = `group()`; **Backend** = the primary crate the workbench drives
(per `crates/valenx-app/Cargo.toml` comments / `src/*_workbench.rs`); **Status**
is the verification to-do (every row starts `⬜ unverified`).

> `TabKind::Blank` ("Untitled") is the empty `➕ New tab` project, not a
> product, so it is excluded — matching `TEMPLATES.len() == 56`.

## Aerospace

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Rocket | `rocket` | Aerospace | `valenx-rocket-demo` (valenx-astro + valenx-fem) | ⬜ unverified |
| Engine | `engine` | Aerospace | `valenx-app/src/engine_workbench.rs` (valenx-astro + valenx-combustion) | ⬜ unverified |
| Astro / Launch | `astro` | Aerospace | `valenx-astro` | ⬜ unverified |
| Aerodynamics | `aero` | Aerospace | `valenx-aero` | ⬜ unverified |
| Gas dynamics | `gasdynamics` | Aerospace | `valenx-gasdynamics` | ⬜ unverified |
| Rotor / Drone (BEMT) | `rotor` | Aerospace | `valenx-rotor` | ⬜ unverified |
| UAS design & counter-UAS | `uas` | Aerospace | `valenx-uas` (valenx-drone / valenx-rotor / valenx-fixedwing) | ⬜ unverified |

## Astrophysics

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Black Hole | `blackhole` | Astrophysics | `valenx-relativity` | ⬜ unverified |

## Simulation

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| CFD | `cfd` | Simulation | `valenx-cfd-native` | ⬜ unverified |
| FEM | `fem` | Simulation | `valenx-fem` | ⬜ unverified |
| Topology optimization (SIMP) | `topopt` | Simulation | `valenx-topopt` (valenx-fem solver) | ⬜ unverified |
| Node graph (visual node editor) | `nodegraph` | Simulation | `valenx-app/src/nodegraph_workbench.rs` (in-house) | ⬜ unverified |
| Bond graph (multi-domain systems) | `bondgraph` | Simulation | `valenx-app/src/bondgraph_workbench.rs` (in-house) | ⬜ unverified |
| Surrogate model (ML solver emulator) | `surrogate` | Simulation | `valenx-app/src/surrogate_workbench.rs` (in-house MLP) | ⬜ unverified |
| Reaction dynamics | `reactdyn` | Simulation | `valenx-reactdyn` (valenx-qchem) | ⬜ unverified |
| Thermodynamics (EoS) | `thermo` | Simulation | `valenx-thermo` | ⬜ unverified |
| Quantum circuit | `quantum` | Simulation | `valenx-quantum` | ⬜ unverified |
| Field statistics | `fields` | Simulation | `valenx-fields` | ⬜ unverified |
| Fluids (SPH particle sim) | `fluids` | Simulation | `valenx-fluids` | ⬜ unverified |
| Ocean (waves + buoyancy) | `ocean` | Simulation | `valenx-ocean` | ⬜ unverified |
| Reduced-order model (POD) | `rom` | Simulation | `valenx-rom` | ⬜ unverified |
| Uncertainty quantification | `uq` | Simulation | `valenx-uq` | ⬜ unverified |
| Mission simulation (constructive) | `missionsim` | Simulation | `valenx-mission-sim` | ⬜ unverified |
| Mission Planner | `missionplanner` | Simulation | `valenx-app/src/mission_planner_workbench.rs` (in-house + walkers map) | ⬜ unverified |
| Survivability / protection | `survivability` | Simulation | `valenx-survivability` (valenx-fem) | ⬜ unverified |
| Co-Simulation (FMI / HELICS) | `cosim` | Simulation | `valenx-adapter-fmi` | ⬜ unverified |
| Multibody dynamics (robot / contact) | `mbd` | Simulation | `valenx-mbd` | ⬜ unverified |
| Optics (thin lens) | `optics` | Simulation | `valenx-optics` | ⬜ unverified |
| Acoustics (radiation) | `acoustics` | Simulation | `valenx-acoustics` | ⬜ unverified |
| Waveform (VCD viewer) | `waveform` | Simulation | `valenx-waveform` | ⬜ unverified |

## CAD & mesh

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Parametric CAD | `cad` | CAD & mesh | `valenx-cad` / `valenx-sketch` | ⬜ unverified |
| Part B-Rep (truck) | `brep` | CAD & mesh | `valenx-truck-cad` | ⬜ unverified |
| Mesh toolbox | `mesh` | CAD & mesh | `valenx-mesh` | ⬜ unverified |
| Sheet metal | `sheetmetal` | CAD & mesh | `valenx-sheet-metal` | ⬜ unverified |
| Reverse engineering | `reverse` | CAD & mesh | `valenx-reverse` | ⬜ unverified |
| Photogrammetry / SfM scan | `photogrammetry` | CAD & mesh | `valenx-photogrammetry` | ⬜ unverified |
| 2D drafting | `draft2d` | CAD & mesh | `valenx-draft` / `valenx-librecad-2d` | ⬜ unverified |
| Path-traced render | `render` | CAD & mesh | `valenx-pathtrace` | ⬜ unverified |
| Animation | `animate` | CAD & mesh | `valenx-animate` | ⬜ unverified |

## Machine design

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Springs | `springs` | Machine design | `valenx-springs` | ⬜ unverified |
| Gears | `gears` | Machine design | `valenx-gears` | ⬜ unverified |
| Fasteners | `fasteners` | Machine design | `valenx-fasteners` | ⬜ unverified |
| Frames / sections | `frames` | Machine design | `valenx-frames` | ⬜ unverified |
| Collision | `collision` | Machine design | `valenx-collision` | ⬜ unverified |

## Civil & AEC

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Piping | `piping` | Civil & AEC | `valenx-piping` | ⬜ unverified |
| HVAC | `hvac` | Civil & AEC | `valenx-hvac` | ⬜ unverified |
| Reinforcement | `reinforcement` | Civil & AEC | `valenx-reinforcement` | ⬜ unverified |
| Interior design | `interior` | Civil & AEC | `valenx-interior` | ⬜ unverified |
| Geomatics | `geomatics` | Civil & AEC | `valenx-geomatics` | ⬜ unverified |

## Life sciences

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Genetics | `genetics` | Life sciences | `valenx-bioseq` / `valenx-align` / `valenx-phylo` / `valenx-md` / `valenx-cheminf` / … | ⬜ unverified |
| Neural interface | `neuro` | Life sciences | `valenx-neuro` (valenx-fem) | ⬜ unverified |
| Variant effect | `variant` | Life sciences | `valenx-variant-effect` | ⬜ unverified |
| Protein interaction (PPI) | `ppi` | Life sciences | `valenx-ppi` | ⬜ unverified |
| Morphogenesis | `morphogenesis` | Life sciences | `valenx-app/src/morphogenesis_workbench.rs` (in-house Gray–Scott) | ⬜ unverified |

## Sensors

| Product | id | Group | Backend | Status |
|---|---|---|---|---|
| Sensors (LiDAR / Radar) | `sensors` | Sensors | `valenx-sensors` | ⬜ unverified |
| Autonomy V&V | `autonomy` | Sensors | `valenx-autonomy-vnv` (valenx-sensors) | ⬜ unverified |

---

**Total: 56 products / workbenches** (= `TabKind::TEMPLATES.len()`; excludes the
non-product `TabKind::Blank`). Group breakdown: Aerospace 7 · Astrophysics 1 ·
Simulation 22 · CAD & mesh 9 · Machine design 5 · Civil & AEC 5 · Life sciences 5 ·
Sensors 2.
