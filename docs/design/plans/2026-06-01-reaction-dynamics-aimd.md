# Reaction Dynamics Workbench (AIMD, Phase 1) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:executing-plans (inline) — the maintainer delegated the full build. Steps use checkbox (`- [ ]`) syntax.

**Goal:** A native 3-D ab-initio MD reaction simulator — small molecules react (bonds form/break) under quantum forces, computed in the background and played back in 3-D with an energy plot.

**Architecture:** New `valenx-reactdyn` crate (engine-agnostic: System/Controls/Frame/Trajectory + `ReactionEngine` trait + `AimdEngine`), forces by central finite-difference of `valenx-qchem`'s single-point energy, velocity-Verlet integrator, **all in atomic units** (qchem-native). A `valenx-app` workbench drives it and plays the trajectory back in 3-D, reusing `genetics::molecule_view` + `valenx_biostruct::geometry::vdw::is_bonded`.

**Tech stack:** Rust, valenx-qchem (energy + Element + geometry), egui/eframe, wgpu viewport.

**Units (locked):** integrator runs in a.u. — position bohr, energy hartree, force hartree/bohr, mass electron-masses, time a.u. Conversions at the boundary: `AMU_TO_AU_MASS = 1822.888486`, `FS_TO_AU_TIME = 41.341_374_575` (dt_au = dt_fs × this), `BOHR_PER_ANGSTROM = 1.889_726_124_625_770_2`.

---

### Task 0: Scaffold the `valenx-reactdyn` crate

**Files:** Create `crates/valenx-reactdyn/Cargo.toml`, `crates/valenx-reactdyn/src/lib.rs`, `crates/valenx-reactdyn/src/error.rs`; Modify root `Cargo.toml` (add to `[workspace] members` + a path dep is not needed yet).

- [ ] Cargo.toml: `name = "valenx-reactdyn"`, deps `valenx-qchem = { path = "../valenx-qchem" }`, `serde` (derive), `thiserror` (match siblings).
- [ ] `error.rs`: `#[derive(Debug, thiserror::Error)] pub enum ReactDynError { Qchem(String), Invalid{reason:String}, GuardExceeded{atoms:usize,steps:usize} }` with a `pub type Result<T> = std::result::Result<T, ReactDynError>;`.
- [ ] `lib.rs`: `#![forbid(unsafe_code)]` + `pub mod {units, integrator, forces, engine, error};` (add modules as tasks land) + re-exports.
- [ ] Add `"crates/valenx-reactdyn"` to root `Cargo.toml` workspace members (keep alphabetical-ish near valenx-qchem).
- [ ] **Verify:** `cargo build -p valenx-reactdyn` → green. **Commit:** `feat(reactdyn): scaffold valenx-reactdyn crate`.

### Task 1: Units module (TDD)

**Files:** Create `crates/valenx-reactdyn/src/units.rs`.

- [ ] **Test first:** `dt 1 fs == 41.3413… au`, `mass C 12.011 amu → 12.011×1822.888 au`, `1 Å == 1.8897 bohr`. Assert the three constants to 1e-6.
- [ ] Implement: `pub const AMU_TO_AU_MASS: f64 = 1822.888486;`, `pub const FS_TO_AU_TIME: f64 = 41.341_374_575;`, `pub const BOHR_PER_ANGSTROM: f64 = 1.889_726_124_625_770_2;` + `pub fn fs_to_au(fs)`, `pub fn amu_to_au_mass(amu)`.
- [ ] **Verify:** `cargo test -p valenx-reactdyn units` → green. **Commit.**

### Task 2: velocity-Verlet integrator (TDD, NO qchem)

**Files:** Create `crates/valenx-reactdyn/src/integrator.rs`.

- [ ] **Test first — harmonic oscillator vs analytic:** a 1-D mass-spring, `F(x) = -k·x`, run `velocity_verlet_step` for many steps via a closure force; assert (a) total energy `½mv² + ½kx²` conserved to < 0.5 % over ~one period, (b) the numerical period ≈ `2π√(m/k)` within a few %.
- [ ] Implement `pub fn velocity_verlet_step(pos: &mut [[f64;3]], vel: &mut [[f64;3]], masses: &[f64], dt: f64, forces_at: impl Fn(&[[f64;3]]) -> Vec<[f64;3]>) -> Vec<[f64;3]>` returning the new forces (so callers reuse them). Half-kick / drift / recompute-force / half-kick. `a = F/m`.
- [ ] `pub fn kinetic_energy(vel, masses) -> f64` = `Σ ½ mᵢ|vᵢ|²`.
- [ ] **Verify:** `cargo test -p valenx-reactdyn integrator` → green. **Commit:** `feat(reactdyn): velocity-Verlet integrator (energy-conserving)`.

### Task 3: Numerical forces from qchem (TDD)

**Files:** Create `crates/valenx-reactdyn/src/forces.rs`.

- [ ] **Test first — H₂ force sign:** at a compressed bond (atoms at ±0.25 bohr on z) the z-force must push the two apart; at a stretched bond (±1.0 bohr) it must pull together; near equilibrium (±0.70 bohr) |F| small. (RHF/STO-3G.)
- [ ] Define `pub enum Method { Rhf, Uhf, Dft }`. Implement `pub fn single_point_energy(elements: &[Element], pos_bohr: &[[f64;3]], charge, mult, method, basis) -> Result<f64>`: build `Atom::new(elem, pos)` (bohr) → `MolecularGeometry::with_charge_multiplicity` → `run_rhf/run_uhf/run_dft(&geom, basis, ScfSettings::default())` → `report.total_energy`; map SCF errors to `ReactDynError::Qchem`.
- [ ] Implement `pub fn numerical_forces(elements, pos_bohr, charge, mult, method, basis, delta_bohr) -> Result<Vec<[f64;3]>>`: central difference `F_iα = -(E(+δ) - E(-δ))/(2δ)` over every coordinate.
- [ ] **Verify:** `cargo test -p valenx-reactdyn forces` → green. **Commit:** `feat(reactdyn): numerical nuclear forces via finite-difference qchem energy`.

### Task 4: AIMD engine end-to-end (TDD)

**Files:** Create `crates/valenx-reactdyn/src/engine.rs`; Modify `lib.rs` re-exports.

- [ ] Types: `System { elements: Vec<Element>, pos_bohr: Vec<[f64;3]>, charge: i32, mult: u32 }`; `enum Thermostat { Nve, Berendsen { target_k: f64, tau_fs: f64 } }`; `Controls { method: Method, basis: String, dt_fs: f64, n_steps: usize, fd_delta_bohr: f64, thermostat: Thermostat, max_cost_guard: usize }` (+ `Default`); `Frame { time_fs, pos_bohr, vel, potential_hartree, kinetic_hartree }`; `Trajectory { system, frames }`; `trait ReactionEngine { fn run(&self, &System, &Controls, progress: &mut dyn FnMut(usize)) -> Result<Trajectory> }`; `struct AimdEngine`.
- [ ] **Test first — H₂ NVE conservation:** H₂ started ~0.1 Å off equilibrium, zero velocity, NVE, ~25 steps, RHF/STO-3G; assert total energy (KE+PE) drift < 2e-3 hartree across the run and `frames.len() == n_steps + 1`.
- [ ] Implement `AimdEngine::run`: guard `elements.len()*n_steps <= max_cost_guard`; masses via `Element::atomic_mass()*AMU_TO_AU_MASS`; `dt_au = fs_to_au(dt_fs)`; init velocities zero (NVE) or Maxwell-Boltzmann (Berendsen target); loop velocity-Verlet using `numerical_forces`; per step record a `Frame` (PE from the step's energy, KE from `kinetic_energy`); apply Berendsen velocity rescale when selected; `progress(step)`.
- [ ] **Verify:** `cargo test -p valenx-reactdyn` → all green; `cargo clippy -p valenx-reactdyn` clean. **Commit:** `feat(reactdyn): AIMD engine (velocity-Verlet + qchem forces, NVE/Berendsen)`.

### Task 5: Workbench wiring (mirror fem/cfd)

**Files:** Modify `crates/valenx-app/Cargo.toml` (+`valenx-reactdyn` dep), `crates/valenx-app/src/lib.rs` (`pub mod reactdyn_workbench;` + `show_reactdyn_workbench: bool` + `reactdyn: ReactdynWorkbenchState` fields + `enable_reactdyn_workbench`), `crates/valenx-app/src/update.rs` (View-menu checkbox + `draw_reactdyn_workbench` dispatch); Create a stub `crates/valenx-app/src/reactdyn_workbench.rs` (state + empty `draw_reactdyn_workbench`).

- [ ] Mirror the exact pattern from `fem_workbench`/`cfd_workbench` (already in the tree). Stub `ReactdynWorkbenchState: Default` + a `draw_reactdyn_workbench(app, ctx)` that shows a titled SidePanel.
- [ ] **Verify:** `cargo build -p valenx-app` green. **Commit:** `feat(app): wire Reaction Dynamics workbench (stub panel)`.

### Task 6: Workbench setup + background run (TDD on the run helper)

**Files:** Modify `crates/valenx-app/src/reactdyn_workbench.rs`.

- [ ] State: preset selector + XYZ text, charge/mult, method/basis, dt_fs/n_steps, thermostat; `Option<Trajectory>` result + `error`; a background `JoinHandle`/channel like aero.
- [ ] **Test first (headless):** a `run_reactdyn(state)`-style helper on the default preset (H₂) with small `n_steps` returns a `Trajectory` with frames and `error.is_none()`; bad XYZ → `error.is_some()` (fail loud).
- [ ] Presets: H₂ dissociation, water, a proton transfer (XYZ literals). XYZ parse via the existing qchem `MolecularGeometry::from_xyz_str` (reuse) → System.
- [ ] **Verify:** `cargo test -p valenx-app reactdyn` green. **Commit:** `feat(app): Reaction Dynamics setup + background AIMD run`.

### Task 7: 3-D trajectory playback + energy plot

**Files:** Modify `crates/valenx-app/src/reactdyn_workbench.rs`.

- [ ] Timeline scrubber (play/pause/step/scrub `frame_idx`).
- [ ] Per displayed frame: build a `molecule_view::ViewMolecule` (atoms at `frame.pos` in Å, **bonds recomputed** via `valenx_biostruct::geometry::vdw::is_bonded(elem_a, elem_b, dist_Å, tol)`) → ball-and-stick mesh → push to the viewport (reuse the genetics molecule→viewport path).
- [ ] Energy plot: potential / kinetic / total vs time (egui plot or a simple line widget already used elsewhere).
- [ ] **Verify:** `cargo build -p valenx-app` green; `cargo clippy -p valenx-app --tests` clean for the new file; run the app, toggle View → Reaction Dynamics, run the H₂ preset, watch the bond stretch/contract + energy line. **Commit:** `feat(app): 3-D AIMD trajectory playback + energy plot`.

---

## Self-review

- **Spec coverage:** crate shell (T0,4) ✓; numerical FD forces (T3) ✓; velocity-Verlet (T2) ✓; NVE+Berendsen (T4) ✓; tiny-system guard (T4) ✓; workbench + presets/XYZ + background run (T5,6) ✓; 3-D dynamic-bond playback + energy plot (T7) ✓; fail-loud (T3,4,6) ✓; tests incl. harmonic-oscillator + force-sign + conservation (T2,3,4) ✓. Phases 2/3 (QM/MM, ReaxFF) intentionally deferred behind the `ReactionEngine` trait.
- **Placeholders:** none — concrete APIs (`Atom::new`, `run_rhf`, `Element::atomic_mass`, `vdw::is_bonded`, `MolecularGeometry::from_xyz_str`) and constants are named.
- **Type consistency:** `Method`, `System`, `Controls`, `Frame`, `Trajectory`, `ReactionEngine`, `AimdEngine` used consistently across T3–T6; units a.u. throughout the engine, Å only at the qchem-XYZ and viewport boundaries.
