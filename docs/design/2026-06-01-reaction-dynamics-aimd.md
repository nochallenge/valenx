# Reaction Dynamics Workbench — AIMD (Phase 1) Design

**Status:** approved design, pre-implementation
**Date:** 2026-06-01
**Author:** brainstormed with the maintainer; review delegated to the implementer.

## Goal

A native 3-D **reaction simulator**: watch small molecules actually *react* —
bonds forming and breaking — driven by first-principles quantum forces, computed
in the background and played back in 3-D with an energy plot that proves the
physics is honest.

This is **Phase 1** of a multi-engine program. The 3-D simulator *shell* is
engine-agnostic; Phase 1 ships the shell plus the **ab-initio MD (AIMD)** backend.
Later phases add **QM/MM** (a molecule reacting in a solvent / pocket / surface)
and **ReaxFF** (reactive force field, materials scale) behind the same trait, with
no shell rework.

## Scope (honest boundaries)

- **Physics:** Born–Oppenheimer AIMD. At each timestep the electronic energy is
  solved by `valenx-qchem` (RHF/UHF/DFT) and the nuclei are advanced by a
  velocity-Verlet integrator.
- **Forces:** **numerical** — central finite differences of the single-point
  energy (see "Force model" below). qchem has *no* analytic nuclear gradient
  (`GeometryOptRequest::run` is an honest `NotYetImplemented` stub), so AIMD must
  not pretend otherwise.
- **Size:** small systems (~≤ 12 atoms in v1), HF or small-basis DFT, tens-to-
  hundreds of steps.
- **Mode:** compute the trajectory on a background thread → **play it back in 3-D**.
  Quantum forces are far too slow to step live; this matches production AIMD
  tooling (VMD-style playback).
- **Environment:** vacuum. (Solvent/pocket = Phase 2 QM/MM. Materials scale =
  Phase 3 ReaxFF.)

## Architecture

A new engine-agnostic crate plus a desktop workbench.

### `valenx-reactdyn` (new crate)

```
struct System    { atoms: Vec<(Element, [f64; 3])>, charge: i32, multiplicity: u32 }
struct Controls  { method: Method, basis: String, dt_fs: f64, n_steps: usize,
                   thermostat: Thermostat, fd_delta: f64, max_cost_guard: usize }
enum   Thermostat { Nve, Berendsen { target_k: f64, tau_fs: f64 } }
struct Frame     { time_fs: f64, positions: Vec<[f64;3]>, velocities: Vec<[f64;3]>,
                   potential_hartree: f64, kinetic_hartree: f64 }
struct Trajectory{ system: System, frames: Vec<Frame> }

trait ReactionEngine {
    fn run(&self, system: &System, controls: &Controls,
           progress: &mut dyn FnMut(usize)) -> Result<Trajectory, ReactDynError>;
}

struct AimdEngine;   // Phase 1. (Phase 2: QmMmEngine; Phase 3: ReaxffEngine.)
```

- **Integrator:** velocity-Verlet (time-reversible, energy-conserving), the
  standard MD integrator.
- **Force model:** `numerical_forces(system, method, basis, delta)` returns
  `Vec<[f64;3]>` by central finite difference of the qchem single-point energy:
  `F_iα = −[E(r + δ·e_iα) − E(r − δ·e_iα)] / (2δ)`. Cost = `6N` single-point
  energies per step; a `max_cost_guard` (atoms × steps) bounds runaway compute.
- **Masses:** atomic masses from the element (reuse the existing periodic-table
  source used by qchem/md).
- **Thermostat:** NVE by default (no coupling); optional Berendsen velocity
  rescaling (reuse `valenx-md`'s thermostat math if exposed, else a few lines).
- **Initial velocities:** zero (NVE) or Maxwell–Boltzmann at `target_k` (NVT).
- **Errors:** SCF non-convergence, bad geometry, or guard-exceeded → a typed
  `ReactDynError`, surfaced loud. Never a fabricated frame.

### `valenx-app/src/reactdyn_workbench.rs` (new workbench)

Mirrors the FEM/CFD/aero workbenches: a resizable right `SidePanel` gated on
`show_reactdyn_workbench`, toggled from the View menu.

- **Setup:** preset-reaction library (e.g. H₂ dissociation, water formation, a
  proton transfer) + XYZ paste; charge/multiplicity; method + basis; `dt`, steps;
  thermostat.
- **Run:** background thread (the aero pattern) producing a `Trajectory`, with a
  progress bar.
- **Playback:** timeline scrubber (play / pause / step / scrub) + an energy plot
  (potential / kinetic / total vs time). For NVE, a flat total-energy line *is*
  the built-in correctness check.

### 3-D rendering (reuse what exists)

Per displayed frame, build a `ViewMolecule` (atoms at `frame.positions`, **bonds
recomputed from covalent-radius distance cutoffs**) and feed
`genetics::molecule_view`'s ball-and-stick mesh builder to the wgpu viewport. The
per-frame bond recompute is what makes bonds visibly appear/disappear — the whole
point of a *reaction* sim.

## Data flow

setup → background AIMD run (`6N` qchem energies × `n_steps`) → `Trajectory` →
3-D playback + energy plot.

## Reuse vs. new

- **Reuse:** qchem single-point energy (RHF/UHF/DFT), `valenx-md` (thermostat +
  atomic masses), `genetics::molecule_view` (3-D mesh), the wgpu viewport, the
  workbench + background-run + View-menu patterns.
- **New:** the `valenx-reactdyn` crate (trait, velocity-Verlet, numerical forces,
  trajectory) + the workbench UI + dynamic-bond playback + the energy plot.

## Testing

- **Integrator (no qchem):** velocity-Verlet on a 1-D harmonic oscillator vs. the
  analytic solution — energy conservation + period. De-risks the integrator alone.
- **Forces:** `numerical_forces` on H₂ — at a compressed bond the force pushes the
  atoms apart, at a stretched bond it pulls them together, ≈ 0 at equilibrium.
- **AIMD end-to-end:** H₂ NVE for a few dozen steps → total energy conserved
  within tolerance; trajectory has the expected length.
- **Bonds:** distance-cutoff bond recompute unit tests (bond present when close,
  absent when far).
- **Fail-loud:** SCF non-convergence / bad input / guard-exceeded → error, no
  panic, no fabricated frame.

## Honesty

- Numerical (not analytic) gradients — stated in the UI and the docs.
- Tiny systems only; `max_cost_guard` prevents runaway compute.
- The energy-conservation plot is the live correctness check the user can see.

## Out of scope for v1 (roadmap)

Analytic gradients; environment/solvent (Phase 2 QM/MM); materials scale
(Phase 3 ReaxFF); enhanced sampling / free energy; ML interatomic potentials;
real-time (live) stepping.
```
