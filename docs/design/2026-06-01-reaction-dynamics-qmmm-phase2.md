# Reaction Dynamics — Phase 2: QM/MM (design + plan)

**Status:** approved, pre-implementation
**Date:** 2026-06-01
**Builds on:** Phase 1 (AIMD) — `valenx-reactdyn` + the Reaction Dynamics workbench.

## Goal

A `QmMmEngine` (behind the existing `ReactionEngine` trait) that simulates a
**quantum reacting region embedded in an explicit classical environment** — a
molecule reacting in solvent / a pocket — and plays it back in 3-D with the
environment visible. Mechanical embedding ships first (working pipeline);
electrostatic embedding is the in-scope accuracy upgrade.

## Architecture (all in `valenx-reactdyn`, behind `ReactionEngine`)

- **QM region:** elements + positions (bohr), charge, multiplicity — forces via
  Phase 1's `forces::numerical_forces` (qchem).
- **MM region:** a `valenx_md::System` + `ForceField` (OPLS-AA) — energy + analytic
  forces via `valenx_md::sim::Simulation::evaluate_forces`.
- **QM–MM coupling (mechanical):** classical Lennard-Jones + Coulomb between every
  QM atom (carrying LJ σ/ε + a fixed charge) and every MM atom. Forces on both.
- **Integrator:** one velocity-Verlet (Phase 1's, reused) over *all* atoms.
- **Units:** the integrator runs in **atomic units**; valenx-md is **GROMACS**
  (nm / ps / u / e / kJ·mol⁻¹). A conversion layer reconciles them — the classic
  QM/MM footgun, so it is built + unit-tested first.

## Honest v1 caveats

- **Mechanical embedding** — the environment does not polarize the QM density.
  Electrostatic embedding (below) is the upgrade.
- **No bonds cross the QM/MM boundary** — QM = a whole molecule, MM = separate
  molecules (solvent). Avoids link atoms.
- Small QM region + modest explicit environment; **cost-guarded + background**
  (same resource-awareness as Phase 1).

## Plan (TDD, each task green + committed)

- **T0 — Unit reconciliation.** Extend `units.rs`: `BOHR_PER_NM`,
  `HARTREE_PER_KJ_MOL`, force factor `hartree·bohr⁻¹ per kJ·mol⁻¹·nm⁻¹`, charge
  1:1 (e == a.u.), mass via the existing `AMU_TO_AU_MASS`. Tests vs known values.
- **T1 — MM forces wrapper.** `mm::mm_forces(system, ff) -> (energy_au, forces_au)`
  wrapping `Simulation::evaluate_forces`, converting GROMACS→a.u. Test on a water
  dimer (finite, sane forces; Newton's third law).
- **T2 — QM–MM mechanical coupling.** `coupling::qmmm_lj_coulomb(qm_pos, qm_lj,
  qm_q, mm_pos, mm_lj, mm_q) -> (energy_au, qm_forces, mm_forces)`. Test: two
  point charges → Coulomb force matches the analytic `k q₁q₂/r²`.
- **T3 — QmMmEngine.** `QmMmSystem` (QM region + MM `valenx_md::System` + QM
  LJ/charges) + `QmMmEngine: ReactionEngine`: assemble QM + MM + coupling forces,
  velocity-Verlet over all atoms, NVE/Berendsen, cost-guarded. Test: a QM atom +
  a few MM atoms → NVE total energy conserved.
- **T4 — Workbench QM/MM mode.** Reaction Dynamics gets a backend selector
  (AIMD | QM/MM); a "solute in water" preset; 3-D playback already renders all
  atoms (QM + MM) since it works off positions + detected bonds.
- **T5 (the "great" upgrade) — Electrostatic embedding.** Add external point
  charges to `valenx-qchem`'s one-electron Hamiltonian (the McMurchie-Davidson
  nuclear-attraction path already computes point-charge integrals — feed it the
  MM charges), expose through the driver, and switch the coupling's electrostatics
  into the QM Hamiltonian. The MM environment then polarizes the QM density.

## Out of scope

Polarizable MM, QM/MM boundary across bonds (link atoms), PME for the QM–MM
electrostatics, free-energy / enhanced sampling.
