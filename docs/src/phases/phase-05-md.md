# Phase 5 — Molecular dynamics

**Status:** 🟢 In progress — LAMMPS adapter live for NVE Lennard-Jones demos; GROMACS still scaffolded.

## Goal

Classical MD is usable from the same app — for both materials
(LAMMPS) and biomolecules (GROMACS) — sharing the same run-queue,
provenance, and visualisation infrastructure as the continuum work.

## Capability inventory

- Force fields: Lennard-Jones, EAM / MEAM, ReaxFF, Stillinger-Weber,
  Tersoff (LAMMPS); AMBER, CHARMM, GROMOS, OPLS (GROMACS).
- Ensembles: NVE, NVT (Nose-Hoover, Langevin), NPT
  (Parrinello-Rahman, MTTK).
- Integration schemes: Verlet, leapfrog, rRESPA.
- Trajectory visualisation with per-atom colouring, time scrubbing,
  RDFs, MSD, autocorrelation.
- Topology / `.pdb` / `.gro` importers.
- Solvation, ion placement, boundary conditions.

## Integrated tools graduating to Implemented

| Tool    | Adapter crate                  | Role                                     |
|---------|--------------------------------|------------------------------------------|
| LAMMPS  | `valenx-adapter-lammps`        | 🟢 Live — NVE / NVT / NPT + LJ / EAM pair styles + procedural or read_data init |
| GROMACS | `valenx-adapter-gromacs`       | Biomolecular MD                          |

## Acceptance checklist

- [ ] LJ argon equilibrium MD, run and visualise trajectory.
- [ ] EAM simulation of FCC copper with defect visualisation.
- [ ] GROMACS: simulate a 10 k-atom protein in explicit water, get
      RMSD plot.
- [ ] MSD + RDF post-processing baked into Valenx (not just shelled
      out).
- [ ] Convert LAMMPS dump → `valenx-fields::Results` losslessly.

## Success metrics

| Metric                                              | Target        |
|-----------------------------------------------------|---------------|
| LAMMPS LJ, 10k atoms, 1M steps                      | < 2 min       |
| GROMACS md.log → residuals panel                    | live streaming|
| Trajectory viewer for 1M+ atoms                     | 60 fps        |

## Leads into

[Phase 6 — Electromagnetics](./phase-06-em.md).
