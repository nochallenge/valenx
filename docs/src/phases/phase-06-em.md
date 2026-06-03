# Phase 6 — Electromagnetics

**Status:** 🟢 In progress — openEMS adapter live for rectangular-box FDTD with Gauss / Sine excitation; Meep still scaffolded.

## Goal

FDTD time-domain electromagnetics — antenna design, EMC, photonics —
land as a peer to CFD and FEA rather than an afterthought.

## Capability inventory

- 3D structured-grid FDTD solvers.
- Near-field → far-field transforms for antenna patterns.
- S-parameter extraction for ports.
- Dispersive-material models (Drude, Lorentz, Debye).
- Perfectly-matched layer (PML) absorbing boundaries.
- Frequency-domain probes via on-the-fly DFT.
- Photonics-oriented workflows (Meep) alongside
  antenna-oriented ones (openEMS).

## Integrated tools graduating to Implemented

| Tool    | Adapter crate                  | Role                                     |
|---------|--------------------------------|------------------------------------------|
| openEMS | `valenx-adapter-openems`       | 🟢 Live — rectangular FDTD + Gauss/Sine excitation + Mur/PML/PEC boundaries, driven via generated Octave script |
| Meep    | `valenx-adapter-meep`          | Photonics / optical FDTD                 |

## Acceptance checklist

- [ ] Patch antenna S11 sweep using openEMS.
- [ ] Far-field radiation pattern visualiser.
- [ ] Photonic crystal bandgap extraction using Meep.
- [ ] Frequency-domain field plot overlayed on geometry.
- [ ] Coupling with FEA for thermal feedback via preCICE (Phase 9+).

## Success metrics

| Metric                                              | Target      |
|-----------------------------------------------------|-------------|
| 10 M-cell FDTD run on consumer laptop               | < 10 min    |
| S-parameter sweep (50 frequency points)             | single pass |

## Leads into

[Phase 7 — Battery modelling](./phase-07-battery.md).
