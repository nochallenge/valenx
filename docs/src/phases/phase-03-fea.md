# Phase 3 — Finite-element analysis

**Status:** 🟢 In progress — CalculiX live for linear-static + linear-dynamic + steady/transient thermal FEA; Elmer live for steady-state + transient heat equation; Code_Aster / OpenRadioss still scaffolded.

## Goal

Structural and thermal FEA becomes a first-class physics in Valenx,
using the same ribbon / browser / viewport shell that CFD already
uses.

## Capability inventory

- Linear and nonlinear static analysis.
- Modal and harmonic response.
- Transient implicit (dynamic) and explicit (crash) solvers.
- Steady and transient heat transfer.
- Contact pairs with friction.
- Material library covering isotropic, orthotropic, elastoplastic,
  and hyperelastic models.
- Post-processing: displacement, stress (von Mises, principal, Cauchy),
  strain, reaction forces, temperature, heat flux.

## Integrated tools graduating to Implemented

| Tool         | Adapter crate                       | Role                              |
|--------------|-------------------------------------|-----------------------------------|
| CalculiX     | `valenx-adapter-calculix`           | 🟢 Live — linear static / linear dynamic (`*DYNAMIC`) / modal / steady + transient thermal deck emission + ccx subprocess + .frd/.dat artifact harvest |
| Code_Aster   | `valenx-adapter-code-aster`         | Industrial-grade nonlinear        |
| Elmer        | `valenx-adapter-elmer`              | 🟢 Live — steady-state + transient (BDF-2) heat equation SIF emission with optional initial-temperature block + ElmerSolver subprocess + .vtu/.result artifact harvest |
| OpenRadioss  | `valenx-adapter-openradioss`        | Explicit crash / impact           |

## Acceptance checklist

- [ ] Linear static study on a cantilever, von Mises contour in the
      viewport.
- [ ] Modal study with mode-shape animation.
- [ ] Thermal steady-state on a heat sink, temperature contour.
- [ ] Contact pair between two solids with friction.
- [ ] Nonlinear material (bilinear elasto-plastic) on a tensile
      specimen.
- [ ] Explicit crash of a thin-walled box against a rigid wall.
- [ ] Compare two solvers on the same problem, show result diff
      overlay.

## Success metrics

| Metric                                              | Target     |
|-----------------------------------------------------|------------|
| Linear static for 1M-DOF model                      | < 60 s     |
| Result-load time for 500 MB `.frd`                  | < 5 s      |
| Cross-solver result diff tool available             | yes        |

## Leads into

[Phase 4 — Reaction kinetics + combustion](./phase-04-chemistry.md).
