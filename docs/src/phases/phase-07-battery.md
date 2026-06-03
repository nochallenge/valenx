# Phase 7 — Battery modelling

**Status:** 🟢 In progress — PyBaMM adapter live for single-protocol discharge / charge on DFN / SPM / SPMe with built-in parameter sets.

## Goal

Electrochemical battery modelling (DFN, SPM, SPMe) inside Valenx so
users can co-simulate pack thermals (CFD/FEA) with cell
electrochemistry on one timeline.

## Capability inventory

- Doyle-Fuller-Newman (DFN) pseudo-2D cell model.
- Single-particle model (SPM) and the electrolyte-enhanced variant
  (SPMe) for speed/accuracy trade-offs.
- Parameterisation library: LFP, NMC, NCA, Si-anode blends.
- Cycling protocols: CC, CCCV, GITT, EIS, drive cycles.
- Degradation models: SEI growth, lithium plating, particle cracking.
- Capacity-fade and power-fade fitting against cycling data.

## Integrated tools graduating to Implemented

| Tool    | Adapter crate                  | Role                                     |
|---------|--------------------------------|------------------------------------------|
| PyBaMM  | `valenx-adapter-pybamm`        | Primary electrochemical solver           |

## Acceptance checklist

- [ ] DFN simulation of a commercial 18650 cell, voltage vs capacity
      plot.
- [ ] Drive-cycle simulation against a WLTP cycle.
- [ ] Degradation study: 500 cycles of CCCV with SEI growth.
- [ ] Thermal coupling with CFD pack simulation via preCICE.
- [ ] Side-by-side DFN vs SPM vs SPMe comparison at the same operating
      point.

## Success metrics

| Metric                                              | Target      |
|-----------------------------------------------------|-------------|
| DFN single-cycle runtime                            | < 10 s      |
| SPM drive-cycle runtime                             | < 2 s       |
| Parameter-sweep UI works for 100-point sweeps       | yes         |

## Leads into

[Phase 8 — Multibody + robotics](./phase-08-dynamics.md).
