# Phase 4 — Reaction kinetics + combustion

**Status:** 🟢 In progress — Cantera adapter live for TP / HP / UV equilibrium.

## Goal

Cantera's thermodynamic / kinetic / transport machinery becomes a
first-class domain, usable on its own (0-D reactor, 1-D flames) and
as the chemistry backend for reactive-flow CFD.

## Capability inventory

- Thermodynamic equilibrium (T-P, H-P, U-V, S-P, …).
- 0-D constant-volume and constant-pressure reactors.
- 1-D freely propagating and burner-stabilised flames.
- Perfectly stirred reactor + plug-flow reactor models.
- Transport: mixture-averaged, multicomponent, Soret.
- Surface chemistry on catalyst beds.
- Large mechanism support (GRI-3.0, USC-II, curated via YAML).

## Integrated tools graduating to Implemented

| Tool    | Adapter crate                  | Role                                     |
|---------|--------------------------------|------------------------------------------|
| Cantera | `valenx-adapter-cantera`       | 🟢 Live — TP/HP/UV equilibrium via Python subprocess + summary JSON |

## Acceptance checklist

- [ ] Load a CHEMKIN mechanism, convert to YAML inside the app.
- [ ] Equilibrium calculator with thermodynamic-state plots.
- [ ] Reactor network editor (reactor + flow-controller + wall icons).
- [ ] 1-D laminar flame with adjustable fresh-gas composition.
- [ ] Reactive-flow coupling: Cantera kinetics inside OpenFOAM
      combustion case (through preCICE, once Phase 9 lands).

## Success metrics

| Metric                                              | Target     |
|-----------------------------------------------------|------------|
| GRI-3.0 equilibrium (C2-C3 fuels)                    | < 200 ms   |
| Freely propagating flame, 53 species                | < 30 s     |
| Mechanism viewer scales to 2000+ species            | yes        |

## Leads into

[Phase 5 — Molecular dynamics](./phase-05-md.md).
