# Phase 9 — Multi-physics coupling

**Status:** 🟢 In progress — preCICE meta-adapter live for config staging + validation + coupling-manifest emission. The concurrent-orchestration contract (launch each participant's solver against the shared interface, aggregate their Results) lives at RFC 0007 and is the Phase 9 tail.

## Goal

Stop treating each physics as an island. Let a fluid solver and a
structural solver exchange surface traction + displacement, or a
thermal solver feed temperatures into a reactive-flow CFD run, on the
same timeline, inside the same project.

## Capability inventory

- Partitioned coupling via preCICE meta-adapter.
- Mesh-to-mesh data mapping: nearest-neighbour, RBF, consistent /
  conservative, on sliding and non-conforming interfaces.
- Implicit and explicit coupling schemes with Aitken / IQN-ILS
  acceleration.
- Sub-cycling: solvers advance at different time-step sizes with a
  convergence check at the coupling window.
- Coupling UI: draw interfaces between participants visually, see
  live iteration-residual plots per exchange.
- Named coupling presets: FSI, CHT, reactive-flow, multi-scale.

## Integrated tools graduating to Implemented

| Tool     | Adapter crate                        | Role                          |
|----------|--------------------------------------|-------------------------------|
| preCICE  | `valenx-adapter-precice`             | Coupling orchestrator         |

See [RFC 0007](../../../rfcs/0007-coupling-adapters.md) for how Valenx
models the coupling semantics internally — separate from any one
tool's details.

## Acceptance checklist

- [ ] FSI: OpenFOAM + CalculiX on a cantilevered flap in cross-flow.
- [ ] CHT: OpenFOAM fluid + CalculiX solid on a fin heat-exchanger
      section.
- [ ] Reactive flow: Cantera kinetics feeding an OpenFOAM compressible
      run.
- [ ] Coupling residual plotted alongside each participant's
      residuals.
- [ ] Failure recovery: if one participant crashes, the other stops
      cleanly with a usable diagnostic.

## Success metrics

| Metric                                              | Target      |
|-----------------------------------------------------|-------------|
| FSI for a 100k-cell fluid + 10k-element solid       | < 1 hr      |
| Reactive-flow coupling overhead vs pure CFD         | < 20 %      |
| Multi-participant dashboard visible while running   | yes         |

## Leads into

[Phase 10 — UX polish + first public release](./phase-10-ux-polish.md).
