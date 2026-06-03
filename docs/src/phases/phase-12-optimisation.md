# Phase 12 — Optimisation + adjoint workflows

**Status:** ⚪ Planned.

## Goal

Turn parametric models into optimisation studies — shape, topology,
and sizing — using gradient information from adjoint-capable solvers.

## Capability inventory

- Parameter sweeps with response-surface interpolation.
- Gradient-based optimisers: BFGS, L-BFGS-B, trust-region.
- Gradient-free: Nelder-Mead, CMA-ES, Bayesian optimisation.
- Shape optimisation using SU2's adjoint + OpenFOAM's adjointFoam.
- Topology optimisation for structures (TopOpt library bridge).
- Design-of-experiments: LHS, Sobol, full-factorial, adaptive.
- Multi-objective: NSGA-II, Pareto-front visualiser.

## Acceptance checklist

- [ ] Shape-optimise an airfoil for drag via SU2 adjoint.
- [ ] Topology-optimise a bracket for stiffness-to-weight.
- [ ] Pareto-front plot for a two-objective optimisation.
- [ ] Kill and resume a study mid-way without losing evaluations.

## Leads into

[Phase 13 — ML-assisted surrogates](./phase-13-ml.md).
