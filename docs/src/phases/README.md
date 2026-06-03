# Per-phase acceptance docs

The 20-year [ROADMAP](../../../ROADMAP.md) plans 16 phases. This
directory holds the per-phase **acceptance criteria** — what
"Phase N ships" actually means, phrased as checklists that anyone
can tick off independently of the roadmap narrative.

Each doc has the same structure:

1. **Goal** — one sentence on why this phase exists.
2. **Capability inventory** — the specific user-visible features the
   phase adds.
3. **Integrated tools** — which external solvers graduate from
   "scaffold" to "implemented" during this phase.
4. **Acceptance checklist** — boxes we tick when the phase is
   declared done.
5. **Success metrics** — numeric targets (install time, first-solve
   latency, test-suite coverage, validation benchmarks passed).
6. **Next phase leads into** — a forward pointer.

These files are living; they drift as we learn what the work
actually takes. When a phase closes, its checklist is frozen and a
short retrospective appended.

## Index

| Phase | Title                                    | Status      |
|-------|------------------------------------------|-------------|
| [1](./phase-01-foundation.md)         | Foundation + first CFD thread        | 🔵 Complete    |
| [2](./phase-02-cad-meshing.md)        | CAD + meshing integration            | 🟢 In progress |
| [3](./phase-03-fea.md)                | Finite-element analysis              | 🟢 In progress |
| [4](./phase-04-chemistry.md)          | Reaction kinetics + combustion       | 🟢 In progress |
| [5](./phase-05-md.md)                 | Molecular dynamics                   | 🟢 In progress |
| [6](./phase-06-em.md)                 | Electromagnetics                     | 🟢 In progress |
| [7](./phase-07-battery.md)            | Battery modelling                    | 🟢 In progress |
| [8](./phase-08-dynamics.md)           | Multibody + robotics                 | 🟢 In progress |
| [9](./phase-09-coupling.md)           | Multi-physics coupling (preCICE)     | 🟢 In progress |
| [10](./phase-10-ux-polish.md)         | UX polish + first public release     | ⚪ Planned    |
| [11](./phase-11-hpc.md)               | HPC / cluster execution              | ⚪ Planned    |
| [12](./phase-12-optimisation.md)      | Optimisation + adjoint workflows     | ⚪ Planned    |
| [13](./phase-13-ml.md)                | ML-assisted surrogates               | ⚪ Planned    |
| [14](./phase-14-plugins.md)           | Plugin marketplace                   | ⚪ Planned    |
| [15](./phase-15-enterprise.md)        | Enterprise deployment                | ⚪ Planned    |
| [16](./phase-16-stewardship.md)       | Stewardship + long-term governance   | ⚪ Planned    |

**Legend:** 🟢 in progress · ⚪ planned/scaffolded · 🔵 done.
