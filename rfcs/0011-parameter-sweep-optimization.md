# RFC 0011 — Parameter sweeps + optimization

| | |
|---|---|
| **Status** | Draft |
| **Phase target** | 12 (optimization) |
| **Type** | Architecture |
| **Author** | Maintainers |
| **Discussion** | open in PR — not yet merged |
| **Related** | RFC 0009 (HPC executors), RFC 0001 (project file format) |

## Summary

Define a small `Optimizer` trait + a per-case `[sweep]` /
`[optimize]` block in `case.toml` that turns a single `case.toml` +
some parameter declarations into N derived runs. Three reference
optimizers ship in the first cut:

- **`grid`** — full Cartesian product over declared parameters.
- **`latin-hypercube`** — N samples maximally spread across the
  parameter space (good for surrogates).
- **`gradient-descent`** — finite-difference gradient + line search
  for differentiable objectives (good for shape opt).

This RFC is a scaffold — the trait + case schema + grid optimizer
are the in-scope deliverables; LHS and gradient descent are sketched
but their implementation is follow-up work.

## Motivation

Engineers running CFD/FEA/etc. don't do single runs. They sweep:

- "What happens to drag if I vary the angle of attack from 0° to
  10° in 1° steps?"
- "What inlet pressure minimises pressure drop subject to a
  cooling-target constraint?"
- "Generate 200 runs with random parameter combinations to train a
  surrogate model" (Phase 13 ML feeds this).

Today every "sweep" requires either hand-editing N copies of
`case.toml` or scripting outside Valenx. Both work; both miss the
point of having a unified workflow loop.

## Design

### Per-case sweep / optimize block in case.toml

```toml
[case]
format  = "1.0"
name    = "airfoil-aoa-sweep"
physics = "cfd"
solver  = "openfoam.simpleFoam"
mesh    = "primary"

# … existing flow / boundaries / solve sections …

[sweep]
optimizer = "grid"

[[sweep.parameter]]
# Path into the case.toml using TOML-pointer syntax.
path = "/boundaries/inlet/velocity/0"
values = [10.0, 20.0, 30.0, 40.0, 50.0]

[[sweep.parameter]]
path = "/flow/turbulence"
values = ["kEpsilon", "kOmegaSST", "SpalartAllmaras"]
```

The grid optimizer takes the Cartesian product of declared
`values` arrays — 5 × 3 = 15 derived runs in this example. Each
derived run is a fresh `case.toml` written to a child directory
under the parent workdir, with the substitution applied.

### `Optimizer` trait

```rust
pub trait Optimizer {
    fn id(&self) -> &str;

    /// Plan the runs this optimizer will execute. Return value is
    /// a list of (run_id, derived_case_toml) pairs the harness
    /// will then submit through the regular Adapter pipeline.
    fn plan(
        &self,
        base_case: &CaseDef,
        sweep: &SweepConfig,
    ) -> Result<Vec<DerivedCase>, OptimizerError>;

    /// Once a batch of derived runs completes, the optimizer can
    /// either declare the sweep done (e.g. grid sweep — all runs
    /// scheduled) or request more runs (gradient descent — pick
    /// the next step based on the gradient observed so far).
    fn step(
        &mut self,
        completed: &[CompletedRun],
    ) -> OptimizerStep;
}

pub enum OptimizerStep {
    Done,
    More(Vec<DerivedCase>),
}
```

### Reference implementations

- **`GridOptimizer`** (in scope): emits the Cartesian product on
  `plan()`, returns `Done` immediately on `step()`.
- **`LatinHypercubeOptimizer`** (sketched): N pseudo-random samples
  via Latin Hypercube Sampling. Same one-shot contract as grid.
- **`GradientDescentOptimizer`** (sketched): emits a single
  baseline + N finite-difference probes on `plan()`. On `step()`
  computes the gradient, picks a step direction, returns
  `More(...)` with the next probe set. Repeats until objective
  stops improving (configurable patience).

### Objective extraction

The optimizer reads each completed run's `Results.scalars` to
extract its objective value:

```toml
[sweep.objective]
# A name from the Results.scalars catalog the adapter populates.
metric = "drag_coefficient"
direction = "minimize"  # or "maximize"
```

This couples nicely with the adapter Results-pipeline work that
shipped in Phases 1-9 — every live adapter now writes scalars or
fields to the catalog, so the same `Cd` extracted by OpenFOAM's
`forceCoeffs` function object becomes the objective for the
optimizer without further plumbing.

### Concurrency

The optimizer's planned runs go through the regular run-pipeline
one at a time today; with RFC 0009's `Executor` trait, derived
runs can be dispatched in parallel to a SLURM cluster — natural
multi-tier scaling.

## Drawbacks

- TOML-pointer syntax for parameter paths is a small DSL we'll need
  to document, support in error messages, and validate. JSON
  Pointer (RFC 6901) is the obvious model.
- "Substitution" is straightforward for scalar values but fiddly
  for nested tables / arrays. Out-of-scope for this RFC.
- Constraint handling (e.g. "minimise drag subject to Cl > 0.8")
  is a real optimization-theory topic. Today's proposal punts:
  the optimizer reads multiple scalars and lets the user combine
  them in a derived "objective" expression. Real constraint
  optimization is Phase 12.2.

## Alternatives considered

- **External Python script driving runs** — works today (call
  `valenx --headless run case.toml` in a loop). Rejected as the
  primary path: defeats the unified-workflow goal and gives users
  no progress / cancellation / results-aggregation.
- **A separate "sweep manager" adapter** — works architecturally
  but loses the static-validation that `[sweep]` lives next to
  `[flow]` / `[structural]` in the same `case.toml`.

## Migration path

Phase 12.0 (this RFC):
- Land `Optimizer` trait + `SweepConfig` parser + `GridOptimizer`
  in a new crate `valenx-optimize`.
- App gains a "Sweep" pane next to "Run" in the Run menu.
- Results aggregation: one `Results` per derived run, the optimizer
  collects them into a per-sweep summary table.

Phase 12.1: LatinHypercubeOptimizer + LhsOptimizerError.

Phase 12.2: GradientDescentOptimizer + finite-difference
infrastructure + line search.

Phase 12.3: Constraint handling + multi-objective Pareto fronts.

## Open questions

1. **Convergence detection across the sweep.** When does the user
   see "the sweep is done" vs "we have an interim best, keep going"?
   The trait's `OptimizerStep::Done` signals it but the UX is
   richer than a single bool.
2. **Per-run cancellation.** Cancel one derived run vs cancel the
   whole sweep. Maps to two cancel buttons in the UI.
3. **Storage.** N derived runs × M timesteps × P fields can grow
   large. Need a "keep only the K best" mode for big sweeps.
4. **Random seed.** LHS / future-stochastic optimizers need a seed
   for reproducibility. Belongs in the `[sweep]` block.

## References

- RFC 0009 — HPC executors (parallel dispatch of derived runs).
- RFC 0001 — project file format (where derived runs live on disk).
- JSON Pointer (RFC 6901) — model for the parameter `path` syntax.
- DAKOTA, OpenMDAO — reference designs for the trait's high-level
  shape (both have a similar plan/step iteration model).
