# valenx Connective + Evaluation Layer — build queue (5 PRs)

**Date:** 2026-06-24
**Status:** queued for the autonomous build loop (branch `feat/capabilities`), **front priority** —
this is the project's biggest readiness gap (turning ~255 isolated solvers into an integrated,
validated, reproducible suite). Source: user-supplied spec, 2026-06-24.

## Reconciliations (how this maps onto the running loop, honestly)

- **PR2 extends the existing `valenx-uq`** (shipped `726d22d` this session: Sobol indices, LHS, MC,
  polynomial surrogate). PR2 adds the missing pieces — the **V&V harness** (ConvergenceStudy /
  observed order, GCI, Method of Manufactured Solutions, error norms), **PCE** (Hermite/Legendre),
  and the **quasi-random Sobol *sequence***. Not a rebuild.
- **Gating is scoped** (`cargo {build,test,clippy,fmt,doc} -p <crate>`) per AGENTS.md "scope first";
  we do **not** run `cargo test --workspace` (valenx-app's rfd file-dialog tests hang headless —
  docs/QA.md). Each new crate here is non-GUI logic, so scoped-green is sound. Workspace `build`/
  `clippy` may be used to "widen" since only the GUI *tests* hang.
- **"PR" = local commit** (GitHub HELD until an explicit "push"). One logical change per commit,
  conventional style, email-safe identity, leak-scan 0.
- **Deps:** `num-complex` (DMD), `rand`/`rand_distr` (if genuinely needed) added **pinned** to root
  `[workspace.dependencies]`, inherited `.workspace = true`. No new linear-algebra stack (use
  `nalgebra`; `ndarray = 0.15` already pinned). Only ONE new-crate/new-dep subagent at a time (root
  `Cargo.toml` + `Cargo.lock` are shared).
- **AGENTS.md non-negotiables apply:** fail-loud on uncalibrated/unsupported input; every numerical
  claim pinned to a published/analytic benchmark; never delete crates/features.

## The 5 PRs (execution order PR1 → PR5; 1–3 are the load-bearing spine)

### PR1 — `valenx-rom` (reduced-order & surrogate modeling)  ← first, new crate
POD (energy-truncated SVD basis; project/reconstruct/error), DMD (standard + exact; complex
eigenvalues, modes, growth rates, frequencies), Operator Inference (LS reduced linear/quadratic
operators + reduced time-stepping), POD-Galerkin projection helper. Pin: rank-2 field → 2 singular
values, recon < 1e-10; DMD on `x_{k+1}=λx_k` → λ mag/phase to < 1e-6; OpInf reproduces a linear
ODE's reduced operator + forecasts held-out steps.

### PR2 — extend `valenx-uq` (UQ + V&V harness)  ← second, existing crate
Add: Sobol' **sequence** (low-discrepancy), **PCE** (Hermite/Legendre, low dim, mean/variance), and
the **V&V harness** — `ConvergenceStudy` (observed order of accuracy via log-error/log-h slope),
**GCI** (Roache, safety factor), **Method of Manufactured Solutions** helper, error norms (L2, L∞,
RMS). Pin: Sobol' on Ishigami (already in); observed_order ≈ 2.0 on a 2nd-order MMS; GCI vs Roache's
worked example; sampler mean convergence with a seeded draw.

### PR3 — `valenx-adapter-fmi` + co-sim master  ← third
Path: `crates/valenx-adapters/coupling/valenx-adapter-fmi` (sibling of `valenx-adapter-precice`).
FMI 2.0/3.0 co-sim import (parse `modelDescription.xml`; native Rust co-sim master, Jacobi +
Gauss-Seidel macro-stepping); binary FMU loading behind a `binary-fmu` feature, native in-process
`Subsystem` path as default. DIS (IEEE 1278.1) Entity State PDU encode/decode (bit-exact). Pin:
two-mass spring-damper split into two Subsystems co-sim'd (Gauss-Seidel) matches monolithic ODE to
< 1e-3; DIS PDU round-trip + reference byte layout.

### PR4 — `valenx-ppi` (PPI / interactome from sequence)
Depends on `valenx-align`, `valenx-binder-score`, `valenx-dock`, `valenx-biostruct`. Interface
contact prediction from a paired MSA (coevolution/MI over inter-chain column pairs); PPI confidence
aggregating coevolution + docked-pose complementarity + interface quality; host×pathogen interactome
screen. **Mirror `valenx-binder-score`: research heuristic only, never a validated "interacts"
verdict, flag for human review, fail-loud.** Pin: ranking AUROC ≥ floor on a toy benchmark;
deterministic score on a fixed pair; contact precision@(L/5) ≥ floor.

### PR5 — provenance + closed-loop campaign (capstone)
Extend `valenx-repro` (content-addressed PROV-style manifest DAG: hash {inputs, params, code
version, outputs} → CID; `record`/`replay`), `valenx-orchestrator` (`Campaign`: propose runs from a
policy/Bayesian-opt via `valenx-optimize` → execute → record provenance → evaluate with `valenx-uq`
→ iterate to a stopping criterion), `valenx-mcp` (expose `campaign.step()` / `replay()` as agent
tools). Pin: replay determinism (identical output hash, pinned CID); 5-iteration campaign on a convex
quadratic → analytic minimum within tol.

## Loop integration
Runs at the front of the queue, interleaved with the in-flight COLMAP-finish + defense M-tracks,
~4 disjoint crates concurrent (one new-crate at a time). Each PR: build subagent (TDD + scoped-gated
green) → orchestrator commits (email-safe, leak 0, real files only) → code + security review to zero
findings → next.
