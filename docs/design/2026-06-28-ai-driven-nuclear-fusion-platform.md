# AI-Driven Nuclear & Fusion Engineering Platform — Design Spec (v3)

**Date:** 2026-06-28
**Status:** **Design only — nothing is built. That's the asset:** the design gets
torn apart *before* concrete is poured. No code until a phase is greenlit.
**Scope:** **CIVILIAN power + research engineering only** — energy generation and the
science behind it. **Explicitly NOT weapons design** (no implosion, weapons-grade
enrichment optimization, criticality-for-a-device, or boosting). Everything below is the
open, university-taught, OpenMC / MOOSE / Grad–Shafranov-class science.

*v3 folds in six review passes — completeness (11 gaps), desktop-infrastructure, agent's-eye
drivability, valenx-native conventions (RFCs), a deep math/numerics pass, and an expert
nuclear/fusion physics-credibility pass — each verified against the actual repo, keeping only
what valenx does not already have.*

---

## 1. Goal

Extend valenx from *"an engineer drives a calculator"* into an **autonomous
engineering-R&D platform**: an AI agent **discovers the design space on its own**,
runs real coupled multiphysics, accelerates it with confidence-aware surrogates,
optimizes across competing objectives (physics + safety + economics), quantifies
uncertainty, and **surfaces validated, re-derivable breakthrough candidates** — first
for civilian **fission microreactors / SMRs** and **fusion energy devices**, then
generalized to *all* of valenx's engineering domains.

The differentiator is the **closed, self-describing, AI-driven discovery loop** wrapped
around validated native physics — not the calculators (those exist elsewhere).

---

## 2. What valenx ALREADY provides — build on, do NOT rebuild

(Verified in-repo. A naive review will recommend building these; they exist.)

- **GUI + render:** `eframe` / `egui` / `egui_plot` / **`wgpu`** desktop stack; `valenx-viz`
  crate; **VTK / VTU export** already wired (OpenFOAM/SU2/Elmer/CalculiX → `Results.fields`).
- **Non-blocking compute:** background worker threads already used in ~13 app sites.
- **Linear algebra:** **`faer`** + `nalgebra` + `nalgebra-sparse` + `ndarray` already deps.
- **Autonomy primitives:** `valenx-orchestrator` (funnel/sweeps), `valenx-uq`, `valenx-rom`,
  `valenx-topopt`; the **agent-command bridge**; the headless **`--self-test`** harness.
- **Provenance store:** **`valenx-audit`** — append-only log with a **SHA-256 chain** (ready
  to become candidate lineage).
- **HPC:** **`valenx-executor-slurm`** — real remote-cluster submission (rsync + sbatch/squeue/
  sacct + GPU `--gres` + multi-rank srun + fetch_results).
- **Physics to reuse:** `valenx-cfd-native`, `-fem`, `-fields`, `-thermo`, `-radioactivity`,
  `-qchem`, `-md`, `-cheminf`.

**Genuinely greenfield (what we build):** `valenx-neutronics`, `valenx-plasma`,
`valenx-autodesign`, `valenx-data`; the **agent interface contract**; the **nuclear-data
pipeline** (incl. self-serve **lattice/MOC group-constant generation — a transport-scale effort**); **GPU-ML surrogates**; **reference-code cross-validation adapters**; the
**provenance + NQA-1/V&V** wiring.

---

## 2a. Alignment with valenx's existing conventions (don't reinvent)

valenx already has RFCs + machinery this module must *follow*, not duplicate (verified in `rfcs/`):

- **RFC 0002 (adapter contract) + RFC 0007 (coupling adapters):** the new physics and the
  multiphysics coupling implement the existing **`prepare() → run() → collect()`** lifecycle and
  the `ProgressSink` contract. A coupled microreactor case is a **concurrent node in the workflow
  DAG**, automatically inheriting valenx's status colors, probes, and UI telemetry. **Do not build
  a custom coupling loop inside `valenx-neutronics` / `-plasma`.**
- **RFC 0004 (results & fields):** structured outputs + field data follow the existing
  `Results.fields` schema (so the renderer + VTK/VTU export work for free).
- **RFC 0009 (HPC job submission):** the compute layer (§8) uses the existing SLURM path.
- **RFC 0011 (parameter-sweep / optimization) + RFC 0012 (ML training-data export):** the autonomy
  spine **extends** these — the sweep + surrogate-training-data plumbing already has a contract.
- **RFC 0013 (audit / RBAC):** provenance (§6) builds on the existing audit chain.
- **Real-time telemetry:** the autodesign loop pipes convergence metrics (k-eff, Q, residuals)
  into **`RunContext`'s `ProgressSink`** (`valenx-core/src/adapter.rs`) so the UI draws live
  optimization-convergence curves **without locking the render thread.**
- **Local-first (manifesto "no cloud, no API keys"):** all surrogates run on the local
  workstation — POD / Gaussian-Process (Kriging) ROMs first (cheap to train locally),
  `burn`/`candle` neural surrogates only when warranted, on the user's CPU/GPU.

---

## 3. The agent interface contract — the headline drivability feature

The single thing that separates *"an agent cold-starts a campaign"* from *"a human tells
the agent what to type."* Every workbench / design object is **self-describing**:

- **`--describe` → machine-readable JSON schema** (OpenAPI-style) of the full design space:
  every parameter with **type, physical units, valid range, and constraints**, plus the
  **output schema**. Fetchable by the agent itself.
- An agent discovers the knobs and objectives with **zero human briefing**, then drives the
  loop. This extends the existing agent-bridge from "set named widget" to "here is the entire
  contract."
- Lives in `valenx-autodesign`; every new + existing workbench gains a `describe()`.

---

## 4. Architecture — four layers

### Layer 1 — Autonomy spine (the "breakthrough engine") · new `valenx-autodesign`

Extends `orchestrator` + `uq` + `rom` + `topopt`. Domain-agnostic; reusable everywhere.

1. **Self-describing parametric design objects** (§3) — the agent reads the contract.
2. **Headless physics evaluation** → structured `key value unit` results (the `--self-test`
   pattern) — **units are first-class and dimensionally checked** (a unit mismatch is the
   worst possible nuclear failure — cf. Mars Climate Orbiter), enforced with the **`uom`** units-of-measure crate (compile-time / runtime dimensional checks across the atoms→thermal-hydraulics multi-scale stack).
3. **AI-simulated surrogates** — fast NN / GP / ROM models, **GPU-accelerated via `burn` or
   `candle`** (WGPU/CUDA backends) — **local-first: POD / Gaussian-Process (Kriging) ROMs first (cheap to train on a workstation), neural surrogates only when warranted, no cloud models.** Each surrogate **reports its own confidence**, and a
   **machine-readable trust gate** decides when to trust it vs. fall back to the real solver
   — so "explore millions cheaply" cannot quietly launder surrogate error into the rankings.
4. **Dry-run cost/time estimate** — before committing a campaign, every eval exposes a cheap
   `estimate()` (compute + wall-time per run) so the agent can plan "200 sweeps overnight"
   without blowing the budget.
5. **Design-space exploration + multi-objective optimization** — Bayesian opt, genetic,
   gradient, topology-opt → Pareto fronts + novel designs.
6. **Active-learning autonomous discovery** — picks the next simulation that maximally
   improves the Pareto front / reduces uncertainty. The real automated-breakthrough engine.
7. **UQ + safety** (`valenx-uq`) — error bars + hard safety-margin gates on every candidate.
8. **Validate the OPTIMIZER, not just the physics** — the discovery loop must re-discover a
   *known-good* design from scratch (an existing SMR; the ITER operating point) as its own
   acceptance test. **Caveat:** re-discovering ITER is **weak validation if the confinement
   scaling was fit to ITER-class data** (circular) — prefer a **held-out reference design the
   model was not calibrated against.**
9. **The AI driver** — reads the contract, launches, reads results, decides, ranks with
   rationale.

### Layer 2 — Fission stack (microreactors / SMRs)

**Prerequisite — the nuclear-data pipeline** (the biggest physics dependency; greenfield):
ingest **ENDF/B cross-section libraries** + a **lattice-physics step** that homogenizes
continuous-energy cross-sections into the **few-group constants** diffusion needs. Without
this, no neutronics runs.

**Two distinct capabilities (don't conflate):** **(a)** the diffusion solver validated using
*provided* few-group constants — Phase 1, tractable; **(b)** *self-serve* cross-section
generation for **arbitrary geometry**, which requires a **lattice / MOC transport** step (later
tier — roughly as hard as the diffusion solver itself; fixed benchmarks ship with their own
constants, arbitrary designs do not). Designing a *novel* microreactor needs (b); reproducing
standard benchmarks needs only (a).

**New crate `valenx-neutronics`:**
- **Multigroup neutron diffusion** (3-D, on `valenx-fem` meshes) → `k-eff`, flux + power.
- **Point kinetics + reactivity feedback** (Doppler, void, temperature) → transients/safety.
- **Burnup / depletion** (Bateman — reuses `valenx-radioactivity` decay chains) → fuel cycle.
- *Later fidelity tier:* Monte Carlo neutron transport.

**Multiphysics coupling — via valenx's RFC 0007 framework + the `prepare→run→collect` lifecycle as a concurrent workflow-DAG node (inherits status/probes/telemetry), NOT a custom loop:** `neutronics ↔ thermal-hydraulics
(cfd-native) ↔ fuel/structure heat + stress (fem) ↔ decay heat + shielding + activation
(radioactivity)`.

**Microreactor Design workbench** — geometry, materials (enrichment, coolant), operating
point → k-eff, power profile, thermal margins, burnup, decay heat, shielding dose. AI-drivable.

**Phase-1 scope (diffusion's limits):** diffusion-friendly **LWR-type SMRs** are the Phase-1
target; **heat-pipe and fast-spectrum microreactors are deferred to the transport (MOC/Sₙ or
Monte Carlo) tier**, because few-group diffusion is unreliable for small, high-leakage, strongly
heterogeneous, fast-spectrum cores. Until transport lands, microreactor results outside the
LWR-like regime are **concept-screening only.**

**Validation** vs *civilian power-reactor* references: analytic k∞, the **IAEA 2D/3D PWR
diffusion benchmark**, the **LRA/LMW kinetics benchmarks**, and TRIGA cores (**C5G7 is deferred
to the transport tier — it is a transport (MOC/Sₙ) benchmark few-group diffusion cannot
reproduce**) — power-reactor benchmarks, not bare-sphere criticality.

### Layer 3 — Fusion stack (compact energy devices)

**New crate `valenx-plasma`:**
- **0-D power balance** — Lawson, triple product *nTτ*, gain **Q** → the viability gate.
- **Grad–Shafranov equilibrium** (2-D MHD) → plasma shape + poloidal field.
- **Confinement scaling** (ITER-class), heating, current drive.
- **Divertor heat loads** + **neutron-damage materials** (DPA / displacements-per-atom,
  transmutation, swelling under 14 MeV neutrons) — two of the hardest problems that *actually
  gate* a real fusion design.

**Coupling:** `plasma ↔ magnets (fields: superconducting coils) ↔ magnet/vessel stress +
thermal (fem) ↔ tritium-breeding blanket (reuses neutronics) ↔ cooling (cfd-native)`.

**Fusion Device workbench** — plasma + magnets + blanket → Q, confinement, magnet stress,
tritium breeding ratio, net power balance. AI-drivable.

**TBR needs transport, not diffusion:** credible **tritium breeding ratio** requires **fast-neutron
transport (multigroup Sₙ or Monte Carlo)** in the blanket — 14 MeV neutrons + (n,2n) multiplication
— and TBR is the *gating* fusion metric. The fusion stack therefore likely needs the **transport
tier the fission stack defers**; blanket TBR is **concept-level until transport lands.**

**Concept breadth:** **tokamak-first** (Grad–Shafranov). Noted but *out of initial scope*:
stellarators (3-D equilibrium — much harder), inertial confinement (laser fusion), and
**field-reversed-configuration / magnetized-target** ("energy cell"-style). Scoped explicitly
rather than implied.

**Honest scope:** 0-D/1-D + equilibrium is tractable and useful for *concept design*. Full
3-D turbulent transport (gyrokinetics) is HPC-scale → **out of scope / BLOCKED**, not faked.

### Layer 4 — "Engineering all kinds" generalization

The Layer-1 spine is domain-agnostic by construction. Once proven on nuclear, the same
machinery wraps every valenx domain. Nuclear is the flagship that *builds the general
capability*; "speed up the field" → "speed up every field valenx touches."

---

## 5. Cross-cutting capabilities

- **Fuel & materials discovery** — AI searches novel fuels / moderators / coolants / alloys;
  properties from `qchem` / `md` / `cheminf` (multi-scale: atoms → component → system).
- **Techno-economics in the loop** — a cost model so optimization targets **$/kWh (LCOE)**. *LCOE for first-of-a-kind designs carries large uncertainty — it is a **soft objective with wide UQ bounds, never a hard optimization target**; optimizing hard against an unreliable cost model drives designs toward whatever the model under-prices.*
- **Generative geometry** — AI *invents shapes* (topology-opt + generative parametric geometry
  on `cad` / `topopt`), not just tunes parameters.
- **Multi-scale coupling** — atoms → component → system → economics in one objective.
- **Literature grounding** — anchor designs to published scaling laws / material DBs (cite, not invent).
- **Safety & licensing pre-screen** — automated defense-in-depth / regulatory-margin checks.

---

## 6. Provenance, reproducibility & units

Every candidate is **re-derivable, not just reported** ("re-run beats recall" — applies to
our own outputs; regulators require it):

- Each result carries: **solver + version**, **which benchmark validated it**, a
  **surrogate-vs-real-solver flag**, the **random seed**, and full **data lineage** (exact
  code version, physics seeds, training data, config).
- **Deterministic regression** — all stochastic solvers (future Monte Carlo) and surrogate
  initializations run from **locked seeds** so results reproduce **deterministically within tolerance** — bit-exact reproduction is guaranteed only on single-thread CPU paths (GPU / parallel floating-point reductions are non-associative); locked seeds give **statistical, not bit-exact**, reproducibility on parallel/GPU runs.
- **Units first-class + dimensionally checked** in every structured output.
- Wire **`valenx-audit`** (SHA-256 chain) as the lineage store.
- **Mandatory human-engineering review gate** before any candidate is called a "breakthrough."

---

## 7. Validation, QA & honesty

- **Cross-validation against reference codes** — the credibility pillar: `valenx-neutronics`
  matches **OpenMC / Serpent / MOOSE** on benchmarks; `valenx-plasma` matches published
  equilibria. *Needs new reference-code adapters (none exist today.)*
- **`--self-test` extension** — every new core gets deep ground-truth checks.
- **NQA-1 awareness** — commercial SMR/microreactor work may eventually need **ASME NQA-1**
  (nuclear-facility software QA); the dev process + lineage are designed to support it.
- **Honesty framing (explicit):** this is a **design / research / screening aid — NOT a
  licensed safety-analysis code.** Real licensing requires NRC-certified codes; we do not
  claim certification.
- **Civilian boundary enforced mechanically, not just declared:** objective functions are
  restricted to **power output, safety margins, and cost — never weapon yield.** No bare-geometry
  criticality optimization, no enrichment-maximization objective, no implosion/hydrodynamics.
  Out-of-scope objectives are **BLOCKED at the contract level.**

---

## 8. Compute & HPC

Surrogates handle *exploration*; the real runs that *train and verify* them go to
**`valenx-executor-slurm`** (HPC). The **dry-run cost/time estimate (§4.4)** drives the
local-vs-cluster split and the campaign budget. Surrogate evals can use the local GPU
(`burn`/`candle`).

---

## 9. Data engineering · new `valenx-data`

Engineering projects evolve for years; saved files must not break on a struct change.

- **Tabular** (optimization logs, surrogate training vectors): **Apache Arrow / Parquet**
  (`arrow` / `polars`), not raw JSON.
- **Large meshes / fields:** **HDF5**.
- **User-saved parametric designs:** **schema-versioned** serialization (`prost` / Protobuf
  or explicitly versioned `serde`) for forward/backward compatibility.
- The **design / result / campaign database** for discovery runs lives here.

---

## 10. Visualization (mostly reuse)

Reuse `wgpu` + `valenx-viz`: 3-D **neutron-flux contours**, **plasma equilibria**, **magnetic
field lines**, mesh deformation, a **Pareto-front explorer**, and design-vs-design compare.
Computed fields (neutronics 3-D volumetric power; plasma 2-D poloidal flux) map directly into `valenx-viz` vertex buffers via the RFC 0004 `Results.fields` schema; VTK/VTU export (already present) feeds ParaView for heavy cases.

---

## 11. Data flow (one loop)

```
agent fetches the design-space CONTRACT (--describe → schema: params, units, ranges, outputs)
   → autodesign samples designs
      → dry-run cost/time estimate → plan the campaign
         → coupled physics eval (neutronics/plasma + cfd + fem + fields + radioactivity)
            → structured results WITH UNITS (k-eff, Q, margins, $/kWh, …)
               → surrogate train/update (GPU) + self-reported confidence
                  → multi-objective optimize + active-learning pick next
                     (trust gate: surrogate vs real-solver fallback)
                        → UQ + safety + licensing pre-screen
                           → ranked candidates, each with FULL PROVENANCE (solver/ver/seed/lineage)
   → human review gate → agent surfaces breakthroughs
```

Every arrow is headless, machine-parseable, dimensionally checked, and exits non-zero on any
failed safety gate (so the whole loop doubles as CI).

---

## 12. New / extended crates

- **New:** `valenx-neutronics`, `valenx-plasma`, `valenx-autodesign` (spine + interface
  contract), `valenx-data` (Arrow/Parquet/HDF5 + schema versioning), reference-code adapters
  (`valenx-adapter-openmc`, …).
- **Extended:** `fields` (superconducting magnets), `cfd-native` (reactor/blanket TH), `fem`
  (multiphysics coupling), `radioactivity` (decay-heat/activation), `uq`/`rom`/`topopt` (into
  autodesign), `audit` (candidate lineage), `executor-slurm` (campaign dispatch), `qchem`/
  `md`/`cheminf` (materials), `viz` (flux/equilibria/Pareto views).
- **Workbenches (`valenx-app`):** Microreactor Design, Fusion Device, Autonomous-Design
  control panel — all AI-drivable, all `--describe`-able, all `--self-test`-covered.

---

## 13. Phasing

- **Phase 0 — Structural skeleton + autonomy spine:** `valenx-autodesign` (incl. the
  **interface contract**, **surrogate trust gate**, **cost estimate**), `valenx-data`, and the
  **provenance skeleton** (`audit` wiring + units + seeds). **Prove the headless→GUI
  non-blocking pipeline + the full provenance/contract loop on a SIMPLE existing problem
  (1-D thermal conduction) before any neutronics.**
- **Phase 1 — Fission:** nuclear-data pipeline → `valenx-neutronics` → coupling → Microreactor
  workbench → AI design loop. Cross-validate vs OpenMC. **LWR-type SMRs first; heat-pipe /
  fast-spectrum microreactors are deferred to the transport tier (diffusion is unreliable there).**
- **Phase 2 — Fusion:** `valenx-plasma` (+ divertor/materials) → coupling (inherits blanket
  neutronics) → Fusion workbench → AI design loop.
- **Phase 3 — Enterprise / generalize:** economics, generative geometry, materials discovery,
  NQA-1 lineage tracking, hybrid local/cluster scaling, and the spine wrapped over the rest of
  valenx's domains.

**Per-phase success criteria + risks are mandatory.** Example — *Phase 1 done when:*
neutronics reproduces the **IAEA PWR diffusion benchmark** within tolerance **AND** matches OpenMC on a shared benchmark
**AND** the coupled loop converges **AND** the AI loop re-discovers a reference microreactor
within tolerance — each with full provenance.

---

## 14. Risks

- **Surrogate error laundering** into ranked candidates (mitigated by the confidence + trust gate, §4.3).
- **Multiphysics coupling convergence** (Picard divergence on stiff feedback).
- **Nuclear-data availability / licensing** (ENDF access + processing).
- **Plasma-model fidelity** (0-D/1-D may mislead on concepts that need transport detail).
- **Agent over-trust** of surrogates / under-specified objectives.
- **Compute cost** of training surrogates at useful fidelity.
- **Lattice-physics / group-constant generation** for arbitrary geometry is a **transport-scale effort, not a simple pipeline** (Layer 2 capability (b)).
- **Fusion TBR requires transport** (Sₙ/MC), not diffusion — the *gating* fusion metric is concept-level until the transport tier lands.
- **Optimizer model-exploitation** (distinct from surrogate laundering, §4.3): a multi-objective optimizer is **adversarial against *any* model, including the real solver** — a "breakthrough" Pareto point is often a **discretization/mesh artifact, not a discovery.** *Mitigation:* mandatory **mesh-convergence check** + an **independent cross-code re-run** (OpenMC/Serpent/published equilibria) on every top candidate **before** the human-review gate.

---

## 15. Out of scope (honest)

- **Weapons design of any kind** — excluded and never assisted.
- Full 3-D plasma turbulence / gyrokinetics (HPC-scale).
- Production Monte Carlo neutron transport at first (diffusion first; MC is a later tier).
- Non-tokamak fusion concepts in v1 (stellarator/ICF/FRC noted, not built initially).
- **NQA-1 certification** itself (designed to *support* it; not claimed).
- Real experimental calibration until real reactor/fusion data exists (digital-twin hook is
  designed-for, not claimed).

---

## 16. Implementation architecture (Phase 0/1 build blueprint)

*Synthesized from a deep math/numerics review + an implementation-architecture review. A
reference for the build — still not code.*

### Crate tree (5 new)

```
valenx-autodesign/  contract.rs (--describe schema) · surrogate/{kriging, neural(DeepONet)}
                    · optimizer/{pareto, active, anderson}
valenx-data/        schema/design.proto (prost) · storage/{columnar (Arrow/Parquet),
                    volume (HDF5 + conservative mesh projection)}
valenx-neutronics/  data/{cross_sections (memmap2 ENDF/B), lattice (few-group)}
                    · solver/{diffusion, kinetics} · adapter.rs
valenx-plasma/      zero_d (Lawson/Q) · mhd/{grad_shafranov, grid} · materials (DPA/divertor)
valenx-adapters-nuclear/  openmc.rs · moose.rs (cross-validation)
```

### Solver formulations (on `faer`)
- **Fission:** 3-D multigroup neutron diffusion `−∇·(D_g ∇φ_g) + Σ_tg φ_g = sources`, groups g∈[1,G].
- **Fusion:** 2-D axisymmetric Grad–Shafranov `Δ*ψ = −μ₀R² dp/dψ − F dF/dψ` on the (R,Z) plane.

### Advanced numerics (the parts that make it converge + scale)
- **Anderson acceleration** (fixed-point history queue, `faer::Mat`) — stops spatial divergence in
  the tight neutronics↔thermal void/temperature feedback coupling; superior to naive Picard.
- **Memory-mapped group constants** (`memmap2`, zero-copy O(1)) — temperature-dependent ENDF/B
  cross-section lookups without ASCII-parse overhead.
- **Conservative mesh-field projection** — map fields across dissimilar meshes (CFD finite-volume ↔
  FEM lattice) while preserving total integrated energy.
- **DeepONet local surrogates** (`burn`) — query distributed field values (flux profiles) without
  discretizing onto a standard mesh.

### Core contracts (bound to valenx's REAL `valenx-core::adapter::Adapter`)
- `DesignSpaceSchema` / `ParameterBounds` / `ObjectiveGoal` / `ConstraintDefinition` — the `--describe` JSON.
- `SelfDescribingDesign { describe() -> DesignSpaceSchema; estimate_cost() -> ComputeEstimate }`.
- `LocalSurrogateModel { train; predict -> PredictionResult { mean, confidence_variance } }` +
  `SurrogateTrustGate { variance_threshold }` (variance < threshold ⇒ trust surrogate, else fall back to solver).
- Provenance `.proto`: `ArchitecturalCandidate` + `ProvenanceRecord { valenx_version_hash, solver_id,
  deterministic_seed, is_surrogate_prediction, surrogate_confidence, validation_benchmark, sha256_audit_link }`.

**⚠ API binding (verified vs `valenx-core/src/adapter.rs`).** The nuclear adapters implement the real
**`trait Adapter`**: `prepare(&self, case, workdir) -> PreparedJob`, `run(&self, job, &mut RunContext)
-> RunReport`, `collect(&self, job) -> Results`; errors are **`AdapterError`** (not `String`); progress
is **`ctx.report_progress(pct: f32 /*0–100*/, msg)`** (not `update_status`). Reviewer drafts assuming a
`PhysicsAdapter`/`update_status` shape must be re-bound to this before they compile.
