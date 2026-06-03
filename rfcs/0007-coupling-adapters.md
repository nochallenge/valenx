# RFC 0007 — Coupling adapters

- **Status:** Draft
- **Authors:** @valenx-maintainers
- **Created:** 2026-04
- **Target phase:** 9

## Summary

Define how Valenx composes two or more `Adapter` implementations into
a single coupled run — concretely, what the orchestrator (`preCICE`
today, potentially others later) sees, and what participants see from
inside their own `prepare() → run() → collect()` lifecycle.

## Motivation

Phase 9 of the ROADMAP adds multi-physics coupling (FSI, CHT,
reactive flow). We already have a one-adapter-per-case model
([RFC 0002](./0002-adapter-contract.md)) and a workflow DAG for
sequencing single-physics steps ([ARCHITECTURE.md § 6]). Neither
handles *concurrent* physics advancing on their own time-steppers,
exchanging fields at a shared interface, with convergence control
driving the whole group.

Ad-hoc scripts written per coupling scenario are what every open
simulation stack already has — and they're exactly what Valenx
promises to replace. A dedicated `CouplingAdapter` surface makes the
orchestrator pluggable, keeps participant adapters oblivious to each
other, and gives the UI one consistent place to show interface
residuals and per-participant progress.

## Design

### The two adapter roles

```rust
pub trait Adapter { /* unchanged */ }

/// A meta-adapter that runs a set of participating adapters in
/// lockstep against a coupling config.
pub trait CouplingAdapter: Adapter {
    /// Participants this coupling wires up. Each maps to a concrete
    /// `Adapter` instance that lives in `AdapterRegistry`.
    fn participants(&self, case: &Case) -> Result<Vec<Participant>, AdapterError>;

    /// Prepare the coupling configuration file (`precice-config.xml`
    /// today) as well as each participant's own workdir.
    fn prepare_coupling(
        &self,
        case: &Case,
        workdir: &Path,
    ) -> Result<CouplingLayout, AdapterError>;

    /// Drive the coupled run: spawn each participant, relay
    /// residuals + logs back to `RunContext`, honour cancellation.
    fn run_coupling(
        &self,
        layout: &CouplingLayout,
        ctx: &mut RunContext,
    ) -> Result<RunReport, AdapterError>;
}

pub struct Participant {
    /// Adapter ID (matches `AdapterInfo::id`).
    pub adapter_id: String,
    /// Case to run under that adapter (typically a sub-case of the
    /// top-level case).
    pub case: Case,
    /// Data exchanged outward. Names are coupling-scheme-local
    /// (e.g. `"Stress"`, `"Displacement"`).
    pub provides: Vec<String>,
    /// Data consumed inward.
    pub requires: Vec<String>,
}

pub struct CouplingLayout {
    pub workdir: PathBuf,
    pub config_path: PathBuf,
    pub participants: Vec<PreparedParticipant>,
}

pub struct PreparedParticipant {
    pub id: String,
    pub job: PreparedJob,
}
```

### Control flow

1. User assembles a coupled case: CFD + FEA + interface mappings.
2. The workflow DAG records one coupling node with the participants
   as sub-children. The node's input / output ports declare the
   exchanged field types.
3. On run, the engine resolves the coupling adapter, calls
   `prepare_coupling` (which calls each participant's `prepare` under
   the hood). Each participant gets a subdirectory under the coupling
   workdir.
4. `run_coupling` spawns the participants (either as subprocesses, or
   as long-running handles when the participant adapter supports it).
   It streams per-participant residuals + the coupling's own
   iteration residual up to `RunContext.progress` / `log`.
5. On completion, each participant's `collect` is called; the
   coupling adapter merges the individual `Results` into one combined
   `Results` with an explicit interface-mapping Provenance record.

### Data mapping model

Valenx's canonical `Mesh` + `Field` already carry the information the
mapper needs:

- `Field::location` tells the orchestrator whether the field lives on
  nodes / cells / faces.
- `Mesh::boundaries` names the surface patches involved.
- `Units` tag every field so the mapper fails loudly if a participant
  tries to hand "pressure in Pa" to one that expects "pressure in
  psi".

The RFC does **not** standardise the low-level mapper algorithms
(nearest-neighbour, RBF, consistent / conservative) — that's the
orchestrator's job. Valenx's contract is about what the orchestrator
sees and what participants see, not what the orchestrator does
internally.

### UI surface

- Coupling is modelled in the browser tree as one node that expands
  into its participants.
- The viewport overlays the exchange interface as a highlighted
  patch with a unit badge ("Pa" / "N·m⁻²") so misconfigurations are
  visible before running.
- The timeline panel shows one row per participant plus one row for
  the coupling iteration loop.
- The run command palette entry is `coupled-solve`, which lists all
  registered `CouplingAdapter` implementations.

### Open questions

1. **Steady-steady coupling** (e.g. CFD steady feeding Aster static):
   should the orchestrator iterate or converge once? Current plan:
   treat as a degenerate case of transient coupling with one
   coupling window.
2. **Mesh evolution during coupling** (e.g. FSI with morphing
   meshes): does the orchestrator see a new `Mesh` each iteration,
   or does it treat remeshing as a separate opaque step? Current
   plan: the participant adapter owns remeshing; it reports a new
   `Mesh` hash in its provenance and the coupler re-derives the
   mapping.
3. **Participant crashes**: the current `AdapterError` taxonomy is
   per-participant. The coupling layer needs a higher-level
   `CouplingError` variant that points at which participant failed
   and preserves its original error.

## Drawbacks

- Adds a second adapter surface, which means another test matrix.
- preCICE is today's only realistic orchestrator; the abstraction
  may be over-engineered for a one-off.

## Alternatives

- Leave coupling to a dedicated `valenx-coupling` crate outside the
  adapter pattern. Loses the registry integration (status colours,
  probe, license-mode tracking).
- Bake coupling directly into `valenx-core::workflow`. Conflates the
  DAG (sequence one-at-a-time) with the concurrent-exchange case.

## Prior art

- preCICE's own tutorials.
- Salome-Meca's YACS for coupled Code_Aster + OpenFOAM.
- Kratos Multiphysics' in-process multi-physics.

## Unresolved

- Naming: `CouplingAdapter` vs `Orchestrator` vs `Composer`?
- Should the coupling's `Results` include per-participant `Results`
  as ancestors in the provenance graph? (Leaning yes.)
