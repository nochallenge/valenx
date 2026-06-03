# RFC 0009 — HPC job submission

| | |
|---|---|
| **Status** | Draft |
| **Phase target** | 11 (HPC + cluster execution) |
| **Type** | Architecture |
| **Author** | Maintainers |
| **Discussion** | open in PR — not yet merged |

## Summary

Define a small, transport-agnostic adapter contract for submitting
prepared cases to remote schedulers (SLURM, PBS, LSF, kubernetes
Jobs, plain SSH-and-nohup). Local execution stays the default; HPC
becomes a sidecar that wraps the existing `Adapter::run()` lifecycle
without forking it.

This RFC is intentionally a scaffold — it sketches the contract and
two reference adapters, but doesn't try to spec the full HPC story
(file staging, result fetch-back, cancellation across the wire,
multi-tenant credentials) which need their own follow-ups.

## Motivation

Today every `Adapter::run()` shells out via
`valenx_core::subprocess::run`, which spawns a child on the local
machine. That's fine for the workflow loop (the user's laptop runs
the solver). It does not work for:

- Cases too large for the user's hardware (a 50M-cell airfoil at LES,
  a 100k-DOF nonlinear FEA, a multi-day battery cycling study).
- Organisations whose policy is "no solver runs on the desktop —
  everything goes to the cluster."
- Reproducibility expectations (cluster runs with pinned MPI ranks,
  fixed CPU pinning, captured environment).

Without an HPC story, Valenx can never replace its commercial
counterparts (Star-CCM+, ANSYS, Abaqus) for any organisation past a
single workstation.

## Design

### Two-tier adapter shape

Today:

```rust
pub trait Adapter {
    fn info(&self) -> AdapterInfo;
    fn probe(&self) -> Result<ProbeReport, AdapterError>;
    fn prepare(&self, case: &Case, workdir: &Path) -> Result<PreparedJob, AdapterError>;
    fn run(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunReport, AdapterError>;
    fn collect(&self, job: &PreparedJob) -> Result<Results, AdapterError>;
    fn capabilities(&self) -> Capabilities;
}
```

After:

```rust
pub trait Adapter { /* unchanged */ }

/// Where an Adapter's `run()` actually executes.
pub trait Executor {
    fn id(&self) -> &str;
    fn submit(&self, job: &PreparedJob, ctx: &mut RunContext) -> Result<RunHandle, ExecutorError>;
    fn poll(&self, handle: &RunHandle) -> Result<RunStatus, ExecutorError>;
    fn cancel(&self, handle: &RunHandle) -> Result<(), ExecutorError>;
    fn fetch_results(&self, handle: &RunHandle, into: &Path) -> Result<(), ExecutorError>;
}
```

The user picks an `Executor` per project (Settings → "Run on…");
the same `Adapter` works against any `Executor`. `LocalExecutor`
(the existing `subprocess::run` flow) is the default; `SlurmExecutor`
/ `PbsExecutor` / `K8sExecutor` are the first remote ones.

### Reference executors (in scope for this RFC)

#### `LocalExecutor` (no behaviour change)

Wraps the current `subprocess::run` path so existing adapters continue
to work without code changes.

#### `SlurmExecutor`

- `submit()`: writes a `submit.sh` script alongside the prepared
  workdir, calls `sbatch submit.sh`, parses the `Submitted batch
  job <N>` reply.
- `poll()`: `squeue -j <N> -h -o %T` → maps `PENDING/RUNNING/COMPLETED/
  FAILED/CANCELLED/TIMEOUT` to `RunStatus`.
- `cancel()`: `scancel <N>`.
- `fetch_results()`: `rsync` (or `scp`, opt-in) the workdir back to
  the local app.
- Per-case `[hpc.slurm]` block in `case.toml` carries `partition`,
  `time_limit`, `nodes`, `cpus_per_task`, `mem_per_cpu`, `account`,
  `qos`. Defaults if omitted: `partition=default`, `time=01:00:00`,
  `nodes=1`, `cpus_per_task=1`.

### Out of scope (follow-up RFCs)

- File staging for huge meshes (the workdir might be 10 GB; we
  shouldn't naively rsync the whole tree every poll).
- Multi-tenant credentials and SSH key management.
- Live log streaming over the wire (today the local executor pipes
  stdout directly; remote SSH-with-tail is more involved).
- Restart/checkpoint integration (resume an interrupted SLURM job).
- GPU resource declaration (`gres=gpu:1`, MIG slices, etc.).
- Hybrid pipelines (run gmsh locally, OpenFOAM on the cluster, post-
  process locally).

## Drawbacks

- A second trait layer adds a real surface area to the canonical
  `Adapter` story. RFC 0002 deliberately kept the adapter contract
  small; this widens it.
- Cancellation semantics are different across schedulers. `scancel`
  is fast; some k8s reconcilers can take seconds.
- "Failed" and "completed-with-non-zero-exit" are different states
  for some schedulers and the same for others. We need a careful
  mapping.

## Alternatives considered

- **Embed scheduler logic into each adapter** — every adapter learns
  to talk to SLURM separately. Rejected: 11 live adapters × 5
  schedulers = 55 implementations, all bound to drift.
- **Wrap the Adapter trait** (`SlurmAdapter<A: Adapter>`) — works
  but composes badly with the existing Arc-based registry; would
  require generic-over-Adapter all the way through `valenx-core`.
- **Pre-build SLURM scripts in `prepare()`, run with `sbatch` from
  the local executor** — confuses "where does the work happen" with
  "what tool does the work." Better to keep them orthogonal.

## Migration path

Phase 11.0 (this RFC):
- Land the `Executor` trait in `valenx-core`.
- `LocalExecutor` is the default and only ships in `valenx-core`.
- `SlurmExecutor` lives in a new crate `valenx-executor-slurm`.

Phase 11.1: PBS / LSF executors as additional crates.

Phase 11.2: file-staging RFC.

Phase 11.3: kubernetes Job executor.

Pre-existing single-executor adapters keep working: the trait is
additive, default-implemented to "use LocalExecutor."

## Open questions

1. **RunHandle shape** — the existing `valenx-app::run::RunHandle`
   carries a thread JoinHandle and an mpsc Receiver. Remote executors
   have neither. Either we add an enum variant or we redesign
   `RunHandle` to be a small ID + state-bag the executor opaquely
   knows how to interpret.
2. **Where does the workdir live remotely?** — the cluster's
   filesystem is not the user's filesystem. The contract probably
   needs a `RemoteWorkdir` newtype distinguishing "the path on this
   machine" from "the path on the executor's machine."
3. **Live progress reporting from remote runs** — today
   `RunContext::progress` and `RunContext::log` route directly to the
   UI thread. With a remote executor, the progress comes back over
   the wire (rsync of log.simpleFoam every N seconds? streaming
   tail-ssh? webhook?). Hard to specify generically; might want
   per-executor strategies.

## References

- RFC 0002 — current `Adapter` contract.
- RFC 0007 — coupling adapters (related: multi-process orchestration).
- SLURM `sbatch(1)` / `squeue(1)` man pages.
- Kubernetes Jobs API.
- HTCondor's `condor_submit` (similar shape; not in scope today
  but the contract should accommodate it).
