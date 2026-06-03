# RFC 0010 — preCICE concurrent participant orchestration

| | |
|---|---|
| **Status** | Draft |
| **Phase target** | 9 (multi-physics coupling) |
| **Type** | Implementation |
| **Author** | Maintainers |
| **Discussion** | open in PR — not yet merged |
| **Supersedes** | — |
| **Related** | RFC 0007 (coupling adapters), RFC 0009 (HPC executors) |

## Summary

Today the preCICE meta-adapter (`valenx-adapter-precice`) parses
`[coupling]` + `[[coupling.participant]]` blocks, stages each
participant's case directory + `precice-config.xml` into the workdir,
and runs `precice-tools check` for config validation. That's the
"compile-time" half of coupling.

This RFC defines the runtime half: launch every participant solver
**concurrently** against the shared coupling interface, stream their
logs into a unified UI feed, watch for the coupling iteration
exchanges in `precice-config.log`, and present a single
`RunReport` summarising all participants.

## Motivation

A coupled simulation isn't a coupled simulation until the
participants are actually running. Today users have to invoke each
participant manually in separate terminals after the staging step.
That's:

- Error-prone (forgetting to start one participant blocks the
  others on the coupling interface and they never time out cleanly).
- Invisible to the UI (no progress bar, no log feed, no convergence
  view across the coupled iteration).
- Impossible to cancel atomically (Ctrl-C in one terminal leaves
  the others orphaned with shared-memory locks held).

To replace commercial coupling tools (preCICE+ANSYS, MpCCI), we need
the runtime story.

## Design

### Adapter contract additions

```rust
impl Adapter for PreciceAdapter {
    // existing prepare() unchanged

    fn run(&self, job: &PreparedJob, ctx: &mut RunContext)
        -> Result<RunReport, AdapterError>
    {
        // 1. Read valenx_coupling.json from workdir (the manifest
        //    prepare() emitted listing every participant's
        //    adapter_id + case_dir + native_command).
        // 2. For each participant, look up its Adapter from the
        //    registry-clone the prepare() step pre-baked into the
        //    manifest.
        // 3. Spawn every participant on a thread:
        //       adapter.run(participant_job, &mut child_ctx)
        //    where child_ctx forwards progress/log into a shared
        //    fan-in channel.
        // 4. Watch precice-config.log for the coupling-iteration
        //    markers ("Iteration N converged" / "Subcycling…").
        //    Map those to ctx.report_progress so the UI shows a
        //    single coupled-iteration counter instead of N
        //    independent solver progress bars.
        // 5. Wait on every participant. RunReport carries:
        //       - aggregate exit code (max of children, 0 if all 0)
        //       - max wall_time
        //       - converged: Some(true) iff every participant
        //         reported converged AND coupling reached the
        //         configured time-window count
        //       - residual_history: empty (per-participant
        //         residuals stay in their own RunHandles)
        //       - warnings: union of children's
        // 6. cancel-token: when ctx.cancel.fire(), forward to
        //    every child's RunContext. Wait up to 5s for graceful
        //    shutdown, then SIGKILL.
    }
}
```

### Spawning model

Every participant is a regular Adapter::run() call on a worker
thread. The fan-in channel uses `std::sync::mpsc::channel` with
each child's `ProgressSink` / `LogSink` cloning the Sender. The
parent thread:

1. Spawns N children.
2. Drains the fan-in channel until all children have signalled
   completion (or cancellation).
3. Joins every child thread.

### Log multiplexing

Each child's log lines are tagged with the participant id
(`[fluid] Solving for Ux…`) before being forwarded to `ctx.log`.
The UI's existing log-panel filter machinery handles the rest;
no new UI surface needed.

### Coupling-iteration tracking

`precice-config.log` contains lines like:

```text
preCICE: Coupling iteration 1
preCICE: it 1 of timestep 0.001 converged
preCICE: timestep 1 of 100 done
```

A small parser sits in the parent thread reading the log file
between participant-progress events. Each `timestep N of M done`
becomes a `ctx.report_progress(N/M, "coupled timestep N of M")`.

### Cancellation semantics

A single `CancellationToken` shared across all children. Parent
flips it on UI cancel; each child's `subprocess::run` honours it
on its next `check_cancel()` poll. After `cancel.fire()`, the
parent waits 5 seconds for children to exit cleanly, then sends
SIGTERM to each. After another 2 seconds, SIGKILL.

## Drawbacks

- Adds ~300 lines of orchestration code to the preCICE adapter.
- Log fan-in creates a single point of contention — a chatty
  participant can starve quieter ones' updates if the channel
  isn't drained fast enough. Mitigation: bounded channel + explicit
  per-participant rate limit.
- The "converged" semantic for a coupled run is ambiguous when
  participants have different convergence definitions (steady RANS
  + transient FEA). The proposal — `Some(true)` only if every
  participant agrees — is conservative. Some users may want
  "the coupling iteration converged regardless of inner solvers."
  Cover via a per-coupling config flag.

## Alternatives considered

- **Run each participant in its own valenx-app run pipeline,
  serialised** — simpler but defeats the whole point of coupling
  (participants exchange data at every coupling timestep).
- **Use preCICE's `precice-run` Bash launcher** — works but
  bypasses the registry / probe / status-badge story and gives the
  UI no progress hooks.
- **Spawn participants as separate OS processes via `mpiexec
  -n 1 …`** — preCICE supports this and it's how their official
  examples do it. Considered, but it would require us to know how
  each adapter wants to be launched in MPI mode, which is
  adapter-specific.

## Migration path

Phase 9.1 (this RFC):
- Implement `PreciceAdapter::run()` per the design above.
- Update `valenx_coupling.json` schema to include enough info per
  participant for the parent to look up its Adapter (it already has
  adapter_id + case_path; it needs the prepared command list too —
  which is what `prepare_participant()` should now persist).

Phase 9.2: error-mode UX. What happens when one participant fails
mid-run? Today's design: cancel all others, surface the failed
participant's error verbatim. Edge cases (participant exited zero
but produced no output, participant deadlocked on the coupling
interface) need their own design pass.

Phase 9.3: per-participant probe. The registry's status badge
should show one row per participant in the case browser (currently
the row only reflects the meta-adapter `precice` itself).

## Open questions

1. **Where does the parent run?** Same machine as participants, or
   on the user's laptop coordinating remote runs (Phase 11 HPC
   integration)? The `Executor` trait from RFC 0009 makes this
   per-participant — needs explicit interaction design.
2. **Mesh staging.** Each participant needs the mesh of the
   surfaces it shares with others. preCICE handles partitioning
   but the file-staging layer is still ours. RFC 0011 (file-
   staging) is implied but not yet drafted.
3. **Live coupling residuals.** preCICE writes coupling iteration
   residuals to its own log; they don't fit the per-solver residual
   chart. Either add a second residual axis or split the chart by
   participant.

## References

- RFC 0007 — coupling adapters (defined the meta-adapter shape).
- RFC 0009 — HPC executors (interacts with this when participants
  run on different machines).
- preCICE docs — `Coupling Schemes` and `Logging` sections.
- preCICE Python tutorials, particularly the Fluent + Code_Aster
  FSI demo which is what we'll reference-test against.
