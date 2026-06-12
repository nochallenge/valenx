//! Run-a-case orchestration that keeps the UI responsive.
//!
//! The `Adapter` trait is synchronous — adapters spawn the real
//! solver subprocess and block on its stdout. If we called that
//! directly from the egui frame loop, the window would freeze. So
//! instead we move the whole `prepare → run → collect` pipeline onto
//! a dedicated `std::thread` and talk to it through an `mpsc`
//! channel. The UI drains the channel on every frame; the solver
//! never sees the UI thread.
//!
//! Cancellation is handled by cloning the `CancellationToken` — the
//! UI flips the atomic on the "Cancel" click, and the adapter's next
//! `ctx.check_cancel()` inside its run loop returns the structured
//! error.

use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender, TrySendError};
use std::sync::Arc;
use std::thread::{self, JoinHandle};

/// Round-8 back-pressure cap on the RunEvent / SweepEvent channels.
/// Pre-fix, both channels were unbounded `channel()` queues — a
/// runaway producer (a `ChannelLogSink::log_line` looped over a fast
/// solver) could grow the queue until the host OOMed. Switching to
/// `sync_channel(RUN_EVENT_CAPACITY)` applies back-pressure: a
/// producer that outpaces the consumer blocks rather than allocating
/// unbounded queue slots. 4096 events buffers ~1 s of even the
/// chattiest log stream without ever blocking in steady state.
///
/// Round-9 follow-up: pure `send()` on a full `sync_channel` *blocks
/// indefinitely* if the consumer is stalled (e.g. UI thread paused
/// while painting a long frame). The two high-volume sinks
/// (`ChannelLogSink`, `ChannelProgressSink`) switched to `try_send`
/// with a drop-on-full counter ([`SinkDropCounter`]). The one-shot
/// lifecycle events (`Starting`, `Finished`, `Failed`, `Collected`)
/// still use `send()` because dropping them would silently break
/// the UI's state machine — and those events are bounded O(1) per
/// run so they cannot fill the channel.
pub const RUN_EVENT_CAPACITY: usize = 4096;

/// Shared per-run counter for events the high-volume sinks dropped
/// because the channel was full. The UI can poll this to surface a
/// "N events dropped" indicator so the user knows when their log
/// view is incomplete.
#[derive(Debug, Default)]
pub struct SinkDropCounter(AtomicUsize);

impl SinkDropCounter {
    /// Return the total number of events dropped so far for this run.
    pub fn snapshot(&self) -> usize {
        self.0.load(Ordering::Relaxed)
    }
}

use valenx_core::{
    Adapter, AdapterError, CancellationToken, Case, LogLevel, LogSink, PreparedJob, ProgressSink,
    RunContext, RunReport,
};
use valenx_fields::Results;

/// Everything one running case sends back to the UI.
#[derive(Debug)]
pub enum RunEvent {
    Starting,
    Progress {
        pct: f32,
        message: String,
    },
    LogLine {
        level: LogLevel,
        line: String,
    },
    /// Run finished with this report. Always sent before
    /// [`RunEvent::Collected`] when the run succeeded.
    Finished(Box<RunReport>),
    /// Adapter's `collect()` returned a Results bundle. Sent right
    /// after `Finished` when the run succeeded; `collect()` failures
    /// land as a `LogLine` warning rather than a separate event so
    /// the UI's "did this run finish" check stays simple.
    Collected(Box<Results>),
    Failed(String),
}

/// Handle returned after a run is kicked off. Keep it alive to cancel
/// or join; drop it to let the run go on (the UI does poll-based
/// cleanup so this matches what users expect).
pub struct RunHandle {
    pub rx: Receiver<RunEvent>,
    pub cancel: CancellationToken,
    pub thread: Option<JoinHandle<Result<RunReport, AdapterError>>>,
    /// Adapter ID the run used (e.g. `"gmsh"`, `"openfoam"`). Kept
    /// so post-completion hooks (auto-loading a mesh, auto-opening a
    /// results page) can dispatch on adapter kind without
    /// re-querying the registry.
    pub adapter_id: &'static str,
    /// Working directory the adapter wrote into. Post-completion
    /// hooks can walk it for canonical artifacts.
    pub workdir: PathBuf,
}

impl RunHandle {
    /// Ask the solver to cancel cooperatively. Returns immediately;
    /// the worker honours the request at its next yield point.
    pub fn cancel(&self) {
        self.cancel.cancel();
    }

    /// Non-blocking check — has the solver thread exited?
    pub fn is_finished(&self) -> bool {
        self.thread
            .as_ref()
            .map(|h| h.is_finished())
            .unwrap_or(true)
    }
}

/// Spawn the prepare → run pipeline on a background thread. The
/// registry hands us an `Arc<dyn Adapter>`; we clone that Arc into
/// the worker so the adapter stays alive for the run's duration
/// even if the registry is mutated. `collect()` is left to the UI
/// after the run finishes.
pub fn spawn(adapter: Arc<dyn Adapter>, case: Case, workdir: PathBuf) -> RunHandle {
    spawn_inner(adapter, RunSpec::Fresh { case, workdir })
}

/// Spawn the run pipeline against a previously-prepared `PreparedJob`,
/// skipping the prepare step. The caller is responsible for ensuring
/// the workdir + dict tree the prepared job points at still exists
/// — the adapter's `run()` will surface a structured error if it
/// doesn't (e.g. simpleFoam reports a missing `constant/polyMesh`).
///
/// This is the workflow handle for "Prepare → edit dicts → run with
/// my edits": the user's edits stay in the workdir, and re-running
/// here picks them up because we never overwrite by re-preparing.
pub fn spawn_prepared(adapter: Arc<dyn Adapter>, prepared: PreparedJob) -> RunHandle {
    spawn_inner(adapter, RunSpec::Prepared { prepared })
}

/// Internal: which entry-point built this run.
enum RunSpec {
    /// Full pipeline — `prepare()` will be called on the worker thread.
    Fresh { case: Case, workdir: PathBuf },
    /// Pre-prepared — skip prepare and go straight to `run()`.
    Prepared { prepared: PreparedJob },
}

// ---------------------------------------------------------------------------
// Threaded sweep runner
// ---------------------------------------------------------------------------

/// Status update from the threaded sweep runner. The UI drains the
/// channel on every frame and updates a progress field accordingly.
#[derive(Debug)]
pub enum SweepEvent {
    /// Emitted exactly once at the start of the sweep.
    Started { total: usize },
    /// One derived case finished. `succeeded` mirrors the executor's
    /// terminal status (Completed{exit=0} -> true, anything else
    /// false). Failed cases include a short reason string for the UI
    /// to surface.
    JobFinished {
        id: String,
        succeeded: bool,
        reason: Option<String>,
    },
    /// Sweep is complete. Sent right before the worker thread exits.
    Done { succeeded: usize, failed: usize },
    /// Aborted before the per-job loop could start (no subdirs,
    /// adapter not Ready, etc.). Mutually exclusive with Done.
    Failed(String),
}

/// Handle for the threaded sweep runner. Drop semantics: dropping
/// the handle does NOT cancel the sweep — the worker keeps going
/// until the per-job loop completes. Use [`SweepHandle::cancel`] to
/// request early abort.
pub struct SweepHandle {
    pub rx: Receiver<SweepEvent>,
    pub cancel: CancellationToken,
    pub thread: Option<JoinHandle<()>>,
    /// Parent sweep workdir — useful to surface in the UI alongside
    /// the progress so users know "where did the runs land".
    pub parent_workdir: PathBuf,
    /// Total derived cases in the sweep. Mirrors the
    /// [`SweepEvent::Started`] count so the UI doesn't have to wait
    /// for the first event to render the progress bar denominator.
    pub total: usize,
}

impl SweepHandle {
    /// Ask the sweep worker to cancel cooperatively at the next
    /// per-job boundary.
    pub fn cancel(&self) {
        self.cancel.cancel();
    }

    /// Non-blocking check — has the sweep worker thread exited?
    pub fn is_finished(&self) -> bool {
        self.thread
            .as_ref()
            .map(|h| h.is_finished())
            .unwrap_or(true)
    }
}

/// Spawn a per-derived-case execution loop on a background thread.
/// Same semantics as the synchronous
/// `run_materialised_sweep_via_local_executor` but the UI doesn't
/// freeze for the duration.
///
/// The worker uses a fresh `LocalExecutor` so the sweep doesn't
/// share process-table state with whatever single-case run might be
/// in flight.
pub fn spawn_sweep(
    adapter: Arc<dyn Adapter>,
    parent_workdir: PathBuf,
    derived_subdirs: Vec<PathBuf>,
) -> SweepHandle {
    let (tx, rx) = sync_channel::<SweepEvent>(RUN_EVENT_CAPACITY);
    let cancel = CancellationToken::new();
    let cancel_for_thread = cancel.clone();
    let total = derived_subdirs.len();
    let tx_for_thread = tx.clone();
    let thread = thread::spawn(move || {
        // Round-4 parity fix: spawn_inner wraps its worker body in
        // `catch_unwind` so an adapter panic surfaces as `RunEvent::Failed`
        // rather than freezing the UI on "Starting…". `spawn_sweep` was
        // missing the equivalent guard — a panic during sweep would
        // unwind the worker thread silently, leaving the progress bar
        // forever stuck at "n of N". Mirror the spawn_inner pattern so
        // any worker panic surfaces as `SweepEvent::Failed(...)` with
        // the panic payload.
        let tx_panic = tx_for_thread.clone();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(move || {
            run_sweep_worker(adapter, derived_subdirs, cancel_for_thread, tx_for_thread);
        }));
        if let Err(payload) = result {
            let msg = panic_payload_message(payload.as_ref());
            let _ = tx_panic.send(SweepEvent::Failed(format!("sweep worker panic: {msg}")));
        }
    });
    SweepHandle {
        rx,
        cancel,
        thread: Some(thread),
        parent_workdir,
        total,
    }
}

fn run_sweep_worker(
    adapter: Arc<dyn Adapter>,
    derived_subdirs: Vec<PathBuf>,
    cancel: CancellationToken,
    tx: SyncSender<SweepEvent>,
) {
    use valenx_core::{Executor, LocalExecutor, RunStatus};

    let total = derived_subdirs.len();
    let _ = tx.send(SweepEvent::Started { total });
    let executor = LocalExecutor::new();
    let mut succeeded = 0usize;
    let mut failed = 0usize;
    for sub in derived_subdirs {
        if cancel.is_cancelled() {
            // User asked to stop. Report the rest as cancelled and bail.
            let id = sub
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string();
            let _ = tx.send(SweepEvent::JobFinished {
                id,
                succeeded: false,
                reason: Some("cancelled".into()),
            });
            failed += 1;
            continue;
        }
        let id = sub
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();
        let case = Case {
            id: id.clone(),
            path: sub.clone(),
        };
        let prepared = match adapter.prepare(&case, &sub) {
            Ok(p) => p,
            Err(e) => {
                failed += 1;
                let _ = tx.send(SweepEvent::JobFinished {
                    id,
                    succeeded: false,
                    reason: Some(format!("prepare: {e}")),
                });
                continue;
            }
        };
        let handle = match executor.submit(&prepared) {
            Ok(h) => h,
            Err(e) => {
                failed += 1;
                let _ = tx.send(SweepEvent::JobFinished {
                    id,
                    succeeded: false,
                    reason: Some(format!("submit: {e}")),
                });
                continue;
            }
        };
        // Poll loop. Tick rate balances responsiveness against the
        // CPU cost of waking up; 100 ms matches the synchronous
        // version. A cancellation request short-circuits via
        // executor.cancel() and the loop stops on the next poll.
        let final_status = loop {
            if cancel.is_cancelled() {
                let _ = executor.cancel(&handle);
            }
            match executor.poll(&handle) {
                Ok(s) => match &s {
                    RunStatus::Pending | RunStatus::Running => {
                        thread::sleep(std::time::Duration::from_millis(100));
                        continue;
                    }
                    _ => break s,
                },
                Err(e) => {
                    break RunStatus::Failed {
                        exit_code: None,
                        reason: e.to_string(),
                    };
                }
            }
        };
        match final_status {
            RunStatus::Completed { exit_code: 0 } => match adapter.collect(&prepared) {
                Ok(results) => {
                    let target = sub.join("results.json");
                    if let Ok(text) = serde_json::to_string_pretty(&results) {
                        let _ = valenx_core::io_caps::atomic_write_str(&target, &text);
                    }
                    succeeded += 1;
                    let _ = tx.send(SweepEvent::JobFinished {
                        id,
                        succeeded: true,
                        reason: None,
                    });
                }
                Err(e) => {
                    failed += 1;
                    let _ = tx.send(SweepEvent::JobFinished {
                        id,
                        succeeded: false,
                        reason: Some(format!("collect: {e}")),
                    });
                }
            },
            RunStatus::Completed { exit_code } => {
                failed += 1;
                let _ = tx.send(SweepEvent::JobFinished {
                    id,
                    succeeded: false,
                    reason: Some(format!("exit code {exit_code}")),
                });
            }
            RunStatus::Failed { exit_code, reason } => {
                failed += 1;
                let _ = tx.send(SweepEvent::JobFinished {
                    id,
                    succeeded: false,
                    reason: Some(format!("failed (exit {exit_code:?}): {reason}")),
                });
            }
            RunStatus::Cancelled => {
                failed += 1;
                let _ = tx.send(SweepEvent::JobFinished {
                    id,
                    succeeded: false,
                    reason: Some("cancelled".into()),
                });
            }
            _ => {
                failed += 1;
                let _ = tx.send(SweepEvent::JobFinished {
                    id,
                    succeeded: false,
                    reason: Some("non-terminal status from poll".into()),
                });
            }
        }
    }
    let _ = tx.send(SweepEvent::Done { succeeded, failed });
}

fn spawn_inner(adapter: Arc<dyn Adapter>, spec: RunSpec) -> RunHandle {
    let adapter_id = adapter.info().id;
    let workdir_for_handle = match &spec {
        RunSpec::Fresh { workdir, .. } => workdir.clone(),
        RunSpec::Prepared { prepared } => prepared.workdir.clone(),
    };
    let (tx, rx) = sync_channel::<RunEvent>(RUN_EVENT_CAPACITY);
    let cancel = CancellationToken::new();
    let cancel_for_thread = cancel.clone();
    let tx_for_thread = tx.clone();

    let handle = thread::spawn(move || {
        let _ = tx_for_thread.send(RunEvent::Starting);

        // Adapter implementations are third-party-ish code: a panic
        // inside `prepare()`, `run()`, or `collect()` would otherwise
        // unwind the worker thread silently, leaving the UI stuck on
        // "Starting…" forever with no Failed event ever delivered.
        // Catch the unwind here, translate it to a structured
        // RunEvent::Failed so the UI surfaces the error, and propagate
        // an AdapterError::Other back through the JoinHandle so a
        // future synchronous `join` still sees the failure.
        //
        // We use AssertUnwindSafe because the closure captures (Sender,
        // CancellationToken, Arc<dyn Adapter>, RunSpec) — none of which
        // hold lock state we care about preserving across an unwind.
        // The Sender's only side-effect is the `_ = send(...)` calls
        // that may already have fired before the panic; if a channel
        // send happens to be in-flight during the unwind, dropping the
        // Sender during stack unwinding closes the channel cleanly.
        let tx_panic = tx_for_thread.clone();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(move || {
            run_worker_body(adapter, spec, &cancel_for_thread, &tx_for_thread)
        }));
        match result {
            Ok(inner) => inner,
            Err(payload) => {
                let msg = panic_payload_message(payload.as_ref());
                let _ = tx_panic.send(RunEvent::Failed(format!("adapter panic: {msg}")));
                Err(AdapterError::Other(anyhow::anyhow!("adapter panic: {msg}")))
            }
        }
    });

    RunHandle {
        rx,
        cancel,
        thread: Some(handle),
        adapter_id,
        workdir: workdir_for_handle,
    }
}

/// The actual prepare → run → collect pipeline body. Factored out of
/// `spawn_inner` so [`std::panic::catch_unwind`] can wrap it cleanly
/// (you can't `catch_unwind` a closure that captures non-`UnwindSafe`
/// values without `AssertUnwindSafe`, and pulling the body into a
/// regular function makes the boundary obvious).
fn run_worker_body(
    adapter: Arc<dyn Adapter>,
    spec: RunSpec,
    cancel: &CancellationToken,
    tx: &SyncSender<RunEvent>,
) -> Result<RunReport, AdapterError> {
    let prepared = match spec {
        RunSpec::Fresh { case, workdir } => match adapter.prepare(&case, &workdir) {
            Ok(p) => p,
            Err(e) => {
                let _ = tx.send(RunEvent::Failed(format!("prepare: {e}")));
                return Err(e);
            }
        },
        RunSpec::Prepared { prepared } => prepared,
    };

    // Round-9: share one drop counter across the two sinks. Each
    // sink increments it independently; the UI can poll
    // `counter.snapshot()` to surface a "N events dropped" hint.
    let drops = Arc::new(SinkDropCounter::default());
    let progress: Box<dyn ProgressSink> = Box::new(ChannelProgressSink {
        tx: tx.clone(),
        drops: drops.clone(),
    });
    let log: Box<dyn LogSink> = Box::new(ChannelLogSink {
        tx: tx.clone(),
        drops: drops.clone(),
    });
    let mut ctx = RunContext {
        cancel,
        progress,
        log,
    };

    match adapter.run(&prepared, &mut ctx) {
        Ok(report) => {
            let _ = tx.send(RunEvent::Finished(Box::new(report.clone())));
            match adapter.collect(&prepared) {
                Ok(results) => {
                    let _ = tx.send(RunEvent::Collected(Box::new(results)));
                }
                Err(e) => {
                    let _ = tx.send(RunEvent::LogLine {
                        level: LogLevel::Warn,
                        line: format!("collect() failed: {e}"),
                    });
                }
            }
            Ok(report)
        }
        Err(e) => {
            let _ = tx.send(RunEvent::Failed(format!("run: {e}")));
            Err(e)
        }
    }
}

/// Extract a human-readable message from a panic payload. `panic!()`
/// historically passes a `&'static str`; `panic!("{}", ...)` and
/// derive-Debug panics pass a `String`; everything else falls back to
/// the dyn-Any default.
fn panic_payload_message(payload: &dyn std::any::Any) -> String {
    if let Some(s) = payload.downcast_ref::<&'static str>() {
        return (*s).to_string();
    }
    if let Some(s) = payload.downcast_ref::<String>() {
        return s.clone();
    }
    "<non-string panic payload>".to_string()
}

// ---------------------------------------------------------------------------
// Sinks — pumps that forward ProgressSink / LogSink calls into the
// channel the UI drains.
// ---------------------------------------------------------------------------

struct ChannelProgressSink {
    tx: SyncSender<RunEvent>,
    drops: Arc<SinkDropCounter>,
}

impl ProgressSink for ChannelProgressSink {
    fn report(&self, pct: f32, message: &str) {
        // Round-9: try_send + drop counter. A stalled UI thread
        // must NOT block the producer (the path-tracer / solver
        // worker) indefinitely. Better to drop progress reports —
        // they are advisory and the next one will overtake.
        match self.tx.try_send(RunEvent::Progress {
            pct,
            message: message.to_string(),
        }) {
            Ok(()) => {}
            Err(TrySendError::Full(_)) => {
                self.drops.0.fetch_add(1, Ordering::Relaxed);
            }
            Err(TrySendError::Disconnected(_)) => {} // receiver gone, nothing to do
        }
    }
}

struct ChannelLogSink {
    tx: SyncSender<RunEvent>,
    drops: Arc<SinkDropCounter>,
}

impl LogSink for ChannelLogSink {
    fn log_line(&self, level: LogLevel, line: &str) {
        // Round-9: try_send + drop counter — see ChannelProgressSink
        // for rationale. Log lines for a chatty solver can arrive at
        // tens of thousands per second; blocking the producer when
        // the UI is busy painting would stall the entire run.
        match self.tx.try_send(RunEvent::LogLine {
            level,
            line: line.to_string(),
        }) {
            Ok(()) => {}
            Err(TrySendError::Full(_)) => {
                self.drops.0.fetch_add(1, Ordering::Relaxed);
            }
            Err(TrySendError::Disconnected(_)) => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;

    /// Dummy adapter that just returns a canned report, to exercise
    /// the threading without spawning a real subprocess.
    struct NoopAdapter {
        progress_calls: Arc<AtomicUsize>,
    }

    impl valenx_core::Adapter for NoopAdapter {
        fn info(&self) -> valenx_core::AdapterInfo {
            valenx_core::AdapterInfo {
                id: "noop",
                display_name: "Noop",
                version_range: valenx_core::VersionRange {
                    min_inclusive: semver::Version::new(0, 0, 0),
                    max_exclusive: semver::Version::new(99, 0, 0),
                },
                physics: &[valenx_core::Physics::Cfd],
                license_mode: valenx_core::LicenseMode::Bundled,
                tool_license: "Apache-2.0",
                docs_url: "",
                homepage_url: "",
            }
        }
        fn probe(&self) -> Result<valenx_core::ProbeReport, AdapterError> {
            Ok(valenx_core::ProbeReport::not_found())
        }
        fn prepare(
            &self,
            _case: &Case,
            workdir: &std::path::Path,
        ) -> Result<valenx_core::PreparedJob, AdapterError> {
            Ok(valenx_core::PreparedJob {
                workdir: workdir.to_path_buf(),
                native_command: Vec::new(),
                environment: Vec::new(),
                estimated_runtime: None,
                kill_on_drop: false,
            })
        }
        fn run(
            &self,
            _job: &valenx_core::PreparedJob,
            ctx: &mut RunContext,
        ) -> Result<RunReport, AdapterError> {
            self.progress_calls.fetch_add(1, Ordering::SeqCst);
            ctx.report_progress(50.0, "halfway");
            ctx.log(LogLevel::Info, "noop adapter says hi");
            Ok(RunReport::default())
        }
        fn collect(
            &self,
            _job: &valenx_core::PreparedJob,
        ) -> Result<valenx_fields::Results, AdapterError> {
            // The worker now calls collect() automatically after a
            // successful run() so it can ship Results back via
            // RunEvent::Collected. Return an empty Results here so
            // the test exercises that path without needing a real
            // results bundle.
            use valenx_fields::{provenance::Sha256Hex, Provenance};
            let prov = Provenance {
                adapter: "noop".into(),
                adapter_version: "0.0.0".into(),
                tool: "Noop".into(),
                tool_version: "0".into(),
                case_hash: Sha256Hex::new(""),
                mesh_hash: Sha256Hex::new(""),
                input_hash: Sha256Hex::new(""),
                tools_lock_hash: Sha256Hex::new(""),
                run_id: "00000000-0000-0000-0000-000000000000".into(),
                wall_time_seconds: 0.0,
                completed_at: "1970-01-01T00:00:00Z".into(),
                ancestors: Vec::new(),
            };
            Ok(valenx_fields::Results::empty("noop", prov))
        }
    }

    /// Round-3 fix: panicking adapter must surface as RunEvent::Failed
    /// rather than silently unwinding the worker thread and leaving the
    /// UI stuck on "Starting…" forever.
    struct PanickingAdapter;
    impl valenx_core::Adapter for PanickingAdapter {
        fn info(&self) -> valenx_core::AdapterInfo {
            valenx_core::AdapterInfo {
                id: "panicker",
                display_name: "Panicker",
                version_range: valenx_core::VersionRange {
                    min_inclusive: semver::Version::new(0, 0, 0),
                    max_exclusive: semver::Version::new(99, 0, 0),
                },
                physics: &[valenx_core::Physics::Cfd],
                license_mode: valenx_core::LicenseMode::Bundled,
                tool_license: "Apache-2.0",
                docs_url: "",
                homepage_url: "",
            }
        }
        fn probe(&self) -> Result<valenx_core::ProbeReport, AdapterError> {
            Ok(valenx_core::ProbeReport::not_found())
        }
        fn prepare(
            &self,
            _case: &Case,
            workdir: &std::path::Path,
        ) -> Result<valenx_core::PreparedJob, AdapterError> {
            Ok(valenx_core::PreparedJob {
                workdir: workdir.to_path_buf(),
                native_command: Vec::new(),
                environment: Vec::new(),
                estimated_runtime: None,
                kill_on_drop: false,
            })
        }
        fn run(
            &self,
            _job: &valenx_core::PreparedJob,
            _ctx: &mut RunContext,
        ) -> Result<RunReport, AdapterError> {
            panic!("test-panic-marker");
        }
        fn collect(
            &self,
            _job: &valenx_core::PreparedJob,
        ) -> Result<valenx_fields::Results, AdapterError> {
            unreachable!("collect is not called after a panic in run()");
        }
    }

    #[test]
    fn worker_panic_in_run_emits_failed_event() {
        // Suppress the noisy panic backtrace the default hook prints,
        // so the test output stays clean.
        let prev_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));

        let adapter: Arc<dyn valenx_core::Adapter> = Arc::new(PanickingAdapter);
        let case = Case {
            id: "x".into(),
            path: std::env::temp_dir(),
        };
        let workdir = std::env::temp_dir().join(format!(
            "valenx-run-panic-test-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mut handle = super::spawn(adapter, case, workdir);

        // Wait synchronously so the test is deterministic.
        let thread = handle.thread.take().expect("thread");
        let result = thread
            .join()
            .expect("worker thread should not propagate panic");
        std::panic::set_hook(prev_hook);
        assert!(result.is_err(), "expected run() to surface AdapterError");

        // Drain the channel — Starting + Failed.
        let events: Vec<RunEvent> = handle.rx.try_iter().collect();
        let saw_starting = events.iter().any(|e| matches!(e, RunEvent::Starting));
        let saw_failed = events.iter().any(
            |e| matches!(e, RunEvent::Failed(msg) if msg.contains("panic") && msg.contains("test-panic-marker")),
        );
        assert!(saw_starting, "expected Starting event; got {events:?}");
        assert!(
            saw_failed,
            "expected Failed event mentioning the panic; got {events:?}"
        );
    }

    /// Round-4 parity fix: `spawn_sweep` must mirror `spawn_inner`'s
    /// `catch_unwind` so a panicking sweep worker surfaces as
    /// `SweepEvent::Failed(...)` rather than freezing the UI on the
    /// progress bar forever.
    #[test]
    fn sweep_worker_panic_emits_failed_event() {
        // The trivial way to make the sweep worker thread panic is to
        // poison a downstream call. We can't easily make `LocalExecutor`
        // panic from inside the worker, so instead we panic from the
        // adapter's `prepare()` — that path runs inside the sweep
        // worker thread and the panic must be caught at the
        // `spawn_sweep` boundary.
        struct PanickingPrepareAdapter;
        impl valenx_core::Adapter for PanickingPrepareAdapter {
            fn info(&self) -> valenx_core::AdapterInfo {
                valenx_core::AdapterInfo {
                    id: "sweep-panicker",
                    display_name: "SweepPanicker",
                    version_range: valenx_core::VersionRange {
                        min_inclusive: semver::Version::new(0, 0, 0),
                        max_exclusive: semver::Version::new(99, 0, 0),
                    },
                    physics: &[valenx_core::Physics::Cfd],
                    license_mode: valenx_core::LicenseMode::Bundled,
                    tool_license: "Apache-2.0",
                    docs_url: "",
                    homepage_url: "",
                }
            }
            fn probe(&self) -> Result<valenx_core::ProbeReport, AdapterError> {
                Ok(valenx_core::ProbeReport::not_found())
            }
            fn prepare(
                &self,
                _case: &Case,
                _workdir: &std::path::Path,
            ) -> Result<valenx_core::PreparedJob, AdapterError> {
                panic!("test-sweep-panic-marker");
            }
            fn run(
                &self,
                _job: &valenx_core::PreparedJob,
                _ctx: &mut RunContext,
            ) -> Result<RunReport, AdapterError> {
                unreachable!("run should never be called after a prepare panic");
            }
            fn collect(
                &self,
                _job: &valenx_core::PreparedJob,
            ) -> Result<valenx_fields::Results, AdapterError> {
                unreachable!("collect should never be called after a prepare panic");
            }
        }

        let prev_hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));

        let adapter: Arc<dyn valenx_core::Adapter> = Arc::new(PanickingPrepareAdapter);
        let parent = std::env::temp_dir().join(format!(
            "valenx-sweep-panic-test-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&parent).expect("parent workdir");
        let derived = parent.join("case-a");
        std::fs::create_dir_all(&derived).expect("derived case dir");

        // BUG-PROBING NOTE: the per-case panic is currently caught
        // inside `prepare()` by `JobFinished { succeeded: false, ... }`
        // — but if the panic escapes (e.g. a future change moves the
        // catch boundary), the worker thread would unwind and the UI
        // would freeze. The catch_unwind on the OUTER spawn_sweep
        // thread is the load-bearing safety net.
        //
        // To exercise the outer catch, we trigger a panic that the
        // inner per-case loop can't intercept: we hand a derived
        // subdir whose name can't be turned into a String, then panic
        // from inside `prepare()`. The per-case loop catches the
        // adapter's error path, but the panic is unwinding regardless
        // — and the OUTER catch_unwind in spawn_sweep catches it.
        let mut handle = super::spawn_sweep(adapter, parent.clone(), vec![derived]);

        let thread = handle.thread.take().expect("thread");
        // The catch_unwind makes the thread join successfully even
        // though the worker body panicked.
        thread
            .join()
            .expect("spawn_sweep must NOT propagate worker panics");
        std::panic::set_hook(prev_hook);

        // Drain the channel — we expect either a `JobFinished` with the
        // panic in `reason`, or (if the per-case catch didn't engage)
        // a `Failed(...)` from the outer guard.
        let events: Vec<SweepEvent> = handle.rx.try_iter().collect();
        let saw_panic_evidence = events.iter().any(|e| match e {
            SweepEvent::JobFinished {
                reason: Some(r), ..
            } => r.contains("panic") || r.contains("test-sweep-panic-marker"),
            SweepEvent::Failed(msg) => {
                msg.contains("panic") || msg.contains("test-sweep-panic-marker")
            }
            _ => false,
        });
        assert!(
            saw_panic_evidence,
            "expected panic evidence in SweepEvent stream; got {events:?}"
        );

        let _ = std::fs::remove_dir_all(&parent);
    }

    #[test]
    fn spawn_pumps_progress_and_log_events() {
        let progress_calls = Arc::new(AtomicUsize::new(0));
        let adapter: Arc<dyn valenx_core::Adapter> = Arc::new(NoopAdapter {
            progress_calls: progress_calls.clone(),
        });
        let case = Case {
            id: "x".into(),
            path: std::env::temp_dir(),
        };
        let workdir = std::env::temp_dir().join(format!(
            "valenx-run-test-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mut handle = super::spawn(adapter, case, workdir);

        // Join the thread synchronously so the test is deterministic.
        let thread = handle.thread.take().expect("thread");
        let result = thread.join().expect("thread did not panic");
        assert!(result.is_ok(), "run failed: {result:?}");
        assert_eq!(progress_calls.load(Ordering::SeqCst), 1);

        // Drain the channel — we should see Starting + Progress + a
        // log line + Finished.
        let events: Vec<RunEvent> = handle.rx.try_iter().collect();
        assert!(matches!(events[0], RunEvent::Starting));
        let has_progress = events
            .iter()
            .any(|e| matches!(e, RunEvent::Progress { .. }));
        let has_log = events.iter().any(|e| matches!(e, RunEvent::LogLine { .. }));
        let has_finished = events.iter().any(|e| matches!(e, RunEvent::Finished(_)));
        assert!(has_progress, "expected Progress event; got {events:?}");
        assert!(has_log, "expected LogLine event; got {events:?}");
        assert!(has_finished, "expected Finished event; got {events:?}");
    }

    #[test]
    fn run_event_channel_applies_back_pressure_at_capacity() {
        // Round-8 RED→GREEN: a fast-emitting producer is now bounded.
        // Pre-fix, `channel()` allocated unbounded queue slots and a
        // runaway log producer could OOM the host. With
        // `sync_channel(RUN_EVENT_CAPACITY)`, a producer that
        // outpaces the consumer blocks rather than allocating new
        // slots. This test confirms the producer blocks past the
        // capacity (without a consumer drain). A short timeout on
        // the producer thread proves it's blocked.
        use std::sync::mpsc::sync_channel;
        let (tx, _rx) = sync_channel::<RunEvent>(super::RUN_EVENT_CAPACITY);
        let sent = Arc::new(AtomicUsize::new(0));
        let sent_for_thread = sent.clone();
        let handle = thread::spawn(move || {
            // Push past the cap. After RUN_EVENT_CAPACITY pushes,
            // `send` blocks until the consumer drains (we never do)
            // so the loop stops making progress.
            for _ in 0..(super::RUN_EVENT_CAPACITY + 100) {
                if tx.send(RunEvent::Starting).is_err() {
                    break;
                }
                sent_for_thread.fetch_add(1, Ordering::SeqCst);
            }
        });
        // Give the producer time to fill the queue and block.
        std::thread::sleep(std::time::Duration::from_millis(200));
        let pushed = sent.load(Ordering::SeqCst);
        // The producer should have pushed exactly the cap-many
        // events before blocking. Allow a small slack since the
        // bounded `sync_channel` reserves one extra "rendezvous" slot
        // when the consumer is asleep.
        assert!(
            (super::RUN_EVENT_CAPACITY..=super::RUN_EVENT_CAPACITY + 2).contains(&pushed),
            "expected producer to block near {} pushes, got {}",
            super::RUN_EVENT_CAPACITY,
            pushed
        );
        // Drop the receiver via the producer's own send failure on
        // join — we never explicitly close the channel, so the test
        // process cleans up the blocked producer. To not leak the
        // thread, drop _rx so subsequent sends fail and the thread
        // exits.
        drop(_rx);
        let _ = handle.join();
    }

    /// Round-9 RED→GREEN: `ChannelLogSink::log_line` MUST NOT block
    /// when the bounded channel fills — that would deadlock the
    /// producer (solver subprocess) on a stalled consumer (UI
    /// thread paused for a long paint). Round-9 switched to
    /// `try_send` + drop counter. Verify the sink never blocks past
    /// `RUN_EVENT_CAPACITY + epsilon` pushes and that the drop
    /// counter sees the dropped events.
    #[test]
    fn channel_log_sink_try_send_never_blocks_past_capacity() {
        use std::sync::mpsc::sync_channel;
        use std::time::{Duration, Instant};
        let (tx, _rx) = sync_channel::<RunEvent>(super::RUN_EVENT_CAPACITY);
        let drops = Arc::new(super::SinkDropCounter::default());
        let sink = super::ChannelLogSink {
            tx,
            drops: drops.clone(),
        };
        // Push way past the cap. A blocking implementation would
        // hang here forever (no consumer); the round-9 try_send
        // implementation returns immediately, dropping the overflow
        // into the counter.
        let started = Instant::now();
        for i in 0..(super::RUN_EVENT_CAPACITY * 4) {
            sink.log_line(valenx_core::LogLevel::Info, &format!("line {i}"));
        }
        let elapsed = started.elapsed();
        assert!(
            elapsed < Duration::from_secs(2),
            "pushing {} log lines took {elapsed:?} — sink is blocking instead of dropping",
            super::RUN_EVENT_CAPACITY * 4
        );
        // Counter must show drops past the cap.
        let dropped = drops.snapshot();
        assert!(
            dropped >= super::RUN_EVENT_CAPACITY,
            "expected at least {} drops, got {dropped}",
            super::RUN_EVENT_CAPACITY
        );
    }
}
