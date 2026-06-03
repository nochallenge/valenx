//! Executor abstraction — where an Adapter's `run()` actually
//! happens. First concrete chunk of
//! [RFC 0009](../../../rfcs/0009-hpc-job-submission.md).
//!
//! The trait is intentionally tiny and orthogonal to
//! [`crate::adapter::Adapter`]. Adapters keep their probe /
//! prepare / run / collect lifecycle unchanged. An Executor
//! decides where `run()` actually executes: the existing
//! in-process subprocess flow (`LocalExecutor`) is the only impl
//! that ships today; SLURM / PBS / Kubernetes executors land as
//! separate crates per RFC 0009 §"Migration path."
//!
//! What's intentionally NOT in this commit:
//!
//! - SlurmExecutor / PbsExecutor / K8sExecutor. Per the RFC, they
//!   live in their own crates (`valenx-executor-slurm` etc.) so
//!   pulling them in is opt-in.
//! - Wiring this into the app's run pipeline. Today
//!   `valenx-app::run::spawn` calls `subprocess::run` directly;
//!   switching it to use a configurable Executor is a follow-up.
//! - Remote workdir handling. Local execution doesn't have the
//!   problem; remote executors will need a `RemoteWorkdir` newtype
//!   per RFC 0009's open-questions section.
//! - Streaming logs. `LocalExecutor::submit` redirects stdout/stderr
//!   into files in the workdir; the existing `subprocess::run` flow
//!   in `valenx-app` is what feeds the live UI residual chart, so
//!   we don't duplicate that pump here.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Mutex;

use thiserror::Error;

use crate::adapter::PreparedJob;
use crate::subprocess::KillOnDropChild;

/// Where an Adapter's `run()` actually executes. Each implementation
/// owns one transport (in-process subprocess, SLURM `sbatch`, k8s
/// Job, plain SSH-and-nohup).
///
/// This trait is deliberately small for v0:
/// - no async (the existing subprocess runner is synchronous-on-
///   a-thread and we don't want to bifurcate the world before HPC
///   actually lands)
/// - no streaming logs across the wire (the local executor pipes
///   stdout directly via `RunContext`; remote executors will need
///   per-implementation strategies — see RFC 0009 §"Open questions")
/// - no result fetch-back (gets its own follow-up RFC for file
///   staging)
///
/// All three of those will land additive — add an async variant,
/// add a `tail_logs(handle)` method, etc. — without breaking the
/// initial trait shape.
pub trait Executor {
    /// Stable executor id for UI display + audit logging.
    /// Convention: lowercase, dot-separated (e.g. `local`, `slurm`,
    /// `k8s`).
    fn id(&self) -> &str;

    /// Submit a prepared job. Returns an opaque handle the caller
    /// uses to poll / cancel / fetch results. For the local
    /// executor the handle is just a process id; for SLURM it's
    /// the batch job id; for k8s it's the Job uid.
    fn submit(&self, job: &PreparedJob) -> Result<ExecutorHandle, ExecutorError>;

    /// Non-blocking status check.
    fn poll(&self, handle: &ExecutorHandle) -> Result<RunStatus, ExecutorError>;

    /// Request cancellation. Best-effort across implementations:
    /// local sends SIGTERM (then SIGKILL after a grace window);
    /// SLURM sends `scancel`; k8s deletes the Job. None of these
    /// guarantee the underlying solver process actually exits
    /// promptly — see RFC 0009 §"Drawbacks."
    fn cancel(&self, handle: &ExecutorHandle) -> Result<(), ExecutorError>;
}

/// Opaque handle returned by [`Executor::submit`]. Implementations
/// stash their own state inside the boxed inner value; callers just
/// keep the handle alive until they're done with it.
#[derive(Debug)]
pub struct ExecutorHandle {
    pub executor_id: String,
    pub native_id: String,
    /// Workdir (local FS — remote executors will switch this to a
    /// `RemoteWorkdir` newtype per RFC 0009; today it's just a
    /// PathBuf).
    pub workdir: PathBuf,
}

/// Coarse status of a running / completed job. Maps cleanly across
/// LocalExecutor (process running / exited) and remote schedulers
/// (PENDING / RUNNING / COMPLETED / FAILED / CANCELLED on SLURM,
/// the Pod conditions on k8s).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RunStatus {
    /// Submitted but not yet started (queued by the scheduler).
    /// Local executor never produces this; SLURM / k8s do.
    Pending,
    /// Actively executing.
    Running,
    /// Exited cleanly with the given exit code (typically 0).
    Completed { exit_code: i32 },
    /// Exited non-zero or was killed by the scheduler.
    Failed {
        exit_code: Option<i32>,
        reason: String,
    },
    /// User-requested cancellation honoured.
    Cancelled,
}

/// Errors surfaced by an [`Executor`] implementation.
#[derive(Debug, Error)]
pub enum ExecutorError {
    /// `submit()` couldn't launch / queue the job (bad command, queue
    /// full, missing credentials, etc.).
    #[error("executor `{executor_id}` couldn't submit: {reason}")]
    SubmitFailed {
        /// The offending executor's [`Executor::id`].
        executor_id: String,
        /// Short human-readable explanation.
        reason: String,
    },
    /// `poll()` failed for an outstanding handle (handle dropped, RPC
    /// failure, etc.).
    #[error("executor `{executor_id}` poll failed for handle `{native_id}`: {reason}")]
    PollFailed {
        /// The offending executor's [`Executor::id`].
        executor_id: String,
        /// Executor-native job id.
        native_id: String,
        /// Short human-readable explanation.
        reason: String,
    },
    /// `cancel()` couldn't be honoured (job already finished, RPC
    /// failure, etc.).
    #[error("executor `{executor_id}` cancel failed for handle `{native_id}`: {reason}")]
    CancelFailed {
        /// The offending executor's [`Executor::id`].
        executor_id: String,
        /// Executor-native job id.
        native_id: String,
        /// Short human-readable explanation.
        reason: String,
    },
    /// Fetching results back to the local workdir failed (used by
    /// remote-cluster executors after a successful poll → Completed).
    /// Distinct from PollFailed because the job itself succeeded —
    /// only the post-completion fetch leg broke.
    #[error("executor `{executor_id}` fetch failed for handle `{native_id}`: {reason}")]
    FetchFailed {
        /// The offending executor's [`Executor::id`].
        executor_id: String,
        /// Executor-native job id.
        native_id: String,
        /// Short human-readable explanation.
        reason: String,
    },
    /// The executor variant doesn't implement this leg of the API yet.
    #[error("not implemented: {0}")]
    NotImplemented(&'static str),
}

// ---------------------------------------------------------------------------
// LocalExecutor
// ---------------------------------------------------------------------------

/// In-process subprocess executor.
///
/// `submit()` forks the prepared job's `native_command` as a child
/// process with stdout/stderr redirected to `workdir/stdout.log` and
/// `workdir/stderr.log`. The `Child` handle is stashed in a
/// process-global side-table keyed by an opaque `native_id`; `poll()`
/// and `cancel()` find the entry by that key.
///
/// Adapters that want live log streaming (the existing residual
/// chart in valenx-app) keep using `subprocess::run` directly —
/// LocalExecutor is the path for fire-and-forget runs that produce
/// their own log files (per the Phase 11 HPC story where remote
/// executors land logs in the workdir for later fetch-back).
///
/// Round-14 M9: every stored child is wrapped in
/// [`KillOnDropChild`] (the same RAII guard the `subprocess::run`
/// path uses). Pre-fix `LocalExecutor` stored bare `Child` values
/// in `children` — when the executor itself was dropped (e.g. the
/// UI exited while a sweep was still running) the `HashMap` was
/// torn down with `std::mem::drop`, the bare `Child` handles got
/// dropped without `kill()`, and every outstanding subprocess
/// orphaned (the OS leaves them running for the rest of their
/// allotted runtime). Now the guard's `Drop` issues `kill()` even
/// when the executor exits abnormally, so a dropped executor
/// guarantees no orphaned children.
pub struct LocalExecutor {
    /// next_id stays monotonic across the executor's lifetime so a
    /// dropped handle's id never gets reused.
    next_id: Mutex<u64>,
    /// Process handles indexed by `native_id`. Entries are cleared
    /// on `poll()` once the child exits, so the table doesn't grow
    /// unboundedly. Each child is wrapped in
    /// [`KillOnDropChild`] so executor-drop = clean teardown (no
    /// orphans).
    children: Mutex<HashMap<String, KillOnDropChild>>,
}

impl LocalExecutor {
    /// New, empty executor — starts with no outstanding child
    /// processes and a fresh monotonic id counter.
    pub fn new() -> Self {
        Self {
            next_id: Mutex::new(1),
            children: Mutex::new(HashMap::new()),
        }
    }
}

impl Default for LocalExecutor {
    fn default() -> Self {
        Self::new()
    }
}

impl Executor for LocalExecutor {
    fn id(&self) -> &str {
        "local"
    }

    fn submit(&self, job: &PreparedJob) -> Result<ExecutorHandle, ExecutorError> {
        if job.native_command.is_empty() {
            return Err(ExecutorError::SubmitFailed {
                executor_id: "local".into(),
                reason: "PreparedJob.native_command is empty".into(),
            });
        }
        std::fs::create_dir_all(&job.workdir).map_err(|e| ExecutorError::SubmitFailed {
            executor_id: "local".into(),
            reason: format!("create workdir {}: {e}", job.workdir.display()),
        })?;
        let stdout = std::fs::File::create(job.workdir.join("stdout.log")).map_err(|e| {
            ExecutorError::SubmitFailed {
                executor_id: "local".into(),
                reason: format!("open stdout.log: {e}"),
            }
        })?;
        let stderr = std::fs::File::create(job.workdir.join("stderr.log")).map_err(|e| {
            ExecutorError::SubmitFailed {
                executor_id: "local".into(),
                reason: format!("open stderr.log: {e}"),
            }
        })?;

        let mut cmd = std::process::Command::new(&job.native_command[0]);
        for arg in &job.native_command[1..] {
            cmd.arg(arg);
        }
        for (k, v) in &job.environment {
            cmd.env(k, v);
        }
        cmd.current_dir(&job.workdir);
        cmd.stdout(stdout);
        cmd.stderr(stderr);
        let child = cmd.spawn().map_err(|e| ExecutorError::SubmitFailed {
            executor_id: "local".into(),
            reason: format!("spawn `{}`: {e}", job.native_command[0].to_string_lossy()),
        })?;

        // Allocate native_id and stash the Child wrapped in
        // KillOnDropChild. Round-14 M9: `enabled = true` means
        // dropping the executor (or removing the entry from the
        // map) issues kill(). Pre-fix the bare Child got dropped
        // without kill on executor teardown — orphaned every
        // outstanding subprocess.
        //
        // Round-24 L2: `.expect("mutex poisoned")` panics if any
        // previous holder of the lock panicked. The executor is a
        // long-lived shared service — a single panic in a background
        // poll thread would poison the lock and tear down EVERY
        // subsequent submit/poll/cancel call workspace-wide.
        // Recover gracefully via `into_inner()` — the poisoned data
        // is still valid (we only mutate the next-id counter and
        // a HashMap), so the worst case is a stale counter that
        // self-corrects on the next non-panicked allocation.
        let mut next = self.next_id.lock().unwrap_or_else(|e| e.into_inner());
        let native_id = format!("local-{}", *next);
        *next += 1;
        drop(next);
        self.children
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .insert(native_id.clone(), KillOnDropChild::new(child, true));

        Ok(ExecutorHandle {
            executor_id: "local".to_string(),
            native_id,
            workdir: job.workdir.clone(),
        })
    }

    fn poll(&self, handle: &ExecutorHandle) -> Result<RunStatus, ExecutorError> {
        // Round-24 L2: graceful recovery — see submit() for rationale.
        let mut children = self.children.lock().unwrap_or_else(|e| e.into_inner());
        let Some(child) = children.get_mut(&handle.native_id) else {
            return Err(ExecutorError::PollFailed {
                executor_id: "local".into(),
                native_id: handle.native_id.clone(),
                reason:
                    "handle not found in process table — already polled-to-completion or never \
                         submitted by this executor"
                        .into(),
            });
        };
        match child.inner_mut().try_wait() {
            Ok(None) => Ok(RunStatus::Running),
            Ok(Some(status)) => {
                // Process exited — clean up the table entry. The
                // KillOnDropChild's Drop sees `try_wait` returns
                // `Some` (because we just polled it) and skips the
                // kill, so this is a clean teardown.
                children.remove(&handle.native_id);
                let exit_code = status.code().unwrap_or(-1);
                if exit_code == 0 {
                    Ok(RunStatus::Completed { exit_code })
                } else {
                    Ok(RunStatus::Failed {
                        exit_code: status.code(),
                        reason: format!("exit code {exit_code}"),
                    })
                }
            }
            Err(e) => Err(ExecutorError::PollFailed {
                executor_id: "local".into(),
                native_id: handle.native_id.clone(),
                reason: format!("try_wait: {e}"),
            }),
        }
    }

    fn cancel(&self, handle: &ExecutorHandle) -> Result<(), ExecutorError> {
        // Round-24 L2: graceful recovery — see submit() for rationale.
        let mut children = self.children.lock().unwrap_or_else(|e| e.into_inner());
        let Some(child) = children.get_mut(&handle.native_id) else {
            return Err(ExecutorError::CancelFailed {
                executor_id: "local".into(),
                native_id: handle.native_id.clone(),
                reason: "handle not found in process table".into(),
            });
        };
        // Best-effort: send SIGKILL on Unix / TerminateProcess on
        // Windows. Real graceful shutdown (SIGTERM + grace window
        // + SIGKILL) is a follow-up — the existing subprocess::run
        // path uses CancellationToken polling for that and we can
        // factor it out later.
        child
            .inner_mut()
            .kill()
            .map_err(|e| ExecutorError::CancelFailed {
                executor_id: "local".into(),
                native_id: handle.native_id.clone(),
                reason: format!("kill: {e}"),
            })?;
        // Drop the entry — `poll()` after this returns
        // PollFailed("handle not found"). Callers should treat
        // cancel as terminal. The KillOnDropChild's Drop re-attempts
        // a kill (no-op if try_wait shows exited).
        children.remove(&handle.native_id);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_job() -> PreparedJob {
        PreparedJob {
            workdir: PathBuf::from("/tmp/x"),
            native_command: vec![],
            environment: vec![],
            estimated_runtime: None,
            kill_on_drop: false,
        }
    }

    #[test]
    fn local_executor_id_is_local() {
        assert_eq!(LocalExecutor::new().id(), "local");
    }

    #[test]
    fn local_executor_rejects_empty_native_command() {
        let exec = LocalExecutor::new();
        let err = exec.submit(&dummy_job()).unwrap_err();
        match err {
            ExecutorError::SubmitFailed { reason, .. } => {
                assert!(reason.contains("native_command is empty"));
            }
            other => panic!("expected SubmitFailed, got {other:?}"),
        }
    }

    /// Spawn a real short-lived subprocess and verify the full
    /// submit → poll-to-completion lifecycle. Uses cmd.exe on
    /// Windows and /bin/sh on Unix; both are universally
    /// available so the test isn't fragile.
    #[test]
    fn local_executor_runs_a_real_command() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-local-exec-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let job = if cfg!(windows) {
            PreparedJob {
                workdir: tmp.clone(),
                native_command: vec![
                    std::ffi::OsString::from("cmd.exe"),
                    std::ffi::OsString::from("/C"),
                    std::ffi::OsString::from("echo hi"),
                ],
                environment: vec![],
                estimated_runtime: None,
                kill_on_drop: false,
            }
        } else {
            PreparedJob {
                workdir: tmp.clone(),
                native_command: vec![
                    std::ffi::OsString::from("/bin/sh"),
                    std::ffi::OsString::from("-c"),
                    std::ffi::OsString::from("echo hi"),
                ],
                environment: vec![],
                estimated_runtime: None,
                kill_on_drop: false,
            }
        };
        let exec = LocalExecutor::new();
        let handle = exec.submit(&job).expect("submit");
        assert!(handle.native_id.starts_with("local-"));

        // Poll until done (up to 5 seconds — short-lived command).
        let mut status = exec.poll(&handle).expect("poll");
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(5);
        while matches!(status, RunStatus::Running) && std::time::Instant::now() < deadline {
            std::thread::sleep(std::time::Duration::from_millis(50));
            status = exec.poll(&handle).expect("poll");
        }
        match status {
            RunStatus::Completed { exit_code } => {
                assert_eq!(exit_code, 0, "echo should exit 0");
            }
            other => panic!("expected Completed, got {other:?}"),
        }

        // stdout.log should exist and contain "hi".
        let stdout = std::fs::read_to_string(tmp.join("stdout.log")).expect("stdout");
        assert!(stdout.contains("hi"), "stdout was: {stdout:?}");

        // Polling again returns "handle not found" — the entry was
        // cleaned up on completion.
        let err = exec.poll(&handle).unwrap_err();
        assert!(matches!(err, ExecutorError::PollFailed { .. }));

        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn local_executor_handle_native_id_is_unique_per_submit() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-local-exec-uniq-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let exec = LocalExecutor::new();
        let mut native_ids = Vec::new();
        for _ in 0..3 {
            let job = if cfg!(windows) {
                PreparedJob {
                    workdir: tmp.clone(),
                    native_command: vec![
                        std::ffi::OsString::from("cmd.exe"),
                        std::ffi::OsString::from("/C"),
                        std::ffi::OsString::from("exit 0"),
                    ],
                    environment: vec![],
                    estimated_runtime: None,
                    kill_on_drop: false,
                }
            } else {
                PreparedJob {
                    workdir: tmp.clone(),
                    native_command: vec![
                        std::ffi::OsString::from("/bin/sh"),
                        std::ffi::OsString::from("-c"),
                        std::ffi::OsString::from("true"),
                    ],
                    environment: vec![],
                    estimated_runtime: None,
                    kill_on_drop: false,
                }
            };
            let handle = exec.submit(&job).expect("submit");
            native_ids.push(handle.native_id);
        }
        // All three ids are distinct + monotonic.
        let mut sorted = native_ids.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(sorted.len(), 3, "ids should be unique: {native_ids:?}");
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn run_status_variants_compare_correctly() {
        assert_eq!(
            RunStatus::Completed { exit_code: 0 },
            RunStatus::Completed { exit_code: 0 }
        );
        assert_ne!(
            RunStatus::Completed { exit_code: 0 },
            RunStatus::Completed { exit_code: 1 }
        );
        assert_ne!(RunStatus::Pending, RunStatus::Running);
        assert_ne!(
            RunStatus::Failed {
                exit_code: Some(1),
                reason: "x".into(),
            },
            RunStatus::Cancelled
        );
    }

    #[test]
    fn executor_error_messages_include_executor_id() {
        let err = ExecutorError::SubmitFailed {
            executor_id: "slurm".into(),
            reason: "sbatch not on PATH".into(),
        };
        let msg = err.to_string();
        assert!(msg.contains("slurm"));
        assert!(msg.contains("sbatch"));
    }

    /// Round-14 M9 RED→GREEN: dropping `LocalExecutor` must kill
    /// every outstanding child. Pre-fix the executor stored bare
    /// `Child` handles, so a drop just torn down the HashMap without
    /// signalling anything — the OS kept the children alive for the
    /// full duration of whatever they were doing.
    ///
    /// We verify by spawning a child that would normally write a
    /// "tombstone" marker file after a 5-second delay. The executor
    /// is dropped well before the 5-second mark; the test then
    /// waits 2 seconds and checks the marker is NOT present.
    /// Without the kill_on_drop guard, the marker file would exist
    /// (the child outlived the executor); with the guard the child
    /// is killed and the marker never gets written.
    #[test]
    fn local_executor_kills_outstanding_children_on_drop() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-local-exec-killondrop-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let marker = tmp.join("tombstone.txt");

        // Cross-platform sleep + write — on Windows we use cmd.exe
        // (`timeout` exists but is interactive-aware and unreliable
        // in non-tty contexts; PowerShell sleep + Out-File works).
        // On Unix we use /bin/sh.
        let job = if cfg!(windows) {
            // PowerShell: Start-Sleep -Seconds 5; New-Item -Path 'marker' -ItemType File
            PreparedJob {
                workdir: tmp.clone(),
                native_command: vec![
                    std::ffi::OsString::from("powershell.exe"),
                    std::ffi::OsString::from("-NoProfile"),
                    std::ffi::OsString::from("-Command"),
                    std::ffi::OsString::from(format!(
                        "Start-Sleep -Seconds 5; New-Item -Path '{}' -ItemType File -Force | Out-Null",
                        marker.display()
                    )),
                ],
                environment: vec![],
                estimated_runtime: None,
                kill_on_drop: true,
            }
        } else {
            PreparedJob {
                workdir: tmp.clone(),
                native_command: vec![
                    std::ffi::OsString::from("/bin/sh"),
                    std::ffi::OsString::from("-c"),
                    std::ffi::OsString::from(format!(
                        "sleep 5 && touch '{}'",
                        marker.display()
                    )),
                ],
                environment: vec![],
                estimated_runtime: None,
                kill_on_drop: true,
            }
        };

        {
            // Scope the executor so it drops before we wait. The
            // submit returns immediately (the spawn is non-blocking);
            // we do NOT poll — the bug we're chasing is exactly "the
            // caller submitted but never polled, then walked away".
            let exec = LocalExecutor::new();
            let _handle = exec.submit(&job).expect("submit");
            // exec drops here → KillOnDropChild::drop fires → child killed.
        }

        // Wait 2 seconds — well short of the 5-second sleep the
        // child would otherwise execute. If the kill fired, the
        // marker never appears.
        std::thread::sleep(std::time::Duration::from_secs(2));
        assert!(
            !marker.exists(),
            "marker file exists at {} — kill_on_drop did NOT fire, child outlived executor",
            marker.display()
        );

        // Cleanup. Note: if for some reason the kill didn't take and
        // the child finishes after our 2-second wait, this cleanup
        // races against the eventual marker write. We give a
        // generous extra few seconds to drain.
        std::thread::sleep(std::time::Duration::from_secs(5));
        let _ = std::fs::remove_dir_all(&tmp);
    }

    /// RED→GREEN (round-24 L2): the executor's `submit`/`poll`/`cancel`
    /// methods must survive a poisoned mutex. Pre-fix any `.expect(
    /// "mutex poisoned")` call would have abort-panicked the calling
    /// thread the first time a previous holder of the lock panicked
    /// — a single bad future would tear down EVERY subsequent
    /// executor call. Post-fix `.unwrap_or_else(|e| e.into_inner())`
    /// recovers the inner guard and lets the call complete.
    ///
    /// We poison the lock from a controlled thread that panics
    /// while holding the lock, then exercise `submit` and `poll`
    /// from the test thread — pre-fix both would panic the test;
    /// post-fix both complete normally.
    #[test]
    fn local_executor_survives_poisoned_mutex() {
        let exec = std::sync::Arc::new(LocalExecutor::new());
        // Poison the `children` lock by panicking while holding it.
        let exec_clone = std::sync::Arc::clone(&exec);
        let poison_thread = std::thread::spawn(move || {
            let _guard = exec_clone
                .children
                .lock()
                .expect("first lock is unpoisoned");
            // Panic with the lock held — the mutex is now poisoned.
            panic!("intentional poison");
        });
        let join_result = poison_thread.join();
        assert!(join_result.is_err(), "poison thread must have panicked");

        // Post-fix: submit() / poll() / cancel() must not panic.
        // Use a no-op job (empty native_command) — `submit` errors
        // synchronously with `SubmitFailed` BEFORE touching either
        // lock, so it can't be used to probe the poisoned lock.
        // Use a working short command instead.
        let tmp = std::env::temp_dir().join(format!(
            "valenx-l2-poison-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&tmp).unwrap();
        let job = if cfg!(windows) {
            PreparedJob {
                workdir: tmp.clone(),
                native_command: vec![
                    std::ffi::OsString::from("cmd.exe"),
                    std::ffi::OsString::from("/C"),
                    std::ffi::OsString::from("exit 0"),
                ],
                environment: vec![],
                estimated_runtime: None,
                kill_on_drop: false,
            }
        } else {
            PreparedJob {
                workdir: tmp.clone(),
                native_command: vec![
                    std::ffi::OsString::from("/bin/sh"),
                    std::ffi::OsString::from("-c"),
                    std::ffi::OsString::from("true"),
                ],
                environment: vec![],
                estimated_runtime: None,
                kill_on_drop: false,
            }
        };
        // submit() touches both locks; if either panics the test
        // thread panics here. Post-fix both .lock() calls use
        // unwrap_or_else(into_inner) so the call completes.
        let handle = exec.submit(&job).expect("submit must succeed despite poison");
        // poll() also touches the children lock — same check.
        let _ = exec.poll(&handle); // status doesn't matter; we're checking it doesn't panic
        let _ = std::fs::remove_dir_all(&tmp);
    }
}
