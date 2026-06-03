//! Concurrent participant orchestration for preCICE meta-cases.
//!
//! Phase 9 tail — closes the loop between `PreciceAdapter::prepare`
//! (which stages the participants + writes a manifest) and the
//! actual run (which needs each participant's solver running
//! simultaneously against the shared preCICE coupling library).
//!
//! ## Orchestration model
//!
//! preCICE participants communicate at run time via either MPI ports
//! or sockets that the preCICE library opens inside each
//! participant's process. From Valenx's perspective each participant
//! is just a subprocess that the user's installed solvers (compiled
//! with `precice` linkage) handle internally — we don't poke into
//! the participant's coupling state at all.
//!
//! What we DO own:
//!
//! 1. **Concurrent submission**: every participant's
//!    [`PreparedJob`] gets submitted via the same
//!    [`Executor`] in quick succession so they're all alive when the
//!    first preCICE handshake happens. Sequential submission would
//!    deadlock — the first participant blocks on the second's
//!    handshake before the second is even running.
//! 2. **Joint polling**: a sweep over all handles each tick;
//!    aggregate status is "Running" while any participant is still
//!    running, "Completed" only when every one finished cleanly,
//!    "Failed" the moment any one fails or the user cancels.
//! 3. **Joint cancellation**: cancelling one participant strands
//!    the others on a coupling read; we propagate the cancel to
//!    every handle so the whole job tears down at once.
//!
//! Per-participant `Results` collection happens after every handle
//! reaches a terminal status; we hand each participant's adapter the
//! relevant subdirectory and aggregate the resulting Results into a
//! shared catalog with a participant-name prefix on every field.
//!
//! ## Why a separate module
//!
//! The orchestrator is independently testable: it only depends on
//! the `valenx_core::Executor` trait, [`PreparedJob`], and an opaque
//! per-participant collect closure. Tests substitute
//! `valenx_core::LocalExecutor` with a no-op subprocess (`echo`) so we
//! can exercise the joint polling / cancellation paths without
//! needing a real preCICE installation.

use std::path::PathBuf;
use std::time::{Duration, Instant};

use thiserror::Error;
use valenx_core::{Executor, ExecutorError, ExecutorHandle, PreparedJob, RunStatus};

/// One participant and its prepared job. The orchestrator submits
/// each one in declaration order, so put the long-running solver
/// first if you care about who establishes the coupling listener.
pub struct ParticipantJob {
    pub name: String,
    pub adapter_id: String,
    pub prepared: PreparedJob,
}

/// Active orchestration session. Returned by [`submit_all`]; the
/// caller polls / cancels via the methods on this struct.
pub struct OrchestratorHandle<'e, E: Executor + ?Sized> {
    /// The executor every participant submits through. Borrowed for
    /// the orchestration session's lifetime so handles + executor
    /// stay paired.
    executor: &'e E,
    /// Per-participant submission record.
    pub participants: Vec<RunningParticipant>,
}

/// One participant's executor handle plus the metadata the caller
/// passed in. Field-public so tests + the run loop can introspect.
pub struct RunningParticipant {
    pub name: String,
    pub adapter_id: String,
    pub handle: ExecutorHandle,
    /// Most recent terminal status, populated when a poll observes
    /// the participant in Completed / Failed / Cancelled. `None`
    /// while still Pending / Running.
    pub final_status: Option<RunStatus>,
}

/// Aggregate status across every participant. Mirrors [`RunStatus`]
/// but with the wider semantics needed for a multi-participant
/// session (`PartialFailure` doesn't exist for a single job).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum OrchestratorStatus {
    /// At least one participant is still running.
    InProgress {
        running: usize,
        completed: usize,
        failed: usize,
    },
    /// Every participant completed cleanly (exit 0).
    AllSucceeded { participants: usize },
    /// Every participant reached a terminal status; at least one
    /// failed (non-zero exit) or was cancelled.
    PartialFailure {
        succeeded: usize,
        failed: usize,
        first_failure: Option<String>,
    },
}

#[derive(Debug, Error)]
pub enum OrchestratorError {
    #[error("submit `{participant}`: {source}")]
    Submit {
        participant: String,
        #[source]
        source: ExecutorError,
    },
    #[error("poll `{participant}`: {source}")]
    Poll {
        participant: String,
        #[source]
        source: ExecutorError,
    },
}

/// Submit every participant via the shared executor in declaration
/// order. The returned [`OrchestratorHandle`] keeps the per-handle
/// state for subsequent polling.
///
/// On submit failure the partial state is rolled back: every
/// participant submitted before the failure gets a best-effort
/// `executor.cancel()` so the user doesn't end up with orphaned
/// half-coupled processes. The original submit error is returned
/// after the rollback completes.
pub fn submit_all<'e, E: Executor + ?Sized>(
    executor: &'e E,
    jobs: Vec<ParticipantJob>,
) -> Result<OrchestratorHandle<'e, E>, OrchestratorError> {
    let mut participants: Vec<RunningParticipant> = Vec::with_capacity(jobs.len());
    for job in jobs {
        match executor.submit(&job.prepared) {
            Ok(handle) => participants.push(RunningParticipant {
                name: job.name,
                adapter_id: job.adapter_id,
                handle,
                final_status: None,
            }),
            Err(source) => {
                // Best-effort rollback: cancel everything submitted so
                // far so partial coupling sessions don't leak processes.
                for already in &participants {
                    let _ = executor.cancel(&already.handle);
                }
                return Err(OrchestratorError::Submit {
                    participant: job.name,
                    source,
                });
            }
        }
    }
    Ok(OrchestratorHandle {
        executor,
        participants,
    })
}

impl<'e, E: Executor + ?Sized> OrchestratorHandle<'e, E> {
    /// Poll every participant once. Updates `final_status` for any
    /// participant that reached a terminal state since the last poll.
    /// Returns the aggregate orchestrator status.
    pub fn poll(&mut self) -> Result<OrchestratorStatus, OrchestratorError> {
        let mut running = 0usize;
        let mut completed = 0usize;
        let mut failed = 0usize;
        let mut first_failure: Option<String> = None;
        for p in &mut self.participants {
            // Skip already-terminal participants — a finished job
            // can't transition back, and re-polling could surface a
            // PollFailed if the executor's table evicted the entry.
            if let Some(prev) = &p.final_status {
                match prev {
                    RunStatus::Completed { exit_code: 0 } => completed += 1,
                    RunStatus::Completed { .. }
                    | RunStatus::Failed { .. }
                    | RunStatus::Cancelled => {
                        failed += 1;
                        if first_failure.is_none() {
                            first_failure = Some(p.name.clone());
                        }
                    }
                    RunStatus::Pending | RunStatus::Running => running += 1,
                }
                continue;
            }
            match self.executor.poll(&p.handle) {
                Ok(status) => match &status {
                    RunStatus::Pending | RunStatus::Running => {
                        running += 1;
                    }
                    RunStatus::Completed { exit_code: 0 } => {
                        p.final_status = Some(status);
                        completed += 1;
                    }
                    RunStatus::Completed { .. }
                    | RunStatus::Failed { .. }
                    | RunStatus::Cancelled => {
                        p.final_status = Some(status);
                        failed += 1;
                        if first_failure.is_none() {
                            first_failure = Some(p.name.clone());
                        }
                    }
                },
                Err(source) => {
                    return Err(OrchestratorError::Poll {
                        participant: p.name.clone(),
                        source,
                    });
                }
            }
        }
        Ok(if running > 0 {
            OrchestratorStatus::InProgress {
                running,
                completed,
                failed,
            }
        } else if failed == 0 {
            OrchestratorStatus::AllSucceeded {
                participants: completed,
            }
        } else {
            OrchestratorStatus::PartialFailure {
                succeeded: completed,
                failed,
                first_failure,
            }
        })
    }

    /// Best-effort cancel every still-running participant. Returns
    /// the number of participants the executor accepted the cancel
    /// for. Does NOT wait for the cancellation to propagate; call
    /// [`Self::wait_until_terminal`] if you need to block.
    pub fn cancel_all(&self) -> usize {
        let mut cancelled = 0;
        for p in &self.participants {
            if p.final_status.is_some() {
                continue;
            }
            if self.executor.cancel(&p.handle).is_ok() {
                cancelled += 1;
            }
        }
        cancelled
    }

    /// Block until every participant reaches a terminal status or
    /// the timeout elapses. Polls every `tick`. Returns the final
    /// aggregate status (which is one of AllSucceeded /
    /// PartialFailure when within the timeout, or InProgress if the
    /// timeout fired with participants still running).
    pub fn wait_until_terminal(
        &mut self,
        tick: Duration,
        timeout: Option<Duration>,
    ) -> Result<OrchestratorStatus, OrchestratorError> {
        let started = Instant::now();
        loop {
            let status = self.poll()?;
            if !matches!(status, OrchestratorStatus::InProgress { .. }) {
                return Ok(status);
            }
            if let Some(t) = timeout {
                if started.elapsed() >= t {
                    return Ok(status);
                }
            }
            std::thread::sleep(tick);
        }
    }

    /// Per-participant workdirs, in declaration order. Useful for
    /// the post-orchestration `collect()` pass that walks each
    /// participant's adapter to gather Results.
    pub fn participant_workdirs(&self) -> Vec<(String, PathBuf)> {
        self.participants
            .iter()
            .map(|p| (p.name.clone(), p.handle.workdir.clone()))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    /// Synthetic Executor that pretends each submitted job
    /// transitions through Pending -> Running -> Completed at a
    /// configurable rate. Avoids the need for real subprocesses in
    /// orchestration tests.
    struct FakeExecutor {
        next_id: Mutex<u64>,
        /// Polls per handle before transitioning to Completed. Each
        /// handle's counter starts at this value and decrements
        /// every poll.
        polls_until_done: u64,
        /// When false, every submit succeeds; when true, the second
        /// submit fails with NotImplemented so we can exercise the
        /// rollback path.
        fail_second_submit: bool,
        /// After a participant finishes, what RunStatus to return.
        final_status: RunStatus,
        /// Per-handle remaining poll count.
        state: Mutex<std::collections::HashMap<String, u64>>,
        /// Per-handle latest status (for cancellation tracking).
        cancelled: Mutex<std::collections::HashSet<String>>,
    }

    impl FakeExecutor {
        fn new(polls_until_done: u64) -> Self {
            Self {
                next_id: Mutex::new(1),
                polls_until_done,
                fail_second_submit: false,
                final_status: RunStatus::Completed { exit_code: 0 },
                state: Mutex::new(std::collections::HashMap::new()),
                cancelled: Mutex::new(std::collections::HashSet::new()),
            }
        }

        fn with_failure_status(mut self, st: RunStatus) -> Self {
            self.final_status = st;
            self
        }

        fn failing_second_submit() -> Self {
            let mut me = Self::new(1);
            me.fail_second_submit = true;
            me
        }
    }

    impl Executor for FakeExecutor {
        fn id(&self) -> &str {
            "fake"
        }
        fn submit(&self, job: &PreparedJob) -> Result<ExecutorHandle, ExecutorError> {
            let mut id = self.next_id.lock().unwrap();
            let n = *id;
            *id += 1;
            if self.fail_second_submit && n == 2 {
                return Err(ExecutorError::SubmitFailed {
                    executor_id: "fake".into(),
                    reason: "synthetic submit failure".into(),
                });
            }
            let native_id = format!("fake-{n}");
            self.state
                .lock()
                .unwrap()
                .insert(native_id.clone(), self.polls_until_done);
            Ok(ExecutorHandle {
                executor_id: "fake".into(),
                native_id,
                workdir: job.workdir.clone(),
            })
        }
        fn poll(&self, handle: &ExecutorHandle) -> Result<RunStatus, ExecutorError> {
            if self.cancelled.lock().unwrap().contains(&handle.native_id) {
                return Ok(RunStatus::Cancelled);
            }
            let mut state = self.state.lock().unwrap();
            let counter = state.entry(handle.native_id.clone()).or_insert(0);
            if *counter == 0 {
                Ok(self.final_status.clone())
            } else {
                *counter -= 1;
                Ok(RunStatus::Running)
            }
        }
        fn cancel(&self, handle: &ExecutorHandle) -> Result<(), ExecutorError> {
            self.cancelled
                .lock()
                .unwrap()
                .insert(handle.native_id.clone());
            Ok(())
        }
    }

    fn dummy_job(name: &str) -> ParticipantJob {
        ParticipantJob {
            name: name.into(),
            adapter_id: "test".into(),
            prepared: PreparedJob {
                workdir: std::env::temp_dir().join(format!("orch-{name}")),
                native_command: vec!["true".into()],
                environment: Vec::new(),
                estimated_runtime: None,
                kill_on_drop: false,
            },
        }
    }

    #[test]
    fn submit_all_assigns_a_handle_to_every_participant() {
        let exec = FakeExecutor::new(0);
        let handle =
            submit_all(&exec, vec![dummy_job("Fluid"), dummy_job("Solid")]).expect("submit");
        assert_eq!(handle.participants.len(), 2);
        assert_eq!(handle.participants[0].name, "Fluid");
        assert_eq!(handle.participants[1].name, "Solid");
    }

    #[test]
    fn submit_failure_rolls_back_already_submitted_handles() {
        let exec = FakeExecutor::failing_second_submit();
        let result = submit_all(&exec, vec![dummy_job("Fluid"), dummy_job("Solid")]);
        match result {
            Err(OrchestratorError::Submit { participant, .. }) => {
                assert_eq!(participant, "Solid");
            }
            Err(other) => panic!("wrong error: {other:?}"),
            Ok(_) => panic!("expected Err, got Ok"),
        }
        // The first participant's submit succeeded; rollback must
        // have requested cancel on it so its handle doesn't leak.
        assert!(exec.cancelled.lock().unwrap().contains("fake-1"));
    }

    #[test]
    fn poll_reports_in_progress_until_every_participant_finishes() {
        // Each fake job needs 3 polls before transitioning to
        // Completed. Two participants => first two polls have
        // running=2; third sees both flip.
        let exec = FakeExecutor::new(3);
        let mut handle = submit_all(&exec, vec![dummy_job("Fluid"), dummy_job("Solid")]).unwrap();
        for _ in 0..3 {
            match handle.poll().unwrap() {
                OrchestratorStatus::InProgress { running, .. } => {
                    assert_eq!(running, 2);
                }
                other => panic!("expected InProgress, got {other:?}"),
            }
        }
        // After enough polls, both finish.
        let final_status = handle.poll().unwrap();
        match final_status {
            OrchestratorStatus::AllSucceeded { participants } => assert_eq!(participants, 2),
            other => panic!("expected AllSucceeded, got {other:?}"),
        }
    }

    #[test]
    fn poll_surfaces_partial_failure_when_any_participant_exits_nonzero() {
        let exec = FakeExecutor::new(0).with_failure_status(RunStatus::Failed {
            exit_code: Some(1),
            reason: "synthetic crash".into(),
        });
        let mut handle = submit_all(&exec, vec![dummy_job("Fluid"), dummy_job("Solid")]).unwrap();
        // First poll has both transition to Failed.
        let status = handle.poll().unwrap();
        match status {
            OrchestratorStatus::PartialFailure {
                succeeded,
                failed,
                first_failure,
            } => {
                assert_eq!(succeeded, 0);
                assert_eq!(failed, 2);
                assert!(first_failure.is_some());
            }
            other => panic!("expected PartialFailure, got {other:?}"),
        }
    }

    #[test]
    fn cancel_all_propagates_to_every_running_participant() {
        let exec = FakeExecutor::new(100); // very long-running
        let handle = submit_all(
            &exec,
            vec![dummy_job("Fluid"), dummy_job("Solid"), dummy_job("Mesh")],
        )
        .unwrap();
        let cancelled = handle.cancel_all();
        assert_eq!(cancelled, 3);
        assert_eq!(exec.cancelled.lock().unwrap().len(), 3);
    }

    #[test]
    fn cancel_all_skips_already_terminal_participants() {
        let exec = FakeExecutor::new(0);
        let mut handle = submit_all(&exec, vec![dummy_job("Fluid"), dummy_job("Solid")]).unwrap();
        // First poll terminates both.
        let _ = handle.poll().unwrap();
        // Cancel after terminal: nothing to cancel.
        let cancelled = handle.cancel_all();
        assert_eq!(cancelled, 0);
    }

    #[test]
    fn wait_until_terminal_returns_when_every_participant_finishes() {
        let exec = FakeExecutor::new(2);
        let mut handle = submit_all(&exec, vec![dummy_job("Fluid"), dummy_job("Solid")]).unwrap();
        let status = handle
            .wait_until_terminal(Duration::from_millis(1), None)
            .unwrap();
        match status {
            OrchestratorStatus::AllSucceeded { participants } => {
                assert_eq!(participants, 2);
            }
            other => panic!("expected AllSucceeded, got {other:?}"),
        }
    }

    #[test]
    fn participant_workdirs_returns_per_participant_paths() {
        let exec = FakeExecutor::new(0);
        let handle = submit_all(&exec, vec![dummy_job("Fluid"), dummy_job("Solid")]).unwrap();
        let dirs = handle.participant_workdirs();
        assert_eq!(dirs.len(), 2);
        assert_eq!(dirs[0].0, "Fluid");
        assert!(dirs[0].1.to_string_lossy().contains("orch-Fluid"));
    }
}
