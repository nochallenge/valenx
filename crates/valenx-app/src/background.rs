//! A one-shot background computation, polled from the UI thread.
//!
//! egui is immediate-mode and runs on a single thread: a heavy solve called
//! straight from a button handler blocks the event loop, so the whole window
//! freezes until it returns. The reactive pattern is to move the solve onto a
//! worker thread and poll its result from the UI thread each frame (requesting
//! a repaint so frames keep ticking while it runs).
//!
//! Several workbenches grew their own copies of that thread + channel
//! boilerplate ([`crate::reactdyn_workbench`], [`crate::aero`], the RNA
//! designer). [`BackgroundJob`] captures the same pattern generically so a
//! workbench needs only a few lines:
//!
//! ```ignore
//! // on the "▶ Solve" click (the disabled-while-running form means one
//! // click == one job):
//! let inputs = snapshot(s);                          // owned, `Send`
//! s.job = Some(BackgroundJob::spawn(move || solve(inputs)));
//!
//! // near the top of the panel's draw, every frame:
//! match s.job.as_mut().map(BackgroundJob::poll) {
//!     Some(JobState::Done(out)) => { s.job = None; apply(s, out); }
//!     Some(JobState::Failed)    => { s.job = None; s.error = Some("…".into()); }
//!     Some(JobState::Pending) | None => {}
//! }
//! if s.job.is_some() {
//!     ui.ctx().request_repaint();                     // keep frames ticking
//! }
//! ```

use std::sync::mpsc::{channel, Receiver, TryRecvError};
use std::thread::{self, JoinHandle};

/// The state of a [`BackgroundJob`] at a poll.
pub enum JobState<T> {
    /// The worker is still running; poll again next frame.
    Pending,
    /// The worker finished; here is its result.
    Done(T),
    /// The worker thread vanished (panicked) without producing a result.
    /// In practice workers run over validated input and don't panic, but a
    /// panel must still handle this so a dead worker can't wedge the UI in a
    /// permanent "running" state.
    Failed,
}

/// A computation running on a worker thread, producing one value of type `T`.
///
/// Poll it each frame with [`BackgroundJob::poll`]: it returns
/// [`JobState::Pending`] while the worker runs and, exactly once, either
/// [`JobState::Done`] with the result or [`JobState::Failed`] if the worker
/// panicked. Drop the job (e.g. set its `Option` to `None`) once you observe a
/// non-`Pending` state; don't keep polling a consumed job.
pub struct BackgroundJob<T: Send + 'static> {
    rx: Receiver<T>,
    handle: Option<JoinHandle<()>>,
}

impl<T: Send + 'static> BackgroundJob<T> {
    /// Spawn `f` on a worker thread. The closure must own everything it needs
    /// (`Send + 'static`) — snapshot any UI state into owned values before
    /// calling this; the worker must not touch `egui`/app state.
    pub fn spawn(f: impl FnOnce() -> T + Send + 'static) -> Self {
        let (tx, rx) = channel();
        let handle = thread::spawn(move || {
            // The receiver is dropped only when the job is dropped; a send
            // error then just means the UI stopped caring, so ignore it.
            let _ = tx.send(f());
        });
        Self {
            rx,
            handle: Some(handle),
        }
    }

    /// Non-blocking. See [`JobState`]. The first non-`Pending` result reaps the
    /// worker thread; don't poll again after that.
    pub fn poll(&mut self) -> JobState<T> {
        match self.rx.try_recv() {
            Ok(value) => {
                // Reap the finished thread so it isn't left as a zombie handle.
                if let Some(h) = self.handle.take() {
                    let _ = h.join();
                }
                JobState::Done(value)
            }
            Err(TryRecvError::Empty) => JobState::Pending,
            Err(TryRecvError::Disconnected) => {
                // Worker dropped its sender without sending — it panicked.
                self.handle = None;
                JobState::Failed
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    /// Bounded poll loop so a test can never hang on a stuck worker.
    fn drive<T: Send + 'static>(job: &mut BackgroundJob<T>) -> JobState<T> {
        for _ in 0..2000 {
            match job.poll() {
                JobState::Pending => thread::sleep(Duration::from_millis(1)),
                done_or_failed => return done_or_failed,
            }
        }
        JobState::Pending
    }

    #[test]
    fn delivers_the_workers_result() {
        let mut job = BackgroundJob::spawn(|| 2 + 2);
        match drive(&mut job) {
            JobState::Done(v) => assert_eq!(v, 4, "worker result must reach the poller"),
            JobState::Pending => panic!("worker never finished within the poll budget"),
            JobState::Failed => panic!("worker should not have failed"),
        }
    }

    #[test]
    fn carries_a_large_owned_value() {
        // Mirror real worker outputs (owned Vecs/Strings sent across the channel).
        let mut job = BackgroundJob::spawn(|| vec![7u8; 4096]);
        match drive(&mut job) {
            JobState::Done(v) => assert!(v.len() == 4096 && v.iter().all(|&b| b == 7)),
            other => panic!("expected Done, got {}", state_name(&other)),
        }
    }

    #[test]
    fn a_panicked_worker_reports_failed_not_pending() {
        // A worker that panics drops its sender without sending; the job must
        // surface `Failed` rather than spin in `Pending` forever.
        let mut job = BackgroundJob::<i32>::spawn(|| panic!("worker blew up"));
        match drive(&mut job) {
            JobState::Failed => {}
            JobState::Done(_) => panic!("a panicked worker must not yield a value"),
            JobState::Pending => panic!("a panicked worker must not spin forever"),
        }
    }

    fn state_name<T>(s: &JobState<T>) -> &'static str {
        match s {
            JobState::Pending => "Pending",
            JobState::Done(_) => "Done",
            JobState::Failed => "Failed",
        }
    }
}
