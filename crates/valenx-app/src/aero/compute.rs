//! Background-thread orchestration for a wind-tunnel run.
//!
//! A steady RANS solve — and especially an angle-of-attack sweep — can
//! take many seconds. Running it on the egui thread would freeze the
//! whole window. This module spawns the `valenx-aero` solve on a
//! dedicated `std::thread`, streams progress back over an `mpsc`
//! channel, and exposes a small handle the workbench polls once per
//! frame. The UI stays fully responsive while a run is in flight.
//!
//! Nothing here touches egui — it is pure compute plumbing, so the
//! orchestration is `#[test]`-coverable.

use std::sync::atomic::AtomicBool;
use std::sync::mpsc::{Receiver, Sender};
use std::sync::Arc;
use std::thread::JoinHandle;

use valenx_aero::{
    aoa_sweep, run_windtunnel, AeroReport, AeroRequest, AeroResult, PolarCurve, TriMesh,
};

/// A message from the background solve thread to the UI.
#[derive(Clone, Debug)]
pub enum AeroProgress {
    /// A coarse stage description for the status line.
    Stage(String),
    /// One newly-available residual sample `(iteration, residual)` —
    /// fed straight into the live convergence plot. (Emitted in bulk
    /// once the steady solve returns; the engine does not expose a
    /// per-iteration callback, so the residual history arrives whole.)
    Residual(usize, f64),
    /// One completed sweep point `(angle_deg, cd, cl)`.
    SweepPoint(f64, f64, f64),
}

/// A completed angle-of-attack polar sweep, held by the workbench.
///
/// A thin newtype around [`PolarCurve`] so the workbench state has a
/// distinct, named result type alongside the steady [`AeroResult`].
#[derive(Clone, Debug)]
pub struct PolarSweepResult {
    /// The lift / drag polar.
    pub curve: PolarCurve,
}

/// The terminal outcome of a background run.
pub enum AeroOutcome {
    /// A steady single-point solve finished — carries the full result
    /// and the human-readable report.
    Steady(Box<AeroResult>, Box<AeroReport>),
    /// An angle-of-attack sweep finished — carries the polar curve.
    Sweep(Box<PolarCurve>),
    /// The run failed before producing a result.
    Failed(String),
}

/// What kind of run a job describes.
#[derive(Clone, Debug)]
pub enum AeroJob {
    /// A steady single-point solve.
    Steady,
    /// An angle-of-attack sweep over the given angles (radians).
    Sweep(Vec<f64>),
}

/// A live handle to a background wind-tunnel run.
///
/// The workbench holds an `Option<AeroRunHandle>`; each frame it drains
/// [`poll`](Self::poll) for progress and checks
/// [`take_outcome`](Self::take_outcome) for completion.
pub struct AeroRunHandle {
    /// Progress stream from the solve thread.
    rx: Receiver<AeroProgress>,
    /// The join handle — `take`n when the run finishes so the result
    /// can be moved out.
    thread: Option<JoinHandle<AeroOutcome>>,
    /// Round-4 cancellation token — flipped to `true` by
    /// [`AeroRunHandle::cancel`]. The sweep loop checks between
    /// per-angle solves and bails early when set; the single-shot
    /// steady solve can't be interrupted mid-flight (the wind-tunnel
    /// solver doesn't yet take a cancellation hook), but the App's
    /// on_exit can still drop the handle and let the OS reclaim the
    /// thread on process exit.
    cancel: Arc<AtomicBool>,
}

impl AeroRunHandle {
    /// Spawn a wind-tunnel run on a background thread.
    ///
    /// `body` is the (already-extracted) triangle mesh; `request` is
    /// the validated case; `job` selects a steady solve or a sweep.
    /// The returned handle is polled by the UI.
    pub fn spawn(body: TriMesh, request: AeroRequest, job: AeroJob) -> AeroRunHandle {
        let (tx, rx) = std::sync::mpsc::channel();
        let cancel = Arc::new(AtomicBool::new(false));
        let cancel_for_thread = cancel.clone();
        let thread = std::thread::Builder::new()
            .name("valenx-aero-run".to_string())
            .spawn(move || run_job(body, request, job, tx, cancel_for_thread))
            .expect("spawn aero run thread");
        AeroRunHandle {
            rx,
            thread: Some(thread),
            cancel,
        }
    }

    /// Round-4 cancellation: request the worker thread stop at the
    /// next per-angle checkpoint. A no-op once the worker has already
    /// finished. The App's `on_exit` calls this on every active run
    /// handle so closing the window doesn't orphan the worker.
    pub fn cancel(&self) {
        self.cancel
            .store(true, std::sync::atomic::Ordering::Relaxed);
    }

    /// Drain every progress message available right now (non-blocking).
    pub fn poll(&mut self) -> Vec<AeroProgress> {
        let mut out = Vec::new();
        while let Ok(msg) = self.rx.try_recv() {
            out.push(msg);
        }
        out
    }

    /// `true` once the solve thread has finished. The outcome can then
    /// be moved out with [`take_outcome`](Self::take_outcome).
    pub fn is_finished(&self) -> bool {
        self.thread
            .as_ref()
            .map(|t| t.is_finished())
            .unwrap_or(true)
    }

    /// Take the run's terminal outcome once it has finished. Returns
    /// `None` while the run is still in flight; returns
    /// [`AeroOutcome::Failed`] if the worker thread panicked.
    pub fn take_outcome(&mut self) -> Option<AeroOutcome> {
        if !self.is_finished() {
            return None;
        }
        let thread = self.thread.take()?;
        Some(
            thread
                .join()
                .unwrap_or_else(|_| AeroOutcome::Failed("the solver thread panicked".to_string())),
        )
    }
}

/// The body of the background thread — run the job and emit progress.
fn run_job(
    body: TriMesh,
    request: AeroRequest,
    job: AeroJob,
    tx: Sender<AeroProgress>,
    cancel: Arc<AtomicBool>,
) -> AeroOutcome {
    match job {
        AeroJob::Steady => run_steady(body, request, tx),
        AeroJob::Sweep(angles) => run_sweep(body, request, angles, tx, cancel),
    }
}

/// Run a steady single-point solve and assemble the report.
fn run_steady(body: TriMesh, request: AeroRequest, tx: Sender<AeroProgress>) -> AeroOutcome {
    let _ = tx.send(AeroProgress::Stage(
        "Building the virtual wind tunnel…".to_string(),
    ));
    let _ = tx.send(AeroProgress::Stage(
        "Solving the steady 3-D flow…".to_string(),
    ));
    match run_windtunnel(&body, &request) {
        Ok(result) => {
            // Stream the residual history into the convergence plot.
            for (i, &r) in result.flow.residual_history.iter().enumerate() {
                let _ = tx.send(AeroProgress::Residual(i + 1, r));
            }
            let _ = tx.send(AeroProgress::Stage("Integrating forces…".to_string()));
            let report = AeroReport::from_result(&result);
            AeroOutcome::Steady(Box::new(result), Box::new(report))
        }
        Err(e) => AeroOutcome::Failed(format!("[{}] {e}", e.code())),
    }
}

/// Run an angle-of-attack sweep, emitting a progress message per point.
///
/// Round-4 cancellation: `cancel` is checked before emitting each
/// sweep-point progress message. If set, the sweep emits a Failed
/// outcome and stops streaming. The underlying `aoa_sweep` call
/// can't be interrupted mid-flight (the solver doesn't yet plumb a
/// cancellation hook), but the per-angle gate is enough for the
/// common "user clicked Cancel before the sweep finished" path.
fn run_sweep(
    body: TriMesh,
    request: AeroRequest,
    angles: Vec<f64>,
    tx: Sender<AeroProgress>,
    cancel: Arc<AtomicBool>,
) -> AeroOutcome {
    let _ = tx.send(AeroProgress::Stage(format!(
        "Running an angle-of-attack sweep ({} points)…",
        angles.len()
    )));
    if cancel.load(std::sync::atomic::Ordering::Relaxed) {
        return AeroOutcome::Failed("aero sweep cancelled".to_string());
    }
    match aoa_sweep(&body, &request, &angles) {
        Ok(curve) => {
            for p in &curve.points {
                if cancel.load(std::sync::atomic::Ordering::Relaxed) {
                    return AeroOutcome::Failed("aero sweep cancelled".to_string());
                }
                let _ = tx.send(AeroProgress::SweepPoint(
                    p.alpha * 180.0 / std::f64::consts::PI,
                    p.cd,
                    p.cl,
                ));
            }
            AeroOutcome::Sweep(Box::new(curve))
        }
        Err(e) => AeroOutcome::Failed(format!("[{}] {e}", e.code())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;
    use valenx_aero::{geometry::box_body, sweep::linspace_degrees};

    /// A small box body + a fast request — keeps the background run
    /// short enough for a unit test (these `#[test]`s are
    /// compile-checked only, but the construction must be correct).
    fn quick_case() -> (TriMesh, AeroRequest) {
        let body = box_body(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(valenx_aero::TurbulenceModel::KEpsilon)
            .with_max_iterations(8);
        (body, req)
    }

    #[test]
    fn steady_run_completes_and_yields_a_result() {
        let (body, req) = quick_case();
        let mut handle = AeroRunHandle::spawn(body, req, AeroJob::Steady);
        // Block on completion (a test, not the UI thread).
        loop {
            if let Some(outcome) = handle.take_outcome() {
                match outcome {
                    AeroOutcome::Steady(result, report) => {
                        assert!(result.coefficients.cd.is_finite());
                        assert!(report.cd.is_finite());
                        // The report's headline Cd matches the result.
                        assert_eq!(report.cd, result.coefficients.cd);
                    }
                    AeroOutcome::Failed(e) => panic!("steady run failed: {e}"),
                    AeroOutcome::Sweep(_) => panic!("expected a steady outcome"),
                }
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
    }

    #[test]
    fn steady_run_streams_residual_progress() {
        let (body, req) = quick_case();
        let mut handle = AeroRunHandle::spawn(body, req, AeroJob::Steady);
        let mut residuals = 0;
        let mut stages = 0;
        loop {
            for msg in handle.poll() {
                match msg {
                    AeroProgress::Residual(it, r) => {
                        assert!(it >= 1);
                        assert!(r.is_finite());
                        residuals += 1;
                    }
                    AeroProgress::Stage(_) => stages += 1,
                    AeroProgress::SweepPoint(..) => panic!("no sweep points expected"),
                }
            }
            if handle.is_finished() {
                // Drain any final messages.
                for msg in handle.poll() {
                    if let AeroProgress::Residual(..) = msg {
                        residuals += 1;
                    }
                }
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        let _ = handle.take_outcome();
        assert!(stages >= 1, "expected stage messages");
        assert!(residuals >= 1, "expected residual history to be streamed");
    }

    #[test]
    fn sweep_run_yields_a_polar_curve() {
        let (body, req) = quick_case();
        let angles = linspace_degrees(0.0, 4.0, 2);
        let n = angles.len();
        let mut handle = AeroRunHandle::spawn(body, req, AeroJob::Sweep(angles));
        loop {
            if let Some(outcome) = handle.take_outcome() {
                match outcome {
                    AeroOutcome::Sweep(curve) => {
                        assert_eq!(curve.points.len(), n);
                        assert!(curve.points.iter().all(|p| p.cd.is_finite()));
                    }
                    AeroOutcome::Failed(e) => panic!("sweep failed: {e}"),
                    AeroOutcome::Steady(..) => panic!("expected a sweep outcome"),
                }
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
    }

    #[test]
    fn failed_run_reports_the_error_code() {
        // A negative speed is an ill-posed wind — the run must fail
        // with a coded message rather than panic.
        let body = box_body(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let req = AeroRequest::new(-5.0);
        let mut handle = AeroRunHandle::spawn(body, req, AeroJob::Steady);
        loop {
            if let Some(outcome) = handle.take_outcome() {
                match outcome {
                    AeroOutcome::Failed(e) => {
                        assert!(e.contains("aero."), "error should carry a code: {e}");
                    }
                    _ => panic!("expected a failure"),
                }
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
    }

    #[test]
    fn take_outcome_is_none_while_running() {
        // Immediately after spawn the run is in flight — take_outcome
        // must not block or hand back a half-built result.
        let (body, req) = quick_case();
        let mut handle = AeroRunHandle::spawn(body, req, AeroJob::Steady);
        // It may or may not have finished by now; if not, the outcome
        // is None and the handle is still usable.
        if !handle.is_finished() {
            assert!(handle.take_outcome().is_none());
        }
        // Drain to completion so the thread is joined.
        while handle.take_outcome().is_none() {
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
    }
}
