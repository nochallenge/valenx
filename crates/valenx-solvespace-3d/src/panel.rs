//! Workbench-panel state envelope for the SolveSpace-3D solver.

use crate::sketch::Sketch3D;
use crate::solver::{SolverReport, SolverStatus};

/// UI-visible solver panel state.
///
/// The panel owns one [`Sketch3D`] plus the most-recent solver report.
/// Bringing this struct into `valenx-app` only requires wrapping it in
/// the usual `*_StateLocal` and registering a menu entry — no FileDialog,
/// no live windowing.
pub struct SolveSpace3DPanelState {
    /// The active sketch.
    pub sketch: Sketch3D,
    /// Most-recent solver report (`None` until `solve_now` runs).
    pub last_report: Option<SolverReport>,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error string.
    pub last_error: Option<String>,
}

impl Default for SolveSpace3DPanelState {
    fn default() -> Self {
        Self {
            sketch: Sketch3D::new(),
            last_report: None,
            last_status: None,
            last_error: None,
        }
    }
}

impl SolveSpace3DPanelState {
    /// Empty state with a fresh `Sketch3D`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Run the solver, recording the report (or error) on `self`.
    pub fn solve_now(&mut self) {
        match self.sketch.solve() {
            Ok(rep) => {
                self.last_status = Some(format!(
                    "{:?} in {} iters — residual {:.3e}",
                    rep.status, rep.iterations, rep.residual_norm
                ));
                self.last_error = None;
                self.last_report = Some(rep);
            }
            Err(e) => {
                self.last_error = Some(e.to_string());
                self.last_status = None;
            }
        }
    }

    /// True if the last solve converged.
    pub fn last_converged(&self) -> bool {
        matches!(
            self.last_report.as_ref().map(|r| r.status),
            Some(SolverStatus::Converged)
        )
    }

    /// Record a status message.
    pub fn set_status(&mut self, msg: impl Into<String>) {
        self.last_status = Some(msg.into());
        self.last_error = None;
    }

    /// Record an error message.
    pub fn set_error(&mut self, msg: impl Into<String>) {
        self.last_error = Some(msg.into());
        self.last_status = None;
    }
}
