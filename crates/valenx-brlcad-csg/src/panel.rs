//! UI panel envelope for the BRL-CAD CSG workbench.

use crate::csg::{CsgNode, SolidHandle};
use crate::error::BrlCadError;

/// Workbench panel.
#[derive(Default)]
pub struct BrlCadPanelState {
    /// Free-form text buffer (Lisp-style MGED form).
    pub text: String,
    /// Last parsed tree (`None` until the user clicks Evaluate).
    pub parsed: Option<CsgNode>,
    /// Last evaluation result.
    pub result: Option<SolidHandle>,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
}

impl BrlCadPanelState {
    /// New empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Parse `self.text` and evaluate. The result is stored in
    /// `self.result`; status / error flags reflect the outcome.
    pub fn evaluate(&mut self) -> Result<(), BrlCadError> {
        let tree = crate::csg::parse_mged(&self.text)?;
        let handle = crate::csg::evaluate(&tree)?;
        self.last_status = Some(format!(
            "evaluated tree -> canonical length {}",
            handle.canonical.len()
        ));
        self.last_error = None;
        self.parsed = Some(tree);
        self.result = Some(handle);
        Ok(())
    }

    /// Status setter.
    pub fn set_status(&mut self, s: impl Into<String>) {
        self.last_status = Some(s.into());
        self.last_error = None;
    }

    /// Error setter.
    pub fn set_error(&mut self, s: impl Into<String>) {
        self.last_error = Some(s.into());
        self.last_status = None;
    }
}
