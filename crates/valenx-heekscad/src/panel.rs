//! UI panel envelope.

use crate::cam::Toolpath;
use crate::drawing::Drawing;

/// Workbench panel state.
pub struct HeeksCadPanelState {
    /// Current drawing.
    pub drawing: Drawing,
    /// Generated toolpaths (most recent first).
    pub toolpaths: Vec<Toolpath>,
    /// Status message.
    pub last_status: Option<String>,
    /// Error message.
    pub last_error: Option<String>,
}

impl Default for HeeksCadPanelState {
    fn default() -> Self {
        Self {
            drawing: Drawing::new(),
            toolpaths: Vec::new(),
            last_status: None,
            last_error: None,
        }
    }
}

impl HeeksCadPanelState {
    /// New empty.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a toolpath.
    pub fn add_toolpath(&mut self, tp: Toolpath) {
        self.last_status = Some(format!(
            "added toolpath `{}` ({} moves)",
            tp.op_name,
            tp.moves.len()
        ));
        self.last_error = None;
        self.toolpaths.push(tp);
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
