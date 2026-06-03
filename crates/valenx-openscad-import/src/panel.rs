//! UI state envelope for the OpenSCAD import panel.
//!
//! The Valenx app drives a text-input path + Import button.  No native
//! file dialog (FileDialog is forbidden in tests; the panel
//! deliberately takes a plain string so it stays unit-testable).

use serde::{Deserialize, Serialize};

/// Workbench-panel state for the OpenSCAD importer.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct OpenScadImportPanelState {
    /// Text input — path to a `.scad` file the user wants to import.
    pub file_path_input: String,
    /// Last status message ("imported `cube.scad`").
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
    /// Last face count of the imported solid (cached for UI display).
    pub last_face_count: usize,
}

impl OpenScadImportPanelState {
    /// New, empty panel state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a successful import.
    pub fn set_status(&mut self, msg: impl Into<String>, faces: usize) {
        self.last_status = Some(msg.into());
        self.last_error = None;
        self.last_face_count = faces;
    }

    /// Record an error.
    pub fn set_error(&mut self, msg: impl Into<String>) {
        self.last_error = Some(msg.into());
        self.last_status = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_state_is_empty() {
        let s = OpenScadImportPanelState::new();
        assert!(s.file_path_input.is_empty());
        assert!(s.last_status.is_none());
        assert!(s.last_error.is_none());
        assert_eq!(s.last_face_count, 0);
    }

    #[test]
    fn set_status_clears_error() {
        let mut s = OpenScadImportPanelState::new();
        s.set_error("bad");
        s.set_status("ok", 42);
        assert!(s.last_error.is_none());
        assert_eq!(s.last_status.as_deref(), Some("ok"));
        assert_eq!(s.last_face_count, 42);
    }

    #[test]
    fn set_error_clears_status() {
        let mut s = OpenScadImportPanelState::new();
        s.set_status("ok", 42);
        s.set_error("oops");
        assert!(s.last_status.is_none());
        assert_eq!(s.last_error.as_deref(), Some("oops"));
    }
}
