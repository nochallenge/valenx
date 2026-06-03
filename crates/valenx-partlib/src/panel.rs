//! UI state envelope for the Parts Library panel.

use serde::{Deserialize, Serialize};

/// Workbench-panel state for the parts library.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct PartLibPanelState {
    /// Text input — library root path the user is browsing.
    pub root_input: String,
    /// Text input — name to register the next install under.
    pub install_name_input: String,
    /// Text input — source file path for the next install.
    pub install_file_input: String,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
    /// Number of parts currently in the loaded library.
    pub last_count: usize,
}

impl PartLibPanelState {
    /// New, empty panel state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a successful action.
    pub fn set_status(&mut self, msg: impl Into<String>, count: usize) {
        self.last_status = Some(msg.into());
        self.last_error = None;
        self.last_count = count;
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
    fn default_is_empty() {
        let s = PartLibPanelState::new();
        assert!(s.last_status.is_none());
        assert_eq!(s.last_count, 0);
    }
}
