//! UI state envelope for the Piping panel.

use serde::{Deserialize, Serialize};

use crate::system::Piping;

/// Workbench-panel state for the piping system.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct PipingPanelState {
    /// Active piping network.
    pub network: Piping,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
}

impl PipingPanelState {
    /// Empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record success.
    pub fn set_status(&mut self, msg: impl Into<String>) {
        self.last_status = Some(msg.into());
        self.last_error = None;
    }

    /// Record failure.
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
        let p = PipingPanelState::new();
        assert!(p.network.sections.is_empty());
        assert!(p.last_status.is_none());
    }
}
