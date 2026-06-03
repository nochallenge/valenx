//! UI state envelope for the HVAC panel.

use serde::{Deserialize, Serialize};

use crate::duct::Duct;
use crate::equipment::Equipment;

/// Workbench-panel state for HVAC.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct HvacPanelState {
    /// All ducts placed by the user.
    pub ducts: Vec<Duct>,
    /// All equipment instances placed by the user.
    pub equipment: Vec<(Equipment, (f64, f64, f64))>,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
}

impl HvacPanelState {
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
    fn default_state_is_empty() {
        let s = HvacPanelState::new();
        assert!(s.ducts.is_empty());
        assert!(s.equipment.is_empty());
    }
}
