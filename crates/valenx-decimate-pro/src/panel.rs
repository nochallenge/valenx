//! UI state envelope for the Decimate-Pro panel.

use serde::{Deserialize, Serialize};

/// Mode dropdown — picks which decimate-pro variant the panel runs.
#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub enum Mode {
    /// Plain valenx-mesh QEM.
    #[default]
    Standard,
    /// Curvature-weighted QEM.
    CurvatureWeighted,
    /// UV-preserving.
    UvPreserving,
    /// Feature-aware.
    FeatureAware,
}

/// Workbench-panel state for decimate-pro.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DecimateProPanelState {
    /// Active mode.
    pub mode: Mode,
    /// Target fraction (0..1) of vertices to keep.
    pub target_fraction: f64,
    /// Curvature weight for the curvature-weighted mode.
    pub curvature_weight: f64,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
    /// Number of nodes in the decimated mesh (for UI display).
    pub last_node_count: usize,
}

impl Default for DecimateProPanelState {
    fn default() -> Self {
        Self {
            mode: Mode::Standard,
            target_fraction: 0.5,
            curvature_weight: 0.5,
            last_status: None,
            last_error: None,
            last_node_count: 0,
        }
    }
}

impl DecimateProPanelState {
    /// Empty default panel state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record success.
    pub fn set_status(&mut self, msg: impl Into<String>, nodes: usize) {
        self.last_status = Some(msg.into());
        self.last_error = None;
        self.last_node_count = nodes;
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
    fn default_panel_has_sane_initial_values() {
        let s = DecimateProPanelState::new();
        assert_eq!(s.mode, Mode::Standard);
        assert!((s.target_fraction - 0.5).abs() < 1e-9);
        assert!((s.curvature_weight - 0.5).abs() < 1e-9);
    }
}
