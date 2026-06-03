//! UI panel state envelope.

use crate::triangle::TriMesh;

/// Workbench-panel state.
#[derive(Default)]
pub struct LibiglPanelState {
    /// Loaded mesh.
    pub mesh: TriMesh,
    /// Last LSCM / ARAP UVs.
    pub last_uvs: Vec<[f64; 2]>,
    /// Last geodesic field.
    pub last_geodesic: Vec<f64>,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
}

impl LibiglPanelState {
    /// Empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record status.
    pub fn set_status(&mut self, msg: impl Into<String>) {
        self.last_status = Some(msg.into());
        self.last_error = None;
    }

    /// Record error.
    pub fn set_error(&mut self, msg: impl Into<String>) {
        self.last_error = Some(msg.into());
        self.last_status = None;
    }
}
