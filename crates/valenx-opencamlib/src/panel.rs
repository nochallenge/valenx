//! UI panel state for the OpenCamLib algorithm demos.

use crate::cutter::Tool;
use crate::triangle::Triangle;

/// Workbench-panel state.
pub struct OpenCamLibPanelState {
    /// Loaded triangle surface.
    pub tris: Vec<Triangle>,
    /// Current tool spec.
    pub tool: Tool,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
    /// Number of XY samples taken in the last AdaptiveDropCutter run.
    pub last_sample_count: usize,
}

impl Default for OpenCamLibPanelState {
    fn default() -> Self {
        Self {
            tris: Vec::new(),
            tool: Tool::new(2.0, 10.0),
            last_status: None,
            last_error: None,
            last_sample_count: 0,
        }
    }
}

impl OpenCamLibPanelState {
    /// Empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Replace the triangle list.
    pub fn set_triangles(&mut self, tris: Vec<Triangle>) {
        self.last_status = Some(format!("loaded {} triangles", tris.len()));
        self.last_error = None;
        self.tris = tris;
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
