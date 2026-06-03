//! UI panel state envelope.

use crate::drawing::Drawing2D;

/// Workbench-panel state.
pub struct LibreCadPanelState {
    /// Active drawing.
    pub drawing: Drawing2D,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
}

impl Default for LibreCadPanelState {
    fn default() -> Self {
        Self {
            drawing: Drawing2D::new(),
            last_status: None,
            last_error: None,
        }
    }
}

impl LibreCadPanelState {
    /// Empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Replace the drawing.
    pub fn set_drawing(&mut self, d: Drawing2D) {
        self.last_status = Some(format!(
            "loaded drawing with {} entities, {} blocks, {} layers",
            d.entities.len(),
            d.blocks.len(),
            d.layers.len()
        ));
        self.last_error = None;
        self.drawing = d;
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
