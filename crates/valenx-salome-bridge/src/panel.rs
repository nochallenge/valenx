//! UI panel envelope for the Salome-style platform — exposes the
//! [`Study`] browser plus diagnostic state.

use crate::study::{NodeId, Study};

/// Workbench panel state.
pub struct SalomePanelState {
    /// The active study.
    pub study: Study,
    /// Optional id of the user-selected node.
    pub selection: Option<NodeId>,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
}

impl Default for SalomePanelState {
    fn default() -> Self {
        Self {
            study: Study::new(),
            selection: None,
            last_status: None,
            last_error: None,
        }
    }
}

impl SalomePanelState {
    /// New empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Select a node id (or clear with `None`).
    pub fn select(&mut self, id: Option<NodeId>) {
        self.selection = id;
        if let Some(id) = id {
            self.last_status = Some(format!("selected node {id}"));
            self.last_error = None;
        }
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
