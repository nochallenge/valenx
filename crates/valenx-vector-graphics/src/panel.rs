//! UI panel envelope.

use crate::entity::VectorEntity;

/// Workbench panel state.
#[derive(Default)]
pub struct VectorPanelState {
    /// Current entity list.
    pub entities: Vec<VectorEntity>,
    /// Selected entity index.
    pub selection: Option<usize>,
    /// Status message.
    pub last_status: Option<String>,
    /// Error message.
    pub last_error: Option<String>,
}

impl VectorPanelState {
    /// New empty.
    pub fn new() -> Self {
        Self::default()
    }

    /// Replace entity list.
    pub fn set_entities(&mut self, e: Vec<VectorEntity>) {
        self.last_status = Some(format!("loaded {} entities", e.len()));
        self.last_error = None;
        self.entities = e;
        self.selection = None;
    }

    /// Add one entity.
    pub fn add(&mut self, e: VectorEntity) {
        self.last_status = Some(format!("added {}", e.kind()));
        self.last_error = None;
        self.entities.push(e);
    }

    /// Select an entity.
    pub fn select(&mut self, i: Option<usize>) {
        self.selection = i;
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
