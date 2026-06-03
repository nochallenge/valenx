//! Trait that lets a macro replay against any host that knows how to
//! perform the actions.
//!
//! Decoupling via this trait keeps the macro crate from depending on
//! the full `valenx-app` crate (which depends on egui + every
//! workbench), so unit tests can ship a tiny in-memory dispatcher.

use thiserror::Error;

use crate::action::MacroAction;

/// Replay-time errors.
#[derive(Debug, Error)]
pub enum DispatchError {
    /// The dispatcher doesn't know how to handle this action.
    #[error("unsupported macro action: {0}")]
    Unsupported(String),
    /// The host returned a runtime error while applying the action.
    #[error("dispatch failed: {0}")]
    Failed(String),
}

/// Receivers of [`MacroAction`] events.
///
/// The desktop shell implements this trait on its `ValenxApp` and
/// forwards each action to the same panel methods the UI calls. Tests
/// can implement it on a tiny struct that just records calls.
pub trait MacroDispatcher {
    /// Apply one action.
    ///
    /// # Errors
    ///
    /// Implementations return [`DispatchError`] if the action can't be
    /// applied. The replay loop aborts on the first failure.
    fn dispatch(&mut self, action: &MacroAction) -> Result<(), DispatchError>;
}

/// An in-memory dispatcher that just records every action it receives.
/// Used by unit tests and by the "Dry run" button in the Macro Library
/// panel.
#[derive(Default)]
pub struct RecordingDispatcher {
    /// Every action this dispatcher received, in order.
    pub received: Vec<MacroAction>,
}

impl RecordingDispatcher {
    /// Build a fresh empty dispatcher.
    pub fn new() -> Self {
        Self::default()
    }
}

impl MacroDispatcher for RecordingDispatcher {
    fn dispatch(&mut self, action: &MacroAction) -> Result<(), DispatchError> {
        self.received.push(action.clone());
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::action::{EntityKind, MacroAction, PanelId};

    #[test]
    fn recording_dispatcher_collects_actions() {
        let mut d = RecordingDispatcher::new();
        d.dispatch(&MacroAction::SwitchPanel {
            panel_id: PanelId::Sketcher,
        })
        .unwrap();
        d.dispatch(&MacroAction::AddSketchEntity {
            entity: EntityKind::Point { x: 1.0, y: 2.0 },
        })
        .unwrap();
        assert_eq!(d.received.len(), 2);
    }
}
