//! Live recording state for the desktop shell.
//!
//! Owned by `valenx-app`'s top-level state; the recorder gets ticked
//! every time a UI panel calls `record_action(...)`. Recording state
//! is held as a [`Recording`] sum type so the recorder cleanly
//! distinguishes "not currently recording" from "actively recording".

use crate::action::MacroAction;
use crate::Macro;

/// Recording state.
#[derive(Clone, Debug, Default)]
pub enum Recording {
    /// Not currently recording — calls to [`MacroRecorder::record`] are
    /// dropped on the floor.
    #[default]
    Off,
    /// Recording in progress. Holds the in-progress macro.
    Active(Macro),
}

/// The recorder itself — usually held as a `Mutex<MacroRecorder>` in
/// the desktop shell so every panel can append actions without
/// threading a `&mut` reference.
#[derive(Clone, Debug, Default)]
pub struct MacroRecorder {
    /// Current recording state.
    pub state: Recording,
    /// All previously stopped recordings, kept for the Macro Library
    /// panel.
    pub library: Vec<Macro>,
}

impl MacroRecorder {
    /// Build a fresh recorder with empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Begin recording. If a recording is already active, the in-
    /// progress macro is moved into the library and a new one is
    /// started.
    pub fn start_recording(&mut self, name: impl Into<String>) {
        if let Recording::Active(prev) = std::mem::take(&mut self.state) {
            self.library.push(prev);
        }
        self.state = Recording::Active(Macro::new(name));
    }

    /// Stop the current recording and return it. Returns `None` if no
    /// recording is active.
    pub fn stop_recording(&mut self) -> Option<Macro> {
        match std::mem::take(&mut self.state) {
            Recording::Active(m) => {
                self.library.push(m.clone());
                Some(m)
            }
            Recording::Off => None,
        }
    }

    /// Record one action. No-op when not recording.
    pub fn record(&mut self, action: MacroAction) {
        if let Recording::Active(m) = &mut self.state {
            m.push(action);
        }
    }

    /// True if currently recording.
    pub fn is_recording(&self) -> bool {
        matches!(self.state, Recording::Active(_))
    }

    /// Number of actions captured by the current recording (0 if not
    /// recording).
    pub fn current_action_count(&self) -> usize {
        match &self.state {
            Recording::Active(m) => m.actions.len(),
            Recording::Off => 0,
        }
    }

    /// Insert a previously-saved macro into the library (e.g. loaded
    /// from `~/.valenx/macros/`).
    pub fn add_to_library(&mut self, m: Macro) {
        self.library.push(m);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::action::{MacroAction, PanelId};

    #[test]
    fn new_recorder_is_off() {
        let r = MacroRecorder::new();
        assert!(!r.is_recording());
        assert_eq!(r.current_action_count(), 0);
    }

    #[test]
    fn start_then_stop_returns_macro() {
        let mut r = MacroRecorder::new();
        r.start_recording("X");
        r.record(MacroAction::SwitchPanel {
            panel_id: PanelId::Sketcher,
        });
        let m = r.stop_recording().unwrap();
        assert_eq!(m.name, "X");
        assert_eq!(m.actions.len(), 1);
        assert!(!r.is_recording());
    }

    #[test]
    fn record_when_off_drops() {
        let mut r = MacroRecorder::new();
        r.record(MacroAction::SwitchPanel {
            panel_id: PanelId::Mesh,
        });
        assert_eq!(r.current_action_count(), 0);
    }

    #[test]
    fn start_while_active_archives_previous() {
        let mut r = MacroRecorder::new();
        r.start_recording("First");
        r.record(MacroAction::SwitchPanel {
            panel_id: PanelId::Mesh,
        });
        r.start_recording("Second");
        // First macro was moved into the library.
        assert_eq!(r.library.len(), 1);
        assert_eq!(r.library[0].name, "First");
        assert_eq!(r.library[0].actions.len(), 1);
        // Second is fresh.
        assert_eq!(r.current_action_count(), 0);
    }

    #[test]
    fn stop_appends_to_library() {
        let mut r = MacroRecorder::new();
        r.start_recording("Lib");
        r.record(MacroAction::SwitchPanel {
            panel_id: PanelId::PartDesign,
        });
        let _ = r.stop_recording();
        assert_eq!(r.library.len(), 1);
    }
}
