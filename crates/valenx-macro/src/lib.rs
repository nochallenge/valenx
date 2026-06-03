//! # valenx-macro
//!
//! Macro recorder + GUI hooks for Valenx.
//!
//! When users want to repeat a multi-step UI workflow — "Open Sketcher,
//! add a 50x50 rectangle, switch to Part Design, Pad 10mm" — they reach
//! for a macro. This crate records those clicks as a structured
//! [`Macro`] of [`MacroAction`] entries that can be:
//!
//! 1. **Replayed** against a fresh `ValenxApp` instance, deterministic.
//! 2. **Exported as Python** — each action turns into a call against
//!    the `valenx-py` module that ships with Phase 11.
//! 3. **Serialized as RON** — saved to `~/.valenx/macros/{name}.ron`
//!    and loaded next session.
//!
//! ## Recording flow
//!
//! The desktop shell (`valenx-app`) owns a [`MacroRecorder`] behind a
//! `parking_lot::Mutex`-free static (we use [`std::sync::OnceLock`] +
//! interior mutability so the recorder lives in the app's runtime
//! without forcing every call site to thread a `&mut` reference).
//!
//! ```no_run
//! use valenx_macro::{MacroRecorder, MacroAction, PanelId};
//! let mut rec = MacroRecorder::new();
//! rec.start_recording("MyWorkflow");
//! rec.record(MacroAction::SwitchPanel { panel_id: PanelId::Sketcher });
//! rec.record(MacroAction::SaveProject { path: "/tmp/p.valenx".into() });
//! let macro_ = rec.stop_recording().expect("recording active");
//! assert_eq!(macro_.actions.len(), 2);
//! ```
//!
//! ## Replay flow
//!
//! Replay is a trait — [`MacroDispatcher`] — so the macro crate doesn't
//! need to know about the concrete `ValenxApp` type. The desktop shell
//! implements the trait on its app and forwards each action to the
//! same panel methods the UI calls.
//!
//! ## Python export
//!
//! `Macro::to_python()` writes one Python statement per action,
//! prefixed with `import valenx`. The full script is:
//!
//! ```python
//! import valenx
//!
//! app = valenx.App.new()
//! app.switch_panel("sketcher")
//! app.save_project("/tmp/p.valenx")
//! ```
//!
//! The script is portable across operating systems and assumes only
//! that the `valenx` wheel (built by `valenx-py` / Phase 11) is
//! importable.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod action;
pub mod dispatcher;
pub mod export;
pub mod persist;
pub mod recorder;

pub use action::{EntityKind, MacroAction, PanelId};
pub use dispatcher::{DispatchError, MacroDispatcher};
pub use recorder::{MacroRecorder, Recording};

/// A captured macro: a name, an optional description, and the
/// ordered list of [`MacroAction`] entries to replay.
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct Macro {
    /// Human-readable macro name; also the filename stem when saved
    /// to `~/.valenx/macros/{name}.ron`.
    pub name: String,
    /// Optional free-form description shown in the Macro Library.
    pub description: String,
    /// Ordered actions to replay.
    pub actions: Vec<MacroAction>,
}

impl Macro {
    /// Build a new empty macro with the given name.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            description: String::new(),
            actions: Vec::new(),
        }
    }

    /// Append an action to the macro.
    pub fn push(&mut self, action: MacroAction) {
        self.actions.push(action);
    }

    /// Replay this macro against `dispatcher`. Returns the first error
    /// encountered (subsequent actions are not attempted).
    ///
    /// # Errors
    ///
    /// Propagates [`DispatchError`] from the dispatcher.
    pub fn replay<D: MacroDispatcher>(&self, dispatcher: &mut D) -> Result<(), DispatchError> {
        for action in &self.actions {
            dispatcher.dispatch(action)?;
        }
        Ok(())
    }

    /// Render this macro as a self-contained Python script that uses
    /// the `valenx` module (the wheel from `valenx-py`).
    pub fn to_python(&self) -> String {
        export::to_python(self)
    }

    /// Render this macro as RON (suitable for writing to disk).
    ///
    /// # Errors
    ///
    /// Bubbles up [`ron::Error`] for any encoding failure (very rare in
    /// practice — the data shape is a fixed enum / struct tree).
    pub fn to_ron(&self) -> Result<String, ron::Error> {
        ron::ser::to_string_pretty(self, ron::ser::PrettyConfig::new())
    }

    /// Parse a macro from RON text.
    ///
    /// # Errors
    ///
    /// Bubbles up [`ron::error::SpannedError`] on malformed input.
    pub fn from_ron(text: &str) -> Result<Self, ron::error::SpannedError> {
        ron::from_str(text)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::action::PanelId;

    #[test]
    fn new_macro_is_empty() {
        let m = Macro::new("Test");
        assert_eq!(m.name, "Test");
        assert!(m.actions.is_empty());
    }

    #[test]
    fn push_appends_actions() {
        let mut m = Macro::new("Push");
        m.push(MacroAction::SwitchPanel {
            panel_id: PanelId::Sketcher,
        });
        m.push(MacroAction::SaveProject {
            path: "/tmp/p.valenx".into(),
        });
        assert_eq!(m.actions.len(), 2);
    }

    #[test]
    fn ron_round_trip_preserves_actions() {
        let mut m = Macro::new("RonTrip");
        m.description = "A test macro".to_string();
        m.push(MacroAction::SwitchPanel {
            panel_id: PanelId::PartDesign,
        });
        let ron = m.to_ron().unwrap();
        let m2 = Macro::from_ron(&ron).unwrap();
        assert_eq!(m.name, m2.name);
        assert_eq!(m.description, m2.description);
        assert_eq!(m.actions.len(), m2.actions.len());
    }

    #[test]
    fn to_python_emits_import() {
        let m = Macro::new("Py");
        let py = m.to_python();
        assert!(py.contains("import valenx"));
    }
}
