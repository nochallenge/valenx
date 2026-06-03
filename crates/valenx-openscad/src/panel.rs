//! UI state envelope for the OpenSCAD CSG panel.
//!
//! The panel holds an [`crate::engine::Engine`] handle plus the
//! current editor text (separate from `engine.source()` so the user
//! can type WITHOUT triggering re-eval on every keystroke), a status
//! line, and the most-recent error message.

use crate::engine::Engine;

/// Workbench-panel state for the OpenSCAD CSG engine.
pub struct OpenScadPanelState {
    /// Multi-line editor buffer.  Sync into the engine via
    /// [`OpenScadPanelState::commit_source`].
    pub editor_text: String,
    /// The live evaluation engine.
    pub engine: Engine,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
}

impl Default for OpenScadPanelState {
    fn default() -> Self {
        Self {
            editor_text: "// OpenSCAD source\ncube([10, 10, 10]);\n".into(),
            engine: Engine::new(),
            last_status: None,
            last_error: None,
        }
    }
}

impl OpenScadPanelState {
    /// Empty default state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Push the editor text into the engine.  Doesn't evaluate.
    pub fn commit_source(&mut self) {
        self.engine.set_source(self.editor_text.clone());
    }

    /// Evaluate the engine, surfacing success or failure into the
    /// status/error fields.
    pub fn evaluate(&mut self) {
        self.commit_source();
        match self.engine.evaluate() {
            Ok(s) => {
                let faces = s.faces();
                self.last_status = Some(format!(
                    "Evaluated — {faces} faces, {} vars",
                    self.engine.variable_bindings().len()
                ));
                self.last_error = None;
            }
            Err(e) => {
                self.last_error = Some(e.to_string());
                self.last_status = None;
            }
        }
    }

    /// Record a successful action message.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_text_evaluates_clean() {
        let mut p = OpenScadPanelState::new();
        p.evaluate();
        assert!(p.last_status.is_some());
        assert!(p.last_error.is_none());
    }

    #[test]
    fn evaluate_propagates_errors() {
        let mut p = OpenScadPanelState::new();
        p.editor_text = "cube((".into(); // garbage
        p.evaluate();
        assert!(p.last_error.is_some());
        assert!(p.last_status.is_none());
    }
}
