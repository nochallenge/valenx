//! Live re-evaluation engine.
//!
//! The engine owns the OpenSCAD source text, the lex/parse cache,
//! the variable bindings (from the most recent eval), and the
//! cached [`valenx_cad::Solid`] result.  Every call to
//! [`Engine::set_source`] marks the cached result dirty; the next
//! [`Engine::evaluate`] re-runs the importer pipeline.
//!
//! This is the analogue of OpenSCAD's "auto-reload" feature.  No
//! native file-watcher: tests + headless callers drive evaluation
//! directly; the UI panel polls [`Engine::is_dirty`] each frame and
//! invokes evaluate on a button press.

use std::collections::HashMap;
use std::time::{SystemTime, UNIX_EPOCH};

use valenx_cad::Solid;
use valenx_openscad_import::{evaluate, lex, parse};

use crate::error::OpenScadCsgError;

/// Live OpenSCAD evaluation engine.
#[derive(Default)]
pub struct Engine {
    /// Current source string.
    source: String,
    /// Cached variable bindings — populated by the last eval.  v1
    /// keeps the importer's `Env` opaque and just records the
    /// *initial* variable assignments the source declares.
    variable_bindings: HashMap<String, f64>,
    /// Cached result from the last successful eval.
    cached_result: Option<Solid>,
    /// Wall-clock epoch (seconds) of the last successful eval.  Used
    /// to display "evaluated 3 s ago" in the panel.
    last_eval_time: Option<f64>,
    /// Set when [`Engine::set_source`] runs; cleared on the next
    /// successful evaluate.
    dirty: bool,
}

impl Engine {
    /// Build an empty engine.
    pub fn new() -> Self {
        Self::default()
    }

    /// Borrow the current source string.
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Replace the source.  Marks the cached result stale.
    pub fn set_source(&mut self, text: impl Into<String>) {
        let text = text.into();
        if text != self.source {
            self.source = text;
            self.dirty = true;
        }
    }

    /// True if the cached result is older than the current source.
    pub fn is_dirty(&self) -> bool {
        self.dirty
    }

    /// Re-lex / re-parse / re-evaluate the source.  Returns a
    /// reference to the cached Solid on success.
    pub fn evaluate(&mut self) -> Result<&Solid, OpenScadCsgError> {
        let toks = lex(&self.source)?;
        let ast = parse(&toks)?;
        // v1: shallow scan for top-level `Assign(name, Number(_))`
        // pairs.  Real binding capture lands in 52.5 when we expose
        // the importer's `Env`.
        self.variable_bindings.clear();
        for node in &ast {
            if let valenx_openscad_import::Ast::Assign(name, expr) = node {
                if let valenx_openscad_import::Ast::Number(v) = expr.as_ref() {
                    self.variable_bindings.insert(name.clone(), *v);
                }
            }
        }
        let solid = evaluate(&ast)?;
        self.cached_result = Some(solid);
        self.dirty = false;
        self.last_eval_time = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .ok()
            .map(|d| d.as_secs_f64());
        Ok(self.cached_result.as_ref().expect("just set"))
    }

    /// Borrow the cached Solid (the one returned by the last
    /// successful eval).
    pub fn cached(&self) -> Option<&Solid> {
        self.cached_result.as_ref()
    }

    /// Variable bindings (top-level numeric assignments) parsed out
    /// of the source on the last eval.
    pub fn variable_bindings(&self) -> &HashMap<String, f64> {
        &self.variable_bindings
    }

    /// Wall-clock seconds since epoch of the last successful eval.
    pub fn last_eval_time(&self) -> Option<f64> {
        self.last_eval_time
    }

    /// Convenience helper — evaluate but discard the borrow so the
    /// caller can immediately mutate the engine again.
    pub fn evaluate_owned(&mut self) -> Result<Solid, OpenScadCsgError> {
        self.evaluate()?;
        Ok(self.cached_result.clone().expect("just evaluated"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_engine_not_dirty() {
        let e = Engine::new();
        assert!(!e.is_dirty());
        assert!(e.cached().is_none());
    }

    #[test]
    fn set_source_dirties() {
        let mut e = Engine::new();
        e.set_source("cube([1, 1, 1]);");
        assert!(e.is_dirty());
        assert_eq!(e.source(), "cube([1, 1, 1]);");
    }

    #[test]
    fn evaluate_clears_dirty() {
        let mut e = Engine::new();
        e.set_source("cube([2, 2, 2]);");
        let s = e.evaluate().expect("ok");
        assert!(s.faces() > 0);
        assert!(!e.is_dirty());
        assert!(e.cached().is_some());
        assert!(e.last_eval_time().is_some());
    }

    #[test]
    fn set_source_same_text_no_dirty() {
        let mut e = Engine::new();
        e.set_source("cube([1, 1, 1]);");
        e.evaluate().expect("ok");
        assert!(!e.is_dirty());
        e.set_source("cube([1, 1, 1]);"); // identical
        assert!(!e.is_dirty()); // no change → no dirty
    }

    #[test]
    fn evaluate_captures_variable_bindings() {
        let mut e = Engine::new();
        e.set_source("x = 5; y = 10; cube([x, y, 1]);");
        e.evaluate().expect("ok");
        assert_eq!(e.variable_bindings().get("x"), Some(&5.0));
        assert_eq!(e.variable_bindings().get("y"), Some(&10.0));
    }

    #[test]
    fn bad_source_surfaces_error() {
        let mut e = Engine::new();
        e.set_source("cube(("); // syntax error
        assert!(e.evaluate().is_err());
        assert!(e.cached().is_none());
    }
}
