//! # valenx-openscad-import
//!
//! OpenSCAD `.scad` importer — lex + parse + evaluate a useful subset
//! of the OpenSCAD language and produce a [`valenx_cad::Solid`].
//!
//! Phase 45 of the FreeCAD-parity roadmap. FreeCAD `OpenSCAD`
//! community workbench equivalent.
//!
//! # Surface
//!
//! - [`lex::lex`] — string → [`lex::Token`] stream, comments skipped.
//! - [`parse::parse`] — token stream → [`ast::Ast`] statement list.
//! - [`eval::evaluate`] — statement list → [`valenx_cad::Solid`].
//! - [`panel::OpenScadImportPanelState`] — UI state envelope used by
//!   the Valenx app.
//! - [`persist::to_ron_string`] / [`persist::from_ron_str`] — versioned
//!   panel-state persistence.
//!
//! # Supported subset
//!
//! Primitives: `cube`, `sphere`, `cylinder`.
//! Booleans: `union`, `difference`, `intersection`.
//! Transforms: `translate`, `rotate`, `scale` (uniform only in v1).
//! Statements: variable assignment, finite `for` loops, bare blocks.
//!
//! Unsupported (raises `OpenScadError::Eval`):
//! `module` / `function` definitions, `if`/`else`, `linear_extrude`,
//! `rotate_extrude`, `mirror`, `intersection_for`, lookup tables, the
//! `*` / `#` / `%` modifier prefixes, and any callee not listed above.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod ast;
pub mod error;
pub mod eval;
pub mod lex;
pub mod panel;
pub mod parse;
pub mod persist;

pub use ast::{Ast, BinOp};
pub use error::{ErrorCategory, OpenScadError};
pub use eval::{evaluate, Env, Value};
pub use lex::{is_keyword, lex, Token};
pub use panel::OpenScadImportPanelState;
pub use parse::parse;
pub use persist::{from_ron_str, to_ron_string, PanelFile, VERSION};

/// One-shot helper: load + lex + parse + evaluate an entire `.scad`
/// source string. Returns the final Solid (the implicit union of
/// every top-level shape).
pub fn import_str(src: &str) -> Result<valenx_cad::Solid, OpenScadError> {
    let toks = lex(src)?;
    let ast = parse(&toks)?;
    evaluate(&ast)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn end_to_end_cube() {
        let s = import_str("cube([4, 4, 4]);").expect("ok");
        assert!(s.faces() > 0);
    }
}
