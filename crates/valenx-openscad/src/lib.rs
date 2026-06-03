//! # valenx-openscad
//!
//! OpenSCAD CSG scripting paradigm — live re-evaluation engine that
//! owns the AST + variable bindings + the cached evaluation result.
//!
//! Phase 52 of the FreeCAD-parity roadmap (the
//! "OpenSCAD as a first-class language" feature).
//!
//! # Difference from `valenx-openscad-import`
//!
//! `valenx-openscad-import` (Phase 45) was the **one-shot importer**:
//! read a `.scad` file, lex + parse + evaluate once, return a
//! [`valenx_cad::Solid`].  This crate wraps that pipeline in an
//! [`engine::Engine`] that:
//!
//! - Holds the current source text and the cached AST.
//! - Tracks a dirty flag — every [`engine::Engine::set_source`] marks
//!   the cached result stale.
//! - Re-evaluates on demand via [`engine::Engine::evaluate`].
//! - Surfaces the most-recent error so a UI panel can show it next to
//!   the editor.
//!
//! # Expanded primitive set
//!
//! [`prims`] re-exports the Phase 45 [`valenx_cad`] primitive helpers
//! and adds the OpenSCAD-named ones (`polyhedron`, `square`, `circle`,
//! `polygon`, `text`).  2-D primitives are emitted as zero-height
//! prisms in the +Z plane so they round-trip through the BRep pipeline
//! the same way as 3-D primitives.
//!
//! # Modifier prefixes
//!
//! [`modifiers`] models the OpenSCAD `*` / `!` / `#` / `%` source
//! prefixes — disable / root / highlight / transparent.  The engine
//! filters / colours children based on the flag set; the lexer is in
//! `valenx-openscad-import` and only reports the prefix.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod engine;
pub mod error;
pub mod modifiers;
pub mod panel;
pub mod persist;
pub mod prims;

pub use engine::Engine;
pub use error::{ErrorCategory, OpenScadCsgError};
pub use modifiers::Modifier;
pub use panel::OpenScadPanelState;
pub use persist::{from_ron_str, to_ron_string, PanelFile, VERSION};
