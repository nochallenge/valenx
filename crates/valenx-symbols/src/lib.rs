//! # valenx-symbols
//!
//! Schematic symbol library — electrical, hydraulic, and pneumatic
//! glyphs as canonical SVG paths, plus a [`Schematic`] composed of
//! placed symbols and wires and an SVG renderer.
//!
//! Phase 37 of the FreeCAD-parity roadmap. FreeCAD `Symbols Library`
//! community workbench equivalent.
//!
//! # Surface
//!
//! - [`SymbolKind`] / [`SymbolFamily`] — 21 standard glyphs across 3
//!   families.
//! - [`SymbolKind::to_svg_path`] — canonical SVG path string per glyph.
//! - [`Schematic`] / [`PlacedSymbol`] / [`Wire`] — schematic data
//!   model.
//! - [`schematic::to_svg`] — render schematic to standalone SVG.
//! - [`persist::to_ron_string`] / [`persist::from_ron_str`] —
//!   round-trip a schematic through the [`persist::SchematicFile`]
//!   RON envelope (versioned).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod persist;
pub mod schematic;
pub mod symbol;

pub use error::{ErrorCategory, SymbolError};
pub use persist::{SchematicFile, VERSION};
pub use schematic::{to_svg, PlacedSymbol, Schematic, Wire};
pub use symbol::{SymbolFamily, SymbolKind};
