//! # valenx-spreadsheet
//!
//! Spreadsheet workbench for Valenx (Phase 16 of the FreeCAD-parity
//! roadmap).
//!
//! ## What's in here
//!
//! - [`Cell`] / [`Sheet`] / [`Spreadsheet`] — data model (HashMap-keyed
//!   cells, named sheets, workbook-level lookup).
//! - [`CellRef`] — canonical cell address, parses `"Sheet.A1"` form
//!   round-trip with `to_string()`.
//! - Formula AST [`formula::Expr`] + hand-rolled recursive-descent
//!   [`parser::parse`] (no `nom`).
//! - [`evaluator::evaluate`] — recursive evaluator with built-in
//!   functions (sin / cos / sqrt / pow / if / min / max / etc.) and
//!   circular-reference detection.
//! - [`persist::SpreadsheetFile`] — RON envelope for round-tripping
//!   workbooks to disk.
//!
//! ## Feature-tree integration
//!
//! `valenx-feature-tree` adds a `Value` enum
//! (`valenx_feature_tree::feature::Value`) so numeric feature
//! parameters (Pad depth, Pocket depth, Revolve angle, ...) can hold
//! either a literal `f64` or a formula source string.
//! `FeatureTree::replay_with_spreadsheet` resolves expressions
//! against a [`Spreadsheet`] before replay.
//!
//! ## End-to-end example
//!
//! ```
//! use valenx_spreadsheet::{Cell, CellRef, Spreadsheet};
//!
//! let mut ss = Spreadsheet::new();
//! ss.add_sheet("Default");
//! let a1 = CellRef::parse("Default.A1").unwrap();
//! ss.set_cell(&a1, Cell::Number(50.0)).unwrap();
//! let a2 = CellRef::parse("Default.A2").unwrap();
//! ss.set_cell(&a2, Cell::Formula("Default.A1 * 2".into())).unwrap();
//! assert_eq!(ss.evaluate_cell(&a2).unwrap(), 100.0);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod cell;
pub mod error;
pub mod evaluator;
pub mod formula;
pub mod parser;
pub mod persist;
pub mod sheet;

pub use cell::{Cell, CellRef};
pub use error::SpreadsheetError;
pub use persist::SpreadsheetFile;
pub use sheet::{Sheet, Spreadsheet};
