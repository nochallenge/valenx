//! # valenx-fasteners
//!
//! Fasteners workbench — standard parts library for bolts, nuts, and
//! washers. The FreeCAD `Fasteners` community workbench equivalent.
//!
//! Phase 36 of the FreeCAD-parity roadmap.
//!
//! # Surface
//!
//! - [`BoltSpec`] / [`bolt::iso4017_hex_table`] /
//!   [`bolt::ansi_b18_2_1_table`] / [`bolt::to_solid`] — bolts.
//! - [`NutSpec`] / [`nut::iso4032_hex_table`] / [`nut::to_solid`] —
//!   nuts.
//! - [`WasherSpec`] / [`washer::iso7089_table`] /
//!   [`washer::to_solid`] — washers.
//!
//! Every `to_solid` returns a [`valenx_cad::Solid::Mesh`] — these are
//! visual parametric placeholders, not full BRep parts. They're fine
//! for assembly viewport rendering, bill-of-materials counts, and
//! TechDraw callouts (the [`valenx_feature_tree::threads::ThreadSpec`]
//! travels along with the bolt for the actual thread callout).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bolt;
pub mod error;
pub mod nut;
pub mod washer;

pub use bolt::{BoltKind, BoltSpec};
pub use error::{ErrorCategory, FastenerError};
pub use nut::NutSpec;
pub use washer::WasherSpec;
