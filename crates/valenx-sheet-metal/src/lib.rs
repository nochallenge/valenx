//! # valenx-sheet-metal
//!
//! Sheet Metal workbench — parametric sheet metal with flanges,
//! bends, and unfold. The FreeCAD `SheetMetal` community workbench
//! equivalent.
//!
//! Phase 34 of the FreeCAD-parity roadmap.
//!
//! # Surface
//!
//! - [`Sheet`] — outline (2D polygon) + thickness + material + k-
//!   factor (neutral-axis fraction for bend allowance).
//! - [`Bend`] — bend line (start + end in sheet-local 2D) + angle +
//!   inside_radius.
//! - [`Flange`] — flange at one outline edge + length + angle.
//! - [`Sheet::add_bend`] / [`Sheet::add_flange`] — recipe-style
//!   accumulators returning a new [`Sheet`].
//! - [`Sheet::to_solid`] — extrude as a thick plate (`Solid::Mesh`),
//!   then apply each bend by rotating downstream vertices around the
//!   bend line (v1 splits the sheet at the bend line; subdivides only
//!   the planar plate before bends are applied).
//! - [`Sheet::unfold`] — flatten the sheet to a 2D pattern using the
//!   k-factor bend allowance.
//! - [`Sheet::cutout`] — subtract a polygonal hole.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bend;
pub mod error;
pub mod sheet;

pub use bend::{Bend, Flange};
pub use error::{ErrorCategory, SheetMetalError};
pub use sheet::{Sheet, SheetMaterial};
