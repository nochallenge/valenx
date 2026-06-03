//! # valenx-reinforcement
//!
//! Reinforcement workbench — parametric concrete rebar + cages. The
//! FreeCAD `Reinforcement` community workbench equivalent.
//!
//! Phase 33 of the FreeCAD-parity roadmap.
//!
//! # Surface
//!
//! - [`Rebar`] — one bar: diameter (mm), length (m), shape
//!   ([`RebarShape::Straight`] / `L` / `U` / `Hook` / `Spiral`),
//!   grade.
//! - [`RebarShape::to_polyline`] — emit the bar's centreline as a
//!   3D polyline.
//! - [`RebarCage`] — longitudinal bars + transverse hoops + spacing +
//!   cover.
//! - [`cage::generate_beam`] / [`cage::generate_column`] /
//!   [`cage::generate_slab`] — common production cage recipes.
//! - [`cage::to_mesh`] — tessellate a cage into one Tri3 mesh
//!   suitable for the viewport.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod cage;
pub mod error;
pub mod rebar;

pub use cage::{generate_beam, generate_column, generate_slab, to_mesh, RebarCage};
pub use error::{ErrorCategory, ReinforcementError};
pub use rebar::{Rebar, RebarGrade, RebarShape};
