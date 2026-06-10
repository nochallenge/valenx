//! # valenx-springs
//!
//! Helical spring generator (compression, extension, torsion) with
//! tessellated swept solid output + canonical shear-modulus
//! stiffness formula.
//!
//! Phase 40 of the FreeCAD-parity roadmap. FreeCAD `Springs`
//! community workbench analogue.
//!
//! # Surface
//!
//! - [`SpringKind`] / [`SpringSpec`] / [`EndTreatment`] —
//!   parametric inputs.
//! - [`springs::compression_centerline`] /
//!   [`springs::extension_centerline`] /
//!   [`springs::torsion_centerline`] — 3D centreline polylines.
//! - [`springs::to_solid`] — wire-circle sweep returning a
//!   [`valenx_cad::Solid::Mesh`].
//! - [`springs::stiffness_n_per_mm`] — axial stiffness from
//!   `k = G d^4 / (8 D^3 n)`.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod spec;
pub mod springs;

pub use error::{ErrorCategory, SpringsError};
pub use spec::{EndTreatment, SpringKind, SpringSpec};
pub use springs::{
    compression_centerline, extension_centerline, spring_index, stiffness_n_per_mm, to_solid,
    torsion_centerline, wahl_factor,
};
