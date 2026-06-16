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
//! - [`springs::spring_index`] / [`springs::wahl_factor`] /
//!   [`springs::shear_stress_mpa`] — the coil index `C = D/d`, its Wahl
//!   curvature-correction factor, and the resulting wire shear stress
//!   `τ = K_w · 8 F D / (π d^3)` under an axial load.
//! - [`springs::deflection_mm`] / [`springs::stored_energy_nmm`] — the
//!   Hooke load response `δ = F / k` and the stored strain energy
//!   `U = ½ F δ = ½ F² / k`.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod spec;
pub mod springs;

pub use error::{ErrorCategory, SpringsError};
pub use spec::{EndTreatment, SpringKind, SpringSpec};
pub use springs::{
    compression_centerline, deflection_mm, extension_centerline, shear_stress_mpa, spring_index,
    stiffness_n_per_mm, stored_energy_nmm, to_solid, torsion_centerline, wahl_factor,
};
