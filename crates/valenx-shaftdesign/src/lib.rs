//! # valenx-shaftdesign
//!
//! Closed-form stress analysis of a solid circular shaft of diameter
//! `d` carrying a steady bending moment `M` and torque `T`.
//!
//! ## What
//!
//! Given a [`ShaftSection`] (a positive diameter plus the applied
//! loads), this crate returns the elementary mechanics-of-materials
//! stresses at the outer fibre, where they are largest:
//!
//! - [`torsion::shear_stress`] — torsional shear stress
//!   `tau = 16 T / (pi d^3)`.
//! - [`bending::bending_stress`] — bending normal stress
//!   `sigma = 32 M / (pi d^3)`.
//! - [`combined::equivalent_torque`] — the equivalent torque
//!   `T_e = sqrt(M^2 + T^2)`, the single twisting moment producing the
//!   same maximum shear stress as the combined `M` and `T`.
//! - [`combined::equivalent_bending_moment`] — the equivalent bending
//!   moment `M_e = (M + sqrt(M^2 + T^2)) / 2`, the single bending
//!   moment producing the same maximum normal stress.
//! - [`combined::max_shear_stress`] — the combined maximum shear
//!   stress `tau_max = sqrt((sigma/2)^2 + tau^2)`, the radius of Mohr's
//!   circle for the uniaxial-bending-plus-torsion stress state, equal
//!   to `16 T_e / (pi d^3)`.
//! - [`combined::max_normal_stress`] — the combined maximum (principal)
//!   normal stress `sigma_1 = sigma/2 + sqrt((sigma/2)^2 + tau^2)`,
//!   equal to `32 M_e / (pi d^3)`.
//!
//! Section properties — [`ShaftSection::polar_section_modulus`] and
//! [`ShaftSection::section_modulus`] — are exposed for reuse.
//!
//! ## Model
//!
//! The bar is a prismatic, solid, circular shaft of homogeneous
//! linear-elastic isotropic material loaded within its elastic range.
//! Saint-Venant torsion gives the linear shear distribution
//! `tau(r) = T r / J` with the polar second moment
//! `J = pi d^4 / 32`; the outer-fibre value is
//! `tau = T (d/2) / J = 16 T / (pi d^3)`. Euler-Bernoulli bending gives
//! `sigma(y) = M y / I` with `I = pi d^4 / 64`; the outer-fibre value
//! is `sigma = M (d/2) / I = 32 M / (pi d^3)`.
//!
//! The combined point on the surface sees uniaxial normal stress
//! `sigma` (from bending) plus shear `tau` (from torsion). Mohr's
//! circle for that state has centre `sigma/2` and radius
//! `sqrt((sigma/2)^2 + tau^2)`, giving the maximum-shear-stress (Tresca)
//! result above. Substituting the fibre formulae shows the radius
//! equals `16 sqrt(M^2 + T^2) / (pi d^3)`, i.e. the same shear an
//! equivalent torque `T_e = sqrt(M^2 + T^2)` would produce alone — the
//! classic shaft-design equivalent-torque shortcut.
//!
//! Every stress scales as `1 / d^3`, so increasing the diameter is a
//! strong lever: doubling `d` cuts every stress to one eighth.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form models
//! for a single idealised section: no stress-concentration factors
//! (keyways, shoulders, holes, fillets), no fatigue / endurance limits,
//! no dynamic, shock, or fluctuating-load factors, no shear-deflection
//! or shaft-whirl effects, no hollow / stepped / non-circular sections,
//! and no material yield check beyond returning the stresses. It is
//! **NOT** a clinical, medical, or production engineering tool. Do not
//! use it to size a real load-bearing shaft: verify every result
//! against a qualified engineer and the governing design code.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bending;
pub mod combined;
pub mod error;
pub mod section;
pub mod torsion;

pub use bending::bending_stress;
pub use combined::{
    equivalent_bending_moment, equivalent_torque, max_normal_stress, max_shear_stress,
};
pub use error::{ErrorCategory, ShaftError};
pub use section::ShaftSection;
pub use torsion::shear_stress;
