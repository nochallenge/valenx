//! # valenx-strainrosette
//!
//! Reduction of a rectangular 0/45/90 strain-gauge rosette to the full
//! in-plane strain and stress state.
//!
//! ## What
//!
//! Bonded strain gauges measure normal strain along their own axis only.
//! A *rectangular rosette* arranges three gauges at 0, 45, and 90
//! degrees so that, by inverting the strain-transformation law, the
//! complete Cartesian strain state can be recovered. This crate
//!
//! - reduces the three readings to `eps_x`, `eps_y`, `gamma_xy`
//!   ([`reduce`]),
//! - computes the principal strains and the orientation of the major
//!   principal axis ([`principal_strains`]), and
//! - maps the strain state to plane stress through 2D Hooke's law for an
//!   isotropic [`ElasticMaterial`] ([`ElasticMaterial::plane_stress`]).
//!
//! [`analyze`] runs all three stages at once.
//!
//! ## Model
//!
//! The normal strain along a line at angle `theta` from the x-axis is
//!
//! ```text
//! eps(theta) = eps_x cos^2(theta) + eps_y sin^2(theta)
//!            + gamma_xy sin(theta) cos(theta).
//! ```
//!
//! Evaluating at the three gauge angles and inverting gives the
//! closed-form reduction
//!
//! ```text
//! eps_x    = eps_0
//! eps_y    = eps_90
//! gamma_xy = 2 eps_45 - eps_0 - eps_90.
//! ```
//!
//! The principal strains are the eigenvalues of the symmetric strain
//! tensor `[[eps_x, gamma_xy/2], [gamma_xy/2, eps_y]]`,
//!
//! ```text
//! eps_1,2 = (eps_x + eps_y)/2
//!         +/- sqrt( ((eps_x - eps_y)/2)^2 + (gamma_xy/2)^2 ),
//! ```
//!
//! with the major axis at `theta_p = 0.5 atan2(gamma_xy, eps_x - eps_y)`.
//! Plane-stress Hooke's law closes the loop:
//!
//! ```text
//! sigma_x = E/(1 - nu^2) (eps_x + nu eps_y)
//! sigma_y = E/(1 - nu^2) (eps_y + nu eps_x)
//! tau_xy  = G gamma_xy,   G = E / (2 (1 + nu)).
//! ```
//!
//! ## Honest scope
//!
//! Research / educational grade. The crate implements the standard
//! textbook closed-form strain-transformation equations and the
//! isotropic linear-elastic (plane-stress) constitutive law, validated
//! against analytic ground truth in the unit tests. It is deliberately
//! small and self-contained:
//!
//! - it covers the *rectangular* (45-degree) rosette only — not delta
//!   (60-degree) or tee rosettes;
//! - it models small-strain, linear, isotropic, homogeneous elasticity
//!   under plane stress only — no plasticity, anisotropy, large strain,
//!   temperature compensation, transverse-sensitivity correction, or
//!   gauge-factor/lead-wire effects;
//! - it performs no measurement-uncertainty propagation.
//!
//! It is NOT a clinical, medical, or production engineering tool and
//! must not be used as the sole basis for any safety-relevant or
//! load-bearing decision.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod analysis;
pub mod error;
pub mod material;
pub mod rosette;

pub use analysis::{analyze, RosetteAnalysis};
pub use error::{ErrorCategory, RosetteError};
pub use material::{ElasticMaterial, PlaneStress};
pub use rosette::{
    principal_from_readings, principal_strains, reduce, CartesianStrain, PrincipalStrain,
    RosetteReadings,
};
