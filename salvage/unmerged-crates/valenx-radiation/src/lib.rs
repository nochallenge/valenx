//! # valenx-radiation
//!
//! Closed-form thermal-radiation primitives: black-body and gray-body
//! emissive power, net surface-to-surroundings exchange, and
//! two-surface view-factor exchange.
//!
//! ## What
//!
//! Three small, dependency-free modules cover the textbook surface-
//! radiation toolkit:
//!
//! [`stefan_boltzmann`] — black-body emissive power `Eb = sigma T^4`,
//! its gray-surface scaling, and the inverse (radiative temperature
//! from a measured flux).
//!
//! [`graybody`] — net radiant exchange between a gray surface and large
//! surroundings, `Q = eps sigma A (T1^4 - T2^4)`, plus the area-free
//! flux form and the linearised radiation heat-transfer coefficient.
//!
//! [`view_factor`] — view-factor algebra (reciprocity
//! `A1 F12 = A2 F21`, the enclosure summation rule) and the
//! diffuse-gray two-surface enclosure exchange via the
//! series-resistance network.
//!
//! Every public function takes validated inputs and returns
//! `Result<f64, RadiationError>`; the [`error`] module's constructors
//! reject non-finite values and out-of-domain inputs (negative absolute
//! temperatures, non-positive areas, emissivities or view factors
//! outside `[0, 1]`).
//!
//! ## Model
//!
//! The governing equations, with `sigma = 5.670374419e-8 W m^-2 K^-4`
//! (CODATA, exact in SI since the 2019 redefinition of the kelvin):
//!
//! ```text
//! Eb   = sigma T^4                                   (black-body emissive power)
//! E    = eps Eb = eps sigma T^4                      (gray-surface emissive power)
//! Q    = eps sigma A (T1^4 - T2^4)                   (surface -> large surroundings)
//! q''  = eps sigma (T1^4 - T2^4)                     (the same, per unit area)
//! h_r  = eps sigma (T1 + T2)(T1^2 + T2^2)            (linearised coefficient)
//! A1 F12 = A2 F21                                    (view-factor reciprocity)
//! sum_j F_ij = 1                                     (enclosure summation rule)
//!
//!            sigma (T1^4 - T2^4)
//! Q12 = ------------------------------------------------------  (two diffuse-gray surfaces)
//!        (1-eps1)/(eps1 A1) + 1/(A1 F12) + (1-eps2)/(eps2 A2)
//! ```
//!
//! ## Honest scope
//!
//! This is research/educational grade: standard textbook closed-form
//! models (Incropera / Modest), validated in-crate against analytic
//! ground truth (the CODATA `sigma`, exact `T^4` scaling, reciprocity,
//! and the black-body and infinite-parallel-plate limits). It is NOT a
//! clinical/medical or production-certified engineering tool. It models
//! diffuse, gray, opaque surfaces only — no spectral or directional
//! dependence, no participating (absorbing/emitting/scattering) media,
//! no conduction or convection coupling, and no component tolerances,
//! safety factors, fatigue/thermal limits, or code compliance.
//!
//! ## Example
//!
//! ```rust
//! use valenx_radiation::{
//!     stefan_boltzmann::black_body_emissive_power,
//!     graybody::net_exchange_to_surroundings,
//!     view_factor::reciprocal_view_factor,
//! };
//!
//! // Black-body emissive power at 1000 K (~56.7 kW/m^2 textbook value).
//! let eb = black_body_emissive_power(1000.0).unwrap();
//! assert!((eb - 56_703.74).abs() < 1.0);
//!
//! // A 2 m^2, eps = 0.8 panel at 500 K losing heat to 300 K walls.
//! let q = net_exchange_to_surroundings(0.8, 2.0, 500.0, 300.0).unwrap();
//! assert!(q > 0.0); // net loss because the panel is hotter
//!
//! // View-factor reciprocity: a 1 m^2 convex body inside a 4 m^2 cavity.
//! let f21 = reciprocal_view_factor(1.0, 1.0, 4.0).unwrap();
//! assert!((f21 - 0.25).abs() < 1e-12);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod graybody;
pub mod stefan_boltzmann;
pub mod view_factor;

pub use error::{ErrorCategory, RadiationError, Result};
pub use graybody::{net_exchange_to_surroundings, net_flux_to_surroundings, radiation_coefficient};
pub use stefan_boltzmann::{
    black_body_emissive_power, gray_emissive_power, temperature_from_emissive_power, SIGMA,
};
pub use view_factor::{
    closing_view_factor, parallel_plates_flux, reciprocal_view_factor,
    self_view_factor_two_surface, two_surface_exchange,
};
