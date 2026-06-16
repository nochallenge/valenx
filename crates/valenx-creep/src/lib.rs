//! # valenx-creep
//!
//! Closed-form **creep and stress-rupture** models for high-temperature
//! materials engineering — the Larson-Miller time-temperature parameter
//! and the Norton-Bailey secondary-creep law.
//!
//! ## What
//!
//! When a metal carries load at a high fraction of its melting
//! temperature it deforms slowly and continuously (creep) and will
//! eventually fail (stress rupture) even below its short-term yield
//! strength. This crate provides the two classic textbook tools used to
//! reason about that behaviour:
//!
//! - [`larson_miller`] — the **Larson-Miller parameter**
//!   `LMP = T * (C + log10(t_r))`, which collapses many
//!   temperature / time rupture tests onto a single master curve, plus
//!   the inverse solve `t_r = 10^(LMP / T - C)` for the time to rupture
//!   at a service temperature (and a temperature solve).
//! - [`norton`] — the **Norton-Bailey** steady-state creep law
//!   `epsilon_dot = A * sigma^n`, with an optional Arrhenius
//!   temperature dependence `A = A0 * exp(-Q / (R T))`, plus the inverse
//!   solve `sigma = (epsilon_dot / A)^(1/n)` for the stress that yields
//!   a target creep rate.
//!
//! Every fallible entry point validates its inputs and returns a typed
//! [`CreepError`] rather than a silent `NaN`.
//!
//! ```
//! use valenx_creep::{larson_miller, norton};
//!
//! // Time to rupture from a master-curve LMP at the operating point.
//! let life_h = larson_miller::rupture_time_hours(27_000.0, 1000.0, 20.0).unwrap();
//! assert!(life_h > 0.0);
//!
//! // Steady-state creep rate under load.
//! let law = norton::NortonLaw::new(1.0e-10, 5.0).unwrap();
//! let rate = law.rate_at(120.0).unwrap();
//! assert!(rate > 0.0);
//! ```
//!
//! ## Model
//!
//! The Larson-Miller relation is the empirical observation that a fixed
//! amount of creep damage (hence rupture) is reached along a locus of
//! constant `T * (C + log10(t_r))`: a part lasts a short time when hot
//! or a long time when cool. With the constant `C` (`~20` for many
//! ferrous alloys, [`larson_miller::DEFAULT_C`]) and an LMP read off the
//! material's master curve at the design stress, the relation predicts
//! the rupture life at any service temperature.
//!
//! Norton's law describes the **secondary** (steady-state) stage of
//! creep, where the strain rate is roughly constant and obeys a power
//! law `epsilon_dot = A * sigma^n` in the applied stress. The exponent
//! `n` (`~3` to `~8` for metals) is the slope of the
//! `log(rate)`-versus-`log(stress)` line; temperature enters through
//! the coefficient `A`, optionally via an Arrhenius factor.
//!
//! ## Honest scope
//!
//! This crate is **research / educational grade**. It implements the
//! standard textbook closed-form and numerical models from the
//! high-temperature-materials literature and nothing more. It ships no
//! material data — you must supply LMP, the constant `C`, and the
//! Norton constants `A`, `n`, `A0`, `Q` from a qualified source — and
//! it does not model scatter, multiaxial stress states, oxidation,
//! microstructural evolution, primary / tertiary creep stages, or the
//! validated extrapolation bounds of a real master curve. It is **not**
//! a clinical / medical tool and **not** a production engineering tool:
//! never use it to qualify, certify, or set the service life of real
//! load-bearing hardware.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod larson_miller;
pub mod norton;

pub use error::{CreepError, ErrorCategory};
pub use larson_miller::{
    larson_miller_parameter, rupture_temperature_k, rupture_time_hours, RupturePoint, DEFAULT_C,
};
pub use norton::{norton_creep_rate, norton_stress_for_rate, NortonLaw, GAS_CONSTANT_J_PER_MOL_K};

#[cfg(test)]
mod tests {
    //! Cross-module integration checks that exercise the public surface
    //! the way a caller would, combining both models.

    use super::*;

    /// Tolerance for floating-point comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn lmp_round_trips_via_reexports() {
        // Build LMP from a rupture point and invert it back, using only
        // the crate-root re-exports.
        let lmp = larson_miller_parameter(1050.0, 3_000.0, DEFAULT_C).unwrap();
        let back = rupture_time_hours(lmp, 1050.0, DEFAULT_C).unwrap();
        assert!((back - 3_000.0).abs() / 3_000.0 < 1e-9, "got {back}");
    }

    #[test]
    fn rupture_point_and_norton_combine() {
        // A hotter point of equal LMP ruptures sooner; an independent
        // Norton law reports a finite steady-state rate at the design
        // stress. The two models are independent but used together.
        let p = RupturePoint::new(1000.0, 1_000.0).unwrap();
        let lmp = p.parameter(DEFAULT_C).unwrap();
        let hotter = rupture_time_hours(lmp, 1100.0, DEFAULT_C).unwrap();
        assert!(hotter < p.time_hours, "hotter ruptures sooner: {hotter}");

        let law = NortonLaw::new(1.0e-12, 5.0).unwrap();
        let rate = law.rate_at(100.0).unwrap();
        // A = 1e-12, sigma = 100, n = 5 → 1e-2.
        assert!((rate - 1.0e-2).abs() < EPS, "got {rate}");
    }

    #[test]
    fn gas_constant_reexport_is_codata_value() {
        assert!((GAS_CONSTANT_J_PER_MOL_K - 8.314_462_618_153_24).abs() < 1e-12);
    }
}
