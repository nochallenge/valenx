//! # valenx-thermistor
//!
//! Resistance-temperature modelling for NTC/PTC thermistors.
//!
//! ## What
//!
//! Two complementary thermistor models, each converting between
//! electrical resistance and absolute temperature in both directions:
//!
//! - [`BetaModel`] — the single-parameter beta (`B`) model
//!   `R = R0 * exp(beta * (1/T - 1/T0))`, the simplest two-point law,
//!   plus [`BetaModel::calibrate_two_point`] to solve `beta` from two
//!   measured `(R, T)` pairs.
//! - [`SteinhartHart`] — the three-coefficient model
//!   `1/T = A + B*ln(R) + C*ln(R)^3`, accurate over a wide span, with an
//!   analytic (Cardano) resistance-from-temperature inverse and
//!   [`SteinhartHart::fit_three_point`] to solve `(A, B, C)` from three
//!   measured `(R, T)` pairs.
//!
//! Both models also report the **temperature coefficient of resistance**
//! `alpha = (1/R) dR/dT` (in `1/K`) at a given temperature
//! ([`BetaModel::temperature_coefficient_at`],
//! [`SteinhartHart::temperature_coefficient_at`]) — the thermistor's
//! sensitivity figure of merit, conventionally quoted as a few percent
//! per kelvin (negative for an NTC).
//!
//! The [`units`] module converts between kelvin (used by every physics
//! function here) and degrees Celsius.
//!
//! ## Model
//!
//! Absolute temperature `T` is in **kelvin** and resistance `R` is in
//! **ohms** throughout, because both governing equations are written in
//! terms of `1/T`. The beta model is exact at its single calibration
//! point and approximate elsewhere; the Steinhart-Hart model is fit to
//! three points and is the standard high-accuracy curve. Both describe
//! NTC behaviour (resistance falls as temperature rises) when their
//! parameters are physical, and the resistance/temperature conversions
//! in each model are exact inverses of one another.
//!
//! The temperature coefficient of resistance is the logarithmic slope of
//! each curve, `alpha = (1/R) dR/dT`: it is `-beta / T^2` for the beta
//! model and `-1 / (T^2 (B + 3 C ln(R)^2))` for Steinhart-Hart. Both are
//! negative for a physical NTC and are validated against a
//! central-difference of the resistance curve.
//!
//! ```
//! use valenx_thermistor::{BetaModel, units::celsius_to_kelvin};
//!
//! // 10 kohm-at-25C NTC, beta = 3950 K.
//! let ntc = BetaModel::new(10_000.0, celsius_to_kelvin(25.0), 3950.0).unwrap();
//! let r_at_50c = ntc.resistance_at(celsius_to_kelvin(50.0)).unwrap();
//! assert!(r_at_50c < 10_000.0); // NTC: hotter means lower resistance
//! let t = ntc.temperature_at(r_at_50c).unwrap();
//! assert!((t - celsius_to_kelvin(50.0)).abs() < 1e-9);
//! ```
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form and
//! numerical models — the beta and Steinhart-Hart equations and their
//! least-squares/exact calibrations — and nothing more. The crate does
//! **not** model self-heating, lead/wire resistance, the
//! voltage-divider readout circuit, dissipation constants, part
//! tolerance, or aging, and it performs no datasheet lookup. It is
//! **not** a clinical, medical, or production engineering tool and is
//! not a substitute for calibrated instrumentation or a manufacturer
//! datasheet. Validate against real hardware before relying on any
//! number it produces.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, ThermistorError>`](error::ThermistorError). The error
//! type exposes a stable [`code`](error::ThermistorError::code) and a
//! coarse [`category`](error::ThermistorError::category) for
//! programmatic handling.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod beta;
pub mod error;
pub mod steinhart;
pub mod units;

pub use beta::BetaModel;
pub use error::{ErrorCategory, ThermistorError};
pub use steinhart::SteinhartHart;

#[cfg(test)]
mod tests {
    use super::*;

    /// The two models, when the Steinhart-Hart coefficients are fit to
    /// the *same* curve the beta model describes near its calibration
    /// point, should agree on resistance close to that point.
    ///
    /// Construction: take a beta model, sample it at three temperatures
    /// to fit a Steinhart-Hart model, then compare predictions at a
    /// temperature near the middle calibration point. The two laws are
    /// different functional forms, so they diverge far from the fit
    /// span, but near it the agreement is tight.
    #[test]
    fn beta_and_steinhart_agree_near_calibration() {
        let beta = BetaModel::new(10_000.0, 298.15, 3950.0).unwrap();

        // Fit Steinhart-Hart to the beta curve at 10, 25, 40 C.
        let temps = [283.15_f64, 298.15, 313.15];
        let mut pts = [(0.0, 0.0); 3];
        for (i, t) in temps.iter().enumerate() {
            pts[i] = (beta.resistance_at(*t).unwrap(), *t);
        }
        let sh = SteinhartHart::fit_three_point(pts).unwrap();

        // At the central calibration point both must reproduce R0.
        let r_beta = beta.resistance_at(298.15).unwrap();
        let r_sh = sh.resistance_at(298.15).unwrap();
        assert!(
            (r_beta - 10_000.0).abs() < 1e-6,
            "beta should give R0 at T0, got {r_beta}"
        );
        assert!(
            (r_beta - r_sh).abs() / r_beta < 1e-9,
            "models disagree at calibration centre: beta={r_beta}, sh={r_sh}"
        );

        // Near the fit span (30 C) they should still be very close,
        // well under 1 ohm on a ~7 kohm value.
        let t_near = 303.15;
        let r_beta_near = beta.resistance_at(t_near).unwrap();
        let r_sh_near = sh.resistance_at(t_near).unwrap();
        assert!(
            (r_beta_near - r_sh_near).abs() < 1.0,
            "models diverge near calibration: beta={r_beta_near}, sh={r_sh_near}"
        );
    }

    /// Cross-model temperature agreement: feed the same resistance to
    /// both models (fit to the same curve) and confirm the recovered
    /// temperature matches near the calibration span.
    #[test]
    fn beta_and_steinhart_agree_on_temperature_near_calibration() {
        let beta = BetaModel::new(10_000.0, 298.15, 3950.0).unwrap();
        let temps = [283.15_f64, 298.15, 313.15];
        let mut pts = [(0.0, 0.0); 3];
        for (i, t) in temps.iter().enumerate() {
            pts[i] = (beta.resistance_at(*t).unwrap(), *t);
        }
        let sh = SteinhartHart::fit_three_point(pts).unwrap();

        let r = 9_000.0; // close to the 10k centre
        let t_beta = beta.temperature_at(r).unwrap();
        let t_sh = sh.temperature_at(r).unwrap();
        assert!(
            (t_beta - t_sh).abs() < 0.05,
            "temperatures disagree: beta={t_beta} K, sh={t_sh} K"
        );
    }

    #[test]
    fn public_reexports_are_usable() {
        // Smoke test the crate-root re-exports compile and resolve.
        let _: ThermistorError = ThermistorError::Degenerate("x");
        let _: ErrorCategory = ErrorCategory::Input;
        let m = BetaModel::new(1.0, 300.0, 3000.0).unwrap();
        assert!(m.resistance_at(300.0).is_ok());
        let sh = SteinhartHart::new(1e-3, 2e-4, 2e-7).unwrap();
        assert!(sh.temperature_at(10_000.0).is_ok());
    }
}
