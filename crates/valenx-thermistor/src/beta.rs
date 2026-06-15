//! Single-parameter beta (B) thermistor model.
//!
//! ## Model
//!
//! The beta model is the simplest useful resistance-temperature law for
//! an NTC thermistor. Given a reference resistance `R0` measured at a
//! reference absolute temperature `T0`, the resistance at any other
//! absolute temperature `T` is
//!
//! ```text
//! R(T) = R0 * exp( beta * (1/T - 1/T0) )
//! ```
//!
//! where `beta` (also written `B`) is a material constant in kelvin,
//! typically a few thousand kelvin for commercial NTC parts. Solving
//! for temperature inverts the relation:
//!
//! ```text
//! 1/T = 1/T0 + (1/beta) * ln(R / R0)
//! ```
//!
//! ## Honest scope
//!
//! This is the two-point (single-`beta`) approximation. Real
//! thermistors deviate from it away from the calibration point; the
//! [`crate::steinhart`] three-coefficient model is more accurate over a
//! wide span. Temperatures are kelvin; resistances are ohms. No
//! self-heating, lead resistance, or tolerance modelling.

use crate::error::{check_resistance, check_temperature, ThermistorError};
use serde::{Deserialize, Serialize};

/// A calibrated single-parameter beta model for an NTC thermistor.
///
/// Construct with [`BetaModel::new`] (validated) and convert in either
/// direction with [`BetaModel::resistance_at`] and
/// [`BetaModel::temperature_at`]. All fields are kept public-read via
/// accessors; the struct is immutable once built.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BetaModel {
    /// Reference resistance `R0` at the reference temperature, in ohms.
    r0: f64,
    /// Reference absolute temperature `T0`, in kelvin.
    t0: f64,
    /// Material constant `beta` (`B`), in kelvin.
    beta: f64,
}

impl BetaModel {
    /// Build a beta model from a reference point and material constant.
    ///
    /// # Parameters
    ///
    /// - `r0_ohms`: reference resistance, strictly positive (ohms).
    /// - `t0_kelvin`: reference absolute temperature, strictly positive
    ///   (kelvin).
    /// - `beta_kelvin`: the material `beta` constant, strictly positive
    ///   and finite (kelvin). A positive `beta` gives NTC behaviour
    ///   (resistance falls as temperature rises).
    ///
    /// # Errors
    ///
    /// Returns [`ThermistorError::BadParameter`] if any argument is
    /// non-finite, or if `r0_ohms`, `t0_kelvin`, or `beta_kelvin` is
    /// not strictly positive.
    pub fn new(r0_ohms: f64, t0_kelvin: f64, beta_kelvin: f64) -> Result<Self, ThermistorError> {
        let r0 = check_resistance("r0_ohms", r0_ohms)?;
        let t0 = check_temperature("t0_kelvin", t0_kelvin)?;
        if !beta_kelvin.is_finite() {
            return Err(ThermistorError::BadParameter {
                name: "beta_kelvin",
                value: beta_kelvin,
                reason: "beta must be finite",
            });
        }
        if beta_kelvin <= 0.0 {
            return Err(ThermistorError::BadParameter {
                name: "beta_kelvin",
                value: beta_kelvin,
                reason: "beta must be strictly positive (kelvin)",
            });
        }
        Ok(BetaModel {
            r0,
            t0,
            beta: beta_kelvin,
        })
    }

    /// Calibrate a beta model from two measured resistance/temperature
    /// pairs.
    ///
    /// Solves the beta equation for `beta` given two points
    /// `(R1, T1)` and `(R2, T2)`:
    ///
    /// ```text
    /// beta = ln(R1 / R2) / (1/T1 - 1/T2)
    /// ```
    ///
    /// The first point `(R1, T1)` is retained as the reference
    /// `(R0, T0)`.
    ///
    /// # Errors
    ///
    /// Returns [`ThermistorError::BadParameter`] if any resistance or
    /// temperature is out of domain, [`ThermistorError::Degenerate`] if
    /// the two temperatures are equal (the denominator vanishes), and
    /// [`ThermistorError::NonFinite`] if the solved `beta` is not a
    /// finite positive number.
    pub fn calibrate_two_point(
        r1_ohms: f64,
        t1_kelvin: f64,
        r2_ohms: f64,
        t2_kelvin: f64,
    ) -> Result<Self, ThermistorError> {
        let r1 = check_resistance("r1_ohms", r1_ohms)?;
        let t1 = check_temperature("t1_kelvin", t1_kelvin)?;
        let r2 = check_resistance("r2_ohms", r2_ohms)?;
        let t2 = check_temperature("t2_kelvin", t2_kelvin)?;

        let inv_dt = 1.0 / t1 - 1.0 / t2;
        if inv_dt == 0.0 {
            return Err(ThermistorError::Degenerate(
                "the two calibration temperatures must differ",
            ));
        }
        let beta = (r1 / r2).ln() / inv_dt;
        if !beta.is_finite() || beta <= 0.0 {
            return Err(ThermistorError::NonFinite(
                "solved beta is not a finite positive number; check that R falls as T rises",
            ));
        }
        Ok(BetaModel {
            r0: r1,
            t0: t1,
            beta,
        })
    }

    /// Reference resistance `R0`, in ohms.
    pub fn r0_ohms(&self) -> f64 {
        self.r0
    }

    /// Reference absolute temperature `T0`, in kelvin.
    pub fn t0_kelvin(&self) -> f64 {
        self.t0
    }

    /// Material constant `beta` (`B`), in kelvin.
    pub fn beta_kelvin(&self) -> f64 {
        self.beta
    }

    /// Resistance predicted at absolute temperature `t_kelvin`, in ohms.
    ///
    /// Evaluates `R = R0 * exp(beta * (1/T - 1/T0))`. By construction
    /// `resistance_at(t0)` returns `r0` exactly (the exponent is zero).
    ///
    /// # Errors
    ///
    /// Returns [`ThermistorError::BadParameter`] if `t_kelvin` is out of
    /// domain, or [`ThermistorError::NonFinite`] if the exponential
    /// overflows to a non-finite value.
    pub fn resistance_at(&self, t_kelvin: f64) -> Result<f64, ThermistorError> {
        let t = check_temperature("t_kelvin", t_kelvin)?;
        let r = self.r0 * (self.beta * (1.0 / t - 1.0 / self.t0)).exp();
        if !r.is_finite() {
            return Err(ThermistorError::NonFinite(
                "beta-model resistance overflowed to a non-finite value",
            ));
        }
        Ok(r)
    }

    /// Absolute temperature predicted at resistance `r_ohms`, in kelvin.
    ///
    /// Inverts the beta law:
    /// `1/T = 1/T0 + (1/beta) * ln(R / R0)`. This is the exact inverse
    /// of [`resistance_at`](BetaModel::resistance_at).
    ///
    /// # Errors
    ///
    /// Returns [`ThermistorError::BadParameter`] if `r_ohms` is out of
    /// domain, or [`ThermistorError::NonFinite`] if the implied `1/T` is
    /// non-positive (which would correspond to a non-physical
    /// temperature).
    pub fn temperature_at(&self, r_ohms: f64) -> Result<f64, ThermistorError> {
        let r = check_resistance("r_ohms", r_ohms)?;
        let inv_t = 1.0 / self.t0 + (r / self.r0).ln() / self.beta;
        if !inv_t.is_finite() || inv_t <= 0.0 {
            return Err(ThermistorError::NonFinite(
                "beta-model inverse implies a non-physical temperature",
            ));
        }
        Ok(1.0 / inv_t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A representative 10 kohm-at-25C NTC with beta = 3950 K, the most
    /// common hobbyist/industrial part.
    fn ntc_10k() -> BetaModel {
        BetaModel::new(10_000.0, 298.15, 3950.0).expect("valid model")
    }

    #[test]
    fn resistance_at_reference_returns_r0_exactly() {
        let m = ntc_10k();
        let r = m.resistance_at(m.t0_kelvin()).unwrap();
        // The exponent is exactly 0 at T == T0, so R == R0 to within
        // rounding of `exp(0.0) == 1.0`.
        assert!((r - 10_000.0).abs() < 1e-6, "got {r}");
    }

    #[test]
    fn ntc_resistance_falls_as_temperature_rises() {
        let m = ntc_10k();
        let cold = m.resistance_at(273.15).unwrap(); // 0 C
        let mid = m.resistance_at(298.15).unwrap(); // 25 C
        let hot = m.resistance_at(323.15).unwrap(); // 50 C
        assert!(cold > mid, "cold {cold} should exceed mid {mid}");
        assert!(mid > hot, "mid {mid} should exceed hot {hot}");
    }

    #[test]
    fn temperature_inverts_resistance_round_trip() {
        let m = ntc_10k();
        for t in [273.15_f64, 280.0, 298.15, 310.15, 333.15, 350.0] {
            let r = m.resistance_at(t).unwrap();
            let back = m.temperature_at(r).unwrap();
            assert!(
                (back - t).abs() < 1e-9,
                "round trip failed at {t}: got {back}"
            );
        }
    }

    #[test]
    fn known_value_against_hand_computation() {
        // R(323.15) = 10000 * exp(3950 * (1/323.15 - 1/298.15)).
        // exponent = 3950 * (0.00309454... - 0.00335402...) = -1.024918...
        // R = 10000 * exp(-1.024918...) = 3588.182582 ohms.
        let m = ntc_10k();
        let r = m.resistance_at(323.15).unwrap();
        let expected: f64 = 10_000.0 * (3950.0_f64 * (1.0 / 323.15 - 1.0 / 298.15)).exp();
        assert!((r - expected).abs() < 1e-9);
        // Cross-check against the independently-computed closed-form
        // value (rounded to 6 decimals).
        assert!((r - 3588.182582).abs() < 1e-3, "got {r}");
    }

    #[test]
    fn two_point_calibration_recovers_beta() {
        // Generate two exact points from a known model, then recover.
        let truth = ntc_10k();
        let r1 = truth.resistance_at(298.15).unwrap();
        let r2 = truth.resistance_at(348.15).unwrap();
        let fitted = BetaModel::calibrate_two_point(r1, 298.15, r2, 348.15).unwrap();
        assert!(
            (fitted.beta_kelvin() - 3950.0).abs() < 1e-6,
            "recovered beta {}",
            fitted.beta_kelvin()
        );
        assert!((fitted.r0_ohms() - r1).abs() < 1e-9);
        assert!((fitted.t0_kelvin() - 298.15).abs() < 1e-12);
    }

    #[test]
    fn calibration_rejects_equal_temperatures() {
        let err = BetaModel::calibrate_two_point(10_000.0, 298.15, 5_000.0, 298.15).unwrap_err();
        assert_eq!(err.code(), "thermistor.degenerate");
    }

    #[test]
    fn constructor_rejects_bad_inputs() {
        assert!(BetaModel::new(-1.0, 298.15, 3950.0).is_err());
        assert!(BetaModel::new(10_000.0, 0.0, 3950.0).is_err());
        assert!(BetaModel::new(10_000.0, 298.15, 0.0).is_err());
        assert!(BetaModel::new(10_000.0, 298.15, f64::NAN).is_err());
    }
}
