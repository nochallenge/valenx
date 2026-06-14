//! Power conservation and efficiency.
//!
//! ## Ideal transformer
//!
//! An ideal transformer is lossless, so all input power appears at the
//! output:
//!
//! ```text
//! Pin = Vp * Ip = Vs * Is = Pout
//! ```
//!
//! This is the direct consequence of the [`crate::ratio`] relations:
//! voltage multiplies by `a` while current divides by `a` (or vice
//! versa), so the product is invariant.
//!
//! ## Real transformer
//!
//! A real transformer dissipates some power in its core (hysteresis and
//! eddy currents) and windings (`I^2 R`), so the output is smaller than
//! the input. The efficiency is the dimensionless ratio
//!
//! ```text
//! eta = Pout / Pin,   with 0 < eta <= 1
//! ```
//!
//! and the dissipated power is `Ploss = Pin - Pout = (1 - eta) * Pin`.
//!
//! ## Honest scope
//!
//! Efficiency here is a single user-supplied number, not a quantity
//! derived from a core-loss curve, a load profile, or a temperature.
//! The relations are exact algebra given that number; they do not model
//! where the loss comes from or how it varies with load. This is a
//! teaching / sizing aid, not a loss-prediction model.

use serde::{Deserialize, Serialize};

use crate::error::TransformerError;

/// Compute apparent power `P = V * I` from a voltage and current.
///
/// Both arguments must be finite. The result carries the sign of the
/// product, so callers passing signed phasor magnitudes get a signed
/// power back; pass magnitudes for an unsigned apparent power.
///
/// # Errors
///
/// Returns [`TransformerError::Invalid`] if either argument is not
/// finite.
pub fn apparent_power(voltage: f64, current: f64) -> Result<f64, TransformerError> {
    if !voltage.is_finite() {
        return Err(TransformerError::invalid(
            "voltage",
            format!("voltage must be finite, got {voltage}"),
        ));
    }
    if !current.is_finite() {
        return Err(TransformerError::invalid(
            "current",
            format!("current must be finite, got {current}"),
        ));
    }
    Ok(voltage * current)
}

/// Output power of an ideal (lossless) transformer given its input
/// power: `Pout = Pin`.
///
/// Provided as a named function so call sites read as the physical
/// statement "ideal output equals input" rather than an identity.
///
/// # Errors
///
/// Returns [`TransformerError::Invalid`] if `power_in` is not finite.
pub fn ideal_output_power(power_in: f64) -> Result<f64, TransformerError> {
    if !power_in.is_finite() {
        return Err(TransformerError::invalid(
            "power_in",
            format!("input power must be finite, got {power_in}"),
        ));
    }
    Ok(power_in)
}

/// An efficiency `eta = Pout / Pin`, validated to lie in `(0, 1]`.
///
/// A standalone newtype so an efficiency can be passed around, stored,
/// and serialised with the domain check guaranteed once at
/// construction. An ideal transformer is [`Efficiency::ideal`]
/// (`eta = 1`).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Efficiency {
    /// The dimensionless efficiency, always in `(0, 1]` by construction.
    eta: f64,
}

impl Efficiency {
    /// Build an efficiency from a dimensionless ratio in `(0, 1]`.
    ///
    /// # Errors
    ///
    /// Returns [`TransformerError::Invalid`] unless `eta` is finite and
    /// in the half-open interval `(0, 1]`. Zero is rejected (a
    /// transformer that delivered no output would not be operating) and
    /// values above one are rejected (a passive transformer cannot
    /// create power).
    pub fn new(eta: f64) -> Result<Self, TransformerError> {
        if !eta.is_finite() || eta <= 0.0 || eta > 1.0 {
            return Err(TransformerError::invalid(
                "efficiency",
                format!("efficiency must be finite and in (0, 1], got {eta}"),
            ));
        }
        Ok(Self { eta })
    }

    /// The efficiency of an ideal lossless transformer, `eta = 1`.
    #[must_use]
    pub fn ideal() -> Self {
        Self { eta: 1.0 }
    }

    /// Derive an efficiency from a measured output and input power,
    /// `eta = Pout / Pin`.
    ///
    /// # Errors
    ///
    /// Returns [`TransformerError::Invalid`] if either power is not
    /// finite, if `power_in` is not strictly positive, or if the
    /// resulting ratio falls outside `(0, 1]` (for example because the
    /// output exceeds the input, which is unphysical for a passive
    /// transformer).
    pub fn from_powers(power_out: f64, power_in: f64) -> Result<Self, TransformerError> {
        if !power_out.is_finite() {
            return Err(TransformerError::invalid(
                "power_out",
                format!("output power must be finite, got {power_out}"),
            ));
        }
        if !power_in.is_finite() || power_in <= 0.0 {
            return Err(TransformerError::invalid(
                "power_in",
                format!("input power must be finite and positive, got {power_in}"),
            ));
        }
        Self::new(power_out / power_in)
    }

    /// The dimensionless efficiency value `eta` in `(0, 1]`.
    #[must_use]
    pub fn value(&self) -> f64 {
        self.eta
    }

    /// Output power for a given input power: `Pout = eta * Pin`.
    ///
    /// # Errors
    ///
    /// Returns [`TransformerError::Invalid`] if `power_in` is not finite.
    pub fn output_power(&self, power_in: f64) -> Result<f64, TransformerError> {
        if !power_in.is_finite() {
            return Err(TransformerError::invalid(
                "power_in",
                format!("input power must be finite, got {power_in}"),
            ));
        }
        Ok(self.eta * power_in)
    }

    /// Input power required to deliver a given output power:
    /// `Pin = Pout / eta`.
    ///
    /// # Errors
    ///
    /// Returns [`TransformerError::Invalid`] if `power_out` is not
    /// finite. The denominator `eta` is guaranteed non-zero by
    /// construction.
    pub fn input_power(&self, power_out: f64) -> Result<f64, TransformerError> {
        if !power_out.is_finite() {
            return Err(TransformerError::invalid(
                "power_out",
                format!("output power must be finite, got {power_out}"),
            ));
        }
        Ok(power_out / self.eta)
    }

    /// Power dissipated as loss for a given input power:
    /// `Ploss = Pin - Pout = (1 - eta) * Pin`.
    ///
    /// Zero for an ideal transformer.
    ///
    /// # Errors
    ///
    /// Returns [`TransformerError::Invalid`] if `power_in` is not finite.
    pub fn power_loss(&self, power_in: f64) -> Result<f64, TransformerError> {
        if !power_in.is_finite() {
            return Err(TransformerError::invalid(
                "power_in",
                format!("input power must be finite, got {power_in}"),
            ));
        }
        Ok((1.0 - self.eta) * power_in)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference tolerance for the analytic float checks.
    const EPS: f64 = 1e-12;

    #[test]
    fn ideal_power_is_conserved() {
        // Vp Ip = Vs Is for a step-down winding: 240 V, 2 A in.
        let vp = 240.0;
        let ip = 2.0;
        let pin = apparent_power(vp, ip).unwrap();
        assert!((pin - 480.0).abs() < EPS, "Pin got {pin}");

        // Step down by a = 10: Vs = 24, Is = 20, Pout = 480 = Pin.
        let vs = 24.0;
        let is = 20.0;
        let pout = apparent_power(vs, is).unwrap();
        assert!(
            (pout - pin).abs() < EPS,
            "ideal transformer must conserve power: Pin {pin}, Pout {pout}"
        );
        assert!((ideal_output_power(pin).unwrap() - pin).abs() < EPS);
    }

    #[test]
    fn efficiency_within_unit_interval() {
        assert!(Efficiency::new(0.97).is_ok());
        assert!(Efficiency::new(1.0).is_ok());
        // Boundary: zero and above-one are rejected; negatives too.
        assert!(Efficiency::new(0.0).is_err());
        assert!(Efficiency::new(-0.1).is_err());
        assert!(Efficiency::new(1.0001).is_err());
        assert!(Efficiency::new(f64::NAN).is_err());

        let eta = Efficiency::new(0.9).unwrap();
        assert!(eta.value() > 0.0 && eta.value() <= 1.0);
    }

    #[test]
    fn efficiency_from_powers_round_trips() {
        // 950 W out of 1000 W in is 95% efficient.
        let eta = Efficiency::from_powers(950.0, 1000.0).unwrap();
        assert!((eta.value() - 0.95).abs() < EPS, "eta got {}", eta.value());
        // Pout = eta * Pin recovers the output.
        assert!(
            (eta.output_power(1000.0).unwrap() - 950.0).abs() < EPS,
            "Pout mismatch"
        );
        // Pin = Pout / eta recovers the input.
        assert!(
            (eta.input_power(950.0).unwrap() - 1000.0).abs() < EPS,
            "Pin mismatch"
        );
        // Loss is the complement.
        assert!(
            (eta.power_loss(1000.0).unwrap() - 50.0).abs() < EPS,
            "Ploss got {}",
            eta.power_loss(1000.0).unwrap()
        );
    }

    #[test]
    fn ideal_efficiency_has_no_loss() {
        let eta = Efficiency::ideal();
        assert!((eta.value() - 1.0).abs() < EPS);
        assert!((eta.output_power(750.0).unwrap() - 750.0).abs() < EPS);
        assert!(
            eta.power_loss(750.0).unwrap().abs() < EPS,
            "ideal transformer must have zero loss"
        );
    }

    #[test]
    fn from_powers_rejects_super_unity_and_bad_input() {
        // Output exceeding input is unphysical for a passive device.
        assert!(Efficiency::from_powers(1100.0, 1000.0).is_err());
        // Non-positive input power.
        assert!(Efficiency::from_powers(500.0, 0.0).is_err());
        assert!(Efficiency::from_powers(500.0, -1.0).is_err());
        // Non-finite arguments.
        assert!(Efficiency::from_powers(f64::NAN, 1000.0).is_err());
        assert!(apparent_power(f64::INFINITY, 1.0).is_err());
    }
}
