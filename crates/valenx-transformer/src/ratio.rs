//! Turns ratio and the ideal voltage / current relations.
//!
//! For an ideal two-winding transformer with `Np` primary turns and
//! `Ns` secondary turns the turns ratio is
//!
//! ```text
//! a = Np / Ns = Vp / Vs = Is / Ip
//! ```
//!
//! Voltage transforms in direct proportion to the turns ratio and
//! current in inverse proportion, so the product `V * I` (the
//! apparent power) is preserved across the winding — see
//! [`crate::power`].
//!
//! ## Honest scope
//!
//! These are the ideal-transformer relations from a circuits textbook.
//! There is no magnetising current, no leakage flux, no winding
//! resistance, and no core loss in this module; those would perturb the
//! exact proportionalities below. Use this for teaching, back-of-the-
//! envelope sizing, and unit-test ground truth, not as a substitute for
//! a measured or finite-element transformer model.

use serde::{Deserialize, Serialize};

use crate::error::TransformerError;

/// An ideal transformer described by its winding turns ratio.
///
/// The single state is the turns ratio
/// `a = Np / Ns = Vp / Vs = Is / Ip`. Construct it either directly from
/// a ratio with [`TurnsRatio::new`] or from the integer (or fractional)
/// turn counts of the two windings with [`TurnsRatio::from_turns`].
///
/// A ratio `a > 1` is a *step-down* transformer (more primary turns, so
/// the secondary voltage is lower than the primary); a ratio `a < 1` is
/// a *step-up* transformer. See [`TurnsRatio::is_step_up`] and
/// [`TurnsRatio::is_step_down`].
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TurnsRatio {
    /// The dimensionless turns ratio `a = Np / Ns`. Always strictly
    /// positive by construction.
    a: f64,
}

impl TurnsRatio {
    /// Build a turns ratio directly from the dimensionless value
    /// `a = Np / Ns`.
    ///
    /// # Errors
    ///
    /// Returns [`TransformerError::Invalid`] unless `a` is finite and
    /// strictly positive.
    pub fn new(a: f64) -> Result<Self, TransformerError> {
        if !a.is_finite() || a <= 0.0 {
            return Err(TransformerError::invalid(
                "turns_ratio",
                format!("turns ratio a must be finite and positive, got {a}"),
            ));
        }
        Ok(Self { a })
    }

    /// Build a turns ratio from the two winding turn counts,
    /// `a = turns_primary / turns_secondary`.
    ///
    /// Turn counts are taken as `f64` so that fractional effective
    /// turns (for example from a tapped winding) are expressible; both
    /// must be finite and strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`TransformerError::Invalid`] if either count is not
    /// finite and strictly positive.
    pub fn from_turns(turns_primary: f64, turns_secondary: f64) -> Result<Self, TransformerError> {
        if !turns_primary.is_finite() || turns_primary <= 0.0 {
            return Err(TransformerError::invalid(
                "turns_primary",
                format!("number of primary turns must be finite and positive, got {turns_primary}"),
            ));
        }
        if !turns_secondary.is_finite() || turns_secondary <= 0.0 {
            return Err(TransformerError::invalid(
                "turns_secondary",
                format!(
                    "number of secondary turns must be finite and positive, got {turns_secondary}"
                ),
            ));
        }
        Self::new(turns_primary / turns_secondary)
    }

    /// The dimensionless turns ratio `a = Np / Ns = Vp / Vs = Is / Ip`.
    #[must_use]
    pub fn ratio(&self) -> f64 {
        self.a
    }

    /// The reciprocal ratio `Ns / Np = Vs / Vp = Ip / Is`.
    ///
    /// Always finite and positive because `a` is validated positive at
    /// construction.
    #[must_use]
    pub fn inverse(&self) -> f64 {
        1.0 / self.a
    }

    /// Secondary voltage from a primary voltage: `Vs = Vp / a`.
    ///
    /// Voltage scales inversely with the turns ratio, so a step-down
    /// ratio (`a > 1`) lowers the voltage.
    ///
    /// # Errors
    ///
    /// Returns [`TransformerError::Invalid`] if `voltage_primary` is not
    /// finite.
    pub fn secondary_voltage(&self, voltage_primary: f64) -> Result<f64, TransformerError> {
        if !voltage_primary.is_finite() {
            return Err(TransformerError::invalid(
                "voltage_primary",
                format!("primary voltage must be finite, got {voltage_primary}"),
            ));
        }
        Ok(voltage_primary / self.a)
    }

    /// Primary voltage from a secondary voltage: `Vp = a * Vs`.
    ///
    /// # Errors
    ///
    /// Returns [`TransformerError::Invalid`] if `voltage_secondary` is
    /// not finite.
    pub fn primary_voltage(&self, voltage_secondary: f64) -> Result<f64, TransformerError> {
        if !voltage_secondary.is_finite() {
            return Err(TransformerError::invalid(
                "voltage_secondary",
                format!("secondary voltage must be finite, got {voltage_secondary}"),
            ));
        }
        Ok(self.a * voltage_secondary)
    }

    /// Secondary current from a primary current: `Is = a * Ip`.
    ///
    /// Current scales in direct proportion to the turns ratio (the
    /// inverse of voltage), so a step-down ratio (`a > 1`) raises the
    /// current.
    ///
    /// # Errors
    ///
    /// Returns [`TransformerError::Invalid`] if `current_primary` is not
    /// finite.
    pub fn secondary_current(&self, current_primary: f64) -> Result<f64, TransformerError> {
        if !current_primary.is_finite() {
            return Err(TransformerError::invalid(
                "current_primary",
                format!("primary current must be finite, got {current_primary}"),
            ));
        }
        Ok(self.a * current_primary)
    }

    /// Primary current from a secondary current: `Ip = Is / a`.
    ///
    /// # Errors
    ///
    /// Returns [`TransformerError::Invalid`] if `current_secondary` is
    /// not finite.
    pub fn primary_current(&self, current_secondary: f64) -> Result<f64, TransformerError> {
        if !current_secondary.is_finite() {
            return Err(TransformerError::invalid(
                "current_secondary",
                format!("secondary current must be finite, got {current_secondary}"),
            ));
        }
        Ok(current_secondary / self.a)
    }

    /// `true` when this is a step-up transformer (`a < 1`): the
    /// secondary voltage exceeds the primary voltage and the secondary
    /// current is lower than the primary current.
    #[must_use]
    pub fn is_step_up(&self) -> bool {
        self.a < 1.0
    }

    /// `true` when this is a step-down transformer (`a > 1`): the
    /// secondary voltage is lower than the primary voltage and the
    /// secondary current exceeds the primary current.
    #[must_use]
    pub fn is_step_down(&self) -> bool {
        self.a > 1.0
    }

    /// `true` when this is an isolation transformer (`a == 1`): equal
    /// turns, so voltage and current pass through unchanged.
    ///
    /// The comparison is exact-equality against `1.0`; for ratios built
    /// from `from_turns` with equal counts this holds exactly because
    /// `x / x == 1.0` for any finite non-zero `x`.
    #[must_use]
    pub fn is_isolation(&self) -> bool {
        self.a == 1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference tolerance for the analytic float checks.
    const EPS: f64 = 1e-12;

    #[test]
    fn ratio_equals_turns_quotient() {
        // a = Np / Ns for 240 primary turns over 24 secondary turns.
        let t = TurnsRatio::from_turns(240.0, 24.0).unwrap();
        assert!((t.ratio() - 10.0).abs() < EPS, "got {}", t.ratio());
        assert!((t.inverse() - 0.1).abs() < EPS, "got {}", t.inverse());
    }

    #[test]
    fn voltage_ratio_equals_turns_ratio() {
        // Vp / Vs must equal Np / Ns. Step down 230 V by a = 23/2 = 11.5.
        let t = TurnsRatio::from_turns(230.0, 20.0).unwrap();
        let vp = 230.0;
        let vs = t.secondary_voltage(vp).unwrap();
        assert!((vp / vs - t.ratio()).abs() < EPS, "Vp/Vs got {}", vp / vs);
        // Inverse direction round-trips to the original primary voltage.
        let vp_back = t.primary_voltage(vs).unwrap();
        assert!((vp_back - vp).abs() < EPS, "round-trip Vp got {vp_back}");
    }

    #[test]
    fn current_ratio_is_inverse_turns_ratio() {
        // Ip / Is must equal Ns / Np = 1/a.
        let t = TurnsRatio::new(4.0).unwrap();
        let ip = 2.5;
        let is = t.secondary_current(ip).unwrap();
        assert!((ip / is - t.inverse()).abs() < EPS, "Ip/Is got {}", ip / is);
        // And Is / Ip == a.
        assert!((is / ip - t.ratio()).abs() < EPS, "Is/Ip got {}", is / ip);
        // Reverse direction round-trips.
        let ip_back = t.primary_current(is).unwrap();
        assert!((ip_back - ip).abs() < EPS, "round-trip Ip got {ip_back}");
    }

    #[test]
    fn step_up_raises_voltage_and_lowers_current() {
        // a = 0.5 is a step-up winding: Vs > Vp and Is < Ip.
        let t = TurnsRatio::new(0.5).unwrap();
        assert!(t.is_step_up());
        assert!(!t.is_step_down());
        let vp = 120.0;
        let ip = 4.0;
        let vs = t.secondary_voltage(vp).unwrap();
        let is = t.secondary_current(ip).unwrap();
        assert!(vs > vp, "expected step-up Vs ({vs}) > Vp ({vp})");
        assert!(is < ip, "expected step-up Is ({is}) < Ip ({ip})");
        // Exact analytic values: Vs = 240, Is = 2.
        assert!((vs - 240.0).abs() < EPS, "Vs got {vs}");
        assert!((is - 2.0).abs() < EPS, "Is got {is}");
    }

    #[test]
    fn step_down_lowers_voltage_and_raises_current() {
        let t = TurnsRatio::new(10.0).unwrap();
        assert!(t.is_step_down());
        assert!(!t.is_step_up());
        let vp = 230.0;
        let ip = 0.5;
        let vs = t.secondary_voltage(vp).unwrap();
        let is = t.secondary_current(ip).unwrap();
        assert!(vs < vp, "expected step-down Vs ({vs}) < Vp ({vp})");
        assert!(is > ip, "expected step-down Is ({is}) > Ip ({ip})");
        assert!((vs - 23.0).abs() < EPS, "Vs got {vs}");
        assert!((is - 5.0).abs() < EPS, "Is got {is}");
    }

    #[test]
    fn isolation_passes_through_unchanged() {
        let t = TurnsRatio::from_turns(100.0, 100.0).unwrap();
        assert!(t.is_isolation());
        assert!(!t.is_step_up());
        assert!(!t.is_step_down());
        let vp = 120.0;
        assert!(
            (t.secondary_voltage(vp).unwrap() - vp).abs() < EPS,
            "isolation should not change voltage"
        );
        let ip = 3.0;
        assert!(
            (t.secondary_current(ip).unwrap() - ip).abs() < EPS,
            "isolation should not change current"
        );
    }

    #[test]
    fn rejects_non_positive_and_non_finite_inputs() {
        assert!(TurnsRatio::new(0.0).is_err());
        assert!(TurnsRatio::new(-2.0).is_err());
        assert!(TurnsRatio::new(f64::NAN).is_err());
        assert!(TurnsRatio::new(f64::INFINITY).is_err());
        assert!(TurnsRatio::from_turns(0.0, 10.0).is_err());
        assert!(TurnsRatio::from_turns(10.0, 0.0).is_err());
        assert!(TurnsRatio::from_turns(-1.0, 10.0).is_err());

        let t = TurnsRatio::new(2.0).unwrap();
        assert!(t.secondary_voltage(f64::NAN).is_err());
        assert!(t.primary_current(f64::INFINITY).is_err());
    }
}
