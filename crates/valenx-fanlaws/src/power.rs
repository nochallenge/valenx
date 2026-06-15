//! Fan power identity and the validated [`Efficiency`] type.
//!
//! The *air power* (useful flow work delivered to the gas per unit
//! time) of a fan is the product of volumetric flow and pressure rise:
//!
//! ```text
//! P_air = Q * dP
//! ```
//!
//! The *shaft power* the motor must supply is larger by the reciprocal
//! of the fan's total efficiency `eta`, which lumps together
//! aerodynamic, volumetric, and mechanical losses:
//!
//! ```text
//! P_shaft = Q * dP / eta
//! ```
//!
//! In consistent SI units (`Q` in m^3/s, `dP` in pascals) the result is
//! in watts. Because `eta` lies in `(0, 1]`, shaft power is always at
//! least the air power and the division never blows up.
//!
//! ## Honest scope
//!
//! `eta` here is a single lumped total-efficiency figure assumed
//! constant; real fan efficiency varies strongly with the operating
//! point along the fan curve. This identity is the textbook definition,
//! not a performance prediction.

use crate::error::{require_non_negative, FanLawError};

/// A fan total efficiency, validated to lie in the physically
/// admissible half-open interval `(0, 1]`.
///
/// The lower bound is open because a zero-efficiency fan would require
/// infinite shaft power for any flow work; the upper bound is closed
/// because unity efficiency (no losses) is the ideal limit.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Efficiency(f64);

impl Efficiency {
    /// Construct a validated efficiency from a fraction in `(0, 1]`.
    ///
    /// # Errors
    ///
    /// Returns [`FanLawError::EfficiencyOutOfRange`] if `value` is not
    /// in `(0, 1]`, or [`FanLawError::NotFinite`] if it is NaN /
    /// infinite.
    pub fn new(value: f64) -> Result<Self, FanLawError> {
        if !value.is_finite() {
            return Err(FanLawError::NotFinite {
                name: "efficiency",
                value,
            });
        }
        if value > 0.0 && value <= 1.0 {
            Ok(Self(value))
        } else {
            Err(FanLawError::EfficiencyOutOfRange { value })
        }
    }

    /// Construct a validated efficiency from a percentage in `(0, 100]`.
    ///
    /// # Errors
    ///
    /// Returns [`FanLawError`] under the same conditions as
    /// [`Efficiency::new`], after dividing by 100.
    pub fn from_percent(percent: f64) -> Result<Self, FanLawError> {
        Self::new(percent / 100.0)
    }

    /// The underlying fraction in `(0, 1]`.
    pub fn fraction(self) -> f64 {
        self.0
    }

    /// The efficiency expressed as a percentage in `(0, 100]`.
    pub fn percent(self) -> f64 {
        self.0 * 100.0
    }
}

/// Air power (useful flow work) `P_air = Q * dP`.
///
/// In SI units (`flow` in m^3/s, `pressure` in Pa) the result is in
/// watts. This is the numerator of the shaft-power identity and is
/// independent of efficiency.
///
/// # Errors
///
/// Returns [`FanLawError`] if `flow` or `pressure` is negative or
/// non-finite.
pub fn air_power(flow: f64, pressure: f64) -> Result<f64, FanLawError> {
    let flow = require_non_negative("flow", flow)?;
    let pressure = require_non_negative("pressure", pressure)?;
    Ok(flow * pressure)
}

/// Shaft power `P_shaft = Q * dP / eta`.
///
/// The efficiency is supplied as a validated [`Efficiency`], so the
/// division is guaranteed finite and the result is never less than the
/// corresponding [`air_power`].
///
/// # Errors
///
/// Returns [`FanLawError`] if `flow` or `pressure` is negative or
/// non-finite.
pub fn shaft_power(flow: f64, pressure: f64, efficiency: Efficiency) -> Result<f64, FanLawError> {
    let p_air = air_power(flow, pressure)?;
    Ok(p_air / efficiency.fraction())
}

/// Back out the total efficiency implied by a measured shaft power for
/// a known flow and pressure rise: `eta = (Q * dP) / P_shaft`.
///
/// This is the inverse of [`shaft_power`]. The result is validated
/// through [`Efficiency::new`], so a measurement that implies an
/// efficiency outside `(0, 1]` (e.g. a shaft power smaller than the air
/// power, which is physically impossible) is reported as an error
/// rather than returned.
///
/// # Errors
///
/// Returns [`FanLawError::NonPositive`] if `shaft` is not strictly
/// positive, [`FanLawError`] if `flow` / `pressure` are negative or
/// non-finite, or [`FanLawError::EfficiencyOutOfRange`] if the implied
/// efficiency falls outside `(0, 1]`.
pub fn implied_efficiency(flow: f64, pressure: f64, shaft: f64) -> Result<Efficiency, FanLawError> {
    use crate::error::require_positive;
    let p_air = air_power(flow, pressure)?;
    let shaft = require_positive("shaft", shaft)?;
    Efficiency::new(p_air / shaft)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for analytic float comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn efficiency_accepts_open_lower_closed_upper_bound() {
        assert!(Efficiency::new(1.0).is_ok()); // closed upper
        assert!(Efficiency::new(0.5).is_ok());
        assert!(Efficiency::new(1e-6).is_ok());
        assert!(Efficiency::new(0.0).is_err()); // open lower
        assert!(Efficiency::new(-0.1).is_err());
        assert!(Efficiency::new(1.0001).is_err()); // above unity
        assert!(Efficiency::new(f64::NAN).is_err());
    }

    #[test]
    fn efficiency_percent_round_trips() {
        let e = Efficiency::from_percent(75.0).unwrap();
        assert!((e.fraction() - 0.75).abs() < EPS);
        assert!((e.percent() - 75.0).abs() < EPS);
        assert!(Efficiency::from_percent(0.0).is_err());
        assert!(Efficiency::from_percent(150.0).is_err());
    }

    #[test]
    fn air_power_is_flow_times_pressure() {
        // 2 m^3/s through a 500 Pa rise = 1000 W of flow work.
        assert!((air_power(2.0, 500.0).unwrap() - 1000.0).abs() < EPS);
        // Zero flow or zero pressure -> zero air power.
        assert!((air_power(0.0, 500.0).unwrap()).abs() < EPS);
        assert!((air_power(2.0, 0.0).unwrap()).abs() < EPS);
    }

    #[test]
    fn shaft_power_equals_air_power_over_efficiency() {
        // 1000 W air power at 50% efficiency -> 2000 W shaft.
        let eta = Efficiency::new(0.5).unwrap();
        assert!((shaft_power(2.0, 500.0, eta).unwrap() - 2000.0).abs() < EPS);

        // Ground-truth worked example: Q=3, dP=400 => air=1200 W,
        // at eta=0.6 => shaft = 2000 W.
        let eta = Efficiency::new(0.6).unwrap();
        assert!((shaft_power(3.0, 400.0, eta).unwrap() - 2000.0).abs() < 1e-7);
    }

    #[test]
    fn unity_efficiency_makes_shaft_equal_air_power() {
        let eta = Efficiency::new(1.0).unwrap();
        let air = air_power(4.0, 250.0).unwrap();
        let shaft = shaft_power(4.0, 250.0, eta).unwrap();
        assert!((air - shaft).abs() < EPS);
    }

    #[test]
    fn shaft_power_never_below_air_power() {
        // For any admissible efficiency, shaft >= air.
        let q = 5.0;
        let dp = 300.0;
        let air = air_power(q, dp).unwrap();
        for frac in [0.1_f64, 0.35, 0.6, 0.85, 1.0] {
            let eta = Efficiency::new(frac).unwrap();
            let shaft = shaft_power(q, dp, eta).unwrap();
            assert!(shaft >= air - EPS, "frac={frac}: shaft {shaft} < air {air}");
        }
    }

    #[test]
    fn implied_efficiency_inverts_shaft_power() {
        // Round trip: pick eta, compute shaft, recover eta.
        let q = 2.5;
        let dp = 480.0;
        let eta0 = Efficiency::new(0.62).unwrap();
        let shaft = shaft_power(q, dp, eta0).unwrap();
        let eta1 = implied_efficiency(q, dp, shaft).unwrap();
        assert!((eta1.fraction() - eta0.fraction()).abs() < 1e-9);
    }

    #[test]
    fn implied_efficiency_rejects_impossible_shaft_power() {
        // Shaft power below the air power implies eta > 1 -> error.
        let air = air_power(2.0, 500.0).unwrap(); // 1000 W
        assert!(implied_efficiency(2.0, 500.0, air / 2.0).is_err());
        // Non-positive shaft power -> error.
        assert!(implied_efficiency(2.0, 500.0, 0.0).is_err());
    }

    #[test]
    fn power_functions_reject_negative_inputs() {
        assert!(air_power(-1.0, 500.0).is_err());
        assert!(air_power(1.0, -500.0).is_err());
        let eta = Efficiency::new(0.5).unwrap();
        assert!(shaft_power(-1.0, 500.0, eta).is_err());
    }
}
