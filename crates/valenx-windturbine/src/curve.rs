//! Idealised wind-turbine power curve with cut-in / rated / cut-out
//! regions.
//!
//! # Model
//!
//! A real turbine does not produce `1/2 rho A v^3 Cp` at every wind
//! speed. Its controller enforces three break-points:
//!
//! - **cut-in** `v_in` — below it the wind cannot overcome friction, so
//!   output is **zero**.
//! - **rated** `v_r` — the speed at which the generator first reaches
//!   its rated electrical power `P_rated`.
//! - **cut-out** `v_out` — above it the machine feathers / brakes to
//!   protect itself, so output drops back to **zero**.
//!
//! This module models the textbook idealisation:
//!
//! ```text
//!            { 0                                      v <  v_in
//!            { P_rated * (v^3 - v_in^3)               v_in <= v < v_r
//!   P(v)  =  {           ---------------
//!            {           (v_r^3 - v_in^3)
//!            { P_rated                                v_r  <= v <= v_out
//!            { 0                                      v >  v_out
//! ```
//!
//! The region between cut-in and rated follows the cube law (the
//! captured power rises with `v^3`), normalised so that `P(v_in) = 0`
//! and `P(v_r) = P_rated`. Between rated and cut-out the output is held
//! **constant** at `P_rated`. Outside `[v_in, v_out]` it is zero.

use crate::error::WindTurbineError;
use serde::{Deserialize, Serialize};

/// Which operating region a wind speed falls into on the power curve.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum Region {
    /// Below cut-in (or above cut-out): the turbine produces no power.
    Idle,
    /// Between cut-in and rated: power rises with the cube of wind
    /// speed.
    Ramp,
    /// Between rated and cut-out: power is held constant at the rated
    /// value.
    Rated,
}

/// A validated idealised power curve.
///
/// Construct with [`PowerCurve::new`], which enforces
/// `0 < cut_in < rated < cut_out` and `rated_power > 0`. Evaluate with
/// [`PowerCurve::power`] (watts at a wind speed) or
/// [`PowerCurve::region`] (which region a speed is in).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PowerCurve {
    cut_in: f64,
    rated: f64,
    cut_out: f64,
    rated_power: f64,
}

impl PowerCurve {
    /// Build a power curve from its four parameters (all SI: wind speeds
    /// in m/s, `rated_power` in W).
    ///
    /// # Errors
    ///
    /// - [`WindTurbineError::BadParameter`] if `rated_power` is not
    ///   strictly positive, or if any speed is negative / non-finite.
    /// - [`WindTurbineError::InconsistentCurve`] unless
    ///   `0 < cut_in < rated < cut_out`.
    ///
    /// ```
    /// use valenx_windturbine::curve::PowerCurve;
    /// let c = PowerCurve::new(3.0, 12.0, 25.0, 2.0e6).unwrap();
    /// assert!((c.rated_power() - 2.0e6).abs() < 1e-6);
    /// ```
    pub fn new(
        cut_in: f64,
        rated: f64,
        cut_out: f64,
        rated_power: f64,
    ) -> Result<Self, WindTurbineError> {
        for (name, v) in [("cut_in", cut_in), ("rated", rated), ("cut_out", cut_out)] {
            if !v.is_finite() || v < 0.0 {
                return Err(WindTurbineError::BadParameter {
                    name,
                    reason: "must be a finite, non-negative wind speed".to_string(),
                });
            }
        }
        if !rated_power.is_finite() || rated_power <= 0.0 {
            return Err(WindTurbineError::BadParameter {
                name: "rated_power",
                reason: "must be > 0".to_string(),
            });
        }
        // All three speeds are finite (validated above), so these direct
        // comparisons are exact — no `partial_cmp` ambiguity.
        if cut_in <= 0.0 {
            return Err(WindTurbineError::InconsistentCurve(
                "cut_in must be > 0".to_string(),
            ));
        }
        if cut_in >= rated {
            return Err(WindTurbineError::InconsistentCurve(format!(
                "cut_in ({cut_in}) must be < rated ({rated})"
            )));
        }
        if rated >= cut_out {
            return Err(WindTurbineError::InconsistentCurve(format!(
                "rated ({rated}) must be < cut_out ({cut_out})"
            )));
        }
        Ok(Self {
            cut_in,
            rated,
            cut_out,
            rated_power,
        })
    }

    /// Cut-in wind speed (m/s): below this, output is zero.
    pub fn cut_in(&self) -> f64 {
        self.cut_in
    }

    /// Rated wind speed (m/s): output first reaches the rated power.
    pub fn rated(&self) -> f64 {
        self.rated
    }

    /// Cut-out wind speed (m/s): above this, the turbine shuts down.
    pub fn cut_out(&self) -> f64 {
        self.cut_out
    }

    /// Rated electrical power (W).
    pub fn rated_power(&self) -> f64 {
        self.rated_power
    }

    /// The [`Region`] a wind speed `v` (m/s) falls into.
    ///
    /// Speeds below cut-in or above cut-out are [`Region::Idle`]; the
    /// boundaries `v_in`, `v_r`, `v_out` themselves are *in-band*
    /// (`v_in` is the first [`Region::Ramp`] point, `v_r` the first
    /// [`Region::Rated`] point, `v_out` the last [`Region::Rated`]
    /// point).
    pub fn region(&self, v: f64) -> Region {
        if !v.is_finite() || v < self.cut_in || v > self.cut_out {
            Region::Idle
        } else if v < self.rated {
            Region::Ramp
        } else {
            Region::Rated
        }
    }

    /// Output power `P(v)` (W) at wind speed `v` (m/s), per the
    /// piecewise model documented at the module level.
    ///
    /// Returns `0.0` for non-finite `v` so callers can feed raw sensor
    /// data without a panic; finite inputs follow the curve exactly.
    ///
    /// ```
    /// use valenx_windturbine::curve::PowerCurve;
    /// let c = PowerCurve::new(3.0, 12.0, 25.0, 2.0e6).unwrap();
    /// assert_eq!(c.power(2.0), 0.0); // below cut-in
    /// assert_eq!(c.power(30.0), 0.0); // above cut-out
    /// // Plateau between rated and cut-out:
    /// assert!((c.power(20.0) - 2.0e6).abs() < 1e-6);
    /// ```
    pub fn power(&self, v: f64) -> f64 {
        match self.region(v) {
            Region::Idle => 0.0,
            Region::Rated => self.rated_power,
            Region::Ramp => {
                // Cube-law interpolation pinned to P(v_in)=0, P(v_r)=P_rated.
                let num = v * v * v - self.cut_in * self.cut_in * self.cut_in;
                let den =
                    self.rated * self.rated * self.rated - self.cut_in * self.cut_in * self.cut_in;
                self.rated_power * num / den
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    fn curve() -> PowerCurve {
        // Representative ~2 MW machine: cut-in 3, rated 12, cut-out 25 m/s.
        PowerCurve::new(3.0, 12.0, 25.0, 2.0e6).unwrap()
    }

    #[test]
    fn rejects_out_of_order_speeds() {
        // cut_in >= rated.
        assert!(PowerCurve::new(12.0, 12.0, 25.0, 1.0).is_err());
        assert!(PowerCurve::new(13.0, 12.0, 25.0, 1.0).is_err());
        // rated >= cut_out.
        assert!(PowerCurve::new(3.0, 25.0, 25.0, 1.0).is_err());
        assert!(PowerCurve::new(3.0, 26.0, 25.0, 1.0).is_err());
        // cut_in == 0.
        assert!(PowerCurve::new(0.0, 12.0, 25.0, 1.0).is_err());
    }

    #[test]
    fn ordering_error_has_curve_code() {
        let e = PowerCurve::new(13.0, 12.0, 25.0, 1.0).unwrap_err();
        assert_eq!(e.code(), "windturbine.inconsistent_curve");
    }

    #[test]
    fn rejects_bad_rated_power() {
        assert!(PowerCurve::new(3.0, 12.0, 25.0, 0.0).is_err());
        assert!(PowerCurve::new(3.0, 12.0, 25.0, -5.0).is_err());
        let e = PowerCurve::new(3.0, 12.0, 25.0, -5.0).unwrap_err();
        assert_eq!(e.code(), "windturbine.bad_parameter");
    }

    #[test]
    fn zero_below_cut_in() {
        let c = curve();
        assert!(c.power(0.0).abs() < EPS);
        assert!(c.power(1.0).abs() < EPS);
        assert!(c.power(2.999).abs() < EPS);
        assert_eq!(c.region(2.0), Region::Idle);
    }

    #[test]
    fn zero_above_cut_out() {
        let c = curve();
        assert!(c.power(25.001).abs() < EPS);
        assert!(c.power(30.0).abs() < EPS);
        assert!(c.power(100.0).abs() < EPS);
        assert_eq!(c.region(40.0), Region::Idle);
    }

    #[test]
    fn zero_at_exactly_cut_in() {
        // The cube-law numerator is zero at v = v_in.
        let c = curve();
        assert!(c.power(3.0).abs() < EPS);
        assert_eq!(c.region(3.0), Region::Ramp);
    }

    #[test]
    fn rated_power_at_exactly_rated_speed() {
        let c = curve();
        assert!((c.power(12.0) - 2.0e6).abs() < EPS);
        assert_eq!(c.region(12.0), Region::Rated);
    }

    #[test]
    fn constant_rated_between_rated_and_cut_out() {
        let c = curve();
        // Output is held flat at P_rated across the whole plateau.
        for &v in &[12.0, 15.0, 18.0, 21.0, 24.0, 25.0] {
            assert!(
                (c.power(v) - 2.0e6).abs() < EPS,
                "expected flat rated power at v = {v}, got {}",
                c.power(v)
            );
            assert_eq!(c.region(v), Region::Rated);
        }
    }

    #[test]
    fn ramp_is_strictly_increasing_and_bounded() {
        let c = curve();
        let mut prev = -1.0;
        let mut v = 3.0;
        while v <= 12.0 {
            let p = c.power(v);
            assert!(p > prev, "power not increasing at v = {v}");
            assert!(
                (0.0..=2.0e6 + EPS).contains(&p),
                "power out of band at v = {v}"
            );
            prev = p;
            v += 0.5;
        }
    }

    #[test]
    fn ramp_follows_cube_law_midpoint() {
        // Hand-computed: at v = 7.5 with v_in = 3, v_r = 12, P_rated = 2e6:
        // P = 2e6 * (7.5^3 - 3^3) / (12^3 - 3^3)
        //   = 2e6 * (421.875 - 27) / (1728 - 27)
        //   = 2e6 * 394.875 / 1701.
        let c = curve();
        let expected = 2.0e6 * (421.875 - 27.0) / (1728.0 - 27.0);
        assert!((c.power(7.5) - expected).abs() < 1e-3);
    }

    #[test]
    fn region_classification_matches_power() {
        let c = curve();
        assert_eq!(c.region(1.0), Region::Idle);
        assert_eq!(c.region(6.0), Region::Ramp);
        assert_eq!(c.region(20.0), Region::Rated);
        assert_eq!(c.region(99.0), Region::Idle);
        assert_eq!(c.region(f64::NAN), Region::Idle);
    }

    #[test]
    fn non_finite_speed_is_zero_power() {
        let c = curve();
        assert!(c.power(f64::NAN).abs() < EPS);
        assert!(c.power(f64::INFINITY).abs() < EPS);
        assert!(c.power(f64::NEG_INFINITY).abs() < EPS);
    }

    #[test]
    fn accessors_round_trip() {
        let c = curve();
        assert!((c.cut_in() - 3.0).abs() < EPS);
        assert!((c.rated() - 12.0).abs() < EPS);
        assert!((c.cut_out() - 25.0).abs() < EPS);
        assert!((c.rated_power() - 2.0e6).abs() < EPS);
    }
}
