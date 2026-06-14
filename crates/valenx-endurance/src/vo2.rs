//! Oxygen uptake (VO2) as a function of mechanical power, and its
//! expression as a percentage of VO2max.
//!
//! ## Model
//!
//! For steady-state cycle ergometry, oxygen uptake rises essentially
//! linearly with external power output. This crate uses the standard
//! linear ramp
//!
//! ```text
//! vo2(power) = resting + slope * power
//! ```
//!
//! where `power` is in watts and `vo2` is in mL O2 per minute. The
//! default slope is `10.8` mL O2/min per watt — the classic ACSM /
//! ergometry oxygen-cost-of-cycling figure — and the default resting
//! offset is `300` mL/min (a representative resting VO2; the textbook
//! "unloaded pedalling" intercept is of similar magnitude). Both
//! defaults are configurable via [`Vo2Model::new`].
//!
//! Exercise intensity relative to an athlete's ceiling is the ratio of
//! the working VO2 to VO2max:
//!
//! ```text
//! percent_vo2max(vo2, vo2max) = 100 * vo2 / vo2max
//! ```

use serde::{Deserialize, Serialize};

use crate::error::{require_non_negative, require_positive, EnduranceError};

/// Default oxygen cost of cycling, in mL O2 per minute per watt.
pub const DEFAULT_O2_COST_ML_MIN_PER_W: f64 = 10.8;

/// Default resting / unloaded VO2 intercept, in mL O2 per minute.
pub const DEFAULT_RESTING_VO2_ML_MIN: f64 = 300.0;

/// A linear VO2-vs-power model for cycle ergometry.
///
/// Construct with [`Vo2Model::standard`] for the textbook defaults or
/// [`Vo2Model::new`] to validate custom coefficients.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Vo2Model {
    /// Resting / unloaded VO2 intercept, mL O2/min. Must be `>= 0`.
    pub resting_vo2_ml_min: f64,
    /// Oxygen cost of power, mL O2/min per watt. Must be `> 0`.
    pub o2_cost_ml_min_per_w: f64,
}

impl Default for Vo2Model {
    fn default() -> Self {
        Self::standard()
    }
}

impl Vo2Model {
    /// The textbook default model: resting `300` mL/min, slope `10.8`
    /// mL/min per watt.
    #[must_use]
    pub fn standard() -> Self {
        Self {
            resting_vo2_ml_min: DEFAULT_RESTING_VO2_ML_MIN,
            o2_cost_ml_min_per_w: DEFAULT_O2_COST_ML_MIN_PER_W,
        }
    }

    /// Construct and validate a custom linear VO2 model.
    ///
    /// # Errors
    ///
    /// Returns [`EnduranceError::OutOfDomain`] if `resting_vo2_ml_min` is
    /// negative or `o2_cost_ml_min_per_w` is not strictly positive, and
    /// [`EnduranceError::NotFinite`] for any non-finite argument.
    pub fn new(resting_vo2_ml_min: f64, o2_cost_ml_min_per_w: f64) -> Result<Self, EnduranceError> {
        require_non_negative(
            "resting_vo2_ml_min",
            resting_vo2_ml_min,
            "resting VO2 must be >= 0",
        )?;
        require_positive(
            "o2_cost_ml_min_per_w",
            o2_cost_ml_min_per_w,
            "oxygen cost must be > 0",
        )?;
        Ok(Self {
            resting_vo2_ml_min,
            o2_cost_ml_min_per_w,
        })
    }

    /// Steady-state oxygen uptake (mL O2/min) at a given mechanical
    /// `power_w` (watts).
    ///
    /// Implements `resting + slope * power`.
    ///
    /// # Errors
    ///
    /// Returns [`EnduranceError::OutOfDomain`] if `power_w` is negative
    /// and [`EnduranceError::NotFinite`] if it is not finite.
    pub fn vo2_at(&self, power_w: f64) -> Result<f64, EnduranceError> {
        require_non_negative("power_w", power_w, "power must be >= 0")?;
        Ok(self.resting_vo2_ml_min + self.o2_cost_ml_min_per_w * power_w)
    }
}

/// Free-function form of [`Vo2Model::vo2_at`] using the standard model.
///
/// Returns oxygen uptake in mL O2/min for a given mechanical power in
/// watts. Equivalent to `Vo2Model::standard().vo2_at(power_w)`.
///
/// # Errors
///
/// Propagates the errors of [`Vo2Model::vo2_at`].
pub fn vo2_at(power_w: f64) -> Result<f64, EnduranceError> {
    Vo2Model::standard().vo2_at(power_w)
}

/// Exercise intensity as a percentage of VO2max.
///
/// Computes `100 * vo2 / vo2max`. Both arguments must use the same units
/// (both mL/min, or both mL/kg/min); the ratio is dimensionless.
///
/// # Arguments
///
/// - `vo2` — current oxygen uptake (`>= 0`).
/// - `vo2max` — the athlete's maximal oxygen uptake (`> 0`).
///
/// # Errors
///
/// Returns [`EnduranceError::OutOfDomain`] if `vo2` is negative or
/// `vo2max` is not strictly positive, and [`EnduranceError::NotFinite`]
/// for any non-finite argument.
pub fn percent_vo2max(vo2: f64, vo2max: f64) -> Result<f64, EnduranceError> {
    require_non_negative("vo2", vo2, "VO2 must be >= 0")?;
    require_positive("vo2max", vo2max, "VO2max must be > 0")?;
    Ok(100.0 * vo2 / vo2max)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn vo2_at_zero_power_is_resting() {
        let m = Vo2Model::standard();
        let v = m.vo2_at(0.0).unwrap();
        assert!(
            (v - m.resting_vo2_ml_min).abs() < EPS,
            "vo2 at 0 W = {v}, expected resting {}",
            m.resting_vo2_ml_min
        );
    }

    #[test]
    fn vo2_is_linear_in_power() {
        let m = Vo2Model::standard();
        let v100 = m.vo2_at(100.0).unwrap();
        let v200 = m.vo2_at(200.0).unwrap();
        let v300 = m.vo2_at(300.0).unwrap();
        // Equal power increments give equal VO2 increments.
        let d_low = v200 - v100;
        let d_high = v300 - v200;
        assert!(
            (d_low - d_high).abs() < EPS,
            "increments differ: {d_low} vs {d_high}"
        );
        // Slope matches the configured oxygen cost.
        assert!(
            (d_low / 100.0 - m.o2_cost_ml_min_per_w).abs() < EPS,
            "recovered slope {} != {}",
            d_low / 100.0,
            m.o2_cost_ml_min_per_w
        );
    }

    #[test]
    fn vo2_known_value() {
        // resting 300 + 10.8 * 250 = 300 + 2700 = 3000 mL/min.
        let v = Vo2Model::standard().vo2_at(250.0).unwrap();
        assert!((v - 3000.0).abs() < 1e-9, "vo2 = {v}, expected 3000");
    }

    #[test]
    fn vo2_strictly_increases_with_power() {
        let m = Vo2Model::standard();
        let mut prev = m.vo2_at(0.0).unwrap();
        let mut p = 10.0;
        while p <= 500.0 {
            let v = m.vo2_at(p).unwrap();
            assert!(v > prev, "vo2 not increasing at {p} W: {v} <= {prev}");
            prev = v;
            p += 10.0;
        }
    }

    #[test]
    fn percent_vo2max_known_value() {
        let pct = percent_vo2max(2000.0, 4000.0).unwrap();
        assert!((pct - 50.0).abs() < EPS, "percent = {pct}, expected 50");
    }

    #[test]
    fn percent_vo2max_at_max_is_hundred() {
        let pct = percent_vo2max(3500.0, 3500.0).unwrap();
        assert!((pct - 100.0).abs() < EPS, "percent = {pct}, expected 100");
    }

    #[test]
    fn percent_vo2max_increases_with_vo2() {
        let lo = percent_vo2max(1000.0, 4000.0).unwrap();
        let hi = percent_vo2max(3000.0, 4000.0).unwrap();
        assert!(hi > lo, "percent should rise with vo2: {hi} <= {lo}");
    }

    #[test]
    fn free_function_matches_standard() {
        let a = vo2_at(175.0).unwrap();
        let b = Vo2Model::standard().vo2_at(175.0).unwrap();
        assert!((a - b).abs() < EPS, "free fn {a} != method {b}");
    }

    #[test]
    fn vo2_rejects_negative_power() {
        assert!(Vo2Model::standard().vo2_at(-1.0).is_err());
    }

    #[test]
    fn new_rejects_bad_parameters() {
        assert!(Vo2Model::new(-1.0, 10.8).is_err());
        assert!(Vo2Model::new(300.0, 0.0).is_err());
        assert!(Vo2Model::new(300.0, -5.0).is_err());
        assert!(Vo2Model::new(f64::NAN, 10.8).is_err());
    }

    #[test]
    fn percent_vo2max_rejects_bad_inputs() {
        assert!(percent_vo2max(-1.0, 4000.0).is_err());
        assert!(percent_vo2max(2000.0, 0.0).is_err());
        assert!(percent_vo2max(2000.0, -4000.0).is_err());
    }
}
