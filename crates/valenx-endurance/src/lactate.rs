//! Blood lactate as a function of exercise intensity.
//!
//! ## Model
//!
//! Blood-lactate concentration is modelled against exercise intensity
//! expressed as a fraction of VO2max in `0.0..=1.0`. Below the **lactate
//! threshold** (default `0.75`, i.e. 75% of VO2max) clearance keeps pace
//! with production and lactate sits near a resting **baseline** (default
//! `1.0` mmol/L). Above the threshold production outstrips clearance and
//! lactate rises sharply; this crate uses an exponential rise in the
//! supra-threshold intensity excess:
//!
//! ```text
//! intensity <= threshold:  lactate = baseline
//! intensity >  threshold:  lactate = baseline
//!                                   + scale * (exp(steepness * x) - 1)
//!                          where x = intensity - threshold
//! ```
//!
//! At `intensity == threshold` the two branches agree (the exponential
//! excess is zero), so the curve is continuous. The result is a flat
//! baseline that turns up into the familiar "lactate hockey-stick".
//!
//! This is a deliberately simple textbook shape, not a metabolic model
//! of any individual athlete.

use serde::{Deserialize, Serialize};

use crate::error::{require_in_closed, require_non_negative, require_positive, EnduranceError};

/// Default lactate threshold as a fraction of VO2max.
pub const DEFAULT_LACTATE_THRESHOLD: f64 = 0.75;

/// Default resting blood-lactate baseline, mmol/L.
pub const DEFAULT_BASELINE_MMOL_L: f64 = 1.0;

/// Default supra-threshold rise scale, mmol/L.
pub const DEFAULT_RISE_SCALE_MMOL_L: f64 = 1.0;

/// Default exponential steepness of the supra-threshold rise
/// (dimensionless multiplier on the intensity excess).
pub const DEFAULT_RISE_STEEPNESS: f64 = 14.0;

/// Parameters of the lactate-vs-intensity curve.
///
/// Construct with [`LactateModel::standard`] for the textbook defaults
/// or [`LactateModel::new`] to validate custom parameters.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LactateModel {
    /// Lactate threshold as a fraction of VO2max, in `0.0..=1.0`.
    pub threshold: f64,
    /// Resting baseline lactate, mmol/L. Must be `>= 0`.
    pub baseline_mmol_l: f64,
    /// Scale of the supra-threshold exponential excess, mmol/L.
    /// Must be `> 0`.
    pub rise_scale_mmol_l: f64,
    /// Exponential steepness of the supra-threshold rise. Must be `> 0`.
    pub rise_steepness: f64,
}

impl Default for LactateModel {
    fn default() -> Self {
        Self::standard()
    }
}

impl LactateModel {
    /// The textbook default curve: threshold `0.75`, baseline `1.0`
    /// mmol/L, rise scale `1.0` mmol/L, steepness `14.0`.
    #[must_use]
    pub fn standard() -> Self {
        Self {
            threshold: DEFAULT_LACTATE_THRESHOLD,
            baseline_mmol_l: DEFAULT_BASELINE_MMOL_L,
            rise_scale_mmol_l: DEFAULT_RISE_SCALE_MMOL_L,
            rise_steepness: DEFAULT_RISE_STEEPNESS,
        }
    }

    /// Construct and validate a custom lactate model.
    ///
    /// # Errors
    ///
    /// Returns [`EnduranceError::OutOfDomain`] if `threshold` is outside
    /// `0.0..=1.0`, `baseline_mmol_l` is negative, or either of
    /// `rise_scale_mmol_l` / `rise_steepness` is not strictly positive;
    /// and [`EnduranceError::NotFinite`] for any non-finite argument.
    pub fn new(
        threshold: f64,
        baseline_mmol_l: f64,
        rise_scale_mmol_l: f64,
        rise_steepness: f64,
    ) -> Result<Self, EnduranceError> {
        require_in_closed(
            "threshold",
            threshold,
            0.0,
            1.0,
            "threshold must be in 0..=1",
        )?;
        require_non_negative("baseline_mmol_l", baseline_mmol_l, "baseline must be >= 0")?;
        require_positive(
            "rise_scale_mmol_l",
            rise_scale_mmol_l,
            "rise scale must be > 0",
        )?;
        require_positive("rise_steepness", rise_steepness, "steepness must be > 0")?;
        Ok(Self {
            threshold,
            baseline_mmol_l,
            rise_scale_mmol_l,
            rise_steepness,
        })
    }

    /// Blood lactate (mmol/L) at a given exercise `intensity`, expressed
    /// as a fraction of VO2max in `0.0..=1.0`.
    ///
    /// Returns the baseline at or below threshold, and the baseline plus
    /// an exponential excess above it. The curve is continuous at the
    /// threshold and monotonically non-decreasing in intensity.
    ///
    /// # Errors
    ///
    /// Returns [`EnduranceError::OutOfDomain`] if `intensity` is outside
    /// `0.0..=1.0` and [`EnduranceError::NotFinite`] if it is not finite.
    pub fn lactate_at(&self, intensity: f64) -> Result<f64, EnduranceError> {
        require_in_closed(
            "intensity",
            intensity,
            0.0,
            1.0,
            "intensity must be in 0..=1",
        )?;
        if intensity <= self.threshold {
            Ok(self.baseline_mmol_l)
        } else {
            let excess = intensity - self.threshold;
            let rise = self.rise_scale_mmol_l * ((self.rise_steepness * excess).exp() - 1.0);
            Ok(self.baseline_mmol_l + rise)
        }
    }
}

/// Free-function form of [`LactateModel::lactate_at`] using the standard
/// curve.
///
/// Equivalent to `LactateModel::standard().lactate_at(intensity)`.
///
/// # Errors
///
/// Propagates the errors of [`LactateModel::lactate_at`].
pub fn lactate_at(intensity: f64) -> Result<f64, EnduranceError> {
    LactateModel::standard().lactate_at(intensity)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn baseline_below_threshold() {
        let m = LactateModel::standard();
        let v = m.lactate_at(0.5).unwrap();
        assert!(
            (v - m.baseline_mmol_l).abs() < EPS,
            "below-threshold lactate = {v}, expected baseline {}",
            m.baseline_mmol_l
        );
    }

    #[test]
    fn equals_baseline_exactly_at_threshold() {
        let m = LactateModel::standard();
        let v = m.lactate_at(m.threshold).unwrap();
        assert!(
            (v - m.baseline_mmol_l).abs() < EPS,
            "at-threshold lactate = {v}, expected baseline {}",
            m.baseline_mmol_l
        );
    }

    #[test]
    fn continuous_at_threshold() {
        let m = LactateModel::standard();
        let below = m.lactate_at(m.threshold - 1e-7).unwrap();
        let at = m.lactate_at(m.threshold).unwrap();
        let above = m.lactate_at(m.threshold + 1e-7).unwrap();
        assert!(
            (below - at).abs() < 1e-5,
            "discontinuity below: {below} vs {at}"
        );
        assert!(
            (above - at).abs() < 1e-4,
            "discontinuity above: {above} vs {at}"
        );
    }

    #[test]
    fn above_threshold_exceeds_below_threshold() {
        let m = LactateModel::standard();
        let below = m.lactate_at(0.5).unwrap();
        let above = m.lactate_at(0.9).unwrap();
        assert!(
            above > below,
            "above-threshold lactate {above} should exceed below-threshold {below}"
        );
    }

    #[test]
    fn rises_with_intensity_above_threshold() {
        let m = LactateModel::standard();
        let mut prev = m.lactate_at(m.threshold).unwrap();
        let mut x = m.threshold + 0.01;
        while x <= 1.0 + EPS {
            let v = m.lactate_at(x.min(1.0)).unwrap();
            assert!(
                v >= prev,
                "lactate should be non-decreasing at intensity {x}: {v} < {prev}"
            );
            prev = v;
            x += 0.01;
        }
        // And the very top should be strictly above the threshold value.
        let top = m.lactate_at(1.0).unwrap();
        let at_thr = m.lactate_at(m.threshold).unwrap();
        assert!(
            top > at_thr,
            "max-intensity lactate {top} not above threshold {at_thr}"
        );
    }

    #[test]
    fn known_supra_threshold_value() {
        // Custom model with steepness chosen so exp() is exact:
        // baseline 1, scale 2, steepness 1, threshold 0.5.
        // At intensity 0.5 + ln(2) is out of range, so use x where
        // steepness*excess = ln(2): excess = ln(2) ~ 0.693 > 0.5 range.
        // Instead pick steepness = 2, threshold 0.5, intensity 0.5+0.5*ln(2)?
        // Simpler: threshold 0.0, steepness ln(2)/0.5 so at 0.5 exponent=ln2.
        let m = LactateModel::new(0.0, 1.0, 2.0, (2.0_f64).ln() / 0.5).unwrap();
        // exponent at intensity 0.5 = (ln2/0.5) * 0.5 = ln2 => exp = 2.
        // lactate = 1 + 2*(2 - 1) = 3.
        let v = m.lactate_at(0.5).unwrap();
        assert!((v - 3.0).abs() < 1e-9, "lactate = {v}, expected 3.0");
    }

    #[test]
    fn free_function_matches_standard() {
        let a = lactate_at(0.85).unwrap();
        let b = LactateModel::standard().lactate_at(0.85).unwrap();
        assert!((a - b).abs() < EPS, "free fn {a} != method {b}");
    }

    #[test]
    fn rejects_out_of_range_intensity() {
        let m = LactateModel::standard();
        assert!(m.lactate_at(-0.01).is_err());
        assert!(m.lactate_at(1.01).is_err());
        assert!(m.lactate_at(f64::INFINITY).is_err());
    }

    #[test]
    fn new_rejects_bad_parameters() {
        assert!(LactateModel::new(1.5, 1.0, 1.0, 14.0).is_err());
        assert!(LactateModel::new(0.75, -1.0, 1.0, 14.0).is_err());
        assert!(LactateModel::new(0.75, 1.0, 0.0, 14.0).is_err());
        assert!(LactateModel::new(0.75, 1.0, 1.0, 0.0).is_err());
    }
}
