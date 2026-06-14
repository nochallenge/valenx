//! Oxygen transport: the O2-hemoglobin dissociation curve, arterial
//! oxygen content, and oxygen delivery.
//!
//! ## Model
//!
//! Hemoglobin oxygen saturation as a function of oxygen partial pressure
//! `po2` (mmHg) is modelled by the **Hill equation**
//!
//! ```text
//! saturation(po2) = po2^n / (p50^n + po2^n)
//! ```
//!
//! where `p50` is the partial pressure at which hemoglobin is half
//! saturated and `n` is the Hill cooperativity coefficient. The function
//! is sigmoidal: it is strictly increasing in `po2`, equals exactly
//! `0.5` at `po2 == p50`, and tends to `1.0` as `po2 -> infinity`.
//! Standard physiological defaults are `p50 == 26.6` mmHg and
//! `n == 2.7` (see [`HillCurve::standard`]).
//!
//! Arterial oxygen **content** (mL O2 per dL of blood) sums the oxygen
//! bound to hemoglobin and the small amount dissolved in plasma:
//!
//! ```text
//! cao2 = 1.34 * hb * saturation + 0.003 * po2
//! ```
//!
//! `1.34` mL O2/g is Hüfner's constant (oxygen-carrying capacity of
//! hemoglobin) and `0.003` mL O2/(dL*mmHg) is the plasma solubility
//! coefficient.
//!
//! Oxygen **delivery** (mL O2 per minute) is the product of cardiac
//! output and arterial content (with a factor of 10 to convert dL to L):
//!
//! ```text
//! do2 = cardiac_output * cao2 * 10
//! ```

use serde::{Deserialize, Serialize};

use crate::error::{require_non_negative, require_positive, EnduranceError};

/// Default `p50`: the oxygen partial pressure (mmHg) at which adult
/// hemoglobin is 50% saturated under standard conditions.
pub const DEFAULT_P50_MMHG: f64 = 26.6;

/// Default Hill cooperativity coefficient `n` for adult hemoglobin.
pub const DEFAULT_HILL_N: f64 = 2.7;

/// Hüfner's constant: oxygen bound per gram of fully saturated
/// hemoglobin, in mL O2 per gram.
pub const HUFNER_CONSTANT_ML_PER_G: f64 = 1.34;

/// Plasma oxygen solubility coefficient, in mL O2 per dL of blood per
/// mmHg of dissolved-oxygen partial pressure.
pub const PLASMA_SOLUBILITY_ML_DL_MMHG: f64 = 0.003;

/// An O2-hemoglobin dissociation curve parameterised by its half-
/// saturation pressure `p50` and Hill coefficient `n`.
///
/// Construct with [`HillCurve::standard`] for the physiological defaults
/// or [`HillCurve::new`] to validate custom parameters (a right- or
/// left-shifted curve, fetal hemoglobin, etc.).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct HillCurve {
    /// Half-saturation oxygen partial pressure, mmHg. Strictly positive.
    pub p50_mmhg: f64,
    /// Hill cooperativity coefficient (dimensionless). Strictly positive.
    pub n: f64,
}

impl Default for HillCurve {
    fn default() -> Self {
        Self::standard()
    }
}

impl HillCurve {
    /// The standard adult curve: `p50 == 26.6` mmHg, `n == 2.7`.
    #[must_use]
    pub fn standard() -> Self {
        Self {
            p50_mmhg: DEFAULT_P50_MMHG,
            n: DEFAULT_HILL_N,
        }
    }

    /// Construct a curve from custom `p50` (mmHg) and Hill coefficient
    /// `n`, validating both.
    ///
    /// # Errors
    ///
    /// Returns [`EnduranceError::OutOfDomain`] if `p50_mmhg <= 0` or
    /// `n <= 0`, and [`EnduranceError::NotFinite`] if either argument is
    /// not finite.
    pub fn new(p50_mmhg: f64, n: f64) -> Result<Self, EnduranceError> {
        require_positive("p50_mmhg", p50_mmhg, "half-saturation pressure must be > 0")?;
        require_positive("n", n, "Hill coefficient must be > 0")?;
        Ok(Self { p50_mmhg, n })
    }

    /// Fractional hemoglobin saturation in `0.0..=1.0` for a given oxygen
    /// partial pressure `po2` (mmHg).
    ///
    /// Implements `po2^n / (p50^n + po2^n)`. At `po2 == 0` the result is
    /// exactly `0.0`; at `po2 == p50` it is exactly `0.5`; as `po2`
    /// grows it approaches `1.0`.
    ///
    /// # Errors
    ///
    /// Returns [`EnduranceError::OutOfDomain`] if `po2_mmhg` is negative
    /// and [`EnduranceError::NotFinite`] if it is not finite.
    pub fn saturation(&self, po2_mmhg: f64) -> Result<f64, EnduranceError> {
        require_non_negative("po2_mmhg", po2_mmhg, "partial pressure must be >= 0")?;
        if po2_mmhg == 0.0 {
            return Ok(0.0);
        }
        let num = po2_mmhg.powf(self.n);
        let den = self.p50_mmhg.powf(self.n) + num;
        Ok(num / den)
    }
}

/// Free-function form of the Hill saturation using the standard adult
/// curve (`p50 == 26.6` mmHg, `n == 2.7`).
///
/// Equivalent to `HillCurve::standard().saturation(po2_mmhg)`.
///
/// # Errors
///
/// Propagates the errors of [`HillCurve::saturation`].
pub fn saturation(po2_mmhg: f64) -> Result<f64, EnduranceError> {
    HillCurve::standard().saturation(po2_mmhg)
}

/// Arterial oxygen content `cao2`, in mL O2 per dL of blood.
///
/// Computes `1.34 * hb * sat + 0.003 * po2`, the sum of hemoglobin-bound
/// and plasma-dissolved oxygen.
///
/// # Arguments
///
/// - `hb_g_dl` — hemoglobin concentration, g/dL (typical adult ~12-17).
/// - `sat` — fractional hemoglobin saturation in `0.0..=1.0`.
/// - `po2_mmhg` — arterial oxygen partial pressure, mmHg.
///
/// # Errors
///
/// Returns [`EnduranceError::OutOfDomain`] if `hb_g_dl` or `po2_mmhg` is
/// negative or if `sat` is outside `0.0..=1.0`, and
/// [`EnduranceError::NotFinite`] for any non-finite argument.
pub fn cao2(hb_g_dl: f64, sat: f64, po2_mmhg: f64) -> Result<f64, EnduranceError> {
    use crate::error::require_in_closed;
    require_non_negative("hb_g_dl", hb_g_dl, "hemoglobin must be >= 0")?;
    require_in_closed("sat", sat, 0.0, 1.0, "saturation fraction must be in 0..=1")?;
    require_non_negative("po2_mmhg", po2_mmhg, "partial pressure must be >= 0")?;
    Ok(HUFNER_CONSTANT_ML_PER_G * hb_g_dl * sat + PLASMA_SOLUBILITY_ML_DL_MMHG * po2_mmhg)
}

/// Oxygen delivery `do2`, in mL O2 per minute.
///
/// Computes `cardiac_output * cao2 * 10`, where the factor of 10
/// converts arterial content from mL/dL to mL/L so it matches a cardiac
/// output expressed in litres per minute.
///
/// # Arguments
///
/// - `cardiac_output_l_min` — cardiac output, L/min (typical resting ~5).
/// - `cao2_ml_dl` — arterial oxygen content, mL O2/dL (see [`cao2`]).
///
/// # Errors
///
/// Returns [`EnduranceError::OutOfDomain`] if either argument is
/// negative and [`EnduranceError::NotFinite`] if either is not finite.
pub fn do2(cardiac_output_l_min: f64, cao2_ml_dl: f64) -> Result<f64, EnduranceError> {
    require_non_negative(
        "cardiac_output_l_min",
        cardiac_output_l_min,
        "cardiac output must be >= 0",
    )?;
    require_non_negative("cao2_ml_dl", cao2_ml_dl, "arterial content must be >= 0")?;
    Ok(cardiac_output_l_min * cao2_ml_dl * 10.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn saturation_at_p50_is_one_half() {
        let curve = HillCurve::standard();
        let s = curve.saturation(curve.p50_mmhg).unwrap();
        assert!((s - 0.5).abs() < EPS, "saturation(p50) = {s}, expected 0.5");
    }

    #[test]
    fn free_function_matches_standard_curve() {
        let p = 40.0;
        let a = saturation(p).unwrap();
        let b = HillCurve::standard().saturation(p).unwrap();
        assert!((a - b).abs() < EPS, "free fn {a} != method {b}");
    }

    #[test]
    fn saturation_at_zero_is_zero() {
        let s = HillCurve::standard().saturation(0.0).unwrap();
        assert!(s.abs() < EPS, "saturation(0) = {s}, expected 0");
    }

    #[test]
    fn saturation_strictly_increases_with_po2() {
        let curve = HillCurve::standard();
        let mut prev = curve.saturation(0.0).unwrap();
        let mut po2 = 1.0_f64;
        while po2 <= 200.0 {
            let s = curve.saturation(po2).unwrap();
            assert!(
                s > prev,
                "saturation not strictly increasing at po2 = {po2}: {s} <= {prev}"
            );
            prev = s;
            po2 += 1.0;
        }
    }

    #[test]
    fn saturation_approaches_one_at_high_po2() {
        let curve = HillCurve::standard();
        let s = curve.saturation(2000.0).unwrap();
        assert!(s > 0.999, "saturation(2000) = {s}, expected near 1.0");
        assert!(s <= 1.0, "saturation must not exceed 1.0, got {s}");
    }

    #[test]
    fn saturation_is_bounded_in_unit_interval() {
        let curve = HillCurve::standard();
        for po2 in [0.0, 5.0, 26.6, 60.0, 100.0, 500.0] {
            let s = curve.saturation(po2).unwrap();
            assert!(
                (0.0..=1.0).contains(&s),
                "saturation {s} out of [0,1] at {po2}"
            );
        }
    }

    #[test]
    fn hand_computed_saturation_value() {
        // po2 = 26.6 * 2^(1/2.7) gives saturation exactly 2/3:
        //   r = 2^(1/n) => (r*p50)^n / (p50^n + (r*p50)^n)
        //              = 2*p50^n / (p50^n + 2*p50^n) = 2/3.
        let curve = HillCurve::standard();
        let po2 = curve.p50_mmhg * 2.0_f64.powf(1.0 / curve.n);
        let s = curve.saturation(po2).unwrap();
        assert!(
            (s - 2.0 / 3.0).abs() < 1e-9,
            "saturation = {s}, expected 2/3"
        );
    }

    #[test]
    fn cao2_increases_with_hemoglobin() {
        let sat = 0.97;
        let po2 = 100.0;
        let low = cao2(10.0, sat, po2).unwrap();
        let high = cao2(15.0, sat, po2).unwrap();
        assert!(high > low, "cao2 should rise with hb: {high} <= {low}");
    }

    #[test]
    fn cao2_known_value() {
        // hb 15 g/dL, sat 1.0, po2 100 mmHg:
        //   1.34*15*1 + 0.003*100 = 20.1 + 0.3 = 20.4 mL/dL.
        let c = cao2(15.0, 1.0, 100.0).unwrap();
        assert!((c - 20.4).abs() < 1e-9, "cao2 = {c}, expected 20.4");
    }

    #[test]
    fn do2_equals_co_times_cao2_times_ten() {
        let co = 5.0;
        let content = 20.4;
        let delivery = do2(co, content).unwrap();
        let expected = co * content * 10.0;
        assert!(
            (delivery - expected).abs() < EPS,
            "do2 = {delivery}, expected {expected}"
        );
        // Resting human reference: ~1020 mL/min.
        assert!(
            (delivery - 1020.0).abs() < EPS,
            "do2 = {delivery}, expected 1020"
        );
    }

    #[test]
    fn new_rejects_non_positive_parameters() {
        assert!(HillCurve::new(0.0, 2.7).is_err());
        assert!(HillCurve::new(-1.0, 2.7).is_err());
        assert!(HillCurve::new(26.6, 0.0).is_err());
        assert!(HillCurve::new(f64::NAN, 2.7).is_err());
    }

    #[test]
    fn saturation_rejects_negative_po2() {
        assert!(HillCurve::standard().saturation(-1.0).is_err());
    }

    #[test]
    fn cao2_rejects_out_of_range_saturation() {
        assert!(cao2(15.0, 1.5, 100.0).is_err());
        assert!(cao2(15.0, -0.1, 100.0).is_err());
        assert!(cao2(-1.0, 0.9, 100.0).is_err());
    }

    #[test]
    fn do2_rejects_negative_inputs() {
        assert!(do2(-1.0, 20.0).is_err());
        assert!(do2(5.0, -20.0).is_err());
    }
}
