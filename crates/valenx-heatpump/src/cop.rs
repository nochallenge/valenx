//! Carnot coefficient of performance (COP) for heat pumps and chillers.
//!
//! A heat pump moving heat from a cold reservoir at `T_c` to a hot
//! reservoir at `T_h` (both *absolute* temperatures, in kelvin) has a
//! thermodynamic upper bound on its efficiency set by the reversible
//! Carnot cycle:
//!
//! ```text
//! COP_heat = T_h / (T_h - T_c)        (heating: useful output is Q_h)
//! COP_cool = T_c / (T_h - T_c)        (cooling: useful output is Q_c)
//! ```
//!
//! Two consequences fall straight out of the algebra and are exposed as
//! first-class, tested relations. **The unity gap** is
//! `COP_heat = COP_cool + 1`: the same compressor work `W` delivers `Q_c`
//! of cooling but `Q_h = Q_c + W` of heating, because the work itself ends
//! up as heat in the hot reservoir. And the **heating COP is always above
//! one** — `T_h / (T_h - T_c) > 1` whenever `T_h > T_c > 0`, so a heat
//! pump always delivers more heat than the electrical work it consumes
//! (unlike a resistive heater, whose COP is exactly 1).
//!
//! Real machines never reach the Carnot limit; the
//! [`derated`](CarnotCop::derated) helper multiplies by a *Carnot
//! fraction* (second-law efficiency) in `(0, 1]` to model irreversibilities.

use serde::{Deserialize, Serialize};

use crate::error::{check_temperature_k, HeatPumpError, Result};

/// The two reversible (Carnot) coefficients of performance for a heat
/// pump operating between a cold and a hot reservoir.
///
/// Construct with [`CarnotCop::new`], which validates both temperatures
/// and the lift. The struct is a plain value type; the fields are the
/// computed COPs and are guaranteed finite and positive on a
/// successfully constructed instance.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CarnotCop {
    /// Cold-reservoir absolute temperature, in kelvin.
    pub t_cold_k: f64,
    /// Hot-reservoir absolute temperature, in kelvin.
    pub t_hot_k: f64,
    /// Reversible heating coefficient of performance,
    /// `T_h / (T_h - T_c)`. Always strictly greater than 1.
    pub cop_heat: f64,
    /// Reversible cooling coefficient of performance,
    /// `T_c / (T_h - T_c)`. Always strictly greater than 0.
    pub cop_cool: f64,
}

impl CarnotCop {
    /// Build the Carnot COPs from two absolute reservoir temperatures.
    ///
    /// Both `t_cold_k` and `t_hot_k` must be finite and strictly above
    /// absolute zero, and the lift `t_hot_k - t_cold_k` must be strictly
    /// positive.
    ///
    /// # Errors
    ///
    /// Returns [`HeatPumpError::Invalid`] if either temperature is not a
    /// finite positive kelvin value, or [`HeatPumpError::DegenerateLift`]
    /// if `t_hot_k <= t_cold_k` (the lift would be zero or negative,
    /// making the COP undefined).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_heatpump::CarnotCop;
    /// // Heating from 0 °C (273.15 K) outside to 35 °C (308.15 K) supply.
    /// let c = CarnotCop::new(273.15, 308.15).unwrap();
    /// assert!(c.cop_heat > 1.0);
    /// assert!((c.cop_heat - c.cop_cool - 1.0).abs() < 1e-9);
    /// ```
    pub fn new(t_cold_k: f64, t_hot_k: f64) -> Result<Self> {
        let t_cold_k = check_temperature_k("t_cold_k", t_cold_k)?;
        let t_hot_k = check_temperature_k("t_hot_k", t_hot_k)?;
        let lift = t_hot_k - t_cold_k;
        if lift <= 0.0 {
            return Err(HeatPumpError::degenerate_lift(t_hot_k, t_cold_k));
        }
        Ok(Self {
            t_cold_k,
            t_hot_k,
            cop_heat: t_hot_k / lift,
            cop_cool: t_cold_k / lift,
        })
    }

    /// The temperature lift `T_h - T_c`, in kelvin. Always strictly
    /// positive on a constructed instance.
    pub fn lift_k(&self) -> f64 {
        self.t_hot_k - self.t_cold_k
    }

    /// The reversible heating COP after applying a *Carnot fraction*
    /// (second-law efficiency) in `(0, 1]` to model real-machine
    /// irreversibilities.
    ///
    /// A `fraction` of `1.0` returns the ideal [`cop_heat`](Self::cop_heat);
    /// a typical real air-source heat pump sits around `0.4`–`0.6`.
    ///
    /// # Errors
    ///
    /// [`HeatPumpError::Invalid`] if `fraction` is not finite or lies
    /// outside the half-open range `(0, 1]`.
    pub fn derated_heat(&self, fraction: f64) -> Result<f64> {
        Ok(self.cop_heat * check_carnot_fraction(fraction)?)
    }

    /// The reversible cooling COP after applying a *Carnot fraction*
    /// (second-law efficiency) in `(0, 1]`. See
    /// [`derated_heat`](Self::derated_heat).
    ///
    /// # Errors
    ///
    /// [`HeatPumpError::Invalid`] if `fraction` is not finite or lies
    /// outside `(0, 1]`.
    pub fn derated_cool(&self, fraction: f64) -> Result<f64> {
        Ok(self.cop_cool * check_carnot_fraction(fraction)?)
    }

    /// Apply a Carnot fraction to *both* COPs at once, returning a new
    /// pair of derated coefficients.
    ///
    /// Note that the unity-gap identity `COP_heat = COP_cool + 1` holds
    /// for the *ideal* COPs but **not** for the derated ones: scaling
    /// both by the same `f < 1` gives `f * cop_heat - f * cop_cool = f`,
    /// not `1`. The returned [`Derated`] therefore reports the gap it
    /// actually has.
    ///
    /// # Errors
    ///
    /// [`HeatPumpError::Invalid`] if `fraction` is not finite or lies
    /// outside `(0, 1]`.
    pub fn derated(&self, fraction: f64) -> Result<Derated> {
        let f = check_carnot_fraction(fraction)?;
        Ok(Derated {
            fraction: f,
            cop_heat: self.cop_heat * f,
            cop_cool: self.cop_cool * f,
        })
    }
}

/// A pair of real-machine COPs obtained from the Carnot limits by a
/// second-law (Carnot) fraction. Returned by [`CarnotCop::derated`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Derated {
    /// The Carnot fraction (second-law efficiency) that was applied,
    /// in `(0, 1]`.
    pub fraction: f64,
    /// Derated heating COP, `fraction * cop_heat_ideal`.
    pub cop_heat: f64,
    /// Derated cooling COP, `fraction * cop_cool_ideal`.
    pub cop_cool: f64,
}

impl Derated {
    /// The heating-minus-cooling COP gap for these *derated* values.
    ///
    /// Equals the applied [`fraction`](Self::fraction) (because
    /// `f * cop_heat - f * cop_cool = f * 1 = f`), which is at most 1 —
    /// so derating shrinks the ideal unity gap.
    pub fn gap(&self) -> f64 {
        self.cop_heat - self.cop_cool
    }
}

/// Validate a Carnot fraction (second-law efficiency): finite and within
/// the half-open range `(0, 1]`.
///
/// # Errors
///
/// [`HeatPumpError::Invalid`] otherwise.
fn check_carnot_fraction(fraction: f64) -> Result<f64> {
    if !fraction.is_finite() {
        return Err(HeatPumpError::invalid(
            "carnot_fraction",
            "must be a finite number",
        ));
    }
    if fraction <= 0.0 || fraction > 1.0 {
        return Err(HeatPumpError::invalid(
            "carnot_fraction",
            format!("must lie in the half-open range (0, 1], got {fraction}"),
        ));
    }
    Ok(fraction)
}

/// The reversible heating coefficient of performance from two absolute
/// temperatures: `COP_heat = T_h / (T_h - T_c)`.
///
/// A free-function shortcut for callers that do not want the
/// [`CarnotCop`] struct. Validation is identical to
/// [`CarnotCop::new`].
///
/// # Errors
///
/// See [`CarnotCop::new`].
///
/// # Examples
///
/// ```
/// use valenx_heatpump::cop::carnot_cop_heat;
/// let cop = carnot_cop_heat(273.15, 308.15).unwrap();
/// assert!(cop > 1.0);
/// ```
pub fn carnot_cop_heat(t_cold_k: f64, t_hot_k: f64) -> Result<f64> {
    Ok(CarnotCop::new(t_cold_k, t_hot_k)?.cop_heat)
}

/// The reversible cooling coefficient of performance from two absolute
/// temperatures: `COP_cool = T_c / (T_h - T_c)`.
///
/// A free-function shortcut for callers that do not want the
/// [`CarnotCop`] struct. Validation is identical to
/// [`CarnotCop::new`].
///
/// # Errors
///
/// See [`CarnotCop::new`].
pub fn carnot_cop_cool(t_cold_k: f64, t_hot_k: f64) -> Result<f64> {
    Ok(CarnotCop::new(t_cold_k, t_hot_k)?.cop_cool)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point comparisons.
    const EPS: f64 = 1e-9;

    /// Ground truth: heating from a 0 °C source to a 35 °C sink.
    /// `T_c = 273.15 K`, `T_h = 308.15 K`, lift = 35 K.
    /// `COP_heat = 308.15 / 35 = 8.804285714...`,
    /// `COP_cool = 273.15 / 35 = 7.804285714...`.
    #[test]
    fn carnot_cop_matches_hand_computed_ground_truth() {
        let c = CarnotCop::new(273.15, 308.15).unwrap();
        assert!((c.cop_heat - 308.15 / 35.0).abs() < EPS);
        assert!((c.cop_cool - 273.15 / 35.0).abs() < EPS);
        assert!((c.cop_heat - 8.804_285_714_285_714).abs() < 1e-12);
        assert!((c.lift_k() - 35.0).abs() < EPS);
    }

    /// VALIDATE: `COP_heat = COP_cool + 1` for any valid pair.
    #[test]
    fn heating_cop_exceeds_cooling_cop_by_exactly_one() {
        for &(tc, th) in &[
            (273.15, 308.15),
            (250.0, 320.0),
            (300.0, 301.0),
            (200.0, 400.0),
            (4.0, 5.0),
        ] {
            let c = CarnotCop::new(tc, th).unwrap();
            assert!(
                (c.cop_heat - (c.cop_cool + 1.0)).abs() < EPS,
                "unity gap broken for T_c={tc}, T_h={th}"
            );
        }
    }

    /// VALIDATE: Carnot `COP_heat = T_h / (T_h - T_c)` exactly.
    #[test]
    fn carnot_heating_formula_is_th_over_lift() {
        let (tc, th) = (260.0, 300.0);
        let c = CarnotCop::new(tc, th).unwrap();
        assert!((c.cop_heat - th / (th - tc)).abs() < EPS);
        assert!((c.cop_cool - tc / (th - tc)).abs() < EPS);
    }

    /// VALIDATE: heating COP is always strictly greater than 1 (a heat
    /// pump beats a resistive heater).
    #[test]
    fn heating_cop_is_always_above_unity() {
        for &(tc, th) in &[(273.15, 308.15), (250.0, 251.0), (1.0, 1000.0)] {
            let c = CarnotCop::new(tc, th).unwrap();
            assert!(c.cop_heat > 1.0, "COP_heat <= 1 for T_c={tc}, T_h={th}");
        }
    }

    /// VALIDATE: COP falls monotonically as the temperature lift grows.
    /// Holding `T_c` fixed and raising `T_h` widens the lift, so both
    /// COPs must decrease.
    #[test]
    fn cop_falls_as_lift_grows() {
        let tc = 273.15;
        let mut prev_heat = f64::INFINITY;
        let mut prev_cool = f64::INFINITY;
        for delta in [5.0, 10.0, 20.0, 40.0, 80.0] {
            let c = CarnotCop::new(tc, tc + delta).unwrap();
            assert!(
                c.cop_heat < prev_heat,
                "COP_heat did not fall at lift {delta}"
            );
            assert!(
                c.cop_cool < prev_cool,
                "COP_cool did not fall at lift {delta}"
            );
            prev_heat = c.cop_heat;
            prev_cool = c.cop_cool;
        }
    }

    /// In the reversible limit the lift `T_h - T_c -> 0` sends both COPs
    /// to infinity; a tiny lift gives a huge but finite COP.
    #[test]
    fn tiny_lift_gives_large_cop() {
        let c = CarnotCop::new(300.0, 300.001).unwrap();
        assert!(c.cop_heat > 100_000.0);
        assert!(c.cop_heat.is_finite());
    }

    #[test]
    fn degenerate_and_inverted_lift_are_rejected() {
        // Equal reservoirs: zero lift.
        assert!(matches!(
            CarnotCop::new(300.0, 300.0),
            Err(HeatPumpError::DegenerateLift { .. })
        ));
        // Inverted: hot below cold.
        assert!(matches!(
            CarnotCop::new(310.0, 300.0),
            Err(HeatPumpError::DegenerateLift { .. })
        ));
    }

    #[test]
    fn non_positive_temperatures_are_rejected() {
        assert!(matches!(
            CarnotCop::new(0.0, 300.0),
            Err(HeatPumpError::Invalid { .. })
        ));
        assert!(matches!(
            CarnotCop::new(273.15, -1.0),
            Err(HeatPumpError::Invalid { .. })
        ));
    }

    /// Derating by 1.0 reproduces the ideal COPs exactly; the ideal
    /// unity gap is preserved only at `fraction == 1`.
    #[test]
    fn derating_by_unity_is_identity_and_preserves_gap() {
        let c = CarnotCop::new(273.15, 308.15).unwrap();
        let d = c.derated(1.0).unwrap();
        assert!((d.cop_heat - c.cop_heat).abs() < EPS);
        assert!((d.cop_cool - c.cop_cool).abs() < EPS);
        assert!((d.gap() - 1.0).abs() < EPS);
    }

    /// A second-law fraction below 1 scales both COPs and shrinks the
    /// gap to exactly the fraction.
    #[test]
    fn derating_below_unity_scales_and_shrinks_gap() {
        let c = CarnotCop::new(273.15, 308.15).unwrap();
        let f = 0.5;
        let d = c.derated(f).unwrap();
        assert!((d.cop_heat - f * c.cop_heat).abs() < EPS);
        assert!((d.cop_cool - f * c.cop_cool).abs() < EPS);
        assert!((d.gap() - f).abs() < EPS);
        assert!((c.derated_heat(f).unwrap() - d.cop_heat).abs() < EPS);
        assert!((c.derated_cool(f).unwrap() - d.cop_cool).abs() < EPS);
    }

    #[test]
    fn out_of_range_carnot_fraction_is_rejected() {
        let c = CarnotCop::new(273.15, 308.15).unwrap();
        assert!(c.derated(0.0).is_err());
        assert!(c.derated(-0.1).is_err());
        assert!(c.derated(1.5).is_err());
        assert!(c.derated(f64::NAN).is_err());
    }

    #[test]
    fn free_function_shortcuts_agree_with_struct() {
        let (tc, th) = (273.15, 308.15);
        let c = CarnotCop::new(tc, th).unwrap();
        assert!((carnot_cop_heat(tc, th).unwrap() - c.cop_heat).abs() < EPS);
        assert!((carnot_cop_cool(tc, th).unwrap() - c.cop_cool).abs() < EPS);
    }

    #[test]
    fn serde_round_trips() {
        let c = CarnotCop::new(273.15, 308.15).unwrap();
        let json = serde_json::to_string(&c).unwrap();
        let back: CarnotCop = serde_json::from_str(&json).unwrap();
        assert_eq!(c, back);
    }
}
