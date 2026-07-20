//! Derived dose-response metrics: `ECx` concentrations, the curve's
//! local slope, and the `EC90 / EC10` dynamic range.
//!
//! These all follow in closed form from the Hill equation
//! `f(C) = C^n / (EC50^n + C^n)`.
//!
//! ## `ECx`
//!
//! The concentration eliciting `x` percent of the maximal effect is
//!
//! ```text
//! ECx = EC50 * (x / (100 - x))^(1/n)
//! ```
//!
//! so `EC50` recovers the half-maximal value (`x = 50`).
//!
//! ## Local slope
//!
//! Differentiating the absolute response `E(C) = Emax * f(C)`:
//!
//! ```text
//! E'(C) = Emax * n * EC50^n * C^(n-1) / (EC50^n + C^n)^2
//! ```
//!
//! At the inflection `C = EC50` this collapses to the clean value
//! `E'(EC50) = n * Emax / (4 * EC50)`, the steepest point of the curve.
//!
//! ## Dynamic range
//!
//! The ratio of the 90-percent to the 10-percent concentration depends
//! only on the Hill slope:
//!
//! ```text
//! EC90 / EC10 = 81^(1/n)
//! ```
//!
//! (`81 = 9 / (1/9)`). A slope of `n = 1` spans roughly two orders of
//! magnitude from 10 to 90 percent; steeper curves compress that span.

use crate::error::DoseResponseError;
use crate::hill::HillCurve;
use crate::Result;

impl HillCurve {
    /// Concentration producing `percent` percent of `Emax`:
    /// `ECx = EC50 * (x / (100 - x))^(1/n)`.
    ///
    /// `percent` must lie in the half-open interval `[0, 100)`; `100`
    /// percent (the full `Emax`) is only reached as `C -> infinity` and
    /// so has no finite concentration.
    ///
    /// # Errors
    ///
    /// Returns [`DoseResponseError::OutOfRange`] unless
    /// `0 <= percent < 100`, or [`DoseResponseError::NonFinite`] for a
    /// `NaN` / infinite argument.
    pub fn ec_percent(&self, percent: f64) -> Result<f64> {
        let p = DoseResponseError::finite("percent", percent)?;
        if !(0.0..100.0).contains(&p) {
            return Err(DoseResponseError::out_of_range("percent", p, 0.0, 100.0));
        }
        // Delegate to the fractional inverse: x percent is fraction
        // x/100. concentration_for_fraction handles the f == 0 case and
        // the (f/(1-f))^(1/n) algebra.
        self.concentration_for_fraction(p / 100.0)
    }

    /// Analytic slope of the absolute response curve at `concentration`:
    /// `E'(C) = Emax * n * EC50^n * C^(n-1) / (EC50^n + C^n)^2`.
    ///
    /// At `C = EC50` this equals `n * Emax / (4 * EC50)`.
    ///
    /// # Errors
    ///
    /// Returns [`DoseResponseError::Negative`] if `concentration < 0`, or
    /// [`DoseResponseError::NonFinite`] for a `NaN` / infinite argument.
    /// The slope at `C = 0` is `0` for `n > 1`, `Emax / EC50` for
    /// `n = 1`, and tends to `+infinity` for `n < 1`; the last case is
    /// reported as the IEEE `f64::INFINITY` value rather than an error.
    pub fn slope_at(&self, concentration: f64) -> Result<f64> {
        let c = DoseResponseError::non_negative("concentration", concentration)?;
        let n = self.hill_slope;
        if c == 0.0 {
            // Limit of E'(C) as C -> 0+: depends on the exponent (n-1).
            return Ok(if n > 1.0 {
                0.0
            } else if n == 1.0 {
                self.emax / self.ec50
            } else {
                f64::INFINITY
            });
        }
        // Work with t = (C/EC50)^n so the denominator (EC50^n + C^n)^2
        // becomes EC50^(2n) * (1 + t)^2 and the whole thing reduces to
        // E'(C) = Emax * n * t / (C * (1 + t)^2), which never forms the
        // overflow-prone bare powers EC50^n or C^n on their own.
        let t = (c / self.ec50).powf(n);
        if t.is_infinite() {
            // Far above EC50 the curve is flat again -> slope 0.
            return Ok(0.0);
        }
        let denom = (1.0 + t) * (1.0 + t);
        Ok(self.emax * n * t / (c * denom))
    }

    /// Ratio `EC90 / EC10 = 81^(1/n)` — the concentration span between 10
    /// and 90 percent of the maximal effect, a slope-only measure of how
    /// graded the response is.
    ///
    /// This never fails for a validly constructed [`HillCurve`] (the
    /// slope is guaranteed positive), so it returns a bare `f64`.
    pub fn dynamic_range_90_10(&self) -> f64 {
        81f64.powf(1.0 / self.hill_slope)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn approx(a: f64, b: f64) {
        assert!((a - b).abs() < EPS, "expected {a} ≈ {b}");
    }

    #[test]
    fn ec50_recovers_half_max() {
        let c = HillCurve::new(100.0, 7.0, 2.0).unwrap();
        approx(c.ec_percent(50.0).unwrap(), 7.0);
    }

    #[test]
    fn ec_percent_textbook_n1() {
        // n = 1, EC50 = 10: EC90 = 10 * (90/10)^1 = 90; EC10 = 10/9.
        let c = HillCurve::new(1.0, 10.0, 1.0).unwrap();
        approx(c.ec_percent(90.0).unwrap(), 90.0);
        approx(c.ec_percent(10.0).unwrap(), 10.0 / 9.0);
        // EC20 = 10 * (20/80) = 10 * 0.25 = 2.5.
        approx(c.ec_percent(20.0).unwrap(), 2.5);
    }

    #[test]
    fn ec_percent_consistent_with_response() {
        // Feeding ECx back into response() must yield x percent of Emax.
        let c = HillCurve::new(200.0, 4.0, 1.8).unwrap();
        for &p in &[5.0, 25.0, 50.0, 75.0, 95.0] {
            let conc = c.ec_percent(p).unwrap();
            approx(c.response(conc).unwrap(), 200.0 * p / 100.0);
        }
    }

    #[test]
    fn ec_percent_at_zero() {
        let c = HillCurve::new(100.0, 5.0, 2.0).unwrap();
        approx(c.ec_percent(0.0).unwrap(), 0.0);
    }

    #[test]
    fn ec_percent_rejects_out_of_range() {
        let c = HillCurve::new(100.0, 5.0, 1.0).unwrap();
        assert!(c.ec_percent(100.0).is_err()); // full Emax unreachable
        assert!(c.ec_percent(150.0).is_err());
        assert!(c.ec_percent(-1.0).is_err());
        assert!(c.ec_percent(f64::NAN).is_err());
    }

    #[test]
    fn slope_at_ec50_closed_form() {
        // GROUND TRUTH: E'(EC50) = n * Emax / (4 * EC50).
        // n = 2, Emax = 100, EC50 = 5 -> 2*100/(4*5) = 200/20 = 10.
        let c = HillCurve::new(100.0, 5.0, 2.0).unwrap();
        approx(c.slope_at(5.0).unwrap(), 10.0);

        // n = 1, Emax = 200, EC50 = 10 -> 1*200/(4*10) = 200/40 = 5.
        let c = HillCurve::new(200.0, 10.0, 1.0).unwrap();
        approx(c.slope_at(10.0).unwrap(), 5.0);
    }

    #[test]
    fn slope_matches_finite_difference() {
        // The analytic slope must match a central finite difference of
        // response() at an off-centre point.
        let c = HillCurve::new(120.0, 8.0, 2.5).unwrap();
        let x = 12.0;
        let h = 1e-5;
        let fd = (c.response(x + h).unwrap() - c.response(x - h).unwrap()) / (2.0 * h);
        let analytic = c.slope_at(x).unwrap();
        assert!(
            (fd - analytic).abs() < 1e-5,
            "fd {fd} vs analytic {analytic}"
        );
    }

    #[test]
    fn slope_at_zero_branches() {
        // n > 1: slope is 0 at the origin.
        let steep = HillCurve::new(100.0, 5.0, 3.0).unwrap();
        approx(steep.slope_at(0.0).unwrap(), 0.0);

        // n == 1: slope is Emax / EC50 at the origin (E = Emax*C/(EC50+C),
        // E'(0) = Emax/EC50). Emax = 100, EC50 = 5 -> 20.
        let flat = HillCurve::new(100.0, 5.0, 1.0).unwrap();
        approx(flat.slope_at(0.0).unwrap(), 20.0);

        // n < 1: slope diverges at the origin.
        let sub = HillCurve::new(100.0, 5.0, 0.5).unwrap();
        assert!(sub.slope_at(0.0).unwrap().is_infinite());
    }

    #[test]
    fn slope_is_nonnegative_and_peaks_at_inflection() {
        // For n > 1 the slope at EC50 must exceed the slope far away.
        let c = HillCurve::new(100.0, 10.0, 2.0).unwrap();
        let at_mid = c.slope_at(10.0).unwrap();
        let far = c.slope_at(1000.0).unwrap();
        let near0 = c.slope_at(0.01).unwrap();
        assert!(at_mid > far, "{at_mid} !> {far}");
        assert!(at_mid > near0, "{at_mid} !> {near0}");
        assert!(far >= 0.0 && near0 >= 0.0);
    }

    #[test]
    fn slope_rejects_bad_input() {
        let c = HillCurve::new(100.0, 5.0, 1.0).unwrap();
        assert!(c.slope_at(-1.0).is_err());
        assert!(c.slope_at(f64::NAN).is_err());
    }

    #[test]
    fn dynamic_range_closed_form() {
        // GROUND TRUTH: EC90/EC10 = 81^(1/n).
        // n = 1 -> 81. n = 2 -> sqrt(81) = 9. n = 4 -> 81^(1/4) = 3.
        approx(
            HillCurve::new(1.0, 1.0, 1.0).unwrap().dynamic_range_90_10(),
            81.0,
        );
        approx(
            HillCurve::new(1.0, 1.0, 2.0).unwrap().dynamic_range_90_10(),
            9.0,
        );
        approx(
            HillCurve::new(1.0, 1.0, 4.0).unwrap().dynamic_range_90_10(),
            3.0,
        );
    }

    #[test]
    fn dynamic_range_matches_ec_ratio() {
        // The closed form must equal the explicit EC90 / EC10 ratio.
        for &n in &[0.7, 1.0, 1.5, 3.0] {
            let c = HillCurve::new(50.0, 13.0, n).unwrap();
            let ec90 = c.ec_percent(90.0).unwrap();
            let ec10 = c.ec_percent(10.0).unwrap();
            approx(c.dynamic_range_90_10(), ec90 / ec10);
        }
    }
}
