//! The Hill equation and its inverse.
//!
//! The Hill (Hill-Langmuir) equation models the response `E` of a system
//! to a drug concentration `C`:
//!
//! ```text
//! E(C) = Emax * C^n / (EC50^n + C^n)
//! ```
//!
//! where `Emax` is the maximal achievable response (the upper asymptote),
//! `EC50` is the concentration producing the half-maximal response, and
//! `n` is the Hill slope (the cooperativity / steepness coefficient).
//!
//! Dividing through by `Emax` gives the dimensionless **fractional
//! response** (a.k.a. fractional occupancy when `Emax = 1`):
//!
//! ```text
//! f(C) = C^n / (EC50^n + C^n)  in  [0, 1)
//! ```
//!
//! Both forms are monotonically increasing in `C`, equal `Emax/2`
//! (respectively `1/2`) exactly at `C = EC50`, tend to `0` as `C -> 0`,
//! and approach `Emax` (respectively `1`) as `C -> infinity`. These four
//! ground-truth facts are what the unit tests check.

use crate::error::DoseResponseError;
use crate::Result;

/// Parameters of a Hill dose-response curve.
///
/// Construct with [`HillCurve::new`], which validates that every field
/// is finite, that `emax`, `ec50`, and `hill_slope` are strictly
/// positive (a concentration and a maximal effect cannot be zero or
/// negative, and a zero / negative slope is not a sigmoid). The fields
/// are public for read-only inspection but the validating constructor is
/// the only way to build one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HillCurve {
    /// Maximal response (upper asymptote), `Emax > 0`.
    pub emax: f64,
    /// Half-maximal-response concentration, `EC50 > 0`, same units as
    /// the concentration passed to [`response`](HillCurve::response).
    pub ec50: f64,
    /// Hill slope `n > 0` (dimensionless cooperativity coefficient).
    pub hill_slope: f64,
}

impl HillCurve {
    /// Build a validated curve.
    ///
    /// # Errors
    ///
    /// Returns [`DoseResponseError::NonFinite`] if any argument is `NaN`
    /// or infinite, or [`DoseResponseError::NonPositive`] if `emax`,
    /// `ec50`, or `hill_slope` is not strictly positive.
    pub fn new(emax: f64, ec50: f64, hill_slope: f64) -> Result<Self> {
        let emax = DoseResponseError::positive("emax", emax)?;
        let ec50 = DoseResponseError::positive("ec50", ec50)?;
        let hill_slope = DoseResponseError::positive("hill_slope", hill_slope)?;
        Ok(HillCurve {
            emax,
            ec50,
            hill_slope,
        })
    }

    /// Absolute response `E(C) = Emax * C^n / (EC50^n + C^n)`.
    ///
    /// At `C = 0` this is exactly `0`; at `C = EC50` it is exactly
    /// `Emax / 2`; as `C -> infinity` it approaches `Emax`.
    ///
    /// # Errors
    ///
    /// Returns [`DoseResponseError::Negative`] if `concentration < 0`
    /// (negative drug concentrations are unphysical), or
    /// [`DoseResponseError::NonFinite`] if it is `NaN` / infinite.
    pub fn response(&self, concentration: f64) -> Result<f64> {
        Ok(self.emax * self.fraction(concentration)?)
    }

    /// Dimensionless fractional response
    /// `f(C) = C^n / (EC50^n + C^n)` in `[0, 1)`.
    ///
    /// This is [`response`](HillCurve::response) divided by `Emax`, i.e.
    /// the fraction of the maximal effect produced at concentration `C`.
    ///
    /// # Errors
    ///
    /// Same validation as [`response`](HillCurve::response).
    pub fn fraction(&self, concentration: f64) -> Result<f64> {
        let c = DoseResponseError::non_negative("concentration", concentration)?;
        if c == 0.0 {
            return Ok(0.0);
        }
        // Compute via the ratio r = (C / EC50)^n so that f = r / (1 + r).
        // This is numerically gentler than forming C^n and EC50^n
        // separately (which overflow for large n) and is algebraically
        // identical: C^n / (EC50^n + C^n) = r / (1 + r).
        let r = (c / self.ec50).powf(self.hill_slope);
        if r.is_infinite() {
            // Saturating limit: an enormous C/EC50 ratio rounds to the
            // upper asymptote fraction of 1.0.
            return Ok(1.0);
        }
        Ok(r / (1.0 + r))
    }

    /// Inverse: the concentration producing a given absolute `response`.
    ///
    /// Solving `E = Emax * C^n / (EC50^n + C^n)` for `C` gives
    /// `C = EC50 * (E / (Emax - E))^(1/n)`.
    ///
    /// # Errors
    ///
    /// Returns [`DoseResponseError::OutOfRange`] unless
    /// `0 <= response < Emax` (the response `Emax` itself is only reached
    /// in the limit `C -> infinity`, so it has no finite inverse), and
    /// [`DoseResponseError::NonFinite`] for a `NaN` / infinite argument.
    pub fn concentration_for_response(&self, response: f64) -> Result<f64> {
        let e = DoseResponseError::finite("response", response)?;
        if !(0.0..self.emax).contains(&e) {
            return Err(DoseResponseError::out_of_range(
                "response", e, 0.0, self.emax,
            ));
        }
        self.concentration_for_fraction(e / self.emax)
    }

    /// Inverse expressed via the fractional response: the concentration
    /// producing fraction `f` of `Emax`.
    ///
    /// `C = EC50 * (f / (1 - f))^(1/n)`.
    ///
    /// # Errors
    ///
    /// Returns [`DoseResponseError::OutOfRange`] unless `0 <= f < 1`, and
    /// [`DoseResponseError::NonFinite`] for a `NaN` / infinite argument.
    pub fn concentration_for_fraction(&self, fraction: f64) -> Result<f64> {
        let f = DoseResponseError::finite("fraction", fraction)?;
        if !(0.0..1.0).contains(&f) {
            return Err(DoseResponseError::out_of_range("fraction", f, 0.0, 1.0));
        }
        if f == 0.0 {
            return Ok(0.0);
        }
        let ratio = f / (1.0 - f);
        Ok(self.ec50 * ratio.powf(1.0 / self.hill_slope))
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
    fn constructor_rejects_bad_params() {
        assert!(HillCurve::new(0.0, 1.0, 1.0).is_err()); // emax = 0
        assert!(HillCurve::new(-1.0, 1.0, 1.0).is_err()); // emax < 0
        assert!(HillCurve::new(100.0, 0.0, 1.0).is_err()); // ec50 = 0
        assert!(HillCurve::new(100.0, -2.0, 1.0).is_err()); // ec50 < 0
        assert!(HillCurve::new(100.0, 1.0, 0.0).is_err()); // slope = 0
        assert!(HillCurve::new(100.0, 1.0, -1.0).is_err()); // slope < 0
        assert!(HillCurve::new(f64::NAN, 1.0, 1.0).is_err()); // non-finite
        assert!(HillCurve::new(100.0, f64::INFINITY, 1.0).is_err());
    }

    #[test]
    fn half_max_at_ec50_n1() {
        // GROUND TRUTH: C = EC50  =>  response = Emax / 2, for any n.
        let c = HillCurve::new(100.0, 5.0, 1.0).unwrap();
        approx(c.response(5.0).unwrap(), 50.0);
        approx(c.fraction(5.0).unwrap(), 0.5);
    }

    #[test]
    fn half_max_at_ec50_independent_of_slope() {
        // The C = EC50 -> Emax/2 identity holds for every Hill slope.
        for &n in &[0.5, 1.0, 2.0, 3.7, 10.0] {
            let c = HillCurve::new(80.0, 12.0, n).unwrap();
            approx(c.response(12.0).unwrap(), 40.0);
            approx(c.fraction(12.0).unwrap(), 0.5);
        }
    }

    #[test]
    fn zero_concentration_gives_zero() {
        // GROUND TRUTH: C -> 0  =>  E -> 0.
        let c = HillCurve::new(100.0, 5.0, 2.0).unwrap();
        approx(c.response(0.0).unwrap(), 0.0);
        approx(c.fraction(0.0).unwrap(), 0.0);
    }

    #[test]
    fn saturates_to_emax_for_large_c() {
        // GROUND TRUTH: C >> EC50  =>  E -> Emax. With C = 1e6 * EC50 and
        // n = 1, f = 1e6 / (1 + 1e6) ≈ 0.999999.
        let c = HillCurve::new(100.0, 1.0, 1.0).unwrap();
        let r = c.response(1.0e6).unwrap();
        assert!(r > 99.9999 && r < 100.0, "got {r}");
        // Hand value: f = 1e6/(1e6+1).
        approx(c.fraction(1.0e6).unwrap(), 1.0e6 / (1.0e6 + 1.0));
    }

    #[test]
    fn n1_closed_form_matches_textbook() {
        // For n = 1 the Hill equation is the rectangular hyperbola
        // E = Emax*C/(EC50+C). Hand-work EC50 = 10, Emax = 200, C = 30:
        // E = 200*30/(10+30) = 6000/40 = 150.
        let c = HillCurve::new(200.0, 10.0, 1.0).unwrap();
        approx(c.response(30.0).unwrap(), 150.0);
        // C = 10/3 gives E = 200*(10/3)/(10 + 10/3) = (2000/3)/(40/3) = 50.
        approx(c.response(10.0 / 3.0).unwrap(), 50.0);
    }

    #[test]
    fn n2_closed_form_matches_textbook() {
        // n = 2, EC50 = 4, Emax = 1, C = 8:
        // f = 8^2 / (4^2 + 8^2) = 64 / (16 + 64) = 64/80 = 0.8.
        let c = HillCurve::new(1.0, 4.0, 2.0).unwrap();
        approx(c.fraction(8.0).unwrap(), 0.8);
        // C = 2: f = 4 / (16 + 4) = 4/20 = 0.2.
        approx(c.fraction(2.0).unwrap(), 0.2);
    }

    #[test]
    fn monotonic_increasing_in_concentration() {
        let c = HillCurve::new(100.0, 7.0, 1.6).unwrap();
        let mut prev = c.response(0.0).unwrap();
        let mut x = 0.1;
        while x < 1000.0 {
            let cur = c.response(x).unwrap();
            assert!(cur > prev, "not increasing at C={x}: {cur} <= {prev}");
            prev = cur;
            x *= 1.5;
        }
    }

    #[test]
    fn inverse_round_trips() {
        // concentration_for_response(response(C)) == C.
        let c = HillCurve::new(120.0, 9.0, 2.3).unwrap();
        for &conc in &[0.5, 3.0, 9.0, 25.0, 200.0] {
            let e = c.response(conc).unwrap();
            let back = c.concentration_for_response(e).unwrap();
            assert!((back - conc).abs() < 1e-6, "C={conc} -> {back}");
        }
    }

    #[test]
    fn inverse_of_half_max_is_ec50() {
        // GROUND TRUTH: the concentration giving Emax/2 is exactly EC50.
        let c = HillCurve::new(50.0, 6.0, 3.0).unwrap();
        approx(c.concentration_for_response(25.0).unwrap(), 6.0);
        approx(c.concentration_for_fraction(0.5).unwrap(), 6.0);
    }

    #[test]
    fn inverse_rejects_out_of_range() {
        let c = HillCurve::new(100.0, 5.0, 1.0).unwrap();
        // response == Emax has no finite solution (open upper bound).
        assert!(c.concentration_for_response(100.0).is_err());
        assert!(c.concentration_for_response(150.0).is_err());
        assert!(c.concentration_for_response(-1.0).is_err());
        assert!(c.concentration_for_fraction(1.0).is_err());
        assert!(c.concentration_for_fraction(1.5).is_err());
        assert!(c.concentration_for_fraction(-0.1).is_err());
        assert!(c.concentration_for_fraction(f64::NAN).is_err());
    }

    #[test]
    fn response_rejects_negative_concentration() {
        let c = HillCurve::new(100.0, 5.0, 1.0).unwrap();
        assert!(c.response(-1.0).is_err());
        assert!(c.fraction(-0.001).is_err());
        assert!(c.response(f64::NAN).is_err());
    }

    #[test]
    fn high_slope_does_not_overflow() {
        // n = 50 with C >> EC50 would overflow C^n if formed directly;
        // the r/(1+r) formulation saturates cleanly to 1.0 instead.
        let c = HillCurve::new(1.0, 1.0, 50.0).unwrap();
        let f = c.fraction(100.0).unwrap();
        assert!((f - 1.0).abs() < EPS, "got {f}");
        // And below EC50 it collapses toward 0.
        let f_lo = c.fraction(0.01).unwrap();
        assert!(f_lo < EPS, "got {f_lo}");
    }
}
