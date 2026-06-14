//! The Hill equation for cooperative binding / catalysis.
//!
//! ```text
//! v(S) = Vmax * S^n / (K^n + S^n)
//! ```
//!
//! `n` is the Hill coefficient (apparent cooperativity): `n = 1` is
//! non-cooperative and reduces *exactly* to Michaelis-Menten with
//! `Km = K`; `n > 1` is positively cooperative (sigmoidal `v`–`S` curve);
//! `0 < n < 1` is negatively cooperative. `K` is the half-saturation
//! constant — the substrate concentration giving `v = Vmax / 2` for any
//! `n` (note `K` equals `Km` only when `n = 1`).
//!
//! ## Key textbook properties (all asserted in the unit tests)
//!
//! - `v(K) = Vmax / 2` exactly, for every `n`.
//! - `n = 1` reproduces [`crate::MichaelisMenten`] with `Km = K` to
//!   floating-point tolerance.
//! - `v -> Vmax` as `S -> ∞`; `v(0) = 0`.
//! - `v` is strictly monotonically increasing in `S` for `Vmax > 0`.
//! - For `n > 1` the curve is sigmoidal: an inflection point exists at
//!   the positive substrate level `S* = K * ((n-1)/(n+1))^(1/n)`.

use crate::error::{require_non_negative, require_positive, Result};
use serde::{Deserialize, Serialize};

/// A validated set of Hill parameters `(Vmax, K, n)`.
///
/// Construct with [`Hill::new`], which enforces `Vmax >= 0`, `K > 0` and
/// `n > 0`. Once built, every method is infallible.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Hill {
    /// Maximal velocity `Vmax`; the asymptote as `S -> ∞`. Non-negative.
    vmax: f64,
    /// Half-saturation constant `K` (same units as `S`); the substrate
    /// level at which `v = Vmax / 2`. Strictly positive. Equals the
    /// Michaelis constant only when `n = 1`.
    k: f64,
    /// Hill coefficient `n` (apparent cooperativity). Strictly positive.
    n: f64,
}

impl Hill {
    /// Build a validated Hill parameter set.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `vmax` is
    /// negative or non-finite, if `k` is not strictly positive, or if `n`
    /// is not strictly positive (a zero or negative Hill coefficient is
    /// not physically meaningful and would invert the saturation curve).
    pub fn new(vmax: f64, k: f64, n: f64) -> Result<Self> {
        let vmax = require_non_negative("vmax", vmax)?;
        let k = require_positive("k", k)?;
        let n = require_positive("n", n)?;
        Ok(Self { vmax, k, n })
    }

    /// The maximal velocity `Vmax`.
    pub fn vmax(&self) -> f64 {
        self.vmax
    }

    /// The half-saturation constant `K`.
    pub fn k(&self) -> f64 {
        self.k
    }

    /// The Hill coefficient `n`.
    pub fn n(&self) -> f64 {
        self.n
    }

    /// Initial velocity `v(S) = Vmax * S^n / (K^n + S^n)`.
    ///
    /// Evaluated in the algebraically equivalent ratio form
    /// `Vmax * x / (1 + x)` with `x = (S / K)^n`, which avoids overflow of
    /// `S^n` and `K^n` separately for large exponents and keeps the
    /// `S = 0` case (`x = 0`, `v = 0`) exact.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `s` is negative
    /// or non-finite. Because `K > 0` by construction the ratio `S / K`
    /// is always well-defined.
    pub fn velocity(&self, s: f64) -> Result<f64> {
        let s = require_non_negative("s", s)?;
        if s == 0.0 {
            return Ok(0.0);
        }
        let x = (s / self.k).powf(self.n);
        // x is finite and >= 0 here; for very large x, x/(1+x) -> 1.
        Ok(self.vmax * x / (1.0 + x))
    }

    /// Fractional saturation `v / Vmax = S^n / (K^n + S^n)`, independent
    /// of `Vmax`; lies in `[0, 1)` and equals `0.5` exactly at `s = K`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `s` is negative
    /// or non-finite.
    pub fn saturation(&self, s: f64) -> Result<f64> {
        let s = require_non_negative("s", s)?;
        if s == 0.0 {
            return Ok(0.0);
        }
        let x = (s / self.k).powf(self.n);
        Ok(x / (1.0 + x))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::michaelis_menten::MichaelisMenten;

    /// Absolute-difference tolerance for the float comparisons below.
    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_bad_parameters() {
        assert!(Hill::new(-1.0, 1.0, 1.0).is_err());
        assert!(Hill::new(1.0, 0.0, 1.0).is_err());
        assert!(Hill::new(1.0, 1.0, 0.0).is_err());
        assert!(Hill::new(1.0, 1.0, -2.0).is_err());
        assert!(Hill::new(1.0, f64::NAN, 1.0).is_err());
    }

    /// v(K) = Vmax / 2 for several Hill coefficients.
    #[test]
    fn half_vmax_at_k_for_all_n() {
        let vmax = 6.5;
        let k = 1.3;
        for &n in &[0.5, 1.0, 2.0, 4.0, 8.0] {
            let h = Hill::new(vmax, k, n).expect("valid");
            let v = h.velocity(k).expect("valid");
            assert!(
                (v - vmax / 2.0).abs() < EPS,
                "n={n}: v(K) = {v}, want {}",
                vmax / 2.0
            );
            let sat = h.saturation(k).expect("valid");
            assert!((sat - 0.5).abs() < EPS, "n={n}: saturation = {sat}");
        }
    }

    /// n = 1 reduces exactly to Michaelis-Menten with Km = K.
    #[test]
    fn n_equals_one_reduces_to_michaelis_menten() {
        let vmax = 3.7;
        let k = 0.9;
        let h = Hill::new(vmax, k, 1.0).expect("valid");
        let mm = MichaelisMenten::new(vmax, k).expect("valid");
        for &s in &[0.0, 0.05, 0.3, 0.9, 2.0, 10.0, 250.0, 5000.0] {
            let vh = h.velocity(s).expect("valid");
            let vm = mm.velocity(s).expect("valid");
            assert!((vh - vm).abs() < EPS, "S={s}: hill {vh} vs mm {vm}");
        }
    }

    /// v(0) = 0 for every n.
    #[test]
    fn zero_velocity_at_zero_substrate() {
        for &n in &[0.5, 1.0, 3.0] {
            let h = Hill::new(5.0, 2.0, n).expect("valid");
            let v = h.velocity(0.0).expect("valid");
            assert!(v.abs() < EPS, "n={n}: v(0) = {v}");
        }
    }

    /// v -> Vmax as S -> ∞ for every n. The exact deficit is
    /// `Vmax - v = Vmax / (1 + (S/K)^n)`; we (a) match it analytically at
    /// a large finite S and (b) confirm it strictly shrinks toward zero as
    /// S grows tenfold. Negative cooperativity (n < 1) approaches Vmax far
    /// more slowly, so a fixed absolute threshold would be wrong — the
    /// limit is verified through the shrinking deficit instead.
    #[test]
    fn approaches_vmax_at_high_substrate() {
        let vmax = 2.0;
        let k = 1.0;
        // Moderate S keeps the deficit representable as a normal f64 for
        // every n (at very large x = (S/K)^n the deficit underflows and v
        // rounds to exactly Vmax, which is a representation artifact, not
        // the mathematics). At these S the deficit is well above ulp(Vmax).
        for &n in &[0.5, 1.0, 2.0, 4.0] {
            let h = Hill::new(vmax, k, n).expect("valid");
            let s1 = 1.0e2;
            let s2 = 1.0e3;
            let v1 = h.velocity(s1).expect("valid");
            let v2 = h.velocity(s2).expect("valid");
            // (a) exact analytic deficit at S = s1.
            let deficit1 = vmax / (1.0 + (s1 / k).powf(n));
            assert!(
                (vmax - v1 - deficit1).abs() < EPS,
                "n={n}: deficit mismatch: {v1}"
            );
            // (b) the velocity stays below Vmax and the deficit shrinks
            //     strictly toward zero as S grows tenfold.
            assert!(v1 < vmax, "n={n}: v should stay below Vmax: {v1}");
            assert!(v2 > v1, "n={n}: deficit should shrink: {v2} <= {v1}");
            assert!(
                (vmax - v2) < (vmax - v1),
                "n={n}: deficit not shrinking toward 0"
            );
        }
    }

    /// v is strictly increasing in S for every n (sampled across decades).
    #[test]
    fn monotonically_increasing_in_substrate() {
        for &n in &[0.5, 1.0, 2.5, 5.0] {
            let h = Hill::new(4.2, 1.7, n).expect("valid");
            let mut prev = h.velocity(0.0).expect("valid");
            for k in 0..300 {
                let s = (k as f64) * 0.1 + 1e-4;
                let v = h.velocity(s).expect("valid");
                assert!(v > prev, "n={n}: not increasing at S={s}: {v} <= {prev}");
                prev = v;
            }
        }
    }

    /// For n > 1 the Hill curve is sigmoidal: it is convex (accelerating)
    /// below the analytic inflection S* and concave above it. We verify
    /// the second derivative changes sign across S* by comparing the
    /// discrete curvature of the saturation curve on each side.
    #[test]
    fn sigmoidal_inflection_for_cooperative_n() {
        let n = 3.0;
        let k = 1.0;
        let h = Hill::new(1.0, k, n).expect("valid");
        // Analytic inflection of S^n/(K^n+S^n): S* = K*((n-1)/(n+1))^(1/n).
        let s_star = k * ((n - 1.0) / (n + 1.0)).powf(1.0 / n);
        assert!(s_star > 0.0, "inflection should be positive: {s_star}");

        // Discrete second difference of saturation, well inside each side.
        let curv = |s: f64| -> f64 {
            let d = 1e-3;
            let a = h.saturation(s - d).expect("valid");
            let b = h.saturation(s).expect("valid");
            let c = h.saturation(s + d).expect("valid");
            (c - 2.0 * b + a) / (d * d)
        };
        let below = curv(s_star * 0.5);
        let above = curv(s_star * 2.0);
        assert!(below > 0.0, "should be convex below S*: {below}");
        assert!(above < 0.0, "should be concave above S*: {above}");
    }

    /// A specific worked value: Vmax=10, K=2, n=2, S=2 → S^n=K^n=4,
    /// v = 10*4/(4+4) = 5 (= Vmax/2 at S=K, cross-check).
    #[test]
    fn known_value() {
        let h = Hill::new(10.0, 2.0, 2.0).expect("valid");
        let v = h.velocity(2.0).expect("valid");
        assert!((v - 5.0).abs() < EPS, "v = {v}");
        // And at S=4 (=2K): x=(4/2)^2=4, v=10*4/5=8.
        let v2 = h.velocity(4.0).expect("valid");
        assert!((v2 - 8.0).abs() < EPS, "v(4) = {v2}");
    }

    #[test]
    fn velocity_rejects_negative_substrate() {
        let h = Hill::new(1.0, 1.0, 2.0).expect("valid");
        assert!(h.velocity(-1.0).is_err());
        assert!(h.velocity(f64::INFINITY).is_err());
    }

    #[test]
    fn parameters_round_trip_through_json() {
        let h = Hill::new(10.0, 2.0, 2.0).expect("valid");
        let json = serde_json::to_string(&h).expect("serialize");
        let back: Hill = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(h, back);
    }
}
