//! The Michaelis-Menten rate law.
//!
//! For a single-substrate enzyme obeying the standard steady-state
//! (Briggs-Haldane) mechanism, the initial reaction velocity is
//!
//! ```text
//! v(S) = Vmax * S / (Km + S)
//! ```
//!
//! where `S` is the substrate concentration, `Vmax` the maximal velocity
//! (reached asymptotically as `S -> ∞`) and `Km` the Michaelis constant —
//! the substrate concentration at which `v = Vmax / 2`.
//!
//! ## Key textbook properties (all asserted in the unit tests)
//!
//! - `v(Km) = Vmax / 2` exactly.
//! - `v -> Vmax` as `S -> ∞` (saturation).
//! - `v` is strictly monotonically increasing in `S` for `Vmax > 0`.
//! - `v(0) = 0`.
//! - At low substrate (`S << Km`) the law linearises to the
//!   first-order form `v ≈ (Vmax / Km) * S`; the slope `Vmax / Km` is the
//!   specificity constant `kcat / Km` scaled by enzyme concentration.

use crate::error::{require_non_negative, require_positive, Result};
use serde::{Deserialize, Serialize};

/// A validated pair of Michaelis-Menten parameters `(Vmax, Km)`.
///
/// Construct with [`MichaelisMenten::new`], which enforces `Vmax >= 0`
/// and `Km > 0`. Once built, every method is infallible.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MichaelisMenten {
    /// Maximal velocity `Vmax` (same units as the returned velocity);
    /// the asymptote approached as `S -> ∞`. Non-negative.
    vmax: f64,
    /// Michaelis constant `Km` (same concentration units as `S`); the
    /// substrate level at which `v = Vmax / 2`. Strictly positive.
    km: f64,
}

impl MichaelisMenten {
    /// Build a validated parameter set.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `vmax` is
    /// negative or non-finite, or if `km` is not strictly positive (a
    /// zero or negative Michaelis constant has no physical meaning and
    /// would divide by zero at `S = 0`).
    pub fn new(vmax: f64, km: f64) -> Result<Self> {
        let vmax = require_non_negative("vmax", vmax)?;
        let km = require_positive("km", km)?;
        Ok(Self { vmax, km })
    }

    /// The maximal velocity `Vmax`.
    pub fn vmax(&self) -> f64 {
        self.vmax
    }

    /// The Michaelis constant `Km`.
    pub fn km(&self) -> f64 {
        self.km
    }

    /// Initial velocity `v(S) = Vmax * S / (Km + S)`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `s` is negative
    /// or non-finite. Because `Km > 0` is guaranteed by construction the
    /// denominator `Km + S` is always strictly positive, so no division
    /// hazard exists for any valid `s`.
    pub fn velocity(&self, s: f64) -> Result<f64> {
        let s = require_non_negative("s", s)?;
        Ok(self.vmax * s / (self.km + s))
    }

    /// Fraction of `Vmax` realised at substrate `s`, i.e. the saturation
    /// `v / Vmax = S / (Km + S)`. Independent of `Vmax`, this lies in
    /// `[0, 1)` and equals `0.5` exactly at `s = Km`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `s` is negative
    /// or non-finite.
    pub fn saturation(&self, s: f64) -> Result<f64> {
        let s = require_non_negative("s", s)?;
        Ok(s / (self.km + s))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference tolerance for the float comparisons below.
    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_bad_parameters() {
        assert!(MichaelisMenten::new(-1.0, 1.0).is_err());
        assert!(MichaelisMenten::new(1.0, 0.0).is_err());
        assert!(MichaelisMenten::new(1.0, -2.0).is_err());
        assert!(MichaelisMenten::new(f64::NAN, 1.0).is_err());
        assert!(MichaelisMenten::new(1.0, f64::INFINITY).is_err());
    }

    #[test]
    fn velocity_rejects_negative_substrate() {
        let mm = MichaelisMenten::new(2.0, 1.0).expect("valid");
        assert!(mm.velocity(-0.5).is_err());
        assert!(mm.velocity(f64::NAN).is_err());
    }

    /// v(Km) = Vmax / 2 — the defining property of the Michaelis constant.
    #[test]
    fn half_vmax_at_km() {
        let vmax = 7.3;
        let km = 0.45;
        let mm = MichaelisMenten::new(vmax, km).expect("valid");
        let v = mm.velocity(km).expect("valid substrate");
        assert!(
            (v - vmax / 2.0).abs() < EPS,
            "v(Km) = {v}, want {}",
            vmax / 2.0
        );
        let sat = mm.saturation(km).expect("valid");
        assert!((sat - 0.5).abs() < EPS, "saturation = {sat}");
    }

    /// v(0) = 0.
    #[test]
    fn zero_velocity_at_zero_substrate() {
        let mm = MichaelisMenten::new(5.0, 2.0).expect("valid");
        let v = mm.velocity(0.0).expect("valid");
        assert!(v.abs() < EPS, "v(0) = {v}");
    }

    /// v -> Vmax as S -> ∞ (checked at a very large but finite S, where
    /// the deficit Vmax - v = Vmax*Km/(Km+S) is tiny).
    #[test]
    fn approaches_vmax_at_high_substrate() {
        let vmax = 3.0;
        let km = 0.8;
        let mm = MichaelisMenten::new(vmax, km).expect("valid");
        let big = 1.0e9;
        let v = mm.velocity(big).expect("valid");
        // Analytic deficit at this S.
        let deficit = vmax * km / (km + big);
        assert!((vmax - v - deficit).abs() < EPS, "deficit mismatch: {v}");
        // And the deficit itself is genuinely small.
        assert!(deficit < 1e-8, "deficit not small: {deficit}");
        assert!(v < vmax, "v should stay below Vmax: {v}");
    }

    /// v is strictly increasing in S (sampled across decades).
    #[test]
    fn monotonically_increasing_in_substrate() {
        let mm = MichaelisMenten::new(4.2, 1.7).expect("valid");
        let mut prev = mm.velocity(0.0).expect("valid");
        for k in 0..200 {
            let s = (k as f64) * 0.5 + 0.01;
            let v = mm.velocity(s).expect("valid");
            assert!(v > prev, "not increasing at S={s}: {v} <= {prev}");
            prev = v;
        }
    }

    /// At S << Km the law linearises: v ≈ (Vmax/Km)*S.
    #[test]
    fn low_substrate_is_first_order() {
        let vmax = 10.0;
        let km = 5.0;
        let mm = MichaelisMenten::new(vmax, km).expect("valid");
        let s = km * 1e-4; // deep in the linear regime
        let v = mm.velocity(s).expect("valid");
        let linear = vmax / km * s;
        // Relative error of the linear approximation is ~ S/Km ≈ 1e-4.
        let rel = (v - linear).abs() / linear;
        assert!(rel < 1e-3, "linear approx rel error too large: {rel}");
    }

    /// A specific worked value: Vmax=2, Km=1, S=3 → v = 2*3/(1+3) = 1.5.
    #[test]
    fn known_value() {
        let mm = MichaelisMenten::new(2.0, 1.0).expect("valid");
        let v = mm.velocity(3.0).expect("valid");
        assert!((v - 1.5).abs() < EPS, "v = {v}");
    }

    #[test]
    fn parameters_round_trip_through_json() {
        let mm = MichaelisMenten::new(2.0, 1.0).expect("valid");
        let json = serde_json::to_string(&mm).expect("serialize");
        let back: MichaelisMenten = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(mm, back);
    }
}
