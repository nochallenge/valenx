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

    /// The substrate concentration needed to reach a target velocity `v`,
    /// inverting the rate law: `S = Km * v / (Vmax - v)`.
    ///
    /// This is the design direction of Michaelis-Menten — "what `[S]`
    /// gives this velocity?". A finite answer only exists for
    /// `0 <= v < Vmax`, since `Vmax` is approached only as `S -> ∞`; a
    /// target at or above `Vmax` is rejected rather than returning an
    /// infinity. By construction `substrate_for_velocity(Vmax/2) == Km`
    /// and feeding the result back through [`velocity`](Self::velocity)
    /// recovers `v`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `v` is negative
    /// or non-finite, or if `v >= Vmax` (unreachable at any finite
    /// substrate; this also rejects every `v` for a dead `Vmax = 0`
    /// enzyme).
    pub fn substrate_for_velocity(&self, v: f64) -> Result<f64> {
        let v = require_non_negative("v", v)?;
        if v >= self.vmax {
            return Err(crate::error::KineticsError::out_of_domain(
                "v",
                v,
                "velocity must be below Vmax; Vmax is only approached as S -> infinity",
            ));
        }
        Ok(self.km * v / (self.vmax - v))
    }

    /// Time for the substrate to fall from an initial level `s0` to a lower
    /// level `s`, from the closed-form **integrated Michaelis-Menten
    /// equation**:
    ///
    /// ```text
    /// Vmax * t = Km * ln(s0 / s) + (s0 - s)
    /// ```
    ///
    /// the exact analytic integral of the depletion ODE `dS/dt = -v(S)`. In
    /// the first-order limit `Km >> S` it reduces to
    /// `t ≈ (Km/Vmax)*ln(s0/s)` (exponential decay); in the zero-order limit
    /// `Km << S` to `t ≈ (s0 - s)/Vmax` (constant-rate turnover at `Vmax`).
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `s0` or `s` is not
    /// finite and strictly positive, if `s > s0` (the substrate cannot rise
    /// during depletion), or if `Vmax = 0` (a dead enzyme never turns over).
    pub fn time_to_deplete(&self, s0: f64, s: f64) -> Result<f64> {
        let s0 = require_positive("s0", s0)?;
        let s = require_positive("s", s)?;
        if s > s0 {
            return Err(crate::error::KineticsError::out_of_domain(
                "s",
                s,
                "must not exceed the initial concentration s0",
            ));
        }
        if self.vmax <= 0.0 {
            return Err(crate::error::KineticsError::out_of_domain(
                "vmax",
                self.vmax,
                "must be > 0 for substrate depletion to occur",
            ));
        }
        Ok((self.km * (s0 / s).ln() + (s0 - s)) / self.vmax)
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

    /// Integrated MM: no elapsed time when the substrate has not moved, and a
    /// specific worked value. Vmax=2, Km=1, s0=10, s=5:
    /// t = (1*ln2 + 5)/2 = (0.693147 + 5)/2 = 2.846574.
    #[test]
    fn integrated_equation_zero_and_known_value() {
        let mm = MichaelisMenten::new(2.0, 1.0).expect("valid");
        assert!(mm.time_to_deplete(10.0, 10.0).expect("valid").abs() < EPS);
        let t = mm.time_to_deplete(10.0, 5.0).expect("valid");
        let want = (1.0 * 2.0_f64.ln() + 5.0) / 2.0;
        assert!((t - want).abs() < EPS, "t = {t}, want {want}");
    }

    /// First-order limit (Km >> S): t -> (Km/Vmax)*ln(s0/s).
    #[test]
    fn integrated_equation_first_order_limit() {
        let mm = MichaelisMenten::new(2.0, 1.0e6).expect("valid");
        let t = mm.time_to_deplete(10.0, 5.0).expect("valid");
        let approx = (1.0e6 / 2.0) * 2.0_f64.ln();
        assert!((t - approx).abs() / approx < 1e-4, "t {t} vs {approx}");
    }

    /// Zero-order limit (Km << S): t -> (s0 - s)/Vmax.
    #[test]
    fn integrated_equation_zero_order_limit() {
        let mm = MichaelisMenten::new(2.0, 1.0e-6).expect("valid");
        let t = mm.time_to_deplete(10.0, 5.0).expect("valid");
        let approx = (10.0 - 5.0) / 2.0;
        assert!((t - approx).abs() / approx < 1e-4, "t {t} vs {approx}");
    }

    /// The initial depletion slope matches the Michaelis-Menten initial rate:
    /// for a tiny depletion delta, t ≈ delta / v(s0).
    #[test]
    fn integrated_equation_initial_slope_is_the_initial_rate() {
        let mm = MichaelisMenten::new(2.0, 1.0).expect("valid");
        let s0 = 10.0;
        let delta = 1.0e-3;
        let t = mm.time_to_deplete(s0, s0 - delta).expect("valid");
        let v0 = mm.velocity(s0).expect("valid"); // 2*10/11 = 1.81818
        let expected = delta / v0;
        assert!(
            (t - expected).abs() / expected < 1e-3,
            "t {t} vs {expected}"
        );
    }

    #[test]
    fn integrated_equation_rejects_bad_inputs() {
        let mm = MichaelisMenten::new(2.0, 1.0).expect("valid");
        assert!(mm.time_to_deplete(10.0, 20.0).is_err()); // s > s0
        assert!(mm.time_to_deplete(0.0, 1.0).is_err()); // s0 <= 0
        assert!(mm.time_to_deplete(10.0, 0.0).is_err()); // s <= 0
                                                         // A dead enzyme (Vmax = 0) never depletes its substrate.
        let dead = MichaelisMenten::new(0.0, 1.0).expect("valid");
        assert!(dead.time_to_deplete(10.0, 5.0).is_err());
    }

    /// substrate_for_velocity inverts velocity exactly (the GOLD
    /// round-trip), across a sweep of target velocities below Vmax.
    #[test]
    fn substrate_for_velocity_inverts_velocity() {
        let mm = MichaelisMenten::new(7.3, 0.45).expect("valid");
        for frac in [0.0, 0.1, 0.25, 0.5, 0.8, 0.99] {
            let v_target = frac * 7.3;
            let s = mm.substrate_for_velocity(v_target).expect("v < Vmax");
            let v_back = mm.velocity(s).expect("valid s");
            assert!(
                (v_back - v_target).abs() < 1e-9 * 7.3_f64.max(1.0),
                "round-trip v {v_back} vs {v_target} at frac={frac}"
            );
        }
    }

    /// The defining property in inverse form: half Vmax needs exactly Km.
    #[test]
    fn half_vmax_needs_km_substrate() {
        let mm = MichaelisMenten::new(7.3, 0.45).expect("valid");
        let s = mm.substrate_for_velocity(7.3 / 2.0).expect("valid");
        assert!((s - 0.45).abs() < EPS, "S(Vmax/2) = {s}, want Km = 0.45");
        // Zero velocity needs zero substrate.
        assert!(mm.substrate_for_velocity(0.0).expect("valid").abs() < EPS);
    }

    /// Closed-form value and monotonic increase in the target velocity.
    #[test]
    fn substrate_for_velocity_closed_form_and_monotonic() {
        // Vmax=2, Km=1, v=1.5 -> S = 1*1.5/(2-1.5) = 3 (inverts known_value).
        let mm = MichaelisMenten::new(2.0, 1.0).expect("valid");
        let s = mm.substrate_for_velocity(1.5).expect("valid");
        assert!((s - 3.0).abs() < EPS, "S = {s}");
        let mut prev = -1.0;
        for k in 0..20 {
            let v = (k as f64) * 0.09 + 0.001; // up to ~1.71 < Vmax = 2
            let cur = mm.substrate_for_velocity(v).expect("valid");
            assert!(
                cur > prev,
                "S not increasing with v at v={v}: {cur} <= {prev}"
            );
            prev = cur;
        }
    }

    #[test]
    fn substrate_for_velocity_rejects_unreachable_and_bad() {
        let mm = MichaelisMenten::new(2.0, 1.0).expect("valid");
        assert!(mm.substrate_for_velocity(2.0).is_err()); // v == Vmax
        assert!(mm.substrate_for_velocity(2.5).is_err()); // v > Vmax
        assert!(mm.substrate_for_velocity(-0.1).is_err());
        assert!(mm.substrate_for_velocity(f64::NAN).is_err());
        // Dead enzyme: no velocity is reachable.
        let dead = MichaelisMenten::new(0.0, 1.0).expect("valid");
        assert!(dead.substrate_for_velocity(0.0).is_err());
    }
}
