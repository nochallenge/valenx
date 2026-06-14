//! Reversible-inhibition rate laws.
//!
//! Each inhibition mode is the Michaelis-Menten law evaluated with
//! *apparent* parameters obtained by scaling `Km` and / or `Vmax` by the
//! dimensionless factors
//!
//! ```text
//! alpha   = 1 + I / Ki    (binds free enzyme  E  — competitive arm)
//! alpha'  = 1 + I / Ki'   (binds the ES complex — uncompetitive arm)
//! ```
//!
//! With those factors the general (mixed) law is
//!
//! ```text
//! v = Vmax * S / (alpha * Km + alpha' * S)
//!   = (Vmax / alpha') * S / ((alpha / alpha') * Km + S)
//! ```
//!
//! so the apparent maximal velocity is `Vmax / alpha'` and the apparent
//! Michaelis constant is `Km * alpha / alpha'`. The three named special
//! cases drop out:
//!
//! | mode            | apparent Km          | apparent Vmax     |
//! |-----------------|----------------------|-------------------|
//! | competitive     | `Km * (1 + I/Ki)`    | `Vmax` (unchanged)|
//! | noncompetitive  | `Km` (unchanged)     | `Vmax / (1 + I/Ki)` |
//! | uncompetitive   | `Km / (1 + I/Ki')`   | `Vmax / (1 + I/Ki')`|
//!
//! (Classical "pure" noncompetitive inhibition is the symmetric mixed
//! case `Ki == Ki'`, so `alpha == alpha'`: the apparent `Km` is then
//! unchanged and only `Vmax` falls.)
//!
//! The diagnostic signatures the unit tests pin down:
//!
//! - a **competitive** inhibitor raises the apparent `Km` but leaves
//!   `Vmax` untouched;
//! - a **noncompetitive** inhibitor lowers the apparent `Vmax` but leaves
//!   `Km` untouched;
//! - an **uncompetitive** inhibitor lowers *both* by the same factor,
//!   leaving their ratio `Vmax/Km` unchanged.

use crate::error::{require_non_negative, require_positive, Result};
use crate::michaelis_menten::MichaelisMenten;
use serde::{Deserialize, Serialize};

/// Competitive inhibition: the inhibitor competes with substrate for the
/// active site, raising the apparent `Km` to `Km * (1 + I/Ki)` while
/// leaving `Vmax` unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Competitive {
    /// Uninhibited Michaelis-Menten parameters.
    base: MichaelisMenten,
    /// Inhibition constant `Ki` (dissociation constant of the
    /// enzyme-inhibitor complex `E·I`), same units as the inhibitor
    /// concentration. Strictly positive.
    ki: f64,
}

impl Competitive {
    /// Build a validated competitive-inhibition model from the
    /// uninhibited parameters and the inhibition constant `ki`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `ki` is not
    /// strictly positive (or non-finite).
    pub fn new(base: MichaelisMenten, ki: f64) -> Result<Self> {
        let ki = require_positive("ki", ki)?;
        Ok(Self { base, ki })
    }

    /// The inhibition constant `Ki`.
    pub fn ki(&self) -> f64 {
        self.ki
    }

    /// Apparent Michaelis constant at inhibitor concentration `i`:
    /// `Km_app = Km * (1 + I / Ki)`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `i` is negative
    /// or non-finite.
    pub fn apparent_km(&self, i: f64) -> Result<f64> {
        let i = require_non_negative("i", i)?;
        Ok(self.base.km() * (1.0 + i / self.ki))
    }

    /// Apparent maximal velocity. For competitive inhibition `Vmax` is
    /// unchanged, so this always equals the uninhibited `Vmax`; the `i`
    /// argument is still validated for symmetry with the other modes.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `i` is negative
    /// or non-finite.
    pub fn apparent_vmax(&self, i: f64) -> Result<f64> {
        let _ = require_non_negative("i", i)?;
        Ok(self.base.vmax())
    }

    /// Initial velocity at substrate `s` and inhibitor `i`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if either `s` or
    /// `i` is negative or non-finite.
    pub fn velocity(&self, s: f64, i: f64) -> Result<f64> {
        let s = require_non_negative("s", s)?;
        let km_app = self.apparent_km(i)?;
        Ok(self.base.vmax() * s / (km_app + s))
    }
}

/// Noncompetitive (here, "pure" / symmetric mixed) inhibition: the
/// inhibitor binds `E` and `ES` with equal affinity, lowering the
/// apparent `Vmax` to `Vmax / (1 + I/Ki)` while leaving `Km` unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Noncompetitive {
    /// Uninhibited Michaelis-Menten parameters.
    base: MichaelisMenten,
    /// Inhibition constant `Ki` (here the common dissociation constant of
    /// both `E·I` and `ES·I`). Strictly positive.
    ki: f64,
}

impl Noncompetitive {
    /// Build a validated noncompetitive-inhibition model.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `ki` is not
    /// strictly positive (or non-finite).
    pub fn new(base: MichaelisMenten, ki: f64) -> Result<Self> {
        let ki = require_positive("ki", ki)?;
        Ok(Self { base, ki })
    }

    /// The inhibition constant `Ki`.
    pub fn ki(&self) -> f64 {
        self.ki
    }

    /// Apparent maximal velocity at inhibitor `i`:
    /// `Vmax_app = Vmax / (1 + I / Ki)`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `i` is negative
    /// or non-finite.
    pub fn apparent_vmax(&self, i: f64) -> Result<f64> {
        let i = require_non_negative("i", i)?;
        Ok(self.base.vmax() / (1.0 + i / self.ki))
    }

    /// Apparent Michaelis constant. For pure noncompetitive inhibition
    /// `Km` is unchanged, so this always equals the uninhibited `Km`; the
    /// `i` argument is still validated for symmetry with the other modes.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `i` is negative
    /// or non-finite.
    pub fn apparent_km(&self, i: f64) -> Result<f64> {
        let _ = require_non_negative("i", i)?;
        Ok(self.base.km())
    }

    /// Initial velocity at substrate `s` and inhibitor `i`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if either `s` or
    /// `i` is negative or non-finite.
    pub fn velocity(&self, s: f64, i: f64) -> Result<f64> {
        let s = require_non_negative("s", s)?;
        let vmax_app = self.apparent_vmax(i)?;
        Ok(vmax_app * s / (self.base.km() + s))
    }
}

/// Uncompetitive inhibition: the inhibitor binds only the `ES` complex,
/// lowering *both* apparent `Vmax` and apparent `Km` by the same factor
/// `1 + I/Ki'`, so the ratio `Vmax/Km` (and hence the low-substrate
/// slope) is unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Uncompetitive {
    /// Uninhibited Michaelis-Menten parameters.
    base: MichaelisMenten,
    /// Inhibition constant `Ki'` (dissociation constant of the `ES·I`
    /// complex). Strictly positive.
    ki_prime: f64,
}

impl Uncompetitive {
    /// Build a validated uncompetitive-inhibition model from the
    /// uninhibited parameters and the `ES·I` inhibition constant
    /// `ki_prime`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `ki_prime` is
    /// not strictly positive (or non-finite).
    pub fn new(base: MichaelisMenten, ki_prime: f64) -> Result<Self> {
        let ki_prime = require_positive("ki_prime", ki_prime)?;
        Ok(Self { base, ki_prime })
    }

    /// The `ES·I` inhibition constant `Ki'`.
    pub fn ki_prime(&self) -> f64 {
        self.ki_prime
    }

    /// Apparent maximal velocity `Vmax / (1 + I / Ki')`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `i` is negative
    /// or non-finite.
    pub fn apparent_vmax(&self, i: f64) -> Result<f64> {
        let i = require_non_negative("i", i)?;
        Ok(self.base.vmax() / (1.0 + i / self.ki_prime))
    }

    /// Apparent Michaelis constant `Km / (1 + I / Ki')`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `i` is negative
    /// or non-finite.
    pub fn apparent_km(&self, i: f64) -> Result<f64> {
        let i = require_non_negative("i", i)?;
        Ok(self.base.km() / (1.0 + i / self.ki_prime))
    }

    /// Initial velocity at substrate `s` and inhibitor `i`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if either `s` or
    /// `i` is negative or non-finite.
    pub fn velocity(&self, s: f64, i: f64) -> Result<f64> {
        let s = require_non_negative("s", s)?;
        let vmax_app = self.apparent_vmax(i)?;
        let km_app = self.apparent_km(i)?;
        Ok(vmax_app * s / (km_app + s))
    }
}

/// General mixed inhibition: independent factors on the competitive
/// (free-enzyme) and uncompetitive (`ES`-complex) arms.
///
/// With `alpha = 1 + I/Ki` and `alpha' = 1 + I/Ki'`,
///
/// ```text
/// v = Vmax * S / (alpha * Km + alpha' * S)
/// ```
///
/// so `Vmax_app = Vmax / alpha'` and `Km_app = Km * alpha / alpha'`. The
/// three named modes are special cases (competitive: `Ki' -> ∞`;
/// uncompetitive: `Ki -> ∞`; pure noncompetitive: `Ki == Ki'`).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Mixed {
    /// Uninhibited Michaelis-Menten parameters.
    base: MichaelisMenten,
    /// Free-enzyme (`E·I`) inhibition constant `Ki`. Strictly positive.
    ki: f64,
    /// `ES`-complex (`ES·I`) inhibition constant `Ki'`. Strictly
    /// positive.
    ki_prime: f64,
}

impl Mixed {
    /// Build a validated mixed-inhibition model.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if either `ki` or
    /// `ki_prime` is not strictly positive (or non-finite).
    pub fn new(base: MichaelisMenten, ki: f64, ki_prime: f64) -> Result<Self> {
        let ki = require_positive("ki", ki)?;
        let ki_prime = require_positive("ki_prime", ki_prime)?;
        Ok(Self { base, ki, ki_prime })
    }

    /// The free-enzyme inhibition constant `Ki`.
    pub fn ki(&self) -> f64 {
        self.ki
    }

    /// The `ES`-complex inhibition constant `Ki'`.
    pub fn ki_prime(&self) -> f64 {
        self.ki_prime
    }

    /// Apparent maximal velocity `Vmax / alpha'` where
    /// `alpha' = 1 + I/Ki'`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `i` is negative
    /// or non-finite.
    pub fn apparent_vmax(&self, i: f64) -> Result<f64> {
        let i = require_non_negative("i", i)?;
        let alpha_prime = 1.0 + i / self.ki_prime;
        Ok(self.base.vmax() / alpha_prime)
    }

    /// Apparent Michaelis constant `Km * alpha / alpha'` where
    /// `alpha = 1 + I/Ki` and `alpha' = 1 + I/Ki'`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if `i` is negative
    /// or non-finite.
    pub fn apparent_km(&self, i: f64) -> Result<f64> {
        let i = require_non_negative("i", i)?;
        let alpha = 1.0 + i / self.ki;
        let alpha_prime = 1.0 + i / self.ki_prime;
        Ok(self.base.km() * alpha / alpha_prime)
    }

    /// Initial velocity `v = Vmax * S / (alpha*Km + alpha'*S)`.
    ///
    /// # Errors
    ///
    /// Returns [`KineticsError`](crate::KineticsError) if either `s` or
    /// `i` is negative or non-finite.
    pub fn velocity(&self, s: f64, i: f64) -> Result<f64> {
        let s = require_non_negative("s", s)?;
        let i = require_non_negative("i", i)?;
        let alpha = 1.0 + i / self.ki;
        let alpha_prime = 1.0 + i / self.ki_prime;
        Ok(self.base.vmax() * s / (alpha * self.base.km() + alpha_prime * s))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference tolerance for the float comparisons below.
    const EPS: f64 = 1e-9;

    fn base() -> MichaelisMenten {
        MichaelisMenten::new(4.0, 2.0).expect("valid base parameters")
    }

    // --- Competitive ---------------------------------------------------

    /// A competitive inhibitor raises the apparent Km but leaves Vmax
    /// untouched — the canonical diagnostic signature.
    #[test]
    fn competitive_raises_km_not_vmax() {
        let c = Competitive::new(base(), 1.0).expect("valid");
        let i = 3.0;
        let km0 = base().km();
        let km_app = c.apparent_km(i).expect("valid");
        // Km_app = Km*(1 + I/Ki) = 2*(1+3) = 8.
        assert!((km_app - km0 * (1.0 + i / 1.0)).abs() < EPS, "{km_app}");
        assert!(km_app > km0, "apparent Km should rise: {km_app} <= {km0}");
        // Vmax is unchanged at every inhibitor level.
        let vmax_app = c.apparent_vmax(i).expect("valid");
        assert!(
            (vmax_app - base().vmax()).abs() < EPS,
            "Vmax moved: {vmax_app}"
        );
    }

    /// Adding competitive inhibitor lowers the velocity at fixed finite S
    /// (because Km_app rises), yet the high-S asymptote is still Vmax.
    #[test]
    fn competitive_lowers_velocity_but_keeps_vmax_asymptote() {
        let c = Competitive::new(base(), 1.0).expect("valid");
        let s = 2.0;
        let v0 = c.velocity(s, 0.0).expect("valid");
        let vi = c.velocity(s, 5.0).expect("valid");
        assert!(vi < v0, "inhibitor should lower v: {vi} >= {v0}");
        // Surmountable: at very large S the velocity recovers toward Vmax.
        let v_big = c.velocity(1e9, 5.0).expect("valid");
        assert!(
            (base().vmax() - v_big).abs() < 1e-4,
            "should approach Vmax: {v_big}"
        );
    }

    /// With no inhibitor (I=0) every mode reduces to bare Michaelis-Menten.
    #[test]
    fn competitive_reduces_to_mm_without_inhibitor() {
        let c = Competitive::new(base(), 1.0).expect("valid");
        for &s in &[0.0, 0.5, 2.0, 10.0, 100.0] {
            let v = c.velocity(s, 0.0).expect("valid");
            let mm = base().velocity(s).expect("valid");
            assert!((v - mm).abs() < EPS, "mismatch at S={s}: {v} vs {mm}");
        }
    }

    // --- Noncompetitive ------------------------------------------------

    /// A noncompetitive inhibitor lowers the apparent Vmax but leaves Km
    /// untouched — its canonical diagnostic signature.
    #[test]
    fn noncompetitive_lowers_vmax_not_km() {
        let n = Noncompetitive::new(base(), 2.0).expect("valid");
        let i = 2.0;
        // Vmax_app = Vmax/(1 + I/Ki) = 4/(1+1) = 2.
        let vmax_app = n.apparent_vmax(i).expect("valid");
        assert!(
            (vmax_app - base().vmax() / (1.0 + i / 2.0)).abs() < EPS,
            "{vmax_app}"
        );
        assert!(vmax_app < base().vmax(), "Vmax should fall: {vmax_app}");
        // Km is unchanged at every inhibitor level.
        let km_app = n.apparent_km(i).expect("valid");
        assert!((km_app - base().km()).abs() < EPS, "Km moved: {km_app}");
    }

    /// Even at saturating substrate the noncompetitive-inhibited velocity
    /// cannot exceed the reduced Vmax_app (insurmountable inhibition).
    #[test]
    fn noncompetitive_caps_velocity_below_reduced_vmax() {
        let n = Noncompetitive::new(base(), 2.0).expect("valid");
        let i = 2.0;
        let vmax_app = n.apparent_vmax(i).expect("valid");
        let v_big = n.velocity(1e9, i).expect("valid");
        assert!(v_big <= vmax_app + EPS, "exceeded Vmax_app: {v_big}");
        assert!(
            (vmax_app - v_big).abs() < 1e-4,
            "should approach Vmax_app: {v_big}"
        );
    }

    /// Half-saturation still occurs at S = Km regardless of inhibitor,
    /// because noncompetitive inhibition does not move Km.
    #[test]
    fn noncompetitive_half_saturation_stays_at_km() {
        let n = Noncompetitive::new(base(), 0.7).expect("valid");
        let i = 4.0;
        let v_at_km = n.velocity(base().km(), i).expect("valid");
        let vmax_app = n.apparent_vmax(i).expect("valid");
        assert!(
            (v_at_km - vmax_app / 2.0).abs() < EPS,
            "v(Km) should be Vmax_app/2: {v_at_km}"
        );
    }

    // --- Uncompetitive -------------------------------------------------

    /// An uncompetitive inhibitor scales Vmax and Km by the *same* factor,
    /// leaving their ratio (the low-substrate slope) unchanged.
    #[test]
    fn uncompetitive_scales_both_equally() {
        let u = Uncompetitive::new(base(), 3.0).expect("valid");
        let i = 6.0; // factor 1 + 6/3 = 3
        let vmax_app = u.apparent_vmax(i).expect("valid");
        let km_app = u.apparent_km(i).expect("valid");
        assert!((vmax_app - base().vmax() / 3.0).abs() < EPS, "{vmax_app}");
        assert!((km_app - base().km() / 3.0).abs() < EPS, "{km_app}");
        // Ratio Vmax/Km is invariant.
        let ratio0 = base().vmax() / base().km();
        let ratio_i = vmax_app / km_app;
        assert!((ratio_i - ratio0).abs() < EPS, "ratio moved: {ratio_i}");
    }

    // --- Mixed (general case + reductions) -----------------------------

    /// Mixed inhibition with Ki' -> ∞ reproduces competitive inhibition.
    #[test]
    fn mixed_reduces_to_competitive_when_ki_prime_huge() {
        let ki = 1.5;
        let c = Competitive::new(base(), ki).expect("valid");
        let m = Mixed::new(base(), ki, 1e15).expect("valid");
        for &s in &[0.3, 2.0, 25.0] {
            for &i in &[0.0, 1.0, 8.0] {
                let vc = c.velocity(s, i).expect("valid");
                let vm = m.velocity(s, i).expect("valid");
                assert!((vc - vm).abs() < 1e-6, "S={s} I={i}: {vc} vs {vm}");
            }
        }
    }

    /// Mixed inhibition with Ki' == Ki reproduces pure noncompetitive
    /// inhibition (the symmetric case).
    #[test]
    fn mixed_reduces_to_noncompetitive_when_constants_equal() {
        let ki = 2.5;
        let n = Noncompetitive::new(base(), ki).expect("valid");
        let m = Mixed::new(base(), ki, ki).expect("valid");
        for &s in &[0.3, 2.0, 25.0] {
            for &i in &[0.0, 1.0, 8.0] {
                let vn = n.velocity(s, i).expect("valid");
                let vm = m.velocity(s, i).expect("valid");
                assert!((vn - vm).abs() < EPS, "S={s} I={i}: {vn} vs {vm}");
            }
        }
    }

    /// Mixed inhibition with Ki -> ∞ reproduces uncompetitive inhibition.
    #[test]
    fn mixed_reduces_to_uncompetitive_when_ki_huge() {
        let ki_prime = 3.0;
        let u = Uncompetitive::new(base(), ki_prime).expect("valid");
        let m = Mixed::new(base(), 1e15, ki_prime).expect("valid");
        for &s in &[0.3, 2.0, 25.0] {
            for &i in &[0.0, 1.0, 8.0] {
                let vu = u.velocity(s, i).expect("valid");
                let vm = m.velocity(s, i).expect("valid");
                assert!((vu - vm).abs() < 1e-6, "S={s} I={i}: {vu} vs {vm}");
            }
        }
    }

    // --- Validation + serde --------------------------------------------

    #[test]
    fn constructors_reject_bad_ki() {
        assert!(Competitive::new(base(), 0.0).is_err());
        assert!(Noncompetitive::new(base(), -1.0).is_err());
        assert!(Uncompetitive::new(base(), f64::NAN).is_err());
        assert!(Mixed::new(base(), 1.0, 0.0).is_err());
        assert!(Mixed::new(base(), -2.0, 1.0).is_err());
    }

    #[test]
    fn velocity_rejects_bad_arguments() {
        let c = Competitive::new(base(), 1.0).expect("valid");
        assert!(c.velocity(-1.0, 0.0).is_err());
        assert!(c.velocity(1.0, -1.0).is_err());
        let m = Mixed::new(base(), 1.0, 1.0).expect("valid");
        assert!(m.velocity(f64::NAN, 0.0).is_err());
        assert!(m.velocity(1.0, f64::INFINITY).is_err());
    }

    #[test]
    fn competitive_round_trips_through_json() {
        let c = Competitive::new(base(), 1.5).expect("valid");
        let json = serde_json::to_string(&c).expect("serialize");
        let back: Competitive = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(c, back);
    }
}
