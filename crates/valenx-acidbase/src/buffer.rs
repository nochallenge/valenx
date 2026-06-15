//! Buffers: the Henderson-Hasselbalch equation and buffer capacity.
//!
//! A buffer is a weak acid `HA` together with its conjugate base `A-`.
//! Taking `-log10` of the equilibrium `Ka = [H+][A-] / [HA]` and
//! rearranging gives the **Henderson-Hasselbalch** equation
//!
//! ```text
//! pH = pKa + log10([A-] / [HA]),
//! ```
//!
//! so the pH equals `pKa` when the acid and conjugate base are present
//! in equal amounts (`[A-] = [HA]`, the log term vanishes), and rises as
//! more conjugate base is added.
//!
//! The resistance of the buffer to pH change on adding strong acid or
//! base is the **buffer capacity** `beta = d(Cb)/d(pH)`. For a single
//! weak-acid / conjugate-base pair the Van Slyke expression is
//!
//! ```text
//! beta = 2.303 * (Ka * [H+] * C_total) / (Ka + [H+])^2,
//! ```
//!
//! where `C_total = [HA] + [A-]` and `[H+] = 10^(-pH)`. `beta` is
//! maximal at `pH = pKa`, where it equals `2.303 * C_total / 4`.

use crate::error::{validate_ka, AcidBaseError, Result};
use serde::{Deserialize, Serialize};

/// `ln(10) ~= 2.302585`, the factor converting natural to base-10 logs
/// that appears in the buffer-capacity derivative.
const LN10: f64 = std::f64::consts::LN_10;

/// A weak-acid buffer: its `pka`, the weak-acid concentration `acid`
/// (`[HA]`) and conjugate-base concentration `base` (`[A-]`), in
/// `mol / L`.
///
/// Construct with [`Buffer::new`], which validates the two
/// concentrations are strictly positive (so the
/// Henderson-Hasselbalch log is well defined), then query
/// [`Buffer::ph`] and [`Buffer::capacity`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Buffer {
    /// `pKa = -log10(Ka)` of the weak acid.
    pub pka: f64,
    /// Weak-acid (`HA`) concentration in `mol / L` (positive).
    pub acid: f64,
    /// Conjugate-base (`A-`) concentration in `mol / L` (positive).
    pub base: f64,
}

impl Buffer {
    /// Build a validated buffer from `pka` and the acid / conjugate-base
    /// concentrations.
    ///
    /// # Errors
    ///
    /// Returns [`AcidBaseError::DegenerateBuffer`] when either `acid` or
    /// `base` is not finite and strictly positive (the
    /// `[A-] / [HA]` ratio would be undefined).
    pub fn new(pka: f64, acid: f64, base: f64) -> Result<Self> {
        if acid.is_finite() && acid > 0.0 && base.is_finite() && base > 0.0 {
            Ok(Self { pka, acid, base })
        } else {
            Err(AcidBaseError::DegenerateBuffer { acid, base })
        }
    }

    /// Build a buffer from a dissociation constant `ka` instead of
    /// `pKa`, converting `pKa = -log10(Ka)`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::AcidBaseError::BadConstant`] when `ka` is
    /// invalid, or [`AcidBaseError::DegenerateBuffer`] when a
    /// concentration is non-positive.
    pub fn from_ka(ka: f64, acid: f64, base: f64) -> Result<Self> {
        let ka = validate_ka("Ka", ka)?;
        Self::new(-ka.log10(), acid, base)
    }

    /// pH via Henderson-Hasselbalch: `pH = pKa + log10([A-] / [HA])`.
    pub fn ph(&self) -> f64 {
        self.pka + (self.base / self.acid).log10()
    }

    /// Van Slyke buffer capacity `beta` at this buffer's own pH, in
    /// `mol / L per pH unit`.
    ///
    /// Uses `beta = 2.303 * Ka * [H+] * C_total / (Ka + [H+])^2` with
    /// `Ka = 10^(-pKa)`, `[H+] = 10^(-pH)`, and
    /// `C_total = [HA] + [A-]`.
    pub fn capacity(&self) -> f64 {
        let ka = 10f64.powf(-self.pka);
        let h = 10f64.powf(-self.ph());
        let c_total = self.acid + self.base;
        let denom = ka + h;
        LN10 * ka * h * c_total / (denom * denom)
    }
}

/// Henderson-Hasselbalch pH from a `pka` and the conjugate-base /
/// weak-acid concentrations: `pH = pKa + log10([A-] / [HA])`.
///
/// Free-function form of [`Buffer::ph`]; validates the concentrations.
///
/// # Errors
///
/// Returns [`AcidBaseError::DegenerateBuffer`] when `acid` or `base` is
/// not finite and strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_acidbase::henderson_hasselbalch;
/// // Equal acid and base -> pH = pKa.
/// let p = henderson_hasselbalch(4.76, 0.1, 0.1).unwrap();
/// assert!((p - 4.76).abs() < 1e-12);
/// ```
pub fn henderson_hasselbalch(pka: f64, acid: f64, base: f64) -> Result<f64> {
    Ok(Buffer::new(pka, acid, base)?.ph())
}

/// Van Slyke buffer capacity `beta` of a weak-acid buffer, in
/// `mol / L per pH unit`.
///
/// Free-function form of [`Buffer::capacity`]; validates the
/// concentrations.
///
/// # Errors
///
/// Returns [`AcidBaseError::DegenerateBuffer`] when `acid` or `base` is
/// not finite and strictly positive.
pub fn buffer_capacity(pka: f64, acid: f64, base: f64) -> Result<f64> {
    Ok(Buffer::new(pka, acid, base)?.capacity())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for analytic float comparisons.
    const EPS: f64 = 1e-12;

    #[test]
    fn ph_equals_pka_at_equal_concentrations() {
        // Ground truth: when [A-] = [HA], the log term is 0 -> pH = pKa.
        let p = henderson_hasselbalch(4.76, 0.1, 0.1).expect("valid");
        assert!((p - 4.76).abs() < EPS, "pH = {p}, expected 4.76");
    }

    #[test]
    fn ph_equals_pka_independent_of_the_shared_concentration() {
        // The equality holds for any equal pair, not just 0.1 M.
        for &c in &[0.001, 0.05, 0.2, 1.0] {
            let p = henderson_hasselbalch(7.0, c, c).expect("valid");
            assert!((p - 7.0).abs() < EPS, "pH = {p} at C = {c}, expected 7");
        }
    }

    #[test]
    fn ph_rises_with_more_conjugate_base() {
        // Ground truth: adding conjugate base raises the pH above pKa.
        let pka = 4.76;
        let equal = henderson_hasselbalch(pka, 0.1, 0.1).expect("valid");
        let more_base = henderson_hasselbalch(pka, 0.1, 0.5).expect("valid");
        assert!(
            more_base > equal,
            "more-base pH {more_base} should exceed equal pH {equal}"
        );
        // 10:1 base:acid -> pH = pKa + log10(10) = pKa + 1 exactly.
        let ten_to_one = henderson_hasselbalch(pka, 0.1, 1.0).expect("valid");
        assert!(
            (ten_to_one - (pka + 1.0)).abs() < EPS,
            "10:1 ratio pH = {ten_to_one}, expected {}",
            pka + 1.0
        );
    }

    #[test]
    fn ph_falls_with_more_acid() {
        // Symmetric check: a 1:10 base:acid ratio gives pKa - 1.
        let pka = 4.76;
        let one_to_ten = henderson_hasselbalch(pka, 1.0, 0.1).expect("valid");
        assert!(
            (one_to_ten - (pka - 1.0)).abs() < EPS,
            "1:10 ratio pH = {one_to_ten}, expected {}",
            pka - 1.0
        );
    }

    #[test]
    fn from_ka_matches_pka_construction() {
        // Building from Ka = 1.8e-5 must equal building from pKa 4.745.
        let via_ka = Buffer::from_ka(1.8e-5, 0.1, 0.1).expect("valid");
        let expected_pka = -(1.8e-5_f64).log10();
        let via_pka = Buffer::new(expected_pka, 0.1, 0.1).expect("valid");
        assert!((via_ka.ph() - via_pka.ph()).abs() < EPS);
        // And at equal concentrations the pH equals that pKa.
        assert!((via_ka.ph() - expected_pka).abs() < EPS);
    }

    #[test]
    fn capacity_is_maximal_at_pka() {
        // Ground truth: beta peaks where pH = pKa (equal acid/base) and
        // there equals 2.303 * C_total / 4. Off-peak ratios give less.
        let pka = 4.76;
        let c_total = 0.2; // 0.1 + 0.1
        let peak = buffer_capacity(pka, 0.1, 0.1).expect("valid");
        let expected_peak = LN10 * c_total / 4.0;
        assert!(
            (peak - expected_peak).abs() < 1e-9,
            "peak beta = {peak}, expected {expected_peak}"
        );

        // Same total concentration but skewed ratio -> lower capacity.
        let skewed = buffer_capacity(pka, 0.18, 0.02).expect("valid");
        assert!(
            skewed < peak,
            "skewed beta {skewed} should be below peak {peak}"
        );
    }

    #[test]
    fn capacity_scales_with_total_concentration() {
        // Doubling both concentrations (same ratio) doubles beta.
        let beta1 = buffer_capacity(4.76, 0.1, 0.1).expect("valid");
        let beta2 = buffer_capacity(4.76, 0.2, 0.2).expect("valid");
        assert!((beta2 - 2.0 * beta1).abs() < 1e-9, "{beta2} vs 2*{beta1}");
    }

    #[test]
    fn capacity_is_positive() {
        let beta = buffer_capacity(7.0, 0.05, 0.15).expect("valid");
        assert!(beta > 0.0, "beta = {beta} must be positive");
    }

    #[test]
    fn degenerate_buffer_rejected() {
        assert!(henderson_hasselbalch(4.76, 0.0, 0.1).is_err());
        assert!(henderson_hasselbalch(4.76, 0.1, 0.0).is_err());
        assert!(henderson_hasselbalch(4.76, -0.1, 0.1).is_err());
        assert!(buffer_capacity(4.76, 0.1, f64::NAN).is_err());
    }

    #[test]
    fn from_ka_rejects_bad_constant() {
        assert!(Buffer::from_ka(0.0, 0.1, 0.1).is_err());
        assert!(Buffer::from_ka(-1.0, 0.1, 0.1).is_err());
    }
}
