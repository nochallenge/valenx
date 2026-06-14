//! Weak monoprotic acid equilibrium.
//!
//! A weak acid `HA` only partially dissociates,
//! `HA <-> H+ + A-`, governed by `Ka = [H+][A-] / [HA]`. Writing the
//! amount dissociated as `x = [H+] = [A-]` and the remaining acid as
//! `C - x`, the equilibrium is the quadratic
//!
//! ```text
//! Ka = x^2 / (C - x)   ==>   x^2 + Ka*x - Ka*C = 0.
//! ```
//!
//! When `Ka << C` the dissociation is slight (`x << C`), so dropping the
//! `-x` term gives the standard approximation `x = sqrt(Ka*C)`, hence
//! `pH = -log10(sqrt(Ka*C))`. [`ph_weak_acid_exact`] keeps the `-x` and
//! takes the positive root of the quadratic; the two agree closely for
//! weak acids and diverge as `Ka` approaches `C`.
//!
//! Self-ionization of water is neglected, valid while
//! `[H+] >> 1.0e-7 M` (the usual regime away from extreme dilution).

use crate::error::{validate_concentration, validate_ka, Result};
use serde::{Deserialize, Serialize};

/// A weak monoprotic acid: its acid-dissociation constant `ka` and
/// formal (analytical) concentration `c` in `mol / L`.
///
/// Construct with [`WeakAcid::new`], which validates both fields, then
/// query [`WeakAcid::ph`], [`WeakAcid::ph_exact`], or
/// [`WeakAcid::fraction_dissociated`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct WeakAcid {
    /// Acid-dissociation constant `Ka = [H+][A-] / [HA]` (positive).
    pub ka: f64,
    /// Formal concentration `C` of the acid in `mol / L` (positive).
    pub c: f64,
}

impl WeakAcid {
    /// Build a validated weak acid from `ka` and concentration `c`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::AcidBaseError::BadConstant`] when `ka` is not
    /// finite and strictly positive, or
    /// [`crate::AcidBaseError::BadConcentration`] when `c` is not finite
    /// and strictly positive.
    pub fn new(ka: f64, c: f64) -> Result<Self> {
        let ka = validate_ka("Ka", ka)?;
        let c = validate_concentration("C", c)?;
        Ok(Self { ka, c })
    }

    /// pH from the small-`x` approximation `[H+] = sqrt(Ka*C)`.
    ///
    /// Accurate when `Ka << C`; see [`WeakAcid::ph_exact`] for the
    /// quadratic that drops this assumption.
    pub fn ph(&self) -> f64 {
        let h = (self.ka * self.c).sqrt();
        -h.log10()
    }

    /// pH from the exact equilibrium quadratic
    /// `x^2 + Ka*x - Ka*C = 0`, taking the positive root for `[H+]`.
    ///
    /// This makes no `x << C` assumption and is the reference the
    /// approximation is compared against.
    pub fn ph_exact(&self) -> f64 {
        let h = self.h_exact();
        -h.log10()
    }

    /// Equilibrium `[H+]` from the exact quadratic (positive root).
    ///
    /// `x = (-Ka + sqrt(Ka^2 + 4*Ka*C)) / 2`. The discriminant is always
    /// positive for positive `Ka`, `C`, so the root is real and the `+`
    /// branch is the physical (positive) one.
    pub fn h_exact(&self) -> f64 {
        let ka = self.ka;
        let disc = ka * ka + 4.0 * ka * self.c;
        (-ka + disc.sqrt()) / 2.0
    }

    /// Fraction of the acid that has dissociated at equilibrium,
    /// `alpha = [H+] / C`, from the exact root. Lies in `(0, 1)`.
    ///
    /// This is the degree of ionization; it rises toward 1 as the
    /// solution is diluted (Ostwald dilution law).
    pub fn fraction_dissociated(&self) -> f64 {
        self.h_exact() / self.c
    }
}

/// pH of a weak monoprotic acid via the approximation
/// `[H+] = sqrt(Ka*C)`.
///
/// Free-function form of [`WeakAcid::ph`]; validates its inputs.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConstant`] /
/// [`crate::AcidBaseError::BadConcentration`] for invalid `ka` / `c`.
///
/// # Examples
///
/// ```
/// use valenx_acidbase::ph_weak_acid;
/// // 0.1 M acetic acid, Ka = 1.8e-5 -> pH ~ 2.87.
/// let p = ph_weak_acid(1.8e-5, 0.1).unwrap();
/// assert!((p - 2.87).abs() < 0.01);
/// ```
pub fn ph_weak_acid(ka: f64, c: f64) -> Result<f64> {
    Ok(WeakAcid::new(ka, c)?.ph())
}

/// pH of a weak monoprotic acid via the exact equilibrium quadratic.
///
/// Free-function form of [`WeakAcid::ph_exact`]; validates its inputs.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConstant`] /
/// [`crate::AcidBaseError::BadConcentration`] for invalid `ka` / `c`.
pub fn ph_weak_acid_exact(ka: f64, c: f64) -> Result<f64> {
    Ok(WeakAcid::new(ka, c)?.ph_exact())
}

/// Degree of dissociation `alpha = [H+] / C` of a weak acid, from the
/// exact root.
///
/// Free-function form of [`WeakAcid::fraction_dissociated`]; validates
/// its inputs.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConstant`] /
/// [`crate::AcidBaseError::BadConcentration`] for invalid `ka` / `c`.
pub fn fraction_dissociated(ka: f64, c: f64) -> Result<f64> {
    Ok(WeakAcid::new(ka, c)?.fraction_dissociated())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ph::ph;
    use crate::strong::ph_strong_acid;

    /// Tolerance for analytic float comparisons.
    const EPS: f64 = 1e-12;

    #[test]
    fn acetic_acid_approximation_matches_textbook() {
        // Ground truth (textbook): 0.1 M acetic acid, Ka = 1.8e-5.
        // [H+] = sqrt(1.8e-5 * 0.1) = sqrt(1.8e-6) = 1.3416e-3 -> pH 2.872.
        let p = ph_weak_acid(1.8e-5, 0.1).expect("valid acid");
        assert!((p - 2.8724).abs() < 1e-3, "pH = {p}, expected ~2.872");
    }

    #[test]
    fn approximation_equals_closed_form_h() {
        // The reported pH must equal -log10(sqrt(Ka*C)) exactly.
        let ka = 1.8e-5;
        let c = 0.1;
        let p = ph_weak_acid(ka, c).expect("valid");
        let expected = -((ka * c).sqrt()).log10();
        assert!((p - expected).abs() < EPS, "{p} vs {expected}");
    }

    #[test]
    fn weak_acid_ph_between_strong_acid_and_neutral() {
        // Ground truth: a weak acid is less acidic than a strong acid of
        // the same concentration but still acidic (pH < 7).
        let c = 0.1;
        let weak = ph_weak_acid(1.8e-5, c).expect("valid");
        let strong = ph_strong_acid(c).expect("valid"); // pH 1
        let neutral = ph(1.0e-7).expect("valid"); // pH 7
        assert!(
            strong < weak && weak < neutral,
            "expected {strong} < {weak} < {neutral}"
        );
    }

    #[test]
    fn exact_and_approx_agree_for_weak_acid() {
        // For Ka << C the quadratic and the sqrt approximation nearly
        // coincide (acetic acid: agreement to < 0.01 pH units).
        let approx = ph_weak_acid(1.8e-5, 0.1).expect("valid");
        let exact = ph_weak_acid_exact(1.8e-5, 0.1).expect("valid");
        assert!(
            (approx - exact).abs() < 0.01,
            "approx {approx} vs exact {exact}"
        );
    }

    #[test]
    fn exact_solver_satisfies_the_equilibrium_quadratic() {
        // Substitute the exact root back: x^2 + Ka*x - Ka*C must be ~0.
        let acid = WeakAcid::new(1.0e-3, 0.05).expect("valid");
        let x = acid.h_exact();
        let residual = x * x + acid.ka * x - acid.ka * acid.c;
        assert!(residual.abs() < 1e-15, "quadratic residual = {residual}");
    }

    #[test]
    fn exact_is_more_acidic_than_approximation() {
        // Keeping the -x term leaves less acid in the denominator, so the
        // exact [H+] is slightly LOWER, i.e. exact pH is slightly HIGHER
        // than the approximation. Use a stronger weak acid to make the gap
        // visible.
        let approx = ph_weak_acid(1.0e-2, 0.1).expect("valid");
        let exact = ph_weak_acid_exact(1.0e-2, 0.1).expect("valid");
        assert!(
            exact > approx,
            "exact {exact} should exceed approx {approx}"
        );
    }

    #[test]
    fn fraction_dissociated_is_between_zero_and_one() {
        let alpha = fraction_dissociated(1.8e-5, 0.1).expect("valid");
        assert!(alpha > 0.0 && alpha < 1.0, "alpha = {alpha} out of range");
        // Acetic acid at 0.1 M is ~1.3% dissociated.
        assert!(
            (alpha - 0.0134).abs() < 1e-3,
            "alpha = {alpha}, expected ~0.0134"
        );
    }

    #[test]
    fn dilution_raises_fraction_dissociated() {
        // Ostwald dilution: diluting a weak acid increases its degree of
        // ionization.
        let concentrated = fraction_dissociated(1.8e-5, 1.0).expect("valid");
        let dilute = fraction_dissociated(1.8e-5, 0.001).expect("valid");
        assert!(
            dilute > concentrated,
            "dilute {dilute} should exceed concentrated {concentrated}"
        );
    }

    #[test]
    fn stronger_ka_lowers_ph() {
        // Larger Ka (stronger acid) at fixed C means lower pH.
        let weaker = ph_weak_acid(1.0e-6, 0.1).expect("valid");
        let stronger = ph_weak_acid(1.0e-4, 0.1).expect("valid");
        assert!(stronger < weaker, "stronger {stronger} < weaker {weaker}");
    }

    #[test]
    fn invalid_inputs_rejected() {
        assert!(WeakAcid::new(0.0, 0.1).is_err());
        assert!(WeakAcid::new(1.8e-5, 0.0).is_err());
        assert!(ph_weak_acid(-1.0, 0.1).is_err());
        assert!(ph_weak_acid_exact(1.8e-5, f64::INFINITY).is_err());
    }
}
