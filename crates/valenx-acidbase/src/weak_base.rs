//! Weak monoprotic base equilibrium.
//!
//! A weak base `B` only partially protonates water,
//! `B + H2O <-> BH+ + OH-`, governed by `Kb = [BH+][OH-] / [B]`. Writing
//! the amount protonated as `x = [OH-] = [BH+]` and the remaining base as
//! `C - x`, the equilibrium is the quadratic
//!
//! ```text
//! Kb = x^2 / (C - x)   ==>   x^2 + Kb*x - Kb*C = 0.
//! ```
//!
//! When `Kb << C` the protonation is slight (`x << C`), so dropping the
//! `-x` term gives the standard approximation `x = sqrt(Kb*C)`, hence
//! `pOH = -log10(sqrt(Kb*C))` and `pH = pKw - pOH`. [`WeakBase::poh_exact`]
//! keeps the `-x` and takes the positive root of the quadratic.
//!
//! This is the exact mirror of the weak-[acid](crate::weak) module with
//! `Ka -> Kb` and `[H+] -> [OH-]`; self-ionization of water is neglected,
//! valid while `[OH-] >> 1.0e-7 M`.

use crate::error::{validate_concentration, validate_ka, Result};
use serde::{Deserialize, Serialize};

/// A weak monoprotic base: its base-dissociation constant `kb` and formal
/// (analytical) concentration `c` in `mol / L`.
///
/// Construct with [`WeakBase::new`], which validates both fields, then
/// query [`WeakBase::poh`], [`WeakBase::poh_exact`], [`WeakBase::ph`], or
/// [`WeakBase::fraction_protonated`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct WeakBase {
    /// Base-dissociation constant `Kb = [BH+][OH-] / [B]` (positive).
    pub kb: f64,
    /// Formal concentration `C` of the base in `mol / L` (positive).
    pub c: f64,
}

impl WeakBase {
    /// Build a validated weak base from `kb` and concentration `c`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::AcidBaseError::BadConstant`] when `kb` is not
    /// finite and strictly positive, or
    /// [`crate::AcidBaseError::BadConcentration`] when `c` is not finite
    /// and strictly positive.
    pub fn new(kb: f64, c: f64) -> Result<Self> {
        let kb = validate_ka("Kb", kb)?;
        let c = validate_concentration("C", c)?;
        Ok(Self { kb, c })
    }

    /// pOH from the small-`x` approximation `[OH-] = sqrt(Kb*C)`.
    pub fn poh(&self) -> f64 {
        let oh = (self.kb * self.c).sqrt();
        -oh.log10()
    }

    /// pOH from the exact equilibrium quadratic `x^2 + Kb*x - Kb*C = 0`,
    /// taking the positive root for `[OH-]`.
    pub fn poh_exact(&self) -> f64 {
        -self.oh_exact().log10()
    }

    /// Equilibrium `[OH-]` from the exact quadratic (positive root)
    /// `x = (-Kb + sqrt(Kb^2 + 4*Kb*C)) / 2`.
    pub fn oh_exact(&self) -> f64 {
        let kb = self.kb;
        let disc = kb * kb + 4.0 * kb * self.c;
        (-kb + disc.sqrt()) / 2.0
    }

    /// pH from the approximate pOH and the water constant `kw`, via
    /// `pH = pKw - pOH`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::AcidBaseError::BadConstant`] if `kw` is not finite
    /// and strictly positive.
    pub fn ph(&self, kw: f64) -> Result<f64> {
        crate::ph::ph_from_poh(self.poh(), kw)
    }

    /// Fraction of the base that has protonated at equilibrium,
    /// `alpha = [OH-] / C`, from the exact root. Lies in `(0, 1)` and rises
    /// toward 1 on dilution (Ostwald dilution law).
    pub fn fraction_protonated(&self) -> f64 {
        self.oh_exact() / self.c
    }
}

/// pOH of a weak monoprotic base via the approximation `[OH-] = sqrt(Kb*C)`.
///
/// Free-function form of [`WeakBase::poh`]; validates its inputs.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConstant`] /
/// [`crate::AcidBaseError::BadConcentration`] for invalid `kb` / `c`.
///
/// # Examples
///
/// ```
/// use valenx_acidbase::poh_weak_base;
/// // 0.1 M ammonia, Kb = 1.8e-5 -> pOH ~ 2.87 (pH ~ 11.13).
/// let p = poh_weak_base(1.8e-5, 0.1).unwrap();
/// assert!((p - 2.87).abs() < 0.01);
/// ```
pub fn poh_weak_base(kb: f64, c: f64) -> Result<f64> {
    Ok(WeakBase::new(kb, c)?.poh())
}

/// pOH of a weak monoprotic base via the exact equilibrium quadratic.
///
/// Free-function form of [`WeakBase::poh_exact`]; validates its inputs.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConstant`] /
/// [`crate::AcidBaseError::BadConcentration`] for invalid `kb` / `c`.
pub fn poh_weak_base_exact(kb: f64, c: f64) -> Result<f64> {
    Ok(WeakBase::new(kb, c)?.poh_exact())
}

/// pH of a weak monoprotic base at water constant `kw`, via
/// `pH = pKw - pOH` with the approximate pOH.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConstant`] /
/// [`crate::AcidBaseError::BadConcentration`] for invalid `kb` / `c`, or a
/// bad `kw`.
pub fn ph_weak_base(kb: f64, c: f64, kw: f64) -> Result<f64> {
    WeakBase::new(kb, c)?.ph(kw)
}

/// Degree of protonation `alpha = [OH-] / C` of a weak base, from the exact
/// root.
///
/// Free-function form of [`WeakBase::fraction_protonated`]; validates its
/// inputs.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConstant`] /
/// [`crate::AcidBaseError::BadConcentration`] for invalid `kb` / `c`.
pub fn fraction_protonated(kb: f64, c: f64) -> Result<f64> {
    Ok(WeakBase::new(kb, c)?.fraction_protonated())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ph::{ph, PKW_25C};
    use crate::strong::ph_strong_base;
    use crate::weak::ph_weak_acid;

    /// Tolerance for analytic float comparisons.
    const EPS: f64 = 1e-12;

    #[test]
    fn ammonia_matches_textbook() {
        // Ground truth: 0.1 M ammonia, Kb = 1.8e-5.
        // [OH-] = sqrt(1.8e-6) = 1.3416e-3 -> pOH 2.872, pH 11.128.
        let poh = poh_weak_base(1.8e-5, 0.1).expect("valid base");
        assert!((poh - 2.8724).abs() < 1e-3, "pOH = {poh}");
        let ph_val = ph_weak_base(1.8e-5, 0.1, crate::ph::KW_25C).expect("valid");
        assert!((ph_val - 11.1276).abs() < 1e-3, "pH = {ph_val}");
    }

    #[test]
    fn poh_equals_closed_form_oh() {
        let kb = 1.8e-5;
        let c = 0.1;
        let p = poh_weak_base(kb, c).expect("valid");
        let expected = -((kb * c).sqrt()).log10();
        assert!((p - expected).abs() < EPS, "{p} vs {expected}");
    }

    /// Acid-base symmetry: the pOH of a weak base equals the pH of a weak
    /// acid with the same constant and concentration (identical math, [OH-]
    /// vs [H+]).
    #[test]
    fn weak_base_poh_mirrors_weak_acid_ph() {
        let poh = poh_weak_base(1.8e-5, 0.1).expect("valid");
        let acid_ph = ph_weak_acid(1.8e-5, 0.1).expect("valid");
        assert!(
            (poh - acid_ph).abs() < EPS,
            "pOH {poh} vs acid pH {acid_ph}"
        );
    }

    #[test]
    fn weak_base_ph_between_neutral_and_strong_base() {
        let c = 0.1;
        let weak = ph_weak_base(1.8e-5, c, crate::ph::KW_25C).expect("valid");
        let strong = ph_strong_base(c, crate::ph::KW_25C).expect("valid"); // pH 13
        let neutral = ph(1.0e-7).expect("valid"); // pH 7
        assert!(
            neutral < weak && weak < strong,
            "expected {neutral} < {weak} < {strong}"
        );
    }

    #[test]
    fn exact_solver_satisfies_the_equilibrium_quadratic() {
        let base = WeakBase::new(1.0e-3, 0.05).expect("valid");
        let x = base.oh_exact();
        let residual = x * x + base.kb * x - base.kb * base.c;
        assert!(residual.abs() < 1e-15, "quadratic residual = {residual}");
    }

    #[test]
    fn exact_and_approx_agree_for_weak_base() {
        let approx = poh_weak_base(1.8e-5, 0.1).expect("valid");
        let exact = poh_weak_base_exact(1.8e-5, 0.1).expect("valid");
        assert!(
            (approx - exact).abs() < 0.01,
            "approx {approx} vs exact {exact}"
        );
    }

    #[test]
    fn fraction_protonated_in_range_and_rises_on_dilution() {
        let alpha = fraction_protonated(1.8e-5, 0.1).expect("valid");
        assert!(alpha > 0.0 && alpha < 1.0, "alpha = {alpha}");
        assert!(
            (alpha - 0.0134).abs() < 1e-3,
            "alpha = {alpha}, expected ~0.0134"
        );
        let dilute = fraction_protonated(1.8e-5, 0.001).expect("valid");
        assert!(dilute > alpha, "dilute {dilute} should exceed {alpha}");
    }

    #[test]
    fn ph_is_pkw_minus_poh() {
        let kb = 1.8e-5;
        let c = 0.1;
        let ph_val = ph_weak_base(kb, c, crate::ph::KW_25C).expect("valid");
        let poh = poh_weak_base(kb, c).expect("valid");
        assert!(
            (ph_val - (PKW_25C - poh)).abs() < 1e-9,
            "{ph_val} vs {}",
            PKW_25C - poh
        );
    }

    #[test]
    fn stronger_kb_raises_ph() {
        // Larger Kb (stronger base) at fixed C means higher pH.
        let weaker = ph_weak_base(1.0e-6, 0.1, crate::ph::KW_25C).expect("valid");
        let stronger = ph_weak_base(1.0e-4, 0.1, crate::ph::KW_25C).expect("valid");
        assert!(stronger > weaker, "stronger {stronger} > weaker {weaker}");
    }

    #[test]
    fn invalid_inputs_rejected() {
        assert!(WeakBase::new(0.0, 0.1).is_err());
        assert!(WeakBase::new(1.8e-5, 0.0).is_err());
        assert!(poh_weak_base(-1.0, 0.1).is_err());
        assert!(poh_weak_base_exact(1.8e-5, f64::INFINITY).is_err());
    }
}
