//! Strong (fully dissociating) monoprotic acids and bases.
//!
//! A strong monoprotic acid `HA` dissociates completely,
//! `HA -> H+ + A-`, so its hydronium concentration equals its formal
//! concentration, `[H+] = C`, giving `pH = -log10(C)`. A strong
//! monoprotic base `BOH -> B+ + OH-` does the same for hydroxide,
//! `[OH-] = C`, with `pH = pKw + log10(C)`.
//!
//! These formulas assume the acid / base concentration dominates the
//! `1.0e-7 M` of protons / hydroxide from water self-ionization — the
//! usual textbook regime for `C` well above `1.0e-6 M`. The validated
//! constructors reject non-positive `C` but do **not** attempt the full
//! charge-balance correction needed in the extremely dilute limit (see
//! the crate-level honest-scope note).

use crate::error::{validate_concentration, Result};
use crate::ph::{poh_from_ph, KW_25C};

/// pH of a strong monoprotic acid at formal concentration `c` (mol/L).
///
/// Complete dissociation gives `[H+] = c`, hence `pH = -log10(c)`.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConcentration`] when `c` is not
/// finite and strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_acidbase::ph_strong_acid;
/// // 0.1 M HCl -> pH 1.
/// let p = ph_strong_acid(0.1).unwrap();
/// assert!((p - 1.0).abs() < 1e-12);
/// ```
pub fn ph_strong_acid(c: f64) -> Result<f64> {
    let c = validate_concentration("C", c)?;
    Ok(-c.log10())
}

/// pH of a strong monoprotic base at formal concentration `c` (mol/L),
/// using a caller-supplied water constant `kw`.
///
/// Complete dissociation gives `[OH-] = c`, so `pOH = -log10(c)` and
/// `pH = pKw - pOH`. Pass [`crate::KW_25C`] for 25 C.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConcentration`] when `c` is not
/// finite and strictly positive, or
/// [`crate::AcidBaseError::BadConstant`] when `kw` is invalid.
///
/// # Examples
///
/// ```
/// use valenx_acidbase::{ph_strong_base, KW_25C};
/// // 0.1 M NaOH at 25 C -> pOH 1 -> pH 13.
/// let p = ph_strong_base(0.1, KW_25C).unwrap();
/// assert!((p - 13.0).abs() < 1e-12);
/// ```
pub fn ph_strong_base(c: f64, kw: f64) -> Result<f64> {
    let c = validate_concentration("C", c)?;
    let poh = -c.log10();
    // pH = pKw - pOH; reuse the validated water relation (also checks kw).
    ph_from_poh_internal(poh, kw)
}

/// Convenience wrapper around [`crate::ph::ph_from_poh`] that keeps the
/// strong-base arithmetic in one place. Re-derives `pH = pKw - pOH`.
fn ph_from_poh_internal(poh: f64, kw: f64) -> Result<f64> {
    // poh_from_ph(0, kw) == pKw; subtract pOH to get pH. Doing it this way
    // funnels Kw validation through a single audited path.
    let pkw = poh_from_ph(0.0, kw)?;
    Ok(pkw - poh)
}

/// Convenience: pH of a strong base at 25 C, where `kw = `[`crate::KW_25C`].
///
/// Equivalent to `ph_strong_base(c, KW_25C)`.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConcentration`] when `c` is not
/// finite and strictly positive.
pub fn ph_strong_base_25c(c: f64) -> Result<f64> {
    ph_strong_base(c, KW_25C)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for analytic float comparisons.
    const EPS: f64 = 1e-12;

    #[test]
    fn decimolar_strong_acid_is_ph_one() {
        // Ground truth: 0.1 M strong acid (HCl) -> pH 1.
        let p = ph_strong_acid(0.1).expect("0.1 M is valid");
        assert!((p - 1.0).abs() < EPS, "pH = {p}, expected 1");
    }

    #[test]
    fn molar_strong_acid_is_ph_zero() {
        // 1.0 M strong acid -> pH = -log10(1) = 0.
        let p = ph_strong_acid(1.0).expect("1 M is valid");
        assert!(p.abs() < EPS, "pH = {p}, expected 0");
    }

    #[test]
    fn millimolar_strong_acid_is_ph_three() {
        let p = ph_strong_acid(1.0e-3).expect("1 mM is valid");
        assert!((p - 3.0).abs() < EPS, "pH = {p}, expected 3");
    }

    #[test]
    fn decimolar_strong_base_is_ph_thirteen() {
        // Ground truth: 0.1 M NaOH at 25 C -> pOH 1 -> pH 13.
        let p = ph_strong_base(0.1, KW_25C).expect("0.1 M / valid Kw");
        assert!((p - 13.0).abs() < EPS, "pH = {p}, expected 13");
    }

    #[test]
    fn molar_strong_base_is_ph_fourteen() {
        let p = ph_strong_base(1.0, KW_25C).expect("1 M / valid Kw");
        assert!((p - 14.0).abs() < EPS, "pH = {p}, expected 14");
    }

    #[test]
    fn strong_base_convenience_matches_explicit_kw() {
        let a = ph_strong_base_25c(0.05).expect("valid");
        let b = ph_strong_base(0.05, KW_25C).expect("valid");
        assert!((a - b).abs() < EPS, "convenience {a} != explicit {b}");
    }

    #[test]
    fn strong_acid_and_base_sum_to_pkw_at_equal_concentration() {
        // For equal C, pH(acid) + pH(base) = pKw, because the acid sits at
        // -log10(C) and the base at pKw + log10(C).
        let c = 0.02_f64;
        let acid = ph_strong_acid(c).expect("valid");
        let base = ph_strong_base(c, KW_25C).expect("valid");
        assert!(
            (acid + base - 14.0).abs() < EPS,
            "acid {acid} + base {base} = {} (expected 14)",
            acid + base
        );
    }

    #[test]
    fn nonpositive_concentration_rejected() {
        assert!(ph_strong_acid(0.0).is_err());
        assert!(ph_strong_acid(-0.1).is_err());
        assert!(ph_strong_base(0.0, KW_25C).is_err());
    }

    #[test]
    fn bad_kw_rejected() {
        assert!(ph_strong_base(0.1, 0.0).is_err());
        assert!(ph_strong_base(0.1, f64::NAN).is_err());
    }
}
