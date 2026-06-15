//! The pH / pOH scale and the water autoionization relation.
//!
//! pH is the negative base-10 logarithm of the hydronium-ion
//! concentration, `pH = -log10([H+])`, and pOH the same for hydroxide.
//! In any aqueous solution the two are tied together by the water
//! self-ionization equilibrium `Kw = [H+][OH-]`, whose logarithmic form
//! is `pH + pOH = pKw`. At 25 C, `Kw = 1.0e-14` and `pKw = 14.00`, so
//! `pH + pOH = 14`.
//!
//! All functions here are pure conversions on already-known
//! concentrations or p-values; the model modules ([`crate::strong`],
//! [`crate::weak`], [`crate::buffer`]) compute the concentrations that
//! feed into them.

use crate::error::{validate_concentration, validate_kw, Result};

/// Water autoionization constant at 25 C: `Kw = 1.0e-14`.
pub const KW_25C: f64 = 1.0e-14;

/// `pKw = -log10(Kw) = 14.00` at 25 C. This is the value such that
/// `pH + pOH = PKW_25C` for any aqueous solution at 25 C.
pub const PKW_25C: f64 = 14.0;

/// pH from a hydronium-ion concentration: `pH = -log10([H+])`.
///
/// `h` is `[H+]` in `mol / L`.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConcentration`] when `h` is not
/// finite and strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_acidbase::ph;
/// // Neutral water at 25 C: [H+] = 1.0e-7 -> pH 7.
/// let p = ph(1.0e-7).unwrap();
/// assert!((p - 7.0).abs() < 1e-12);
/// ```
pub fn ph(h: f64) -> Result<f64> {
    let h = validate_concentration("[H+]", h)?;
    Ok(-h.log10())
}

/// pOH from a hydroxide-ion concentration: `pOH = -log10([OH-])`.
///
/// `oh` is `[OH-]` in `mol / L`.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConcentration`] when `oh` is not
/// finite and strictly positive.
pub fn poh(oh: f64) -> Result<f64> {
    let oh = validate_concentration("[OH-]", oh)?;
    Ok(-oh.log10())
}

/// Hydronium-ion concentration from pH: `[H+] = 10^(-pH)`.
///
/// The inverse of [`ph`]. Any finite `p` is accepted (a strongly basic
/// solution has a large positive pH and a vanishingly small `[H+]`).
///
/// # Examples
///
/// ```
/// use valenx_acidbase::h_from_ph;
/// let h = h_from_ph(3.0);
/// assert!((h - 1.0e-3).abs() < 1e-15);
/// ```
pub fn h_from_ph(p: f64) -> f64 {
    10f64.powf(-p)
}

/// Hydroxide-ion concentration from pOH: `[OH-] = 10^(-pOH)`.
///
/// The inverse of [`poh`].
pub fn oh_from_poh(poh_value: f64) -> f64 {
    10f64.powf(-poh_value)
}

/// pOH from pH using a caller-supplied `Kw`: `pOH = pKw - pH`.
///
/// Pass [`KW_25C`] for 25 C. Supplying `Kw` lets the relation track
/// temperature, since `Kw` (and therefore `pKw`) is temperature
/// dependent.
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConstant`] when `kw` is not finite
/// and strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_acidbase::{poh_from_ph, KW_25C};
/// // At 25 C, pH 7 implies pOH 7 (they sum to 14).
/// let q = poh_from_ph(7.0, KW_25C).unwrap();
/// assert!((q - 7.0).abs() < 1e-12);
/// ```
pub fn poh_from_ph(ph_value: f64, kw: f64) -> Result<f64> {
    let kw = validate_kw(kw)?;
    let pkw = -kw.log10();
    Ok(pkw - ph_value)
}

/// pH from pOH using a caller-supplied `Kw`: `pH = pKw - pOH`.
///
/// The symmetric counterpart of [`poh_from_ph`].
///
/// # Errors
///
/// Returns [`crate::AcidBaseError::BadConstant`] when `kw` is not finite
/// and strictly positive.
pub fn ph_from_poh(poh_value: f64, kw: f64) -> Result<f64> {
    let kw = validate_kw(kw)?;
    let pkw = -kw.log10();
    Ok(pkw - poh_value)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for analytic float comparisons.
    const EPS: f64 = 1e-12;

    #[test]
    fn neutral_water_is_ph_seven() {
        // Ground truth: pure water at 25 C has [H+] = 1.0e-7 M -> pH 7.0.
        let p = ph(1.0e-7).expect("neutral [H+] is valid");
        assert!((p - 7.0).abs() < EPS, "pH(1e-7) = {p}, expected 7");
    }

    #[test]
    fn ph_of_decimolar_proton_is_one() {
        // [H+] = 0.1 M -> pH = -log10(0.1) = 1.
        let p = ph(0.1).expect("0.1 M is valid");
        assert!((p - 1.0).abs() < EPS, "pH(0.1) = {p}, expected 1");
    }

    #[test]
    fn ph_and_h_round_trip() {
        // ph and h_from_ph are exact inverses over a wide range.
        for &target in &[0.5, 1.0, 3.7, 7.0, 10.2, 13.5] {
            let h = h_from_ph(target);
            let back = ph(h).expect("reconstructed [H+] is valid");
            assert!(
                (back - target).abs() < EPS,
                "round trip pH {target} -> [H+] {h} -> pH {back}"
            );
        }
    }

    #[test]
    fn poh_round_trip() {
        for &target in &[0.5, 7.0, 12.3] {
            let oh = oh_from_poh(target);
            let back = poh(oh).expect("reconstructed [OH-] is valid");
            assert!(
                (back - target).abs() < EPS,
                "pOH round trip {target} -> {back}"
            );
        }
    }

    #[test]
    fn ph_plus_poh_equals_fourteen_at_25c() {
        // Ground truth: pH + pOH = pKw = 14 at 25 C, for ANY solution.
        for &ph_value in &[0.0, 1.0, 3.5, 7.0, 9.2, 13.0, 14.0] {
            let q = poh_from_ph(ph_value, KW_25C).expect("Kw is valid");
            let sum = ph_value + q;
            assert!(
                (sum - PKW_25C).abs() < EPS,
                "pH {ph_value} + pOH {q} = {sum}, expected 14"
            );
        }
    }

    #[test]
    fn ph_from_poh_is_symmetric_inverse() {
        // pH -> pOH -> pH must return the original at fixed Kw.
        let ph_value = 4.2_f64;
        let q = poh_from_ph(ph_value, KW_25C).expect("Kw is valid");
        let back = ph_from_poh(q, KW_25C).expect("Kw is valid");
        assert!(
            (back - ph_value).abs() < EPS,
            "symmetry: {ph_value} -> {back}"
        );
    }

    #[test]
    fn pkw_constant_matches_kw_constant() {
        // The two published constants are mutually consistent.
        assert!((-KW_25C.log10() - PKW_25C).abs() < EPS);
    }

    #[test]
    fn nonpositive_concentration_rejected() {
        assert!(ph(0.0).is_err());
        assert!(ph(-1.0).is_err());
        assert!(poh(f64::NAN).is_err());
    }

    #[test]
    fn bad_kw_rejected_in_conversions() {
        assert!(poh_from_ph(7.0, 0.0).is_err());
        assert!(ph_from_poh(7.0, -1.0).is_err());
    }
}
