//! Sound-pressure level and decibel arithmetic.
//!
//! ## Model
//!
//! Sound-pressure level (SPL) expresses an RMS sound pressure `p` on a
//! logarithmic decibel scale relative to the standard reference pressure
//! [`P_REF`] = 20 micropascals (the nominal threshold of human hearing in
//! air):
//!
//! ```text
//! L = 20 * log10(p / p_ref)        [dB]
//! ```
//!
//! Because power is proportional to pressure squared, levels combine on a
//! *power* basis. For `N` mutually **incoherent** sources (random relative
//! phase — the usual assumption for unrelated noise sources) the combined
//! level is the decibel sum of their mean-square pressures:
//!
//! ```text
//! L_total = 10 * log10( sum_i 10^(L_i / 10) )        [dB]
//! ```
//!
//! Two consequences fall straight out of these definitions and are the
//! everyday rules of thumb of room acoustics:
//!
//! - Doubling the pressure adds `20*log10(2) ≈ 6.02` dB.
//! - Two **equal** incoherent sources add `10*log10(2) ≈ 3.01` dB.
//!
//! ## Honest scope
//!
//! These are the textbook closed-form definitions. They assume a single
//! RMS pressure already in hand, free-field combination on a mean-square
//! basis, and the standard 20 micropascal air reference. There is no
//! frequency weighting (A / C), no octave-band machinery, and no
//! coherent-interference term. This is research/educational grade, not a
//! calibrated metrology or regulatory measurement chain.

use crate::error::{AcousticsError, Result};

/// Reference sound pressure for airborne SPL: 20 micropascals
/// (`2e-5` Pa), the standard 0 dB datum.
pub const P_REF: f64 = 20e-6;

/// `20*log10(2)` — the decibel increase from doubling sound *pressure*.
/// Provided as a named constant so callers can assert against the exact
/// value rather than re-deriving the `≈ 6.0206` dB rule of thumb.
pub const DOUBLE_PRESSURE_DB: f64 = 6.020_599_913_279_624;

/// `10*log10(2)` — the decibel increase from summing two *equal*
/// incoherent sources (a doubling of acoustic *power*).
pub const DOUBLE_POWER_DB: f64 = 3.010_299_956_639_812;

fn check_pressure(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(AcousticsError::InvalidPressure { name, value })
    }
}

/// Sound-pressure level of an RMS pressure `pressure` (pascals) relative
/// to the standard [`P_REF`] = 20 micropascal reference.
///
/// `L = 20 * log10(pressure / P_REF)` in decibels.
///
/// At the reference pressure this is exactly `0` dB.
///
/// # Errors
///
/// Returns [`AcousticsError::InvalidPressure`] if `pressure` is not finite
/// and strictly positive.
///
/// # Examples
///
/// ```
/// use valenx_acoustics::spl::{spl, P_REF};
/// // The reference pressure itself is the 0 dB datum.
/// assert!((spl(P_REF).unwrap()).abs() < 1e-9);
/// ```
pub fn spl(pressure: f64) -> Result<f64> {
    spl_ref(pressure, P_REF)
}

/// Sound-pressure level of `pressure` relative to an explicit reference
/// pressure `p_ref`, both in pascals.
///
/// `L = 20 * log10(pressure / p_ref)` in decibels. Use [`spl`] for the
/// standard 20 micropascal air reference.
///
/// # Errors
///
/// Returns [`AcousticsError::InvalidPressure`] if either `pressure` or
/// `p_ref` is not finite and strictly positive.
pub fn spl_ref(pressure: f64, p_ref: f64) -> Result<f64> {
    check_pressure("pressure", pressure)?;
    check_pressure("p_ref", p_ref)?;
    Ok(20.0 * (pressure / p_ref).log10())
}

/// Inverse of [`spl`]: the RMS pressure (pascals) corresponding to a
/// sound-pressure level `level_db` measured against the standard
/// [`P_REF`] reference.
///
/// `pressure = P_REF * 10^(level_db / 20)`.
///
/// This is the exact round-trip partner of [`spl`]; `0` dB maps back to
/// [`P_REF`].
///
/// # Examples
///
/// ```
/// use valenx_acoustics::spl::{pressure_from_spl, spl};
/// let p = pressure_from_spl(94.0);
/// assert!((spl(p).unwrap() - 94.0).abs() < 1e-9);
/// ```
pub fn pressure_from_spl(level_db: f64) -> f64 {
    pressure_from_spl_ref(level_db, P_REF)
}

/// Inverse of [`spl_ref`]: the RMS pressure (pascals) corresponding to a
/// level `level_db` measured against an explicit reference `p_ref`.
///
/// `pressure = p_ref * 10^(level_db / 20)`.
///
/// `p_ref` is taken on trust here (no validation) so the function is
/// total; pass a positive reference such as [`P_REF`].
pub fn pressure_from_spl_ref(level_db: f64, p_ref: f64) -> f64 {
    p_ref * 10f64.powf(level_db / 20.0)
}

/// Combine several sound-pressure levels from mutually **incoherent**
/// sources into a single overall level.
///
/// `L_total = 10 * log10( sum_i 10^(L_i / 10) )` in decibels — the
/// decibel (mean-square pressure / power) sum.
///
/// An empty slice yields negative infinity (`log10(0)`), the level of
/// silence; callers that prefer to reject the empty case should check the
/// slice length first.
///
/// # Examples
///
/// ```
/// use valenx_acoustics::spl::combine_incoherent_levels;
/// // Two equal 60 dB sources add ~3 dB.
/// let total = combine_incoherent_levels(&[60.0, 60.0]);
/// assert!((total - 63.010_3).abs() < 1e-3);
/// ```
pub fn combine_incoherent_levels(levels_db: &[f64]) -> f64 {
    let power_sum: f64 = levels_db.iter().map(|&l| 10f64.powf(l / 10.0)).sum();
    10.0 * power_sum.log10()
}

/// Increase in level when two incoherent levels `a_db` and `b_db` are
/// summed, expressed as the *difference* above the louder of the two.
///
/// Convenience wrapper over [`combine_incoherent_levels`]; equal sources
/// give [`DOUBLE_POWER_DB`] (≈ 3.01 dB).
pub fn incoherent_sum_excess(a_db: f64, b_db: f64) -> f64 {
    let total = combine_incoherent_levels(&[a_db, b_db]);
    total - a_db.max(b_db)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    /// SPL at the reference pressure is exactly 0 dB.
    #[test]
    fn spl_at_reference_is_zero() {
        let l = spl(P_REF).unwrap();
        assert!(l.abs() < EPS, "expected 0 dB, got {l}");
    }

    /// Doubling the pressure adds ~6.02 dB (20*log10(2)).
    #[test]
    fn doubling_pressure_adds_six_db() {
        let p = 0.5; // arbitrary positive reference pressure
        let delta = spl(2.0 * p).unwrap() - spl(p).unwrap();
        assert!(
            (delta - DOUBLE_PRESSURE_DB).abs() < EPS,
            "expected {DOUBLE_PRESSURE_DB} dB, got {delta}"
        );
        // Sanity: the rule-of-thumb is ~6 dB.
        assert!((delta - 6.0).abs() < 0.05, "delta = {delta}");
    }

    /// A 10x pressure increase is exactly 20 dB.
    #[test]
    fn decade_of_pressure_is_twenty_db() {
        let delta = spl(10.0 * P_REF).unwrap() - spl(P_REF).unwrap();
        assert!((delta - 20.0).abs() < EPS, "got {delta}");
    }

    /// Two equal incoherent sources add ~3.01 dB (10*log10(2)).
    #[test]
    fn two_equal_incoherent_sources_add_three_db() {
        let total = combine_incoherent_levels(&[70.0, 70.0]);
        assert!(
            (total - (70.0 + DOUBLE_POWER_DB)).abs() < EPS,
            "got {total}"
        );
        assert!((total - 73.0103).abs() < 1e-3, "got {total}");
    }

    /// The excess of two equal sources over the louder is ~3.01 dB.
    #[test]
    fn equal_sources_excess_is_three_db() {
        let excess = incoherent_sum_excess(55.0, 55.0);
        assert!((excess - DOUBLE_POWER_DB).abs() < EPS, "got {excess}");
    }

    /// A source 10 dB below another adds well under 0.5 dB.
    #[test]
    fn ten_db_quieter_source_barely_adds() {
        let excess = incoherent_sum_excess(80.0, 70.0);
        assert!(excess > 0.0 && excess < 0.5, "got {excess}");
    }

    /// Combining a single level returns that level unchanged.
    #[test]
    fn single_level_is_identity() {
        let total = combine_incoherent_levels(&[63.7]);
        assert!((total - 63.7).abs() < EPS, "got {total}");
    }

    /// Combination is order-independent.
    #[test]
    fn combination_is_commutative() {
        let a = combine_incoherent_levels(&[40.0, 55.0, 60.0]);
        let b = combine_incoherent_levels(&[60.0, 40.0, 55.0]);
        assert!((a - b).abs() < EPS, "{a} vs {b}");
    }

    /// `pressure_from_spl` is the exact inverse of `spl`.
    #[test]
    fn pressure_round_trips_through_spl() {
        for &level in &[0.0, 20.0, 60.0, 94.0, 120.0] {
            let p = pressure_from_spl(level);
            let back = spl(p).unwrap();
            assert!((back - level).abs() < 1e-9, "level {level} -> {back}");
        }
        // 0 dB maps back to the reference pressure exactly.
        assert!((pressure_from_spl(0.0) - P_REF).abs() < EPS);
    }

    /// 94 dB SPL corresponds to ~1 Pa (the classic calibrator level).
    #[test]
    fn ninety_four_db_is_about_one_pascal() {
        let p = pressure_from_spl(94.0);
        assert!((p - 1.0).abs() < 0.01, "got {p} Pa");
    }

    /// Non-positive / non-finite pressures are rejected.
    #[test]
    fn invalid_pressures_rejected() {
        assert!(spl(0.0).is_err());
        assert!(spl(-1.0).is_err());
        assert!(spl(f64::NAN).is_err());
        assert!(spl_ref(1.0, 0.0).is_err());
        assert!(spl_ref(1.0, f64::INFINITY).is_err());
    }
}
