//! Aliasing and frequency folding.
//!
//! ## Model
//!
//! When a continuous tone of frequency `f` is sampled at rate `fs`, the
//! sampled sequence is indistinguishable from a tone at the *aliased*
//! (apparent) frequency
//!
//! ```text
//! f_alias = | f - round(f / fs) * fs |,
//! ```
//!
//! where `round` is round-half-away-from-zero. This maps any input
//! frequency onto the baseband interval `[0, fs / 2]`: it subtracts the
//! nearest integer multiple of `fs` and takes the magnitude, which is
//! exactly the textbook spectral-folding picture. A component already
//! below the Nyquist frequency `fs / 2` is returned unchanged; a
//! component above it folds back down.
//!
//! [`alias_frequency`] computes the apparent frequency; [`is_aliased`]
//! reports whether folding actually moved the component (i.e. whether
//! the original tone was at or above the Nyquist frequency).
//!
//! ## Honest scope
//!
//! This is the ideal point-spectrum folding rule for a single real
//! sinusoid under uniform ideal sampling. It says nothing about the
//! amplitude of the alias, anti-alias-filter attenuation, or
//! wide-band / random signals. Research / educational grade only.

use crate::error::Result;
use crate::nyquist::{validate_frequency, validate_sample_rate};

/// The aliased (apparent) frequency of a tone at `f` hertz sampled at
/// `fs` hertz.
///
/// Implements `f_alias = |f - round(f / fs) * fs|`, folding `f` onto the
/// baseband interval `[0, fs / 2]`. A tone strictly below the Nyquist
/// frequency is returned essentially unchanged (to within floating-point
/// rounding); a tone above it is folded back down.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`](crate::error::SamplingError) if
/// `fs` is not a finite positive value, or if `f` is negative or not
/// finite.
pub fn alias_frequency(f: f64, fs: f64) -> Result<f64> {
    let fs = validate_sample_rate(fs)?;
    let f = validate_frequency(f)?;
    let folded = (f - (f / fs).round() * fs).abs();
    Ok(folded)
}

/// Whether sampling a tone at `f` hertz at rate `fs` produces aliasing,
/// i.e. whether `f` is at or above the Nyquist frequency `fs / 2`.
///
/// Returns `false` for a tone strictly below the Nyquist frequency
/// (captured faithfully) and `true` once `f >= fs / 2` (the component
/// folds to a different apparent frequency, or sits exactly on the
/// folding frequency).
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`](crate::error::SamplingError) if
/// `fs` is not a finite positive value, or if `f` is negative or not
/// finite.
pub fn is_aliased(f: f64, fs: f64) -> Result<bool> {
    let fs = validate_sample_rate(fs)?;
    let f = validate_frequency(f)?;
    Ok(f >= fs / 2.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point comparisons throughout the tests.
    const EPS: f64 = 1e-9;

    #[test]
    fn below_nyquist_returns_input_unchanged() {
        // 30 Hz tone, 100 Hz sampling: 30 < 50, captured faithfully.
        let a = alias_frequency(30.0, 100.0).unwrap();
        assert!((a - 30.0).abs() < EPS, "got {a}");
        assert!(!is_aliased(30.0, 100.0).unwrap());
    }

    #[test]
    fn above_nyquist_folds_back() {
        // 70 Hz tone, 100 Hz sampling: folds to |70 - 100| = 30 Hz.
        let a = alias_frequency(70.0, 100.0).unwrap();
        assert!((a - 30.0).abs() < EPS, "got {a}");
        assert!(is_aliased(70.0, 100.0).unwrap());
    }

    #[test]
    fn classic_textbook_case_folds() {
        // Oppenheim & Schafer staple: a 6 kHz tone sampled at 8 kHz
        // aliases to |6000 - 8000| = 2000 Hz.
        let a = alias_frequency(6_000.0, 8_000.0).unwrap();
        assert!((a - 2_000.0).abs() < EPS, "got {a}");
    }

    #[test]
    fn frequency_above_sampling_rate_folds_into_baseband() {
        // 120 Hz tone, 100 Hz sampling: round(1.2) = 1, |120 - 100| = 20.
        let a = alias_frequency(120.0, 100.0).unwrap();
        assert!((a - 20.0).abs() < EPS, "got {a}");
        // The fold always lands in [0, fs/2].
        assert!(a <= 100.0 / 2.0 + EPS, "got {a}");
    }

    #[test]
    fn far_above_sampling_rate_still_in_baseband() {
        // 250 Hz tone, 100 Hz sampling: round(2.5) = 3 (half away from
        // zero), |250 - 300| = 50 Hz, which is the Nyquist frequency.
        let a = alias_frequency(250.0, 100.0).unwrap();
        assert!((a - 50.0).abs() < EPS, "got {a}");
        assert!(a <= 100.0 / 2.0 + EPS, "got {a}");
    }

    #[test]
    fn integer_multiple_of_fs_aliases_to_dc() {
        // A tone exactly at fs (and at 2*fs) folds to 0 Hz.
        let a1 = alias_frequency(100.0, 100.0).unwrap();
        let a2 = alias_frequency(200.0, 100.0).unwrap();
        assert!(a1.abs() < EPS, "got {a1}");
        assert!(a2.abs() < EPS, "got {a2}");
    }

    #[test]
    fn exactly_at_nyquist_stays_at_nyquist() {
        // f == fs/2: folds to itself and counts as aliased.
        let a = alias_frequency(50.0, 100.0).unwrap();
        assert!((a - 50.0).abs() < EPS, "got {a}");
        assert!(is_aliased(50.0, 100.0).unwrap());
    }

    #[test]
    fn fold_never_exceeds_nyquist_frequency_for_a_sweep() {
        // Sweep f across several periods of fs; the alias must always
        // remain within the baseband [0, fs/2].
        let fs = 1_000.0;
        let fnyq = fs / 2.0;
        let mut f = 0.0;
        while f <= 5_000.0 {
            let a = alias_frequency(f, fs).unwrap();
            assert!(
                a >= -EPS && a <= fnyq + EPS,
                "f={f} aliased to {a}, outside [0, {fnyq}]"
            );
            f += 7.5;
        }
    }

    #[test]
    fn dc_is_never_aliased() {
        assert!(!is_aliased(0.0, 100.0).unwrap());
        let a = alias_frequency(0.0, 100.0).unwrap();
        assert!(a.abs() < EPS, "got {a}");
    }

    #[test]
    fn rejects_bad_inputs() {
        assert!(alias_frequency(10.0, 0.0).is_err());
        assert!(alias_frequency(-10.0, 100.0).is_err());
        assert!(alias_frequency(f64::NAN, 100.0).is_err());
        assert!(is_aliased(10.0, -1.0).is_err());
    }
}
