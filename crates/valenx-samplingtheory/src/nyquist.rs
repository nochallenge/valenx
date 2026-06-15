//! Nyquist-Shannon sampling criterion.
//!
//! ## Model
//!
//! A continuous real signal whose spectrum is zero outside
//! `[-fmax, +fmax]` (a strictly band-limited signal) is uniquely
//! determined by, and can be reconstructed exactly from, samples taken
//! at a uniform rate `fs` provided
//!
//! ```text
//! fs >= 2 * fmax.
//! ```
//!
//! The threshold rate `2 * fmax` is the *Nyquist rate*; half the
//! sampling rate, `fs / 2`, is the *Nyquist frequency* (also called the
//! folding frequency). A signal component strictly below the Nyquist
//! frequency is captured without ambiguity; a component at or above it
//! aliases (see [`crate::aliasing`]).
//!
//! This module exposes the criterion as a predicate
//! ([`satisfies_nyquist`]), the two derived frequencies
//! ([`nyquist_rate`], [`nyquist_frequency`]), and the oversampling
//! ratio ([`oversampling_ratio`]).
//!
//! ## Honest scope
//!
//! These are the exact textbook relations for an ideal, strictly
//! band-limited signal sampled uniformly with an ideal sampler. Real
//! signals are never perfectly band-limited and real anti-alias filters
//! have finite roll-off, so practical systems oversample (choose
//! `fs > 2 * fmax`) by a guard margin this module does not prescribe.
//! Research / educational grade only.

use crate::error::{Result, SamplingError};

/// Validate that a sampling rate is a usable positive, finite value.
///
/// Returns the rate unchanged on success. Used internally by the other
/// functions in this module and re-exported for callers that want to
/// validate a rate up front.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `fs` is not strictly positive
/// or is not finite (NaN / infinite).
pub fn validate_sample_rate(fs: f64) -> Result<f64> {
    if !fs.is_finite() {
        return Err(SamplingError::invalid(
            "sample_rate_hz",
            "must be a finite number",
        ));
    }
    if fs <= 0.0 {
        return Err(SamplingError::invalid(
            "sample_rate_hz",
            "must be strictly positive",
        ));
    }
    Ok(fs)
}

/// Validate that a frequency is non-negative and finite.
///
/// Returns the frequency unchanged on success.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `f` is negative or not finite.
pub fn validate_frequency(f: f64) -> Result<f64> {
    if !f.is_finite() {
        return Err(SamplingError::invalid(
            "frequency_hz",
            "must be a finite number",
        ));
    }
    if f < 0.0 {
        return Err(SamplingError::invalid(
            "frequency_hz",
            "must be non-negative",
        ));
    }
    Ok(f)
}

/// The Nyquist frequency (folding frequency) `fs / 2` for a sampling
/// rate `fs`.
///
/// Signal energy strictly below this frequency is represented without
/// aliasing; energy at or above it folds back into the baseband.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `fs` is not a finite positive
/// value.
pub fn nyquist_frequency(fs: f64) -> Result<f64> {
    let fs = validate_sample_rate(fs)?;
    Ok(fs / 2.0)
}

/// The Nyquist rate `2 * fmax`: the minimum sampling rate that captures
/// a signal band-limited to `fmax` without aliasing.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `fmax` is negative or not
/// finite.
pub fn nyquist_rate(fmax: f64) -> Result<f64> {
    let fmax = validate_frequency(fmax)?;
    Ok(2.0 * fmax)
}

/// Whether sampling rate `fs` satisfies the Nyquist-Shannon criterion
/// for a signal band-limited to `fmax`, i.e. whether `fs >= 2 * fmax`.
///
/// The comparison is inclusive at the exact Nyquist rate: a 100 Hz tone
/// sampled at exactly 200 Hz returns `true` (the mathematical boundary
/// of the theorem), even though practical systems leave a guard margin.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `fs` is not a finite positive
/// value, or if `fmax` is negative or not finite.
pub fn satisfies_nyquist(fs: f64, fmax: f64) -> Result<bool> {
    let fs = validate_sample_rate(fs)?;
    let rate = nyquist_rate(fmax)?;
    Ok(fs >= rate)
}

/// The oversampling ratio `fs / (2 * fmax)`.
///
/// A value of `1.0` means sampling exactly at the Nyquist rate; values
/// greater than `1.0` quantify how much headroom there is above the
/// theorem's minimum; values below `1.0` mean the signal is
/// undersampled and will alias.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `fs` is not a finite positive
/// value, or if `fmax` is not strictly positive (the ratio is undefined
/// for `fmax == 0`).
pub fn oversampling_ratio(fs: f64, fmax: f64) -> Result<f64> {
    let fs = validate_sample_rate(fs)?;
    let fmax = validate_frequency(fmax)?;
    if fmax == 0.0 {
        return Err(SamplingError::invalid(
            "fmax_hz",
            "must be strictly positive to define an oversampling ratio",
        ));
    }
    Ok(fs / (2.0 * fmax))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point comparisons throughout the tests.
    const EPS: f64 = 1e-12;

    #[test]
    fn nyquist_frequency_is_half_the_rate() {
        // Ground truth: CD audio fs = 44_100 Hz, folding freq = 22_050 Hz.
        let fnyq = nyquist_frequency(44_100.0).unwrap();
        assert!((fnyq - 22_050.0).abs() < EPS, "got {fnyq}");
    }

    #[test]
    fn nyquist_rate_is_twice_fmax() {
        // A 20 kHz audio bandwidth needs a 40 kHz Nyquist rate.
        let rate = nyquist_rate(20_000.0).unwrap();
        assert!((rate - 40_000.0).abs() < EPS, "got {rate}");
    }

    #[test]
    fn threshold_is_inclusive_at_exactly_nyquist() {
        // fs == 2*fmax: the boundary of the theorem holds (returns true).
        assert!(satisfies_nyquist(200.0, 100.0).unwrap());
    }

    #[test]
    fn just_above_threshold_satisfies() {
        assert!(satisfies_nyquist(200.000_001, 100.0).unwrap());
    }

    #[test]
    fn just_below_threshold_violates() {
        // Undersampling by a hair fails the criterion.
        assert!(!satisfies_nyquist(199.999_999, 100.0).unwrap());
    }

    #[test]
    fn clearly_oversampled_satisfies() {
        // 48 kHz comfortably captures a 20 kHz bandwidth.
        assert!(satisfies_nyquist(48_000.0, 20_000.0).unwrap());
    }

    #[test]
    fn clearly_undersampled_violates() {
        // 30 kHz cannot capture a 20 kHz bandwidth (needs >= 40 kHz).
        assert!(!satisfies_nyquist(30_000.0, 20_000.0).unwrap());
    }

    #[test]
    fn oversampling_ratio_unity_at_nyquist_rate() {
        let r = oversampling_ratio(200.0, 100.0).unwrap();
        assert!((r - 1.0).abs() < EPS, "got {r}");
    }

    #[test]
    fn oversampling_ratio_greater_than_one_when_oversampled() {
        // 48 kHz over a 20 kHz band: ratio = 48000 / 40000 = 1.2.
        let r = oversampling_ratio(48_000.0, 20_000.0).unwrap();
        assert!((r - 1.2).abs() < EPS, "got {r}");
        assert!(r > 1.0);
    }

    #[test]
    fn oversampling_ratio_less_than_one_when_undersampled() {
        // 30 kHz over a 20 kHz band: ratio = 30000 / 40000 = 0.75 < 1.
        let r = oversampling_ratio(30_000.0, 20_000.0).unwrap();
        assert!((r - 0.75).abs() < EPS, "got {r}");
        assert!(r < 1.0);
    }

    #[test]
    fn higher_fmax_needs_higher_rate_monotonic() {
        // The Nyquist rate increases monotonically with bandwidth.
        let r1 = nyquist_rate(1_000.0).unwrap();
        let r2 = nyquist_rate(2_000.0).unwrap();
        let r3 = nyquist_rate(4_000.0).unwrap();
        assert!(r1 < r2 && r2 < r3, "got {r1}, {r2}, {r3}");
    }

    #[test]
    fn dc_signal_satisfies_any_positive_rate() {
        // fmax == 0 means a DC / constant signal: any positive rate works.
        assert!(satisfies_nyquist(1.0, 0.0).unwrap());
    }

    #[test]
    fn rejects_non_positive_rate() {
        assert!(nyquist_frequency(0.0).is_err());
        assert!(nyquist_frequency(-5.0).is_err());
        assert!(satisfies_nyquist(-1.0, 10.0).is_err());
    }

    #[test]
    fn rejects_negative_frequency() {
        assert!(nyquist_rate(-1.0).is_err());
        assert!(satisfies_nyquist(100.0, -1.0).is_err());
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert!(nyquist_frequency(f64::NAN).is_err());
        assert!(nyquist_frequency(f64::INFINITY).is_err());
        assert!(nyquist_rate(f64::NAN).is_err());
    }

    #[test]
    fn oversampling_ratio_rejects_zero_fmax() {
        // Ratio is undefined for a DC signal (division by zero).
        assert!(oversampling_ratio(1_000.0, 0.0).is_err());
    }
}
