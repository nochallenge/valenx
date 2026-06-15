//! Uniform quantization and the ideal-quantizer SNR.
//!
//! ## Model
//!
//! An ideal uniform quantizer with `N` bits divides a full-scale input
//! range `R` into `2^N` equal steps, so the quantization step size
//! (the value of one least-significant bit) is
//!
//! ```text
//! q = R / 2^N.
//! ```
//!
//! Modelling the quantization error as a uniform random variable on
//! `[-q/2, +q/2]` (the standard additive-noise approximation) and the
//! input as a full-scale sinusoid gives the classic ideal
//! signal-to-quantization-noise ratio
//!
//! ```text
//! SNR_dB = 6.02 * N + 1.76.
//! ```
//!
//! Each extra bit adds about `6.02` dB (a factor of two in amplitude
//! resolution); the `1.76` dB offset is `10 * log10(3 / 2)`, the gain
//! of the full-scale-sinusoid case over the raw step quantization.
//!
//! This module exposes the step size ([`quant_step`]), the ideal SNR
//! ([`ideal_snr_db`]), the number of effective bits implied by a
//! measured SNR ([`enob_from_snr_db`]), a single-sample mid-tread
//! quantizer ([`quantize`]), and a small serializable
//! [`QuantizationAnalysis`] bundle.
//!
//! ## Honest scope
//!
//! `6.02 N + 1.76` is the textbook ideal-quantizer figure: it assumes a
//! perfect uniform quantizer, a full-scale sinusoidal input, and the
//! white-noise quantization-error model. Real converters fall short of
//! it (their datasheets quote a lower *effective number of bits*) and
//! this module models none of those non-idealities. Research /
//! educational grade only.

use crate::error::{Result, SamplingError};
use serde::{Deserialize, Serialize};

/// Decibels of SNR gained per quantizer bit in the ideal model.
pub const DB_PER_BIT: f64 = 6.02;

/// The constant dB offset in the ideal full-scale-sinusoid SNR formula,
/// equal to `10 * log10(3 / 2) ≈ 1.7609`.
pub const SNR_OFFSET_DB: f64 = 1.76;

/// Validate a quantizer bit depth: it must be at least one bit and small
/// enough that `2^bits` is representable without overflow.
///
/// Returns the bit depth unchanged on success.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `bits == 0` or if `bits` is so
/// large that `2u64.pow(bits)` would overflow (`bits >= 64`).
pub fn validate_bits(bits: u32) -> Result<u32> {
    if bits == 0 {
        return Err(SamplingError::invalid("bits", "must be at least 1"));
    }
    if bits >= 64 {
        return Err(SamplingError::invalid(
            "bits",
            "must be < 64 so that 2^bits is representable",
        ));
    }
    Ok(bits)
}

/// Validate that a full-scale range is strictly positive and finite.
///
/// Returns the range unchanged on success.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `range` is not strictly
/// positive or is not finite.
pub fn validate_range(range: f64) -> Result<f64> {
    if !range.is_finite() {
        return Err(SamplingError::invalid("range", "must be a finite number"));
    }
    if range <= 0.0 {
        return Err(SamplingError::invalid("range", "must be strictly positive"));
    }
    Ok(range)
}

/// The number of quantization levels `2^bits` for a `bits`-bit
/// quantizer, returned as an exact integer.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `bits` is zero or `>= 64`.
pub fn levels(bits: u32) -> Result<u64> {
    let bits = validate_bits(bits)?;
    Ok(1u64 << bits)
}

/// The quantization step size `q = range / 2^bits` (the value of one
/// least-significant bit) for a `bits`-bit quantizer spanning a
/// full-scale input `range`.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `range` is not a finite
/// positive value, or if `bits` is zero or `>= 64`.
pub fn quant_step(range: f64, bits: u32) -> Result<f64> {
    let range = validate_range(range)?;
    let n_levels = levels(bits)?;
    Ok(range / n_levels as f64)
}

/// The ideal signal-to-quantization-noise ratio in decibels,
/// `SNR_dB = 6.02 * bits + 1.76`, for a full-scale sinusoid into an
/// ideal `bits`-bit quantizer.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `bits` is zero or `>= 64`.
pub fn ideal_snr_db(bits: u32) -> Result<f64> {
    let bits = validate_bits(bits)?;
    Ok(DB_PER_BIT * bits as f64 + SNR_OFFSET_DB)
}

/// The effective number of bits (ENOB) implied by a measured SNR, the
/// inverse of [`ideal_snr_db`]:
///
/// ```text
/// ENOB = (SNR_dB - 1.76) / 6.02.
/// ```
///
/// The result is real-valued (not rounded) and may be fractional or, for
/// an SNR below the offset, negative; callers decide how to interpret
/// it.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `snr_db` is not finite.
pub fn enob_from_snr_db(snr_db: f64) -> Result<f64> {
    if !snr_db.is_finite() {
        return Err(SamplingError::invalid("snr_db", "must be a finite number"));
    }
    Ok((snr_db - SNR_OFFSET_DB) / DB_PER_BIT)
}

/// Quantize a single sample with an ideal uniform mid-tread quantizer.
///
/// The input `sample` must lie within the symmetric full-scale interval
/// `[-range/2, +range/2]`. The returned value is `sample` rounded to the
/// nearest multiple of the step size `q = range / 2^bits`, i.e. the
/// reconstructed analog value the quantizer would emit.
///
/// # Errors
///
/// Returns [`SamplingError::Invalid`] if `range` is not a finite
/// positive value, if `bits` is zero or `>= 64`, or if `sample` is not
/// finite; returns [`SamplingError::OutOfRange`] if `sample` falls
/// outside `[-range/2, +range/2]`.
pub fn quantize(sample: f64, range: f64, bits: u32) -> Result<f64> {
    let q = quant_step(range, bits)?;
    if !sample.is_finite() {
        return Err(SamplingError::invalid("sample", "must be a finite number"));
    }
    let half = range / 2.0;
    if sample < -half || sample > half {
        return Err(SamplingError::out_of_range("sample", sample, -half, half));
    }
    Ok((sample / q).round() * q)
}

/// A bundle of the standard quantizer figures for a given bit depth and
/// full-scale range, convenient to compute once and serialize.
///
/// All fields are derived from `bits` and `range_full_scale` via the
/// formulas in this module; see [`QuantizationAnalysis::new`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct QuantizationAnalysis {
    /// Quantizer resolution in bits.
    pub bits: u32,
    /// Full-scale input range `R`.
    pub range_full_scale: f64,
    /// Number of quantization levels `2^bits`.
    pub levels: u64,
    /// Step size `q = R / 2^bits`.
    pub step: f64,
    /// Ideal SNR in decibels, `6.02 * bits + 1.76`.
    pub ideal_snr_db: f64,
}

impl QuantizationAnalysis {
    /// Compute the analysis bundle for a `bits`-bit quantizer over a
    /// full-scale input `range`.
    ///
    /// # Errors
    ///
    /// Returns [`SamplingError::Invalid`] if `range` is not a finite
    /// positive value or if `bits` is zero or `>= 64`.
    pub fn new(range: f64, bits: u32) -> Result<Self> {
        Ok(Self {
            bits,
            range_full_scale: validate_range(range)?,
            levels: levels(bits)?,
            step: quant_step(range, bits)?,
            ideal_snr_db: ideal_snr_db(bits)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point comparisons throughout the tests.
    const EPS: f64 = 1e-9;

    #[test]
    fn levels_is_two_to_the_bits() {
        assert_eq!(levels(1).unwrap(), 2);
        assert_eq!(levels(8).unwrap(), 256);
        assert_eq!(levels(16).unwrap(), 65_536);
    }

    #[test]
    fn step_is_range_over_levels() {
        // 8-bit converter over a 1.0 V range: q = 1 / 256.
        let q = quant_step(1.0, 8).unwrap();
        assert!((q - 1.0 / 256.0).abs() < EPS, "got {q}");
    }

    #[test]
    fn step_ground_truth_unit_range_per_bit_depth() {
        // 10 V full-scale, 12-bit: q = 10 / 4096 ≈ 2.4414 mV.
        let q = quant_step(10.0, 12).unwrap();
        assert!((q - 10.0 / 4096.0).abs() < EPS, "got {q}");
        assert!((q - 0.002_441_406_25).abs() < EPS, "got {q}");
    }

    #[test]
    fn more_bits_means_smaller_step() {
        // Halving the step per added bit is the defining property.
        let q8 = quant_step(1.0, 8).unwrap();
        let q9 = quant_step(1.0, 9).unwrap();
        assert!(q9 < q8, "q8={q8}, q9={q9}");
        assert!((q8 / q9 - 2.0).abs() < EPS, "ratio {}", q8 / q9);
    }

    #[test]
    fn ideal_snr_matches_formula_at_landmark_depths() {
        // 8-bit: 6.02*8 + 1.76 = 49.92 dB.
        let s8 = ideal_snr_db(8).unwrap();
        assert!((s8 - 49.92).abs() < EPS, "got {s8}");

        // 16-bit (CD audio): 6.02*16 + 1.76 = 98.08 dB.
        let s16 = ideal_snr_db(16).unwrap();
        assert!((s16 - 98.08).abs() < EPS, "got {s16}");

        // 1-bit: 6.02 + 1.76 = 7.78 dB.
        let s1 = ideal_snr_db(1).unwrap();
        assert!((s1 - 7.78).abs() < EPS, "got {s1}");
    }

    #[test]
    fn more_bits_means_more_snr() {
        // SNR is strictly increasing in bit depth.
        let mut prev = ideal_snr_db(1).unwrap();
        for b in 2..=24u32 {
            let cur = ideal_snr_db(b).unwrap();
            assert!(cur > prev, "snr({b})={cur} not > {prev}");
            prev = cur;
        }
    }

    #[test]
    fn each_added_bit_adds_six_point_oh_two_db() {
        // The incremental SNR gain per bit is exactly DB_PER_BIT.
        let a = ideal_snr_db(10).unwrap();
        let b = ideal_snr_db(11).unwrap();
        assert!((b - a - DB_PER_BIT).abs() < EPS, "delta {}", b - a);
    }

    #[test]
    fn enob_inverts_ideal_snr() {
        // Round-trip: ENOB(SNR(N)) == N for integer N.
        for n in 1..=20u32 {
            let snr = ideal_snr_db(n).unwrap();
            let enob = enob_from_snr_db(snr).unwrap();
            assert!((enob - n as f64).abs() < EPS, "n={n}, enob={enob}");
        }
    }

    #[test]
    fn enob_of_offset_is_zero() {
        // An SNR equal to the bare offset implies zero effective bits.
        let enob = enob_from_snr_db(SNR_OFFSET_DB).unwrap();
        assert!(enob.abs() < EPS, "got {enob}");
    }

    #[test]
    fn quantize_snaps_to_nearest_step() {
        // 2-bit over [-1, 1]: q = 2/4 = 0.5; levels at ..., -0.5, 0, 0.5,
        // 1.0. 0.3 rounds to 0.5; 0.2 rounds to 0.0.
        let q_up = quantize(0.3, 2.0, 2).unwrap();
        assert!((q_up - 0.5).abs() < EPS, "got {q_up}");

        let q_dn = quantize(0.2, 2.0, 2).unwrap();
        assert!(q_dn.abs() < EPS, "got {q_dn}");
    }

    #[test]
    fn quantize_is_exact_on_a_grid_point() {
        // A sample already on a step boundary is unchanged.
        let v = quantize(0.5, 2.0, 2).unwrap();
        assert!((v - 0.5).abs() < EPS, "got {v}");
    }

    #[test]
    fn quantize_error_never_exceeds_half_step() {
        // The defining bound of mid-tread quantization: |error| <= q/2.
        let range = 2.0;
        let bits = 6;
        let q = quant_step(range, bits).unwrap();
        let half = range / 2.0;
        let mut x = -half;
        while x <= half {
            let y = quantize(x, range, bits).unwrap();
            assert!(
                (y - x).abs() <= q / 2.0 + EPS,
                "x={x} quantized to {y}, error {} exceeds q/2={}",
                (y - x).abs(),
                q / 2.0
            );
            x += q / 17.0;
        }
    }

    #[test]
    fn quantize_rejects_out_of_range_sample() {
        // 1.5 lies outside the [-1, 1] full-scale interval.
        let err = quantize(1.5, 2.0, 8).unwrap_err();
        assert_eq!(err.code(), "samplingtheory.out_of_range");
    }

    #[test]
    fn analysis_bundle_is_internally_consistent() {
        let a = QuantizationAnalysis::new(2.0, 10).unwrap();
        assert_eq!(a.bits, 10);
        assert_eq!(a.levels, 1024);
        assert!((a.step - 2.0 / 1024.0).abs() < EPS, "got {}", a.step);
        assert!(
            (a.ideal_snr_db - (DB_PER_BIT * 10.0 + SNR_OFFSET_DB)).abs() < EPS,
            "got {}",
            a.ideal_snr_db
        );
    }

    #[test]
    fn rejects_zero_bits_and_bad_range() {
        assert!(quant_step(1.0, 0).is_err());
        assert!(ideal_snr_db(0).is_err());
        assert!(levels(64).is_err());
        assert!(quant_step(0.0, 8).is_err());
        assert!(quant_step(-1.0, 8).is_err());
        assert!(quant_step(f64::INFINITY, 8).is_err());
    }
}
