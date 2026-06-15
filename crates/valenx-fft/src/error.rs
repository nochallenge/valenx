//! Error taxonomy for the discrete Fourier transform routines.

use thiserror::Error;

/// Errors raised while validating DFT inputs.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum FftError {
    /// The input signal had length zero. A DFT is only defined for at
    /// least one sample (`N >= 1`).
    #[error("empty signal: a DFT needs at least one sample (N >= 1)")]
    EmptySignal,

    /// The sample rate was not a finite, strictly positive value.
    /// Stored as a human-readable string because `f64` is not `Eq`.
    #[error("invalid sample rate `{0}`: must be finite and > 0")]
    InvalidSampleRate(String),

    /// A requested bin index `k` was outside the valid range `0..N`.
    #[error("bin index {index} out of range for a length-{len} spectrum (valid 0..{len})")]
    BinOutOfRange {
        /// The offending bin index.
        index: usize,
        /// The spectrum length `N`.
        len: usize,
    },
}

/// Coarse classification of a [`FftError`], useful for routing /
/// telemetry without matching every variant.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied bad data (empty signal, out-of-range bin).
    Input,
    /// A tunable parameter was invalid (sample rate).
    Config,
}

impl FftError {
    /// Stable, kebab-cased identifier for this error, suitable for logs
    /// and assertions that should not depend on the `Display` wording.
    pub fn code(&self) -> &'static str {
        match self {
            FftError::EmptySignal => "fft.empty_signal",
            FftError::InvalidSampleRate(_) => "fft.invalid_sample_rate",
            FftError::BinOutOfRange { .. } => "fft.bin_out_of_range",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            FftError::EmptySignal => ErrorCategory::Input,
            FftError::BinOutOfRange { .. } => ErrorCategory::Input,
            FftError::InvalidSampleRate(_) => ErrorCategory::Config,
        }
    }
}

/// Validate that a sample rate is finite and strictly positive,
/// returning it unchanged on success.
///
/// # Errors
///
/// Returns [`FftError::InvalidSampleRate`] when `fs` is `NaN`,
/// infinite, zero, or negative.
pub(crate) fn validate_sample_rate(fs: f64) -> Result<f64, FftError> {
    if fs.is_finite() && fs > 0.0 {
        Ok(fs)
    } else {
        Err(FftError::InvalidSampleRate(format!("{fs}")))
    }
}
