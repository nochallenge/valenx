//! # valenx-fft
//!
//! ## What
//!
//! A small, dependency-light discrete Fourier transform (DFT) library
//! over `f64` complex numbers. It provides the forward transform, the
//! inverse transform, and the frequency-domain helpers most analyses
//! need: per-bin frequencies, the Nyquist frequency, the magnitude
//! spectrum, and a Parseval energy check.
//!
//! Complex numbers are represented by a tiny in-house [`Complex`]
//! `(re, im)` pair — there is deliberately no dependency on an external
//! complex-number crate.
//!
//! ## Model
//!
//! For an `N`-point signal `x[0..N]` the forward DFT is the textbook
//! direct summation
//!
//! ```text
//! X[k] = sum_{n=0}^{N-1} x[n] * exp(-2 pi i k n / N),   k = 0..N-1
//! ```
//!
//! and the inverse DFT is
//!
//! ```text
//! x[n] = (1 / N) * sum_{k=0}^{N-1} X[k] * exp(+2 pi i k n / N).
//! ```
//!
//! This is the NumPy `numpy.fft` normalisation: the `1/N` factor sits
//! entirely on the inverse side, so `idft(dft(x))` reproduces `x`. The
//! `k`-th bin maps to the frequency `f_k = k * fs / N` for a sample
//! rate `fs`, and the Nyquist frequency is `fs / 2`. The implementation
//! is the naive `O(N^2)` algorithm — it is the literal definition, not
//! a fast (radix) FFT, and works for any length `N >= 1`.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form and
//! numerical models, intended for learning, prototyping, and
//! cross-checking. This crate is NOT a clinical / medical instrument
//! and NOT a production digital-signal-processing engine: the `O(N^2)`
//! transform is too slow for large signals, no windowing /
//! anti-aliasing / detrending is applied for you, and the numerics are
//! plain `f64` with no error-bound guarantees. Do not rely on it for
//! diagnostic, safety-critical, or regulated use.
//!
//! ## Quick example
//!
//! ```
//! use valenx_fft::{dft_real, idft, magnitude_spectrum, nyquist_frequency, Complex};
//!
//! // A unit impulse has a perfectly flat spectrum.
//! let spectrum = dft_real(&[1.0, 0.0, 0.0, 0.0]).unwrap();
//! for m in magnitude_spectrum(&spectrum) {
//!     assert!((m - 1.0).abs() < 1e-12);
//! }
//!
//! // The forward transform inverts back to the original samples.
//! let recovered = idft(&spectrum).unwrap();
//! assert!((recovered[0].re - 1.0).abs() < 1e-12);
//!
//! // Nyquist is half the sample rate.
//! assert!((nyquist_frequency(48_000.0).unwrap() - 24_000.0).abs() < 1e-9);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod complex;
pub mod dft;
pub mod error;

pub use complex::Complex;
pub use dft::{
    bin_frequencies, bin_frequency, dft, dft_real, idft, magnitude_spectrum, nyquist_frequency,
    spectral_energy, time_energy,
};
pub use error::{ErrorCategory, FftError};
