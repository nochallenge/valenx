//! # valenx-samplingtheory — DSP sampling theory
//!
//! Closed-form models for the four cornerstone results of
//! digital-signal-processing sampling: the Nyquist-Shannon criterion,
//! aliasing / spectral folding, uniform quantization step size, and the
//! ideal-quantizer signal-to-noise ratio. Everything here is a small,
//! exact, well-validated formula from a standard DSP textbook
//! (Oppenheim & Schafer, *Discrete-Time Signal Processing*; Proakis &
//! Manolakis, *Digital Signal Processing*).
//!
//! The crate is pure: it computes numbers from numbers, performs no I/O,
//! and pulls in no heavy or platform-specific dependencies.
//!
//! ## What
//!
//! - **Nyquist criterion** ([`nyquist`]) — the Nyquist frequency
//!   `fs / 2` ([`nyquist::nyquist_frequency`]), the Nyquist rate
//!   `2 * fmax` ([`nyquist::nyquist_rate`]), the predicate
//!   `fs >= 2 * fmax` ([`nyquist::satisfies_nyquist`]), and the
//!   oversampling ratio ([`nyquist::oversampling_ratio`]).
//! - **Aliasing** ([`aliasing`]) — the folded apparent frequency
//!   `|f - round(f / fs) * fs|` ([`aliasing::alias_frequency`]) and a
//!   test for whether folding occurs ([`aliasing::is_aliased`]).
//! - **Quantization** ([`quantization`]) — the step size
//!   `q = range / 2^bits` ([`quantization::quant_step`]), the ideal SNR
//!   `6.02 * N + 1.76` ([`quantization::ideal_snr_db`]), its inverse
//!   ([`quantization::enob_from_snr_db`]), a single-sample mid-tread
//!   quantizer ([`quantization::quantize`]), and a serializable
//!   [`quantization::QuantizationAnalysis`] bundle.
//! - **Errors** ([`error`]) — a validated [`error::SamplingError`] enum
//!   with stable [`code`](error::SamplingError::code) and
//!   [`category`](error::SamplingError::category) accessors.
//!
//! ## Model
//!
//! The relations implemented, for a real signal band-limited to `fmax`
//! sampled uniformly at rate `fs` and quantized with `N` bits over a
//! full-scale range `R`:
//!
//! ```text
//! Nyquist criterion : fs >= 2 * fmax
//! Nyquist frequency : f_nyq = fs / 2
//! aliased frequency : f_alias = | f - round(f / fs) * fs |
//! quantization step : q = R / 2^N
//! ideal SNR (dB)    : SNR_dB = 6.02 * N + 1.76
//! ```
//!
//! Each is the idealized closed form: a strictly band-limited signal, an
//! ideal uniform sampler, an ideal uniform quantizer, the full-scale
//! sinusoid and white-noise quantization-error assumptions behind
//! `6.02 N + 1.76`.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form and
//! numerical models, not measurements of a real system. Physical signals
//! are never perfectly band-limited, real anti-alias filters have finite
//! roll-off, and real converters achieve fewer effective bits than the
//! `6.02 N + 1.76` ideal. This crate is NOT a clinical, medical, or
//! production signal-processing / instrumentation tool, and its outputs
//! must not be relied upon for safety-critical, diagnostic, or
//! certification purposes.
//!
//! ## Example
//!
//! ```
//! use valenx_samplingtheory::nyquist::satisfies_nyquist;
//! use valenx_samplingtheory::aliasing::alias_frequency;
//! use valenx_samplingtheory::quantization::{ideal_snr_db, quant_step};
//!
//! // 48 kHz comfortably captures a 20 kHz audio bandwidth.
//! assert!(satisfies_nyquist(48_000.0, 20_000.0).unwrap());
//!
//! // A 6 kHz tone sampled at 8 kHz aliases down to 2 kHz.
//! let a = alias_frequency(6_000.0, 8_000.0).unwrap();
//! assert!((a - 2_000.0).abs() < 1e-6);
//!
//! // A 16-bit converter over a 2.0 V range: q and the ideal SNR.
//! let q = quant_step(2.0, 16).unwrap();
//! assert!((q - 2.0 / 65_536.0).abs() < 1e-12);
//! assert!((ideal_snr_db(16).unwrap() - 98.08).abs() < 1e-9);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod aliasing;
pub mod error;
pub mod nyquist;
pub mod quantization;

pub use error::{ErrorCategory, Result, SamplingError};
