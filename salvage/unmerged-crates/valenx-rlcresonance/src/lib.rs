//! # valenx-rlcresonance
//!
//! Closed-form resonance, quality-factor, bandwidth and transient-damping
//! models for **series** and **parallel** RLC circuits, in SI units.
//!
//! ## What
//!
//! Given an inductance `L` (henries), a capacitance `C` (farads) and a
//! resistance `R` (ohms), this crate computes the small set of numbers an
//! engineer reaches for when reasoning about a resonant tank or a
//! second-order filter: the resonant frequency, the loaded quality
//! factor, the `-3 dB` bandwidth, and the transient damping regime
//! (underdamped, critically damped, overdamped). Every public function
//! validates its inputs and returns a [`Result`], so a non-finite or
//! out-of-domain argument is rejected up front instead of silently
//! producing `NaN`/`Inf`.
//!
//! ## Model
//!
//! The governing relations (standard textbook AC-circuit and
//! second-order-ODE results) are
//!
//! ```text
//! f0     = 1 / (2 pi sqrt(L C))        resonant frequency        (Hz)
//! omega0 = 1 / sqrt(L C) = 2 pi f0     resonant angular freq.    (rad/s)
//! Q_s    = (1 / R) sqrt(L / C)         series quality factor     (-)
//! Q_p    = R sqrt(C / L)               parallel quality factor   (-)
//! BW     = f0 / Q = R / (2 pi L)       -3 dB bandwidth           (Hz)
//! alpha  = R / (2 L)                   neper / damping freq.     (rad/s)
//! zeta   = alpha / omega0 = 1 / (2 Q)  damping ratio             (-)
//! R_crit = 2 sqrt(L / C)               critical resistance       (ohm)
//! omega_d = omega0 sqrt(1 - zeta^2)    damped ringing freq.      (rad/s)
//! ```
//!
//! The headline ground-truth identity `Q * BW = f0` is exact (it is the
//! definition `BW = f0 / Q` rearranged) and the test-suite pins it,
//! together with hand-worked textbook numbers for each formula and the
//! known limiting cases (lossless `R = 0`, the `zeta = 1` critical
//! boundary, the undamped `omega_d = omega0` limit).
//!
//! ## Honest scope
//!
//! Research/educational grade: these are standard textbook closed-form
//! models, validated against analytic ground truth. This is **NOT** a
//! clinical/medical or production-certified engineering tool. It models
//! the ideal lumped `R`/`L`/`C` only — there are no component tolerances,
//! safety factors, fatigue/thermal limits, parasitics (lead inductance,
//! ESR, self-resonance), core saturation, skin/proximity effect, or any
//! code-compliance checks. Do not size real hardware from these numbers.
//!
//! ## Example
//!
//! ```
//! use valenx_rlcresonance::{
//!     resonant_frequency_hz, series_quality_factor, series_bandwidth_hz,
//!     classify_regime, DampingRegime,
//! };
//!
//! // A series tank: R = 10 ohm, L = 1 mH, C = 1 uF.
//! let (r, l, c) = (10.0, 1e-3, 1e-6);
//!
//! let f0 = resonant_frequency_hz(l, c).unwrap();
//! let q = series_quality_factor(r, l, c).unwrap();
//! let bw = series_bandwidth_hz(r, l).unwrap();
//!
//! // Ground-truth identity: Q * BW == f0.
//! assert!((q * bw - f0).abs() < 1e-6);
//!
//! // Light damping -> the circuit rings (underdamped).
//! assert_eq!(classify_regime(r, l, c).unwrap(), DampingRegime::Underdamped);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod damping;
pub mod error;
pub mod resonance;

pub use damping::{
    classify_regime, critical_resistance_ohm, damped_frequency_rad_s, damping_ratio,
    neper_frequency_rad_s, DampingRegime,
};
pub use error::{require_non_negative, require_positive, Result, RlcError};
pub use resonance::{
    bandwidth_from_f0_quality, parallel_quality_factor, quality_from_f0_bandwidth,
    resonant_angular_frequency_rad_s, resonant_frequency_hz, series_bandwidth_hz,
    series_quality_factor,
};
