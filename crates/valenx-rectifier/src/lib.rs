//! # valenx-rectifier
//!
//! Closed-form analysis of single-phase diode rectifiers.
//!
//! ## What
//!
//! Pure-function helpers for the textbook figures of an ideal half-wave
//! and full-wave rectifier driven by a sinusoidal mains:
//!
//! - average (DC) output [`rectifier::half_wave_vdc`] = `Vpeak/pi` and
//!   [`rectifier::full_wave_vdc`] = `2*Vpeak/pi`,
//! - RMS output [`rectifier::half_wave_vrms`] = `Vpeak/2` and
//!   [`rectifier::full_wave_vrms`] = `Vpeak/sqrt(2)`,
//! - the dimensionless [`rectifier::ripple_factor`]
//!   `r = sqrt((Vrms/Vdc)^2 - 1)`,
//! - the capacitor-input-filter peak-to-peak ripple
//!   [`rectifier::capacitor_ripple_pp`] `Vr = I/(f*C)`.
//!
//! All entry points validate their inputs and return [`RectifierError`]
//! rather than panicking or producing `NaN`.
//!
//! ## Model
//!
//! The input is the ideal sinusoid `v(t) = Vpeak * sin(omega t)`. Diodes
//! are ideal switches (zero forward voltage drop, zero reverse leakage,
//! instantaneous turn-on). Averages and RMS values come from analytic
//! integration over one mains period; the ripple factor is their exact
//! algebraic combination. The capacitor-filter formula `Vr = I/(f*C)` is
//! the standard linear-discharge approximation: a constant-current load
//! `I` discharges the reservoir capacitor `C` between conduction pulses
//! arriving at frequency `f` (the mains frequency for a half-wave
//! rectifier, twice the mains frequency for a full-wave one).
//!
//! ## Honest scope
//!
//! Research/educational grade. These are first-order, idealized
//! textbook formulae: real diodes drop ~0.6 to 1.0 V, transformers and
//! wiring add series resistance, the constant-current load assumption is
//! only approximate, and the linear-discharge ripple estimate ignores
//! the finite conduction angle and capacitor ESR. The crate is for
//! learning and quick hand-calculation cross-checks, **not** a
//! clinical/medical/production engineering tool. Do not size real
//! power-supply components from these numbers without proper derating
//! and SPICE/bench verification.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod rectifier;

pub use error::{ErrorCategory, RectifierError};
pub use rectifier::{
    capacitor_ripple_pp, capacitor_ripple_pp_for, full_wave_vdc, full_wave_vrms, half_wave_vdc,
    half_wave_vrms, ripple_factor, ripple_factor_from, vdc, vrms, Topology,
};
