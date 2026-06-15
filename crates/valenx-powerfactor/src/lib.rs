//! # valenx-powerfactor
//!
//! Single-phase AC power-triangle arithmetic and shunt-capacitor
//! power-factor correction.
//!
//! ## What
//!
//! Given two of the quantities that describe a sinusoidal AC load —
//! RMS voltage and current, a phase angle or power factor, or the real
//! and reactive powers directly — this crate resolves the full power
//! triangle and classifies the load as leading, lagging, or unity. It
//! also sizes the shunt capacitor needed to raise a lagging load's
//! power factor to a target value.
//!
//! ## Model
//!
//! For RMS voltage `V`, RMS current `I`, and the phase angle `phi` by
//! which the current lags the voltage:
//!
//! ```text
//!   S  = V * I                 apparent power   [VA]
//!   P  = S * cos(phi)          real power       [W]
//!   Q  = S * sin(phi)          reactive power   [var]
//!   PF = P / S = cos(phi)      power factor     [-]
//! ```
//!
//! These obey the Pythagorean identity `S^2 = P^2 + Q^2`. A positive
//! `phi` (and hence `Q`) is an inductive, lagging load; a negative one
//! is a capacitive, leading load; `phi = 0` is a resistive, unity load.
//!
//! Shunt-capacitor correction moves a lagging load from angle `phi1` to
//! a smaller `phi2` (higher power factor) at constant real power by
//! supplying capacitor reactive power
//!
//! ```text
//!   Qc = P * (tan(phi1) - tan(phi2))     [var]
//! ```
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form,
//! steady-state, single-frequency sinusoidal models that assume purely
//! sinusoidal voltage and current and an ideal capacitor. They ignore
//! harmonics, transients, three-phase imbalance, capacitor tolerance,
//! resonance, and any real-world derating. This crate is NOT a
//! clinical/medical or production engineering tool — do not use it to
//! size protective devices, capacitor banks, conductors, or any
//! safety-relevant electrical installation.
//!
//! ## Surface
//!
//! - [`Phase`] — leading / lagging / unity classification.
//! - [`PowerTriangle`] — the resolved triangle with constructors
//!   [`PowerTriangle::from_vi_phase`], [`PowerTriangle::from_vi_pf`] and
//!   [`PowerTriangle::from_p_q`].
//! - [`Correction`] — shunt-capacitor sizing via
//!   [`Correction::for_target_pf`] and [`Correction::for_triangle`].
//! - [`PowerError`] — the validation error taxonomy.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod correction;
pub mod error;
pub mod triangle;

pub use correction::Correction;
pub use error::PowerError;
pub use triangle::{Phase, PowerTriangle};
