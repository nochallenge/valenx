//! # valenx-thermocouple
//!
//! A small, self-contained model of a **thermocouple** — the
//! thermoelectric voltage produced by two dissimilar conductors when
//! their two junctions sit at different temperatures (the Seebeck
//! effect), the **cold-junction compensation** needed to turn that
//! voltage back into an absolute temperature, and the inverse map from
//! a measured voltage to the sensed temperature.
//!
//! ## What
//!
//! A thermocouple reports the temperature *difference* between a
//! measurement junction (the "hot" junction) and a reference junction
//! (the "cold" junction). This crate provides:
//!
//! - [`Thermocouple`] — a thermocouple characterised by a single
//!   Seebeck sensitivity `S` (in volts per kelvin/degree-Celsius) plus
//!   the named [`TcType`] presets (type K, J, T, E) that pre-fill `S`
//!   with a representative near-room-temperature value.
//! - [`Thermocouple::emf`] — the open-circuit thermoelectric voltage
//!   `EMF = S * (T_hot - T_cold)` for a junction pair.
//! - [`Thermocouple::emf_compensated`] — the voltage a real
//!   instrument reads relative to a `0` reference, recovered by
//!   **cold-junction compensation**: the measured EMF across the leads
//!   plus the EMF the reference junction would itself generate against
//!   `0`.
//! - [`Thermocouple::temperature_from_emf`] — invert a measured EMF
//!   back to the hot-junction temperature given the known cold-junction
//!   temperature.
//!
//! ```
//! use valenx_thermocouple::{TcType, Thermocouple};
//!
//! let tc = Thermocouple::of_type(TcType::K);
//! // A 100 degree-C difference on a ~41 uV/C type-K gives ~4.1 mV.
//! let emf = tc.emf(125.0, 25.0).expect("valid junctions");
//! assert!((emf - 0.0041).abs() < 1e-4);
//!
//! // Recover the hot-junction temperature from that EMF.
//! let t_hot = tc.temperature_from_emf(emf, 25.0).expect("valid emf");
//! assert!((t_hot - 125.0).abs() < 1e-6);
//! ```
//!
//! ## Model
//!
//! The whole crate is the **linear Seebeck approximation**. Over a
//! modest span around the calibration point the thermoelectric voltage
//! of a junction pair is very nearly proportional to the temperature
//! difference:
//!
//! ```text
//! EMF = S * (T_hot - T_cold)
//! ```
//!
//! where `S` is the (assumed constant) Seebeck coefficient of the pair.
//! The sign convention here is that a hotter measurement junction than
//! reference junction yields a positive EMF, and the magnitude grows
//! with the temperature difference.
//!
//! **Cold-junction compensation.** A real voltmeter cannot sit at
//! absolute zero, so the raw reading only encodes `T_hot - T_cold`. To
//! recover an absolute reading the instrument adds back the voltage the
//! cold junction *would* produce against a `0` reference:
//!
//! ```text
//! V_compensated = S * (T_hot - T_cold) + S * (T_cold - 0)
//!               = S *  T_hot
//! ```
//!
//! which is exactly what [`Thermocouple::emf_compensated`] returns, and
//! is the basis for the inverse [`Thermocouple::temperature_from_emf`]:
//!
//! ```text
//! T_hot = T_cold + EMF / S
//! ```
//!
//! All temperatures are in degrees Celsius and all voltages in volts;
//! because the model is linear in temperature *difference*, the same
//! arithmetic holds unchanged if every temperature is read in kelvin.
//!
//! ## Honest scope
//!
//! This is **research/educational grade**. It implements the textbook
//! closed-form linear thermoelectric relations and nothing more; it is
//! **NOT a clinical, medical, or production engineering instrument**
//! and must not be used to make safety, diagnostic, or process-control
//! decisions.
//!
//! Specifically, a real thermocouple is **not** perfectly linear: the
//! Seebeck coefficient varies with temperature, and accurate
//! instruments use the NIST ITS-90 reference polynomials (and
//! type-specific inverse polynomials) rather than a single constant
//! `S`. The presets here use one representative near-room-temperature
//! sensitivity per type, so error grows as you move away from that
//! point. Lead resistance, amplifier offset, junction
//! inhomogeneity, drift, and noise are all out of scope.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod thermocouple;

pub use error::{ErrorCategory, ThermocoupleError};
pub use thermocouple::{TcType, Thermocouple};
