//! # valenx-zener
//!
//! Closed-form analysis of the textbook zener-diode shunt voltage
//! regulator: size the series resistor, find the regulated output, the
//! diode and resistor currents, the diode power, and the minimum supply
//! voltage that keeps the diode in breakdown.
//!
//! ## What
//!
//! A zener regulator drops the difference between an unregulated supply
//! and a fixed output across a series resistor `Rs`, while a reverse-
//! biased zener diode in parallel with the load pins the output node at
//! its breakdown voltage `Vz`. This crate turns the standard hand-
//! analysis equations into validated functions:
//!
//! [`size_series_resistor`] solves `Rs = (Vin - Vz) / (Iz + Il)`;
//! [`ZenerRegulator`] then evaluates a chosen design point —
//! [`ZenerRegulator::resistor_current_a`] gives `I_Rs = (Vin - Vz)/Rs`,
//! [`ZenerRegulator::output_voltage_v`] clamps the output at `Vz`,
//! [`ZenerRegulator::zener_current_a`] applies `Iz = I_Rs - I_load`,
//! [`ZenerRegulator::zener_power_w`] applies `Pz = Vz * Iz`, and
//! [`ZenerRegulator::min_input_voltage_v`] gives the lowest supply that
//! still keeps the diode conducting. [`ZenerRegulator::operating_point`]
//! bundles the whole state into one [`OperatingPoint`].
//!
//! ## Model
//!
//! The diode is the ideal piecewise model: open below `Vz`, and a
//! perfect constant-voltage source of exactly `Vz` once in reverse
//! breakdown. Consequences of that idealisation:
//!
//! Output regulation is treated as a hard clamp at `Vz` (zero dynamic
//! impedance `r_z`, so line- and load-regulation ripple are not
//! modelled). The series resistor is ideal and ohmic. Kirchhoff's
//! current law splits the resistor current between diode and load,
//! `I_Rs = Iz + Il`, and energy balances as `Vin·I_Rs = P_Rs + Pz + P_load`
//! in regulation. The diode is "in regulation" exactly when the supply
//! has headroom (`Vin > Vz`) and the load has not stolen the entire
//! resistor current (`Iz > 0`); outside that window the output droops
//! and the reported signed `Iz` goes non-positive.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form / ideal-
//! diode models — constant `Vz`, no dynamic impedance, no temperature
//! coefficient, no leakage or knee behaviour near breakdown, no part
//! tolerances, and no resistor power- or diode power-rating checks
//! against a real datasheet. This crate is NOT a clinical, medical, or
//! production engineering tool; do not use it to qualify a power supply
//! or to set component ratings for hardware that ships. Validate any
//! real design against measured device curves and manufacturer
//! datasheets.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod regulator;

pub use error::{ErrorCategory, ZenerError};
pub use regulator::{size_series_resistor, OperatingPoint, ZenerRegulator};
