//! # valenx-electrochem
//!
//! Textbook closed-form electrochemistry: the Nernst equation, cell
//! potential, and Faraday's law of electrolysis, as pure `f64` algorithms
//! with validated constructors.
//!
//! ## What
//!
//! Three classical results, each with a free-function form and a small
//! validated struct. The [`nernst`] module gives the Nernst equation
//! `E = E0 - (R T / (n F)) ln Q` for the reduction potential of a
//! half-reaction at arbitrary reaction quotient `Q` and temperature `T`.
//! The [`cell`] module gives the cell potential
//! `E_cell = E_cathode - E_anode` from two reduction potentials, plus the
//! spontaneity sign that follows. The [`faraday`] module gives Faraday's
//! law of electrolysis `m = (Q M) / (n F)` for the mass moved by a charge,
//! with `Q = I t` for a constant current.
//!
//! The shared physical constants (the molar gas constant `R`, the Faraday
//! constant `F`, and the standard temperature) live in [`constants`], and
//! every fallible constructor reports through [`error::ElectrochemError`].
//!
//! ## Model
//!
//! All quantities are SI-flavoured: potentials in volts, temperature in
//! kelvin, charge in coulombs, current in amperes, time in seconds, molar
//! mass in grams per mole, and the reaction quotient `Q` dimensionless.
//! Reduction potentials are referenced to the standard hydrogen electrode.
//! The reaction quotient is written reduced-over-oxidised, in the same
//! direction as the half-reaction, so that `E == E0` exactly at `Q = 1`.
//!
//! Useful sanity checks the code is tested against: the thermal voltage
//! `R T / F` is about `0.0257 V` at 298.15 K, giving a base-10 Nernst slope
//! of about `0.0592 V` per decade of `Q` for `n = 1` at 25 C; with `n > 0`
//! the potential falls as `Q` rises; the cell potential is exactly the
//! cathode reduction potential minus the anode reduction potential; and the
//! Faraday mass scales linearly with charge and inversely with `n`.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form analytic
//! models intended for learning and back-of-envelope estimates. The
//! reaction quotient uses ideal activities only — no activity coefficients,
//! no liquid-junction potentials, no electrode kinetics or overpotential,
//! and (for electrolysis) an assumed 100 % current efficiency. The results
//! are equilibrium / open-circuit values, not voltages or yields measured
//! under load. This crate is NOT a clinical, medical, or production
//! engineering tool and is not validated for safety-critical, regulatory,
//! or industrial use.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod cell;
pub mod constants;
pub mod error;
pub mod faraday;
pub mod nernst;

pub use cell::{cell_potential, spontaneity, Cell, Spontaneity};
pub use constants::{
    FARADAY_C_PER_MOL, GAS_CONSTANT_J_PER_MOL_K, STANDARD_TEMPERATURE_K, ZERO_CELSIUS_IN_KELVIN,
};
pub use error::{ElectrochemError, ErrorCategory};
pub use faraday::{
    charge_from_current, mass_from_charge, mass_from_current, moles_from_charge, Electrolysis,
};
pub use nernst::{
    nernst_potential, nernst_slope_per_decade, thermal_voltage, thermal_voltage_standard,
    HalfReaction,
};
