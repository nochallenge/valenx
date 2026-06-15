//! # valenx-coil — inductor / solenoid coil models
//!
//! Closed-form textbook electromagnetics for a long, uniformly wound
//! solenoid. Pure scalar algorithms, no external processes.
//!
//! ## What
//!
//! - [`Solenoid`] — a validated `N` turns / area `A` / length `l` /
//!   relative-permeability `mu_r` coil. Construct with
//!   [`Solenoid::new`] or [`Solenoid::air_core`].
//! - [`Solenoid::inductance_henries`] — the self-inductance
//!   `L = mu0 * mu_r * N^2 * A / l`.
//! - [`energy_joules`] — the stored magnetic energy `E = 0.5 * L * I^2`.
//! - [`reactance_ohms`] — the inductive reactance `X_L = 2 * pi * f * L`.
//! - [`time_constant_seconds`] — the series-RL time constant
//!   `tau = L / R`.
//! - [`VACUUM_PERMEABILITY`] — the constant `mu0 = 4 * pi * 1e-7 H/m`.
//!
//! Every fallible function returns [`Result<_, CoilError>`](CoilError),
//! whose [`code`](CoilError::code) / [`category`](CoilError::category)
//! accessors give stable handles for tests and telemetry.
//!
//! ## Model
//!
//! Under Ampère's law, an ideal long solenoid carries a uniform axial
//! field inside and (approximately) none outside, giving the four
//! canonical relations:
//!
//! ```text
//! L   = mu0 * mu_r * N^2 * A / l   (henries)
//! E   = 0.5 * L * I^2              (joules)
//! X_L = 2 * pi * f * L             (ohms)
//! tau = L / R                      (seconds)
//! ```
//!
//! Inductance is quadratic in the turn count and linear in `A`, `mu_r`
//! and `1 / l`, so doubling the turns quadruples `L`.
//!
//! ```
//! use valenx_coil::Solenoid;
//!
//! // 100 turns, 1 cm^2 cross-section, 10 cm long, air core.
//! let coil = Solenoid::air_core(100.0, 1.0e-4, 0.1).unwrap();
//! let l = coil.inductance_henries();
//!
//! // Doubling the turn count gives 4x the inductance.
//! let bigger = Solenoid::air_core(200.0, 1.0e-4, 0.1).unwrap();
//! assert!((bigger.inductance_henries() / l - 4.0).abs() < 1.0e-9);
//! ```
//!
//! ## Honest scope
//!
//! Research / educational grade: textbook closed-form / numerical
//! models, NOT a clinical / medical / production engineering tool. The
//! long-solenoid formula assumes `l` is large compared with the coil
//! diameter, so it over-estimates `L` for short, fat coils (it omits the
//! Nagaoka end-correction). It also ignores the proximity / skin effect,
//! winding self-capacitance and self-resonance, core saturation and
//! hysteresis, and frequency-dependent core loss. Use the outputs as
//! first-order estimates, not certified magnetics-design values.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod solenoid;

pub use error::{CoilError, ErrorCategory};
pub use solenoid::{
    energy_joules, reactance_ohms, time_constant_seconds, Solenoid, VACUUM_PERMEABILITY,
};
