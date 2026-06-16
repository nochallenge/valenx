//! # valenx-inductionmotor
//!
//! Steady-state rotating-field kinematics of a three-phase induction
//! (asynchronous) motor.
//!
//! ## What
//!
//! Given a supply frequency `f` (Hz) and a pole count, this crate
//! computes the four canonical kinematic quantities of a poly-phase
//! induction machine and bundles them into a single validated
//! [`InductionMotor`] operating point:
//!
//! - synchronous (field) speed `Ns` in rev/min,
//! - fractional slip `s` (dimensionless),
//! - rotor electrical frequency `f_r` (the slip frequency) in Hz, and
//! - rotor mechanical speed `Nr` in rev/min.
//!
//! Each relation is also exposed as a free function
//! ([`sync_speed_rpm`], [`slip`], [`rotor_frequency_hz`],
//! [`rotor_speed_rpm`]) for ad-hoc use without building an aggregate.
//!
//! From a caller-supplied **air-gap power** it also gives the slip-derived
//! **power split** — the rotor copper loss [`rotor_copper_loss_w`]
//! (`s * Pag`) and the developed mechanical power
//! [`developed_mechanical_power_w`] (`(1 - s) * Pag`) — the energy
//! partition that gives slip its physical meaning.
//!
//! ## Model
//!
//! All four quantities follow directly from the speed of the rotating
//! stator field. With `f` the line frequency in hertz and `poles` the
//! (even) number of magnetic poles:
//!
//! ```text
//! Ns  = 120 * f / poles      (synchronous speed, rev/min)
//! s   = (Ns - Nr) / Ns       (fractional slip, dimensionless)
//! f_r = s * f                (rotor / slip frequency, Hz)
//! Nr  = Ns * (1 - s)         (rotor mechanical speed, rev/min)
//! P_cu   = s * Pag           (rotor copper loss, W)
//! P_mech = (1 - s) * Pag     (developed mechanical power, W)
//! ```
//!
//! The `120` is `60 s/min` (hertz to rev/min) times the two poles per
//! pole pair. Slip is `0` exactly at synchronous speed (`Nr = Ns`) and
//! `1` at standstill (`Nr = 0`); the validated [`InductionMotor::new`]
//! constructor restricts it to the motoring range `[0, 1]`, while the
//! free [`slip`] function returns the raw arithmetic ratio so that
//! generating (`s < 0`) and plugging (`s > 1`) regimes can still be
//! inspected.
//!
//! ## Honest scope
//!
//! Research/educational grade. This crate implements only the
//! **textbook closed-form** rotating-field speed and slip relations of
//! an idealised machine. It is **NOT** a clinical/medical/production
//! engineering tool and must not be used to size, select, certify, or
//! protect real motors or drives. From a **caller-supplied** air-gap
//! power it does give the slip-derived power split (rotor copper loss
//! `s * Pag`, developed mechanical power `(1 - s) * Pag`), which follows
//! from slip and energy conservation alone. It still deliberately omits
//! the per-phase equivalent circuit — so it cannot compute the air-gap
//! power, the torque-slip, or the current-slip curves itself — along with
//! stator copper / iron / friction-windage losses, **overall** efficiency
//! and power factor, magnetic saturation and skin effect, starting and
//! locked-rotor behaviour, and any thermal, mechanical, or
//! variable-frequency-drive dynamics. Results are exact for the
//! algebraic model above and nothing more.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod motor;

pub use error::{ErrorCategory, InductionMotorError};
pub use motor::{
    developed_mechanical_power_w, rotor_copper_loss_w, rotor_frequency_hz, rotor_speed_rpm, slip,
    sync_speed_rpm, InductionMotor, SYNC_SPEED_CONSTANT,
};
