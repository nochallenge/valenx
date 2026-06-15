//! # valenx-dcmotor
//!
//! Steady-state model of a brushed permanent-magnet DC motor.
//!
//! ## What
//!
//! Closed-form armature relations for a brushed PM DC motor, evaluated
//! at a single steady operating point. The model provides the back-EMF
//! `E = Ke * omega`, the electromagnetic torque `T = Kt * I`, the
//! armature voltage balance `V = I*R + Ke*omega`, the stall current
//! `I = V/R` and stall torque `T = Kt*V/R` (`omega = 0`), the ideal
//! no-load speed `omega = V/Ke` (`T = 0`), the linear torque-speed line
//! joining the no-load and stall points, and the electrical / mechanical
//! power split with efficiency `T*omega/(V*I)`.
//!
//! The public surface is [`DcMotor`] (the validated parameter set and
//! all evaluators), [`OperatingPoint`] (a solved point with its power
//! split), and [`DcMotorError`] for the validated constructors.
//!
//! ## Model
//!
//! A single-loop armature with constant magnetic flux. In steady state
//! the inductive term `L di/dt` vanishes, leaving the algebraic voltage
//! balance `V = I*R + E` with `E = Ke*omega`. Torque is proportional to
//! current, `T = Kt*I`. Eliminating `I` between the two yields the
//! straight torque-speed characteristic
//!
//! `omega = V/Ke - (R / (Kt*Ke)) * T`
//!
//! whose intercepts are the no-load speed `V/Ke` and the stall torque
//! `Kt*V/R`. In coherent SI units `Kt` and `Ke` are numerically equal,
//! so the back-EMF power `E*I` equals the mechanical power `T*omega` and
//! the only loss is armature copper loss `I^2 R`. These are the
//! relations in any introductory electric-machines text (e.g. Fitzgerald
//! and Kingsley, *Electric Machinery*).
//!
//! ## Honest scope
//!
//! Research/educational grade; textbook closed-form/numerical models;
//! NOT a clinical/medical/production engineering tool. The model is the
//! idealised steady-state line only. It deliberately ignores armature
//! inductance and all electrical transients, brush/contact voltage drop,
//! magnetic saturation and armature reaction, temperature-dependent
//! winding resistance, and every mechanical loss (friction, windage,
//! iron/eddy losses) — so the no-load speed is the ideal `V/Ke` rather
//! than the slightly lower real value, and efficiency reflects copper
//! loss only. Do not use it to size, qualify, or certify a real motor,
//! drive, or product; cross-check against measured data and a vendor
//! datasheet for any engineering decision.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod motor;

pub use error::{DcMotorError, ErrorCategory};
pub use motor::{DcMotor, OperatingPoint};
