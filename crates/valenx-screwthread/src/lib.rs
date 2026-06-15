//! # valenx-screwthread
//!
//! Closed-form analysis of a square-thread power / lead screw: lead,
//! lead angle, raising and lowering torque, mechanical efficiency, and
//! the self-locking criterion.
//!
//! ## What
//!
//! Given a screw's mean diameter, pitch, number of thread starts, and
//! the screw/nut friction coefficient, this crate computes:
//!
//! - the lead `l = pitch * starts` and the lead (helix) angle
//!   `lambda = atan(l / (pi * dm))`;
//! - the raising torque `T_R` (driving a load against its direction,
//!   e.g. lifting) and the lowering torque `T_L`;
//! - the thread-friction mechanical efficiency of the raising motion;
//! - whether the screw is self-locking (cannot be back-driven by the
//!   axial load alone).
//!
//! ## Model
//!
//! The square-thread power-screw equations from a standard machine
//! design text (Shigley / Budynas, *Mechanical Engineering Design*)
//! are used throughout, with `F` the axial load, `dm` the mean
//! (pitch-line) diameter, `l` the lead, and `mu` the coefficient of
//! friction:
//!
//! - raise: `T_R = (F*dm/2) * (l + pi*mu*dm) / (pi*dm - mu*l)`;
//! - lower: `T_L = (F*dm/2) * (pi*mu*dm - l) / (pi*dm + mu*l)`;
//! - efficiency: `e = (F*l) / (2*pi*T_R)`, equal to
//!   `tan(lambda) * (1 - mu*tan(lambda)) / (tan(lambda) + mu)`;
//! - self-locking when the friction angle `phi = atan(mu)` is at least
//!   the lead angle `lambda`, equivalently `mu >= tan(lambda)`,
//!   equivalently `T_L >= 0`.
//!
//! Acme / trapezoidal threads are approximated by the square-thread
//! relations; the thread half-angle correction that divides the
//! friction terms by `cos(alpha_n)` is intentionally omitted, as is
//! collar friction. The returned torques are the thread-only
//! components.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form models
//! built on the idealised Coulomb-friction power-screw derivation; they
//! ignore thread-form half-angle effects, collar/thrust-bearing
//! friction, elastic and thermal deformation, manufacturing tolerance,
//! lubrication regime, wear, and dynamic / fatigue behaviour. This is
//! NOT a clinical/medical or production engineering tool — do not use
//! it to certify or size a safety-critical jack, actuator, clamp, or
//! fastener. Validate against a qualified reference and physical test
//! before any real-world reliance.
//!
//! ## Example
//!
//! ```
//! use valenx_screwthread::ScrewThread;
//!
//! // Single-start square-thread screw: dm = 25 mm, p = 5 mm, mu = 0.08.
//! let screw = ScrewThread::new(25.0, 5.0, 1, 0.08).unwrap();
//! assert!((screw.lead() - 5.0).abs() < 1e-12);
//! assert!(screw.is_self_locking());
//!
//! let raise = screw.raise_torque(6500.0).unwrap(); // N·mm
//! let lower = screw.lower_torque(6500.0);
//! assert!(raise > lower);
//!
//! let e = screw.efficiency().unwrap();
//! assert!(e > 0.0 && e < 1.0);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod screw;

pub use error::{ErrorCategory, ScrewThreadError};
pub use screw::ScrewThread;
