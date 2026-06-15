//! # valenx-camdynamics
//!
//! Closed-form cam-follower motion laws for the rise (and, by symmetry,
//! the return) segment of a disk-cam profile.
//!
//! ## What
//!
//! Given a follower `lift` (the total rise displacement) and a cam
//! rotation interval `beta` (the rise angle, in radians) over which that
//! rise occurs, this crate evaluates the four kinematic quantities of the
//! follower as a function of the cam angle `theta` measured from the start
//! of the rise:
//!
//! - displacement `s(theta)`,
//! - velocity `v(theta) = ds/dtheta`,
//! - acceleration `a(theta) = d2s/dtheta2`,
//! - jerk `j(theta) = d3s/dtheta3`.
//!
//! Two textbook motion laws are provided:
//!
//! - [`MotionLaw::SimpleHarmonic`] — simple-harmonic motion (SHM), and
//! - [`MotionLaw::Cycloidal`] — cycloidal motion.
//!
//! All derivatives are taken with respect to the cam angle `theta`, so
//! velocity carries units of `length / radian`, acceleration
//! `length / radian^2`, and jerk `length / radian^3`. To convert to a
//! time basis at a constant cam angular speed `omega` (rad/s), multiply
//! by the corresponding power of `omega`: `v_t = omega * v`,
//! `a_t = omega^2 * a`, `j_t = omega^3 * j`.
//!
//! ## Model
//!
//! Let `x = theta / beta` be the normalised position in the interval
//! `[0, 1]`. The two laws use the standard normalised forms (see Norton,
//! *Design of Machinery*, or Shigley, *Mechanical Engineering Design*):
//!
//! Simple harmonic (SHM):
//!
//! - `s = (L/2) (1 - cos(pi x))`
//! - `v = (pi L / (2 beta)) sin(pi x)`
//! - `a = (pi^2 L / (2 beta^2)) cos(pi x)`
//! - `j = -(pi^3 L / (2 beta^3)) sin(pi x)`
//!
//! Cycloidal:
//!
//! - `s = L (x - sin(2 pi x) / (2 pi))`
//! - `v = (L / beta) (1 - cos(2 pi x))`
//! - `a = (2 pi L / beta^2) sin(2 pi x)`
//! - `j = (4 pi^2 L / beta^3) cos(2 pi x)`
//!
//! where `L` is the lift and `beta` the rise angle.
//!
//! The two laws differ in smoothness at the segment ends. SHM leaves a
//! finite, non-zero acceleration at both ends (a step in acceleration
//! when joined to a dwell), whereas the cycloidal law brings both the
//! velocity and the acceleration to zero at each end, leaving only a
//! finite jerk discontinuity. The cycloidal law is therefore the
//! smoother of the two and is the classic choice for high-speed cams.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are idealised, rigid-body,
//! closed-form motion laws drawn directly from textbook kinematics. They
//! assume a perfectly rigid follower train, zero clearance, no friction,
//! and no dynamic (inertial or vibratory) response of the real mechanism.
//! Pressure angle, radius of curvature, contact stress, follower jump,
//! and manufacturing tolerances are all out of scope. This crate is NOT a
//! clinical, medical, or production engineering tool and must not be used
//! to certify a physical cam design; validate any real mechanism with
//! measurement and appropriate engineering analysis.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod motion;

pub use error::{CamError, ErrorCategory};
pub use motion::{FollowerState, MotionLaw, RiseProfile};
