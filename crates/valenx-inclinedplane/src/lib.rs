//! # valenx-inclinedplane
//!
//! Closed-form statics of the inclined plane (ramp) — the simplest of
//! the six classical simple machines.
//!
//! ## What
//!
//! Given a block of weight `W` resting on a ramp inclined at angle
//! `theta` from the horizontal, with a Coulomb friction coefficient
//! `mu` between block and ramp, this crate computes:
//!
//! - the normal reaction `N = W cos(theta)`;
//! - the down-slope gravity component `W sin(theta)`;
//! - the available Coulomb friction force `mu * N`;
//! - the ideal (frictionless) mechanical advantage `MA = 1 / sin(theta)`;
//! - the slope-parallel effort to raise the load,
//!   `F_up = W (sin(theta) + mu cos(theta))`;
//! - the slope-parallel effort to lower / hold the load,
//!   `F_down = W (sin(theta) - mu cos(theta))`;
//! - whether the ramp is *self-locking* — the block will not slide
//!   under gravity alone — which happens exactly when the friction
//!   angle `phi = atan(mu)` is at least the slope angle `theta`.
//!
//! ## Model
//!
//! A single rigid block in static equilibrium on a planar ramp, in a
//! uniform gravity field `g`. The effort is assumed to act parallel to
//! the inclined surface (the classic textbook configuration that gives
//! the clean `1 / sin(theta)` ideal mechanical advantage). Friction
//! follows the Amontons-Coulomb dry-friction law with a single
//! coefficient used for both the static threshold and the kinetic
//! resistance; the normal force is taken as `W cos(theta)` (no applied
//! load component pressing into or pulling off the surface beyond the
//! block's own weight). Angles are in radians; all forces share
//! whatever consistent unit the caller's weight `W` is expressed in
//! (newtons, pounds-force, etc.).
//!
//! The ideal mechanical advantage is the *velocity ratio*: sliding the
//! block a distance `L` along the ramp raises it by `L sin(theta)`, so
//! a frictionless effort `W sin(theta)` over distance `L` does the same
//! work as lifting `W` straight up by `L sin(theta)`, giving
//! `MA = W / F = 1 / sin(theta)`. A steeper ramp (larger `theta`) has a
//! smaller `MA` — less force multiplication but a shorter push.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are the textbook closed-form
//! rigid-body equilibrium relations for a point/block on an ideal
//! planar ramp: no deformation, no rolling, no toppling about an edge,
//! no distinction between static and kinetic friction coefficients, no
//! dynamics (the block is assumed at rest or on the verge of motion),
//! and no second body, pulley, or wedge interaction. It is NOT a
//! clinical/medical/production engineering tool and must not be used to
//! certify a real ramp, lift, conveyor, or restraint. Always validate
//! against measured data and a qualified engineer before relying on any
//! of these numbers in the physical world.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod geometry;
pub mod statics;

pub use error::{ErrorCategory, InclinedPlaneError};
pub use geometry::{IdealRamp, GRAVITY_STANDARD};
pub use statics::{Ramp, RampForces};
