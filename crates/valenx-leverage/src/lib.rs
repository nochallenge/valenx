//! # valenx-leverage
//!
//! Ideal rigid-lever mechanics: lever classification, mechanical
//! advantage, and the static moment-balance law.
//!
//! ## What
//!
//! A small, dependency-light core for reasoning about simple machines of
//! the lever family. It models a lever by its two arm lengths and offers
//! the textbook relations between effort, load, arm geometry, and the
//! resulting mechanical advantage:
//!
//! - [`Lever`] — a beam defined by `effort_arm` and `load_arm`.
//! - [`LeverClass`] — first / second / third class.
//! - [`Lever::mechanical_advantage`] — the ratio `effort_arm / load_arm`.
//! - [`Lever::balance_load`] / [`Lever::balance_effort`] — solve the
//!   balance law for the unknown force.
//! - [`Lever::net_moment`] / [`Lever::is_balanced`] — equilibrium check
//!   for an arbitrary effort/load pair.
//!
//! ## Model
//!
//! The lever is treated as an *ideal rigid lever*: a massless,
//! perfectly rigid beam pivoting on a frictionless fulcrum, with point
//! forces applied perpendicular to their moment arms. Under those
//! assumptions the system is in static equilibrium when the moments
//! about the fulcrum cancel,
//!
//! `effort * effort_arm = load * load_arm`
//!
//! and the dimensionless ideal mechanical advantage is the arm ratio,
//!
//! `MA = effort_arm / load_arm = load / effort`   (at balance).
//!
//! `MA > 1` multiplies force (second-class levers always; first-class
//! when the fulcrum sits nearer the load); `MA < 1` divides it, trading
//! force for distance and speed (third-class levers always).
//!
//! ## Honest scope
//!
//! Research/educational grade. These are closed-form, idealized
//! statics formulas straight out of an introductory mechanics text:
//! no beam mass or flexure, no friction or bearing losses, no dynamics,
//! no material strength, and forces assumed perpendicular to the arms.
//! Real levers deviate (the beam has weight and bends, the fulcrum has
//! friction, forces act at an angle, the actual mechanical advantage is
//! lower than the ideal). This crate is **NOT** a clinical/medical or
//! production engineering tool — do not use it to size real structural
//! members, machinery, prosthetics, or any load-bearing component.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod lever;

pub use error::{ErrorCategory, LeverError};
pub use lever::{Balance, Lever, LeverClass};
