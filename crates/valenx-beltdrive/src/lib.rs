//! # valenx-beltdrive
//!
//! Closed-form analysis of flat- and V-belt drives: speed ratio, belt
//! linear speed, capstan friction limit, transmitted power, and
//! centrifugal tension.
//!
//! ## What
//!
//! Given pulley diameters, rotational speeds, a belt/pulley coefficient
//! of friction, an angle of wrap, and tensions, this crate computes:
//!
//! - the transmission (velocity) ratio `i = D_driven / D_driver`
//!   ([`geometry::speed_ratio`]),
//! - the belt linear speed `v = pi * D * N` ([`geometry::belt_speed`]),
//! - open-belt wrap angles and exact belt length
//!   ([`geometry::wrap_angles_open`], [`geometry::belt_length_open`]),
//!   and the crossed-belt (counter-rotating) geometry
//!   ([`geometry::wrap_angle_crossed`], [`geometry::belt_length_crossed`]),
//! - the Euler / capstan tension ratio `T1 / T2 = exp(mu * theta)`
//!   ([`friction::tension_ratio`]), including the V-belt
//!   `mu_eff = mu / sin(beta)` amplification
//!   ([`friction::v_belt_tension_ratio`]),
//! - transmitted power `P = (T1 - T2) * v`
//!   ([`power::transmitted_power`]),
//! - centrifugal tension `Tc = m * v^2`
//!   ([`power::centrifugal_tension`]),
//! - and the slipping-limited power capacity and the speed that
//!   maximises it ([`power::max_power`], [`power::speed_for_max_power`]).
//!
//! A whole configuration can be bundled in a serialisable
//! [`spec::DriveSpec`] and analysed in one call via
//! [`spec::DriveSpec::analyze`], which returns a [`spec::DriveAnalysis`]
//! whose fields come from the free functions above.
//!
//! ## Model
//!
//! Every formula is the standard machine-design idealisation: rigid
//! pulleys, an inextensible belt of uniform linear mass density, and no
//! slip (the belt is on the verge of slipping when the capstan limit is
//! applied). Surface speeds at the two pulleys are equal, so the speed
//! ratio follows directly from the diameters. The tension ratio is the
//! Euler capstan equation; for a V-belt the groove half-angle `beta`
//! wedges the belt and raises the effective friction to
//! `mu / sin(beta)`. Power is the net driving force `T1 - T2` times the
//! belt speed, and the centrifugal contribution `Tc = m * v^2` is
//! carried equally on both sides, reducing the tension available to
//! drive the load. All quantities are SI unless a function's
//! documentation states otherwise.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form and
//! numerical models for learning and first-cut estimation — they are
//! **not** a clinical/medical/production engineering tool. The crate
//! deliberately omits everything a real drive-selection workflow needs:
//! belt fatigue and creep, elastic slip, service factors, manufacturer
//! power-rating and length tables, groove-and-belt fit tolerances,
//! thermal effects, dynamic / resonance behaviour, and arc-of-contact
//! correction factors beyond the bare capstan relation. Do not size,
//! certify, or operate physical machinery from its output. Validate any
//! result against the relevant standard (e.g. an ISO / manufacturer
//! catalogue) before use.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod friction;
pub mod geometry;
pub mod power;
pub mod spec;

pub use error::{BeltError, ErrorCategory};
pub use friction::{slack_tension, tension_ratio, v_belt_effective_mu, v_belt_tension_ratio};
pub use geometry::{
    belt_length_open, belt_speed, driven_speed, rpm_to_rev_per_sec, speed_ratio, wrap_angles_open,
};
pub use power::{centrifugal_tension, max_power, speed_for_max_power, transmitted_power};
pub use spec::{DriveAnalysis, DriveSpec};
