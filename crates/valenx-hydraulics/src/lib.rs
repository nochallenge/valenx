//! # valenx-hydraulics
//!
//! Closed-form analysis of steady-state hydraulic-power components:
//! linear actuators (cylinders), control valves, and the power they
//! transmit.
//!
//! ## What
//!
//! Given pressures, areas, velocities, and valve flow coefficients,
//! this crate computes:
//!
//! - cylinder thrust on the bore (cap) side and the rod (annulus)
//!   side, with the extend vs. retract asymmetry that the rod
//!   area introduces;
//! - volumetric flow from the continuity relation `Q = A v`, and the
//!   piston speed it implies for a commanded flow;
//! - control-valve flow from the flow coefficient relation
//!   `Q = Cv sqrt(dP / SG)`;
//! - hydraulic power `Power = p Q`.
//!
//! Every quantity is SI-coherent inside the formulas; the public API
//! documents the unit attached to each argument and result.
//!
//! ## Model
//!
//! The relations are the standard incompressible, steady-flow
//! textbook forms:
//!
//! - Pascal force balance on a piston: `F = p A`. The bore side sees
//!   the full piston area `A_bore = pi D^2 / 4`; the rod side sees the
//!   annulus `A_rod = A_bore - pi d^2 / 4`, where `d` is the rod
//!   diameter. Because `A_rod < A_bore`, the extend stroke develops
//!   more force than the retract stroke at equal pressure.
//! - Continuity for an incompressible fluid through a constant area:
//!   `Q = A v`, with `v = Q / A` the implied piston speed.
//! - Orifice / control-valve sizing in the flow-coefficient form
//!   `Q = Cv sqrt(dP / SG)`. Here `SG` is the fluid specific gravity
//!   (dimensionless, water = 1). Flow scales with the square root of
//!   the pressure drop, so quadrupling `dP` doubles `Q`.
//! - Transmitted hydraulic power `Power = p Q` (pressure times
//!   volumetric flow).
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form and
//! numerical models, NOT a clinical, medical, or production
//! engineering tool. The model is steady-state and incompressible and
//! deliberately ignores: fluid compressibility and entrained air,
//! distributed line and fitting losses, valve hysteresis and dynamic
//! response, cavitation and aeration, viscous heating and temperature
//! dependence, seal friction, and structural limits of the hardware.
//! `Cv` is treated as a pure constant rather than an
//! opening-dependent characteristic. Do not size, certify, or operate
//! real hydraulic equipment from these numbers; use validated
//! engineering tooling and physical testing for that.
//!
//! # Surface
//!
//! - [`Cylinder`] — bore / rod geometry with area, force, and
//!   speed accessors.
//! - [`Stroke`] — extend vs. retract selector.
//! - [`valve::valve_flow`] — flow-coefficient valve flow.
//! - [`power::hydraulic_power`] — `Power = p Q`.
//! - [`HydraulicsError`] — validated-constructor error type.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod cylinder;
pub mod error;
pub mod power;
pub mod valve;

pub use cylinder::{Cylinder, Stroke};
pub use error::{ErrorCategory, HydraulicsError};
pub use power::hydraulic_power;
pub use valve::valve_flow;
