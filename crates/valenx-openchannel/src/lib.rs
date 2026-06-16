//! # valenx-openchannel
//!
//! Open-channel (free-surface) hydraulics for the Valenx workspace.
//!
//! ## What
//!
//! Textbook steady-uniform-flow models for prismatic open channels:
//!
//! - Cross-section geometry ([`Channel`]) for rectangular and
//!   trapezoidal shapes — flow area `A`, wetted perimeter `P`, top
//!   width `T`, hydraulic radius `R = A/P` and hydraulic depth
//!   `D = A/T`.
//! - Manning's equation ([`manning`]) — mean velocity
//!   `v = (1/n) R^(2/3) S^(1/2)` and discharge `Q = v A`, plus a
//!   normal-depth solver that inverts `Q(y)`.
//! - The Froude number ([`froude`]) `Fr = v / sqrt(g D)`, flow-regime
//!   classification, specific energy `E = y + v^2 / (2 g)`, the
//!   critical-depth solution (`Fr = 1`), and the hydraulic-jump sequent
//!   (conjugate) depth from the Bélanger equation.
//!
//! ## Model
//!
//! Everything assumes steady, gradually-varied-free, **uniform**
//! (normal) flow in a prismatic channel with a mild bed slope and a
//! single roughness around the wetted perimeter. The governing relations
//! are the classical empirical / energy formulas:
//!
//! - Manning (SI form): `Q = (1/n) A R^(2/3) S^(1/2)`, `R = A / P`.
//! - Froude: `Fr = v / sqrt(g D)`, with `D = A / T` the hydraulic depth.
//! - Specific energy: `E = y + v^2 / (2 g)`, minimised at critical depth.
//! - Critical flow: `Q^2 T = g A^3` (i.e. `Fr = 1`).
//! - Hydraulic jump (rectangular): the Bélanger sequent-depth relation
//!   `y2 = (y1/2)(sqrt(1 + 8 Fr1^2) - 1)`, conjugate to `y1`.
//!
//! Inputs and outputs are SI throughout (metres, m³/s, m/s); the
//! US-customary `1.49` Manning factor is intentionally not used. The
//! iterative solvers (normal depth, critical depth) use bisection on
//! monotone residuals, so the roots are unique and bracketed.
//!
//! ## Honest scope
//!
//! Research/educational grade: textbook closed-form / numerical models,
//! NOT a clinical/medical/production engineering tool. The Manning
//! relation is an empirical correlation valid for fully-turbulent rough
//! flow; real channels involve compound sections, varying roughness,
//! sediment transport, backwater and unsteady effects that this crate
//! does not model. Do not use these results for real hydraulic design,
//! flood modelling, or any safety-critical decision — they are intended
//! for learning and first-order estimates only.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod froude;
pub mod geometry;
pub mod manning;

pub use error::OpenChannelError;
pub use froude::{
    classify_regime, critical_depth, froude_for_discharge, froude_number, sequent_depth,
    specific_energy_for_discharge_m, specific_energy_m, FlowRegime, GRAVITY_M_S2,
};
pub use geometry::Channel;
pub use manning::{discharge, normal_depth, velocity, velocity_from_radius};
