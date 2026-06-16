//! # valenx-rcbeam — reinforced-concrete beam flexure
//!
//! Closed-form ultimate-strength flexure for a singly-reinforced
//! rectangular concrete section, using the Whitney equivalent
//! rectangular stress block.
//!
//! ## What
//!
//! Given a rectangular beam section (width `b`, effective depth `d`),
//! its materials (concrete strength `fc'`, steel yield `fy`) and the
//! tension-steel area `As`, this crate computes the equivalent
//! stress-block depth `a`, the internal lever arm `d - a/2`, the
//! nominal moment capacity `Mn`, the reinforcement ratio `rho` and
//! balanced ratio `rho_b`, whether the section is under-reinforced
//! (ductile), and the design strength `phi*Mn`.
//!
//! It also runs the **design direction**:
//! [`BeamSection::required_steel_area`] inverts the capacity equation to
//! the tension-steel area `As` needed for a target nominal moment, and
//! [`BeamSection::for_nominal_moment`] returns the sized section directly.
//!
//! The single public type is [`BeamSection`]; the [`error`] module
//! carries the [`RcBeamError`] taxonomy.
//!
//! ## Model
//!
//! Horizontal-force equilibrium on the section (concrete compression
//! block `C = 0.85*fc'*b*a` balancing steel tension `T = As*fy`)
//! yields the four governing equations:
//!
//! ```text
//! a      = As * fy / (0.85 * fc' * b)
//! Mn     = As * fy * (d - a/2)
//! rho    = As / (b * d)
//! phi_Mn = phi * Mn
//! ```
//!
//! All quantities are evaluated in a single *consistent* unit system
//! (the docs and tests use SI: mm, MPa, mm^2 -> N·mm). The model
//! assumes the tension steel yields at nominal capacity — the
//! under-reinforced, tension-controlled regime that code-conforming
//! flexural members are designed to. A grossly over-reinforced section
//! (stress block `a` reaching the effective depth `d`) is rejected as a
//! degenerate input rather than silently returning a negative lever
//! arm.
//!
//! ## Honest scope
//!
//! This is a **research / educational grade** implementation of the
//! textbook closed-form flexure equations (the ACI-318-style
//! equivalent-stress-block method). It is **NOT** a clinical, medical,
//! or production structural-engineering tool. It models only
//! singly-reinforced rectangular sections in pure flexure: there is no
//! compression steel, no T-/L-/non-rectangular geometry, no shear,
//! torsion, deflection, crack-width, development-length, fatigue,
//! seismic detailing, load combinations, or strain-compatibility
//! check; `phi` is supplied by the caller rather than derived from the
//! net tensile strain. Do **not** use it for real construction or any
//! safety-critical decision — every design must be verified against the
//! governing building code and signed off by a licensed professional
//! engineer.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod beam;
pub mod error;

pub use beam::{BeamSection, PHI_TENSION_CONTROLLED, STRESS_BLOCK_INTENSITY};
pub use error::{ErrorCategory, RcBeamError};
