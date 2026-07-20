//! # valenx-thevenin
//!
//! Thevenin / Norton equivalent two-terminal DC sources, the source
//! transformation between them, and the maximum power transfer theorem.
//!
//! ## What
//!
//! Given a linear two-terminal DC network reduced to its Thevenin form
//! (open-circuit voltage `V_th` in series with internal resistance `R_th`)
//! or its Norton form (short-circuit current `I_n` in parallel with `R_n`),
//! this crate provides:
//!
//! Validated [`Thevenin`] and [`Norton`] source types whose constructors
//! reject non-finite and out-of-domain inputs.
//!
//! The source transformation [`Thevenin::to_norton`] /
//! [`Norton::to_thevenin`], which round-trips exactly.
//!
//! Load analysis: [`Thevenin::load_current`], [`Thevenin::load_voltage`],
//! [`power::load_power`], and [`power::efficiency`].
//!
//! The maximum power transfer results [`power::max_power`] and
//! [`power::matched_load`].
//!
//! ## Model
//!
//! The governing equations (all standard textbook DC circuit theory):
//!
//! ```text
//! Source transform : I_n = V_th / R_th ,  R_n = R_th
//!                    V_th = I_n  * R_n  ,  R_th = R_n
//!
//! Load (divider)   : i      = V_th / (R_th + R_load)
//!                    v_load  = V_th * R_load / (R_th + R_load)
//!                    P(R_load) = V_th^2 * R_load / (R_th + R_load)^2
//!
//! Max power xfer   : dP/dR_load = 0  =>  R_load = R_th
//!                    P_max = V_th^2 / (4 * R_th)
//!                    efficiency at match = 50%
//! ```
//!
//! ## Honest scope
//!
//! Research / educational grade. These are standard textbook closed-form
//! DC circuit models, validated against analytic ground truth (source
//! transformation round-trips, the maximum power transfer theorem, and a
//! numeric power sweep). This is NOT a clinical/medical or
//! production-certified engineering tool: it models ideal lumped linear
//! elements only and accounts for no component tolerances, safety factors,
//! fatigue/thermal limits, power ratings, or code compliance.
//!
//! ## Example
//!
//! ```rust
//! use valenx_thevenin::{Thevenin, power};
//!
//! // A source with open-circuit voltage 12 V and internal resistance 6 ohm.
//! let src = Thevenin::new(12.0, 6.0).unwrap();
//!
//! // Source transformation to the Norton equivalent: I_n = 12 / 6 = 2 A.
//! let norton = src.to_norton();
//! assert!((norton.i_n - 2.0).abs() < 1e-9);
//! assert!((norton.r_n - 6.0).abs() < 1e-9);
//!
//! // Maximum power transfer: matched load equals R_th, P_max = 144 / 24.
//! assert!((power::matched_load(src) - 6.0).abs() < 1e-9);
//! assert!((power::max_power(src) - 6.0).abs() < 1e-9);
//!
//! // The matched load really does deliver the most power.
//! let p_match = power::load_power(src, 6.0).unwrap();
//! let p_off   = power::load_power(src, 20.0).unwrap();
//! assert!(p_off < p_match);
//!
//! // Efficiency at the matched load is exactly 50%.
//! assert!((power::efficiency(src, 6.0).unwrap() - 0.5).abs() < 1e-9);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod equivalent;
pub mod error;
pub mod power;

pub use equivalent::{Norton, Thevenin};
pub use error::TheveninError;
