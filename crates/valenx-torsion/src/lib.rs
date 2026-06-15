//! # valenx-torsion
//!
//! Closed-form linear-elastic **torsion of circular shafts**: the polar
//! second moment of area, the shear-stress distribution, the angle of
//! twist, and the power a rotating shaft transmits.
//!
//! ## What
//!
//! Describe a round shaft as a [`Shaft`] (solid bar or hollow tube),
//! apply a torque, and read back the engineering answers:
//!
//! - the polar second moment of area [`Shaft::polar_moment`] `J`,
//! - the shear stress [`shear_stress_at`] anywhere across the section and
//!   its maximum [`max_shear_stress`] at the surface,
//! - the [`angle_of_twist`] over a length of shaft, and
//! - the [`power`] transmitted at a given angular speed.
//!
//! [`TorsionCase`] bundles a shaft with a load case and
//! [`TorsionCase::analyse`] returns every quantity at once in a
//! [`TorsionResult`].
//!
//! ```
//! use valenx_torsion::{Shaft, TorsionCase};
//!
//! // A 30 mm solid steel shaft, 1.5 m long, carrying 250 N·mm at 12 rad/s.
//! let shaft = Shaft::solid(30.0).expect("positive diameter");
//! let case = TorsionCase::new(shaft, 250.0, 1_500.0, 79_300.0, 12.0)
//!     .expect("valid load case");
//! let r = case.analyse().expect("closed-form evaluation");
//! println!(
//!     "J = {:.0}, tau_max = {:.3}, theta = {:.4} rad, P = {:.0}",
//!     r.polar_moment, r.max_shear_stress, r.angle_of_twist, r.power,
//! );
//! ```
//!
//! ## Model
//!
//! These are the textbook results of **St. Venant torsion of a prismatic
//! circular bar** in the linear-elastic regime. With torque `T`, length
//! `L`, shear modulus `G`, radius `r`, diameter `d` (outer `D` and bore
//! `d` for a tube), and angular speed `omega`:
//!
//! ```text
//! J        = pi * d^4 / 32                (solid)
//! J        = pi * (D^4 - d^4) / 32        (hollow)
//! tau(r)   = T * r / J                    (shear stress, linear in r)
//! tau_max  = T * (d / 2) / J             (at the outer surface)
//! theta    = T * L / (G * J)             (angle of twist)
//! P        = T * omega                    (transmitted power)
//! ```
//!
//! For a circular section the polar second moment of area equals the
//! torsion constant, so the same `J` appears in both the stress and the
//! twist. Units are the caller's responsibility but must be consistent;
//! SI (`T` in N·m, lengths in m, `G` in Pa, `omega` in rad/s) yields
//! stress in Pa, twist in rad and power in W.
//!
//! ## Honest scope
//!
//! Research/educational grade. Every formula here is the exact
//! closed-form / numerical model from a mechanics-of-materials textbook
//! and is validated against analytic ground truth in the test suite, but
//! the crate deliberately models only the idealised case:
//!
//! - **Linear-elastic, prismatic, circular** shafts only — no plasticity,
//!   creep, or material non-linearity, and no non-circular sections
//!   (which warp and need a different torsion constant).
//! - **Uniform pure torsion** — no stress-concentration factors at
//!   shoulders, keyways, holes or fillets, no combined axial/bending
//!   loads, and no fatigue or failure-criterion checks.
//! - **Static / quasi-static** — the power relation is the kinematic
//!   `P = T omega`; there is no rotordynamics, vibration, or shaft-whirl
//!   modelling.
//!
//! It is **not** a clinical/medical tool and **not** a production
//! engineering or code-compliance tool. Do not use it as the sole basis
//! for a safety-critical design decision; verify against an applicable
//! standard and qualified review.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod analysis;
pub mod error;
pub mod response;
pub mod shaft;

pub use analysis::{TorsionCase, TorsionResult};
pub use error::{ErrorCategory, TorsionError};
pub use response::{angle_of_twist, max_shear_stress, power, shear_stress_at, torsional_rigidity};
pub use shaft::Shaft;
