//! # valenx-statics — 2D rigid-body statics
//!
//! Force and moment equilibrium for coplanar systems, and the support
//! reactions of a simply-supported beam. The whole crate is the
//! textbook "statics of rigid bodies" chapter rendered as a small,
//! validated Rust library: planar forces, the scalar moment
//! `M = r x F`, the three equilibrium equations, and the closed-form
//! lever-rule solution for a pin + roller beam.
//!
//! ## What
//!
//! - **Forces & moments** ([`force`]) — [`AppliedForce`] (a force at a
//!   point) and the scalar out-of-plane moment
//!   [`moment_z`] / [`AppliedForce::moment_about`], the `z`-component of
//!   `r x F`.
//! - **Equilibrium** ([`equilibrium`]) — a [`ForceSystem`] of coplanar
//!   forces with resultant force `(sum Fx, sum Fy)`, resultant moment
//!   `sum M` about any pivot, and an [`ForceSystem::is_in_equilibrium`]
//!   test of `sum Fx = sum Fy = sum M = 0`.
//! - **Beam reactions** ([`beam`]) — a [`SimpleBeam`] (pin + roller)
//!   carrying [`PointLoad`]s, solved in closed form to its
//!   [`Reactions`], plus the single equivalent [`beam::LoadResultant`]
//!   of its vertical loads (magnitude `sum P` at the load centroid).
//!
//! ```
//! use valenx_statics::SimpleBeam;
//!
//! // A 10 m simply-supported beam with a 100 N load at mid-span.
//! let mut beam = SimpleBeam::new(0.0, 10.0).expect("valid span");
//! beam.add_vertical(5.0, 100.0).expect("finite load");
//! let r = beam.reactions();
//!
//! // Central load -> equal reactions of P/2 = 50 N each.
//! assert!((r.pin_vertical - 50.0).abs() < 1e-9);
//! assert!((r.roller_vertical - 50.0).abs() < 1e-9);
//! ```
//!
//! ## Model
//!
//! A **rigid body in the plane** has three degrees of freedom, so it is
//! governed by exactly three scalar equilibrium equations:
//!
//! ```text
//! sum Fx = 0,    sum Fy = 0,    sum M = 0.
//! ```
//!
//! The moment of a force about a point is the `z`-component of the cross
//! product of the position vector (pivot to point of application) with
//! the force:
//!
//! ```text
//! M = r_x * F_y - r_y * F_x      (counter-clockwise positive).
//! ```
//!
//! A **simply-supported beam** carries a pin (2 vertical + horizontal
//! reaction unknowns) and a roller (1 vertical reaction unknown) — three
//! unknowns for three equations, so it is **statically determinate**.
//! Taking moments about the pin isolates the roller reaction; vertical
//! and horizontal balance recover the pin's two components. The result
//! is the classic **lever rule**: a load at distance `d` from one
//! support throws the fraction `1 - d/L` of its weight onto that
//! support, so a mid-span load splits `P/2 : P/2` and a load near a
//! support loads that support more.
//!
//! Errors are reported through [`StaticsError`], which carries stable
//! [`code`](StaticsError::code) and [`category`](StaticsError::category)
//! accessors for telemetry. Geometry is validated at construction: a
//! beam whose supports coincide is rejected as
//! [`StaticsError::Degenerate`] rather than dividing by a zero span.
//!
//! ## Honest scope
//!
//! Research / educational grade. Every formula here is the genuine
//! textbook article — the scalar moment, the three planar equilibrium
//! equations, and the determinate pin + roller reaction solution are
//! exact closed form, validated against analytic known values (central
//! load gives `P/2`; the solved beam satisfies all three equilibrium
//! equations; the lever rule is symmetric). It is deliberately a
//! focused v1 and is **not** a clinical/medical or production
//! engineering tool:
//!
//! - **Rigid bodies only.** No material deformation, stress, strain,
//!   bending moment / shear-force diagrams, or beam deflection — those
//!   belong to mechanics-of-materials and finite-element analysis
//!   (`valenx-fem`), not to rigid-body statics.
//! - **Planar (2-D) and statically determinate.** Three reactions for
//!   three equations. Statically indeterminate structures (extra
//!   supports / redundant members) and full 3-D (6-equation) statics
//!   are out of scope; the [`StaticsError::Indeterminate`] variant
//!   names that boundary.
//! - **Idealised supports and point loads.** Frictionless pin and
//!   roller, point forces at exact positions. Distributed loads must be
//!   reduced to equivalent point resultants by the caller; there is no
//!   truss method-of-joints / sections solver, no friction, and no
//!   self-weight unless added explicitly as a load.
//!
//! None of those omissions makes a result wrong — the reactions, the
//! resultant force and the resultant moment are real engineering
//! numbers for the idealised rigid body, each a documented,
//! well-understood building block on the way to a fuller
//! structural-mechanics suite.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod beam;
pub mod equilibrium;
pub mod error;
pub mod force;

pub use beam::{LoadResultant, PointLoad, Reactions, SimpleBeam};
pub use equilibrium::ForceSystem;
pub use error::{ErrorCategory, Result, StaticsError};
pub use force::{moment_z, AppliedForce};

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    /// End-to-end: solve a beam, then independently confirm the solved
    /// load + reaction system is in equilibrium under all three planar
    /// equations.
    #[test]
    fn beam_solution_is_an_equilibrium_force_system() {
        let mut beam = SimpleBeam::new(0.0, 6.0).expect("valid span");
        beam.add_vertical(2.0, 30.0).expect("finite");
        beam.add_vertical(4.0, 50.0).expect("finite");

        let r = beam.reactions();
        // Reactions sum to the applied load (sum Fy = 0).
        assert!((r.total_vertical() - 80.0).abs() < EPS);

        // The assembled system passes the independent equilibrium check.
        let sys = beam.equilibrium_system();
        assert!(sys.is_in_equilibrium(1e-6));
    }

    /// The lever rule, end-to-end through the public surface: a load at
    /// 1/4 span puts 3/4 of its weight on the near support.
    #[test]
    fn quarter_point_load_lever_rule() {
        let mut beam = SimpleBeam::new(0.0, 4.0).expect("valid span");
        beam.add_vertical(1.0, 80.0).expect("finite"); // 1/4 span from pin
        let r = beam.reactions();
        // Near (pin) support carries 3/4; far (roller) carries 1/4.
        assert!(
            (r.pin_vertical - 60.0).abs() < EPS,
            "pin {}",
            r.pin_vertical
        );
        assert!(
            (r.roller_vertical - 20.0).abs() < EPS,
            "roller {}",
            r.roller_vertical
        );
    }

    /// A standalone moment computation through the re-exported helper.
    #[test]
    fn reexported_moment_helper_matches_definition() {
        use nalgebra::Vector2;
        let m = moment_z(Vector2::new(2.0, 0.0), Vector2::new(0.0, 5.0));
        assert!((m - 10.0).abs() < EPS, "got {m}");
    }
}
