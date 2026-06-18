//! # valenx-pulley
//!
//! Closed-form mechanics of rope-and-pulley machines: fixed, movable, and
//! block-and-tackle (compound) arrangements, plus the effort, velocity
//! ratio, work, and efficiency relations that connect a load to the force
//! an operator must apply.
//!
//! ## What
//!
//! - [`PulleyKind`] / [`PulleySystem`] ([`system`]) — the three textbook
//!   arrangements, each reduced to the one number that fixes its ideal
//!   behaviour: the count of rope segments supporting the load.
//! - [`system::PulleySystem::ideal_mechanical_advantage`] /
//!   [`system::PulleySystem::velocity_ratio`] — the ideal mechanical
//!   advantage `MA = n` and the equal velocity ratio `VR = n`.
//! - [`mechanics`] — [`mechanics::ideal_effort`],
//!   [`mechanics::real_effort`], their load inverses
//!   [`mechanics::ideal_load_from_effort`] /
//!   [`mechanics::load_from_effort`],
//!   [`mechanics::actual_mechanical_advantage`],
//!   [`mechanics::efficiency_from_effort`],
//!   [`mechanics::effort_distance`], [`mechanics::output_work`],
//!   [`mechanics::input_work`] and [`mechanics::work_lost`].
//! - [`PulleyError`] ([`error`]) — a validated `thiserror` taxonomy with
//!   stable [`code`](error::PulleyError::code) /
//!   [`category`](error::PulleyError::category) accessors.
//!
//! ## Model
//!
//! The ideal mechanical advantage of a rope machine is the number of rope
//! segments directly supporting the movable load:
//!
//! ```text
//! MA = n_supporting_ropes.
//! ```
//!
//! A fixed pulley has `n = 1` (it only changes the rope's direction); a
//! movable pulley has `n = 2`; a block and tackle has `n` equal to its
//! supporting-segment count. For an ideal, friction-free machine the
//! effort is `F = W / MA`, the velocity ratio equals the mechanical
//! advantage (`VR = MA`, because the inextensible rope makes the effort
//! end travel `MA` times the load's travel), and energy is conserved so
//! input work equals output work.
//!
//! A real machine wastes part of the input work to friction. A single
//! lumped efficiency `eta` in `(0, 1]` captures this: the actual
//! mechanical advantage is `AMA = MA * eta`, the required effort rises to
//! `F = W / (MA * eta)`, and the work lost over one lift is
//! `W_out (1 / eta - 1)`.
//!
//! Read backwards, a given effort raises `W = F * MA * eta` (`F * MA` for
//! the ideal machine) — with friction the same pull lifts less.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are idealized rigid-body textbook
//! closed-form models — inextensible massless rope, point load, and either
//! friction-free sheaves (ideal) or a single scalar efficiency (real).
//! They reproduce the standard physics-class pulley results but are NOT a
//! clinical, medical, or production engineering tool. Do not size real
//! rigging, hoists, lifting tackle, or any load-bearing hardware from
//! them; real systems require rated components, dynamic and shock loads,
//! safety factors, and rope-mechanics this crate deliberately omits.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod mechanics;
pub mod system;

pub use error::{ErrorCategory, PulleyError, Result};
pub use mechanics::{
    actual_mechanical_advantage, efficiency_from_effort, effort_distance, ideal_effort,
    ideal_load_from_effort, input_work, load_from_effort, output_work, real_effort, work_lost,
};
pub use system::{PulleyKind, PulleySystem};

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    /// End-to-end worked example: lift a 1200 N load with a 4-rope
    /// block-and-tackle at 90% efficiency through 1.5 m.
    ///
    /// Ground-truth hand computation:
    ///
    /// - `MA = VR = 4`.
    /// - ideal effort `= 1200 / 4 = 300 N`.
    /// - real effort `= 1200 / (4 * 0.9) = 333.333... N`.
    /// - actual MA `= 4 * 0.9 = 3.6`.
    /// - effort travel `= 4 * 1.5 = 6 m`.
    /// - output work `= 1200 * 1.5 = 1800 J`.
    /// - input work `= 1800 / 0.9 = 2000 J`.
    /// - lost work `= 2000 - 1800 = 200 J`.
    #[test]
    fn worked_block_and_tackle_example() {
        let sys = PulleySystem::block_and_tackle(4).unwrap();
        let load = 1200.0;
        let eta = 0.9;
        let s_load = 1.5;

        assert!((sys.ideal_mechanical_advantage() - 4.0).abs() < EPS);
        assert!((sys.velocity_ratio() - 4.0).abs() < EPS);

        assert!((ideal_effort(sys, load).unwrap() - 300.0).abs() < EPS);
        assert!((real_effort(sys, load, eta).unwrap() - 1200.0 / 3.6).abs() < EPS);
        assert!((actual_mechanical_advantage(sys, eta).unwrap() - 3.6).abs() < EPS);
        assert!((effort_distance(sys, s_load).unwrap() - 6.0).abs() < EPS);
        assert!((output_work(load, s_load).unwrap() - 1800.0).abs() < EPS);
        assert!((input_work(load, s_load, eta).unwrap() - 2000.0).abs() < EPS);
        assert!((work_lost(load, s_load, eta).unwrap() - 200.0).abs() < EPS);
    }

    /// The three named arrangements have mechanical advantages 1, 2, n.
    #[test]
    fn canonical_advantages() {
        assert!((PulleySystem::fixed().ideal_mechanical_advantage() - 1.0).abs() < EPS);
        assert!((PulleySystem::movable().ideal_mechanical_advantage() - 2.0).abs() < EPS);
        assert!(
            (PulleySystem::block_and_tackle(7)
                .unwrap()
                .ideal_mechanical_advantage()
                - 7.0)
                .abs()
                < EPS
        );
    }

    /// The convenience re-exports resolve and the error type is usable.
    #[test]
    fn reexports_resolve() {
        let err: PulleyError = PulleyError::degenerate("x");
        assert_eq!(err.category_enum(), ErrorCategory::Geometry);
    }
}
