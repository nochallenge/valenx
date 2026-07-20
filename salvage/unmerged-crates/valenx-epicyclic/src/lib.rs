//! # valenx-epicyclic
//!
//! ## What
//!
//! Closed-form kinematics of a single-stage planetary (epicyclic) gear
//! train made of a central sun, a set of orbiting planets carried on an
//! arm (the carrier), and an internal ring gear (annulus). The crate
//! computes:
//!
//! the tooth-count meshing constraint `R = S + 2 P` ([`TeethSet`]);
//! the general Willis equation relating sun, ring and carrier speeds
//! ([`willis`]);
//! the three named single-input reductions (ring fixed, sun fixed,
//! carrier fixed); and
//! the tabular (superposition) method that resolves every member,
//! including the planet's own spin ([`tabular`]).
//!
//! All speeds are angular velocities in a single arbitrary unit (rad/s,
//! rev/min, ...); every relation is homogeneous, so the output is in
//! whatever unit the inputs used.
//!
//! ## Model
//!
//! With the carrier held fixed the set is an ordinary gear train, giving
//! the signed basic train ratio
//!
//! ```text
//! e = omega_ring / omega_sun  |_(carrier fixed)  =  -S / R
//! ```
//!
//! Measuring the sun and ring speeds in the carrier's rotating frame
//! gives the Willis equation
//!
//! ```text
//! (omega_ring - omega_carrier) / (omega_sun - omega_carrier) = e
//! ```
//!
//! which rearranges to the linear constraint solved throughout:
//!
//! ```text
//! omega_ring - e*omega_sun + (e - 1)*omega_carrier = 0
//! ```
//!
//! Grounding one member yields the standard reductions
//!
//! ```text
//! ring fixed:    omega_carrier / omega_sun  =  S / (S + R)
//! sun  fixed:    omega_carrier / omega_ring =  R / (S + R)
//! carrier fixed: omega_ring    / omega_sun  = -S / R
//! ```
//!
//! ## Honest scope
//!
//! This is research/educational grade: standard textbook closed-form
//! kinematics (Willis equation, superposition table) validated against
//! analytic ground truth. It is NOT a clinical/medical or
//! production-certified engineering tool. It models only rigid-body
//! speed ratios — there is no treatment of gear module, addendum,
//! pressure angle, profile/tip interference, backlash, transmission
//! error, efficiency/losses, torque/power load ratings, component
//! tolerances, safety factors, fatigue/thermal limits, or design-code
//! compliance. Do not use it to certify a physical gearbox.
//!
//! ## Example
//!
//! ```
//! use valenx_epicyclic::{willis, TeethSet};
//!
//! // A 24-tooth sun, 24-tooth planets, 72-tooth ring (R = S + 2P).
//! let teeth = TeethSet::new(24, 24, 72).unwrap();
//!
//! // Ground the ring and drive the sun at 100 rev/min; the carrier is
//! // the output of this classic reducer.
//! let state = willis::solve_carrier_ring_fixed(&teeth, 100.0).unwrap();
//!
//! // Reduction ratio carrier/sun = S/(S+R) = 24/96 = 1/4, so 25 rev/min.
//! assert!((state.omega_carrier - 25.0).abs() < 1e-9);
//! assert!(state.omega_ring.abs() < 1e-12);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod tabular;
pub mod teeth;
pub mod willis;

pub use error::EpicyclicError;
pub use tabular::{Driver, TableRows};
pub use teeth::TeethSet;
pub use willis::TrainState;

#[cfg(test)]
mod tests {
    use super::*;

    /// Float comparison tolerance for exact-fraction ground truth.
    const EPS: f64 = 1e-9;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS
    }

    // ----------------------------------------------------------------
    // TeethSet construction + meshing constraint
    // ----------------------------------------------------------------

    #[test]
    fn teethset_accepts_valid_set() {
        let t = TeethSet::new(24, 24, 72).unwrap();
        assert_eq!((t.sun(), t.planet(), t.ring()), (24, 24, 72));
    }

    #[test]
    fn teethset_derives_ring_from_constraint() {
        // R = S + 2P = 18 + 2*21 = 60. Ground truth: a 18/21/60 set.
        let t = TeethSet::from_sun_planet(18, 21).unwrap();
        assert_eq!(t.ring(), 60);
        // Re-deriving via new() with that ring must agree.
        assert_eq!(t, TeethSet::new(18, 21, 60).unwrap());
    }

    #[test]
    fn teethset_rejects_non_positive() {
        let err = TeethSet::new(0, 24, 72).unwrap_err();
        assert_eq!(err.code(), "epicyclic.non_positive_teeth");
        assert!(TeethSet::new(24, -1, 72).is_err());
        assert!(TeethSet::from_sun_planet(24, 0).is_err());
    }

    #[test]
    fn teethset_rejects_bad_mesh() {
        // 24 + 2*24 = 72, so a ring of 70 must be rejected.
        let err = TeethSet::new(24, 24, 70).unwrap_err();
        assert_eq!(err.code(), "epicyclic.meshing_constraint");
        match err {
            EpicyclicError::MeshingConstraint { expected, .. } => assert_eq!(expected, 72),
            other => panic!("wrong variant: {other:?}"),
        }
    }

    // ----------------------------------------------------------------
    // Basic train ratio (carrier-fixed) ground truth
    // ----------------------------------------------------------------

    #[test]
    fn basic_train_ratio_is_minus_s_over_r() {
        // S=24, R=72 -> e = -1/3 (independent hand value).
        let t = TeethSet::new(24, 24, 72).unwrap();
        assert!(close(t.basic_train_ratio(), -1.0 / 3.0));

        // S=18, R=60 -> e = -0.30 exactly.
        let t2 = TeethSet::from_sun_planet(18, 21).unwrap();
        assert!(close(t2.basic_train_ratio(), -0.3));
    }

    #[test]
    fn planet_sun_ratio_carrier_fixed_ground_truth() {
        // S=24, P=24 -> planet/sun = -1 (equal pinions, opposite sense).
        let t = TeethSet::new(24, 24, 72).unwrap();
        assert!(close(t.planet_sun_ratio_carrier_fixed(), -1.0));
        // S=20, P=16 -> -20/16 = -1.25.
        let t2 = TeethSet::from_sun_planet(20, 16).unwrap();
        assert!(close(t2.planet_sun_ratio_carrier_fixed(), -1.25));
    }

    // ----------------------------------------------------------------
    // Named single-input reductions vs analytic closed forms
    // ----------------------------------------------------------------

    #[test]
    fn ring_fixed_reduction_quarter() {
        // Ring grounded, sun in, carrier out: ratio S/(S+R)=24/96=1/4.
        let t = TeethSet::new(24, 24, 72).unwrap();
        let s = willis::solve_carrier_ring_fixed(&t, 100.0).unwrap();
        assert!(close(s.omega_carrier, 25.0));
        assert!(close(s.omega_ring, 0.0));
        assert!(close(s.omega_sun, 100.0));
        // Independent closed form.
        let expected = 100.0 * (t.sun() as f64) / ((t.sun() + t.ring()) as f64);
        assert!(close(s.omega_carrier, expected));
    }

    #[test]
    fn sun_fixed_reduction_three_quarter() {
        // Sun grounded, ring in, carrier out: ratio R/(S+R)=72/96=3/4.
        let t = TeethSet::new(24, 24, 72).unwrap();
        let s = willis::solve_carrier_sun_fixed(&t, 100.0).unwrap();
        assert!(close(s.omega_carrier, 75.0));
        assert!(close(s.omega_sun, 0.0));
        let expected = 100.0 * (t.ring() as f64) / ((t.sun() + t.ring()) as f64);
        assert!(close(s.omega_carrier, expected));
    }

    #[test]
    fn carrier_fixed_reverses_and_steps_down() {
        // Carrier grounded, sun in, ring out: ratio -S/R = -1/3.
        let t = TeethSet::new(24, 24, 72).unwrap();
        let s = willis::solve_ring_carrier_fixed(&t, 90.0).unwrap();
        assert!(close(s.omega_ring, -30.0));
        assert!(close(s.omega_carrier, 0.0));
        assert!(close(s.omega_ring, 90.0 * t.basic_train_ratio()));
    }

    // ----------------------------------------------------------------
    // General Willis solver: round-trip + limiting cases
    // ----------------------------------------------------------------

    #[test]
    fn willis_solvers_are_mutually_consistent() {
        // Pick an arbitrary sun and ring speed, derive the carrier that
        // satisfies Willis, then require the other two solvers to recover
        // the original sun and ring from that consistent triple.
        let t = TeethSet::from_sun_planet(20, 16).unwrap(); // R=52
        let (ws, wr) = (130.0_f64, -10.0_f64);
        let wc = willis::solve_carrier(&t, ws, wr).unwrap();

        let ws_back = willis::solve_sun(&t, wr, wc).unwrap();
        let wr_back = willis::solve_ring(&t, ws, wc).unwrap();
        assert!(close(ws_back, ws));
        assert!(close(wr_back, wr));
    }

    #[test]
    fn willis_locked_train_spins_rigidly() {
        // If all three speeds are equal the relative motion is zero and
        // the Willis constraint is satisfied for any e: a direct-drive
        // (locked) planetary. Solve the ring from sun=carrier=55.
        let t = TeethSet::new(24, 24, 72).unwrap();
        let wr = willis::solve_ring(&t, 55.0, 55.0).unwrap();
        assert!(close(wr, 55.0));
    }

    #[test]
    fn willis_satisfies_its_own_constraint() {
        // Resolve a carrier, then plug back into the raw linear form
        // omega_ring - e*omega_sun + (e-1)*omega_carrier = 0.
        let t = TeethSet::from_sun_planet(30, 15).unwrap(); // R=60
        let e = t.basic_train_ratio();
        let (ws, wr) = (200.0, 50.0);
        let wc = willis::solve_carrier(&t, ws, wr).unwrap();
        let residual = wr - e * ws + (e - 1.0) * wc;
        assert!(residual.abs() < 1e-9);
    }

    // ----------------------------------------------------------------
    // Planet spin
    // ----------------------------------------------------------------

    #[test]
    fn planet_spin_ground_truth_ring_fixed() {
        // Ring fixed, sun=96 -> carrier=24. Planet carrier-frame spin =
        // (96-24)*(-S/P)=(72)*(-1)=-72; absolute = -72 + carrier 24 = -48.
        let t = TeethSet::new(24, 24, 72).unwrap();
        let s = willis::solve_carrier_ring_fixed(&t, 96.0).unwrap();
        assert!(close(s.omega_carrier, 24.0));
        assert!(close(s.planet_spin(&t), -48.0));
    }

    // ----------------------------------------------------------------
    // Tabular (superposition) method
    // ----------------------------------------------------------------

    #[test]
    fn tabular_matches_ring_fixed_closed_form() {
        // Drive sun=100 with carrier unknown is not a tabular input; but
        // ring-fixed means we can supply carrier=25 (the answer) and the
        // sun driver=100, and the table must reproduce ring=0.
        let t = TeethSet::new(24, 24, 72).unwrap();
        let rows = tabular::resolve(&t, 25.0, Driver::Sun, 100.0).unwrap();
        assert!(close(rows.state.omega_ring, 0.0));
        assert!(close(rows.state.omega_sun, 100.0));
        assert!(close(rows.state.omega_carrier, 25.0));
        // x is the carrier; y is sun - carrier.
        assert!(close(rows.x, 25.0));
        assert!(close(rows.y, 75.0));
    }

    #[test]
    fn tabular_equals_willis_for_random_state() {
        // For an arbitrary (carrier, sun) pair the tabular ring must match
        // the Willis ring solver exactly.
        let t = TeethSet::from_sun_planet(20, 16).unwrap(); // R=52
        let (wc, ws) = (40.0, 130.0);
        let rows = tabular::resolve(&t, wc, Driver::Sun, ws).unwrap();
        let wr_willis = willis::solve_ring(&t, ws, wc).unwrap();
        assert!(close(rows.state.omega_ring, wr_willis));
        // And the planet column matches TrainState::planet_spin.
        let st = TrainState {
            omega_sun: ws,
            omega_ring: wr_willis,
            omega_carrier: wc,
        };
        assert!(close(rows.omega_planet, st.planet_spin(&t)));
    }

    #[test]
    fn tabular_ring_driver_inverts_consistently() {
        // Drive from the ring instead of the sun: carrier=40, ring=10.
        // Solve y=(ring-x)/e, then the sun must match Willis solve_sun.
        let t = TeethSet::from_sun_planet(24, 18).unwrap(); // R=60
        let (wc, wr) = (40.0, 10.0);
        let rows = tabular::resolve(&t, wc, Driver::Ring, wr).unwrap();
        let ws_willis = willis::solve_sun(&t, wr, wc).unwrap();
        assert!(close(rows.state.omega_sun, ws_willis));
        assert!(close(rows.state.omega_ring, wr));
    }

    // ----------------------------------------------------------------
    // Input validation on the solvers
    // ----------------------------------------------------------------

    #[test]
    fn solvers_reject_non_finite() {
        let t = TeethSet::new(24, 24, 72).unwrap();
        assert_eq!(
            willis::solve_carrier(&t, f64::NAN, 0.0).unwrap_err().code(),
            "epicyclic.non_finite"
        );
        assert!(willis::solve_sun(&t, f64::INFINITY, 0.0).is_err());
        assert!(willis::solve_ring(&t, 0.0, f64::NEG_INFINITY).is_err());
        assert!(willis::solve_carrier_ring_fixed(&t, f64::NAN).is_err());
        assert!(tabular::resolve(&t, f64::NAN, Driver::Sun, 1.0).is_err());
        assert!(tabular::resolve(&t, 1.0, Driver::Ring, f64::INFINITY).is_err());
    }

    #[test]
    fn error_codes_are_stable_and_distinct() {
        let a = EpicyclicError::NonPositiveTeeth {
            name: "sun",
            value: 0,
        };
        let b = EpicyclicError::NonFinite { name: "omega_sun" };
        let c = EpicyclicError::Singular { reason: "x" };
        assert_eq!(a.code(), "epicyclic.non_positive_teeth");
        assert_eq!(b.code(), "epicyclic.non_finite");
        assert_eq!(c.code(), "epicyclic.singular");
        assert_ne!(a.code(), b.code());
    }
}
