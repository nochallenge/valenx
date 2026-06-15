//! # valenx-rivet
//!
//! Closed-form **riveted-joint strength** calculator. Describe a row (or
//! rows) of rivets connecting two plates loaded in tension, and read back
//! the three classic failure loads — rivet shear, plate bearing, and
//! plate tension on the net section — the **governing** (smallest) joint
//! strength, and the **efficiency** of the joint against the un-drilled
//! solid plate.
//!
//! ## What
//!
//! The entry point is [`Joint`], built from a [`RivetGroup`] (diameter,
//! how many rivets per row, how many rows, how many shear planes), a
//! [`Plate`] (gross width and thickness of the critical section) and a
//! set of [`Allowables`] (permissible shear, bearing and tensile
//! stress). [`Joint::analyze`] returns a [`JointResult`] with every
//! failure load, the governing [`FailureMode`], the governing strength
//! and the efficiency.
//!
//! ```
//! use valenx_rivet::{Allowables, Joint, Plate, RivetGroup};
//!
//! // One row of three 20 mm rivets in single shear, joining 10 mm
//! // plates 150 mm wide; mild-steel working stresses (MPa).
//! let group = RivetGroup::new(0.020, 3, 1, 1).unwrap();
//! let plate = Plate::new(0.150, 0.010).unwrap();
//! let allow = Allowables::new(80.0e6, 160.0e6, 100.0e6).unwrap();
//!
//! let r = Joint::new(group, plate, allow).analyze().unwrap();
//! println!(
//!     "governs in {:?} at {:.1} kN, efficiency {:.1} %",
//!     r.mode,
//!     r.strength / 1.0e3,
//!     r.efficiency * 100.0,
//! );
//! assert!(r.efficiency < 1.0);
//! ```
//!
//! ## Model
//!
//! Each strength is an allowable stress times the area that resists it,
//! the standard strength-of-materials treatment of a riveted joint (all
//! lengths metres, stresses pascals, forces newtons):
//!
//! - Rivet **shear**: `P_shear = n · s · (π d² / 4) · τ`, summing the
//!   shank cross-sections over `n` rivets and `s` shear planes (`s = 1`
//!   lap / single-cover butt, `s = 2` double-cover butt).
//! - Plate / rivet **bearing**: `P_bearing = n · d · t · σ_b`, the
//!   allowable crushing stress on the projected contact area `d · t`.
//! - Plate **tension** on the net section:
//!   `P_tension = (w − n_row · d) · t · σ_t`, the gross width less the
//!   holes punched in the critical row.
//!
//! The joint fails in whichever mode reaches its allowable first, so the
//! **joint strength** is `min(P_shear, P_bearing, P_tension)`, and the
//! **efficiency** is `η = P / P_solid` against the solid plate
//! `P_solid = w · t · σ_t`. Because drilling holes can only remove
//! material, `η` is always strictly less than one.
//!
//! ## Honest scope
//!
//! This is a **research / educational-grade** calculator implementing the
//! textbook closed-form rivet-joint formulae. It is deliberately a first-
//! principles strength model and is **not** a clinical/medical tool nor a
//! production, code-of-record structural-design tool. In particular it
//! does **not**:
//!
//! - apply any design code's safety factors, partial factors, edge- and
//!   pitch-distance rules, or workmanship/fabrication allowances (e.g.
//!   Eurocode 3, AISC, or a boiler-code joint-efficiency table);
//! - model **friction-grip / slip-critical** behaviour, pre-tension, or
//!   the difference between rivets and high-strength bolts;
//! - account for **eccentric** (moment-loaded) rivet groups, secondary
//!   bending, prying, fatigue, stress concentration around the hole, or
//!   load sharing that is not uniform across the rivets;
//! - distinguish driven-rivet hole clearance, hole vs. shank diameter, or
//!   tearing/shear-out at the plate edge as separate modes.
//!
//! Every number it returns is a genuine first-order strength from the
//! governing formula, suitable for learning, sizing studies and sanity
//! checks — not for certifying a real joint.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod joint;
pub mod material;

pub use error::{Result, RivetError};
pub use joint::{FailureMode, Joint, JointResult, Plate, RivetGroup};
pub use material::Allowables;

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    /// Absolute tolerance for force comparisons (newtons). The loads here
    /// are tens to hundreds of kN, so 1e-6 N is ~1e-11 relative — well
    /// inside f64 round-off but far tighter than any physical meaning.
    const EPS_N: f64 = 1.0e-6;
    /// Absolute tolerance for dimensionless ratios (efficiency).
    const EPS: f64 = 1.0e-12;

    /// A reusable mild-steel-ish allowable set: τ = 80, σ_b = 160,
    /// σ_t = 100 MPa.
    fn allow() -> Allowables {
        Allowables::new(80.0e6, 160.0e6, 100.0e6).unwrap()
    }

    // -- Ground-truth: the shear formula ----------------------------------

    #[test]
    fn shear_matches_closed_form() {
        // n rivets, s planes, area π d²/4, stress τ.
        let d = 0.020;
        let group = RivetGroup::new(d, 4, 1, 1).unwrap();
        let plate = Plate::new(0.300, 0.012).unwrap();
        let j = Joint::new(group, plate, allow());

        let area = PI * d * d / 4.0;
        let expected = 4.0 * 1.0 * area * 80.0e6;
        assert!(
            (j.shear_strength() - expected).abs() < EPS_N,
            "shear {} vs expected {expected}",
            j.shear_strength()
        );
    }

    #[test]
    fn double_shear_is_twice_single_shear() {
        // s = 2 (double-cover butt) carries exactly twice s = 1.
        let single = RivetGroup::new(0.018, 3, 1, 1).unwrap();
        let double = RivetGroup::new(0.018, 3, 1, 2).unwrap();
        let plate = Plate::new(0.200, 0.010).unwrap();

        let ps = Joint::new(single, plate, allow()).shear_strength();
        let pd = Joint::new(double, plate, allow()).shear_strength();
        assert!(
            (pd - 2.0 * ps).abs() < EPS_N,
            "double {pd} should be 2x single {ps}"
        );
    }

    #[test]
    fn shear_scales_with_diameter_squared() {
        // Doubling d quadruples the shear area, hence the shear load.
        let plate = Plate::new(0.500, 0.010).unwrap();
        let g1 = RivetGroup::new(0.010, 2, 1, 1).unwrap();
        let g2 = RivetGroup::new(0.020, 2, 1, 1).unwrap();
        let p1 = Joint::new(g1, plate, allow()).shear_strength();
        let p2 = Joint::new(g2, plate, allow()).shear_strength();
        assert!(
            (p2 - 4.0 * p1).abs() < EPS_N,
            "p2 {p2} should be 4x p1 {p1}"
        );
    }

    // -- Ground-truth: the bearing formula --------------------------------

    #[test]
    fn bearing_matches_closed_form() {
        // n · d · t · σ_b.
        let d = 0.020;
        let t = 0.012;
        let group = RivetGroup::new(d, 4, 1, 1).unwrap();
        let plate = Plate::new(0.300, t).unwrap();
        let j = Joint::new(group, plate, allow());

        let expected = 4.0 * d * t * 160.0e6;
        assert!(
            (j.bearing_strength() - expected).abs() < EPS_N,
            "bearing {} vs expected {expected}",
            j.bearing_strength()
        );
    }

    #[test]
    fn bearing_is_independent_of_shear_planes() {
        // Bearing depends on projected area d·t only, not on s.
        let plate = Plate::new(0.250, 0.011).unwrap();
        let g1 = RivetGroup::new(0.016, 3, 1, 1).unwrap();
        let g2 = RivetGroup::new(0.016, 3, 1, 2).unwrap();
        let b1 = Joint::new(g1, plate, allow()).bearing_strength();
        let b2 = Joint::new(g2, plate, allow()).bearing_strength();
        assert!((b1 - b2).abs() < EPS_N, "bearing must ignore shear planes");
    }

    // -- Ground-truth: tension on the net section -------------------------

    #[test]
    fn tension_matches_net_section_closed_form() {
        // (w − n_row·d) · t · σ_t.
        let d = 0.020;
        let w = 0.150;
        let t = 0.010;
        let n_row = 3;
        let group = RivetGroup::new(d, n_row, 1, 1).unwrap();
        let plate = Plate::new(w, t).unwrap();
        let j = Joint::new(group, plate, allow());

        let net_w = w - n_row as f64 * d; // 0.150 - 0.060 = 0.090 m
        let expected = net_w * t * 100.0e6;
        let got = j.tension_strength().unwrap();
        assert!(
            (got - expected).abs() < EPS_N,
            "tension {got} vs expected {expected}"
        );
    }

    #[test]
    fn net_section_non_positive_is_rejected() {
        // 5 holes of 30 mm across a 150 mm plate removes 150 mm → none left.
        let group = RivetGroup::new(0.030, 5, 1, 1).unwrap();
        let plate = Plate::new(0.150, 0.010).unwrap();
        let j = Joint::new(group, plate, allow());
        let err = j.tension_strength().unwrap_err();
        assert!(matches!(err, RivetError::NetSectionNonPositive { .. }));
        // analyze() must surface the same error, not silently produce NaN.
        assert!(j.analyze().is_err());
    }

    // -- Ground-truth: joint strength = minimum failure mode --------------

    #[test]
    fn joint_strength_is_minimum_of_the_three() {
        let group = RivetGroup::new(0.020, 3, 1, 1).unwrap();
        let plate = Plate::new(0.150, 0.010).unwrap();
        let j = Joint::new(group, plate, allow());
        let r = j.analyze().unwrap();

        let lo = r.shear.min(r.bearing).min(r.tension);
        assert!(
            (r.strength - lo).abs() < EPS_N,
            "strength {} should equal min {lo}",
            r.strength
        );
        // And the reported mode must actually own that minimum value.
        let reported = match r.mode {
            FailureMode::Shear => r.shear,
            FailureMode::Bearing => r.bearing,
            FailureMode::Tension => r.tension,
        };
        assert!((reported - r.strength).abs() < EPS_N);
    }

    #[test]
    fn shear_can_govern() {
        // Tiny rivets, generous plate + bearing allowable → shear is least.
        // d = 6 mm, single shear, but a wide thick strong plate.
        let group = RivetGroup::new(0.006, 2, 1, 1).unwrap();
        let plate = Plate::new(0.300, 0.030).unwrap();
        let allow = Allowables::new(50.0e6, 400.0e6, 250.0e6).unwrap();
        let r = Joint::new(group, plate, allow).analyze().unwrap();
        assert_eq!(r.mode, FailureMode::Shear);
    }

    #[test]
    fn tension_can_govern() {
        // Many large rivets in one row gut the net section, while fat
        // rivets and a high shear allowable keep shear/bearing high.
        let group = RivetGroup::new(0.040, 3, 1, 2).unwrap();
        let plate = Plate::new(0.150, 0.020).unwrap();
        let allow = Allowables::new(300.0e6, 600.0e6, 80.0e6).unwrap();
        let r = Joint::new(group, plate, allow).analyze().unwrap();
        assert_eq!(r.mode, FailureMode::Tension);
    }

    #[test]
    fn bearing_can_govern() {
        // Thin plate + low bearing allowable, but a high tensile and
        // shear allowable, drives bearing to be the weakest link.
        let group = RivetGroup::new(0.020, 2, 1, 2).unwrap();
        let plate = Plate::new(0.300, 0.004).unwrap();
        let allow = Allowables::new(400.0e6, 60.0e6, 400.0e6).unwrap();
        let r = Joint::new(group, plate, allow).analyze().unwrap();
        assert_eq!(r.mode, FailureMode::Bearing);
    }

    // -- Ground-truth: efficiency = joint / solid, strictly < 1 -----------

    #[test]
    fn efficiency_is_joint_over_solid() {
        let group = RivetGroup::new(0.020, 3, 1, 1).unwrap();
        let plate = Plate::new(0.150, 0.010).unwrap();
        let j = Joint::new(group, plate, allow());
        let r = j.analyze().unwrap();

        let expected = r.strength / (0.150 * 0.010 * 100.0e6);
        assert!(
            (r.efficiency - expected).abs() < EPS,
            "efficiency {} vs expected {expected}",
            r.efficiency
        );
    }

    #[test]
    fn efficiency_strictly_below_one() {
        // Across a spread of geometries the efficiency must stay in (0,1):
        // the net/bearing/shear strength can never beat the solid plate.
        for &d in &[0.010, 0.016, 0.020, 0.024] {
            for &n in &[1u32, 2, 3, 4] {
                let group = RivetGroup::new(d, n, 1, 2).unwrap();
                let plate = Plate::new(0.300, 0.012).unwrap();
                let r = Joint::new(group, plate, allow()).analyze().unwrap();
                assert!(
                    r.efficiency > 0.0 && r.efficiency < 1.0,
                    "d={d} n={n}: efficiency {} out of (0,1)",
                    r.efficiency
                );
            }
        }
    }

    #[test]
    fn tension_limited_efficiency_equals_net_over_gross() {
        // When tension governs, η reduces to the pure section ratio
        // (w − n·d)/w, independent of σ_t. Force tension to govern with a
        // low tensile allowable and strong rivets/bearing.
        let d = 0.020;
        let w = 0.150;
        let n_row = 3;
        let group = RivetGroup::new(d, n_row, 1, 2).unwrap();
        let plate = Plate::new(w, 0.020).unwrap();
        let allow = Allowables::new(300.0e6, 600.0e6, 90.0e6).unwrap();
        let r = Joint::new(group, plate, allow).analyze().unwrap();
        assert_eq!(r.mode, FailureMode::Tension);

        let section_ratio = (w - n_row as f64 * d) / w; // 0.090/0.150 = 0.6
        assert!(
            (r.efficiency - section_ratio).abs() < EPS,
            "η {} should equal section ratio {section_ratio}",
            r.efficiency
        );
    }

    // -- Ground-truth: more rivets -> higher capacity ---------------------

    #[test]
    fn more_rivets_raise_shear_and_bearing_capacity() {
        // Adding rows (no extra holes in the critical tension row) lifts
        // both shear and bearing strength monotonically.
        let plate = Plate::new(0.400, 0.012).unwrap();
        let mut last_shear = 0.0;
        let mut last_bearing = 0.0;
        for rows in 1u32..=5 {
            let group = RivetGroup::new(0.020, 3, rows, 1).unwrap();
            let j = Joint::new(group, plate, allow());
            let s = j.shear_strength();
            let b = j.bearing_strength();
            assert!(s > last_shear, "shear should grow with rows={rows}");
            assert!(b > last_bearing, "bearing should grow with rows={rows}");
            last_shear = s;
            last_bearing = b;
        }
    }

    #[test]
    fn adding_a_shear_governed_rivet_raises_joint_strength() {
        // In a shear-governed joint, adding a row raises the governing
        // strength (extra rivets share the load). Keep the tension row
        // fixed so the net section — and thus tension strength — is
        // unchanged, leaving shear/bearing to climb.
        let plate = Plate::new(0.400, 0.030).unwrap();
        let allow = Allowables::new(50.0e6, 400.0e6, 250.0e6).unwrap();

        let one_row = RivetGroup::new(0.010, 2, 1, 1).unwrap();
        let two_row = RivetGroup::new(0.010, 2, 2, 1).unwrap();
        let r1 = Joint::new(one_row, plate, allow).analyze().unwrap();
        let r2 = Joint::new(two_row, plate, allow).analyze().unwrap();
        assert_eq!(r1.mode, FailureMode::Shear);
        assert!(
            r2.strength > r1.strength,
            "two rows {} should beat one row {}",
            r2.strength,
            r1.strength
        );
    }

    // -- Validated constructors reject bad input --------------------------

    #[test]
    fn constructors_reject_non_positive() {
        assert!(matches!(
            RivetGroup::new(0.0, 1, 1, 1),
            Err(RivetError::NotPositive {
                name: "diameter",
                ..
            })
        ));
        assert!(matches!(
            RivetGroup::new(f64::NAN, 1, 1, 1),
            Err(RivetError::NotPositive { .. })
        ));
        assert!(matches!(
            Plate::new(-1.0, 0.01),
            Err(RivetError::NotPositive { name: "width", .. })
        ));
        assert!(matches!(
            Plate::new(0.1, f64::INFINITY),
            Err(RivetError::NotPositive {
                name: "thickness",
                ..
            })
        ));
        assert!(matches!(
            Allowables::new(0.0, 1.0, 1.0),
            Err(RivetError::NotPositive { name: "shear", .. })
        ));
    }

    #[test]
    fn constructors_reject_zero_counts() {
        assert!(matches!(
            RivetGroup::new(0.02, 0, 1, 1),
            Err(RivetError::ZeroCount {
                name: "rivets_per_row"
            })
        ));
        assert!(matches!(
            RivetGroup::new(0.02, 1, 0, 1),
            Err(RivetError::ZeroCount { name: "rows" })
        ));
        assert!(matches!(
            RivetGroup::new(0.02, 1, 1, 0),
            Err(RivetError::ZeroCount {
                name: "shear_planes"
            })
        ));
    }

    #[test]
    fn group_helpers_are_consistent() {
        let g = RivetGroup::new(0.020, 3, 2, 1).unwrap();
        assert_eq!(g.total_rivets(), 6);
        let area = PI * 0.020 * 0.020 / 4.0;
        assert!((g.rivet_area() - area).abs() < EPS);
    }

    #[test]
    fn result_round_trips_through_json() {
        // The public result struct serializes and deserializes through
        // JSON. Floats are compared with a tolerance, not `==`, because a
        // text round-trip need not reproduce the exact same bit pattern.
        let group = RivetGroup::new(0.020, 3, 1, 1).unwrap();
        let plate = Plate::new(0.150, 0.010).unwrap();
        let r = Joint::new(group, plate, allow()).analyze().unwrap();
        let text = serde_json::to_string(&r).unwrap();
        let back: JointResult = serde_json::from_str(&text).unwrap();

        assert!((r.shear - back.shear).abs() < EPS_N);
        assert!((r.bearing - back.bearing).abs() < EPS_N);
        assert!((r.tension - back.tension).abs() < EPS_N);
        assert!((r.solid - back.solid).abs() < EPS_N);
        assert!((r.strength - back.strength).abs() < EPS_N);
        assert!((r.efficiency - back.efficiency).abs() < EPS);
        assert_eq!(r.mode, back.mode);
    }
}
