//! # valenx-bolt — bolted-joint analysis
//!
//! A small, dependency-light core for the classic textbook analysis of a
//! **preloaded bolted joint**: how hard you have to tighten a bolt to
//! reach a target clamp force, how an external service load is shared
//! between the bolt and the clamped members, and when the joint lets go.
//!
//! ## What
//!
//! Given a bolt diameter, a tightening torque (or a directly specified
//! preload), a nut factor and the bolt/member stiffness split, this
//! crate answers:
//!
//! - **Preload from torque** — [`BoltedJoint::from_torque`] inverts
//!   `T = K F d` to the achieved preload `F`, and
//!   [`BoltedJoint::required_torque_nm`] runs it forward.
//! - **Load sharing** — [`BoltedJoint::bolt_load_increment_n`] (`C P`),
//!   [`BoltedJoint::bolt_load_n`] (`F + C P`) and
//!   [`BoltedJoint::clamping_force_n`] (`F - (1 - C) P`).
//! - **Separation** — [`BoltedJoint::separation_load_n`] (`F / (1 - C)`),
//!   [`BoltedJoint::stays_clamped`] and
//!   [`BoltedJoint::separation_safety_factor`].
//! - **Strength** — [`stress`] turns an ISO-metric thread into a
//!   tensile-stress area `A_t` and, with a [`material::BoltMaterial`]
//!   grade, into proof / tensile loads and axial stress.
//!
//! ## Model
//!
//! The mechanics are the standard linear elastic bolted-joint model
//! (Shigley, VDI 2230 in its first-principles form):
//!
//! - **Torque ↔ preload:** `T = K F d`, with the nut factor `K` lumping
//!   thread and under-head friction and the helix geometry into one
//!   dimensionless coefficient (`K ≈ 0.2` for plain steel). A larger `K`
//!   needs a larger torque for the same preload.
//! - **Stiffness split:** the bolt (`kb`) and members (`km`) are springs
//!   in parallel; the joint constant `C = kb / (kb + km)` lies in
//!   `(0, 1)`. An external tensile load `P` raises the bolt tension by
//!   `C P` and relieves the member clamp by `(1 - C) P`.
//! - **Separation:** the members stay together while `F - (1 - C) P > 0`,
//!   i.e. up to `P_sep = F / (1 - C)`.
//! - **Stress area:** `A_t = (pi/4)(d - 0.938_194 P)^2` for ISO metric
//!   threads; loads are strength × `A_t`, stress is force / `A_t`.
//!
//! Everything is unit-agnostic provided the inputs are consistent. The
//! ergonomic path uses SI: metres, newtons, newton-metres and pascals.
//!
//! ```
//! use valenx_bolt::{BoltGrade, BoltedJoint, NutFactor, StiffnessRatio, stress};
//!
//! // An M10 class-8.8 bolt, P = 1.5 mm pitch, tightened to 40 N·m.
//! let d = 0.010_f64; // 10 mm
//! let area = stress::tensile_stress_area(0.010, 0.0015).unwrap();
//! let material = BoltGrade::Class8_8.material().unwrap();
//!
//! let k = NutFactor::new(NutFactor::STEEL_AS_RECEIVED).unwrap();
//! let c = StiffnessRatio::new(0.25).unwrap();
//! let joint = BoltedJoint::from_torque(40.0, k, d, c).unwrap();
//!
//! // Preload F = T / (K d) = 40 / (0.2 * 0.01) = 20 000 N.
//! assert!((joint.preload_n() - 20_000.0).abs() < 1e-6);
//!
//! // Stays clamped well below the separation load, and below proof.
//! assert!(joint.stays_clamped(5_000.0).unwrap());
//! assert!(joint.bolt_load_n(5_000.0).unwrap()
//!     < stress::proof_load(&material, area).unwrap());
//! ```
//!
//! ## Honest scope
//!
//! This is a **research/educational-grade** library of closed-form,
//! textbook bolted-joint formulas. It is **NOT a clinical/medical tool
//! and NOT a production engineering certification tool.** In particular:
//!
//! - The nut factor `K` and the stiffness constant `C` are **inputs**.
//!   This crate does not derive `K` from a friction model, nor `C` from
//!   a member-stiffness frustum / `km` calculation — you supply them (or
//!   compute `C` from stiffnesses you bring yourself via
//!   [`StiffnessRatio::from_stiffnesses`]).
//! - The model is **static and linear elastic**. There is no fatigue /
//!   endurance (Goodman) analysis, no thread stripping or pull-out
//!   check, no bending or eccentric-load / prying, no thermal
//!   relaxation, embedding or creep, and no gasketed-joint behaviour.
//! - Material grades are **nominal table values** (ISO 898-1), not
//!   measured lot certificates.
//!
//! Use it to learn, prototype and sanity-check — not to release hardware
//! without an independent, standards-based verification.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod joint;
pub mod material;
pub mod stress;

pub use error::{BoltError, ErrorCategory};
pub use joint::{BoltedJoint, NutFactor, StiffnessRatio};
pub use material::{BoltGrade, BoltMaterial};

#[cfg(test)]
mod tests {
    use super::*;

    /// Tight tolerance for exact closed-form identities.
    const EPS: f64 = 1e-9;

    // --- Preload from torque: F = T / (K d) -------------------------

    #[test]
    fn preload_matches_torque_relation() {
        // T = 50 N·m, K = 0.2, d = 12 mm → F = 50 / (0.2 * 0.012).
        let k = NutFactor::new(0.2).unwrap();
        let c = StiffnessRatio::new(0.3).unwrap();
        let joint = BoltedJoint::from_torque(50.0, k, 0.012, c).unwrap();
        let expected = 50.0 / (0.2 * 0.012);
        assert!((joint.preload_n() - expected).abs() < EPS * expected);
    }

    #[test]
    fn torque_roundtrips_through_preload() {
        // from_torque then required_torque_nm must recover the torque.
        let k = NutFactor::new(0.18).unwrap();
        let c = StiffnessRatio::new(0.25).unwrap();
        let torque = 33.7_f64;
        let joint = BoltedJoint::from_torque(torque, k, 0.01, c).unwrap();
        let recovered = joint.required_torque_nm(k);
        assert!((recovered - torque).abs() < EPS * torque);
    }

    #[test]
    fn larger_nut_factor_needs_more_torque_for_same_preload() {
        // Fix the preload F; the torque to reach it scales linearly with
        // K: T = K F d. A bigger K (more friction lost) ⇒ bigger torque.
        let c = StiffnessRatio::new(0.25).unwrap();
        let joint = BoltedJoint::with_preload(20_000.0, 0.01, c).unwrap();

        let k_low = NutFactor::new(0.12).unwrap();
        let k_high = NutFactor::new(0.30).unwrap();
        let t_low = joint.required_torque_nm(k_low);
        let t_high = joint.required_torque_nm(k_high);

        assert!(t_high > t_low);
        // Ratio of torques equals the ratio of nut factors exactly.
        assert!((t_high / t_low - 0.30 / 0.12).abs() < EPS);
    }

    #[test]
    fn lower_nut_factor_gives_more_preload_for_same_torque() {
        // Conversely, at fixed torque a lower K yields a higher preload.
        let c = StiffnessRatio::new(0.25).unwrap();
        let d = 0.01;
        let torque = 40.0;
        let f_low_k = BoltedJoint::from_torque(torque, NutFactor::new(0.15).unwrap(), d, c)
            .unwrap()
            .preload_n();
        let f_high_k = BoltedJoint::from_torque(torque, NutFactor::new(0.25).unwrap(), d, c)
            .unwrap()
            .preload_n();
        assert!(f_low_k > f_high_k);
    }

    // --- Stiffness ratio C lives in (0, 1) --------------------------

    #[test]
    fn stiffness_ratio_from_stiffnesses_is_in_open_unit_interval() {
        // For any positive kb, km, C = kb/(kb+km) ∈ (0, 1).
        for &(kb, km) in &[(1.0, 1.0), (1.0, 9.0), (7.0, 0.001), (0.001, 7.0)] {
            let c = StiffnessRatio::from_stiffnesses(kb, km).unwrap();
            assert!(
                c.value() > 0.0 && c.value() < 1.0,
                "C out of range for {kb},{km}"
            );
            let expected = kb / (kb + km);
            assert!((c.value() - expected).abs() < EPS);
        }
    }

    #[test]
    fn equal_stiffnesses_split_load_in_half() {
        let c = StiffnessRatio::from_stiffnesses(5.0, 5.0).unwrap();
        assert!((c.value() - 0.5).abs() < EPS);
        assert!((c.member_fraction() - 0.5).abs() < EPS);
    }

    #[test]
    fn member_fraction_complements_c() {
        let c = StiffnessRatio::new(0.27).unwrap();
        assert!((c.value() + c.member_fraction() - 1.0).abs() < EPS);
    }

    // --- Bolt load increment = C * P --------------------------------

    #[test]
    fn bolt_load_increment_is_c_times_p() {
        let c = StiffnessRatio::new(0.25).unwrap();
        let joint = BoltedJoint::with_preload(20_000.0, 0.01, c).unwrap();
        let p = 8_000.0;
        let inc = joint.bolt_load_increment_n(p).unwrap();
        assert!((inc - 0.25 * p).abs() < EPS * p);
    }

    #[test]
    fn total_bolt_load_is_preload_plus_share() {
        let c = StiffnessRatio::new(0.25).unwrap();
        let f = 20_000.0;
        let joint = BoltedJoint::with_preload(f, 0.01, c).unwrap();
        let p = 8_000.0;
        let total = joint.bolt_load_n(p).unwrap();
        assert!((total - (f + 0.25 * p)).abs() < EPS * f);
    }

    #[test]
    fn zero_external_load_leaves_bolt_at_preload() {
        let c = StiffnessRatio::new(0.4).unwrap();
        let f = 12_345.0;
        let joint = BoltedJoint::with_preload(f, 0.01, c).unwrap();
        assert!((joint.bolt_load_n(0.0).unwrap() - f).abs() < EPS * f);
        assert!((joint.clamping_force_n(0.0).unwrap() - f).abs() < EPS * f);
    }

    // --- Clamping force and separation ------------------------------

    #[test]
    fn clamping_force_drops_by_member_share() {
        let c = StiffnessRatio::new(0.25).unwrap();
        let f = 20_000.0;
        let joint = BoltedJoint::with_preload(f, 0.01, c).unwrap();
        let p = 8_000.0;
        let clamp = joint.clamping_force_n(p).unwrap();
        // F - (1 - C) P = 20000 - 0.75 * 8000 = 14000.
        assert!((clamp - (f - 0.75 * p)).abs() < EPS * f);
        assert!((clamp - 14_000.0).abs() < EPS * f);
    }

    #[test]
    fn separation_load_is_preload_over_one_minus_c() {
        let c = StiffnessRatio::new(0.25).unwrap();
        let f = 20_000.0;
        let joint = BoltedJoint::with_preload(f, 0.01, c).unwrap();
        let p_sep = joint.separation_load_n();
        // 20000 / 0.75 = 26666.67.
        assert!((p_sep - f / 0.75).abs() < EPS * p_sep);
    }

    #[test]
    fn clamping_force_is_zero_at_separation_load() {
        // At P = P_sep the residual clamp must be exactly zero.
        let c = StiffnessRatio::new(0.3).unwrap();
        let f = 15_000.0;
        let joint = BoltedJoint::with_preload(f, 0.01, c).unwrap();
        let p_sep = joint.separation_load_n();
        let clamp = joint.clamping_force_n(p_sep).unwrap();
        assert!(clamp.abs() < 1e-6, "residual clamp at P_sep was {clamp}");
    }

    #[test]
    fn joint_separates_above_separation_load() {
        let c = StiffnessRatio::new(0.25).unwrap();
        let f = 20_000.0;
        let joint = BoltedJoint::with_preload(f, 0.01, c).unwrap();
        let p_sep = joint.separation_load_n();

        // Just below: clamped, positive residual.
        let below = p_sep * 0.999;
        assert!(joint.stays_clamped(below).unwrap());
        assert!(joint.clamping_force_n(below).unwrap() > 0.0);

        // Just above: separated, negative residual.
        let above = p_sep * 1.001;
        assert!(!joint.stays_clamped(above).unwrap());
        assert!(joint.clamping_force_n(above).unwrap() < 0.0);
    }

    #[test]
    fn separation_safety_factor_is_psep_over_p() {
        let c = StiffnessRatio::new(0.25).unwrap();
        let f = 20_000.0;
        let joint = BoltedJoint::with_preload(f, 0.01, c).unwrap();
        let p = 10_000.0;
        let n = joint.separation_safety_factor(p).unwrap();
        assert!((n - joint.separation_load_n() / p).abs() < EPS * n);
        // n > 1 ⇔ stays clamped.
        assert!(n > 1.0);
        assert!(joint.stays_clamped(p).unwrap());
    }

    // --- Tensile-stress area & strength -----------------------------

    #[test]
    fn tensile_stress_area_matches_iso_table_m10() {
        // M10 coarse, P = 1.5 mm. ISO 898-1 tabulates A_t = 58.0 mm².
        // Work in millimetres so the result is in mm².
        let at = stress::tensile_stress_area(10.0, 1.5).unwrap();
        assert!(
            (at - 58.0).abs() < 0.3,
            "A_t(M10) = {at} mm², expected ≈ 58"
        );
    }

    #[test]
    fn tensile_stress_area_matches_iso_table_m8() {
        // M8 coarse, P = 1.25 mm. ISO 898-1: A_t = 36.6 mm².
        let at = stress::tensile_stress_area(8.0, 1.25).unwrap();
        assert!(
            (at - 36.6).abs() < 0.3,
            "A_t(M8) = {at} mm², expected ≈ 36.6"
        );
    }

    #[test]
    fn axial_stress_is_force_over_area() {
        let s = stress::axial_stress(10_000.0, 58.0e-6).unwrap();
        assert!((s - 10_000.0 / 58.0e-6).abs() < 1e-3 * s);
    }

    #[test]
    fn proof_load_is_strength_times_area() {
        // Class 8.8: S_p = 600 MPa. A_t = 58 mm² = 58e-6 m².
        // F_p = 600e6 * 58e-6 = 34 800 N.
        let m = BoltGrade::Class8_8.material().unwrap();
        let fp = stress::proof_load(&m, 58.0e-6).unwrap();
        assert!((fp - 34_800.0).abs() < 1e-3 * fp, "F_p = {fp} N");
    }

    #[test]
    fn tensile_load_exceeds_proof_load() {
        let m = BoltGrade::Class10_9.material().unwrap();
        let area = 58.0e-6;
        let fp = stress::proof_load(&m, area).unwrap();
        let fu = stress::tensile_load(&m, area).unwrap();
        assert!(fu > fp);
    }

    #[test]
    fn recommended_preload_is_75_percent_of_proof() {
        let m = BoltGrade::Class8_8.material().unwrap();
        let area = 58.0e-6;
        let rp = stress::recommended_preload(&m, area).unwrap();
        let fp = stress::proof_load(&m, area).unwrap();
        assert!((rp - 0.75 * fp).abs() < 1e-6 * fp);
    }

    // --- Material grades --------------------------------------------

    #[test]
    fn named_grades_have_proof_below_tensile() {
        for grade in [
            BoltGrade::Class4_6,
            BoltGrade::Class5_8,
            BoltGrade::Class8_8,
            BoltGrade::Class10_9,
            BoltGrade::Class12_9,
        ] {
            let m = grade.material().unwrap();
            assert!(
                m.proof_strength_pa < m.tensile_strength_pa,
                "grade {} proof !< tensile",
                grade.label()
            );
        }
    }

    #[test]
    fn custom_grade_validates_proof_le_tensile() {
        // Proof above ultimate is impossible and must be rejected.
        let bad = BoltMaterial::new(900.0e6, 800.0e6);
        assert!(bad.is_err());
        // A sane custom grade is accepted.
        let good = BoltMaterial::new(640.0e6, 800.0e6);
        assert!(good.is_ok());
    }

    // --- Error / validation paths -----------------------------------

    #[test]
    fn nut_factor_rejects_out_of_range() {
        assert!(matches!(
            NutFactor::new(0.0),
            Err(BoltError::NutFactorRange { .. })
        ));
        assert!(matches!(
            NutFactor::new(1.0),
            Err(BoltError::NutFactorRange { .. })
        ));
        assert!(matches!(
            NutFactor::new(1.5),
            Err(BoltError::NutFactorRange { .. })
        ));
        assert!(matches!(
            NutFactor::new(f64::NAN),
            Err(BoltError::NotFinite { .. })
        ));
        assert!(NutFactor::new(0.2).is_ok());
    }

    #[test]
    fn stiffness_ratio_rejects_out_of_range() {
        assert!(matches!(
            StiffnessRatio::new(0.0),
            Err(BoltError::StiffnessRatioRange { .. })
        ));
        assert!(matches!(
            StiffnessRatio::new(1.0),
            Err(BoltError::StiffnessRatioRange { .. })
        ));
        assert!(matches!(
            StiffnessRatio::from_stiffnesses(-1.0, 5.0),
            Err(BoltError::NonPositive { .. })
        ));
    }

    #[test]
    fn negative_torque_and_diameter_rejected() {
        let k = NutFactor::new(0.2).unwrap();
        let c = StiffnessRatio::new(0.25).unwrap();
        assert!(matches!(
            BoltedJoint::from_torque(-1.0, k, 0.01, c),
            Err(BoltError::NegativeLoad { .. })
        ));
        assert!(matches!(
            BoltedJoint::from_torque(10.0, k, 0.0, c),
            Err(BoltError::NonPositive { .. })
        ));
    }

    #[test]
    fn negative_external_load_rejected() {
        let c = StiffnessRatio::new(0.25).unwrap();
        let joint = BoltedJoint::with_preload(20_000.0, 0.01, c).unwrap();
        assert!(matches!(
            joint.bolt_load_n(-5.0),
            Err(BoltError::NegativeLoad { .. })
        ));
        assert!(matches!(
            joint.clamping_force_n(-5.0),
            Err(BoltError::NegativeLoad { .. })
        ));
    }

    #[test]
    fn degenerate_thread_rejected() {
        // Pitch larger than diameter ⇒ effective diameter non-positive.
        assert!(matches!(
            stress::tensile_stress_area(1.0, 2.0),
            Err(BoltError::NonPositive { .. })
        ));
    }

    #[test]
    fn error_code_and_category_are_stable() {
        let e = NutFactor::new(2.0).unwrap_err();
        assert_eq!(e.code(), "bolt.nut-factor-range");
        assert_eq!(e.category(), ErrorCategory::Input);
        let nf = BoltError::NotFinite {
            name: "x",
            value: f64::NAN,
        };
        assert_eq!(nf.category(), ErrorCategory::Numeric);
    }

    // --- Worked end-to-end example ----------------------------------

    #[test]
    fn worked_m10_class88_joint() {
        // M10 8.8, P=1.5mm, K=0.2, tightened to 40 N·m, C=0.25.
        let d = 0.010;
        let area_m2 = stress::tensile_stress_area(0.010, 0.0015).unwrap();
        let material = BoltGrade::Class8_8.material().unwrap();
        let k = NutFactor::new(0.2).unwrap();
        let c = StiffnessRatio::new(0.25).unwrap();
        let joint = BoltedJoint::from_torque(40.0, k, d, c).unwrap();

        // Preload = 40 / (0.2 * 0.01) = 20 000 N.
        assert!((joint.preload_n() - 20_000.0).abs() < 1e-6);

        // Proof load = 600 MPa * A_t.
        let fp = stress::proof_load(&material, area_m2).unwrap();
        // Preload is below proof (good practice).
        assert!(joint.preload_n() < fp);

        // Under a 6 kN service load the bolt picks up 0.25*6000 = 1500 N.
        let inc = joint.bolt_load_increment_n(6_000.0).unwrap();
        assert!((inc - 1_500.0).abs() < 1e-6);
        // Total bolt force still below proof.
        assert!(joint.bolt_load_n(6_000.0).unwrap() < fp);
        // And the joint is nowhere near separation.
        assert!(joint.stays_clamped(6_000.0).unwrap());
        assert!(joint.separation_safety_factor(6_000.0).unwrap() > 4.0);
    }
}
