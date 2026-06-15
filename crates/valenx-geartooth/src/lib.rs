//! # valenx-geartooth
//!
//! Gear-tooth bending strength: the classic **Lewis** equation plus
//! **AGMA** bending-stress basics.
//!
//! ## What
//!
//! Given a gear's transmitted load and geometry, estimate the bending
//! stress at the tooth root — the quantity that governs tooth-breakage
//! (bending-fatigue) failure. Two model tiers are provided:
//!
//! - [`lewis`] — the bare Lewis cantilever-beam equation
//!   `sigma = Wt / (F m Y)`, plus the kinematics that supply `Wt`
//!   (pitch-line velocity, load from power or torque).
//! - [`agma`] — the AGMA refinement that multiplies in overload,
//!   dynamic, size, load-distribution, and rim-thickness factors and
//!   replaces `Y` with the geometry factor `J = Y / Kf`.
//! - [`lewis_factor`] — the tabulated Lewis form factor `Y(N)` for
//!   20-degree full-depth involute teeth, with interpolation.
//!
//! ## Model
//!
//! The Lewis equation treats a tooth as a uniform-strength cantilever
//! loaded by the tangential force `Wt`. In SI metric-module form, with
//! `Wt` in newtons and the face width `F` and module `m` in
//! millimetres, the root bending stress comes out in megapascals:
//!
//! ```text
//! sigma = Wt / (F * m * Y)
//! ```
//!
//! The dimensionless form factor `Y` rises with tooth count (a tooth
//! with more, smaller teeth has a stockier root), so larger gears carry
//! more load at the same stress. The AGMA stress layers correction
//! factors on top:
//!
//! ```text
//! sigma_agma = Wt * Ko * Kv * Ks * (1 / (b m)) * (Kh * Kb / J)
//! ```
//!
//! ## Honest scope
//!
//! Research / educational grade. These are **textbook closed-form
//! models** (Shigley's *Mechanical Engineering Design*; AGMA 2001-D04
//! factor forms). They are first-pass sizing estimates only and ignore
//! effects a real rating must include (measured fillet geometry,
//! material S-N data, lubrication regime, thermal and dynamic coupling,
//! statistical load spectra). This crate is **NOT a clinical, medical,
//! or production engineering certification tool** and must not be used
//! to rate life-critical gearing. Validate against a qualified analysis
//! and the governing standard before relying on any number.
//!
//! ## Example
//!
//! ```
//! use valenx_geartooth::{lewis_bending_stress_for_teeth, LewisResult};
//!
//! // 20-tooth pinion, module 5 mm, 50 mm face width, 3500 N tangential.
//! let LewisResult { form_factor_y, bending_stress_mpa } =
//!     lewis_bending_stress_for_teeth(3500.0, 50.0, 5.0, 20).unwrap();
//!
//! // Shigley Table 14-2 gives Y = 0.322 for N = 20.
//! assert!((form_factor_y - 0.322).abs() < 1e-12);
//! // sigma = 3500 / (50 * 5 * 0.322) ~= 43.48 MPa.
//! assert!((bending_stress_mpa - 43.4783).abs() < 1e-3);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod agma;
pub mod error;
pub mod lewis;
pub mod lewis_factor;
pub mod spec;

pub use agma::{
    agma_bending_stress, dynamic_factor_kv, geometry_factor_j, AgmaFactors, QV_MAX, QV_MIN,
};
pub use error::{ErrorCategory, GearToothError};
pub use lewis::{
    lewis_bending_stress, lewis_bending_stress_for_teeth, lewis_bending_stress_of,
    pitch_line_velocity_m_per_s, tangential_load_from_power_n, tangential_load_from_torque_n,
    LewisResult,
};
pub use lewis_factor::{
    lewis_form_factor, table_len, table_row, MAX_TABULATED_TEETH, MIN_TABULATED_TEETH,
};
pub use spec::ToothLoad;

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute-difference tolerance for tight closed-form checks.
    const EPS: f64 = 1e-9;

    // ---- Lewis form factor: ground-truth table values -----------------

    #[test]
    fn lewis_factor_matches_shigley_table_exactly() {
        // Spot-check several rows of Shigley Table 14-2 (20 deg full
        // depth, load at highest point of single-tooth contact).
        let cases = [
            (12u32, 0.245),
            (20, 0.322),
            (30, 0.359),
            (50, 0.409),
            (100, 0.447),
            (400, 0.480),
        ];
        for (teeth, expected) in cases {
            let y = lewis_form_factor(teeth).unwrap();
            assert!(
                (y - expected).abs() < EPS,
                "Y({teeth}) = {y}, expected {expected}"
            );
        }
    }

    #[test]
    fn lewis_factor_rises_with_tooth_count() {
        // Monotone non-decreasing across the whole table.
        let mut prev = lewis_form_factor(MIN_TABULATED_TEETH).unwrap();
        for teeth in (MIN_TABULATED_TEETH + 1)..=MAX_TABULATED_TEETH {
            let y = lewis_form_factor(teeth).unwrap();
            assert!(
                y + EPS >= prev,
                "Y not monotone: Y({teeth}) = {y} < previous {prev}"
            );
            prev = y;
        }
        // Strictly larger at the extremes.
        let small = lewis_form_factor(12).unwrap();
        let large = lewis_form_factor(400).unwrap();
        assert!(large > small, "large {large} should exceed small {small}");
    }

    #[test]
    fn lewis_factor_interpolates_between_rows() {
        // 23 teeth sits midway between 22 (Y=0.331) and 24 (Y=0.337):
        // expect 0.334.
        let y = lewis_form_factor(23).unwrap();
        assert!((y - 0.334).abs() < EPS, "interpolated Y(23) = {y}");
        // The interpolated value lies strictly between its neighbours.
        let lo = lewis_form_factor(22).unwrap();
        let hi = lewis_form_factor(24).unwrap();
        assert!(lo < y && y < hi);
    }

    #[test]
    fn lewis_factor_saturates_above_table() {
        let rack = lewis_form_factor(MAX_TABULATED_TEETH).unwrap();
        let beyond = lewis_form_factor(10_000).unwrap();
        assert!(
            (beyond - rack).abs() < EPS,
            "saturation: {beyond} vs {rack}"
        );
    }

    #[test]
    fn lewis_factor_rejects_undercut_pinion() {
        let err = lewis_form_factor(MIN_TABULATED_TEETH - 1).unwrap_err();
        assert_eq!(err.code(), "geartooth.out_of_domain");
        assert_eq!(err.category(), ErrorCategory::Domain);
    }

    #[test]
    fn lewis_table_is_strictly_increasing() {
        // Both columns of the embedded table strictly ascend.
        for i in 1..table_len() {
            let (t0, y0) = table_row(i - 1).unwrap();
            let (t1, y1) = table_row(i).unwrap();
            assert!(t1 > t0, "teeth not ascending at row {i}: {t0} then {t1}");
            assert!(y1 > y0, "Y not ascending at row {i}: {y0} then {y1}");
        }
        assert!(table_row(table_len()).is_none());
    }

    // ---- Lewis bending equation: closed-form ground truth -------------

    #[test]
    fn lewis_stress_matches_hand_calculation() {
        // sigma = Wt / (F m Y). Pick round numbers: Wt=1000 N, F=20 mm,
        // m=2 mm, Y=0.25  ->  1000 / (20*2*0.25) = 1000/10 = 100 MPa.
        let sigma = lewis_bending_stress(1000.0, 20.0, 2.0, 0.25).unwrap();
        assert!((sigma - 100.0).abs() < EPS, "sigma = {sigma}");
    }

    #[test]
    fn lewis_stress_for_teeth_reports_consistent_factor() {
        let r = lewis_bending_stress_for_teeth(3500.0, 50.0, 5.0, 20).unwrap();
        assert!((r.form_factor_y - 0.322).abs() < EPS);
        // sigma = 3500 / (50 * 5 * 0.322).
        let expected = 3500.0 / (50.0 * 5.0 * 0.322);
        assert!((r.bending_stress_mpa - expected).abs() < EPS);
    }

    #[test]
    fn lewis_stress_rises_with_load() {
        let base = lewis_bending_stress(1000.0, 20.0, 2.0, 0.3).unwrap();
        let heavier = lewis_bending_stress(2000.0, 20.0, 2.0, 0.3).unwrap();
        assert!(heavier > base);
        // Stress is exactly linear in load: doubling Wt doubles sigma.
        assert!((heavier - 2.0 * base).abs() < EPS);
    }

    #[test]
    fn lewis_stress_falls_with_face_width() {
        let narrow = lewis_bending_stress(1000.0, 10.0, 2.0, 0.3).unwrap();
        let wide = lewis_bending_stress(1000.0, 40.0, 2.0, 0.3).unwrap();
        assert!(wide < narrow);
        // Inverse proportionality: 4x face width -> 1/4 stress.
        assert!((wide - narrow / 4.0).abs() < EPS);
    }

    #[test]
    fn lewis_stress_falls_with_module() {
        let small_m = lewis_bending_stress(1000.0, 20.0, 1.0, 0.3).unwrap();
        let big_m = lewis_bending_stress(1000.0, 20.0, 5.0, 0.3).unwrap();
        assert!(big_m < small_m);
        // Inverse proportionality in module too.
        assert!((big_m - small_m / 5.0).abs() < EPS);
    }

    #[test]
    fn lewis_stress_falls_as_more_teeth_raise_y() {
        // For fixed load/face/module, more teeth -> larger Y -> lower
        // stress (the form factor is in the denominator).
        let few = lewis_bending_stress_for_teeth(1000.0, 20.0, 2.0, 14).unwrap();
        let many = lewis_bending_stress_for_teeth(1000.0, 20.0, 2.0, 100).unwrap();
        assert!(many.form_factor_y > few.form_factor_y);
        assert!(many.bending_stress_mpa < few.bending_stress_mpa);
    }

    #[test]
    fn lewis_stress_rejects_nonpositive_inputs() {
        assert!(lewis_bending_stress(-1.0, 20.0, 2.0, 0.3).is_err());
        assert!(lewis_bending_stress(1000.0, 0.0, 2.0, 0.3).is_err());
        assert!(lewis_bending_stress(1000.0, 20.0, f64::NAN, 0.3).is_err());
        assert!(lewis_bending_stress(1000.0, 20.0, 2.0, f64::INFINITY).is_err());
    }

    #[test]
    fn tooth_load_path_agrees_with_direct_call() {
        let load = ToothLoad::new(3500.0, 50.0, 5.0, 20).unwrap();
        let via_bundle = lewis_bending_stress_of(&load).unwrap();
        let direct = lewis_bending_stress_for_teeth(3500.0, 50.0, 5.0, 20).unwrap();
        assert!((via_bundle.bending_stress_mpa - direct.bending_stress_mpa).abs() < EPS);
        // Pitch diameter convenience: module * teeth.
        assert!((load.pitch_diameter_mm() - 100.0).abs() < EPS);
    }

    #[test]
    fn tooth_load_rejects_bad_inputs() {
        assert!(ToothLoad::new(0.0, 50.0, 5.0, 20).is_err());
        assert!(ToothLoad::new(3500.0, 50.0, 5.0, 0).is_err());
    }

    // ---- Kinematics: load from power / torque, pitch-line velocity ----

    #[test]
    fn pitch_line_velocity_ground_truth() {
        // d = 100 mm, n = 1000 rpm.
        // V = pi * 0.1 m * (1000/60) rev/s = 5 pi / 3 m/s ~= 5.235988.
        let v = pitch_line_velocity_m_per_s(100.0, 1000.0).unwrap();
        let expected = 5.0 * std::f64::consts::PI / 3.0;
        assert!((v - expected).abs() < EPS, "V = {v}");
    }

    #[test]
    fn load_from_power_ground_truth() {
        // P = V * Wt. With V = 5 m/s, Wt should recover from P = 5000 W
        // as 1000 N.
        let wt = tangential_load_from_power_n(5000.0, 5.0).unwrap();
        assert!((wt - 1000.0).abs() < EPS, "Wt = {wt}");
    }

    #[test]
    fn load_from_torque_ground_truth() {
        // T = 100 N·m, d = 100 mm -> r = 0.05 m -> Wt = 100/0.05 = 2000 N.
        let wt = tangential_load_from_torque_n(100.0, 100.0).unwrap();
        assert!((wt - 2000.0).abs() < EPS, "Wt = {wt}");
    }

    #[test]
    fn power_torque_load_paths_are_self_consistent() {
        // A torque T at diameter d running at n rpm transmits
        // P = T * omega = T * (2 pi n / 60) watts, and the tangential
        // load from power-at-V must equal the load from torque-at-r.
        let d = 80.0;
        let n = 1500.0;
        let torque_nm = 42.0;
        let omega = 2.0 * std::f64::consts::PI * n / 60.0;
        let power_w = torque_nm * omega;
        let v = pitch_line_velocity_m_per_s(d, n).unwrap();

        let wt_power = tangential_load_from_power_n(power_w, v).unwrap();
        let wt_torque = tangential_load_from_torque_n(torque_nm, d).unwrap();
        assert!(
            (wt_power - wt_torque).abs() < 1e-6,
            "power path {wt_power} vs torque path {wt_torque}"
        );
    }

    // ---- AGMA dynamic factor, geometry factor, bending stress ---------

    #[test]
    fn dynamic_factor_is_at_least_one_and_grows_with_speed() {
        let slow = dynamic_factor_kv(7.0, 1.0).unwrap();
        let fast = dynamic_factor_kv(7.0, 20.0).unwrap();
        assert!(slow >= 1.0, "Kv(slow) = {slow}");
        assert!(fast > slow, "Kv should grow with speed: {fast} > {slow}");
    }

    #[test]
    fn dynamic_factor_approaches_one_for_higher_quality() {
        // At a fixed speed, a higher quality number Qv (better teeth)
        // gives a Kv closer to 1.
        let coarse = dynamic_factor_kv(6.0, 10.0).unwrap();
        let fine = dynamic_factor_kv(11.0, 10.0).unwrap();
        assert!(fine < coarse, "finer quality Kv {fine} < coarser {coarse}");
        assert!(fine >= 1.0);
    }

    #[test]
    fn dynamic_factor_matches_shigley_curve_fit() {
        // Recompute the closed form independently: Qv = 6, V = 10 m/s.
        let qv: f64 = 6.0;
        let v: f64 = 10.0;
        let b = 0.25 * (12.0 - qv).powf(2.0 / 3.0);
        let a = 50.0 + 56.0 * (1.0 - b);
        let expected = ((a + (200.0 * v).sqrt()) / a).powf(b);
        let kv = dynamic_factor_kv(qv, v).unwrap();
        assert!(
            (kv - expected).abs() < EPS,
            "Kv = {kv}, expected {expected}"
        );
    }

    #[test]
    fn dynamic_factor_rejects_out_of_range_quality() {
        assert!(dynamic_factor_kv(5.0, 10.0).is_err());
        assert!(dynamic_factor_kv(12.0, 10.0).is_err());
        assert!(dynamic_factor_kv(7.0, -1.0).is_err());
    }

    #[test]
    fn geometry_factor_is_y_over_kf() {
        // J = Y / Kf. Y=0.4, Kf=1.6 -> J = 0.25.
        let j = geometry_factor_j(0.4, 1.6).unwrap();
        assert!((j - 0.25).abs() < EPS, "J = {j}");
        // With Kf = 1 the geometry factor equals Y itself.
        let j_unity = geometry_factor_j(0.322, 1.0).unwrap();
        assert!((j_unity - 0.322).abs() < EPS);
        // A larger fillet concentration lowers J.
        assert!(geometry_factor_j(0.4, 2.0).unwrap() < geometry_factor_j(0.4, 1.5).unwrap());
    }

    #[test]
    fn geometry_factor_rejects_relieving_kf() {
        // Kf < 1 would mean the fillet relieves stress — physically
        // impossible, so it is rejected.
        assert!(geometry_factor_j(0.4, 0.9).is_err());
        assert!(geometry_factor_j(-0.1, 1.5).is_err());
    }

    #[test]
    fn agma_unity_factors_reduce_to_lewis_with_j() {
        // With every correction factor = 1, the AGMA stress is exactly
        // Wt / (b m J) — the Lewis equation with Y replaced by J.
        let wt = 3500.0;
        let b = 50.0;
        let m = 5.0;
        let j = 0.25;
        let agma = agma_bending_stress(wt, b, m, j, &AgmaFactors::unity()).unwrap();
        let lewis_like = lewis_bending_stress(wt, b, m, j).unwrap();
        assert!(
            (agma - lewis_like).abs() < EPS,
            "AGMA-unity {agma} should equal Lewis-with-J {lewis_like}"
        );
    }

    #[test]
    fn agma_stress_scales_with_each_factor() {
        let wt = 3500.0;
        let b = 50.0;
        let m = 5.0;
        let j = 0.3;
        let base = agma_bending_stress(wt, b, m, j, &AgmaFactors::unity()).unwrap();

        // Doubling the dynamic factor doubles the stress; same for the
        // overload, size, load-distribution, and rim factors.
        let kv2 = AgmaFactors::new(1.0, 2.0, 1.0, 1.0, 1.0).unwrap();
        let ko2 = AgmaFactors::new(2.0, 1.0, 1.0, 1.0, 1.0).unwrap();
        let ks2 = AgmaFactors::new(1.0, 1.0, 2.0, 1.0, 1.0).unwrap();
        let kh2 = AgmaFactors::new(1.0, 1.0, 1.0, 2.0, 1.0).unwrap();
        let kb2 = AgmaFactors::new(1.0, 1.0, 1.0, 1.0, 2.0).unwrap();
        for f in [kv2, ko2, ks2, kh2, kb2] {
            let doubled = agma_bending_stress(wt, b, m, j, &f).unwrap();
            assert!(
                (doubled - 2.0 * base).abs() < EPS,
                "factor doubling should double stress: {doubled} vs {}",
                2.0 * base
            );
        }
    }

    #[test]
    fn agma_stress_exceeds_lewis_when_factors_above_one() {
        // Real correction factors are >= 1, so the AGMA stress is never
        // below the equivalent Lewis-with-J stress.
        let wt = 3500.0;
        let b = 50.0;
        let m = 5.0;
        let j = 0.3;
        let factors = AgmaFactors::new(1.25, 1.4, 1.1, 1.3, 1.0).unwrap();
        let agma = agma_bending_stress(wt, b, m, j, &factors).unwrap();
        let lewis_like = lewis_bending_stress(wt, b, m, j).unwrap();
        assert!(agma > lewis_like, "AGMA {agma} should exceed {lewis_like}");

        // And the ratio equals the product of the factors.
        let product = 1.25 * 1.4 * 1.1 * 1.3 * 1.0;
        assert!((agma / lewis_like - product).abs() < 1e-6);
    }

    #[test]
    fn agma_factors_reject_below_unity() {
        assert!(AgmaFactors::new(0.9, 1.0, 1.0, 1.0, 1.0).is_err());
        assert!(AgmaFactors::new(1.0, 1.0, 1.0, 1.0, 0.5).is_err());
        assert!(AgmaFactors::new(1.0, f64::NAN, 1.0, 1.0, 1.0).is_err());
    }

    #[test]
    fn agma_stress_rejects_nonpositive_geometry() {
        let f = AgmaFactors::unity();
        assert!(agma_bending_stress(3500.0, 50.0, 5.0, 0.0, &f).is_err());
        assert!(agma_bending_stress(3500.0, 0.0, 5.0, 0.3, &f).is_err());
    }

    // ---- Error metadata -----------------------------------------------

    #[test]
    fn error_codes_and_categories_are_stable() {
        let bad = GearToothError::bad_parameter("x", "nope");
        assert_eq!(bad.code(), "geartooth.bad_parameter");
        assert_eq!(bad.category(), ErrorCategory::Input);

        let dom = GearToothError::OutOfDomain("nope".into());
        assert_eq!(dom.code(), "geartooth.out_of_domain");
        assert_eq!(dom.category(), ErrorCategory::Domain);
    }

    // ---- Serde round-trip ---------------------------------------------

    #[test]
    fn result_serde_round_trips() {
        let r = lewis_bending_stress_for_teeth(3500.0, 50.0, 5.0, 20).unwrap();
        let json = serde_json::to_string(&r).unwrap();
        let back: LewisResult = serde_json::from_str(&json).unwrap();
        assert!((back.bending_stress_mpa - r.bending_stress_mpa).abs() < EPS);
        assert!((back.form_factor_y - r.form_factor_y).abs() < EPS);
    }
}
