//! # valenx-thickcylinder
//!
//! ## What
//!
//! A small, dependency-light solver for the stress field in a
//! thick-walled circular cylinder (a pressure-vessel wall, gun barrel,
//! hydraulic cylinder, pipe, etc.) loaded by a uniform internal and/or
//! external pressure. It computes the radial and hoop (circumferential)
//! stress at any radius through the wall, the location and value of the
//! peak hoop stress, and the elementary thin-wall membrane estimate used
//! to sanity-check the exact result.
//!
//! ## Model
//!
//! For an axisymmetric, long (plane-strain or open-ended) cylinder of
//! isotropic linear-elastic material with inner radius `a` and outer
//! radius `b`, the Lame solution gives the stresses as a constant term
//! plus an inverse-square term:
//!
//! ```text
//! sigma_r(r)     = A - B / r^2     (radial)
//! sigma_theta(r) = A + B / r^2     (hoop / circumferential)
//! ```
//!
//! The two integration constants follow from the pressure boundary
//! conditions `sigma_r(a) = -p_i`, `sigma_r(b) = -p_o`:
//!
//! ```text
//! A = (p_i a^2 - p_o b^2) / (b^2 - a^2)
//! B = a^2 b^2 (p_i - p_o) / (b^2 - a^2)
//! ```
//!
//! With internal pressure only (`p_o = 0`) the peak hoop stress is at the
//! bore `r = a`:
//!
//! ```text
//! sigma_theta,max = p_i (b^2 + a^2) / (b^2 - a^2)
//! ```
//!
//! As the wall becomes thin (`t = b - a`, `t/r -> 0`) this converges to
//! the membrane formula `sigma = p r_m / t` with mean radius
//! `r_m = (a + b)/2`, which the crate provides for comparison.
//!
//! ## Honest scope
//!
//! This is research/educational grade: standard textbook closed-form
//! linear-elasticity, validated against analytic ground truth. It is
//! NOT a clinical/medical or production-certified engineering tool. It
//! models none of the things a real pressure-vessel design must address:
//! no component tolerances, no safety factors, no yield/von-Mises check,
//! no fatigue or creep, no thermal stress, no autofrettage residual
//! stress, no end-cap / longitudinal stress, no material plasticity, and
//! no code compliance (ASME BPVC, EN 13445, PD 5500, ...). Do not size a
//! real vessel with it.
//!
//! ## Example
//!
//! ```rust
//! use valenx_thickcylinder::{Cylinder, LameConstants, hoop_stress, max_hoop_stress};
//!
//! // Bore 50 mm, outer 100 mm, 30 MPa internal pressure, no external.
//! let cyl = Cylinder::new(50.0, 100.0)?;
//! let k = LameConstants::solve(&cyl, 30.0, 0.0)?;
//!
//! // Peak hoop stress is at the bore and equals 50 MPa here.
//! let peak = max_hoop_stress(&cyl, &k)?;
//! assert!((peak - 50.0).abs() < 1e-9);
//!
//! // The same value via a direct evaluation at r = a.
//! let at_bore = hoop_stress(&cyl, &k, 50.0)?;
//! assert!((at_bore - peak).abs() < 1e-12);
//! # Ok::<(), valenx_thickcylinder::ThickCylinderError>(())
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod cylinder;
pub mod error;
pub mod stress;

pub use cylinder::{Cylinder, LameConstants};
pub use error::ThickCylinderError;
pub use stress::{hoop_stress, max_hoop_stress, radial_stress, thin_wall_hoop_stress};

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // --- geometry validation -------------------------------------------

    #[test]
    fn cylinder_accessors_match_inputs() {
        let cyl = Cylinder::new(50.0, 100.0).unwrap();
        assert!((cyl.inner() - 50.0).abs() < EPS);
        assert!((cyl.outer() - 100.0).abs() < EPS);
        assert!((cyl.thickness() - 50.0).abs() < EPS);
        assert!((cyl.mean_radius() - 75.0).abs() < EPS);
        assert!((cyl.diameter_ratio() - 2.0).abs() < EPS);
    }

    #[test]
    fn rejects_non_positive_radius() {
        let err = Cylinder::new(0.0, 100.0).unwrap_err();
        assert_eq!(err.code(), "thickcylinder.non-positive");
        let err = Cylinder::new(-5.0, 100.0).unwrap_err();
        assert_eq!(err.code(), "thickcylinder.non-positive");
    }

    #[test]
    fn rejects_outer_not_greater_than_inner() {
        // equal radii -> zero-thickness wall
        let err = Cylinder::new(100.0, 100.0).unwrap_err();
        assert_eq!(err.code(), "thickcylinder.bad-radii");
        // inverted
        let err = Cylinder::new(100.0, 50.0).unwrap_err();
        assert_eq!(err.code(), "thickcylinder.bad-radii");
    }

    #[test]
    fn rejects_non_finite_radius() {
        let err = Cylinder::new(f64::NAN, 100.0).unwrap_err();
        assert_eq!(err.code(), "thickcylinder.not-finite");
        let err = Cylinder::new(50.0, f64::INFINITY).unwrap_err();
        assert_eq!(err.code(), "thickcylinder.not-finite");
    }

    #[test]
    fn rejects_radius_outside_wall() {
        let cyl = Cylinder::new(50.0, 100.0).unwrap();
        let k = LameConstants::solve(&cyl, 30.0, 0.0).unwrap();
        let err = hoop_stress(&cyl, &k, 49.0).unwrap_err();
        assert_eq!(err.code(), "thickcylinder.outside-wall");
        let err = radial_stress(&cyl, &k, 101.0).unwrap_err();
        assert_eq!(err.code(), "thickcylinder.outside-wall");
        // NaN radius is rejected as not-finite, not outside-wall.
        let err = hoop_stress(&cyl, &k, f64::NAN).unwrap_err();
        assert_eq!(err.code(), "thickcylinder.not-finite");
    }

    #[test]
    fn rejects_negative_pressure() {
        let cyl = Cylinder::new(50.0, 100.0).unwrap();
        let err = LameConstants::solve(&cyl, -1.0, 0.0).unwrap_err();
        assert_eq!(err.code(), "thickcylinder.non-positive");
        let err = thin_wall_hoop_stress(&cyl, 0.0, -2.0).unwrap_err();
        assert_eq!(err.code(), "thickcylinder.non-positive");
    }

    // --- Lame constants vs hand-worked algebra --------------------------

    #[test]
    fn lame_constants_internal_pressure_only() {
        // a=50, b=100, p_i=30, p_o=0.
        // A = p_i a^2 / (b^2 - a^2) = 30*2500/7500 = 10.
        // B = a^2 b^2 p_i / (b^2 - a^2) = 2500*10000*30/7500 = 100000.
        let cyl = Cylinder::new(50.0, 100.0).unwrap();
        let k = LameConstants::solve(&cyl, 30.0, 0.0).unwrap();
        assert!((k.a - 10.0).abs() < EPS);
        assert!((k.b - 100_000.0).abs() < 1e-6);
    }

    // --- boundary conditions reproduced exactly ------------------------

    #[test]
    fn radial_stress_matches_pressure_at_boundaries() {
        // sigma_r(a) = -p_i, sigma_r(b) = -p_o.
        let cyl = Cylinder::new(50.0, 100.0).unwrap();
        let k = LameConstants::solve(&cyl, 30.0, 7.0).unwrap();
        let at_a = radial_stress(&cyl, &k, 50.0).unwrap();
        let at_b = radial_stress(&cyl, &k, 100.0).unwrap();
        assert!((at_a - (-30.0)).abs() < 1e-9);
        assert!((at_b - (-7.0)).abs() < 1e-9);
    }

    #[test]
    fn hoop_stress_max_at_bore_internal_only() {
        // Classic result: sigma_theta,max = p_i (b^2+a^2)/(b^2-a^2)
        //               = 30 * 12500 / 7500 = 50 MPa.
        let cyl = Cylinder::new(50.0, 100.0).unwrap();
        let k = LameConstants::solve(&cyl, 30.0, 0.0).unwrap();
        let peak = max_hoop_stress(&cyl, &k).unwrap();
        assert!((peak - 50.0).abs() < 1e-9);

        // At the outer surface: sigma_theta(b) = 20 MPa (A + B/b^2 = 10 + 10).
        let outer = hoop_stress(&cyl, &k, 100.0).unwrap();
        assert!((outer - 20.0).abs() < 1e-9);

        // Hoop stress decreases monotonically a -> b.
        let mid = hoop_stress(&cyl, &k, 75.0).unwrap();
        assert!(peak > mid && mid > outer);
    }

    #[test]
    fn stress_invariant_sum_is_constant_through_wall() {
        // sigma_r + sigma_theta = 2A is independent of r.
        let cyl = Cylinder::new(40.0, 90.0).unwrap();
        let k = LameConstants::solve(&cyl, 25.0, 3.0).unwrap();
        for &r in &[40.0_f64, 55.0, 70.0, 90.0] {
            let sr = radial_stress(&cyl, &k, r).unwrap();
            let st = hoop_stress(&cyl, &k, r).unwrap();
            assert!((sr + st - 2.0 * k.a).abs() < 1e-9);
        }
    }

    #[test]
    fn external_pressure_only_gives_compressive_hoop() {
        // External pressure only: peak compressive hoop at the bore,
        // sigma_theta(a) = -2 p_o b^2 / (b^2 - a^2).
        // a=50,b=100,p_o=12 -> -2*12*10000/7500 = -32 MPa.
        let cyl = Cylinder::new(50.0, 100.0).unwrap();
        let k = LameConstants::solve(&cyl, 0.0, 12.0).unwrap();
        let bore = hoop_stress(&cyl, &k, 50.0).unwrap();
        assert!((bore - (-32.0)).abs() < 1e-9);
        // radial at outer surface equals -p_o.
        let at_b = radial_stress(&cyl, &k, 100.0).unwrap();
        assert!((at_b - (-12.0)).abs() < 1e-9);
    }

    #[test]
    fn uniform_external_equals_internal_is_hydrostatic() {
        // p_i = p_o = p: B = 0, A = -p, so both stresses equal -p
        // everywhere (pure hydrostatic compression).
        let cyl = Cylinder::new(30.0, 80.0).unwrap();
        let k = LameConstants::solve(&cyl, 15.0, 15.0).unwrap();
        assert!(k.b.abs() < EPS);
        for &r in &[30.0_f64, 50.0, 80.0] {
            let sr = radial_stress(&cyl, &k, r).unwrap();
            let st = hoop_stress(&cyl, &k, r).unwrap();
            assert!((sr - (-15.0)).abs() < 1e-9);
            assert!((st - (-15.0)).abs() < 1e-9);
        }
    }

    // --- thin-wall limit recovery (the headline ground truth) ----------

    #[test]
    fn thin_wall_membrane_formula_value() {
        // sigma = p r_m / t. a=1000,b=1001 -> r_m=1000.5, t=1.
        // 10 * 1000.5 / 1 = 10005.0.
        let cyl = Cylinder::new(1000.0, 1001.0).unwrap();
        let membrane = thin_wall_hoop_stress(&cyl, 10.0, 0.0).unwrap();
        assert!((membrane - 10005.0).abs() < 1e-9);
    }

    #[test]
    fn exact_hoop_converges_to_thin_wall_as_wall_thins() {
        // As t/r -> 0 the bore hoop stress -> membrane estimate.
        // a=1000,b=1001,p_i=10: exact = 10*(b^2+a^2)/(b^2-a^2)
        //   = 10*2002001/2001 ~ 10005.0025; membrane = 10005.0.
        let cyl = Cylinder::new(1000.0, 1001.0).unwrap();
        let k = LameConstants::solve(&cyl, 10.0, 0.0).unwrap();
        let exact = max_hoop_stress(&cyl, &k).unwrap();
        let membrane = thin_wall_hoop_stress(&cyl, 10.0, 0.0).unwrap();
        // Relative error must be tiny for this slender wall.
        let rel = (exact - membrane).abs() / membrane;
        assert!(rel < 1e-6, "rel error {rel} too large");
    }

    #[test]
    fn thin_wall_error_grows_for_thick_wall() {
        // For a genuinely thick wall the membrane formula is a poor
        // approximation: the exact bore hoop stress should differ by
        // well over 10%, demonstrating why the full Lame model matters.
        let cyl = Cylinder::new(50.0, 100.0).unwrap();
        let k = LameConstants::solve(&cyl, 30.0, 0.0).unwrap();
        let exact = max_hoop_stress(&cyl, &k).unwrap(); // 50.0
        let membrane = thin_wall_hoop_stress(&cyl, 30.0, 0.0).unwrap();
        // membrane = 30 * 75 / 50 = 45.0  -> ~10% under-prediction.
        assert!((membrane - 45.0).abs() < 1e-9);
        let rel = (exact - membrane).abs() / exact;
        assert!(rel > 0.05, "thick-wall membrane error unexpectedly small");
    }

    // --- error code stability / display --------------------------------

    #[test]
    fn error_codes_are_stable_strings() {
        let e = Cylinder::new(0.0, 1.0).unwrap_err();
        assert_eq!(e.code(), "thickcylinder.non-positive");
        // Display is non-empty and mentions the parameter.
        assert!(e.to_string().contains("inner_radius"));
    }
}
