//! # valenx-beam
//!
//! Closed-form **Euler-Bernoulli beam bending** for the two classic
//! support cases — the **cantilever** and the **simply-supported** beam
//! — under a concentrated **point** load or a **uniformly distributed**
//! load (UDL).
//!
//! ## What
//!
//! Describe a beam as a [`Beam`] (span `L`, Young's modulus `E`, and a
//! cross-[`Section`]), pick a [`Support`] and a [`Load`], and
//! [`Beam::analyze`] returns the three engineering headline numbers:
//!
//! - the **maximum transverse deflection**,
//! - the **maximum bending moment**, and
//! - the **peak bending stress** they imply.
//!
//! Cross-sections are reduced to the only two scalars bending cares
//! about: the second moment of area `I` and the extreme-fibre distance
//! `c`. Two textbook closed forms ship out of the box — a solid
//! **rectangle** (`I = b*h^3 / 12`) and a solid **circle**
//! (`I = pi*d^4 / 64`).
//!
//! ```
//! use valenx_beam::{Beam, Load, Section, Support};
//!
//! // A 1 m steel bar (E = 200 GPa) of 20 x 30 mm rectangular section,
//! // 1 kN hung off the free tip of a cantilever. Units: N, mm, MPa.
//! let section = Section::rectangular(20.0, 30.0).unwrap();
//! let beam = Beam::new(1000.0, 200_000.0, section).unwrap();
//! let r = beam
//!     .analyze(Support::Cantilever, Load::Point { force: 1000.0 })
//!     .unwrap();
//!
//! // delta = P L^3 / (3 E I), with I = b h^3 / 12 = 45000 mm^4
//! assert!((r.max_deflection - 1000.0 * 1000.0_f64.powi(3)
//!         / (3.0 * 200_000.0 * 45_000.0)).abs() < 1e-9);
//! // M = P L = 1e6 N·mm, sigma = M c / I with c = 15 mm
//! assert!((r.max_moment - 1.0e6).abs() < 1e-6);
//! ```
//!
//! ## Model
//!
//! Each result is the standard **Euler-Bernoulli** closed form for a
//! prismatic, homogeneous, linear-elastic beam in small-deflection
//! bending:
//!
//! | Support          | Load           | Max deflection        | Max moment  |
//! |------------------|----------------|-----------------------|-------------|
//! | Cantilever       | tip point `P`  | `P L^3 / (3 E I)`     | `P L`       |
//! | Cantilever       | UDL `w`        | `w L^4 / (8 E I)`     | `w L^2 / 2` |
//! | Simply-supported | centre `P`     | `P L^3 / (48 E I)`    | `P L / 4`   |
//! | Simply-supported | UDL `w`        | `5 w L^4 / (384 E I)` | `w L^2 / 8` |
//!
//! and the peak bending stress in every case is the flexure formula
//! `sigma = M c / I`. Use one consistent unit system throughout (for
//! example N, mm and MPa).
//!
//! ## Honest scope
//!
//! Research / educational grade. These are **textbook closed-form**
//! Euler-Bernoulli expressions — the well-established, exact solutions
//! for the four idealised cases tabulated above — not a finite-element
//! solver and **NOT a clinical/medical or production structural-
//! engineering tool**. The underlying theory assumes:
//!
//! - small deflections (linear, geometrically undeformed equilibrium),
//! - a linear-elastic, homogeneous, isotropic, prismatic beam,
//! - slender geometry where **transverse shear is neglected**
//!   (pure bending; no Timoshenko shear correction), and
//! - the point load applied at the characteristic location (cantilever
//!   tip / simply-supported mid-span).
//!
//! Stress concentrations, buckling, plasticity, dynamic / fatigue
//! effects, support settlement, self-weight and arbitrary load
//! positions are out of scope. Do not use the output for the design of
//! a load-bearing structure without independent qualified review.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod beam;
pub mod error;
pub mod section;

pub use beam::{Beam, BeamResult, Load, Support};
pub use error::{BeamError, ErrorCategory};
pub use section::Section;

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Absolute tolerance for floating-point comparisons throughout the
    /// suite. The formulas are exact closed forms, so the only error is
    /// IEEE-754 rounding; a tight epsilon scaled to magnitude is used.
    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    // ---------------------------------------------------------------
    // Section: second moment of area, extreme fibre, modulus, area.
    // ---------------------------------------------------------------

    #[test]
    fn rectangular_second_moment_matches_formula() {
        // b = 20, h = 30 -> I = b h^3 / 12 = 20 * 27000 / 12 = 45000.
        let s = Section::rectangular(20.0, 30.0).unwrap();
        assert!(close(s.second_moment_area(), 45_000.0, 1e-9));
        assert!(close(s.extreme_fibre(), 15.0, 1e-12));
        assert!(close(s.area(), 600.0, 1e-12));
    }

    #[test]
    fn circular_second_moment_matches_formula() {
        // d = 10 -> I = pi d^4 / 64 = pi * 10000 / 64.
        let s = Section::circular(10.0).unwrap();
        let expected = PI * 10_000.0 / 64.0;
        assert!(close(s.second_moment_area(), expected, 1e-9));
        assert!(close(s.extreme_fibre(), 5.0, 1e-12));
        // Area = pi d^2 / 4 = 25 pi.
        assert!(close(s.area(), 25.0 * PI, 1e-9));
    }

    #[test]
    fn section_modulus_is_i_over_c() {
        // Rectangle S = b h^2 / 6 = 20 * 900 / 6 = 3000.
        let s = Section::rectangular(20.0, 30.0).unwrap();
        assert!(close(s.section_modulus().unwrap(), 3_000.0, 1e-9));
        assert!(close(
            s.section_modulus().unwrap(),
            s.second_moment_area() / s.extreme_fibre(),
            1e-9
        ));
    }

    #[test]
    fn rectangular_second_moment_scales_with_height_cubed() {
        // Doubling height multiplies I by 2^3 = 8.
        let base = Section::rectangular(20.0, 30.0).unwrap();
        let tall = Section::rectangular(20.0, 60.0).unwrap();
        assert!(close(
            tall.second_moment_area(),
            8.0 * base.second_moment_area(),
            1e-6
        ));
    }

    #[test]
    fn circular_second_moment_scales_with_diameter_to_the_fourth() {
        // Doubling diameter multiplies I by 2^4 = 16.
        let base = Section::circular(10.0).unwrap();
        let big = Section::circular(20.0).unwrap();
        assert!(close(
            big.second_moment_area(),
            16.0 * base.second_moment_area(),
            1e-6
        ));
    }

    #[test]
    fn section_constructors_reject_bad_input() {
        assert!(Section::rectangular(0.0, 10.0).is_err());
        assert!(Section::rectangular(10.0, -1.0).is_err());
        assert!(Section::circular(0.0).is_err());
        assert!(Section::circular(f64::INFINITY).is_err());
        assert!(Section::circular(f64::NAN).is_err());
    }

    // ---------------------------------------------------------------
    // Cantilever, tip point load: delta = P L^3 / (3 E I), M = P L.
    // ---------------------------------------------------------------

    #[test]
    fn cantilever_point_matches_textbook() {
        let section = Section::rectangular(20.0, 30.0).unwrap();
        let i = section.second_moment_area();
        let c = section.extreme_fibre();
        let (p, l, e) = (1000.0, 1000.0, 200_000.0);
        let beam = Beam::new(l, e, section).unwrap();
        let r = beam
            .analyze(Support::Cantilever, Load::Point { force: p })
            .unwrap();

        let expected_defl = p * l.powi(3) / (3.0 * e * i);
        let expected_moment = p * l;
        let expected_stress = expected_moment * c / i;

        assert!(close(r.max_deflection, expected_defl, 1e-9));
        assert!(close(r.max_moment, expected_moment, 1e-6));
        assert!(close(r.max_stress, expected_stress, 1e-9));
    }

    // ---------------------------------------------------------------
    // Simply-supported, centre point load: delta = P L^3 / (48 E I),
    // M = P L / 4.
    // ---------------------------------------------------------------

    #[test]
    fn simply_supported_point_matches_textbook() {
        let section = Section::circular(10.0).unwrap();
        let i = section.second_moment_area();
        let c = section.extreme_fibre();
        let (p, l, e) = (500.0, 800.0, 70_000.0);
        let beam = Beam::new(l, e, section).unwrap();
        let r = beam
            .analyze(Support::SimplySupported, Load::Point { force: p })
            .unwrap();

        let expected_defl = p * l.powi(3) / (48.0 * e * i);
        let expected_moment = p * l / 4.0;
        let expected_stress = expected_moment * c / i;

        assert!(close(r.max_deflection, expected_defl, 1e-9));
        assert!(close(r.max_moment, expected_moment, 1e-9));
        assert!(close(r.max_stress, expected_stress, 1e-9));
    }

    // ---------------------------------------------------------------
    // Simply-supported, UDL: delta = 5 w L^4 / (384 E I), M = w L^2 / 8.
    // ---------------------------------------------------------------

    #[test]
    fn simply_supported_udl_matches_textbook() {
        let section = Section::rectangular(50.0, 100.0).unwrap();
        let i = section.second_moment_area();
        let c = section.extreme_fibre();
        let (w, l, e) = (2.0, 3000.0, 12_000.0);
        let beam = Beam::new(l, e, section).unwrap();
        let r = beam
            .analyze(Support::SimplySupported, Load::Udl { intensity: w })
            .unwrap();

        let expected_defl = 5.0 * w * l.powi(4) / (384.0 * e * i);
        let expected_moment = w * l * l / 8.0;
        let expected_stress = expected_moment * c / i;

        assert!(close(r.max_deflection, expected_defl, 1e-6));
        assert!(close(r.max_moment, expected_moment, 1e-6));
        assert!(close(r.max_stress, expected_stress, 1e-9));
    }

    // ---------------------------------------------------------------
    // Cantilever, UDL: delta = w L^4 / (8 E I), M = w L^2 / 2.
    // ---------------------------------------------------------------

    #[test]
    fn cantilever_udl_matches_textbook() {
        let section = Section::rectangular(40.0, 40.0).unwrap();
        let i = section.second_moment_area();
        let (w, l, e) = (3.0, 1500.0, 100_000.0);
        let beam = Beam::new(l, e, section).unwrap();
        let r = beam
            .analyze(Support::Cantilever, Load::Udl { intensity: w })
            .unwrap();

        let expected_defl = w * l.powi(4) / (8.0 * e * i);
        let expected_moment = w * l * l / 2.0;

        assert!(close(r.max_deflection, expected_defl, 1e-6));
        assert!(close(r.max_moment, expected_moment, 1e-6));
    }

    // ---------------------------------------------------------------
    // Scaling laws (the prompt's explicit validations).
    // ---------------------------------------------------------------

    #[test]
    fn point_deflection_scales_with_length_cubed() {
        // Doubling L multiplies a point-load deflection by 2^3 = 8.
        let section = Section::rectangular(20.0, 30.0).unwrap();
        let load = Load::Point { force: 1000.0 };

        let short = Beam::new(1000.0, 200_000.0, section).unwrap();
        let long = Beam::new(2000.0, 200_000.0, section).unwrap();

        let d_short = short.max_deflection(Support::Cantilever, load).unwrap();
        let d_long = long.max_deflection(Support::Cantilever, load).unwrap();
        assert!(close(d_long, 8.0 * d_short, 1e-6));

        // Same cube law for the simply-supported centre-load case.
        let s_short = short
            .max_deflection(Support::SimplySupported, load)
            .unwrap();
        let s_long = long.max_deflection(Support::SimplySupported, load).unwrap();
        assert!(close(s_long, 8.0 * s_short, 1e-6));
    }

    #[test]
    fn udl_deflection_scales_with_length_to_the_fourth() {
        // Doubling L multiplies a UDL deflection by 2^4 = 16.
        let section = Section::rectangular(20.0, 30.0).unwrap();
        let load = Load::Udl { intensity: 2.0 };

        let short = Beam::new(1000.0, 200_000.0, section).unwrap();
        let long = Beam::new(2000.0, 200_000.0, section).unwrap();

        let c_short = short.max_deflection(Support::Cantilever, load).unwrap();
        let c_long = long.max_deflection(Support::Cantilever, load).unwrap();
        assert!(close(c_long, 16.0 * c_short, 1e-4));

        let s_short = short
            .max_deflection(Support::SimplySupported, load)
            .unwrap();
        let s_long = long.max_deflection(Support::SimplySupported, load).unwrap();
        assert!(close(s_long, 16.0 * s_short, 1e-4));
    }

    #[test]
    fn doubling_modulus_halves_deflection() {
        let section = Section::circular(12.0).unwrap();
        let load = Load::Point { force: 800.0 };

        let soft = Beam::new(1000.0, 100_000.0, section).unwrap();
        let stiff = Beam::new(1000.0, 200_000.0, section).unwrap();

        let d_soft = soft.max_deflection(Support::Cantilever, load).unwrap();
        let d_stiff = stiff.max_deflection(Support::Cantilever, load).unwrap();
        assert!(close(d_stiff, d_soft / 2.0, 1e-9));
        // Stiffness change must NOT move the bending moment (statics).
        let m_soft = soft.max_moment(Support::Cantilever, load).unwrap();
        let m_stiff = stiff.max_moment(Support::Cantilever, load).unwrap();
        assert!(close(m_soft, m_stiff, 1e-9));
    }

    #[test]
    fn deflection_scales_linearly_with_load() {
        // delta is linear in P / w (superposition).
        let section = Section::rectangular(20.0, 30.0).unwrap();
        let beam = Beam::new(1000.0, 200_000.0, section).unwrap();

        let d1 = beam
            .max_deflection(Support::Cantilever, Load::Point { force: 1000.0 })
            .unwrap();
        let d3 = beam
            .max_deflection(Support::Cantilever, Load::Point { force: 3000.0 })
            .unwrap();
        assert!(close(d3, 3.0 * d1, 1e-9));
    }

    // ---------------------------------------------------------------
    // Stress: sigma = M c / I, and the section-modulus equivalence.
    // ---------------------------------------------------------------

    #[test]
    fn stress_equals_moment_c_over_i() {
        let section = Section::rectangular(20.0, 30.0).unwrap();
        let beam = Beam::new(1000.0, 200_000.0, section).unwrap();
        let moment = 1.0e6; // N·mm
        let stress = beam.stress_from_moment(moment).unwrap();
        let expected = moment * section.extreme_fibre() / section.second_moment_area();
        assert!(close(stress, expected, 1e-9));
        // Equivalent via section modulus: sigma = M / S.
        assert!(close(
            stress,
            moment / section.section_modulus().unwrap(),
            1e-9
        ));
    }

    #[test]
    fn deeper_section_lowers_both_deflection_and_stress() {
        // A taller rectangle (more I, more S) deflects less and is less
        // stressed under the same load.
        let shallow =
            Beam::new(1000.0, 200_000.0, Section::rectangular(20.0, 30.0).unwrap()).unwrap();
        let deep = Beam::new(1000.0, 200_000.0, Section::rectangular(20.0, 60.0).unwrap()).unwrap();
        let load = Load::Point { force: 1000.0 };

        let r_shallow = shallow.analyze(Support::Cantilever, load).unwrap();
        let r_deep = deep.analyze(Support::Cantilever, load).unwrap();

        assert!(r_deep.max_deflection < r_shallow.max_deflection);
        assert!(r_deep.max_stress < r_shallow.max_stress);
        // Moment depends only on statics, so it is unchanged.
        assert!(close(r_deep.max_moment, r_shallow.max_moment, 1e-9));
    }

    #[test]
    fn flexural_rigidity_is_e_times_i() {
        let section = Section::rectangular(20.0, 30.0).unwrap();
        let beam = Beam::new(1000.0, 200_000.0, section).unwrap();
        assert!(close(beam.flexural_rigidity(), 200_000.0 * 45_000.0, 1e-3));
    }

    // ---------------------------------------------------------------
    // Relative ordering of the standard cases (known inequalities).
    // ---------------------------------------------------------------

    #[test]
    fn cantilever_point_is_softest_of_the_point_cases() {
        // For equal P, L, E, I the cantilever tip deflection
        // P L^3/(3 E I) is 16x the simply-supported centre deflection
        // P L^3/(48 E I).
        let section = Section::rectangular(20.0, 30.0).unwrap();
        let beam = Beam::new(1000.0, 200_000.0, section).unwrap();
        let load = Load::Point { force: 1000.0 };
        let cant = beam.max_deflection(Support::Cantilever, load).unwrap();
        let ss = beam.max_deflection(Support::SimplySupported, load).unwrap();
        assert!(close(cant, 16.0 * ss, 1e-6));
    }

    #[test]
    fn udl_case_ratio_is_known_constant() {
        // Cantilever UDL deflection w L^4/(8 E I) over simply-supported
        // UDL 5 w L^4/(384 E I) = (1/8)/(5/384) = 384/40 = 9.6.
        let section = Section::rectangular(20.0, 30.0).unwrap();
        let beam = Beam::new(1000.0, 200_000.0, section).unwrap();
        let load = Load::Udl { intensity: 2.0 };
        let cant = beam.max_deflection(Support::Cantilever, load).unwrap();
        let ss = beam.max_deflection(Support::SimplySupported, load).unwrap();
        assert!(close(cant / ss, 9.6, 1e-9));
    }

    // ---------------------------------------------------------------
    // Error paths.
    // ---------------------------------------------------------------

    #[test]
    fn beam_constructor_rejects_bad_input() {
        let section = Section::rectangular(20.0, 30.0).unwrap();
        assert!(Beam::new(0.0, 200_000.0, section).is_err());
        assert!(Beam::new(-1.0, 200_000.0, section).is_err());
        assert!(Beam::new(1000.0, 0.0, section).is_err());
        assert!(Beam::new(1000.0, f64::NAN, section).is_err());
    }

    #[test]
    fn analyze_rejects_nonpositive_load() {
        let section = Section::rectangular(20.0, 30.0).unwrap();
        let beam = Beam::new(1000.0, 200_000.0, section).unwrap();
        assert!(beam
            .analyze(Support::Cantilever, Load::Point { force: 0.0 })
            .is_err());
        assert!(beam
            .analyze(Support::Cantilever, Load::Point { force: -5.0 })
            .is_err());
        assert!(beam
            .analyze(Support::SimplySupported, Load::Udl { intensity: -1.0 })
            .is_err());
    }

    #[test]
    fn stress_rejects_negative_moment() {
        let section = Section::rectangular(20.0, 30.0).unwrap();
        let beam = Beam::new(1000.0, 200_000.0, section).unwrap();
        assert!(beam.stress_from_moment(-1.0).is_err());
        assert!(beam.stress_from_moment(f64::NAN).is_err());
        // Zero moment is a valid (unloaded) case -> zero stress.
        assert!(close(beam.stress_from_moment(0.0).unwrap(), 0.0, 1e-12));
    }

    #[test]
    fn error_metadata_is_stable() {
        let e = BeamError::bad_parameter("length", -1.0);
        assert_eq!(e.code(), "beam.bad-parameter");
        assert_eq!(e.category(), ErrorCategory::Input);

        let g = BeamError::DegenerateSection {
            reason: "second moment of area is zero",
        };
        assert_eq!(g.code(), "beam.degenerate-section");
        assert_eq!(g.category(), ErrorCategory::Geometry);
    }

    #[test]
    fn require_positive_guards_correctly() {
        assert!(close(
            BeamError::require_positive("x", 3.5).unwrap(),
            3.5,
            1e-12
        ));
        assert!(BeamError::require_positive("x", 0.0).is_err());
        assert!(BeamError::require_positive("x", -2.0).is_err());
        assert!(BeamError::require_positive("x", f64::INFINITY).is_err());
        assert!(BeamError::require_positive("x", f64::NAN).is_err());
    }

    // ---------------------------------------------------------------
    // Serialization round-trips for the public structs.
    // ---------------------------------------------------------------

    #[test]
    fn result_round_trips_through_json() {
        let section = Section::circular(10.0).unwrap();
        let beam = Beam::new(800.0, 70_000.0, section).unwrap();
        let r = beam
            .analyze(Support::SimplySupported, Load::Point { force: 500.0 })
            .unwrap();
        let json = serde_json::to_string(&r).unwrap();
        let back: BeamResult = serde_json::from_str(&json).unwrap();
        assert!(close(back.max_deflection, r.max_deflection, 1e-12));
        assert!(close(back.max_moment, r.max_moment, 1e-12));
        assert!(close(back.max_stress, r.max_stress, 1e-12));
    }

    #[test]
    fn beam_and_section_round_trip_through_json() {
        let section = Section::rectangular(20.0, 30.0).unwrap();
        let beam = Beam::new(1000.0, 200_000.0, section).unwrap();
        let json = serde_json::to_string(&beam).unwrap();
        let back: Beam = serde_json::from_str(&json).unwrap();
        assert_eq!(beam, back);
    }
}
