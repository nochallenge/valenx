//! # valenx-plate
//!
//! Closed-form small-deflection bending of thin flat plates.
//!
//! ## What
//!
//! This crate evaluates the classical textbook formulas for the transverse
//! bending of a thin flat plate of uniform thickness made of a linear-
//! elastic isotropic material:
//!
//! - the **flexural rigidity** `D = E t^3 / (12 (1 - nu^2))`
//!   ([`PlateMaterial::flexural_rigidity`]);
//! - the **centre deflection** `w = k p a^4 / D` of a uniformly-loaded
//!   circular plate ([`CircularPlate::center_deflection`]);
//! - the **maximum bending stress** in that plate
//!   ([`CircularPlate::max_bending_stress`]);
//!
//! for both clamped and simply-supported rims ([`EdgeSupport`]).
//!
//! ## Model
//!
//! Kirchhoff-Love thin-plate theory in the small-deflection (linear)
//! regime. The plate is flat, of uniform thickness, isotropic and linear-
//! elastic; deflections are assumed small compared with the thickness, so
//! the governing biharmonic equation `D nabla^4 w = p` is linear and the
//! centre deflection is proportional to the load. The dimensionless
//! coefficient `k` and the peak-moment expressions are the standard results
//! of Timoshenko & Woinowsky-Krieger, *Theory of Plates and Shells*
//! (2nd ed.):
//!
//! ```text
//! D       = E t^3 / (12 (1 - nu^2))
//!
//! clamped:           w = (1/64) p a^4 / D
//!                    sigma_max = 0.75 p (a/t)^2          (radial, at the edge)
//! simply-supported:  w = (5 + nu) / (64 (1 + nu)) p a^4 / D
//!                    sigma_max = 3 (3 + nu) / 8 p (a/t)^2   (at the centre)
//! ```
//!
//! Because `D` is proportional to `t^3`, the deflection scales as `1/t^3`
//! and as `a^4`; both stresses scale linearly with `p` and quadratically
//! with the radius-to-thickness ratio `a / t`. For any admissible Poisson
//! ratio the simply-supported deflection coefficient exceeds `1/64`, so a
//! clamped plate is always the stiffer of the two.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form /
//! numerical models, NOT a clinical / medical / production engineering
//! tool. The crate covers only the small-deflection, linear-elastic,
//! isotropic, uniformly-loaded *circular* case: there is no large-
//! deflection membrane stiffening, no plasticity, no orthotropy, no
//! transverse-shear (Mindlin) correction, no rectangular / annular / point-
//! loaded plates, and no buckling or dynamics. Do not size a real
//! structure from these numbers without an independent, validated analysis.
//!
//! ## Example
//!
//! ```
//! use valenx_plate::{CircularPlate, EdgeSupport, PlateMaterial};
//!
//! // Steel disc: E = 200 GPa, nu = 0.3, t = 5 mm, radius 250 mm,
//! // uniform pressure 20 kPa, clamped rim. SI units (Pa, m).
//! let steel = PlateMaterial::new(200.0e9, 0.3, 0.005).unwrap();
//! let plate = CircularPlate::new(steel, 0.25, 20.0e3, EdgeSupport::Clamped).unwrap();
//!
//! let w = plate.center_deflection();      // metres
//! let s = plate.max_bending_stress();     // pascals
//! assert!(w > 0.0 && s > 0.0);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod circular;
pub mod error;
pub mod material;

pub use circular::{CircularPlate, EdgeSupport};
pub use error::{ErrorCategory, PlateError};
pub use material::PlateMaterial;

#[cfg(test)]
mod tests {
    use super::*;

    /// Relative tolerance for comparing two finite floats.
    fn approx_rel(a: f64, b: f64, rel: f64) -> bool {
        let diff = (a - b).abs();
        let scale = a.abs().max(b.abs()).max(1.0);
        diff < rel * scale
    }

    // -- flexural rigidity ------------------------------------------------

    #[test]
    fn flexural_rigidity_matches_closed_form() {
        let e = 200.0e9;
        let nu = 0.3;
        let t = 0.01;
        let m = PlateMaterial::new(e, nu, t).unwrap();

        let expected = e * t.powi(3) / (12.0 * (1.0 - nu * nu));
        assert!(
            approx_rel(m.flexural_rigidity(), expected, 1e-12),
            "D = {got}, expected {expected}",
            got = m.flexural_rigidity()
        );
    }

    #[test]
    fn flexural_rigidity_unit_inputs() {
        // E = 1, nu = 0, t = 1 -> D = 1/12 exactly.
        let m = PlateMaterial::new(1.0, 0.0, 1.0).unwrap();
        let expected = 1.0 / 12.0;
        let diff = (m.flexural_rigidity() - expected).abs();
        assert!(diff < 1e-15, "D = {d}", d = m.flexural_rigidity());
    }

    #[test]
    fn flexural_rigidity_scales_thickness_cubed() {
        // Doubling thickness multiplies D by 2^3 = 8.
        let base = PlateMaterial::new(70.0e9, 0.33, 0.004).unwrap();
        let thick = PlateMaterial::new(70.0e9, 0.33, 0.008).unwrap();
        let ratio = thick.flexural_rigidity() / base.flexural_rigidity();
        assert!(approx_rel(ratio, 8.0, 1e-12), "ratio = {ratio}");
    }

    #[test]
    fn flexural_rigidity_poisson_factor() {
        // D(nu) / D(0) = 1 / (1 - nu^2).
        let nu = 0.3;
        let d0 = PlateMaterial::new(1.0, 0.0, 1.0)
            .unwrap()
            .flexural_rigidity();
        let dn = PlateMaterial::new(1.0, nu, 1.0)
            .unwrap()
            .flexural_rigidity();
        let expected = 1.0 / (1.0 - nu * nu);
        assert!(
            approx_rel(dn / d0, expected, 1e-12),
            "ratio = {r}",
            r = dn / d0
        );
    }

    // -- circular plate centre deflection ---------------------------------

    #[test]
    fn clamped_center_deflection_closed_form() {
        // w_clamped = p a^4 / (64 D).
        let m = PlateMaterial::new(200.0e9, 0.3, 0.005).unwrap();
        let a = 0.25;
        let p = 20.0e3;
        let plate = CircularPlate::new(m, a, p, EdgeSupport::Clamped).unwrap();

        let expected = p * a.powi(4) / (64.0 * m.flexural_rigidity());
        assert!(
            approx_rel(plate.center_deflection(), expected, 1e-12),
            "w = {w}, expected {expected}",
            w = plate.center_deflection()
        );
    }

    #[test]
    fn simply_supported_center_deflection_closed_form() {
        // w_ss = (5 + nu)/(1 + nu) * p a^4 / (64 D).
        let nu = 0.3;
        let m = PlateMaterial::new(200.0e9, nu, 0.005).unwrap();
        let a = 0.25;
        let p = 20.0e3;
        let plate = CircularPlate::new(m, a, p, EdgeSupport::SimplySupported).unwrap();

        let expected = (5.0 + nu) / (1.0 + nu) * p * a.powi(4) / (64.0 * m.flexural_rigidity());
        assert!(
            approx_rel(plate.center_deflection(), expected, 1e-12),
            "w = {w}, expected {expected}",
            w = plate.center_deflection()
        );
    }

    #[test]
    fn clamped_deflection_coefficient_is_one_over_64() {
        let plate = sample_plate(EdgeSupport::Clamped);
        let diff = (plate.deflection_coefficient() - 1.0 / 64.0).abs();
        assert!(diff < 1e-15, "k = {k}", k = plate.deflection_coefficient());
    }

    #[test]
    fn deflection_scales_radius_to_the_fourth() {
        // Same plate, radius x2 -> deflection x16 (k and D unchanged).
        let m = PlateMaterial::new(70.0e9, 0.33, 0.003).unwrap();
        let p = 5.0e3;
        let small = CircularPlate::new(m, 0.1, p, EdgeSupport::Clamped).unwrap();
        let big = CircularPlate::new(m, 0.2, p, EdgeSupport::Clamped).unwrap();
        let ratio = big.center_deflection() / small.center_deflection();
        assert!(approx_rel(ratio, 16.0, 1e-12), "ratio = {ratio}");
    }

    #[test]
    fn deflection_scales_inverse_thickness_cubed() {
        // Halving thickness multiplies deflection by 8 (D ~ t^3), with
        // radius/thickness kept admissible in both cases.
        let thick = PlateMaterial::new(70.0e9, 0.33, 0.010).unwrap();
        let thin = PlateMaterial::new(70.0e9, 0.33, 0.005).unwrap();
        let p = 5.0e3;
        let a = 0.20;
        let plate_thick = CircularPlate::new(thick, a, p, EdgeSupport::Clamped).unwrap();
        let plate_thin = CircularPlate::new(thin, a, p, EdgeSupport::Clamped).unwrap();
        let ratio = plate_thin.center_deflection() / plate_thick.center_deflection();
        assert!(approx_rel(ratio, 8.0, 1e-12), "ratio = {ratio}");
    }

    #[test]
    fn deflection_scales_linearly_with_pressure() {
        let m = PlateMaterial::new(70.0e9, 0.33, 0.004).unwrap();
        let a = 0.15;
        let low = CircularPlate::new(m, a, 1.0e3, EdgeSupport::Clamped).unwrap();
        let high = CircularPlate::new(m, a, 3.0e3, EdgeSupport::Clamped).unwrap();
        let ratio = high.center_deflection() / low.center_deflection();
        assert!(approx_rel(ratio, 3.0, 1e-12), "ratio = {ratio}");
    }

    #[test]
    fn clamped_is_stiffer_than_simply_supported() {
        // For identical inputs the clamped plate deflects less.
        let m = PlateMaterial::new(200.0e9, 0.3, 0.005).unwrap();
        let a = 0.25;
        let p = 20.0e3;
        let clamped = CircularPlate::new(m, a, p, EdgeSupport::Clamped).unwrap();
        let ss = CircularPlate::new(m, a, p, EdgeSupport::SimplySupported).unwrap();
        assert!(
            clamped.center_deflection() < ss.center_deflection(),
            "clamped {c} should be < simply-supported {s}",
            c = clamped.center_deflection(),
            s = ss.center_deflection()
        );

        // The ratio is exactly (1 + nu) / (5 + nu) for nu = 0.3 -> 1.3/5.3.
        let nu = 0.3;
        let expected = (1.0 + nu) / (5.0 + nu);
        let ratio = clamped.center_deflection() / ss.center_deflection();
        assert!(approx_rel(ratio, expected, 1e-12), "ratio = {ratio}");
    }

    // -- maximum bending stress -------------------------------------------

    #[test]
    fn clamped_max_stress_closed_form() {
        // sigma_max = 3 p a^2 / (4 t^2) = 0.75 p (a/t)^2.
        let t = 0.005;
        let m = PlateMaterial::new(200.0e9, 0.3, t).unwrap();
        let a = 0.25;
        let p = 20.0e3;
        let plate = CircularPlate::new(m, a, p, EdgeSupport::Clamped).unwrap();

        let expected = 0.75 * p * (a / t).powi(2);
        assert!(
            approx_rel(plate.max_bending_stress(), expected, 1e-12),
            "sigma = {s}, expected {expected}",
            s = plate.max_bending_stress()
        );
    }

    #[test]
    fn simply_supported_max_stress_closed_form() {
        // sigma_max = 3 (3 + nu) p a^2 / (8 t^2).
        let nu = 0.3;
        let t = 0.005;
        let m = PlateMaterial::new(200.0e9, nu, t).unwrap();
        let a = 0.25;
        let p = 20.0e3;
        let plate = CircularPlate::new(m, a, p, EdgeSupport::SimplySupported).unwrap();

        let expected = 3.0 * (3.0 + nu) / 8.0 * p * (a / t).powi(2);
        assert!(
            approx_rel(plate.max_bending_stress(), expected, 1e-12),
            "sigma = {s}, expected {expected}",
            s = plate.max_bending_stress()
        );
    }

    #[test]
    fn stress_scales_linearly_with_pressure() {
        let m = PlateMaterial::new(70.0e9, 0.33, 0.004).unwrap();
        let a = 0.15;
        let low = CircularPlate::new(m, a, 2.0e3, EdgeSupport::Clamped).unwrap();
        let high = CircularPlate::new(m, a, 10.0e3, EdgeSupport::Clamped).unwrap();
        let ratio = high.max_bending_stress() / low.max_bending_stress();
        assert!(approx_rel(ratio, 5.0, 1e-12), "ratio = {ratio}");
    }

    #[test]
    fn stress_scales_radius_thickness_ratio_squared() {
        // Doubling a/t (here via radius) quadruples the stress.
        let m = PlateMaterial::new(70.0e9, 0.33, 0.004).unwrap();
        let p = 5.0e3;
        let small = CircularPlate::new(m, 0.10, p, EdgeSupport::Clamped).unwrap();
        let big = CircularPlate::new(m, 0.20, p, EdgeSupport::Clamped).unwrap();
        let ratio = big.max_bending_stress() / small.max_bending_stress();
        assert!(approx_rel(ratio, 4.0, 1e-12), "ratio = {ratio}");

        // And via thickness: halving t (doubling a/t) also quadruples stress.
        let thick = PlateMaterial::new(70.0e9, 0.33, 0.008).unwrap();
        let thin = PlateMaterial::new(70.0e9, 0.33, 0.004).unwrap();
        let a = 0.20;
        let s_thick = CircularPlate::new(thick, a, p, EdgeSupport::Clamped)
            .unwrap()
            .max_bending_stress();
        let s_thin = CircularPlate::new(thin, a, p, EdgeSupport::Clamped)
            .unwrap()
            .max_bending_stress();
        assert!(
            approx_rel(s_thin / s_thick, 4.0, 1e-12),
            "ratio = {r}",
            r = s_thin / s_thick
        );
    }

    #[test]
    fn stress_is_independent_of_modulus() {
        // The thin-plate bending stress does not depend on E.
        let soft = PlateMaterial::new(10.0e9, 0.3, 0.005).unwrap();
        let stiff = PlateMaterial::new(400.0e9, 0.3, 0.005).unwrap();
        let a = 0.25;
        let p = 20.0e3;
        let s_soft = CircularPlate::new(soft, a, p, EdgeSupport::Clamped)
            .unwrap()
            .max_bending_stress();
        let s_stiff = CircularPlate::new(stiff, a, p, EdgeSupport::Clamped)
            .unwrap()
            .max_bending_stress();
        assert!(approx_rel(s_soft, s_stiff, 1e-12), "{s_soft} vs {s_stiff}");
    }

    // -- validation / errors ----------------------------------------------

    #[test]
    fn rejects_nonpositive_thickness() {
        let err = PlateMaterial::new(200.0e9, 0.3, 0.0).unwrap_err();
        assert_eq!(err.code(), "plate.invalid_parameter");
        assert_eq!(err.category(), ErrorCategory::Input);
        assert!(matches!(
            err,
            PlateError::InvalidParameter {
                name: "thickness",
                ..
            }
        ));
    }

    #[test]
    fn rejects_nonfinite_modulus() {
        let err = PlateMaterial::new(f64::NAN, 0.3, 0.01).unwrap_err();
        assert!(matches!(
            err,
            PlateError::InvalidParameter {
                name: "youngs_modulus",
                ..
            }
        ));
        let err2 = PlateMaterial::new(f64::INFINITY, 0.3, 0.01).unwrap_err();
        assert!(matches!(
            err2,
            PlateError::InvalidParameter {
                name: "youngs_modulus",
                ..
            }
        ));
    }

    #[test]
    fn rejects_poisson_at_or_beyond_bounds() {
        // nu = 0.5 (incompressible upper bound) is rejected: 1 - nu^2 would
        // not blow up, but the open-interval contract excludes it.
        assert!(PlateMaterial::new(1.0, 0.5, 1.0).is_err());
        // nu = 1.0 makes 1 - nu^2 = 0 (singular D).
        assert!(PlateMaterial::new(1.0, 1.0, 1.0).is_err());
        // nu = -1.0 (lower thermodynamic bound) is rejected.
        assert!(PlateMaterial::new(1.0, -1.0, 1.0).is_err());
        // A value just inside the interval is accepted.
        assert!(PlateMaterial::new(1.0, 0.49, 1.0).is_ok());
        assert!(PlateMaterial::new(1.0, -0.99, 1.0).is_ok());
    }

    #[test]
    fn rejects_nonpositive_pressure_and_radius() {
        let m = PlateMaterial::new(200.0e9, 0.3, 0.005).unwrap();
        let bad_p = CircularPlate::new(m, 0.25, 0.0, EdgeSupport::Clamped).unwrap_err();
        assert!(matches!(
            bad_p,
            PlateError::InvalidParameter {
                name: "pressure",
                ..
            }
        ));
        let bad_a = CircularPlate::new(m, -1.0, 1.0e3, EdgeSupport::Clamped).unwrap_err();
        assert!(matches!(
            bad_a,
            PlateError::InvalidParameter { name: "radius", .. }
        ));
    }

    #[test]
    fn rejects_thick_plate() {
        // radius / thickness = 0.02 / 0.01 = 2 < 5 -> NotThin.
        let m = PlateMaterial::new(200.0e9, 0.3, 0.01).unwrap();
        let err = CircularPlate::new(m, 0.02, 1.0e3, EdgeSupport::Clamped).unwrap_err();
        assert_eq!(err.code(), "plate.not_thin");
        assert_eq!(err.category(), ErrorCategory::Algorithm);
        match err {
            PlateError::NotThin { ratio, min } => {
                assert!(approx_rel(ratio, 2.0, 1e-12), "ratio = {ratio}");
                let diff = (min - 5.0).abs();
                assert!(diff < 1e-15, "min = {min}");
            }
            other => panic!("expected NotThin, got {other:?}"),
        }
    }

    #[test]
    fn thin_boundary_is_accepted() {
        // Exactly at the threshold a/t = 5 is admissible.
        let m = PlateMaterial::new(200.0e9, 0.3, 0.01).unwrap();
        assert!(CircularPlate::new(m, 0.05, 1.0e3, EdgeSupport::Clamped).is_ok());
    }

    // -- serde round-trip -------------------------------------------------

    #[test]
    fn material_serde_round_trip() {
        let m = PlateMaterial::new(200.0e9, 0.3, 0.005).unwrap();
        let json = serde_json::to_string(&m).unwrap();
        let back: PlateMaterial = serde_json::from_str(&json).unwrap();
        assert_eq!(m, back);
    }

    #[test]
    fn plate_serde_round_trip() {
        let plate = sample_plate(EdgeSupport::SimplySupported);
        let json = serde_json::to_string(&plate).unwrap();
        let back: CircularPlate = serde_json::from_str(&json).unwrap();
        assert_eq!(plate, back);
        let diff = (plate.center_deflection() - back.center_deflection()).abs();
        assert!(diff < 1e-15);
    }

    /// A representative steel-disc plate for the structural-property tests.
    fn sample_plate(support: EdgeSupport) -> CircularPlate {
        let m = PlateMaterial::new(200.0e9, 0.3, 0.005).unwrap();
        CircularPlate::new(m, 0.25, 20.0e3, support).unwrap()
    }
}
