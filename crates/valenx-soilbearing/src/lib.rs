//! # valenx-soilbearing
//!
//! Shallow-foundation bearing capacity by the classical Terzaghi
//! general-shear equation.
//!
//! ## What
//!
//! Given soil strength parameters and a strip-footing geometry, this
//! crate computes:
//!
//! - the dimensionless bearing-capacity factors `Nc`, `Nq`, `Ngamma`
//!   ([`BearingFactors::from_friction_angle`]);
//! - the ultimate bearing capacity `qult`
//!   ([`ultimate_bearing_capacity`]); and
//! - the allowable bearing pressure `qall = qult / FS` together with a
//!   term-by-term breakdown ([`bearing_capacity`] returning a
//!   [`BearingResult`]).
//!
//! Inputs are wrapped in validated value types ([`SoilProperties`],
//! [`Footing`]) whose constructors reject out-of-domain numbers, so the
//! solver itself is total.
//!
//! ## Model
//!
//! The Terzaghi (1943) ultimate bearing capacity of a continuous strip
//! footing under general shear failure is the sum of three terms:
//!
//! `qult = c * Nc + q * Nq + 0.5 * gamma * B * Ngamma`
//!
//! where:
//!
//! - `c` is the soil cohesion and `Nc` the cohesion factor;
//! - `q = gamma * Df` is the surcharge from soil of unit weight `gamma`
//!   over the founding depth `Df`, weighted by the surcharge factor
//!   `Nq`;
//! - `B` is the footing width and `Ngamma` the self-weight factor.
//!
//! The bearing-capacity factors are closed-form functions of the
//! drained friction angle `phi` (radians):
//!
//! `Nq = exp(pi * tan(phi)) * tan^2(pi/4 + phi/2)`
//!
//! `Nc = (Nq - 1) * cot(phi)`  (Prandtl/Reissner form; `Nc -> pi + 2`
//! as `phi -> 0`)
//!
//! `Ngamma = 2 * (Nq + 1) * tan(phi)`  (Vesic 1973 form)
//!
//! The allowable pressure applies a global factor of safety:
//! `qall = qult / FS`.
//!
//! Units are the caller's responsibility but must be self-consistent. A
//! convenient SI system is metres for `B` and `Df`, kN/m^3 for `gamma`,
//! and kPa for `c`, which yields `qult` and `qall` in kPa.
//!
//! ## Honest scope
//!
//! Research/educational grade. This crate implements textbook
//! closed-form geotechnical models (Terzaghi 1943; Prandtl/Reissner
//! `Nc`; Vesic 1973 `Ngamma`) for a single, homogeneous, drained
//! stratum under a vertically and concentrically loaded continuous
//! strip footing on level ground. It deliberately omits effects that a
//! real design must consider, including but not limited to: shape,
//! depth, and load-inclination factors (Meyerhof/Hansen/Vesic
//! corrections); the groundwater table and effective-stress reduction
//! of `gamma`; local- versus general-shear failure transition;
//! eccentric loading and effective-width reduction; layered or
//! anisotropic soils; settlement (this is a strength check only, not a
//! serviceability check); and seismic or dynamic loading. The choice of
//! `Ngamma` expression alone changes results materially between
//! published correlations. It is NOT a clinical/medical/production
//! engineering tool and is not a substitute for a licensed
//! geotechnical engineer, a site-specific investigation, laboratory
//! testing, or the governing building code. Do not use it for design,
//! construction, or any decision affecting safety.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod capacity;
pub mod error;
pub mod factors;
pub mod footing;
pub mod soil;

pub use capacity::{
    allowable_from_ultimate, bearing_capacity, ultimate_bearing_capacity, BearingResult,
};
pub use error::SoilBearingError;
pub use factors::BearingFactors;
pub use footing::Footing;
pub use soil::SoilProperties;

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Absolute tolerance for analytic float comparisons.
    const EPS: f64 = 1e-9;

    // --- Bearing-capacity factors: known anchor values -----------------

    #[test]
    fn factors_at_phi_zero_are_prandtl_limits() {
        let soil = SoilProperties::new(0.0, 10.0, 18.0).unwrap();
        let f = BearingFactors::from_friction_angle(&soil);
        assert!(
            (f.nq - 1.0).abs() < EPS,
            "Nq(0) should be 1, got {nq}",
            nq = f.nq
        );
        assert!(
            f.ngamma.abs() < EPS,
            "Ngamma(0) should be 0, got {ng}",
            ng = f.ngamma
        );
        assert!(
            (f.nc - (PI + 2.0)).abs() < EPS,
            "Nc(0) should be pi+2 = {expected}, got {nc}",
            expected = PI + 2.0,
            nc = f.nc
        );
    }

    #[test]
    fn factors_at_phi_30_match_textbook() {
        // Standard tabulated values (Das, Principles of Foundation
        // Engineering) for phi = 30 deg: Nq = 18.401, Nc = 30.140,
        // and Vesic Ngamma = 22.402. Check to 2 decimal places.
        let soil = SoilProperties::new(30.0, 0.0, 18.0).unwrap();
        let f = BearingFactors::from_friction_angle(&soil);
        assert!((f.nq - 18.401).abs() < 1e-2, "Nq(30) = {nq}", nq = f.nq);
        assert!((f.nc - 30.140).abs() < 1e-2, "Nc(30) = {nc}", nc = f.nc);
        assert!(
            (f.ngamma - 22.402).abs() < 1e-2,
            "Ngamma(30) = {ng}",
            ng = f.ngamma
        );
    }

    #[test]
    fn nq_closed_form_independent_recompute() {
        // Re-derive Nq from the closed form independently and compare.
        for phi_deg in [5.0_f64, 12.5, 20.0, 35.0, 45.0] {
            let soil = SoilProperties::new(phi_deg, 0.0, 18.0).unwrap();
            let f = BearingFactors::from_friction_angle(&soil);
            let phi = phi_deg.to_radians();
            let expected = (PI * phi.tan()).exp() * (PI / 4.0 + phi / 2.0).tan().powi(2);
            assert!(
                (f.nq - expected).abs() < EPS,
                "Nq mismatch at {phi_deg} deg: got {got}, expected {expected}",
                got = f.nq
            );
        }
    }

    #[test]
    fn nc_satisfies_defining_identity() {
        // Nc = (Nq - 1) cot(phi) must hold exactly for phi > 0.
        for phi_deg in [1.0_f64, 10.0, 25.0, 40.0] {
            let soil = SoilProperties::new(phi_deg, 0.0, 18.0).unwrap();
            let f = BearingFactors::from_friction_angle(&soil);
            let phi = phi_deg.to_radians();
            let from_identity = (f.nq - 1.0) / phi.tan();
            assert!(
                (f.nc - from_identity).abs() < EPS,
                "Nc identity broken at {phi_deg} deg"
            );
        }
    }

    // --- Monotonicity: factors rise with friction angle ----------------

    #[test]
    fn all_factors_increase_with_friction_angle() {
        let mut prev: Option<BearingFactors> = None;
        for phi_deg in [0.0_f64, 5.0, 10.0, 20.0, 30.0, 40.0, 45.0] {
            let soil = SoilProperties::new(phi_deg, 0.0, 18.0).unwrap();
            let f = BearingFactors::from_friction_angle(&soil);
            if let Some(p) = prev {
                assert!(f.nq > p.nq, "Nq not increasing at {phi_deg} deg");
                assert!(f.nc > p.nc, "Nc not increasing at {phi_deg} deg");
                // Ngamma is 0 at phi = 0 and strictly positive after, so
                // it is non-decreasing across the first step and strictly
                // increasing thereafter.
                assert!(f.ngamma >= p.ngamma, "Ngamma decreased at {phi_deg} deg");
            }
            prev = Some(f);
        }
    }

    // --- Ultimate-capacity equation ------------------------------------

    #[test]
    fn qult_equals_sum_of_three_terms() {
        let soil = SoilProperties::new(25.0, 20.0, 19.0).unwrap();
        let footing = Footing::new(2.5, 1.2).unwrap();
        let f = BearingFactors::from_friction_angle(&soil);
        let q = soil.unit_weight() * footing.depth();
        let expected = soil.cohesion() * f.nc
            + q * f.nq
            + 0.5 * soil.unit_weight() * footing.width() * f.ngamma;
        let got = ultimate_bearing_capacity(&soil, &footing);
        assert!(
            (got - expected).abs() < EPS,
            "qult mismatch: got {got}, expected {expected}"
        );
    }

    #[test]
    fn result_terms_sum_to_qult() {
        let soil = SoilProperties::new(28.0, 15.0, 18.0).unwrap();
        let footing = Footing::new(1.8, 0.9).unwrap();
        let r = bearing_capacity(&soil, &footing, 2.5).unwrap();
        let summed = r.cohesion_term + r.surcharge_term + r.self_weight_term;
        assert!(
            (summed - r.q_ultimate).abs() < EPS,
            "terms {summed} should sum to qult {qult}",
            qult = r.q_ultimate
        );
    }

    #[test]
    fn phi_zero_clay_reduces_to_c_nc_plus_surcharge() {
        // Undrained clay: qult = c*(pi+2) + gamma*Df  (Nq=1, Ngamma=0).
        let c = 60.0;
        let gamma = 18.0;
        let df = 1.5;
        let soil = SoilProperties::new(0.0, c, gamma).unwrap();
        let footing = Footing::new(3.0, df).unwrap();
        let expected = c * (PI + 2.0) + gamma * df;
        let got = ultimate_bearing_capacity(&soil, &footing);
        assert!(
            (got - expected).abs() < 1e-7,
            "phi=0 qult: got {got}, expected {expected}"
        );
    }

    // --- Allowable = qult / FS -----------------------------------------

    #[test]
    fn allowable_is_ultimate_over_fs() {
        let soil = SoilProperties::new(32.0, 5.0, 18.5).unwrap();
        let footing = Footing::new(2.0, 1.0).unwrap();
        let fs = 3.0;
        let r = bearing_capacity(&soil, &footing, fs).unwrap();
        assert!(
            (r.q_allowable - r.q_ultimate / fs).abs() < EPS,
            "qall {qall} != qult/FS {ratio}",
            qall = r.q_allowable,
            ratio = r.q_ultimate / fs
        );
        assert!((r.factor_of_safety - fs).abs() < EPS);
    }

    #[test]
    fn allowable_from_ultimate_divides_correctly() {
        let qall = allowable_from_ultimate(450.0, 3.0).unwrap();
        assert!((qall - 150.0).abs() < EPS, "got {qall}");
    }

    #[test]
    fn larger_fs_gives_smaller_allowable() {
        let soil = SoilProperties::new(30.0, 10.0, 18.0).unwrap();
        let footing = Footing::new(2.0, 1.0).unwrap();
        let low = bearing_capacity(&soil, &footing, 2.0).unwrap();
        let high = bearing_capacity(&soil, &footing, 4.0).unwrap();
        // Same qult, larger divisor -> smaller allowable.
        assert!((low.q_ultimate - high.q_ultimate).abs() < EPS);
        assert!(high.q_allowable < low.q_allowable);
    }

    // --- Cohesionless soil drops the c term ----------------------------

    #[test]
    fn cohesionless_soil_drops_cohesion_term() {
        let footing = Footing::new(2.0, 1.0).unwrap();
        let sand = SoilProperties::new(34.0, 0.0, 18.5).unwrap();
        assert!(sand.is_cohesionless());
        let r = bearing_capacity(&sand, &footing, 3.0).unwrap();
        assert!(r.cohesion_term.abs() < EPS, "cohesion term should vanish");
        // qult is then exactly surcharge + self-weight.
        assert!((r.q_ultimate - (r.surcharge_term + r.self_weight_term)).abs() < EPS);
    }

    #[test]
    fn adding_cohesion_only_changes_cohesion_term() {
        // Holding phi, gamma, B, Df fixed, switching c from 0 to c0
        // raises qult by exactly c0 * Nc and leaves the other two terms
        // untouched.
        let footing = Footing::new(2.2, 1.1).unwrap();
        let phi = 26.0;
        let gamma = 18.0;
        let dry = SoilProperties::new(phi, 0.0, gamma).unwrap();
        let cohesive = SoilProperties::new(phi, 25.0, gamma).unwrap();
        let r0 = bearing_capacity(&dry, &footing, 3.0).unwrap();
        let r1 = bearing_capacity(&cohesive, &footing, 3.0).unwrap();
        assert!((r0.surcharge_term - r1.surcharge_term).abs() < EPS);
        assert!((r0.self_weight_term - r1.self_weight_term).abs() < EPS);
        let delta = r1.q_ultimate - r0.q_ultimate;
        assert!(
            (delta - 25.0 * r1.factors.nc).abs() < EPS,
            "delta qult {delta} should equal c*Nc"
        );
    }

    // --- Deeper footing (q up) raises capacity -------------------------

    #[test]
    fn deeper_footing_raises_capacity() {
        let soil = SoilProperties::new(30.0, 5.0, 18.0).unwrap();
        let shallow = Footing::new(2.0, 0.5).unwrap();
        let deep = Footing::new(2.0, 2.0).unwrap();
        let q_shallow = ultimate_bearing_capacity(&soil, &shallow);
        let q_deep = ultimate_bearing_capacity(&soil, &deep);
        assert!(
            q_deep > q_shallow,
            "deeper footing should raise qult: deep {q_deep} <= shallow {q_shallow}"
        );
    }

    #[test]
    fn surface_footing_has_no_surcharge_term() {
        let soil = SoilProperties::new(30.0, 5.0, 18.0).unwrap();
        let surface = Footing::new(2.0, 0.0).unwrap();
        let r = bearing_capacity(&soil, &surface, 3.0).unwrap();
        assert!(
            r.surcharge_term.abs() < EPS,
            "Df=0 -> q=0 -> surcharge term 0"
        );
    }

    #[test]
    fn deeper_footing_raises_only_via_surcharge_term() {
        // Increasing Df changes only q = gamma*Df, hence only the
        // surcharge term; the increment equals gamma*dDf*Nq.
        let soil = SoilProperties::new(28.0, 10.0, 18.0).unwrap();
        let f0 = Footing::new(2.0, 1.0).unwrap();
        let f1 = Footing::new(2.0, 2.5).unwrap();
        let r0 = bearing_capacity(&soil, &f0, 3.0).unwrap();
        let r1 = bearing_capacity(&soil, &f1, 3.0).unwrap();
        assert!((r0.cohesion_term - r1.cohesion_term).abs() < EPS);
        assert!((r0.self_weight_term - r1.self_weight_term).abs() < EPS);
        let increment = r1.surcharge_term - r0.surcharge_term;
        let expected = soil.unit_weight() * (2.5 - 1.0) * r0.factors.nq;
        assert!(
            (increment - expected).abs() < EPS,
            "surcharge increment {increment} != gamma*dDf*Nq {expected}"
        );
    }

    // --- Wider footing raises the self-weight term ---------------------

    #[test]
    fn wider_footing_raises_self_weight_term() {
        // For a frictional soil (Ngamma > 0), a wider B raises qult only
        // through the 0.5*gamma*B*Ngamma term.
        let soil = SoilProperties::new(33.0, 0.0, 18.5).unwrap();
        let narrow = Footing::new(1.0, 1.0).unwrap();
        let wide = Footing::new(3.0, 1.0).unwrap();
        let qn = ultimate_bearing_capacity(&soil, &narrow);
        let qw = ultimate_bearing_capacity(&soil, &wide);
        assert!(qw > qn, "wider footing should raise qult");
    }

    // --- Validation / error paths --------------------------------------

    #[test]
    fn rejects_friction_angle_at_or_above_90() {
        assert!(SoilProperties::new(90.0, 0.0, 18.0).is_err());
        assert!(SoilProperties::new(95.0, 0.0, 18.0).is_err());
        // Just under 90 is accepted.
        assert!(SoilProperties::new(89.999, 0.0, 18.0).is_ok());
    }

    #[test]
    fn rejects_negative_inputs() {
        assert!(SoilProperties::new(-1.0, 0.0, 18.0).is_err());
        assert!(SoilProperties::new(30.0, -5.0, 18.0).is_err());
        assert!(SoilProperties::new(30.0, 0.0, -1.0).is_err());
        assert!(Footing::new(-1.0, 1.0).is_err());
        assert!(Footing::new(0.0, 1.0).is_err());
        assert!(Footing::new(2.0, -0.1).is_err());
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert!(SoilProperties::new(f64::NAN, 0.0, 18.0).is_err());
        assert!(SoilProperties::new(30.0, f64::INFINITY, 18.0).is_err());
        assert!(Footing::new(f64::NAN, 1.0).is_err());
    }

    #[test]
    fn rejects_factor_of_safety_at_or_below_one() {
        let soil = SoilProperties::new(30.0, 5.0, 18.0).unwrap();
        let footing = Footing::new(2.0, 1.0).unwrap();
        assert!(bearing_capacity(&soil, &footing, 1.0).is_err());
        assert!(bearing_capacity(&soil, &footing, 0.5).is_err());
        assert!(bearing_capacity(&soil, &footing, f64::NAN).is_err());
        assert!(bearing_capacity(&soil, &footing, 1.0001).is_ok());
        assert!(allowable_from_ultimate(300.0, 1.0).is_err());
        assert!(allowable_from_ultimate(f64::INFINITY, 3.0).is_err());
    }

    #[test]
    fn error_codes_are_stable() {
        let err = SoilProperties::new(120.0, 0.0, 18.0).unwrap_err();
        assert_eq!(err.code(), "soilbearing.invalid_parameter");
        let nf = Footing::new(f64::NAN, 1.0).unwrap_err();
        assert_eq!(nf.code(), "soilbearing.not_finite");
    }

    // --- Serde round-trip ----------------------------------------------

    #[test]
    fn soil_properties_serde_round_trip() {
        let soil = SoilProperties::new(31.5, 12.0, 18.7).unwrap();
        let json = serde_json::to_string(&soil).unwrap();
        let back: SoilProperties = serde_json::from_str(&json).unwrap();
        assert_eq!(soil, back);
    }
}
