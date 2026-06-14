//! Steady laminar flow through a single rigid cylindrical vessel.
//!
//! These are the textbook Hagen-Poiseuille relations for the steady,
//! fully-developed, laminar flow of an incompressible Newtonian fluid
//! through a straight rigid tube of circular cross-section. They are
//! exact for that idealised case and are the standard first-order model
//! of resistance in a vascular segment.
//!
//! All quantities are in SI base units (metres, pascals, seconds,
//! pascal-seconds), so a result in [`poiseuille_flow`] is in m^3/s. The
//! relations are scale-free, so any consistent unit system works; the
//! `_si` helpers below just document the SI choice. The caller is
//! responsible for converting to physiological units (mmHg, mL/min,
//! centipoise) at the boundary.
//!
//! # Models
//!
//! For a tube of radius `r`, length `L`, fluid dynamic viscosity `mu`
//! and pressure drop `dP = P_in - P_out` along it:
//!
//! - Volumetric flow ([`poiseuille_flow`]):
//!   `Q = pi * r^4 * dP / (8 * mu * L)`.
//! - Hydraulic resistance ([`vascular_resistance`]):
//!   `R = 8 * mu * L / (pi * r^4)`, so that `Q = dP / R`.
//! - Wall shear stress ([`wall_shear_stress`]):
//!   `tau = 4 * mu * Q / (pi * r^3)`.
//!
//! The fourth-power dependence on radius is the headline physiological
//! fact: halving a vessel's radius cuts its flow (at fixed `dP`) by a
//! factor of sixteen and raises its resistance by the same factor.

use std::f64::consts::PI;

use crate::error::{require_non_negative, require_positive};
use crate::HemodynamicsError;

/// Hagen-Poiseuille volumetric flow rate through a rigid cylindrical
/// vessel.
///
/// Computes `Q = pi * r^4 * dP / (8 * mu * L)`, the steady laminar flow
/// driven by a pressure drop `dP` across a tube of radius `radius` and
/// length `length` for a fluid of dynamic viscosity `viscosity`.
///
/// # Units
///
/// With SI inputs (radius m, pressure drop Pa, viscosity Pa·s, length
/// m) the result is in m^3/s.
///
/// # Errors
///
/// Returns [`HemodynamicsError::NonPositive`] if `radius`, `viscosity`
/// or `length` is not strictly positive, or [`HemodynamicsError::NotFinite`]
/// for a non-finite input. The pressure drop may be any finite value
/// (a negative `dP` simply yields a negative — reversed — flow).
pub fn poiseuille_flow(
    radius: f64,
    pressure_drop: f64,
    viscosity: f64,
    length: f64,
) -> Result<f64, HemodynamicsError> {
    let r = require_positive("radius", radius)?;
    let mu = require_positive("viscosity", viscosity)?;
    let l = require_positive("length", length)?;
    if !pressure_drop.is_finite() {
        return Err(HemodynamicsError::NotFinite {
            name: "pressure_drop",
            value: pressure_drop,
        });
    }
    Ok(PI * r.powi(4) * pressure_drop / (8.0 * mu * l))
}

/// Hydraulic (vascular) resistance of a rigid cylindrical vessel.
///
/// Computes `R = 8 * mu * L / (pi * r^4)`. This is the proportionality
/// constant in the Ohm-analogue `dP = Q * R`: dividing a pressure drop
/// by this resistance reproduces [`poiseuille_flow`] exactly.
///
/// # Units
///
/// With SI inputs the result is in Pa·s/m^3.
///
/// # Errors
///
/// Returns [`HemodynamicsError::NonPositive`] / [`HemodynamicsError::NotFinite`]
/// if `viscosity`, `length` or `radius` is not a strictly-positive
/// finite number.
pub fn vascular_resistance(
    viscosity: f64,
    length: f64,
    radius: f64,
) -> Result<f64, HemodynamicsError> {
    let mu = require_positive("viscosity", viscosity)?;
    let l = require_positive("length", length)?;
    let r = require_positive("radius", radius)?;
    Ok(8.0 * mu * l / (PI * r.powi(4)))
}

/// Flow rate from a pressure drop and a resistance via the Ohm analogue
/// `Q = dP / R`.
///
/// This is the resistance-form of [`poiseuille_flow`]; the two agree
/// when `R` is the [`vascular_resistance`] of the same vessel.
///
/// # Errors
///
/// Returns [`HemodynamicsError::NonPositive`] / [`HemodynamicsError::NotFinite`]
/// if `resistance` is not a strictly-positive finite number, or
/// [`HemodynamicsError::NotFinite`] if `pressure_drop` is non-finite.
pub fn flow_from_resistance(pressure_drop: f64, resistance: f64) -> Result<f64, HemodynamicsError> {
    let r = require_positive("resistance", resistance)?;
    if !pressure_drop.is_finite() {
        return Err(HemodynamicsError::NotFinite {
            name: "pressure_drop",
            value: pressure_drop,
        });
    }
    Ok(pressure_drop / r)
}

/// Wall shear stress exerted by a Poiseuille flow on the vessel wall.
///
/// Computes `tau = 4 * mu * Q / (pi * r^3)`, the magnitude of the
/// viscous traction the flowing fluid applies tangentially to the tube
/// wall. Equivalently `tau = mu * gamma_w` with the wall shear rate
/// `gamma_w = 4 Q / (pi r^3)`. Wall shear stress is the mechanical
/// stimulus that endothelium senses, so it rises with flow and falls
/// steeply as the vessel dilates.
///
/// # Units
///
/// With SI inputs (viscosity Pa·s, flow m^3/s, radius m) the result is
/// in pascals.
///
/// # Errors
///
/// Returns [`HemodynamicsError::NonPositive`] / [`HemodynamicsError::NotFinite`]
/// if `viscosity` or `radius` is not a strictly-positive finite number,
/// or [`HemodynamicsError::Negative`] / [`HemodynamicsError::NotFinite`]
/// if `flow` is negative or non-finite.
pub fn wall_shear_stress(viscosity: f64, flow: f64, radius: f64) -> Result<f64, HemodynamicsError> {
    let mu = require_positive("viscosity", viscosity)?;
    let q = require_non_negative("flow", flow)?;
    let r = require_positive("radius", radius)?;
    Ok(4.0 * mu * q / (PI * r.powi(3)))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference physiological-ish numbers (arbitrary but consistent SI):
    /// a 2 mm-radius arteriole, 0.05 m long, blood viscosity 3.5 mPa·s,
    /// driven by 1000 Pa (~7.5 mmHg).
    const R0: f64 = 2.0e-3;
    const L0: f64 = 0.05;
    const MU0: f64 = 3.5e-3;
    const DP0: f64 = 1000.0;

    #[test]
    fn poiseuille_matches_closed_form() {
        let q = poiseuille_flow(R0, DP0, MU0, L0).expect("valid");
        let expected = PI * R0.powi(4) * DP0 / (8.0 * MU0 * L0);
        assert!((q - expected).abs() < 1e-15 * expected.abs().max(1.0));
        // Sanity: positive driving pressure gives positive flow.
        assert!(q > 0.0);
    }

    #[test]
    fn flow_scales_with_radius_to_the_fourth() {
        // VALIDATE: doubling the radius multiplies flow by 2^4 = 16.
        let q1 = poiseuille_flow(R0, DP0, MU0, L0).expect("valid");
        let q2 = poiseuille_flow(2.0 * R0, DP0, MU0, L0).expect("valid");
        let ratio = q2 / q1;
        assert!(
            (ratio - 16.0).abs() < 1e-9,
            "expected 16x, got {ratio} (q1={q1}, q2={q2})"
        );

        // Tripling the radius multiplies flow by 3^4 = 81.
        let q3 = poiseuille_flow(3.0 * R0, DP0, MU0, L0).expect("valid");
        let ratio3 = q3 / q1;
        assert!((ratio3 - 81.0).abs() < 1e-9, "expected 81x, got {ratio3}");
    }

    #[test]
    fn flow_is_linear_in_pressure_drop() {
        let q1 = poiseuille_flow(R0, DP0, MU0, L0).expect("valid");
        let q2 = poiseuille_flow(R0, 3.0 * DP0, MU0, L0).expect("valid");
        assert!((q2 / q1 - 3.0).abs() < 1e-12);
    }

    #[test]
    fn flow_reverses_with_negative_pressure() {
        let qf = poiseuille_flow(R0, DP0, MU0, L0).expect("valid");
        let qr = poiseuille_flow(R0, -DP0, MU0, L0).expect("valid");
        assert!((qf + qr).abs() < 1e-15 * qf.abs().max(1.0));
        assert!(qr < 0.0);
    }

    #[test]
    fn resistance_equals_pressure_over_flow() {
        // VALIDATE: R = dP / Q.
        let q = poiseuille_flow(R0, DP0, MU0, L0).expect("valid");
        let r = vascular_resistance(MU0, L0, R0).expect("valid");
        let r_from_ratio = DP0 / q;
        assert!(
            (r - r_from_ratio).abs() < 1e-6 * r,
            "R={r}, dP/Q={r_from_ratio}"
        );
    }

    #[test]
    fn flow_from_resistance_round_trips_poiseuille() {
        let q_direct = poiseuille_flow(R0, DP0, MU0, L0).expect("valid");
        let r = vascular_resistance(MU0, L0, R0).expect("valid");
        let q_ohm = flow_from_resistance(DP0, r).expect("valid");
        assert!((q_direct - q_ohm).abs() < 1e-12 * q_direct);
    }

    #[test]
    fn resistance_scales_inverse_fourth_power_of_radius() {
        // Doubling radius cuts resistance to 1/16.
        let r1 = vascular_resistance(MU0, L0, R0).expect("valid");
        let r2 = vascular_resistance(MU0, L0, 2.0 * R0).expect("valid");
        assert!((r1 / r2 - 16.0).abs() < 1e-9, "expected 16x, got {}", r1 / r2);
    }

    #[test]
    fn wall_shear_stress_increases_with_flow() {
        // VALIDATE: tau rises with Q (linearly, at fixed geometry).
        let q = poiseuille_flow(R0, DP0, MU0, L0).expect("valid");
        let tau1 = wall_shear_stress(MU0, q, R0).expect("valid");
        let tau2 = wall_shear_stress(MU0, 2.0 * q, R0).expect("valid");
        assert!(tau2 > tau1);
        assert!(
            (tau2 / tau1 - 2.0).abs() < 1e-12,
            "expected 2x, got {}",
            tau2 / tau1
        );
    }

    #[test]
    fn wall_shear_stress_matches_closed_form() {
        let q = 1.0e-6;
        let tau = wall_shear_stress(MU0, q, R0).expect("valid");
        let expected = 4.0 * MU0 * q / (PI * R0.powi(3));
        assert!((tau - expected).abs() < 1e-18 * expected.abs().max(1.0));
    }

    #[test]
    fn wall_shear_stress_zero_flow_is_zero() {
        let tau = wall_shear_stress(MU0, 0.0, R0).expect("valid");
        assert!(tau.abs() < 1e-18);
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        assert!(poiseuille_flow(0.0, DP0, MU0, L0).is_err());
        assert!(poiseuille_flow(R0, DP0, -MU0, L0).is_err());
        assert!(poiseuille_flow(R0, DP0, MU0, 0.0).is_err());
        assert!(poiseuille_flow(R0, f64::NAN, MU0, L0).is_err());
        assert!(vascular_resistance(MU0, L0, 0.0).is_err());
        assert!(flow_from_resistance(DP0, 0.0).is_err());
        assert!(wall_shear_stress(MU0, -1.0, R0).is_err());
    }
}
