//! Mass-estimating / vehicle-sizing relations.
//!
//! The closed-form algebra that turns a `Δv` budget and a stage's
//! propulsive + structural figures into the numbers that size a rocket:
//! the **mass ratio** from the Tsiolkovsky rocket equation, the
//! **structural mass fraction** `ε`, the **deliverable payload fraction**
//! of a single stage, and **propellant-tank volumes** from the
//! propellant densities.
//!
//! Every relation here is an exact algebraic / Tsiolkovsky result — no
//! integration — so the unit tests pin the values against the textbook
//! formulas directly.
//!
//! ## The single-stage payload relation
//!
//! Split the lift-off mass `m₀` into structure `m_s`, propellant `m_p`
//! and payload `m_pl`. The **structural mass fraction** is taken on the
//! *stage* (structure + propellant, excluding payload):
//!
//! ```text
//!   ε = m_s / (m_s + m_p)
//! ```
//!
//! Burning all of `m_p` gives a mass ratio `MR = m₀ / m_f` with
//! `m_f = m_s + m_pl`, and the rocket equation ties `MR` to the `Δv`:
//!
//! ```text
//!   MR = exp( Δv / (Isp · g₀) )
//! ```
//!
//! Eliminating the masses (normalising `m₀ = 1`) gives the deliverable
//! **payload fraction**
//!
//! ```text
//!   λ = m_pl / m₀ = (1 − ε·MR) / (MR · (1 − ε))
//! ```
//!
//! which is positive only when `ε·MR < 1` (equivalently `MR < 1/ε`). A
//! `Δv` budget too large for the structural fraction therefore has **no
//! positive payload solution** — a real, testable infeasibility that this
//! module reports as an [`AstroError::InvalidParameter`] rather than a
//! meaningless negative fraction.

use crate::constants::G0;
use crate::error::{AstroError, Result};

/// Bulk density of liquid oxygen (LOX) at its boiling point (kg/m³).
pub const RHO_LOX: f64 = 1_141.0;
/// Bulk density of RP-1 (refined kerosene) propellant (kg/m³).
pub const RHO_RP1: f64 = 810.0;
/// Bulk density of liquid hydrogen (LH2) at its boiling point (kg/m³).
pub const RHO_LH2: f64 = 71.0;
/// Bulk density of liquid methane (LCH4) at its boiling point (kg/m³).
pub const RHO_LCH4: f64 = 423.0;

/// Tsiolkovsky **mass ratio** `MR = m₀/m_f = exp(Δv / (Isp · g₀))`.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if `delta_v` is non-finite or
/// negative, or `isp` is non-finite or non-positive — either would make
/// the `Δv/(Isp·g₀)` term or its `exp` a silent `NaN`/`Inf`.
pub fn mass_ratio(delta_v: f64, isp: f64) -> Result<f64> {
    if !delta_v.is_finite() || delta_v < 0.0 {
        return Err(AstroError::InvalidParameter(
            "delta_v must be finite and >= 0",
        ));
    }
    if !isp.is_finite() || isp <= 0.0 {
        return Err(AstroError::InvalidParameter("isp must be finite and > 0"));
    }
    Ok((delta_v / (isp * G0)).exp())
}

/// The `Δv` (m/s) a stage of mass ratio `MR` delivers at the given `isp`
/// (s): the inverse of [`mass_ratio`], `Δv = Isp · g₀ · ln(MR)`.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if `mass_ratio` is non-finite
/// or `< 1` (a mass ratio below one is non-physical and makes `ln` give a
/// negative `Δv`), or `isp` is non-finite or non-positive.
pub fn delta_v_from_mass_ratio(mass_ratio: f64, isp: f64) -> Result<f64> {
    if !mass_ratio.is_finite() || mass_ratio < 1.0 {
        return Err(AstroError::InvalidParameter("mass_ratio must be >= 1"));
    }
    if !isp.is_finite() || isp <= 0.0 {
        return Err(AstroError::InvalidParameter("isp must be finite and > 0"));
    }
    Ok(isp * G0 * mass_ratio.ln())
}

/// Structural mass fraction `ε = m_dry / (m_dry + m_prop)` of a stage.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if either mass is non-finite
/// or negative, or both are zero (a zero total mass makes the ratio
/// `0/0 = NaN`).
pub fn structural_fraction(dry_mass: f64, propellant_mass: f64) -> Result<f64> {
    if !dry_mass.is_finite() || dry_mass < 0.0 {
        return Err(AstroError::InvalidParameter(
            "dry_mass must be finite and >= 0",
        ));
    }
    if !propellant_mass.is_finite() || propellant_mass < 0.0 {
        return Err(AstroError::InvalidParameter(
            "propellant_mass must be finite and >= 0",
        ));
    }
    let total = dry_mass + propellant_mass;
    if total <= 0.0 {
        return Err(AstroError::InvalidParameter(
            "total stage mass must be > 0",
        ));
    }
    Ok(dry_mass / total)
}

/// Deliverable **payload fraction** `λ = m_pl/m₀` of a single stage that
/// burns its full propellant load to achieve `delta_v`, given the engine
/// `isp` (s) and the stage structural fraction `epsilon`
/// (`ε = m_dry/(m_dry+m_prop)`, excluding payload):
///
/// ```text
///   λ = (1 − ε·MR) / (MR · (1 − ε)),   MR = exp(Δv/(Isp·g₀))
/// ```
///
/// # Errors
///
/// - [`AstroError::InvalidParameter`] if `epsilon` is non-finite or not in
///   the open interval `(0, 1)`, or for a bad `delta_v`/`isp` (see
///   [`mass_ratio`]).
/// - [`AstroError::InvalidParameter`] (`"infeasible"`) when the budget has
///   **no positive payload** — i.e. `ε·MR ≥ 1` (`MR ≥ 1/ε`). The stage's
///   structure alone already exceeds what the mass ratio allows, so no
///   payload can be carried. This is a correct, deliberate result, not a
///   numerical failure.
pub fn payload_fraction(delta_v: f64, isp: f64, epsilon: f64) -> Result<f64> {
    if !epsilon.is_finite() || epsilon <= 0.0 || epsilon >= 1.0 {
        return Err(AstroError::InvalidParameter(
            "structural fraction epsilon must be in (0, 1)",
        ));
    }
    let mr = mass_ratio(delta_v, isp)?;
    // Positive payload requires ε·MR < 1; at/above that the structure
    // alone exhausts the mass-ratio allowance.
    if epsilon * mr >= 1.0 {
        return Err(AstroError::InvalidParameter(
            "infeasible: delta-v too high for this structural fraction (no positive payload)",
        ));
    }
    Ok((1.0 - epsilon * mr) / (mr * (1.0 - epsilon)))
}

/// Propellant-tank volume (m³) for a propellant mass (kg) at bulk density
/// `rho` (kg/m³): `V = m_prop / ρ`.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if `propellant_mass` is
/// non-finite or negative, or `rho` is non-finite or non-positive (which
/// would make the division a silent `NaN`/`Inf`).
pub fn tank_volume(propellant_mass: f64, rho: f64) -> Result<f64> {
    if !propellant_mass.is_finite() || propellant_mass < 0.0 {
        return Err(AstroError::InvalidParameter(
            "propellant_mass must be finite and >= 0",
        ));
    }
    if !rho.is_finite() || rho <= 0.0 {
        return Err(AstroError::InvalidParameter(
            "propellant density must be finite and > 0",
        ));
    }
    Ok(propellant_mass / rho)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mass_ratio_matches_tsiolkovsky() {
        // Δv = 9000 m/s, Isp = 350 s. v_e = 350·9.80665 = 3432.3275 m/s.
        // MR = exp(9000/3432.3275) = 13.764976134...
        let mr = mass_ratio(9_000.0, 350.0).expect("valid");
        let expected = (9_000.0 / (350.0 * G0)).exp();
        assert!((mr - expected).abs() < 1e-12);
        assert!((mr - 13.764_976_134_343_75).abs() < 1e-9, "MR = {mr}");
    }

    #[test]
    fn mass_ratio_and_delta_v_are_inverses() {
        let mr = mass_ratio(7_800.0, 311.0).expect("valid");
        let dv = delta_v_from_mass_ratio(mr, 311.0).expect("valid");
        assert!((dv - 7_800.0).abs() < 1e-6, "round-trip dv = {dv}");
    }

    #[test]
    fn delta_v_known_value() {
        // MR = 10, Isp = 300 -> Δv = 300·g₀·ln(10).
        let dv = delta_v_from_mass_ratio(10.0, 300.0).expect("valid");
        let expected = 300.0 * G0 * 10.0_f64.ln();
        assert!((dv - expected).abs() < 1e-9, "dv = {dv}");
    }

    #[test]
    fn structural_fraction_basic() {
        // m_dry = 8, m_prop = 92 -> ε = 0.08.
        let e = structural_fraction(8.0, 92.0).expect("valid");
        assert!((e - 0.08).abs() < 1e-12, "eps = {e}");
    }

    #[test]
    fn pinned_single_stage_9000_350_008_is_infeasible() {
        // ORACLE: Δv=9000, Isp=350, ε=0.08 -> MR=13.765 > 1/ε=12.5, so
        // the closed-form payload fraction is negative: NO positive
        // payload. The function must report this as an Err, not return a
        // nonsense negative number.
        let r = payload_fraction(9_000.0, 350.0, 0.08);
        assert!(
            matches!(r, Err(AstroError::InvalidParameter(_))),
            "expected infeasible Err, got {r:?}"
        );
    }

    #[test]
    fn payload_fraction_feasible_pinned() {
        // A feasible budget: Δv=3000, Isp=350, ε=0.08.
        // MR = exp(3000/3432.3275) = 2.3965793941...
        // λ = (1 − 0.08·MR)/(MR·0.92) = 0.36658844506...
        let lam = payload_fraction(3_000.0, 350.0, 0.08).expect("feasible");
        let mr = (3_000.0 / (350.0 * G0)).exp();
        let expected = (1.0 - 0.08 * mr) / (mr * 0.92);
        assert!((lam - expected).abs() < 1e-12);
        assert!((lam - 0.366_588_445_060_05).abs() < 1e-9, "lambda = {lam}");
        assert!((0.0..1.0).contains(&lam));
    }

    #[test]
    fn payload_fraction_at_feasibility_boundary_is_err() {
        // Exactly at ε·MR = 1 the payload is zero; the strict `>= 1`
        // guard treats the boundary itself as infeasible (no *positive*
        // payload). Construct Δv so MR = 1/ε exactly.
        let eps = 0.1;
        let mr_boundary = 1.0 / eps; // = 10
        let dv = delta_v_from_mass_ratio(mr_boundary, 300.0).expect("valid");
        assert!(payload_fraction(dv, 300.0, eps).is_err());
        // Just under the boundary is feasible and gives a tiny payload.
        let lam = payload_fraction(dv * 0.99, 300.0, eps).expect("just feasible");
        assert!(lam > 0.0 && lam < 0.01, "near-boundary lambda = {lam}");
    }

    #[test]
    fn tank_volumes_match_densities() {
        // 100 t of each propellant.
        assert!((tank_volume(100_000.0, RHO_LOX).expect("ok") - 100_000.0 / 1_141.0).abs() < 1e-9);
        assert!((tank_volume(100_000.0, RHO_RP1).expect("ok") - 100_000.0 / 810.0).abs() < 1e-9);
        // LH2 is very low density -> a much larger tank than LOX.
        let v_lh2 = tank_volume(100_000.0, RHO_LH2).expect("ok");
        let v_lox = tank_volume(100_000.0, RHO_LOX).expect("ok");
        assert!(v_lh2 > v_lox * 10.0, "LH2 {v_lh2} vs LOX {v_lox}");
        // Methane sits between RP-1 and LH2.
        let v_ch4 = tank_volume(100_000.0, RHO_LCH4).expect("ok");
        assert!(v_ch4 > tank_volume(100_000.0, RHO_RP1).expect("ok"));
    }

    #[test]
    fn rejects_non_physical_inputs() {
        assert!(mass_ratio(f64::NAN, 300.0).is_err());
        assert!(mass_ratio(-1.0, 300.0).is_err());
        assert!(mass_ratio(1_000.0, 0.0).is_err());
        assert!(delta_v_from_mass_ratio(0.5, 300.0).is_err()); // MR < 1
        assert!(structural_fraction(0.0, 0.0).is_err()); // 0/0
        assert!(structural_fraction(-1.0, 10.0).is_err());
        assert!(payload_fraction(3_000.0, 350.0, 0.0).is_err()); // ε not in (0,1)
        assert!(payload_fraction(3_000.0, 350.0, 1.0).is_err());
        assert!(tank_volume(100.0, 0.0).is_err());
        assert!(tank_volume(-1.0, 800.0).is_err());
    }
}
