//! Archimedes' principle: buoyant force, displaced volume, and the
//! neutral-buoyancy density of a fully submerged body.
//!
//! ## Governing equations
//!
//! For a body that displaces a volume `V` of fluid of density `rho`
//! under gravitational acceleration `g`, the upward buoyant force is
//!
//! ```text
//! F_b = rho * g * V
//! ```
//!
//! A *fully submerged* solid body of mass `m` and volume `V` is in
//! vertical equilibrium (neutrally buoyant) when its weight equals the
//! buoyant force, i.e. when its mean density equals the fluid density:
//!
//! ```text
//! rho_body = m / V = rho_fluid   (neutral)
//! ```
//!
//! It sinks when `rho_body > rho_fluid` and rises when
//! `rho_body < rho_fluid`. All inputs are validated to be finite and
//! within their physical domain before any arithmetic is performed.

use crate::error::{BuoyancyError, Result};

/// Standard gravitational acceleration at Earth's surface, in
/// metres per second squared. Convenience default for callers that do
/// not need a site-specific value.
pub const STANDARD_GRAVITY: f64 = 9.80665;

/// Density of fresh water at roughly 4 degrees Celsius, in kilograms
/// per cubic metre — a common reference fluid density.
pub const RHO_FRESH_WATER: f64 = 1000.0;

/// Representative density of standard sea water, in kilograms per cubic
/// metre.
pub const RHO_SEA_WATER: f64 = 1025.0;

/// Buoyant (upward) force on a body that displaces volume `volume_m3`
/// of fluid of density `rho_fluid` (kg/m^3) under acceleration
/// `gravity` (m/s^2).
///
/// Implements `F_b = rho * g * V`. Returns the force in newtons.
///
/// All three arguments must be finite and strictly positive.
///
/// ```
/// use valenx_buoyancy::archimedes::{buoyant_force, STANDARD_GRAVITY};
/// // 1 m^3 of fresh water displaced -> about 9806.65 N of lift.
/// let f = buoyant_force(1000.0, 1.0, STANDARD_GRAVITY).unwrap();
/// assert!((f - 9806.65).abs() < 1e-9);
/// ```
pub fn buoyant_force(rho_fluid: f64, volume_m3: f64, gravity: f64) -> Result<f64> {
    let rho = BuoyancyError::positive("rho_fluid", rho_fluid)?;
    let v = BuoyancyError::positive("volume_m3", volume_m3)?;
    let g = BuoyancyError::positive("gravity", gravity)?;
    Ok(rho * g * v)
}

/// Volume of fluid (m^3) that must be displaced for the buoyant force
/// to support a given `weight_n` (newtons), in a fluid of density
/// `rho_fluid` (kg/m^3) under acceleration `gravity` (m/s^2).
///
/// Inverts `F_b = rho * g * V` to give `V = F_b / (rho * g)`. This is
/// the submerged volume of a freely floating body whose weight is
/// `weight_n`.
///
/// All arguments must be finite and strictly positive.
///
/// ```
/// use valenx_buoyancy::archimedes::{displaced_volume_for_weight, STANDARD_GRAVITY};
/// // A 9806.65 N weight floating in fresh water displaces ~1 m^3.
/// let v = displaced_volume_for_weight(9806.65, 1000.0, STANDARD_GRAVITY).unwrap();
/// assert!((v - 1.0).abs() < 1e-9);
/// ```
pub fn displaced_volume_for_weight(weight_n: f64, rho_fluid: f64, gravity: f64) -> Result<f64> {
    let w = BuoyancyError::positive("weight_n", weight_n)?;
    let rho = BuoyancyError::positive("rho_fluid", rho_fluid)?;
    let g = BuoyancyError::positive("gravity", gravity)?;
    Ok(w / (rho * g))
}

/// Mean density (kg/m^3) of a solid body of mass `mass_kg` and volume
/// `volume_m3`. This is `rho_body = m / V`, the quantity compared
/// against the fluid density to decide whether a fully submerged body
/// sinks, rises, or hovers.
///
/// Both arguments must be finite and strictly positive.
///
/// ```
/// use valenx_buoyancy::archimedes::mean_density;
/// // 2 kg in 0.002 m^3 -> 1000 kg/m^3.
/// let rho = mean_density(2.0, 0.002).unwrap();
/// assert!((rho - 1000.0).abs() < 1e-9);
/// ```
pub fn mean_density(mass_kg: f64, volume_m3: f64) -> Result<f64> {
    let m = BuoyancyError::positive("mass_kg", mass_kg)?;
    let v = BuoyancyError::positive("volume_m3", volume_m3)?;
    Ok(m / v)
}

/// Vertical-equilibrium verdict for a *fully submerged* solid body,
/// returned by [`submerged_verdict`].
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SubmergedVerdict {
    /// Body density exceeds fluid density: net downward force, body
    /// sinks.
    Sinks,
    /// Body density equals fluid density (within tolerance): neutral
    /// buoyancy, body hovers at any depth.
    Neutral,
    /// Body density is below fluid density: net upward force, body
    /// rises toward the surface.
    Rises,
}

/// Classify a fully submerged body of mean density `rho_body` in a
/// fluid of density `rho_fluid` (both kg/m^3) as
/// [`SubmergedVerdict::Sinks`], [`SubmergedVerdict::Neutral`], or
/// [`SubmergedVerdict::Rises`].
///
/// The densities are treated as equal when they differ by at most
/// `tol_rel` *relative to the fluid density*, i.e. when
/// `|rho_body - rho_fluid| <= tol_rel * rho_fluid`. This is the
/// fully-submerged neutral-buoyancy ground-truth case: equal densities
/// give exact neutral equilibrium.
///
/// `rho_body`, `rho_fluid` and `tol_rel` must be finite; the densities
/// must be strictly positive and `tol_rel` must be non-negative.
///
/// ```
/// use valenx_buoyancy::archimedes::{submerged_verdict, SubmergedVerdict};
/// let v = submerged_verdict(1000.0, 1000.0, 1e-9).unwrap();
/// assert_eq!(v, SubmergedVerdict::Neutral);
/// assert_eq!(submerged_verdict(7850.0, 1000.0, 1e-9).unwrap(), SubmergedVerdict::Sinks);
/// assert_eq!(submerged_verdict(500.0, 1000.0, 1e-9).unwrap(), SubmergedVerdict::Rises);
/// ```
pub fn submerged_verdict(rho_body: f64, rho_fluid: f64, tol_rel: f64) -> Result<SubmergedVerdict> {
    let rho_b = BuoyancyError::positive("rho_body", rho_body)?;
    let rho_f = BuoyancyError::positive("rho_fluid", rho_fluid)?;
    let tol = BuoyancyError::non_negative("tol_rel", tol_rel)?;
    let abs_tol = tol * rho_f;
    if (rho_b - rho_f).abs() <= abs_tol {
        Ok(SubmergedVerdict::Neutral)
    } else if rho_b > rho_f {
        Ok(SubmergedVerdict::Sinks)
    } else {
        Ok(SubmergedVerdict::Rises)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Tolerance for analytic float comparisons.
    const EPS: f64 = 1e-9;

    #[test]
    fn buoyant_force_matches_rho_g_v() {
        // 2 m^3 of sea water (1025) at g = 9.81 -> 2 * 1025 * 9.81.
        let f = buoyant_force(1025.0, 2.0, 9.81).unwrap();
        let expected = 1025.0 * 9.81 * 2.0; // 20110.5
        assert!((f - expected).abs() < EPS);
        assert!((f - 20110.5).abs() < 1e-6);
    }

    #[test]
    fn buoyant_force_scales_linearly_in_volume() {
        let f1 = buoyant_force(1000.0, 1.0, 9.80665).unwrap();
        let f3 = buoyant_force(1000.0, 3.0, 9.80665).unwrap();
        assert!((f3 - 3.0 * f1).abs() < 1e-6);
    }

    #[test]
    fn buoyant_force_is_inverse_of_displaced_volume() {
        // Round-trip: weight -> volume -> force should recover the weight.
        let weight = 5000.0;
        let v = displaced_volume_for_weight(weight, 1000.0, 9.80665).unwrap();
        let f = buoyant_force(1000.0, v, 9.80665).unwrap();
        assert!((f - weight).abs() < 1e-6);
    }

    #[test]
    fn displaced_volume_known_value() {
        // 9806.65 N in fresh water at standard g -> exactly 1 m^3.
        let v = displaced_volume_for_weight(9806.65, RHO_FRESH_WATER, STANDARD_GRAVITY).unwrap();
        assert!((v - 1.0).abs() < EPS);
    }

    #[test]
    fn mean_density_known_value() {
        let rho = mean_density(7.85, 0.001).unwrap();
        assert!((rho - 7850.0).abs() < 1e-6); // steel-ish
    }

    #[test]
    fn fully_submerged_equal_density_is_neutral_ground_truth() {
        // Ground truth: a fully submerged body with rho_body == rho_fluid
        // is in exact neutral equilibrium (zero net vertical force).
        let v = submerged_verdict(1025.0, 1025.0, 1e-12).unwrap();
        assert_eq!(v, SubmergedVerdict::Neutral);
        // And the net force is zero: weight (rho*g*V) minus buoyancy
        // (rho*g*V) cancels exactly for any V.
        let vol = 0.5;
        let weight = 1025.0 * STANDARD_GRAVITY * vol;
        let buoyancy = buoyant_force(1025.0, vol, STANDARD_GRAVITY).unwrap();
        assert!((weight - buoyancy).abs() < EPS);
    }

    #[test]
    fn dense_body_sinks_light_body_rises() {
        assert_eq!(
            submerged_verdict(7850.0, 1000.0, 1e-9).unwrap(),
            SubmergedVerdict::Sinks
        );
        assert_eq!(
            submerged_verdict(600.0, 1000.0, 1e-9).unwrap(),
            SubmergedVerdict::Rises
        );
    }

    #[test]
    fn submerged_verdict_respects_relative_tolerance() {
        // 0.05% denser than the fluid: neutral under a 0.1% tolerance,
        // but sinks under a tight tolerance.
        assert_eq!(
            submerged_verdict(1000.5, 1000.0, 1e-3).unwrap(),
            SubmergedVerdict::Neutral
        );
        assert_eq!(
            submerged_verdict(1000.5, 1000.0, 1e-9).unwrap(),
            SubmergedVerdict::Sinks
        );
    }

    #[test]
    fn rejects_nonpositive_and_nonfinite_inputs() {
        assert!(buoyant_force(0.0, 1.0, 9.81).is_err());
        assert!(buoyant_force(1000.0, -1.0, 9.81).is_err());
        assert!(buoyant_force(1000.0, 1.0, f64::NAN).is_err());
        assert!(displaced_volume_for_weight(-5.0, 1000.0, 9.81).is_err());
        assert!(mean_density(1.0, 0.0).is_err());
        assert!(submerged_verdict(1000.0, 0.0, 1e-9).is_err());
        assert!(submerged_verdict(1000.0, 1000.0, -1e-9).is_err());
    }
}
