//! **Armor / impact protection sizing** — the minimum protection needed to
//! *stop* a threat, framed strictly defensively.
//!
//! The single question this module answers is: *given a projectile that arrives
//! with kinetic energy `E_k`, what is the **minimum plate thickness** (and the
//! corresponding **areal density**) of a protective plate that just defeats it?*
//! This is the protective-design inverse of a perforation model — it sizes the
//! armor that keeps the threat out. It is **not**, and must not be read as, a
//! tool for defeating armor or optimizing penetration: there is no model here
//! for how to get *through* a plate, only how thick a plate must be to be
//! *safe*.
//!
//! ## Energy-balance ballistic limit (cited)
//!
//! For a normal-impact projectile defeated by **ductile plugging / shear** of a
//! plate (the dominant mechanism for compact projectiles on metallic armor),
//! the plate is at its **ballistic limit** when the energy needed to shear out
//! a plug equals the projectile's kinetic energy. The plug-shear work is
//!
//! ```text
//!   W_plug = (π · d · t) · t · τ_d  =  π · d · t² · τ_d
//! ```
//!
//! the cylindrical shear area `π·d·t` times the plate thickness `t` times the
//! material's dynamic shear strength `τ_d` (the average shear stress over the
//! plug travel). Setting `W_plug ≥ E_k` and solving for `t` gives the **minimum
//! protective thickness**
//!
//! ```text
//!   t_min = √( E_k / (π · d · τ_d) )
//! ```
//!
//! This plug-shear energy balance is the standard textbook ballistic-limit form
//! (Recht & Ipson, "Ballistic Perforation Dynamics," *J. Appl. Mech.*
//! 30:384–390, 1963; Backman & Goldsmith, "The mechanics of penetration of
//! projectiles into targets," *Int. J. Eng. Sci.* 16:1–99, 1978). The dynamic
//! shear strength `τ_d` is taken as an input (≈ 0.5–0.6 × dynamic ultimate
//! tensile strength for ductile metals). Research/educational, validation-
//! pending — the absolute thickness depends strongly on the (uncertain) dynamic
//! material data and on the impact mechanism, so treat the output as a sizing
//! *screen*, not a certified protection rating.
//!
//! Monotonicity the benchmark pins: `t_min` (and hence areal density) **increases
//! with threat kinetic energy** — a faster/heavier threat needs more protection.

use crate::error::SurvivabilityError;
use serde::{Deserialize, Serialize};

/// Validate a finite, strictly-positive quantity, naming it on error.
fn require_pos(name: &str, v: f64) -> Result<f64, SurvivabilityError> {
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(SurvivabilityError::InvalidParameter(format!(
            "{name} must be finite and > 0, got {v}"
        )))
    }
}

/// The threat a protective plate must defeat, expressed by its kinetic energy.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Threat {
    /// Projectile mass `m_p` (kg).
    pub mass_kg: f64,
    /// Projectile impact speed `v` (m/s).
    pub velocity_m_s: f64,
    /// Projectile (or contact) diameter `d` (m) — sets the sheared plug area.
    pub diameter_m: f64,
}

impl Threat {
    /// Build a threat, validating that all three quantities are
    /// finite-and-positive.
    ///
    /// # Errors
    ///
    /// [`SurvivabilityError::InvalidParameter`] on any non-physical input.
    pub fn new(
        mass_kg: f64,
        velocity_m_s: f64,
        diameter_m: f64,
    ) -> Result<Threat, SurvivabilityError> {
        Ok(Threat {
            mass_kg: require_pos("projectile mass", mass_kg)?,
            velocity_m_s: require_pos("impact velocity", velocity_m_s)?,
            diameter_m: require_pos("projectile diameter", diameter_m)?,
        })
    }

    /// The threat's kinetic energy `E_k = ½ m v²` (joules).
    pub fn kinetic_energy_j(&self) -> f64 {
        0.5 * self.mass_kg * self.velocity_m_s * self.velocity_m_s
    }
}

/// The sized protection: minimum plate thickness and areal density that just
/// defeats a [`Threat`] by the plug-shear energy balance.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ArmorSizing {
    /// Threat kinetic energy `E_k` (J) that this plate is sized to stop.
    pub threat_energy_j: f64,
    /// Minimum protective plate thickness `t_min` (m).
    pub min_thickness_m: f64,
    /// Areal density `ρ · t_min` (kg/m²) — the mass-per-area cost of the
    /// protection, the figure of merit for vehicle/structure weight budgets.
    pub areal_density_kg_m2: f64,
}

/// Size the **minimum protective plate** (thickness and areal density) that just
/// defeats `threat`, for a plate of density `plate_density_kg_m3` (kg/m³) and
/// dynamic shear strength `shear_strength_pa` (Pa), via the plug-shear energy
/// balance `t_min = √(E_k / (π d τ_d))`.
///
/// # Errors
///
/// [`SurvivabilityError::InvalidParameter`] if the density or shear strength is
/// not finite-and-positive (the [`Threat`] is already validated at construction;
/// `d > 0` guarantees the divide and the `sqrt` argument are well-defined).
pub fn minimum_protection(
    threat: &Threat,
    plate_density_kg_m3: f64,
    shear_strength_pa: f64,
) -> Result<ArmorSizing, SurvivabilityError> {
    let rho = require_pos("plate density", plate_density_kg_m3)?;
    let tau = require_pos("dynamic shear strength", shear_strength_pa)?;
    let d = threat.diameter_m; // > 0 by Threat invariant
    let e_k = threat.kinetic_energy_j();

    // d, tau, π all > 0 ⇒ denominator > 0; e_k ≥ 0 ⇒ sqrt argument ≥ 0.
    let denom = std::f64::consts::PI * d * tau;
    let t_min = (e_k / denom).sqrt();
    let areal = rho * t_min;

    Ok(ArmorSizing {
        threat_energy_j: e_k,
        min_thickness_m: t_min,
        areal_density_kg_m2: areal,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kinetic_energy_formula() {
        let t = Threat::new(0.01, 800.0, 0.008).unwrap(); // 10 g at 800 m/s
        let e = t.kinetic_energy_j();
        assert!((e - 0.5 * 0.01 * 800.0 * 800.0).abs() < 1e-9); // 3200 J
    }

    #[test]
    fn thickness_increases_with_threat_energy() {
        // PIN: more threat KE ⇒ more protection (monotone increasing).
        let rho = 7850.0; // steel
        let tau = 5.0e8; // ~500 MPa dynamic shear
        let d = 0.008;
        let mut prev = 0.0;
        for v in [200.0, 400.0, 600.0, 800.0, 1000.0, 1500.0] {
            let threat = Threat::new(0.01, v, d).unwrap();
            let s = minimum_protection(&threat, rho, tau).unwrap();
            assert!(
                s.min_thickness_m > prev,
                "thickness not increasing at v={v}"
            );
            assert!(s.areal_density_kg_m2 > 0.0);
            prev = s.min_thickness_m;
        }
    }

    #[test]
    fn areal_density_increases_with_threat_energy() {
        let rho = 2700.0; // aluminium
        let tau = 2.0e8;
        let d = 0.0127;
        let light = minimum_protection(&Threat::new(0.005, 300.0, d).unwrap(), rho, tau).unwrap();
        let heavy = minimum_protection(&Threat::new(0.05, 1000.0, d).unwrap(), rho, tau).unwrap();
        assert!(heavy.areal_density_kg_m2 > light.areal_density_kg_m2);
    }

    #[test]
    fn thickness_matches_closed_form() {
        let rho = 7850.0;
        let tau = 5.0e8;
        let threat = Threat::new(0.01, 800.0, 0.008).unwrap();
        let s = minimum_protection(&threat, rho, tau).unwrap();
        let e_k = threat.kinetic_energy_j();
        let expect = (e_k / (std::f64::consts::PI * 0.008 * tau)).sqrt();
        assert!((s.min_thickness_m - expect).abs() < 1e-12);
        assert!((s.areal_density_kg_m2 - rho * expect).abs() < 1e-9);
    }

    #[test]
    fn degenerate_inputs_error_not_panic() {
        assert!(Threat::new(0.0, 800.0, 0.008).is_err()); // zero mass
        assert!(Threat::new(0.01, 0.0, 0.008).is_err()); // zero velocity
        assert!(Threat::new(0.01, 800.0, 0.0).is_err()); // zero diameter
        assert!(Threat::new(f64::NAN, 800.0, 0.008).is_err());
        let threat = Threat::new(0.01, 800.0, 0.008).unwrap();
        assert!(minimum_protection(&threat, 0.0, 5.0e8).is_err()); // zero density
        assert!(minimum_protection(&threat, 7850.0, 0.0).is_err()); // zero strength
        assert!(minimum_protection(&threat, -1.0, 5.0e8).is_err());
    }

    #[test]
    fn serde_round_trip() {
        let threat = Threat::new(0.01, 800.0, 0.008).unwrap();
        let s = minimum_protection(&threat, 7850.0, 5.0e8).unwrap();
        let json = serde_json::to_string(&s).unwrap();
        let back: ArmorSizing = serde_json::from_str(&json).unwrap();
        assert_eq!(s, back);
    }
}
