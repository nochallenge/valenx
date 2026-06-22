//! Black-hole thermodynamics for the Kerr–Newman family.
//!
//! Geometrized units `G = c = 1`, and additionally `ħ = k_B = 1` for the
//! quantum quantities (Hawking temperature, entropy). Use [`crate::units`] to
//! express these in SI for a real black hole.

use std::f64::consts::PI;

use crate::observables::horizons;
use crate::spacetimes::KerrNewman;
use crate::Result;

/// Horizon thermodynamic quantities of a black hole.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Thermodynamics {
    /// Surface gravity `κ = (r+ − r−) / (2(r+² + a²))` at the outer horizon.
    pub surface_gravity: f64,
    /// Hawking temperature `T_H = κ / 2π` (with `ħ = k_B = 1`).
    pub hawking_temperature: f64,
    /// Outer-horizon area `A = 4π(r+² + a²)`.
    pub horizon_area: f64,
    /// Bekenstein–Hawking entropy `S = A / 4`.
    pub entropy: f64,
    /// Angular velocity of the horizon `Ω_H = a / (r+² + a²)`.
    pub horizon_angular_velocity: f64,
}

/// Compute the outer-horizon thermodynamics of a (sub-extremal) black hole.
///
/// For Schwarzschild this gives `κ = 1/4M`, `T_H = 1/8πM`, `S = 4πM²`,
/// `Ω_H = 0`; at the extremal limit (`a² + Q² = M²`) the surface gravity and
/// Hawking temperature vanish.
///
/// # Errors
/// Propagates [`crate::RelativityError::InvalidParameter`] from [`horizons`]
/// for a non-positive-mass or super-extremal hole.
pub fn thermodynamics(bh: &KerrNewman) -> Result<Thermodynamics> {
    let h = horizons(bh)?;
    let a = bh.spin;
    let denom = h.outer * h.outer + a * a;
    let kappa = (h.outer - h.inner) / (2.0 * denom);
    let area = 4.0 * PI * denom;
    Ok(Thermodynamics {
        surface_gravity: kappa,
        hawking_temperature: kappa / (2.0 * PI),
        horizon_area: area,
        entropy: area / 4.0,
        horizon_angular_velocity: a / denom,
    })
}
