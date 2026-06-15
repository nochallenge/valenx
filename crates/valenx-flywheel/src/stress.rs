//! First-order rim (hoop) stress for a thin rotating ring.
//!
//! A thin ring of material density `rho` whose particles move at rim
//! speed `v = omega r` carries a tangential (hoop) tensile stress
//! `sigma = rho v^2 = rho (omega r)^2`, in pascals. This leading-order,
//! radius-independent result sets the burst speed of rim-type flywheels:
//! the energy density a rim can store is capped by `sigma_allow / rho`,
//! which is why high-speed flywheels favour low-density, high-strength
//! materials.
//!
//! See the crate-level docs for the assumptions this elementary model
//! omits (radial stress, Poisson coupling, thickness gradients, hub
//! interface, anisotropy, fatigue).

use crate::error::FlywheelError;

/// Rim (tangential) speed `v = omega r`, in metres per second.
///
/// `omega` is the angular speed (rad/s) and `radius` the rim radius (m).
///
/// # Errors
///
/// Returns [`FlywheelError::InvalidParameter`] if `omega` is negative /
/// non-finite or `radius` is not strictly positive.
pub fn rim_speed(omega: f64, radius: f64) -> Result<f64, FlywheelError> {
    let omega = FlywheelError::require_non_negative("omega", omega)?;
    let radius = FlywheelError::require_positive("radius", radius)?;
    Ok(omega * radius)
}

/// Thin-ring hoop stress `sigma = rho (omega r)^2 = rho v^2`, in pascals.
///
/// `density` is the material density `rho` (kg/m^3), `omega` the angular
/// speed (rad/s), and `radius` the rim radius (m).
///
/// # Errors
///
/// Returns [`FlywheelError::InvalidParameter`] if `density` or `radius`
/// is not strictly positive, or `omega` is negative / non-finite.
pub fn rim_stress(density: f64, omega: f64, radius: f64) -> Result<f64, FlywheelError> {
    let density = FlywheelError::require_positive("density", density)?;
    let v = rim_speed(omega, radius)?;
    Ok(density * v * v)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn rim_speed_is_omega_times_radius() {
        // omega = 50 rad/s, r = 0.4 m -> v = 20 m/s
        let v = rim_speed(50.0, 0.4).unwrap();
        assert!((v - 20.0).abs() < EPS);
    }

    #[test]
    fn rim_stress_equals_rho_v_squared() {
        // rho = 7800 (steel), omega = 100, r = 0.5 -> v = 50,
        // sigma = 7800 * 2500 = 19_500_000 Pa
        let sigma = rim_stress(7800.0, 100.0, 0.5).unwrap();
        let v = rim_speed(100.0, 0.5).unwrap();
        assert!((sigma - 7800.0 * v * v).abs() < 1e-3);
        assert!((sigma - 19_500_000.0).abs() < 1e-3);
    }

    #[test]
    fn doubling_speed_quadruples_rim_stress() {
        // sigma ~ omega^2, mirroring the energy scaling.
        let rho = 2700.0;
        let r = 0.3;
        let omega = 80.0;
        let s1 = rim_stress(rho, omega, r).unwrap();
        let s2 = rim_stress(rho, 2.0 * omega, r).unwrap();
        assert!((s2 / s1 - 4.0).abs() < 1e-9);
    }

    #[test]
    fn stress_is_independent_of_a_thin_ring_radius_at_fixed_rim_speed() {
        // Two rings with the same rim speed v but different radii carry
        // the same hoop stress, since sigma = rho v^2.
        let rho = 1600.0; // a composite-ish density
                          // Pair A: r = 0.2, choose omega so v = 60.
        let sa = rim_stress(rho, 300.0, 0.2).unwrap();
        // Pair B: r = 0.6, omega = 100 -> v = 60 as well.
        let sb = rim_stress(rho, 100.0, 0.6).unwrap();
        assert!((sa - sb).abs() < 1e-6);
    }

    #[test]
    fn stress_rejects_bad_inputs() {
        assert!(rim_stress(0.0, 10.0, 1.0).is_err());
        assert!(rim_stress(1000.0, -1.0, 1.0).is_err());
        assert!(rim_stress(1000.0, 10.0, 0.0).is_err());
        assert!(rim_speed(10.0, -1.0).is_err());
    }
}
