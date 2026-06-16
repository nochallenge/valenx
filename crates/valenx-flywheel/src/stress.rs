//! First-order rim (hoop) stress for a thin rotating ring.
//!
//! A thin ring of material density `rho` whose particles move at rim
//! speed `v = omega r` carries a tangential (hoop) tensile stress
//! `sigma = rho v^2 = rho (omega r)^2`, in pascals. This leading-order,
//! radius-independent result sets the burst speed of rim-type flywheels:
//! the energy density a rim can store is capped by `sigma_allow / rho`,
//! which is why high-speed flywheels favour low-density, high-strength
//! materials. Inverting it gives that burst speed,
//! `v_max = sqrt(sigma_allow / rho)` ([`max_rim_speed`]), and the
//! corresponding [`burst_angular_speed`].
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

/// The maximum (burst) rim speed for a given allowable hoop stress,
/// inverting [`rim_stress`]: `v_max = sqrt(sigma_allow / rho)`, in m/s.
///
/// This is the rim speed at which a thin ring of density `density`
/// reaches the allowable tangential stress `allowable_stress` — the
/// leading-order burst-speed limit that caps a rim-type flywheel's
/// energy density. It is radius-independent (a direct consequence of
/// `sigma = rho v^2`), and falls as the square root of density, which is
/// why low-density, high-strength rims spin fastest.
///
/// # Errors
///
/// Returns [`FlywheelError::InvalidParameter`] if `allowable_stress` is
/// negative / non-finite or `density` is not strictly positive.
pub fn max_rim_speed(allowable_stress: f64, density: f64) -> Result<f64, FlywheelError> {
    let sigma = FlywheelError::require_non_negative("allowable_stress", allowable_stress)?;
    let density = FlywheelError::require_positive("density", density)?;
    Ok((sigma / density).sqrt())
}

/// The maximum (burst) angular speed of a thin ring of the given radius
/// and density at an allowable hoop stress, in rad/s:
/// `omega_max = sqrt(sigma_allow / rho) / r = v_max / r`.
///
/// Combines [`max_rim_speed`] with `v = omega r` to give the spin speed
/// at which the rim reaches `allowable_stress`.
///
/// # Errors
///
/// Returns [`FlywheelError::InvalidParameter`] if `allowable_stress` is
/// negative / non-finite, or `density` or `radius` is not strictly
/// positive.
pub fn burst_angular_speed(
    allowable_stress: f64,
    density: f64,
    radius: f64,
) -> Result<f64, FlywheelError> {
    let radius = FlywheelError::require_positive("radius", radius)?;
    let v_max = max_rim_speed(allowable_stress, density)?;
    Ok(v_max / radius)
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

    #[test]
    fn max_rim_speed_inverts_rim_stress() {
        // The burst angular speed, fed back through rim_stress, reproduces
        // the allowable stress exactly.
        let (rho, r, sigma_allow) = (7800.0, 0.5, 19_500_000.0);
        let omega_max = burst_angular_speed(sigma_allow, rho, r).unwrap();
        let sigma_back = rim_stress(rho, omega_max, r).unwrap();
        assert!(
            (sigma_back - sigma_allow).abs() / sigma_allow < 1e-9,
            "got {sigma_back}"
        );
    }

    #[test]
    fn max_rim_speed_closed_form() {
        // sqrt(sigma/rho): sigma = 7800*2500 = 19.5e6, rho = 7800 -> 50 m/s.
        let v = max_rim_speed(19_500_000.0, 7800.0).unwrap();
        assert!((v - 50.0).abs() < 1e-6, "got {v}");
    }

    #[test]
    fn burst_omega_equals_vmax_over_radius() {
        let (sigma, rho, r) = (3.0e8, 2700.0, 0.25);
        let v_max = max_rim_speed(sigma, rho).unwrap();
        let omega_max = burst_angular_speed(sigma, rho, r).unwrap();
        assert!((omega_max - v_max / r).abs() < 1e-9);
        // rim_speed at omega_max recovers v_max.
        assert!((rim_speed(omega_max, r).unwrap() - v_max).abs() < 1e-9);
    }

    #[test]
    fn lower_density_bursts_at_higher_speed() {
        // v_max ~ 1/sqrt(rho): a lighter rim spins faster before bursting.
        let sigma = 5.0e8;
        let steel = max_rim_speed(sigma, 7800.0).unwrap();
        let composite = max_rim_speed(sigma, 1600.0).unwrap();
        assert!(composite > steel, "composite {composite} vs steel {steel}");
        assert!((composite / steel - (7800.0_f64 / 1600.0).sqrt()).abs() < 1e-9);
    }

    #[test]
    fn zero_stress_gives_zero_speed() {
        assert!(max_rim_speed(0.0, 7800.0).unwrap().abs() < 1e-12);
    }

    #[test]
    fn burst_speed_rejects_bad_inputs() {
        assert!(max_rim_speed(-1.0, 7800.0).is_err()); // negative stress
        assert!(max_rim_speed(1.0e6, 0.0).is_err()); // density <= 0
        assert!(max_rim_speed(f64::NAN, 7800.0).is_err()); // non-finite
        assert!(burst_angular_speed(1.0e6, 7800.0, 0.0).is_err()); // radius <= 0
    }
}
