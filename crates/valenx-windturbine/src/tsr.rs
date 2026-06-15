//! Tip-speed ratio — the ratio of blade-tip speed to wind speed.
//!
//! # Model
//!
//! The tip-speed ratio is the dimensionless
//!
//! ```text
//! lambda = omega * R / v
//! ```
//!
//! where `omega` is the rotor angular velocity (rad/s), `R` the rotor
//! (tip) radius (m), and `v` the free-stream wind speed (m/s). It sets
//! where on its `Cp(lambda)` curve a turbine operates; modern
//! three-blade machines peak near `lambda ~ 6..8`.
//!
//! Two convenience converters relate angular velocity in rad/s to the
//! rev/min figure quoted on a nameplate, via `omega = 2 pi * rpm / 60`.

use crate::error::WindTurbineError;

/// Convert rotor speed from rev/min to rad/s: `omega = 2 pi * rpm / 60`.
///
/// # Errors
///
/// [`WindTurbineError::BadParameter`] if `rpm` is negative or
/// non-finite.
///
/// ```
/// use valenx_windturbine::tsr::rpm_to_rad_per_s;
/// // 60 rpm = 1 rev/s = 2*pi rad/s
/// let w = rpm_to_rad_per_s(60.0).unwrap();
/// assert!((w - std::f64::consts::TAU).abs() < 1e-12);
/// ```
pub fn rpm_to_rad_per_s(rpm: f64) -> Result<f64, WindTurbineError> {
    if !rpm.is_finite() || rpm < 0.0 {
        return Err(WindTurbineError::BadParameter {
            name: "rpm",
            reason: "must be >= 0".to_string(),
        });
    }
    Ok(rpm * std::f64::consts::TAU / 60.0)
}

/// Convert rotor speed from rad/s to rev/min: `rpm = 60 * omega / 2 pi`.
///
/// The inverse of [`rpm_to_rad_per_s`].
///
/// # Errors
///
/// [`WindTurbineError::BadParameter`] if `omega` is negative or
/// non-finite.
pub fn rad_per_s_to_rpm(omega: f64) -> Result<f64, WindTurbineError> {
    if !omega.is_finite() || omega < 0.0 {
        return Err(WindTurbineError::BadParameter {
            name: "omega",
            reason: "must be >= 0".to_string(),
        });
    }
    Ok(omega * 60.0 / std::f64::consts::TAU)
}

/// Tip-speed ratio `lambda = omega * R / v`.
///
/// `omega` is the rotor angular velocity (rad/s), `radius` the tip
/// radius (m), and `wind_speed` the free-stream speed (m/s).
///
/// # Errors
///
/// [`WindTurbineError::BadParameter`] if `omega` is negative, if
/// `radius` or `wind_speed` is not strictly positive (a zero wind speed
/// would divide by zero), or if any argument is non-finite.
///
/// ```
/// use valenx_windturbine::tsr::tip_speed_ratio;
/// // omega = 2 rad/s, R = 30 m, v = 12 m/s  ->  lambda = 5
/// let lambda = tip_speed_ratio(2.0, 30.0, 12.0).unwrap();
/// assert!((lambda - 5.0).abs() < 1e-12);
/// ```
pub fn tip_speed_ratio(omega: f64, radius: f64, wind_speed: f64) -> Result<f64, WindTurbineError> {
    if !omega.is_finite() || omega < 0.0 {
        return Err(WindTurbineError::BadParameter {
            name: "omega",
            reason: "must be >= 0".to_string(),
        });
    }
    if !radius.is_finite() || radius <= 0.0 {
        return Err(WindTurbineError::BadParameter {
            name: "radius",
            reason: "must be > 0".to_string(),
        });
    }
    if !wind_speed.is_finite() || wind_speed <= 0.0 {
        return Err(WindTurbineError::BadParameter {
            name: "wind_speed",
            reason: "must be > 0".to_string(),
        });
    }
    Ok(omega * radius / wind_speed)
}

/// Tip speed `v_tip = omega * R` (m/s) — the linear speed of a blade
/// tip.
///
/// # Errors
///
/// [`WindTurbineError::BadParameter`] if `omega` is negative, if
/// `radius` is not strictly positive, or if either is non-finite.
pub fn tip_speed(omega: f64, radius: f64) -> Result<f64, WindTurbineError> {
    if !omega.is_finite() || omega < 0.0 {
        return Err(WindTurbineError::BadParameter {
            name: "omega",
            reason: "must be >= 0".to_string(),
        });
    }
    if !radius.is_finite() || radius <= 0.0 {
        return Err(WindTurbineError::BadParameter {
            name: "radius",
            reason: "must be > 0".to_string(),
        });
    }
    Ok(omega * radius)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn tip_speed_ratio_known_value() {
        // omega=2, R=30, v=12 -> lambda = 60/12 = 5.
        let lambda = tip_speed_ratio(2.0, 30.0, 12.0).unwrap();
        assert!((lambda - 5.0).abs() < EPS);
    }

    #[test]
    fn tip_speed_ratio_inversely_proportional_to_wind() {
        let fast_wind = tip_speed_ratio(3.0, 40.0, 24.0).unwrap();
        let slow_wind = tip_speed_ratio(3.0, 40.0, 12.0).unwrap();
        // Halving the wind speed doubles lambda.
        assert!((slow_wind / fast_wind - 2.0).abs() < EPS);
    }

    #[test]
    fn tip_speed_ratio_proportional_to_omega_and_radius() {
        let base = tip_speed_ratio(2.0, 30.0, 10.0).unwrap();
        let dbl_omega = tip_speed_ratio(4.0, 30.0, 10.0).unwrap();
        let dbl_radius = tip_speed_ratio(2.0, 60.0, 10.0).unwrap();
        assert!((dbl_omega / base - 2.0).abs() < EPS);
        assert!((dbl_radius / base - 2.0).abs() < EPS);
    }

    #[test]
    fn tip_speed_ratio_zero_omega_is_zero() {
        let lambda = tip_speed_ratio(0.0, 30.0, 12.0).unwrap();
        assert!(lambda.abs() < EPS);
    }

    #[test]
    fn tip_speed_ratio_rejects_zero_wind() {
        // Would divide by zero.
        let e = tip_speed_ratio(2.0, 30.0, 0.0).unwrap_err();
        assert_eq!(e.code(), "windturbine.bad_parameter");
    }

    #[test]
    fn tip_speed_ratio_rejects_bad_inputs() {
        assert!(tip_speed_ratio(-1.0, 30.0, 12.0).is_err()); // omega < 0
        assert!(tip_speed_ratio(2.0, 0.0, 12.0).is_err()); // radius 0
        assert!(tip_speed_ratio(2.0, -5.0, 12.0).is_err()); // radius < 0
        assert!(tip_speed_ratio(f64::NAN, 30.0, 12.0).is_err());
    }

    #[test]
    fn tip_speed_known_value() {
        // omega=2 rad/s, R=30 m -> 60 m/s.
        let v = tip_speed(2.0, 30.0).unwrap();
        assert!((v - 60.0).abs() < EPS);
    }

    #[test]
    fn lambda_equals_tip_speed_over_wind() {
        // Cross-check: lambda = v_tip / v.
        let omega = 2.5;
        let r = 45.0;
        let wind = 11.0;
        let lambda = tip_speed_ratio(omega, r, wind).unwrap();
        let vtip = tip_speed(omega, r).unwrap();
        assert!((lambda - vtip / wind).abs() < EPS);
    }

    #[test]
    fn rpm_to_rad_per_s_known_value() {
        // 60 rpm = 1 rev/s = 2*pi rad/s.
        let w = rpm_to_rad_per_s(60.0).unwrap();
        assert!((w - std::f64::consts::TAU).abs() < EPS);
        // 0 rpm = 0 rad/s.
        assert!(rpm_to_rad_per_s(0.0).unwrap().abs() < EPS);
    }

    #[test]
    fn rad_per_s_round_trips_with_rpm() {
        for &rpm in &[1.0, 15.0, 60.0, 1500.0] {
            let w = rpm_to_rad_per_s(rpm).unwrap();
            let back = rad_per_s_to_rpm(w).unwrap();
            assert!((back - rpm).abs() < 1e-9, "round-trip failed for {rpm}");
        }
    }

    #[test]
    fn rpm_converters_reject_negative() {
        assert!(rpm_to_rad_per_s(-1.0).is_err());
        assert!(rad_per_s_to_rpm(-1.0).is_err());
    }
}
