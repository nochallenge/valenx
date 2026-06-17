//! Actuator-disc power extraction: available power, the Betz limit, and
//! the power actually captured for a given power coefficient.
//!
//! # Model
//!
//! A wind turbine sweeps a disc of area `A` (m²). Air of density `rho`
//! (kg/m³) flows through it at free-stream speed `v` (m/s). The
//! kinetic-energy flux through that disc — the *available* power — is
//!
//! ```text
//! P_avail = 1/2 * rho * A * v^3
//! ```
//!
//! the classic cube-in-wind-speed law. Betz's momentum theory shows an
//! ideal, loss-free actuator disc can convert at most a fraction
//! `Cp_max = 16/27 ~ 0.593` of that flux into shaft power, because
//! slowing the wind too much chokes the mass flow. The captured power is
//!
//! ```text
//! P = 1/2 * rho * A * v^3 * Cp,   0 <= Cp <= 16/27
//! ```

use crate::error::WindTurbineError;

/// The Betz limit, `16 / 27 ~ 0.5926`.
///
/// The maximum fraction of the wind's kinetic-energy flux that an ideal
/// actuator disc can extract, from one-dimensional momentum theory.
pub const BETZ_LIMIT: f64 = 16.0 / 27.0;

/// Standard sea-level air density, `1.225` kg/m³ (ISA, 15 °C).
///
/// A convenient default for [`available_power`]; real density varies
/// with altitude, temperature, and humidity.
pub const AIR_DENSITY_SEA_LEVEL: f64 = 1.225;

/// Swept disc area `A = pi * r^2` (m²) for a rotor of radius `radius`
/// (m).
///
/// # Errors
///
/// [`WindTurbineError::BadParameter`] if `radius` is not strictly
/// positive or is non-finite.
///
/// ```
/// use valenx_windturbine::power::swept_area;
/// let a = swept_area(2.0).unwrap();
/// assert!((a - std::f64::consts::PI * 4.0).abs() < 1e-12);
/// ```
pub fn swept_area(radius: f64) -> Result<f64, WindTurbineError> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(WindTurbineError::BadParameter {
            name: "radius",
            reason: "must be > 0".to_string(),
        });
    }
    Ok(std::f64::consts::PI * radius * radius)
}

/// Available wind power `P = 1/2 * rho * A * v^3` (W) through a disc of
/// area `area` (m²) for air density `air_density` (kg/m³) and wind speed
/// `wind_speed` (m/s).
///
/// This is the *total* kinetic-energy flux; the extractable shaft power
/// is at most [`BETZ_LIMIT`] times this. See [`extracted_power`].
///
/// # Errors
///
/// [`WindTurbineError::BadParameter`] if `air_density` or `area` is not
/// strictly positive, if `wind_speed` is negative, or if any argument is
/// non-finite.
///
/// ```
/// use valenx_windturbine::power::available_power;
/// // 1/2 * 1.0 * 2.0 * 3^3 = 27 W
/// let p = available_power(1.0, 2.0, 3.0).unwrap();
/// assert!((p - 27.0).abs() < 1e-12);
/// ```
pub fn available_power(
    air_density: f64,
    area: f64,
    wind_speed: f64,
) -> Result<f64, WindTurbineError> {
    if !air_density.is_finite() || air_density <= 0.0 {
        return Err(WindTurbineError::BadParameter {
            name: "air_density",
            reason: "must be > 0".to_string(),
        });
    }
    if !area.is_finite() || area <= 0.0 {
        return Err(WindTurbineError::BadParameter {
            name: "area",
            reason: "must be > 0".to_string(),
        });
    }
    if !wind_speed.is_finite() || wind_speed < 0.0 {
        return Err(WindTurbineError::BadParameter {
            name: "wind_speed",
            reason: "must be >= 0".to_string(),
        });
    }
    Ok(0.5 * air_density * area * wind_speed * wind_speed * wind_speed)
}

/// Validate that a power coefficient `cp` is a physical value in
/// `[0, 16/27]`.
///
/// # Errors
///
/// - [`WindTurbineError::BadParameter`] if `cp` is negative or
///   non-finite.
/// - [`WindTurbineError::AboveBetz`] if `cp` exceeds [`BETZ_LIMIT`].
pub fn validate_cp(cp: f64) -> Result<(), WindTurbineError> {
    if !cp.is_finite() || cp < 0.0 {
        return Err(WindTurbineError::BadParameter {
            name: "cp",
            reason: "must be >= 0".to_string(),
        });
    }
    // Allow a tiny epsilon so a Cp of exactly 16/27 (the Betz limit
    // itself) is accepted despite floating-point round-off.
    if cp > BETZ_LIMIT + 1e-12 {
        return Err(WindTurbineError::AboveBetz {
            cp,
            betz: BETZ_LIMIT,
        });
    }
    Ok(())
}

/// Extracted shaft power `P = 1/2 * rho * A * v^3 * Cp` (W).
///
/// The available power ([`available_power`]) scaled by the power
/// coefficient `cp`. Because `cp` is validated against the Betz limit,
/// the result never exceeds `BETZ_LIMIT * available_power(..)`.
///
/// # Errors
///
/// Propagates the errors of [`available_power`] for the physical inputs,
/// and of [`validate_cp`] for `cp` (including
/// [`WindTurbineError::AboveBetz`]).
///
/// ```
/// use valenx_windturbine::power::{extracted_power, BETZ_LIMIT};
/// // At the Betz limit, extracted = 16/27 of the 27 W available.
/// let p = extracted_power(1.0, 2.0, 3.0, BETZ_LIMIT).unwrap();
/// assert!((p - 27.0 * BETZ_LIMIT).abs() < 1e-12);
/// ```
pub fn extracted_power(
    air_density: f64,
    area: f64,
    wind_speed: f64,
    cp: f64,
) -> Result<f64, WindTurbineError> {
    validate_cp(cp)?;
    let avail = available_power(air_density, area, wind_speed)?;
    Ok(avail * cp)
}

/// Maximum (Betz-limited) extractable power
/// `P = 1/2 * rho * A * v^3 * 16/27` (W).
///
/// A convenience wrapper for `extracted_power(.., BETZ_LIMIT)`.
///
/// # Errors
///
/// Propagates the errors of [`available_power`].
pub fn betz_power(air_density: f64, area: f64, wind_speed: f64) -> Result<f64, WindTurbineError> {
    let avail = available_power(air_density, area, wind_speed)?;
    Ok(avail * BETZ_LIMIT)
}

/// The rotor radius `R` (m) whose swept disc captures a target shaft power
/// at a design wind speed and power coefficient — the turbine-sizing
/// inverse of [`extracted_power`] (with [`swept_area`]).
///
/// Inverting `P = 1/2 * rho * (pi R^2) * v^3 * Cp` for the radius gives
///
/// ```text
/// R = sqrt( 2 P / (rho * pi * v^3 * Cp) ).
/// ```
///
/// Building a [`swept_area`] at this radius and feeding it to
/// [`extracted_power`] (same `rho`, `v`, `Cp`) reproduces `target_power`.
///
/// # Errors
///
/// [`WindTurbineError::AboveBetz`] if `cp` exceeds the Betz limit;
/// [`WindTurbineError::BadParameter`] if `target_power`, `air_density`,
/// `wind_speed`, or `cp` is not strictly positive and finite (a zero `Cp`
/// or wind speed would demand an infinite rotor).
pub fn rotor_radius_for_power(
    target_power: f64,
    air_density: f64,
    wind_speed: f64,
    cp: f64,
) -> Result<f64, WindTurbineError> {
    validate_cp(cp)?;
    for (name, value) in [
        ("target_power", target_power),
        ("air_density", air_density),
        ("wind_speed", wind_speed),
        ("cp", cp),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(WindTurbineError::BadParameter {
                name,
                reason: "must be > 0".to_string(),
            });
        }
    }
    let area = 2.0 * target_power / (air_density * wind_speed.powi(3) * cp);
    Ok((area / std::f64::consts::PI).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    const EPS: f64 = 1e-9;

    #[test]
    fn betz_limit_is_sixteen_over_twentyseven() {
        // The headline ground-truth number: 16/27 ~ 0.5926.
        assert!((BETZ_LIMIT - 16.0 / 27.0).abs() < EPS);
        assert!((BETZ_LIMIT - 0.592_592_592_592_592_6).abs() < 1e-12);
    }

    #[test]
    fn swept_area_matches_pi_r_squared() {
        let a = swept_area(2.0).unwrap();
        assert!((a - std::f64::consts::PI * 4.0).abs() < EPS);
        // Doubling radius quadruples area.
        let a2 = swept_area(4.0).unwrap();
        assert!((a2 / a - 4.0).abs() < EPS);
    }

    #[test]
    fn swept_area_rejects_non_positive() {
        for r in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let e = swept_area(r).unwrap_err();
            assert_eq!(e.code(), "windturbine.bad_parameter");
            assert_eq!(e.category(), ErrorCategory::Input);
        }
    }

    #[test]
    fn available_power_known_value() {
        // 1/2 * 1.0 * 2.0 * 3^3 = 27 W (hand-computed).
        let p = available_power(1.0, 2.0, 3.0).unwrap();
        assert!((p - 27.0).abs() < EPS);
    }

    #[test]
    fn available_power_scales_with_cube_of_speed() {
        let base = available_power(1.225, 10.0, 4.0).unwrap();
        let doubled = available_power(1.225, 10.0, 8.0).unwrap();
        // v -> 2v  =>  P -> 8 P.
        assert!((doubled / base - 8.0).abs() < 1e-6);
        let tripled = available_power(1.225, 10.0, 12.0).unwrap();
        assert!((tripled / base - 27.0).abs() < 1e-6);
    }

    #[test]
    fn available_power_scales_linearly_with_area() {
        let a1 = available_power(1.225, 5.0, 7.0).unwrap();
        let a2 = available_power(1.225, 15.0, 7.0).unwrap();
        // 3x area => 3x power.
        assert!((a2 / a1 - 3.0).abs() < 1e-9);
    }

    #[test]
    fn available_power_scales_linearly_with_density() {
        let d1 = available_power(1.0, 5.0, 7.0).unwrap();
        let d2 = available_power(2.5, 5.0, 7.0).unwrap();
        assert!((d2 / d1 - 2.5).abs() < 1e-9);
    }

    #[test]
    fn available_power_zero_at_zero_wind() {
        let p = available_power(1.225, 10.0, 0.0).unwrap();
        assert!(p.abs() < EPS);
    }

    #[test]
    fn available_power_rejects_bad_inputs() {
        assert!(available_power(0.0, 1.0, 1.0).is_err()); // density 0
        assert!(available_power(-1.0, 1.0, 1.0).is_err()); // density < 0
        assert!(available_power(1.0, 0.0, 1.0).is_err()); // area 0
        assert!(available_power(1.0, 1.0, -1.0).is_err()); // speed < 0
        assert!(available_power(f64::NAN, 1.0, 1.0).is_err());
    }

    #[test]
    fn extracted_never_exceeds_betz_never_exceeds_available() {
        let rho = 1.225;
        let area = 30.0;
        for &v in &[1.0, 5.0, 9.0, 13.0, 25.0] {
            let avail = available_power(rho, area, v).unwrap();
            let betz = betz_power(rho, area, v).unwrap();
            // A realistic operating Cp ~ 0.45 < Betz.
            let ext = extracted_power(rho, area, v, 0.45).unwrap();
            assert!(ext <= betz + EPS, "extracted {ext} > betz {betz}");
            assert!(betz <= avail + EPS, "betz {betz} > available {avail}");
            // Betz is exactly 16/27 of available.
            assert!((betz / avail - BETZ_LIMIT).abs() < 1e-9);
        }
    }

    #[test]
    fn extracted_at_betz_equals_betz_power() {
        let a = extracted_power(1.0, 2.0, 3.0, BETZ_LIMIT).unwrap();
        let b = betz_power(1.0, 2.0, 3.0).unwrap();
        assert!((a - b).abs() < EPS);
        // == 27 * 16/27 = 16 exactly.
        assert!((a - 16.0).abs() < 1e-9);
    }

    #[test]
    fn extracted_scales_linearly_with_cp() {
        let p_low = extracted_power(1.225, 20.0, 8.0, 0.2).unwrap();
        let p_high = extracted_power(1.225, 20.0, 8.0, 0.4).unwrap();
        assert!((p_high / p_low - 2.0).abs() < 1e-9);
    }

    #[test]
    fn cp_zero_gives_zero_power() {
        let p = extracted_power(1.225, 20.0, 8.0, 0.0).unwrap();
        assert!(p.abs() < EPS);
    }

    #[test]
    fn validate_cp_accepts_exactly_betz() {
        assert!(validate_cp(BETZ_LIMIT).is_ok());
        assert!(validate_cp(0.0).is_ok());
        assert!(validate_cp(0.45).is_ok());
    }

    #[test]
    fn validate_cp_rejects_above_betz() {
        // 0.6 > 16/27 ~ 0.593 — physically impossible.
        let e = validate_cp(0.6).unwrap_err();
        assert_eq!(e.code(), "windturbine.above_betz");
        assert_eq!(e.category(), ErrorCategory::Algorithm);
        // 1.0 (more than the wind carries) likewise.
        assert!(validate_cp(1.0).is_err());
    }

    #[test]
    fn validate_cp_rejects_negative() {
        let e = validate_cp(-0.01).unwrap_err();
        assert_eq!(e.code(), "windturbine.bad_parameter");
    }

    #[test]
    fn extracted_propagates_above_betz() {
        let e = extracted_power(1.225, 20.0, 8.0, 0.7).unwrap_err();
        assert_eq!(e.code(), "windturbine.above_betz");
    }

    #[test]
    fn rotor_radius_for_power_round_trips_through_extracted() {
        // Size R for a target power, build the disc, recover the power.
        let (rho, v, cp) = (1.225, 12.0, 0.45);
        let target = 2.0e6; // 2 MW
        let r = rotor_radius_for_power(target, rho, v, cp).unwrap();
        let area = swept_area(r).unwrap();
        let p = extracted_power(rho, area, v, cp).unwrap();
        assert!((p - target).abs() < 1e-6 * target, "r={r} p={p}");
    }

    #[test]
    fn rotor_radius_for_power_hand_value() {
        // P = 0.5 * 1 * A * 3^3 * 0.5 with A = 2 gives P = 13.5 W, so the
        // radius is sqrt(2/pi).
        let r = rotor_radius_for_power(13.5, 1.0, 3.0, 0.5).unwrap();
        let expected = (2.0_f64 / std::f64::consts::PI).sqrt();
        assert!((r - expected).abs() < 1e-9, "r={r}");
    }

    #[test]
    fn rotor_radius_scales_and_rejects_bad() {
        let base = rotor_radius_for_power(1.0e6, 1.225, 10.0, 0.4).unwrap();
        // R ~ sqrt(P): 4x power -> 2x radius.
        let quad = rotor_radius_for_power(4.0e6, 1.225, 10.0, 0.4).unwrap();
        assert!((quad / base - 2.0).abs() < 1e-9, "base={base} quad={quad}");
        assert!(rotor_radius_for_power(0.0, 1.225, 10.0, 0.4).is_err()); // P
        assert!(rotor_radius_for_power(1.0e6, 1.225, 0.0, 0.4).is_err()); // v
        assert!(rotor_radius_for_power(1.0e6, 1.225, 10.0, 0.0).is_err()); // cp = 0
        assert!(rotor_radius_for_power(1.0e6, 1.225, 10.0, 0.7).is_err()); // > Betz
        assert!(rotor_radius_for_power(f64::NAN, 1.225, 10.0, 0.4).is_err());
    }
}
