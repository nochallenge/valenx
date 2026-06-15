//! Braking kinetic energy and stopping distance.
//!
//! A mass `m` translating at speed `v` carries kinetic energy
//!
//! ```text
//! E = 0.5 * m * v^2
//! ```
//!
//! Ignoring all losses, this is the heat the brakes must absorb to
//! bring the mass to rest. Under a constant deceleration `a > 0` the
//! work–energy theorem (`E = m * a * d`) gives the stopping distance
//! `d = v^2 / (2 a)` and stopping time `t = v / a`.

use crate::error::{check_non_negative, check_positive, BrakeError};

/// Translational kinetic energy `E = 0.5 * m * v^2`, in joules.
///
/// This is the energy the brakes must dissipate to stop the mass
/// (ignoring rotational inertia, rolling resistance, aero drag, etc.).
///
/// # Parameters
/// - `mass_kg` — translating mass `m`, in kilograms (> 0).
/// - `speed_mps` — speed `v`, in metres per second (>= 0).
///
/// # Errors
/// [`BrakeError`] if `mass_kg` is non-finite or non-positive, or
/// `speed_mps` is non-finite or negative.
///
/// # Examples
/// ```
/// use valenx_brake::energy::kinetic_energy;
/// // 1500 kg at 30 m/s -> 0.5*1500*900 = 675 kJ.
/// let e = kinetic_energy(1500.0, 30.0).unwrap();
/// assert!((e - 675_000.0).abs() < 1e-6);
/// ```
pub fn kinetic_energy(mass_kg: f64, speed_mps: f64) -> Result<f64, BrakeError> {
    let m = check_positive("mass_kg", mass_kg)?;
    let v = check_non_negative("speed_mps", speed_mps)?;
    Ok(0.5 * m * v * v)
}

/// Kinetic energy dissipated decelerating from `v_initial` to
/// `v_final`, in joules.
///
/// `ΔE = 0.5 * m * (v_initial^2 - v_final^2)`. The full-stop energy is
/// the `v_final = 0` case, which equals [`kinetic_energy`].
///
/// # Errors
/// [`BrakeError`] if `mass_kg` is non-finite or non-positive, or either
/// speed is non-finite or negative. (No ordering is required; the
/// signed difference is returned, so a speed-up yields a negative
/// value.)
///
/// # Examples
/// ```
/// use valenx_brake::energy::energy_dissipated;
/// // 1000 kg from 20 to 10 m/s -> 0.5*1000*(400-100) = 150 kJ.
/// let e = energy_dissipated(1000.0, 20.0, 10.0).unwrap();
/// assert!((e - 150_000.0).abs() < 1e-6);
/// ```
pub fn energy_dissipated(
    mass_kg: f64,
    v_initial_mps: f64,
    v_final_mps: f64,
) -> Result<f64, BrakeError> {
    let m = check_positive("mass_kg", mass_kg)?;
    let vi = check_non_negative("v_initial_mps", v_initial_mps)?;
    let vf = check_non_negative("v_final_mps", v_final_mps)?;
    Ok(0.5 * m * (vi * vi - vf * vf))
}

/// Stopping distance under a constant deceleration, in metres.
///
/// From `E = m * a * d` with `E = 0.5 * m * v^2`, the mass cancels and
/// `d = v^2 / (2 * a)`.
///
/// # Parameters
/// - `speed_mps` — initial speed `v`, in metres per second (>= 0).
/// - `decel_mps2` — constant deceleration magnitude `a`, in metres per
///   second squared (> 0).
///
/// # Errors
/// [`BrakeError`] if `speed_mps` is non-finite or negative, or
/// `decel_mps2` is non-finite or non-positive.
///
/// # Examples
/// ```
/// use valenx_brake::energy::stopping_distance;
/// // 30 m/s at 7.5 m/s^2 -> 900 / 15 = 60 m.
/// let d = stopping_distance(30.0, 7.5).unwrap();
/// assert!((d - 60.0).abs() < 1e-9);
/// ```
pub fn stopping_distance(speed_mps: f64, decel_mps2: f64) -> Result<f64, BrakeError> {
    let v = check_non_negative("speed_mps", speed_mps)?;
    let a = check_positive("decel_mps2", decel_mps2)?;
    Ok(v * v / (2.0 * a))
}

/// Stopping time under a constant deceleration, in seconds.
///
/// `t = v / a`.
///
/// # Errors
/// [`BrakeError`] if `speed_mps` is non-finite or negative, or
/// `decel_mps2` is non-finite or non-positive.
///
/// # Examples
/// ```
/// use valenx_brake::energy::stopping_time;
/// // 30 m/s at 7.5 m/s^2 -> 4 s.
/// let t = stopping_time(30.0, 7.5).unwrap();
/// assert!((t - 4.0).abs() < 1e-9);
/// ```
pub fn stopping_time(speed_mps: f64, decel_mps2: f64) -> Result<f64, BrakeError> {
    let v = check_non_negative("speed_mps", speed_mps)?;
    let a = check_positive("decel_mps2", decel_mps2)?;
    Ok(v / a)
}

/// Average braking power needed to stop in a given time, in watts.
///
/// `P_avg = E / t = (0.5 * m * v^2) / t`. This is the mean rate the
/// brakes must dissipate the kinetic energy over the stop.
///
/// # Errors
/// [`BrakeError`] if `mass_kg` is non-finite or non-positive,
/// `speed_mps` is non-finite or negative, or `time_s` is non-finite or
/// non-positive.
///
/// # Examples
/// ```
/// use valenx_brake::energy::average_braking_power;
/// // 675 kJ dissipated in 4 s -> 168.75 kW.
/// let p = average_braking_power(1500.0, 30.0, 4.0).unwrap();
/// assert!((p - 168_750.0).abs() < 1e-6);
/// ```
pub fn average_braking_power(mass_kg: f64, speed_mps: f64, time_s: f64) -> Result<f64, BrakeError> {
    let e = kinetic_energy(mass_kg, speed_mps)?;
    let t = check_positive("time_s", time_s)?;
    Ok(e / t)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    #[test]
    fn kinetic_energy_matches_half_m_v_squared() {
        // 0.5 * 1500 * 30^2 = 675_000 J.
        let e = kinetic_energy(1500.0, 30.0).unwrap();
        assert!((e - 675_000.0).abs() < EPS, "got {e}");
    }

    #[test]
    fn kinetic_energy_zero_at_rest() {
        let e = kinetic_energy(1500.0, 0.0).unwrap();
        assert!(e.abs() < EPS, "got {e}");
    }

    #[test]
    fn energy_scales_with_speed_squared() {
        // Doubling v quadruples E (quadratic dependence).
        let base = kinetic_energy(1000.0, 10.0).unwrap();
        let doubled = kinetic_energy(1000.0, 20.0).unwrap();
        assert!(
            (doubled - 4.0 * base).abs() < EPS,
            "base {base} doubled {doubled}"
        );
    }

    #[test]
    fn energy_scales_linearly_with_mass() {
        // Doubling m doubles E (linear in mass).
        let base = kinetic_energy(1000.0, 25.0).unwrap();
        let doubled = kinetic_energy(2000.0, 25.0).unwrap();
        assert!(
            (doubled - 2.0 * base).abs() < EPS,
            "base {base} doubled {doubled}"
        );
    }

    #[test]
    fn dissipated_difference_of_squares() {
        // 0.5 * 1000 * (20^2 - 10^2) = 150_000 J.
        let e = energy_dissipated(1000.0, 20.0, 10.0).unwrap();
        assert!((e - 150_000.0).abs() < EPS, "got {e}");
    }

    #[test]
    fn dissipated_full_stop_equals_kinetic_energy() {
        // v_final = 0 must equal the standalone kinetic energy.
        let m = 1234.0;
        let v = 27.5;
        let full = energy_dissipated(m, v, 0.0).unwrap();
        let ke = kinetic_energy(m, v).unwrap();
        assert!((full - ke).abs() < EPS, "full {full} ke {ke}");
    }

    #[test]
    fn stopping_distance_closed_form() {
        // d = v^2 / (2a) = 900 / 15 = 60 m.
        let d = stopping_distance(30.0, 7.5).unwrap();
        assert!((d - 60.0).abs() < 1e-9, "got {d}");
    }

    #[test]
    fn stopping_distance_consistent_with_work_energy() {
        // E = m*a*d must reproduce 0.5*m*v^2 for the same stop.
        let m = 1500.0;
        let v = 30.0;
        let a = 7.5;
        let d = stopping_distance(v, a).unwrap();
        let work = m * a * d;
        let ke = kinetic_energy(m, v).unwrap();
        assert!((work - ke).abs() < 1e-3, "work {work} ke {ke}");
    }

    #[test]
    fn stopping_time_closed_form() {
        // t = v / a = 30 / 7.5 = 4 s.
        let t = stopping_time(30.0, 7.5).unwrap();
        assert!((t - 4.0).abs() < 1e-9, "got {t}");
    }

    #[test]
    fn average_power_is_energy_over_time() {
        // 675 kJ over 4 s = 168.75 kW.
        let p = average_braking_power(1500.0, 30.0, 4.0).unwrap();
        assert!((p - 168_750.0).abs() < EPS, "got {p}");
        // Cross-check: E / t.
        let e = kinetic_energy(1500.0, 30.0).unwrap();
        assert!((p - e / 4.0).abs() < EPS, "p {p} e {e}");
    }

    #[test]
    fn rejects_bad_inputs() {
        assert_eq!(
            kinetic_energy(0.0, 30.0).unwrap_err().code(),
            "brake.non_positive"
        );
        assert_eq!(
            kinetic_energy(1500.0, -1.0).unwrap_err().code(),
            "brake.negative"
        );
        assert_eq!(
            stopping_distance(30.0, 0.0).unwrap_err().code(),
            "brake.non_positive"
        );
        assert_eq!(
            stopping_time(30.0, -2.0).unwrap_err().code(),
            "brake.non_positive"
        );
        assert_eq!(
            average_braking_power(1500.0, 30.0, 0.0).unwrap_err().code(),
            "brake.non_positive"
        );
        assert_eq!(
            kinetic_energy(f64::NAN, 30.0).unwrap_err().code(),
            "brake.not_finite"
        );
    }
}
