//! Rotational energy, usable-energy, and governor speed-fluctuation math.
//!
//! These free functions express the core energy relations independently
//! of any rotor geometry, taking a moment of inertia `I` (kg.m^2) and one
//! or more angular speeds `omega` (rad/s) directly. The [`crate::Flywheel`]
//! type wraps them with a concrete rotor.
//!
//! Unit helpers ([`rpm_to_rad_s`], [`rad_s_to_rpm`]) bridge the
//! revolutions-per-minute speeds quoted on real machines and the radians
//! per second the formulas require.

use std::f64::consts::PI;

use crate::error::FlywheelError;

/// Rotational kinetic energy `E = 1/2 I omega^2`, in joules.
///
/// `inertia` is the moment of inertia about the spin axis (kg.m^2) and
/// `omega` the angular speed (rad/s).
///
/// # Errors
///
/// Returns [`FlywheelError::InvalidParameter`] if `inertia` is not
/// strictly positive or `omega` is negative / non-finite (a speed of
/// exactly zero is allowed and yields zero energy).
pub fn kinetic_energy(inertia: f64, omega: f64) -> Result<f64, FlywheelError> {
    let inertia = FlywheelError::require_positive("inertia", inertia)?;
    let omega = FlywheelError::require_non_negative("omega", omega)?;
    Ok(0.5 * inertia * omega * omega)
}

/// Usable (extractable) energy as the rotor decelerates from
/// `omega_max` down to `omega_min`:
///
/// `dE = 1/2 I (omega_max^2 - omega_min^2)`, in joules.
///
/// This is simply the difference of the kinetic energies at the two
/// speeds and is the energy a flywheel can deliver while staying within
/// an allowed speed band.
///
/// # Errors
///
/// Returns [`FlywheelError::InvalidParameter`] for a non-positive
/// inertia or a negative / non-finite speed, and
/// [`FlywheelError::Inconsistent`] if `omega_min > omega_max`.
pub fn usable_energy(inertia: f64, omega_min: f64, omega_max: f64) -> Result<f64, FlywheelError> {
    let inertia = FlywheelError::require_positive("inertia", inertia)?;
    let omega_min = FlywheelError::require_non_negative("omega_min", omega_min)?;
    let omega_max = FlywheelError::require_non_negative("omega_max", omega_max)?;
    if omega_min > omega_max {
        return Err(FlywheelError::Inconsistent(
            "omega_min must be <= omega_max",
        ));
    }
    Ok(0.5 * inertia * (omega_max * omega_max - omega_min * omega_min))
}

/// Coefficient of fluctuation `Cs = (omega_max - omega_min) / omega_avg`,
/// where `omega_avg = (omega_max + omega_min) / 2`.
///
/// A dimensionless measure of how much the speed of a flywheel-smoothed
/// shaft swings about its mean over one cycle; smaller is steadier.
///
/// # Errors
///
/// Returns [`FlywheelError::InvalidParameter`] for a negative /
/// non-finite speed, [`FlywheelError::Inconsistent`] if
/// `omega_min > omega_max`, and [`FlywheelError::InvalidParameter`] if
/// the mean speed is zero (both speeds zero), for which the coefficient
/// is undefined.
pub fn coefficient_of_fluctuation(omega_min: f64, omega_max: f64) -> Result<f64, FlywheelError> {
    let omega_min = FlywheelError::require_non_negative("omega_min", omega_min)?;
    let omega_max = FlywheelError::require_non_negative("omega_max", omega_max)?;
    if omega_min > omega_max {
        return Err(FlywheelError::Inconsistent(
            "omega_min must be <= omega_max",
        ));
    }
    let omega_avg = 0.5 * (omega_max + omega_min);
    if omega_avg <= 0.0 {
        return Err(FlywheelError::invalid(
            "omega_avg",
            "mean speed must be > 0",
            omega_avg,
        ));
    }
    Ok((omega_max - omega_min) / omega_avg)
}

/// Energy fluctuation in the governor-sizing form `dE = I omega_avg^2 Cs`,
/// in joules.
///
/// This is algebraically identical to [`usable_energy`] between the band
/// edges: with `Cs = (omega_max - omega_min)/omega_avg` and
/// `omega_avg = (omega_max + omega_min)/2`,
/// `I omega_avg^2 Cs = 1/2 I (omega_max^2 - omega_min^2)`. It is the form
/// used when sizing an engine flywheel from a target coefficient of
/// fluctuation and a mean operating speed.
///
/// # Errors
///
/// Returns [`FlywheelError::InvalidParameter`] for a non-positive
/// inertia, a non-positive mean speed, or a negative / non-finite
/// coefficient of fluctuation.
pub fn energy_fluctuation(inertia: f64, omega_avg: f64, cs: f64) -> Result<f64, FlywheelError> {
    let inertia = FlywheelError::require_positive("inertia", inertia)?;
    let omega_avg = FlywheelError::require_positive("omega_avg", omega_avg)?;
    let cs = FlywheelError::require_non_negative("cs", cs)?;
    Ok(inertia * omega_avg * omega_avg * cs)
}

/// Solve the governor-sizing relation for the moment of inertia required
/// to supply a target energy fluctuation `dE` at mean speed `omega_avg`
/// with coefficient of fluctuation `Cs`:
///
/// `I = dE / (omega_avg^2 Cs)`, in kg.m^2.
///
/// The inverse of [`energy_fluctuation`]; answers "how big a flywheel do
/// I need?".
///
/// # Errors
///
/// Returns [`FlywheelError::InvalidParameter`] for a negative /
/// non-finite `energy`, or a non-positive `omega_avg` or `cs` (each must
/// be strictly positive to divide by).
pub fn flywheel_inertia_for_energy(
    energy: f64,
    omega_avg: f64,
    cs: f64,
) -> Result<f64, FlywheelError> {
    let energy = FlywheelError::require_non_negative("energy", energy)?;
    let omega_avg = FlywheelError::require_positive("omega_avg", omega_avg)?;
    let cs = FlywheelError::require_positive("cs", cs)?;
    Ok(energy / (omega_avg * omega_avg * cs))
}

/// Convert revolutions per minute to radians per second:
/// `omega = rpm * 2 PI / 60`.
///
/// # Errors
///
/// Returns [`FlywheelError::InvalidParameter`] if `rpm` is negative or
/// non-finite.
pub fn rpm_to_rad_s(rpm: f64) -> Result<f64, FlywheelError> {
    let rpm = FlywheelError::require_non_negative("rpm", rpm)?;
    Ok(rpm * 2.0 * PI / 60.0)
}

/// Convert radians per second to revolutions per minute:
/// `rpm = omega * 60 / (2 PI)`.
///
/// # Errors
///
/// Returns [`FlywheelError::InvalidParameter`] if `omega` is negative or
/// non-finite.
pub fn rad_s_to_rpm(omega: f64) -> Result<f64, FlywheelError> {
    let omega = FlywheelError::require_non_negative("omega", omega)?;
    Ok(omega * 60.0 / (2.0 * PI))
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn kinetic_energy_matches_half_i_omega_squared() {
        // I = 2, omega = 3 -> E = 0.5 * 2 * 9 = 9
        let e = kinetic_energy(2.0, 3.0).unwrap();
        assert!((e - 9.0).abs() < EPS);
    }

    #[test]
    fn kinetic_energy_is_zero_at_rest() {
        let e = kinetic_energy(5.0, 0.0).unwrap();
        assert!(e.abs() < EPS);
    }

    #[test]
    fn doubling_speed_quadruples_energy() {
        // E ~ omega^2, so E(2 omega) / E(omega) = 4 exactly.
        let inertia = 3.7;
        let omega = 11.0;
        let e1 = kinetic_energy(inertia, omega).unwrap();
        let e2 = kinetic_energy(inertia, 2.0 * omega).unwrap();
        assert!((e2 - 4.0 * e1).abs() < 1e-6);
        assert!((e2 / e1 - 4.0).abs() < EPS);
    }

    #[test]
    fn tripling_speed_gives_nine_times_energy() {
        let inertia = 1.25;
        let omega = 4.0;
        let e1 = kinetic_energy(inertia, omega).unwrap();
        let e3 = kinetic_energy(inertia, 3.0 * omega).unwrap();
        assert!((e3 / e1 - 9.0).abs() < EPS);
    }

    #[test]
    fn usable_energy_is_difference_of_kinetic_energies() {
        let inertia = 2.0;
        let lo = 10.0;
        let hi = 20.0;
        let du = usable_energy(inertia, lo, hi).unwrap();
        let manual = kinetic_energy(inertia, hi).unwrap() - kinetic_energy(inertia, lo).unwrap();
        assert!((du - manual).abs() < 1e-6);
        // 0.5 * 2 * (400 - 100) = 300
        assert!((du - 300.0).abs() < 1e-6);
    }

    #[test]
    fn usable_energy_zero_when_band_collapses() {
        let du = usable_energy(4.0, 15.0, 15.0).unwrap();
        assert!(du.abs() < EPS);
    }

    #[test]
    fn coefficient_of_fluctuation_matches_definition() {
        // omega_max = 105, omega_min = 95 -> avg = 100, Cs = 10/100 = 0.1
        let cs = coefficient_of_fluctuation(95.0, 105.0).unwrap();
        assert!((cs - 0.1).abs() < EPS);
    }

    #[test]
    fn energy_fluctuation_matches_usable_energy_form() {
        // The governor form dE = I omega_avg^2 Cs must equal the
        // difference-of-energies form 0.5 I (wmax^2 - wmin^2).
        let inertia = 3.3;
        let wmin = 90.0;
        let wmax = 110.0;
        let omega_avg = 0.5 * (wmin + wmax);
        let cs = coefficient_of_fluctuation(wmin, wmax).unwrap();

        let governor = energy_fluctuation(inertia, omega_avg, cs).unwrap();
        let direct = usable_energy(inertia, wmin, wmax).unwrap();
        assert!((governor - direct).abs() < 1e-6);
    }

    #[test]
    fn inertia_for_energy_inverts_energy_fluctuation() {
        let omega_avg = 157.08; // ~1500 rpm
        let cs = 0.02;
        let target = 5_000.0; // joules
        let inertia = flywheel_inertia_for_energy(target, omega_avg, cs).unwrap();
        let back = energy_fluctuation(inertia, omega_avg, cs).unwrap();
        assert!((back - target).abs() < 1e-6);
    }

    #[test]
    fn rpm_round_trips_through_rad_s() {
        let rpm = 3000.0;
        let omega = rpm_to_rad_s(rpm).unwrap();
        // 3000 rpm = 50 rev/s = 100 PI rad/s.
        assert!((omega - 100.0 * PI).abs() < 1e-9);
        let back = rad_s_to_rpm(omega).unwrap();
        assert!((back - rpm).abs() < 1e-9);
    }

    #[test]
    fn energy_rejects_bad_inputs() {
        assert!(kinetic_energy(0.0, 1.0).is_err());
        assert!(kinetic_energy(1.0, -1.0).is_err());
        assert!(usable_energy(1.0, 20.0, 10.0).is_err());
        assert!(coefficient_of_fluctuation(0.0, 0.0).is_err());
        assert!(flywheel_inertia_for_energy(1.0, 0.0, 0.1).is_err());
        assert!(flywheel_inertia_for_energy(1.0, 10.0, 0.0).is_err());
    }
}
