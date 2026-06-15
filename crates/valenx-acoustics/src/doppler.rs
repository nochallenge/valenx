//! Speed of sound and the classical Doppler frequency shift.
//!
//! ## Model
//!
//! The speed of sound in dry air is taken from the standard first-order
//! temperature relation
//!
//! ```text
//! c(T) = 331.3 * sqrt(1 + T / 273.15)        [m/s]
//! ```
//!
//! where `T` is the air temperature in degrees Celsius and `331.3` m/s is
//! the 0 degrees Celsius value. At 20 degrees Celsius this recovers the
//! familiar `≈ 343` m/s.
//!
//! The classical (non-relativistic) Doppler shift for a source and an
//! observer moving along the line joining them, through a still medium of
//! sound speed `c`, is
//!
//! ```text
//! f_obs = f_src * (c + v_o) / (c - v_s)
//! ```
//!
//! with the **sign convention** used throughout this module:
//!
//! - `v_o` (observer velocity) is **positive when the observer moves
//!   toward the source**.
//! - `v_s` (source velocity) is **positive when the source moves toward
//!   the observer**.
//!
//! So an approaching source (`v_s > 0`) raises the pitch, a receding
//! source (`v_s < 0`) lowers it, and a stationary geometry returns the
//! emitted frequency unchanged.
//!
//! ## Honest scope
//!
//! This is the textbook 1-D closed form. It assumes motion strictly along
//! the source–observer line (no transverse / angle-of-arrival term), a
//! still homogeneous medium, and subsonic source motion (`v_s < c`). The
//! `c - v_s <= 0` sonic / supersonic case is reported as an error rather
//! than returned as a non-physical value. Research/educational grade.

use crate::error::{AcousticsError, Result};

/// Speed of sound at 0 degrees Celsius in dry air, in metres per second —
/// the leading coefficient of the [`speed_of_sound`] relation.
pub const C0: f64 = 331.3;

fn check_velocity(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(AcousticsError::InvalidVelocity { name, value })
    }
}

fn check_speed_of_sound(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(AcousticsError::InvalidSpeedOfSound { name, value })
    }
}

/// Speed of sound in dry air at a temperature `temperature_c` (degrees
/// Celsius), in metres per second.
///
/// `c = 331.3 * sqrt(1 + T / 273.15)`.
///
/// # Errors
///
/// Returns [`AcousticsError::InvalidTemperature`] if `temperature_c` is
/// not finite or is at/below absolute zero (`<= -273.15` degrees Celsius),
/// where the radicand `1 + T/273.15` is non-positive.
///
/// # Examples
///
/// ```
/// use valenx_acoustics::doppler::speed_of_sound;
/// // ~343 m/s at 20 degrees Celsius.
/// assert!((speed_of_sound(20.0).unwrap() - 343.0).abs() < 1.0);
/// ```
pub fn speed_of_sound(temperature_c: f64) -> Result<f64> {
    let radicand = 1.0 + temperature_c / 273.15;
    if !temperature_c.is_finite() || radicand <= 0.0 {
        return Err(AcousticsError::InvalidTemperature {
            name: "temperature_c",
            value: temperature_c,
        });
    }
    Ok(C0 * radicand.sqrt())
}

/// Classical Doppler-shifted frequency observed when a source emitting at
/// `source_hz` and an observer move along the line joining them through a
/// still medium of sound speed `speed_of_sound`.
///
/// `f_obs = source_hz * (c + v_o) / (c - v_s)`, with `observer_velocity`
/// (`v_o`) positive **toward** the source and `source_velocity` (`v_s`)
/// positive **toward** the observer (see the module docs for the full
/// sign convention).
///
/// # Errors
///
/// - [`AcousticsError::InvalidFrequency`] if `source_hz` is negative or
///   non-finite.
/// - [`AcousticsError::InvalidSpeedOfSound`] if `speed_of_sound` is not
///   finite and strictly positive.
/// - [`AcousticsError::InvalidVelocity`] if either velocity is non-finite.
/// - [`AcousticsError::DopplerSingularity`] if `c - v_s <= 0` (the source
///   reaches or exceeds the speed of sound toward the observer), where the
///   subsonic formula has no finite value.
///
/// # Examples
///
/// ```
/// use valenx_acoustics::doppler::doppler_shift;
/// let c = 343.0;
/// // Source approaching at 30 m/s, observer still: pitch rises.
/// let f = doppler_shift(1000.0, 0.0, 30.0, c).unwrap();
/// assert!(f > 1000.0);
/// ```
pub fn doppler_shift(
    source_hz: f64,
    observer_velocity: f64,
    source_velocity: f64,
    speed_of_sound: f64,
) -> Result<f64> {
    if !source_hz.is_finite() || source_hz < 0.0 {
        return Err(AcousticsError::InvalidFrequency {
            name: "source_hz",
            value: source_hz,
        });
    }
    check_speed_of_sound("speed_of_sound", speed_of_sound)?;
    check_velocity("observer_velocity", observer_velocity)?;
    check_velocity("source_velocity", source_velocity)?;

    let denom = speed_of_sound - source_velocity;
    if denom <= 0.0 {
        return Err(AcousticsError::DopplerSingularity {
            source_velocity,
            speed_of_sound,
        });
    }
    Ok(source_hz * (speed_of_sound + observer_velocity) / denom)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    /// c(20 degrees C) is about 343 m/s.
    #[test]
    fn speed_of_sound_at_twenty_c() {
        let c = speed_of_sound(20.0).unwrap();
        assert!((c - 343.0).abs() < 1.0, "got {c} m/s");
    }

    /// c(0 degrees C) is exactly the 331.3 m/s coefficient.
    #[test]
    fn speed_of_sound_at_zero_c() {
        let c = speed_of_sound(0.0).unwrap();
        assert!((c - C0).abs() < EPS, "got {c} m/s");
    }

    /// Closed-form check against an independent evaluation of the formula.
    #[test]
    fn speed_of_sound_matches_formula() {
        for &t in &[-40.0f64, -10.0, 0.0, 15.0, 25.0, 100.0] {
            let expected = C0 * (1.0 + t / 273.15).sqrt();
            let got = speed_of_sound(t).unwrap();
            assert!((got - expected).abs() < EPS, "T={t}: {got} vs {expected}");
        }
    }

    /// Temperatures at or below absolute zero are rejected.
    #[test]
    fn speed_of_sound_rejects_below_absolute_zero() {
        assert!(speed_of_sound(-273.15).is_err());
        assert!(speed_of_sound(-300.0).is_err());
        assert!(speed_of_sound(f64::NAN).is_err());
    }

    /// A stationary source and observer return the emitted frequency.
    #[test]
    fn stationary_geometry_is_unshifted() {
        let f = doppler_shift(440.0, 0.0, 0.0, 343.0).unwrap();
        assert!((f - 440.0).abs() < EPS, "got {f}");
    }

    /// Approaching > stationary-source > receding (monotone ordering).
    #[test]
    fn approaching_exceeds_source_exceeds_receding() {
        let c = 343.0;
        let f0 = 1000.0;
        let approaching = doppler_shift(f0, 0.0, 30.0, c).unwrap();
        let receding = doppler_shift(f0, 0.0, -30.0, c).unwrap();
        assert!(
            approaching > f0,
            "approaching {approaching} should exceed source {f0}"
        );
        assert!(
            f0 > receding,
            "source {f0} should exceed receding {receding}"
        );
        assert!(approaching > receding, "{approaching} !> {receding}");
    }

    /// Closed-form value: source approaching at 34.3 m/s (= c/10) gives
    /// f_obs = f * c/(c - c/10) = f * 10/9.
    #[test]
    fn approaching_source_exact_value() {
        let c = 343.0;
        let f = doppler_shift(900.0, 0.0, c / 10.0, c).unwrap();
        let expected = 900.0 * 10.0 / 9.0;
        assert!((f - expected).abs() < 1e-6, "got {f}, expected {expected}");
    }

    /// Observer moving toward the source also raises the pitch, by the
    /// numerator factor (c + v_o)/c.
    #[test]
    fn observer_approaching_raises_pitch_by_known_factor() {
        let c = 343.0;
        let f = doppler_shift(500.0, c / 100.0, 0.0, c).unwrap();
        let expected = 500.0 * (c + c / 100.0) / c;
        assert!((f - expected).abs() < 1e-9, "got {f}, expected {expected}");
        assert!(f > 500.0);
    }

    /// A sonic / supersonic source (v_s >= c) is a reported singularity.
    #[test]
    fn sonic_source_is_singular() {
        let err = doppler_shift(1000.0, 0.0, 343.0, 343.0).unwrap_err();
        assert_eq!(err.code(), "acoustics.doppler_singularity");
        assert!(doppler_shift(1000.0, 0.0, 400.0, 343.0).is_err());
    }

    /// Invalid scalar arguments are rejected.
    #[test]
    fn invalid_arguments_rejected() {
        assert!(doppler_shift(-1.0, 0.0, 0.0, 343.0).is_err());
        assert!(doppler_shift(440.0, 0.0, 0.0, 0.0).is_err());
        assert!(doppler_shift(440.0, f64::NAN, 0.0, 343.0).is_err());
        assert!(doppler_shift(440.0, 0.0, f64::INFINITY, 343.0).is_err());
    }
}
