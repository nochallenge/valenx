//! Classic textbook gyroscopic-couple cases: an aeroplane making a turn
//! and a ship steering / pitching. Each helper reduces to the core
//! relation `C = I*w*wp` after computing the appropriate precession rate
//! `wp` from the vehicle motion, and reports the sense of the **reactive**
//! couple the rotor exerts on the frame.

use crate::couple::gyroscopic_couple;
use crate::error::{positive, Result};

/// Sense of a gyroscopic couple relative to the rotor and the imposed
/// precession. The **active** couple is what the structure applies to the
/// rotor to force the precession; the **reactive** couple is the equal
/// and opposite reaction the rotor applies back on its bearings / frame.
///
/// `C_active = +I*w*wp` and `C_reactive = -C_active`, so they have equal
/// magnitude and opposite sign along the couple axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoupleSense {
    /// The couple supplied to the rotor to sustain the precession.
    Active,
    /// The reaction couple the rotor exerts on the frame (opposite to
    /// [`CoupleSense::Active`]).
    Reactive,
}

impl CoupleSense {
    /// Multiplier applied to the active-couple magnitude to get the
    /// signed couple for this sense: `+1` for active, `-1` for reactive.
    #[inline]
    pub fn sign(self) -> f64 {
        match self {
            CoupleSense::Active => 1.0,
            CoupleSense::Reactive => -1.0,
        }
    }

    /// Signed couple for this sense given an active-couple magnitude.
    ///
    /// # Errors
    ///
    /// Returns [`crate::GyroError::Negative`] if `magnitude` is `< 0`, or
    /// [`crate::GyroError::NonFinite`] if it is `NaN`/`Inf`.
    #[inline]
    pub fn signed(self, magnitude: f64) -> Result<f64> {
        let magnitude = crate::error::non_negative("magnitude", magnitude)?;
        Ok(self.sign() * magnitude)
    }
}

/// Reactive gyroscopic couple on the bearings of an aeroplane engine /
/// propeller while the aeroplane makes a level turn.
///
/// The turn imposes a yaw (precession) rate `wp = v / R` on the spin
/// axis, where `v` is the forward speed (m/s) and `R` the turn radius
/// (m). The reactive couple magnitude is then `C = I * w * wp`.
///
/// # Errors
///
/// Returns an error if `inertia`, `spin`, `speed`, or `radius` is not
/// strictly positive, or if any input is non-finite.
pub fn aeroplane_turn_couple(inertia: f64, spin: f64, speed: f64, radius: f64) -> Result<f64> {
    let speed = positive("speed", speed)?;
    let radius = positive("radius", radius)?;
    let precession = speed / radius;
    gyroscopic_couple(inertia, spin, precession)
}

/// Reactive gyroscopic couple on a ship's turbine rotor while the ship
/// steers (changes heading) at a steady rate `wp` (rad/s).
///
/// This is just the core relation with the steering rate supplied
/// directly as the precession rate; provided as a named entry point so
/// the steering case reads clearly at the call site.
///
/// # Errors
///
/// Returns an error if `inertia` or `spin` is not strictly positive,
/// if `steer_rate` is negative, or if any input is non-finite.
#[inline]
pub fn ship_steering_couple(inertia: f64, spin: f64, steer_rate: f64) -> Result<f64> {
    gyroscopic_couple(inertia, spin, steer_rate)
}

/// Peak reactive gyroscopic couple on a ship's turbine rotor due to
/// simple-harmonic pitching.
///
/// The bow pitches as `theta = phi0 * sin(2*pi*t / period)`, so the
/// angular velocity of the pitch is `dtheta/dt`, whose maximum magnitude
/// is the angular amplitude times the circular frequency:
///
/// `wp_max = phi0 * (2*pi / period)`
///
/// The peak couple is `C_max = I * w * wp_max`.
///
/// # Arguments
///
/// `inertia` rotor polar moment (kg*m^2), `spin` rotor spin rate (rad/s),
/// `amplitude` pitch angular amplitude `phi0` (rad), `period` pitch period
/// (s).
///
/// # Errors
///
/// Returns an error if `inertia`, `spin`, `amplitude`, or `period` is not
/// strictly positive, or if any input is non-finite.
pub fn ship_pitching_peak_couple(
    inertia: f64,
    spin: f64,
    amplitude: f64,
    period: f64,
) -> Result<f64> {
    let amplitude = positive("amplitude", amplitude)?;
    let period = positive("period", period)?;
    let precession_max = amplitude * (2.0 * core::f64::consts::PI / period);
    gyroscopic_couple(inertia, spin, precession_max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::couple::rpm_to_rad_s;
    use crate::error::GyroError;

    const EPS: f64 = 1e-9;

    #[test]
    fn couple_sense_signs() {
        assert!((CoupleSense::Active.sign() - 1.0).abs() < EPS);
        assert!((CoupleSense::Reactive.sign() + 1.0).abs() < EPS);
        // Active and reactive are equal in magnitude, opposite in sign.
        let a = CoupleSense::Active.signed(1047.197755).unwrap();
        let r = CoupleSense::Reactive.signed(1047.197755).unwrap();
        assert!((a + r).abs() < EPS);
        assert!((a - 1047.197755).abs() < EPS);
    }

    #[test]
    fn couple_sense_rejects_negative_magnitude() {
        assert!(matches!(
            CoupleSense::Active.signed(-1.0),
            Err(GyroError::Negative { .. })
        ));
    }

    #[test]
    fn aeroplane_turn_uses_v_over_r() {
        // I = 20, N = 1000 rpm, v = 50 m/s, R = 100 m => wp = 0.5 rad/s.
        // C = 20 * 104.71975511965977 * 0.5 = 1047.1975512 N*m (hand-worked).
        let w = rpm_to_rad_s(1000.0).unwrap();
        let c = aeroplane_turn_couple(20.0, w, 50.0, 100.0).unwrap();
        assert!((c - 1047.1975511965977).abs() < 1e-4);

        // It must equal the bare relation evaluated at wp = v/R.
        let direct = crate::gyroscopic_couple(20.0, w, 50.0 / 100.0).unwrap();
        assert!((c - direct).abs() < EPS);
    }

    #[test]
    fn ship_steering_matches_hand_worked_number() {
        // Rotor: m = 750 kg, k = 0.3 m => I = 750 * 0.3^2 = 67.5 kg*m^2.
        // N = 1800 rpm => w = 188.49555921538757 rad/s.
        // Steering rate wp = 0.1 rad/s.
        // C = 67.5 * 188.49555921538757 * 0.1 = 1272.3450247 N*m.
        let inertia = 750.0_f64 * 0.3 * 0.3;
        assert!((inertia - 67.5).abs() < EPS);
        let w = rpm_to_rad_s(1800.0).unwrap();
        let c = ship_steering_couple(inertia, w, 0.1).unwrap();
        // Independent recomputation of the same product.
        let reference = inertia * w * 0.1;
        assert!((c - reference).abs() < EPS);
        assert!((c - 1272.345024703807).abs() < 1e-3);
    }

    #[test]
    fn ship_pitching_peak_uses_amplitude_times_omega() {
        // phi0 = 0.1 rad, T = 20 s => wp_max = 0.1 * 2*pi/20 = 0.0314159265 rad/s.
        let inertia = 67.5;
        let w = rpm_to_rad_s(1800.0).unwrap();
        let c = ship_pitching_peak_couple(inertia, w, 0.1, 20.0).unwrap();

        // Independent recomputation of wp_max and the couple.
        let wp_max = 0.1 * (2.0 * core::f64::consts::PI / 20.0);
        assert!((wp_max - 0.031_415_926_535_897_93).abs() < EPS);
        let reference = inertia * w * wp_max;
        assert!((c - reference).abs() < EPS);
        assert!((c - 399.718_972_675_984).abs() < 1e-3);
    }

    #[test]
    fn ship_pitching_zero_period_is_rejected() {
        let w = rpm_to_rad_s(1800.0).unwrap();
        assert_eq!(
            ship_pitching_peak_couple(67.5, w, 0.1, 0.0)
                .unwrap_err()
                .code(),
            "gyroscope.non_positive"
        );
    }

    #[test]
    fn aeroplane_turn_rejects_non_positive_radius() {
        let w = rpm_to_rad_s(1000.0).unwrap();
        assert_eq!(
            aeroplane_turn_couple(20.0, w, 50.0, 0.0)
                .unwrap_err()
                .code(),
            "gyroscope.non_positive"
        );
        assert!(matches!(
            aeroplane_turn_couple(20.0, w, f64::NAN, 100.0),
            Err(GyroError::NonFinite { .. })
        ));
    }

    #[test]
    fn ship_steering_couple_vanishes_at_zero_steer_rate() {
        // Holding a straight course (wp = 0) produces no gyroscopic couple.
        let w = rpm_to_rad_s(1800.0).unwrap();
        let c = ship_steering_couple(67.5, w, 0.0).unwrap();
        assert!(c.abs() < EPS);
    }
}
