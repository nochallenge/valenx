//! Launch-site geometry: the relationship between launch latitude,
//! launch azimuth, and the orbital inclination that results.
//!
//! From a launch site at geodetic latitude `φ`, an instantaneous launch
//! on inertial azimuth `β` (measured clockwise from north) injects into
//! an orbit of inclination `i` given by the spherical-triangle relation
//!
//! ```text
//!   cos i = cos φ · sin β
//! ```
//!
//! Two consequences fall out and are provided here: the **minimum
//! inclination** reachable directly from a site is `|φ|` (you cannot
//! launch into an orbit less inclined than your latitude without a
//! plane change), and the eastward **Earth-rotation velocity** the site
//! carries (`ω·R·cos φ`) is a free contribution to a prograde launch.

use crate::constants::{OMEGA_EARTH, R_EARTH};

/// Inertial launch azimuth (rad, clockwise from north) that injects into
/// an orbit of the given `inclination` from the given `latitude`
/// (both rad).
///
/// Returns `None` when the inclination is unreachable directly
/// (`inclination < |latitude|`). For a reachable target there are two
/// azimuths (ascending / descending); this returns the prograde
/// ascending solution in `[0, π/2]` for `i ≤ π/2`.
pub fn azimuth_for_inclination(latitude: f64, inclination: f64) -> Option<f64> {
    let ratio = inclination.cos() / latitude.cos();
    if !ratio.is_finite() || ratio.abs() > 1.0 {
        return None;
    }
    Some(ratio.asin())
}

/// Minimum orbital inclination (rad) reachable from a launch site at the
/// given `latitude` (rad) without an out-of-plane maneuver: `|φ|`.
pub fn min_inclination(latitude: f64) -> f64 {
    latitude.abs()
}

/// Eastward surface speed (m/s) due to Earth's rotation at the given
/// `latitude` (rad): `ω·R·cos φ`. This is the free prograde velocity a
/// due-east launch starts with — maximal at the equator, zero at the
/// poles.
pub fn earth_rotation_speed(latitude: f64) -> f64 {
    OMEGA_EARTH * R_EARTH * latitude.cos()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equator_due_east_is_zero_inclination() {
        // From the equator, a due-east (β = 90°) launch -> equatorial
        // (i = 0) orbit.
        let beta = azimuth_for_inclination(0.0, 0.0).unwrap();
        assert!(
            (beta.to_degrees() - 90.0).abs() < 1e-6,
            "β {}",
            beta.to_degrees()
        );
    }

    #[test]
    fn ksc_to_iss_inclination() {
        // KSC latitude 28.5°, ISS inclination 51.6° -> azimuth ~45°.
        let beta = azimuth_for_inclination(28.5_f64.to_radians(), 51.6_f64.to_radians()).unwrap();
        assert!(
            (beta.to_degrees() - 44.96).abs() < 0.2,
            "β {}",
            beta.to_degrees()
        );
    }

    #[test]
    fn inclination_below_latitude_is_unreachable() {
        // Can't reach a 20° orbit directly from a 28.5° site.
        assert!(azimuth_for_inclination(28.5_f64.to_radians(), 20.0_f64.to_radians()).is_none());
        assert!((min_inclination(28.5_f64.to_radians()).to_degrees() - 28.5).abs() < 1e-9);
    }

    #[test]
    fn polar_launch_is_due_north() {
        // A 90° (polar) orbit needs a due-north/south launch from the
        // equator: cos i = 0 -> sin β = 0 -> β = 0.
        let beta = azimuth_for_inclination(0.0, 90.0_f64.to_radians()).unwrap();
        assert!(beta.abs() < 1e-9);
    }

    #[test]
    fn rotation_speed_peaks_at_equator() {
        let eq = earth_rotation_speed(0.0);
        assert!((eq - 465.1).abs() < 1.0, "{eq}");
        // At 28.5° it is reduced by cos(28.5°).
        let ksc = earth_rotation_speed(28.5_f64.to_radians());
        assert!((ksc - eq * 28.5_f64.to_radians().cos()).abs() < 1e-6);
        // Zero at the pole.
        assert!(earth_rotation_speed(90.0_f64.to_radians()).abs() < 1e-9);
    }
}
