//! Launch-window geometry and orbital-plane targeting.
//!
//! Extends [`crate::launch`] from "what inclination does an azimuth give"
//! to "how do I hit a *target orbital plane*, and what does it cost if I
//! cannot reach it directly". From a launch site at geodetic latitude
//! `φ`, the reachable inclinations are bounded below by `|φ|` (the
//! spherical-triangle relation `cos i = cos φ · sin β` has no solution for
//! `i < |φ|`). For a reachable target this returns the launch azimuth;
//! for an unreachable one it gives the **residual plane-change `Δv`** an
//! on-orbit maneuver must supply, from the exact impulsive relation
//!
//! ```text
//!   Δv = 2·v·sin(Δi/2)
//! ```
//!
//! at orbital speed `v`, where `Δi` is the plane-change angle. These are
//! the numbers that decide whether a payload's target plane is a
//! free launch-azimuth choice or an expensive dogleg / on-orbit burn.
//!
//! All relations are exact closed forms, pinned directly by the tests
//! (and cross-checked against [`crate::launch::azimuth_for_inclination`]
//! for self-consistency of the `cos i = cos φ · sin β` round-trip).

use serde::{Deserialize, Serialize};

use crate::constants::MU_EARTH;
use crate::error::{AstroError, Result};
use crate::launch;

/// The launch-azimuth solution for a target orbital plane from a site.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LaunchAzimuth {
    /// Inertial launch azimuth (rad, clockwise from north) for the
    /// **ascending** (prograde) pass.
    pub ascending: f64,
    /// Inertial launch azimuth (rad) for the **descending** pass — the
    /// other azimuth that reaches the same inclination (`π − β`).
    pub descending: f64,
}

/// Inertial launch azimuth(s) (rad, clockwise from north) that inject into
/// an orbit of the given `inclination` (rad) from a site at `latitude`
/// (rad), via `cos i = cos φ · sin β`.
///
/// Returns both the ascending solution (`β ∈ [0, π/2]` for a prograde
/// target) and its descending counterpart (`π − β`).
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if either angle is non-finite,
/// or [`AstroError::OutOfRange`]-style **`InvalidParameter`** ("inclination
/// unreachable") when `inclination < |latitude|` — the plane cannot be hit
/// directly and needs a dogleg or an on-orbit plane change (see
/// [`plane_change_delta_v`]).
pub fn launch_azimuth(latitude: f64, inclination: f64) -> Result<LaunchAzimuth> {
    if !latitude.is_finite() || !inclination.is_finite() {
        return Err(AstroError::InvalidParameter(
            "latitude and inclination must be finite",
        ));
    }
    match launch::azimuth_for_inclination(latitude, inclination) {
        Some(beta) => Ok(LaunchAzimuth {
            ascending: beta,
            descending: std::f64::consts::PI - beta,
        }),
        None => Err(AstroError::InvalidParameter(
            "inclination unreachable directly from this latitude (i < |phi|)",
        )),
    }
}

/// Impulsive plane-change `Δv` (m/s) for a plane rotation of
/// `delta_inclination` (rad) at orbital `speed` (m/s):
/// `Δv = 2·v·sin(Δi/2)`.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if `speed` is non-finite or
/// negative, or `delta_inclination` is non-finite or negative.
pub fn plane_change_delta_v(speed: f64, delta_inclination: f64) -> Result<f64> {
    if !speed.is_finite() || speed < 0.0 {
        return Err(AstroError::InvalidParameter(
            "speed must be finite and >= 0",
        ));
    }
    if !delta_inclination.is_finite() || delta_inclination < 0.0 {
        return Err(AstroError::InvalidParameter(
            "delta_inclination must be finite and >= 0",
        ));
    }
    Ok(2.0 * speed * (0.5 * delta_inclination).sin())
}

/// The residual plane-change `Δv` (m/s) to reach a `target_inclination`
/// (rad) that lies **below** a site `latitude` (rad), at circular orbital
/// `speed` (m/s). The residual angle is `|φ| − i`; a reachable target
/// (`i ≥ |φ|`) needs no plane change and returns `0`.
///
/// # Errors
///
/// As [`plane_change_delta_v`], plus a non-finite `latitude` /
/// `target_inclination`.
pub fn residual_plane_change_delta_v(
    latitude: f64,
    target_inclination: f64,
    speed: f64,
) -> Result<f64> {
    if !latitude.is_finite() || !target_inclination.is_finite() {
        return Err(AstroError::InvalidParameter(
            "latitude and target_inclination must be finite",
        ));
    }
    let residual = (latitude.abs() - target_inclination).max(0.0);
    plane_change_delta_v(speed, residual)
}

/// Circular orbital speed (m/s) at altitude `altitude_m` (m) above the
/// equatorial radius — a convenience for sizing the plane-change `Δv`.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if `altitude_m` is non-finite
/// or `<= −R⊕` (a radius at or below the centre).
pub fn circular_speed_at_altitude(altitude_m: f64) -> Result<f64> {
    if !altitude_m.is_finite() {
        return Err(AstroError::InvalidParameter("altitude must be finite"));
    }
    let r = crate::constants::R_EARTH + altitude_m;
    if r <= 0.0 {
        return Err(AstroError::InvalidParameter(
            "radius must be > 0 (altitude above the centre)",
        ));
    }
    Ok((MU_EARTH / r).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::R_EARTH;

    #[test]
    fn azimuth_round_trips_with_launch_relation() {
        // ORACLE: cos i = cos φ · sin β. Solve β, then recover i.
        let phi = 28.5_f64.to_radians();
        let i = 51.6_f64.to_radians();
        let az = launch_azimuth(phi, i).expect("reachable");
        // KSC->ISS azimuth ≈ 45°.
        assert!(
            (az.ascending.to_degrees() - 44.975).abs() < 0.1,
            "β {}",
            az.ascending.to_degrees()
        );
        // Recover the inclination from the azimuth: must round-trip.
        let i_rec = (phi.cos() * az.ascending.sin()).acos();
        assert!(
            (i_rec - i).abs() < 1e-12,
            "recovered i {}",
            i_rec.to_degrees()
        );
        // The descending azimuth is π − β, reaching the same inclination.
        let i_rec2 = (phi.cos() * az.descending.sin()).acos();
        assert!((i_rec2 - i).abs() < 1e-12);
    }

    #[test]
    fn equatorial_due_east_from_equator() {
        // φ=0, i=0 -> due east (β=90°); descending is 90° as well (π−π/2).
        let az = launch_azimuth(0.0, 0.0).expect("reachable");
        assert!((az.ascending.to_degrees() - 90.0).abs() < 1e-9);
    }

    #[test]
    fn unreachable_inclination_below_latitude_errors() {
        // Can't reach a 20° plane directly from a 28.5° site.
        assert!(launch_azimuth(28.5_f64.to_radians(), 20.0_f64.to_radians()).is_err());
    }

    #[test]
    fn plane_change_delta_v_is_exact() {
        // Δv = 2·v·sin(Δi/2). v=7668.56 (≈400 km LEO), Δi=1°:
        // 2·7668.56·sin(0.5°) = 133.84 m/s.
        let v = circular_speed_at_altitude(400_000.0).expect("ok");
        assert!((v - 7_668.56).abs() < 0.1, "v = {v}");
        let dv = plane_change_delta_v(v, 1.0_f64.to_radians()).expect("ok");
        let expected = 2.0 * v * (0.5_f64.to_radians()).sin();
        assert!((dv - expected).abs() < 1e-9);
        assert!((dv - 133.84).abs() < 0.5, "dv = {dv}");
        // A zero plane change costs nothing.
        assert!(plane_change_delta_v(v, 0.0).expect("ok").abs() < 1e-12);
    }

    #[test]
    fn residual_plane_change_only_when_below_latitude() {
        let v = circular_speed_at_altitude(400_000.0).expect("ok");
        let phi = 28.5_f64.to_radians();
        // Reachable (i >= φ): no residual.
        assert!(
            residual_plane_change_delta_v(phi, 51.6_f64.to_radians(), v)
                .expect("ok")
                .abs()
                < 1e-12
        );
        // Below latitude (i = 0, equatorial from KSC): residual = φ = 28.5°.
        let dv = residual_plane_change_delta_v(phi, 0.0, v).expect("ok");
        let expected = plane_change_delta_v(v, phi).expect("ok");
        assert!((dv - expected).abs() < 1e-12);
        // ~3775 m/s — the well-known cost of reaching GEO's plane from KSC.
        assert!((dv - 3_775.0).abs() < 5.0, "residual dv = {dv}");
    }

    #[test]
    fn circular_speed_matches_vis_viva() {
        let v = circular_speed_at_altitude(0.0).expect("ok");
        let expected = (MU_EARTH / R_EARTH).sqrt();
        assert!((v - expected).abs() < 1e-9);
        // Surface circular speed ≈ 7.9 km/s.
        assert!((v - 7_905.0).abs() < 5.0, "v_surface = {v}");
    }

    #[test]
    fn rejects_non_physical_inputs() {
        assert!(launch_azimuth(f64::NAN, 1.0).is_err());
        assert!(plane_change_delta_v(-1.0, 0.1).is_err());
        assert!(plane_change_delta_v(7_000.0, -0.1).is_err());
        assert!(circular_speed_at_altitude(f64::NAN).is_err());
        assert!(circular_speed_at_altitude(-R_EARTH - 1.0).is_err());
    }
}
