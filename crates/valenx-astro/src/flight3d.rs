//! 3-D powered ascent: orient a planar ascent into a real inclined
//! orbit from a launch site at a given latitude and azimuth.
//!
//! The in-plane ascent dynamics (gravity turn, staging, drag, thrust)
//! are rotationally symmetric, so the planar [`crate::sim`] result *is*
//! the in-plane trajectory. This module embeds that 2-D launch-plane
//! state into the 3-D Earth-centred inertial frame using the launch
//! geometry, yielding the full [`crate::orbit3d::ClassicalElements`] —
//! in particular the **inclination**, which comes out exactly as
//! `cos i = cos φ · sin β` (latitude φ, azimuth β).
//!
//! Scope: the orbital plane orientation and the in-plane orbit shape are
//! exact; the Earth-rotation launch bonus is inherited from the planar
//! sim (its equatorial value), a small approximation at high latitude.
//! A fully native 3-D integrator (with latitude-correct rotation and
//! out-of-plane steering) is the next step beyond this v1.

use nalgebra::Vector3;

use crate::config::AscentConfig;
use crate::error::AstroError;
use crate::orbit3d::{self, ClassicalElements, StateVector};
use crate::result::AscentResult;
use crate::sim::simulate_ascent;
use crate::vehicle::Vehicle;

/// A 3-D ascent result: the planar flight record plus the embedded 3-D
/// insertion state and its orbital elements.
#[derive(Debug, Clone)]
pub struct Ascent3d {
    /// The underlying planar ascent record (Δv, max-Q, events, …).
    pub planar: AscentResult,
    /// Insertion state in the 3-D ECI frame.
    pub insertion: StateVector,
    /// Classical orbital elements of the inserted orbit.
    pub elements: ClassicalElements,
}

/// Orthonormal basis of the orbital plane in ECI: the launch radial
/// `P̂` (planar +radial axis) and the downrange direction `D̂` (planar
/// +downrange axis) for a site at `latitude` launching on `azimuth`
/// (both rad, azimuth clockwise from north).
fn plane_basis(latitude: f64, azimuth: f64) -> (Vector3<f64>, Vector3<f64>) {
    let (sphi, cphi) = latitude.sin_cos();
    let (sb, cb) = azimuth.sin_cos();
    // Launch radial (local up) and local ENU horizontal directions at
    // longitude 0.
    let up = Vector3::new(cphi, 0.0, sphi);
    let east = Vector3::new(0.0, 1.0, 0.0);
    let north = Vector3::new(-sphi, 0.0, cphi);
    // Downrange direction on the given azimuth (from north toward east).
    let downrange = sb * east + cb * north;
    (up, downrange)
}

/// Fly a planar ascent and embed it into a 3-D orbit launched from
/// `latitude` on `azimuth` (both rad).
pub fn ascent_to_orbit(
    vehicle: &Vehicle,
    config: &AscentConfig,
    latitude: f64,
    azimuth: f64,
) -> Result<Ascent3d, AstroError> {
    let planar = simulate_ascent(vehicle, config)?;
    let (e_radial, e_downrange) = plane_basis(latitude, azimuth);

    let p = planar.final_position_m;
    let v = planar.final_velocity_ms;
    let position = p[0] * e_radial + p[1] * e_downrange;
    let velocity = v[0] * e_radial + v[1] * e_downrange;
    let insertion = StateVector { position, velocity };
    let elements = orbit3d::rv_to_coe(&insertion)?;

    Ok(Ascent3d {
        planar,
        insertion,
        elements,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::launch::azimuth_for_inclination;
    use crate::presets;

    #[test]
    fn equatorial_due_east_launch_is_low_inclination() {
        let v = presets::two_stage_medium_lift();
        let c = presets::leo_insertion_config();
        // Equator, due east (azimuth 90°).
        let a = ascent_to_orbit(&v, &c, 0.0, 90.0_f64.to_radians()).unwrap();
        assert!(
            a.elements.inclination.to_degrees() < 0.5,
            "i {}",
            a.elements.inclination.to_degrees()
        );
    }

    #[test]
    fn launch_hits_target_inclination() {
        let v = presets::two_stage_medium_lift();
        let c = presets::leo_insertion_config();
        let lat = 28.5_f64.to_radians();
        let target = 51.6_f64.to_radians();
        let beta = azimuth_for_inclination(lat, target).unwrap();
        let a = ascent_to_orbit(&v, &c, lat, beta).unwrap();
        // Inclination comes out exactly as cos i = cos φ sin β.
        assert!(
            (a.elements.inclination.to_degrees() - 51.6).abs() < 0.2,
            "i {}",
            a.elements.inclination.to_degrees()
        );
    }

    #[test]
    fn three_d_orbit_shape_matches_planar() {
        // Embedding preserves the in-plane orbit: 3-D apoapsis/periapsis
        // radii match the planar result.
        let v = presets::two_stage_medium_lift();
        let c = presets::leo_insertion_config();
        let lat = 28.5_f64.to_radians();
        let beta = azimuth_for_inclination(lat, 51.6_f64.to_radians()).unwrap();
        let a = ascent_to_orbit(&v, &c, lat, beta).unwrap();
        assert!(
            (a.elements.apoapsis_radius() - a.planar.orbit.apoapsis_radius).abs() < 5_000.0,
            "apo 3d {} vs planar {}",
            a.elements.apoapsis_radius(),
            a.planar.orbit.apoapsis_radius
        );
        // Near-circular insertion stays near-circular in 3-D.
        assert!(
            a.elements.eccentricity < 0.05,
            "ecc {}",
            a.elements.eccentricity
        );
    }
}
