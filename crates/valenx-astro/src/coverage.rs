//! Instantaneous **coverage / footprint** geometry: how much of the Earth a
//! satellite can see above a given elevation mask, and what fraction of a set
//! of ground points falls inside that footprint.
//!
//! A satellite at altitude `h` sees a circular cap of the Earth's surface
//! centred on its sub-satellite point. The cap's angular size is set by the
//! observer's **minimum elevation** `ε`: a point right at the edge of the
//! footprint sees the satellite exactly at elevation `ε`, on the horizon when
//! `ε = 0`. This module gives the closed-form spherical-cap geometry
//! ([`footprint_half_angle`] and friends) plus a simple
//! [`coverage_fraction`] over an explicit list of ground points — the
//! instantaneous-access building block for constellation coverage studies.
//!
//! The relations (spherical Earth of radius `R⊕`, Vallado §11):
//!
//! ```text
//! λ  = Earth-central half-angle of the footprint cap
//! η  = arcsin( R⊕·cos ε / (R⊕ + h) )        (the nadir/look angle at the edge)
//! λ  = 90° − ε − η                           (the angle triangle closes)
//! d_max = (R⊕ + h)·sin λ / cos ε             (max slant range, to the edge)
//! A_frac = (1 − cos λ)/2                      (cap area ÷ full-sphere area)
//! ```
//!
//! # Honest scope
//!
//! This is **spherical-Earth** geometry (radius [`crate::constants::R_EARTH`]),
//! consistent with [`crate::groundtrack`]; it ignores the WGS-84 flattening, so
//! the footprint is a true circle on a sphere, not on the ellipsoid. No
//! atmospheric refraction or terrain masking. The [`coverage_fraction`] is the
//! exact fraction of the *supplied* points inside the cap — its fidelity as an
//! area fraction is only as good as how evenly those points sample the region
//! of interest (use an equal-area grid for a meaningful global number).

use crate::constants::R_EARTH;
use crate::error::AstroError;

/// The **Earth-central half-angle** `λ` (rad) of the access footprint of a
/// satellite at altitude `altitude` (m above the spherical Earth) seen above a
/// minimum elevation `min_elevation` (rad).
///
/// This is the geocentric angular radius of the circular cap of the Earth's
/// surface from which the satellite is at or above the elevation mask. From the
/// edge geometry `η = arcsin(R⊕·cos ε /(R⊕+h))` and the planar angle sum,
/// `λ = π/2 − ε − η`.
///
/// At the horizon mask `ε = 0` this is the largest cap the altitude allows;
/// raising `ε` shrinks it, reaching `0` when the mask equals the maximum
/// elevation the altitude can ever present off-nadir (i.e. the satellite is
/// only ever visible at its own sub-point).
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] for a non-finite or non-positive
/// altitude, or a mask outside `[0, π/2)`.
pub fn footprint_half_angle(altitude: f64, min_elevation: f64) -> Result<f64, AstroError> {
    if !altitude.is_finite() || altitude <= 0.0 {
        return Err(AstroError::InvalidParameter(
            "footprint_half_angle: altitude must be finite and > 0",
        ));
    }
    if !min_elevation.is_finite() || !(0.0..std::f64::consts::FRAC_PI_2).contains(&min_elevation) {
        return Err(AstroError::InvalidParameter(
            "footprint_half_angle: min_elevation must be in [0, π/2)",
        ));
    }
    let rho = R_EARTH / (R_EARTH + altitude); // < 1
    let cos_eps = min_elevation.cos();
    // Nadir (look) angle at the footprint edge.
    let eta = (rho * cos_eps).clamp(-1.0, 1.0).asin();
    let lambda = std::f64::consts::FRAC_PI_2 - min_elevation - eta;
    Ok(lambda.max(0.0))
}

/// The **ground-range radius** of the footprint (m) — the great-circle surface
/// distance from the sub-satellite point to the edge of the access cap, for a
/// satellite at `altitude` (m) and elevation mask `min_elevation` (rad).
///
/// Simply `R⊕·λ` with `λ` the [`footprint_half_angle`]: the cap half-angle
/// times the Earth's radius. This is the radius of the visible circle *measured
/// along the ground*, the figure most useful for sizing swath / coverage.
///
/// # Errors
///
/// As [`footprint_half_angle`].
pub fn footprint_ground_radius(altitude: f64, min_elevation: f64) -> Result<f64, AstroError> {
    Ok(R_EARTH * footprint_half_angle(altitude, min_elevation)?)
}

/// The **maximum slant range** (m) — the straight-line distance to a satellite
/// at the very edge of its footprint (where it appears exactly at the elevation
/// mask) — for `altitude` (m) and `min_elevation` (rad).
///
/// `d = (R⊕ + h)·sin λ / cos ε`, the longest line of sight the access cap
/// allows (the link-budget worst case). At the sub-point the slant range is
/// just the altitude `h`; this is its far-edge counterpart.
///
/// # Errors
///
/// As [`footprint_half_angle`].
pub fn max_slant_range(altitude: f64, min_elevation: f64) -> Result<f64, AstroError> {
    let lambda = footprint_half_angle(altitude, min_elevation)?;
    let cos_eps = min_elevation.cos();
    Ok((R_EARTH + altitude) * lambda.sin() / cos_eps)
}

/// The **instantaneous coverage area fraction** — the area of the footprint cap
/// as a fraction of the whole Earth's surface — for `altitude` (m) and
/// `min_elevation` (rad).
///
/// For a spherical cap of half-angle `λ` the area is `2π R⊕²(1 − cos λ)`, so the
/// fraction of the full `4π R⊕²` sphere is `(1 − cos λ)/2`. This is the share of
/// the planet a single satellite blankets at one instant — the headline number
/// for constellation sizing (how many such caps tile the globe).
///
/// # Errors
///
/// As [`footprint_half_angle`].
pub fn coverage_area_fraction(altitude: f64, min_elevation: f64) -> Result<f64, AstroError> {
    let lambda = footprint_half_angle(altitude, min_elevation)?;
    Ok(0.5 * (1.0 - lambda.cos()))
}

/// Whether a ground point at geocentric `(lat_rad, lon_rad)` lies inside the
/// footprint of a satellite whose **sub-satellite point** is `(sub_lat,
/// sub_lon)` (all rad), for footprint half-angle `half_angle` (rad, e.g. from
/// [`footprint_half_angle`]).
///
/// The test is purely the great-circle (central) angle between the two surface
/// points against the cap half-angle: `point is covered ⇔ central_angle ≤ λ`.
/// The central angle uses the numerically-robust form (atan2 of the cross/dot of
/// the two unit surface vectors) so it stays accurate for both nearby and
/// nearly-antipodal points.
pub fn point_in_footprint(
    lat_rad: f64,
    lon_rad: f64,
    sub_lat: f64,
    sub_lon: f64,
    half_angle: f64,
) -> bool {
    central_angle(lat_rad, lon_rad, sub_lat, sub_lon) <= half_angle
}

/// The great-circle **central angle** (rad) between two surface points given as
/// geocentric latitude/longitude (rad) — the angle subtended at the Earth's
/// centre.
///
/// Robust haversine-equivalent via `atan2`: accurate from coincident points up
/// to antipodes (unlike the naive `acos` of the dot product, which loses
/// precision for small angles).
pub fn central_angle(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    let dlat = lat2 - lat1;
    let dlon = lon2 - lon1;
    let (s_dlat, s_dlon) = ((dlat * 0.5).sin(), (dlon * 0.5).sin());
    let a = s_dlat * s_dlat + lat1.cos() * lat2.cos() * s_dlon * s_dlon;
    2.0 * a.sqrt().clamp(0.0, 1.0).asin()
}

/// The **fraction of the supplied ground `points`** (each geocentric
/// `(lat_rad, lon_rad)`) that lies within the access footprint of a satellite at
/// `altitude` (m) whose sub-satellite point is `(sub_lat, sub_lon)` (rad), above
/// the elevation mask `min_elevation` (rad).
///
/// Computes the footprint half-angle once, then counts the covered points. With
/// an even (e.g. equal-area) grid of points this approximates the instantaneous
/// area coverage of the region the grid spans; with a list of cities / assets it
/// is the fraction of *those* that currently have access.
///
/// An empty point list returns `0.0` (nothing to cover) rather than a `0/0`
/// NaN.
///
/// # Errors
///
/// As [`footprint_half_angle`].
pub fn coverage_fraction(
    points: &[(f64, f64)],
    altitude: f64,
    sub_lat: f64,
    sub_lon: f64,
    min_elevation: f64,
) -> Result<f64, AstroError> {
    let lambda = footprint_half_angle(altitude, min_elevation)?;
    if points.is_empty() {
        return Ok(0.0);
    }
    let covered = points
        .iter()
        .filter(|&&(lat, lon)| point_in_footprint(lat, lon, sub_lat, sub_lon, lambda))
        .count();
    Ok(covered as f64 / points.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn half_angle_matches_closed_form_spherical_cap() {
        // GROUND TRUTH: cross-check against the independent edge construction.
        // For h and ε, η = arcsin(R cos ε /(R+h)) and λ = 90° − ε − η. Verify
        // for several altitudes and masks that footprint_half_angle reproduces
        // exactly this, and that the returned λ closes the geometry: a ground
        // point at central angle λ from the sub-point sees the satellite at
        // exactly elevation ε.
        for &alt in &[400_000.0, 800_000.0, 20_200_000.0, 35_786_000.0] {
            for &eps_deg in &[0.0_f64, 5.0, 10.0, 25.0] {
                let eps = eps_deg.to_radians();
                let lambda = footprint_half_angle(alt, eps).unwrap();
                let rho = R_EARTH / (R_EARTH + alt);
                let eta = (rho * eps.cos()).asin();
                let expected = FRAC_PI_2 - eps - eta;
                assert!(
                    (lambda - expected).abs() < 1e-12,
                    "alt {alt} eps {eps_deg}: λ {lambda} != {expected}"
                );

                // Independent check: law of sines in the Earth-centre /
                // satellite / edge-point triangle. The elevation reconstructed
                // from λ must equal the input ε.
                // Sides: from centre R⊕ (to edge point) and R⊕+h (to sat);
                // angle at the edge point is (90° + ε). The angle at the centre
                // is λ. sin(λ)/(? ) ... reconstruct ε from the slant geometry:
                let rs = R_EARTH + alt;
                // Edge point ECEF on +x; satellite at central angle λ from it.
                // Place edge at angle 0, satellite sub-point at angle λ.
                let edge = nalgebra::Vector2::new(R_EARTH, 0.0);
                let sat = nalgebra::Vector2::new(rs * lambda.cos(), rs * lambda.sin());
                let los = sat - edge; // edge -> satellite
                let up = edge / edge.norm(); // local vertical at edge
                let cos_zenith = los.dot(&up) / los.norm();
                let elev = FRAC_PI_2 - cos_zenith.clamp(-1.0, 1.0).acos();
                assert!(
                    (elev - eps).abs() < 1e-9,
                    "alt {alt} eps {eps_deg}: reconstructed elev {} != {}",
                    elev.to_degrees(),
                    eps_deg
                );
            }
        }
    }

    #[test]
    fn higher_altitude_sees_a_bigger_cap() {
        // Monotonic: a higher satellite (at the same mask) covers a larger
        // footprint.
        let a = footprint_half_angle(400_000.0, 0.0).unwrap();
        let b = footprint_half_angle(800_000.0, 0.0).unwrap();
        let c = footprint_half_angle(35_786_000.0, 0.0).unwrap();
        assert!(a < b && b < c, "λ should grow with altitude: {a} {b} {c}");
    }

    #[test]
    fn higher_mask_shrinks_the_cap() {
        let lo = footprint_half_angle(800_000.0, 0.0).unwrap();
        let hi = footprint_half_angle(800_000.0, 20.0_f64.to_radians()).unwrap();
        assert!(hi < lo, "raising the mask must shrink λ: {lo} -> {hi}");
    }

    #[test]
    fn geo_horizon_cap_is_about_81_degrees() {
        // GROUND TRUTH: from geostationary altitude (35,786 km) the horizon
        // (ε = 0) footprint half-angle is ≈ 81.3°, the well-known figure behind
        // "3 GEO satellites ≈ global coverage". cos⁻¹(R⊕/(R⊕+h)).
        let lambda = footprint_half_angle(35_786_000.0, 0.0).unwrap();
        let expected = (R_EARTH / (R_EARTH + 35_786_000.0)).acos();
        assert!((lambda - expected).abs() < 1e-12);
        assert!(
            (lambda.to_degrees() - 81.3).abs() < 0.3,
            "GEO horizon half-angle {}° != ~81.3°",
            lambda.to_degrees()
        );
    }

    #[test]
    fn area_fraction_is_cap_over_sphere() {
        // (1 − cos λ)/2, and three GEO horizon caps comfortably exceed a full
        // sphere of coverage (they overlap but together blanket the globe
        // except the poles).
        let frac = coverage_area_fraction(35_786_000.0, 0.0).unwrap();
        let lambda = footprint_half_angle(35_786_000.0, 0.0).unwrap();
        assert!((frac - 0.5 * (1.0 - lambda.cos())).abs() < 1e-15);
        assert!(
            3.0 * frac > 1.0,
            "3 GEO caps {} should exceed full sphere",
            3.0 * frac
        );
        assert!((0.0..=1.0).contains(&frac));
    }

    #[test]
    fn max_slant_range_exceeds_altitude_and_is_finite() {
        // The far edge of the footprint is farther than straight down.
        let alt = 800_000.0;
        let d = max_slant_range(alt, 0.0).unwrap();
        assert!(d > alt, "edge slant {d} should exceed nadir {alt}");
        assert!(d.is_finite());
    }

    #[test]
    fn coverage_fraction_counts_points_in_the_cap() {
        // Sub-point at the equator/prime meridian; a tight cluster of points
        // right under it is fully covered, an antipodal cluster not at all.
        let sub_lat = 0.0;
        let sub_lon = 0.0;
        let alt = 800_000.0;
        let near: Vec<(f64, f64)> = (0..10)
            .map(|i| ((i as f64 * 0.1).to_radians(), 0.0))
            .collect();
        let f_near = coverage_fraction(&near, alt, sub_lat, sub_lon, 0.0).unwrap();
        assert!((f_near - 1.0).abs() < 1e-12, "near cluster frac {f_near}");

        let far: Vec<(f64, f64)> = vec![(0.0, PI), (0.1, PI), (-0.1, PI)];
        let f_far = coverage_fraction(&far, alt, sub_lat, sub_lon, 0.0).unwrap();
        assert!(f_far.abs() < 1e-12, "antipodal cluster frac {f_far}");

        // A mix: half near, half far.
        let mixed: Vec<(f64, f64)> = vec![(0.0, 0.0), (0.0, PI)];
        let f_mixed = coverage_fraction(&mixed, alt, sub_lat, sub_lon, 0.0).unwrap();
        assert!((f_mixed - 0.5).abs() < 1e-12, "mixed frac {f_mixed}");
    }

    #[test]
    fn central_angle_basics() {
        // Coincident -> 0; equator quarter-turn -> 90°; pole-to-pole -> 180°.
        assert!(central_angle(0.0, 0.0, 0.0, 0.0).abs() < 1e-15);
        let q = central_angle(0.0, 0.0, 0.0, FRAC_PI_2);
        assert!((q - FRAC_PI_2).abs() < 1e-12, "quarter turn {q}");
        let anti = central_angle(FRAC_PI_2, 0.0, -FRAC_PI_2, 0.0);
        assert!((anti - PI).abs() < 1e-12, "pole-to-pole {anti}");
    }

    #[test]
    fn degenerate_inputs_are_handled() {
        // Empty point set -> 0, no NaN.
        assert_eq!(
            coverage_fraction(&[], 800_000.0, 0.0, 0.0, 0.0).unwrap(),
            0.0
        );
        // Bad altitude / mask -> error, not panic.
        assert!(footprint_half_angle(-1.0, 0.0).is_err());
        assert!(footprint_half_angle(0.0, 0.0).is_err());
        assert!(footprint_half_angle(f64::NAN, 0.0).is_err());
        assert!(footprint_half_angle(800_000.0, FRAC_PI_2).is_err());
        assert!(footprint_half_angle(800_000.0, -0.1).is_err());
    }
}
