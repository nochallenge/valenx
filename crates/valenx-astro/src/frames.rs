//! Earth reference frames and coordinate transforms for ground-relative
//! geometry: Greenwich sidereal time, geodetic ↔ Earth-fixed (ECEF), and
//! inertial (ECI) ↔ Earth-fixed.
//!
//! The on-orbit machinery in [`crate::orbit3d`] lives entirely in the
//! Earth-centred **inertial** (ECI) frame. To relate a satellite to a point on
//! the ground — a station, a target, a sensor footprint — you must rotate that
//! inertial state into the **Earth-fixed** (ECEF) frame that turns with the
//! planet, and place the ground point on the WGS-84 reference ellipsoid. This
//! module supplies those primitives:
//!
//! 1. [`gmst`] — Greenwich Mean Sidereal Time (rad): the rotation angle of the
//!    Earth-fixed frame about the pole, from the IAU-1982 series in `UT1`
//!    (approximated by `UTC`).
//! 2. [`geodetic_to_ecef`] / [`ecef_to_geodetic`] — a point's geodetic
//!    latitude/longitude/altitude ↔ its ECEF Cartesian position, on the WGS-84
//!    ellipsoid (round-trip stable).
//! 3. [`eci_to_ecef`] / [`ecef_to_eci`] — rotate a position vector between the
//!    inertial and Earth-fixed frames by the Greenwich angle.
//!
//! # Honest scope
//!
//! This is a **mean-rotation** frame model, the natural companion to the rest of
//! the crate's point-mass / spherical-Earth fidelity:
//!
//! - The ECI ↔ ECEF rotation is a *simple spin about the pole* by GMST. It uses
//!   **mean** sidereal time (no equation of the equinoxes / nutation term, so it
//!   is GMST not GAST) and omits **precession, nutation, polar motion and the
//!   length-of-day (`UT1 − UTC`) correction**. The resulting pointing error is
//!   well under a degree — fine for access / coverage / visibility planning, not
//!   a precision-ephemeris / geolocation product.
//! - The geodetic conversion is the full WGS-84 *ellipsoid* (not the spherical
//!   geocentric latitude used by [`crate::groundtrack`]), so station altitudes
//!   and the latitude flattening are handled correctly.

use nalgebra::Vector3;

use crate::constants::{OMEGA_EARTH, R_EARTH, WGS84_ECC_SQ};

/// Julian Date of the J2000.0 epoch (2000-01-01 12:00 TT). Re-exported value
/// shared with [`crate::eclipse::J2000`]; kept here so the frame transforms are
/// self-contained.
pub const J2000: f64 = 2_451_545.0;

/// **Greenwich Mean Sidereal Time** (radians, reduced to `[0, 2π)`) at Julian
/// Date `julian_date` (UT1, approximated by UTC) — the rotation angle of the
/// Earth-fixed frame about the celestial pole, i.e. the right ascension of the
/// Greenwich meridian.
///
/// Uses the IAU-1982 polynomial (Vallado, *Fundamentals of Astrodynamics and
/// Applications*, eq. 3-47): with `T` the number of Julian centuries of UT1
/// since J2000.0,
///
/// ```text
/// GMST = 67310.54841 + (876600ʰ·3600 + 8640184.812866)·T
///                    + 0.093104·T² − 6.2e-6·T³   [seconds of time]
/// ```
///
/// converted to an angle at `360°/86400ˢ` and wrapped. The value is what the
/// ECI → ECEF rotation [`eci_to_ecef`] needs and what [`crate::groundtrack`]
/// takes as its Greenwich angle `θ`.
pub fn gmst(julian_date: f64) -> f64 {
    // Julian centuries of UT1 since J2000.
    let t = (julian_date - J2000) / 36_525.0;
    // GMST in seconds of time (IAU 1982).
    let gmst_sec =
        67_310.548_41 + (876_600.0 * 3_600.0 + 8_640_184.812_866) * t + 0.093_104 * t * t
            - 6.2e-6 * t * t * t;
    // Seconds of time -> radians: 86400 s = 2π. Wrap into [0, 2π).
    let two_pi = std::f64::consts::TAU;
    let mut theta = (gmst_sec * two_pi / 86_400.0) % two_pi;
    if theta < 0.0 {
        theta += two_pi;
    }
    theta
}

/// Convert a geodetic position — `lat_rad`/`lon_rad` (geodetic latitude and
/// east longitude, radians) and `alt_m` (height above the WGS-84 ellipsoid,
/// metres) — to an **ECEF** Cartesian position vector (m).
///
/// Standard closed-form WGS-84 transform (Vallado eq. 3-7) with the
/// prime-vertical radius of curvature `N = a/√(1 − e²·sin²φ)`:
///
/// ```text
/// x = (N + h)·cosφ·cosλ
/// y = (N + h)·cosφ·sinλ
/// z = (N·(1 − e²) + h)·sinφ
/// ```
///
/// The inverse is [`ecef_to_geodetic`].
pub fn geodetic_to_ecef(lat_rad: f64, lon_rad: f64, alt_m: f64) -> Vector3<f64> {
    let (sin_lat, cos_lat) = lat_rad.sin_cos();
    let (sin_lon, cos_lon) = lon_rad.sin_cos();
    let n = R_EARTH / (1.0 - WGS84_ECC_SQ * sin_lat * sin_lat).sqrt();
    let x = (n + alt_m) * cos_lat * cos_lon;
    let y = (n + alt_m) * cos_lat * sin_lon;
    let z = (n * (1.0 - WGS84_ECC_SQ) + alt_m) * sin_lat;
    Vector3::new(x, y, z)
}

/// Convert an **ECEF** Cartesian position (m) to geodetic
/// `(lat_rad, lon_rad, alt_m)` on the WGS-84 ellipsoid — the inverse of
/// [`geodetic_to_ecef`].
///
/// Longitude is exact (`λ = atan2(y, x)`, east-positive, in `(−π, π]`).
/// Latitude/altitude use **Bowring's** iteration on the reduced latitude, which
/// converges to full `f64` precision in a few steps for all near-surface and
/// orbital radii. Returns latitude in `[−π/2, π/2]`.
///
/// Degenerate inputs are handled without a NaN: a point on the polar axis
/// (`x = y = 0`) returns longitude `0` and latitude `±π/2`; the geocentre
/// (the zero vector) returns all zeros.
pub fn ecef_to_geodetic(ecef: Vector3<f64>) -> (f64, f64, f64) {
    let (x, y, z) = (ecef.x, ecef.y, ecef.z);
    let p = (x * x + y * y).sqrt(); // distance from the spin axis

    // On (or extremely near) the polar axis: longitude is undefined, latitude
    // is ±90°, altitude is |z| − b. Guard the 0/0 in atan2 / the iteration.
    if p < 1e-9 {
        let b = R_EARTH * (1.0 - WGS84_ECC_SQ).sqrt();
        if z.abs() < 1e-9 {
            return (0.0, 0.0, 0.0); // the geocentre itself
        }
        let lat = std::f64::consts::FRAC_PI_2 * z.signum();
        return (lat, 0.0, z.abs() - b);
    }

    let lon = y.atan2(x);
    let a = R_EARTH;
    let b = a * (1.0 - WGS84_ECC_SQ).sqrt();
    let ep2 = (a * a - b * b) / (b * b); // second eccentricity squared

    // Bowring's seed: tan(reduced latitude) then geodetic latitude.
    let theta = (z * a).atan2(p * b);
    let (st, ct) = theta.sin_cos();
    let mut lat = (z + ep2 * b * st * st * st).atan2(p - WGS84_ECC_SQ * a * ct * ct * ct);

    // One Newton-style refinement pass tightens it to machine precision.
    for _ in 0..5 {
        let sin_lat = lat.sin();
        let n = a / (1.0 - WGS84_ECC_SQ * sin_lat * sin_lat).sqrt();
        let new_lat = (z + WGS84_ECC_SQ * n * sin_lat).atan2(p);
        if (new_lat - lat).abs() < 1e-14 {
            lat = new_lat;
            break;
        }
        lat = new_lat;
    }

    let sin_lat = lat.sin();
    let n = a / (1.0 - WGS84_ECC_SQ * sin_lat * sin_lat).sqrt();
    // Altitude: away from the poles use p/cosφ − N; near the poles fall back to
    // the z-based form to avoid dividing by a tiny cosφ.
    let cos_lat = lat.cos();
    let alt = if cos_lat.abs() > 1e-3 {
        p / cos_lat - n
    } else {
        z.abs() / sin_lat.abs() - n * (1.0 - WGS84_ECC_SQ)
    };
    (lat, lon, alt)
}

/// Rotate an **ECI** position vector (m) into the **ECEF** frame, given the
/// Greenwich sidereal angle `gmst_rad` (e.g. from [`gmst`]).
///
/// The Earth-fixed frame is the inertial frame turned by `+θ` about the pole, so
/// a fixed inertial vector expressed in Earth-fixed coordinates is rotated by
/// `−θ`:
///
/// ```text
/// x_ecef =  cosθ·x_eci + sinθ·y_eci
/// y_ecef = −sinθ·x_eci + cosθ·y_eci
/// z_ecef =  z_eci
/// ```
///
/// This matches the rotation used by [`crate::groundtrack::subpoint`]. The
/// inverse is [`ecef_to_eci`].
pub fn eci_to_ecef(eci: Vector3<f64>, gmst_rad: f64) -> Vector3<f64> {
    let (s, c) = gmst_rad.sin_cos();
    Vector3::new(c * eci.x + s * eci.y, -s * eci.x + c * eci.y, eci.z)
}

/// Rotate an **ECEF** position vector (m) into the **ECI** frame, given the
/// Greenwich sidereal angle `gmst_rad` — the inverse of [`eci_to_ecef`].
pub fn ecef_to_eci(ecef: Vector3<f64>, gmst_rad: f64) -> Vector3<f64> {
    let (s, c) = gmst_rad.sin_cos();
    Vector3::new(c * ecef.x - s * ecef.y, s * ecef.x + c * ecef.y, ecef.z)
}

/// The Earth-fixed position of a point fixed in the inertial frame drifts west
/// at the sidereal rate; this returns the GMST angle at `t` seconds after an
/// epoch whose GMST is `gmst0` — `θ(t) = θ₀ + ω⊕·t`, wrapped to `[0, 2π)`.
///
/// A convenience for stepping ground geometry through a pass without
/// re-evaluating the full [`gmst`] polynomial at every sample (the secular
/// `ω⊕·t` term dominates over a single pass).
pub fn gmst_after(gmst0: f64, t_seconds: f64) -> f64 {
    let two_pi = std::f64::consts::TAU;
    let mut theta = (gmst0 + OMEGA_EARTH * t_seconds) % two_pi;
    if theta < 0.0 {
        theta += two_pi;
    }
    theta
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::R_EARTH;

    #[test]
    fn geodetic_ecef_round_trips() {
        // A spread of stations: equator, mid-latitude, high-latitude, with and
        // without altitude, both hemispheres and both longitudes.
        let cases = [
            (0.0_f64, 0.0_f64, 0.0_f64),
            (45.0, -75.0, 100.0),  // Ottawa-ish
            (-33.9, 151.2, 58.0),  // Sydney-ish
            (51.48, 0.0, 47.0),    // Greenwich
            (78.0, 15.0, 2_000.0), // Svalbard-ish, high alt
            (-45.0, -170.0, 0.0),  // southern, far east-of-dateline
            (10.0, 179.9, 500.0),  // near the antimeridian
        ];
        for (lat_deg, lon_deg, alt) in cases {
            let lat = lat_deg.to_radians();
            let lon = lon_deg.to_radians();
            let ecef = geodetic_to_ecef(lat, lon, alt);
            let (lat2, lon2, alt2) = ecef_to_geodetic(ecef);
            assert!(
                (lat2 - lat).abs() < 1e-9,
                "lat {lat_deg}: {} != {}",
                lat2.to_degrees(),
                lat_deg
            );
            assert!(
                (lon2 - lon).abs() < 1e-9,
                "lon {lon_deg}: {} != {}",
                lon2.to_degrees(),
                lon_deg
            );
            assert!((alt2 - alt).abs() < 1e-6, "alt {alt}: {alt2} != {alt}");
        }
    }

    #[test]
    fn equator_prime_meridian_is_on_the_equatorial_radius() {
        // (0,0,0) geodetic sits at exactly the WGS-84 equatorial radius on +x.
        let ecef = geodetic_to_ecef(0.0, 0.0, 0.0);
        assert!((ecef.x - R_EARTH).abs() < 1e-6, "x = {}", ecef.x);
        assert!(ecef.y.abs() < 1e-6 && ecef.z.abs() < 1e-6);
    }

    #[test]
    fn pole_is_polar_radius_and_no_nan() {
        // The geographic North Pole at sea level: ECEF should be (0,0,b) with
        // b the polar radius, and the inverse must return +90° latitude (not a
        // NaN from the p = 0 singularity).
        let b = R_EARTH * (1.0 - WGS84_ECC_SQ).sqrt();
        let ecef = geodetic_to_ecef(std::f64::consts::FRAC_PI_2, 0.3, 0.0);
        assert!(ecef.x.abs() < 1e-6 && ecef.y.abs() < 1e-6);
        assert!((ecef.z - b).abs() < 1e-6, "z {} != b {}", ecef.z, b);
        let (lat, _, alt) = ecef_to_geodetic(Vector3::new(0.0, 0.0, b));
        assert!(
            (lat.to_degrees() - 90.0).abs() < 1e-9,
            "lat {}",
            lat.to_degrees()
        );
        assert!(alt.abs() < 1e-6, "alt {alt}");
        // The geocentre maps to all-zeros, no NaN.
        let (lat0, lon0, alt0) = ecef_to_geodetic(Vector3::zeros());
        assert!(lat0 == 0.0 && lon0 == 0.0);
        assert!((alt0).abs() < 1.0 || alt0.is_finite());
    }

    #[test]
    fn gmst_at_j2000_matches_the_known_value() {
        // GROUND TRUTH: GMST at J2000.0 (2000-01-01 12:00 UT1) is
        // 18h 41m 50.548s ≈ 280.4606°. (Vallado example / IAU 1982 series.)
        let theta = gmst(J2000).to_degrees();
        assert!(
            (theta - 280.4606).abs() < 1e-3,
            "GMST(J2000) = {theta}° != 280.4606°"
        );
    }

    #[test]
    fn gmst_advances_one_sidereal_turn_per_sidereal_day() {
        // Over one mean solar day (86400 s of UT) GMST advances slightly more
        // than a full turn (~360.9856°), the well-known sidereal-vs-solar gap.
        let two_pi = std::f64::consts::TAU;
        let g0 = gmst(J2000);
        let g1 = gmst(J2000 + 1.0); // +1 day
        let mut advance = (g1 - g0).rem_euclid(two_pi).to_degrees();
        // It went round once plus the gap; recover the total advance.
        advance += 360.0;
        assert!(
            (advance - 360.9856).abs() < 1e-2,
            "one-day GMST advance {advance}° != 360.9856°"
        );
    }

    #[test]
    fn eci_ecef_rotation_is_an_inverse_pair_and_preserves_length() {
        let v = Vector3::new(7.0e6, -2.0e6, 3.0e6);
        for &theta in &[0.0, 0.3, 1.7, 4.5, 6.0] {
            let ecef = eci_to_ecef(v, theta);
            assert!(
                (ecef.norm() - v.norm()).abs() < 1e-6,
                "rotation must be isometric"
            );
            let back = ecef_to_eci(ecef, theta);
            assert!((back - v).norm() < 1e-6, "round-trip at θ={theta}");
        }
        // At θ = 0 the frames coincide.
        assert!((eci_to_ecef(v, 0.0) - v).norm() < 1e-9);
    }

    #[test]
    fn gmst_after_matches_full_gmst_over_a_short_span() {
        // Over a single pass the secular ω⊕·t advance reproduces the full
        // polynomial to high accuracy.
        let jd0 = J2000 + 1234.5;
        let g0 = gmst(jd0);
        let dt = 600.0; // 10 minutes
        let approx = gmst_after(g0, dt);
        let exact = gmst(jd0 + dt / 86_400.0);
        let two_pi = std::f64::consts::TAU;
        let diff = (approx - exact).rem_euclid(two_pi);
        let diff = diff.min(two_pi - diff);
        assert!(diff < 1e-6, "gmst_after vs gmst diff {diff} rad");
    }
}
