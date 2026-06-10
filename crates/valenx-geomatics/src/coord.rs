//! WGS84 ↔ UTM coordinate conversion.
//!
//! v1 uses the standard 4-term Krüger series for the transverse-
//! Mercator forward/inverse, accurate to ≈cm within any 6° UTM
//! zone — plenty for the engineering-surveying use cases this
//! workbench targets. Karney's higher-order series can drop into
//! [`wgs84_to_utm`] / [`utm_to_wgs84`] as a Phase 35.5 follow-up if
//! sub-mm fidelity becomes a need.

use serde::{Deserialize, Serialize};

use crate::error::GeomaticsError;

/// Geographic coordinates on the WGS84 ellipsoid.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct LatLon {
    /// Latitude in decimal degrees, [-90, 90].
    pub latitude_deg: f64,
    /// Longitude in decimal degrees, [-180, 180].
    pub longitude_deg: f64,
    /// Elevation above the ellipsoid (m). Kept around so round-
    /// tripping through UTM doesn't lose it.
    pub elevation_m: f64,
}

/// Northern or southern hemisphere flag for a UTM coord (the false-
/// northing offset differs).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Hemisphere {
    /// Northern hemisphere (false northing = 0).
    North,
    /// Southern hemisphere (false northing = 10 000 000 m).
    South,
}

/// Universal Transverse Mercator coordinates.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Utm {
    /// Zone number 1..=60.
    pub zone: u32,
    /// N or S.
    pub hemisphere: Hemisphere,
    /// Easting (m) from the zone's central meridian + 500 000 false.
    pub easting_m: f64,
    /// Northing (m).
    pub northing_m: f64,
    /// Elevation above the ellipsoid (m) — copied across from
    /// LatLon so round-trip preserves it.
    pub elevation_m: f64,
}

// WGS84 ellipsoid constants.
const A: f64 = 6_378_137.0;
const F_INV: f64 = 298.257_223_563;
const K0: f64 = 0.999_6;
const FALSE_EASTING: f64 = 500_000.0;
const FALSE_NORTHING_S: f64 = 10_000_000.0;

/// Convert WGS84 (lat, lon) to UTM (zone, hem, E, N) using the
/// 4-term Krüger series.
pub fn wgs84_to_utm(p: LatLon) -> Result<Utm, GeomaticsError> {
    if !(-80.0..=84.0).contains(&p.latitude_deg) {
        return Err(GeomaticsError::BadParameter {
            name: "latitude_deg",
            reason: format!("must be in [-80, 84] for UTM, got {}", p.latitude_deg),
        });
    }
    let zone = utm_zone_for_longitude(p.longitude_deg);
    let central_meridian = (zone as f64 - 1.0) * 6.0 - 180.0 + 3.0;
    let lat = p.latitude_deg.to_radians();
    let lon = p.longitude_deg.to_radians();
    let lon0 = central_meridian.to_radians();

    let f = 1.0 / F_INV;
    let n = f / (2.0 - f);
    let n2 = n * n;
    let n3 = n2 * n;
    let n4 = n3 * n;

    // Meridional arc terms.
    let a_hat = A / (1.0 + n) * (1.0 + n2 / 4.0 + n4 / 64.0);

    // Conformal latitude.
    let e2 = f * (2.0 - f);
    let e = e2.sqrt();
    let sin_lat = lat.sin();
    let tau = sin_lat / lat.cos();
    // Karney's conformal-latitude step: σ = sinh(e·atanh(e·τ/√(1+τ²)))
    // and τ/√(1+τ²) = sinφ, so the atanh argument is e·sinφ. (The
    // inverse below uses the same e·τ/√(1+τ²) form — the two must
    // match, and this is the analytically-exact conformal latitude.)
    let sigma = (e * (e * sin_lat).atanh()).sinh();
    let tau_prime = tau * (1.0 + sigma * sigma).sqrt() - sigma * (1.0 + tau * tau).sqrt();
    let xi_prime = tau_prime.atan2((lon - lon0).cos());
    let eta_prime = (((lon - lon0).sin()) / (tau_prime * tau_prime + (lon - lon0).cos().powi(2)).sqrt()).asinh();

    // Krüger series coefficients (alpha_1 .. alpha_4).
    let alpha_1 = 0.5 * n - 2.0 / 3.0 * n2 + 5.0 / 16.0 * n3 + 41.0 / 180.0 * n4;
    let alpha_2 = 13.0 / 48.0 * n2 - 0.6 * n3 + 557.0 / 1440.0 * n4;
    let alpha_3 = 61.0 / 240.0 * n3 - 103.0 / 140.0 * n4;
    let alpha_4 = 49561.0 / 161_280.0 * n4;

    let xi = xi_prime
        + alpha_1 * (2.0 * xi_prime).sin() * (2.0 * eta_prime).cosh()
        + alpha_2 * (4.0 * xi_prime).sin() * (4.0 * eta_prime).cosh()
        + alpha_3 * (6.0 * xi_prime).sin() * (6.0 * eta_prime).cosh()
        + alpha_4 * (8.0 * xi_prime).sin() * (8.0 * eta_prime).cosh();
    let eta = eta_prime
        + alpha_1 * (2.0 * xi_prime).cos() * (2.0 * eta_prime).sinh()
        + alpha_2 * (4.0 * xi_prime).cos() * (4.0 * eta_prime).sinh()
        + alpha_3 * (6.0 * xi_prime).cos() * (6.0 * eta_prime).sinh()
        + alpha_4 * (8.0 * xi_prime).cos() * (8.0 * eta_prime).sinh();

    let hem = if p.latitude_deg >= 0.0 {
        Hemisphere::North
    } else {
        Hemisphere::South
    };
    let easting = K0 * a_hat * eta + FALSE_EASTING;
    let northing_base = K0 * a_hat * xi;
    let northing = if matches!(hem, Hemisphere::South) {
        northing_base + FALSE_NORTHING_S
    } else {
        northing_base
    };
    Ok(Utm {
        zone,
        hemisphere: hem,
        easting_m: easting,
        northing_m: northing,
        elevation_m: p.elevation_m,
    })
}

/// Convert UTM to WGS84 (lat, lon) using the inverse 4-term Krüger
/// series.
pub fn utm_to_wgs84(u: Utm) -> Result<LatLon, GeomaticsError> {
    if u.zone == 0 || u.zone > 60 {
        return Err(GeomaticsError::BadParameter {
            name: "zone",
            reason: format!("must be in 1..=60, got {}", u.zone),
        });
    }
    let central_meridian = (u.zone as f64 - 1.0) * 6.0 - 180.0 + 3.0;
    let lon0 = central_meridian.to_radians();
    let f = 1.0 / F_INV;
    let n = f / (2.0 - f);
    let n2 = n * n;
    let n3 = n2 * n;
    let n4 = n3 * n;
    let a_hat = A / (1.0 + n) * (1.0 + n2 / 4.0 + n4 / 64.0);

    let northing_eff = if matches!(u.hemisphere, Hemisphere::South) {
        u.northing_m - FALSE_NORTHING_S
    } else {
        u.northing_m
    };
    let xi = northing_eff / (a_hat * K0);
    let eta = (u.easting_m - FALSE_EASTING) / (a_hat * K0);

    let beta_1 = 0.5 * n - 2.0 / 3.0 * n2 + 37.0 / 96.0 * n3 - n4 / 360.0;
    let beta_2 = n2 / 48.0 + n3 / 15.0 - 437.0 / 1440.0 * n4;
    let beta_3 = 17.0 / 480.0 * n3 - 37.0 / 840.0 * n4;
    let beta_4 = 4397.0 / 161_280.0 * n4;

    let xi_prime = xi
        - beta_1 * (2.0 * xi).sin() * (2.0 * eta).cosh()
        - beta_2 * (4.0 * xi).sin() * (4.0 * eta).cosh()
        - beta_3 * (6.0 * xi).sin() * (6.0 * eta).cosh()
        - beta_4 * (8.0 * xi).sin() * (8.0 * eta).cosh();
    let eta_prime = eta
        - beta_1 * (2.0 * xi).cos() * (2.0 * eta).sinh()
        - beta_2 * (4.0 * xi).cos() * (4.0 * eta).sinh()
        - beta_3 * (6.0 * xi).cos() * (6.0 * eta).sinh()
        - beta_4 * (8.0 * xi).cos() * (8.0 * eta).sinh();

    let tau_prime = xi_prime.sin() / (eta_prime.sinh().powi(2) + xi_prime.cos().powi(2)).sqrt();
    let e2 = f * (2.0 - f);
    let e = e2.sqrt();
    // Recover the geodetic τ = tanφ from the conformal τ′ by Newton's
    // method on the forward relation F(τ) = τ·√(1+σ²) − σ·√(1+τ²)
    // (Karney 2011). A bare fixed-point on that relation does NOT
    // invert F — the forward map must be Newton-solved for its root.
    let mut tau = tau_prime;
    for _ in 0..8 {
        let sigma = (e * (e * tau / (1.0 + tau * tau).sqrt()).atanh()).sinh();
        let tau_p_i = tau * (1.0 + sigma * sigma).sqrt() - sigma * (1.0 + tau * tau).sqrt();
        let d_tau = (tau_prime - tau_p_i) * (1.0 + (1.0 - e2) * tau * tau)
            / ((1.0 - e2) * (1.0 + tau * tau).sqrt() * (1.0 + tau_p_i * tau_p_i).sqrt());
        tau += d_tau;
        if d_tau.abs() < 1e-14 {
            break;
        }
    }
    let lat = tau.atan();
    let lon = lon0 + eta_prime.sinh().atan2(xi_prime.cos());
    Ok(LatLon {
        latitude_deg: lat.to_degrees(),
        longitude_deg: lon.to_degrees(),
        elevation_m: u.elevation_m,
    })
}

/// UTM zone for a longitude (1..=60). Norway / Svalbard exceptions
/// are not applied in v1.
pub fn utm_zone_for_longitude(lon_deg: f64) -> u32 {
    let z = ((lon_deg + 180.0) / 6.0).floor() as i32 + 1;
    z.clamp(1, 60) as u32
}

/// Great-circle distance (m) between two WGS84 points via the haversine formula,
/// `d = 2R·asin(√(sin²(Δφ/2) + cosφ₁·cosφ₂·sin²(Δλ/2)))`, with `R = 6_371_008.8 m` (the WGS84
/// mean radius). Elevation is ignored. This is the shortest surface path on a sphere; for short
/// spans (<100 km) it is within ~0.1 % of the true ellipsoidal geodesic. Identical points → `0.0`.
pub fn haversine_distance(a: LatLon, b: LatLon) -> f64 {
    const EARTH_RADIUS_M: f64 = 6_371_008.8;
    let lat1 = a.latitude_deg.to_radians();
    let lat2 = b.latitude_deg.to_radians();
    let dlat = (b.latitude_deg - a.latitude_deg).to_radians();
    let dlon = (b.longitude_deg - a.longitude_deg).to_radians();
    let h = (dlat * 0.5).sin().powi(2) + lat1.cos() * lat2.cos() * (dlon * 0.5).sin().powi(2);
    2.0 * EARTH_RADIUS_M * h.sqrt().asin()
}

/// Initial great-circle bearing (forward azimuth) from `a` to `b`, in degrees clockwise from
/// north and normalised to `[0, 360)`:
/// `θ = atan2(sinΔλ·cosφ₂, cosφ₁·sinφ₂ − sinφ₁·cosφ₂·cosΔλ)`. Distinct from
/// [`haversine_distance`] (a distance, not a heading). Identical points → `0.0`.
pub fn initial_bearing(a: LatLon, b: LatLon) -> f64 {
    let lat1 = a.latitude_deg.to_radians();
    let lat2 = b.latitude_deg.to_radians();
    let dlon = (b.longitude_deg - a.longitude_deg).to_radians();
    let y = dlon.sin() * lat2.cos();
    let x = lat1.cos() * lat2.sin() - lat1.sin() * lat2.cos() * dlon.cos();
    (y.atan2(x).to_degrees() + 360.0) % 360.0
}

/// Cross-track distance (m) — the signed perpendicular offset of point `p` from the great-circle
/// path `path_start` → `path_end`: `d_xt = asin(sin(d₁₃/R)·sin(θ₁₃−θ₁₂))·R`, where d₁₃ =
/// [`haversine_distance`]`(path_start, p)`, θ₁₃/θ₁₂ are the [`initial_bearing`]s from `path_start`
/// to `p` and to `path_end`, and R = 6_371_008.8 m. Positive = right of the path, negative = left.
/// Returns `0.0` for a degenerate path (`path_start == path_end`) or non-finite input.
pub fn cross_track_distance(p: LatLon, path_start: LatLon, path_end: LatLon) -> f64 {
    const EARTH_RADIUS_M: f64 = 6_371_008.8;
    let degenerate = path_start.latitude_deg == path_end.latitude_deg
        && path_start.longitude_deg == path_end.longitude_deg;
    let all_finite = p.latitude_deg.is_finite()
        && p.longitude_deg.is_finite()
        && path_start.latitude_deg.is_finite()
        && path_start.longitude_deg.is_finite()
        && path_end.latitude_deg.is_finite()
        && path_end.longitude_deg.is_finite();
    if degenerate || !all_finite {
        return 0.0;
    }
    let d13 = haversine_distance(path_start, p);
    let theta13 = initial_bearing(path_start, p).to_radians();
    let theta12 = initial_bearing(path_start, path_end).to_radians();
    let arg = (d13 / EARTH_RADIUS_M).sin() * (theta13 - theta12).sin();
    arg.clamp(-1.0, 1.0).asin() * EARTH_RADIUS_M
}

/// Along-track distance (m) — the distance from `path_start` to the foot of the perpendicular
/// dropped from `p` onto the great-circle path `path_start` → `path_end`: the parallel projection
/// of `p` along the path, complementing [`cross_track_distance`] (the perpendicular offset).
/// `d_at = acos(cos(d₁₃/R)/cos(d_xt/R))·R`, with d₁₃ = [`haversine_distance`]`(path_start, p)`,
/// d_xt = [`cross_track_distance`]`(p, path_start, path_end)`, R = 6_371_008.8 m. Returns `0.0` for
/// a degenerate path (`path_start == path_end`) or non-finite input.
pub fn along_track_distance(p: LatLon, path_start: LatLon, path_end: LatLon) -> f64 {
    const EARTH_RADIUS_M: f64 = 6_371_008.8;
    let degenerate = path_start.latitude_deg == path_end.latitude_deg
        && path_start.longitude_deg == path_end.longitude_deg;
    let all_finite = p.latitude_deg.is_finite()
        && p.longitude_deg.is_finite()
        && path_start.latitude_deg.is_finite()
        && path_start.longitude_deg.is_finite()
        && path_end.latitude_deg.is_finite()
        && path_end.longitude_deg.is_finite();
    if degenerate || !all_finite {
        return 0.0;
    }
    let d13 = haversine_distance(path_start, p);
    let d_xt = cross_track_distance(p, path_start, path_end);
    let arg = (d13 / EARTH_RADIUS_M).cos() / (d_xt / EARTH_RADIUS_M).cos();
    arg.clamp(-1.0, 1.0).acos() * EARTH_RADIUS_M
}

/// Final bearing (the arrival/back azimuth) at `b` when travelling the great circle from `a`, in
/// degrees clockwise from north and normalised to `[0, 360)`. Computed as
/// `(initial_bearing(b, a) + 180) mod 360`. Distinct from [`initial_bearing`] (the departure
/// heading): on a sphere the arrival heading differs from the departure heading by ≠ 180° in
/// general because the meridians converge (they coincide only on the equator or a meridian).
pub fn final_bearing(a: LatLon, b: LatLon) -> f64 {
    (initial_bearing(b, a) + 180.0) % 360.0
}

/// Rhumb-line (loxodrome) distance (m) between two WGS84 points — the length of the path of
/// constant bearing, via the Mercator loxodrome `d = √(Δφ² + q²·Δλ²)·R`, R = 6_371_008.8 m.
/// Distinct from [`haversine_distance`] (the great circle): the rhumb line is never shorter, and
/// is equal only along a meridian or the equator. Returns `0.0` for identical or non-finite points.
pub fn rhumb_distance(a: LatLon, b: LatLon) -> f64 {
    use std::f64::consts::PI;
    const EARTH_RADIUS_M: f64 = 6_371_008.8;
    let finite = a.latitude_deg.is_finite()
        && a.longitude_deg.is_finite()
        && b.latitude_deg.is_finite()
        && b.longitude_deg.is_finite();
    if !finite || (a.latitude_deg == b.latitude_deg && a.longitude_deg == b.longitude_deg) {
        return 0.0;
    }
    let lat1 = a.latitude_deg.to_radians();
    let lat2 = b.latitude_deg.to_radians();
    let dlat = lat2 - lat1;
    let mut dlon = (b.longitude_deg - a.longitude_deg).to_radians();
    if dlon.abs() > PI {
        dlon -= dlon.signum() * 2.0 * PI;
    }
    // Mercator stretched-latitude difference; for an E–W course (Δψ → 0) fall back to cos φ₁.
    let dpsi = ((lat2 / 2.0 + PI / 4.0).tan() / (lat1 / 2.0 + PI / 4.0).tan()).ln();
    let q = if dpsi.abs() > 1e-12 { dlat / dpsi } else { lat1.cos() };
    (dlat * dlat + q * q * dlon * dlon).sqrt() * EARTH_RADIUS_M
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zone_for_londons_prime_meridian() {
        assert_eq!(utm_zone_for_longitude(0.0), 31);
    }

    #[test]
    fn round_trip_near_equator() {
        let p = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 3.0,
            elevation_m: 0.0,
        };
        let u = wgs84_to_utm(p).unwrap();
        assert_eq!(u.zone, 31);
        let back = utm_to_wgs84(u).unwrap();
        assert!((back.latitude_deg - p.latitude_deg).abs() < 1e-6);
        assert!((back.longitude_deg - p.longitude_deg).abs() < 1e-6);
    }

    #[test]
    fn round_trip_mid_lat_north() {
        let p = LatLon {
            latitude_deg: 51.5,
            longitude_deg: -0.13,
            elevation_m: 12.0,
        };
        let u = wgs84_to_utm(p).unwrap();
        assert_eq!(u.hemisphere, Hemisphere::North);
        let back = utm_to_wgs84(u).unwrap();
        assert!((back.latitude_deg - p.latitude_deg).abs() < 1e-5);
        assert!((back.longitude_deg - p.longitude_deg).abs() < 1e-5);
        assert_eq!(back.elevation_m, p.elevation_m);
    }

    #[test]
    fn south_hemisphere_has_false_northing() {
        let p = LatLon {
            latitude_deg: -33.86,
            longitude_deg: 151.21,
            elevation_m: 0.0,
        };
        let u = wgs84_to_utm(p).unwrap();
        assert_eq!(u.hemisphere, Hemisphere::South);
        assert!(u.northing_m > 6_000_000.0);
    }

    #[test]
    fn bad_zone_errors() {
        let bad = Utm {
            zone: 0,
            hemisphere: Hemisphere::North,
            easting_m: 0.0,
            northing_m: 0.0,
            elevation_m: 0.0,
        };
        assert!(matches!(
            utm_to_wgs84(bad),
            Err(GeomaticsError::BadParameter { .. })
        ));
    }

    #[test]
    fn polar_lat_errors() {
        let p = LatLon {
            latitude_deg: 85.0,
            longitude_deg: 0.0,
            elevation_m: 0.0,
        };
        assert!(matches!(
            wgs84_to_utm(p),
            Err(GeomaticsError::BadParameter { .. })
        ));
    }

    #[test]
    fn haversine_great_circle_distance() {
        let london = LatLon {
            latitude_deg: 51.5074,
            longitude_deg: -0.1278,
            elevation_m: 0.0,
        };
        let paris = LatLon {
            latitude_deg: 48.8566,
            longitude_deg: 2.3522,
            elevation_m: 0.0,
        };
        // London–Paris ≈ 343.6 km.
        assert!((haversine_distance(london, paris) - 343_600.0).abs() < 1500.0);
        // Symmetric.
        assert!(
            (haversine_distance(london, paris) - haversine_distance(paris, london)).abs() < 1e-6
        );
        // Identical point → 0.
        assert_eq!(haversine_distance(london, london), 0.0);
        // Antipodal points → half the circumference πR.
        let p0 = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 0.0,
            elevation_m: 0.0,
        };
        let p180 = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 180.0,
            elevation_m: 0.0,
        };
        assert!((haversine_distance(p0, p180) - std::f64::consts::PI * 6_371_008.8).abs() < 1.0);
    }

    #[test]
    fn initial_bearing_forward_azimuth() {
        let origin = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 0.0,
            elevation_m: 0.0,
        };
        // Due east along the equator → 90°; due west → 270°.
        let east = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 90.0,
            elevation_m: 0.0,
        };
        let west = LatLon {
            latitude_deg: 0.0,
            longitude_deg: -90.0,
            elevation_m: 0.0,
        };
        assert!((initial_bearing(origin, east) - 90.0).abs() < 1e-6);
        assert!((initial_bearing(origin, west) - 270.0).abs() < 1e-6);
        // Along a meridian the forward/back bearings are exactly 0° (north) and 180° (south).
        let up = LatLon {
            latitude_deg: 10.0,
            longitude_deg: 0.0,
            elevation_m: 0.0,
        };
        assert!(initial_bearing(origin, up).abs() < 1e-9);
        assert!((initial_bearing(up, origin) - 180.0).abs() < 1e-9);
        // Identical points → 0.0 (no NaN); result always in [0, 360).
        assert_eq!(initial_bearing(origin, origin), 0.0);
        assert!((0.0..360.0).contains(&initial_bearing(origin, east)));
    }

    #[test]
    fn cross_track_distance_perpendicular_offset() {
        // Equator path (0,0)→(0,10); point (1,5) lies ~1° of latitude north of it.
        let start = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 0.0,
            elevation_m: 0.0,
        };
        let end = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 10.0,
            elevation_m: 0.0,
        };
        let off = LatLon {
            latitude_deg: 1.0,
            longitude_deg: 5.0,
            elevation_m: 0.0,
        };
        // ~1° of latitude in metres ≈ 111_195 m.
        assert!((cross_track_distance(off, start, end).abs() - 111_195.0).abs() < 200.0);
        // A point on the path → ~0.
        let on = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 5.0,
            elevation_m: 0.0,
        };
        assert!(cross_track_distance(on, start, end).abs() < 1.0);
        // Degenerate path (start == end) → 0.0.
        assert_eq!(cross_track_distance(off, start, start), 0.0);
    }

    #[test]
    fn along_track_distance_parallel_projection() {
        // Equator path (0,0)→(0,10); point (1,5)'s foot of perpendicular ≈ (0,5).
        let start = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 0.0,
            elevation_m: 0.0,
        };
        let end = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 10.0,
            elevation_m: 0.0,
        };
        let off = LatLon {
            latitude_deg: 1.0,
            longitude_deg: 5.0,
            elevation_m: 0.0,
        };
        // Along-track ≈ haversine to (0,5) ≈ 5° longitude at the equator ≈ 555_975 m.
        assert!((along_track_distance(off, start, end) - 555_975.0).abs() < 500.0);
        // A point at path_start → 0.
        assert!(along_track_distance(start, start, end).abs() < 1.0);
        // Degenerate path → 0.0.
        assert_eq!(along_track_distance(off, start, start), 0.0);
    }

    #[test]
    fn final_bearing_is_arrival_heading() {
        let origin = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 0.0,
            elevation_m: 0.0,
        };
        // Meridian (0,0)→(10,0): arrive heading due north → 0°.
        let north = LatLon {
            latitude_deg: 10.0,
            longitude_deg: 0.0,
            elevation_m: 0.0,
        };
        assert!(final_bearing(origin, north).abs() < 1e-9);
        // Equator (0,0)→(0,10): arrive heading due east → 90°.
        let east = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 10.0,
            elevation_m: 0.0,
        };
        assert!((final_bearing(origin, east) - 90.0).abs() < 1e-6);
        // Diagonal: arrival heading differs from the departure heading (meridian convergence).
        let diag = LatLon {
            latitude_deg: 10.0,
            longitude_deg: 10.0,
            elevation_m: 0.0,
        };
        let f = final_bearing(origin, diag);
        assert!((0.0..360.0).contains(&f));
        assert!((f - initial_bearing(origin, diag)).abs() > 0.1);
    }

    #[test]
    fn rhumb_distance_loxodrome() {
        let origin = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 0.0,
            elevation_m: 0.0,
        };
        // Equator (0,0)→(0,10): rhumb = great-circle = 10° of longitude ≈ 1_111_949 m.
        let east = LatLon {
            latitude_deg: 0.0,
            longitude_deg: 10.0,
            elevation_m: 0.0,
        };
        assert!((rhumb_distance(origin, east) - 1_111_949.0).abs() < 100.0);
        // Meridian (0,0)→(10,0): rhumb = great-circle ≈ 1_111_949 m.
        let north = LatLon {
            latitude_deg: 10.0,
            longitude_deg: 0.0,
            elevation_m: 0.0,
        };
        assert!((rhumb_distance(origin, north) - 1_111_949.0).abs() < 100.0);
        // Symmetric; identical → 0.
        assert!((rhumb_distance(origin, north) - rhumb_distance(north, origin)).abs() < 1e-6);
        assert_eq!(rhumb_distance(origin, origin), 0.0);
        // The rhumb line is never shorter than the great circle (a mid-latitude diagonal).
        let a = LatLon {
            latitude_deg: 45.0,
            longitude_deg: 0.0,
            elevation_m: 0.0,
        };
        let b = LatLon {
            latitude_deg: 50.0,
            longitude_deg: 30.0,
            elevation_m: 0.0,
        };
        assert!(rhumb_distance(a, b) >= haversine_distance(a, b) - 1.0);
    }
}
