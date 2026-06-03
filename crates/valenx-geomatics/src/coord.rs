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
}
