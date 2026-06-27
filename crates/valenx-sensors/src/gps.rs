//! Simulated **GPS** (GNSS) position sensor.
//!
//! A [`Gps`] turns a vehicle's **local ENU** position — east / north / up metres
//! relative to a fixed geodetic **datum** (the origin of the local tangent
//! plane) — into a geodetic fix `(latitude, longitude, altitude)`, with optional
//! additive position noise. This is the inverse of the "flat-Earth local frame"
//! that the autonomy harness integrates in.
//!
//! The conversion goes ENU → ECEF → geodetic on the WGS-84 ellipsoid:
//!
//! 1. The datum's geodetic `(lat₀, lon₀, alt₀)` fixes an ECEF origin and the
//!    local east/north/up basis (the rows of the standard ENU rotation).
//! 2. A local offset `(e, n, u)` becomes an ECEF point
//!    `r = r₀ + e·ê + n·n̂ + u·û`.
//! 3. That ECEF point is converted back to geodetic by Bowring's method.
//!
//! Noise (`horizontal_std`, `vertical_std`, metres) is added in the local ENU
//! frame before the conversion, drawn from the crate's seeded
//! [`crate::SplitMix64`], so a run reproduces. With both standard deviations
//! zero the fix is the exact geodetic position of the ENU offset, which is what
//! the tests pin.
//!
//! Honest scope: this is a **geometry-grade** position source — it models *where*
//! a perfect receiver would place a point, plus a simple Gaussian error. It does
//! **not** model GNSS signal physics: no satellite geometry / dilution of
//! precision, no multipath, no ionospheric / tropospheric delay, no clock bias,
//! and no correlated/biased error process. The WGS-84 conversion itself is exact
//! to floating point.

use nalgebra::Vector3;

use crate::error::SensorError;
use crate::rng::SplitMix64;

/// WGS-84 semi-major axis (m).
const WGS84_A: f64 = 6_378_137.0;
/// WGS-84 first eccentricity squared.
const WGS84_E2: f64 = 6.694_379_990_141_316e-3;

/// A geodetic position on the WGS-84 ellipsoid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Geodetic {
    /// Latitude (radians, geodetic).
    pub lat: f64,
    /// Longitude (radians, east-positive).
    pub lon: f64,
    /// Altitude above the ellipsoid (m).
    pub alt: f64,
}

impl Geodetic {
    /// Construct from **degrees** of latitude / longitude and metres of
    /// altitude.
    #[must_use]
    pub fn from_degrees(lat_deg: f64, lon_deg: f64, alt_m: f64) -> Self {
        Self {
            lat: lat_deg.to_radians(),
            lon: lon_deg.to_radians(),
            alt: alt_m,
        }
    }

    /// Latitude in degrees.
    #[must_use]
    pub fn lat_deg(&self) -> f64 {
        self.lat.to_degrees()
    }

    /// Longitude in degrees.
    #[must_use]
    pub fn lon_deg(&self) -> f64 {
        self.lon.to_degrees()
    }
}

/// Convert geodetic `(lat, lon, alt)` to an ECEF position (m).
fn geodetic_to_ecef(g: &Geodetic) -> Vector3<f64> {
    let (sin_lat, cos_lat) = g.lat.sin_cos();
    let (sin_lon, cos_lon) = g.lon.sin_cos();
    let n = WGS84_A / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();
    Vector3::new(
        (n + g.alt) * cos_lat * cos_lon,
        (n + g.alt) * cos_lat * sin_lon,
        (n * (1.0 - WGS84_E2) + g.alt) * sin_lat,
    )
}

/// Convert an ECEF position (m) to geodetic, via Bowring's iteration.
fn ecef_to_geodetic(r: Vector3<f64>) -> Geodetic {
    let (x, y, z) = (r.x, r.y, r.z);
    let p = (x * x + y * y).sqrt();

    // Polar-axis guard (p ≈ 0): longitude undefined, latitude ±90°.
    if p < 1e-9 {
        let b = WGS84_A * (1.0 - WGS84_E2).sqrt();
        if z.abs() < 1e-9 {
            return Geodetic {
                lat: 0.0,
                lon: 0.0,
                alt: -b,
            };
        }
        return Geodetic {
            lat: std::f64::consts::FRAC_PI_2 * z.signum(),
            lon: 0.0,
            alt: z.abs() - b,
        };
    }

    let lon = y.atan2(x);
    let a = WGS84_A;
    let b = a * (1.0 - WGS84_E2).sqrt();
    let ep2 = (a * a - b * b) / (b * b);

    let theta = (z * a).atan2(p * b);
    let (st, ct) = theta.sin_cos();
    let mut lat = (z + ep2 * b * st * st * st).atan2(p - WGS84_E2 * a * ct * ct * ct);

    for _ in 0..5 {
        let sin_lat = lat.sin();
        let n = a / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();
        let new_lat = (z + WGS84_E2 * n * sin_lat).atan2(p);
        if (new_lat - lat).abs() < 1e-14 {
            lat = new_lat;
            break;
        }
        lat = new_lat;
    }

    let sin_lat = lat.sin();
    let n = a / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();
    let cos_lat = lat.cos();
    let alt = if cos_lat.abs() > 1e-3 {
        p / cos_lat - n
    } else {
        z.abs() / sin_lat.abs() - n * (1.0 - WGS84_E2)
    };
    Geodetic { lat, lon, alt }
}

/// A simulated GPS receiver referenced to a fixed geodetic datum.
#[derive(Debug, Clone)]
pub struct Gps {
    datum: Geodetic,
    datum_ecef: Vector3<f64>,
    /// Local ENU basis as columns: [east | north | up].
    east: Vector3<f64>,
    north: Vector3<f64>,
    up: Vector3<f64>,
    horizontal_std: f64,
    vertical_std: f64,
    rng: SplitMix64,
}

/// A GPS fix: the geodetic position plus the (noisy) ENU offset it came from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpsFix {
    /// The reported geodetic position.
    pub position: Geodetic,
    /// The (noise-perturbed) local ENU offset that produced it (m).
    pub enu: Vector3<f64>,
}

impl Gps {
    /// Build a GPS at a geodetic `datum` with horizontal/vertical noise standard
    /// deviations (m) and a noise seed.
    ///
    /// # Errors
    /// - [`SensorError::InvalidConfig`] if the datum latitude is outside
    ///   `[−π/2, π/2]` or any datum component is non-finite.
    /// - [`SensorError::InvalidNoise`] if either standard deviation is negative
    ///   or non-finite.
    pub fn new(
        datum: Geodetic,
        horizontal_std: f64,
        vertical_std: f64,
        seed: u64,
    ) -> Result<Self, SensorError> {
        if !(datum.lat.is_finite() && datum.lon.is_finite() && datum.alt.is_finite()) {
            return Err(SensorError::InvalidConfig("datum must be finite".into()));
        }
        if datum.lat.abs() > std::f64::consts::FRAC_PI_2 + 1e-9 {
            return Err(SensorError::InvalidConfig(format!(
                "datum latitude must be in [−π/2, π/2], got {} rad",
                datum.lat
            )));
        }
        if !(horizontal_std.is_finite() && horizontal_std >= 0.0) {
            return Err(SensorError::InvalidNoise(format!(
                "horizontal_std must be finite and ≥ 0, got {horizontal_std}"
            )));
        }
        if !(vertical_std.is_finite() && vertical_std >= 0.0) {
            return Err(SensorError::InvalidNoise(format!(
                "vertical_std must be finite and ≥ 0, got {vertical_std}"
            )));
        }

        let (sin_lat, cos_lat) = datum.lat.sin_cos();
        let (sin_lon, cos_lon) = datum.lon.sin_cos();
        // Standard local ENU basis expressed in ECEF.
        let east = Vector3::new(-sin_lon, cos_lon, 0.0);
        let north = Vector3::new(-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat);
        let up = Vector3::new(cos_lat * cos_lon, cos_lat * sin_lon, sin_lat);

        Ok(Self {
            datum,
            datum_ecef: geodetic_to_ecef(&datum),
            east,
            north,
            up,
            horizontal_std,
            vertical_std,
            rng: SplitMix64::new(seed),
        })
    }

    /// The reference datum.
    #[must_use]
    pub fn datum(&self) -> Geodetic {
        self.datum
    }

    /// The **exact** geodetic position of a local ENU offset (no noise).
    #[must_use]
    pub fn enu_to_geodetic(&self, enu: Vector3<f64>) -> Geodetic {
        let ecef = self.datum_ecef + self.east * enu.x + self.north * enu.y + self.up * enu.z;
        ecef_to_geodetic(ecef)
    }

    /// Sample a fix at local ENU offset `enu` (east, north, up; m), adding
    /// seeded Gaussian noise (horizontal on east+north, vertical on up).
    pub fn sample(&mut self, enu: Vector3<f64>) -> GpsFix {
        let noisy = Vector3::new(
            enu.x + self.rng.next_normal(0.0, self.horizontal_std),
            enu.y + self.rng.next_normal(0.0, self.horizontal_std),
            enu.z + self.rng.next_normal(0.0, self.vertical_std),
        );
        GpsFix {
            position: self.enu_to_geodetic(noisy),
            enu: noisy,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    fn datum() -> Geodetic {
        // ~Seattle-ish, an arbitrary mid-latitude datum.
        Geodetic::from_degrees(47.6062, -122.3321, 50.0)
    }

    #[test]
    fn zero_offset_returns_the_datum() {
        let gps = Gps::new(datum(), 0.0, 0.0, 0).unwrap();
        let g = gps.enu_to_geodetic(v(0.0, 0.0, 0.0));
        assert!((g.lat - datum().lat).abs() < 1e-12, "lat");
        assert!((g.lon - datum().lon).abs() < 1e-12, "lon");
        assert!((g.alt - datum().alt).abs() < 1e-6, "alt = {}", g.alt);
    }

    #[test]
    fn pure_up_offset_only_raises_altitude() {
        // 100 m straight up: altitude +100, lat/lon (essentially) unchanged.
        let gps = Gps::new(datum(), 0.0, 0.0, 0).unwrap();
        let g = gps.enu_to_geodetic(v(0.0, 0.0, 100.0));
        assert!(
            (g.alt - (datum().alt + 100.0)).abs() < 1e-4,
            "alt = {}",
            g.alt
        );
        assert!((g.lat - datum().lat).abs() < 1e-9, "lat moved");
        assert!((g.lon - datum().lon).abs() < 1e-9, "lon moved");
    }

    #[test]
    fn east_offset_increases_longitude_north_increases_latitude() {
        let gps = Gps::new(datum(), 0.0, 0.0, 0).unwrap();
        let east = gps.enu_to_geodetic(v(500.0, 0.0, 0.0));
        let north = gps.enu_to_geodetic(v(0.0, 500.0, 0.0));
        // East ⇒ more longitude (we're at negative lon, so it grows toward 0).
        assert!(east.lon > datum().lon, "east should raise lon");
        assert!(
            (east.lat - datum().lat).abs() < 1e-7,
            "east shouldn't move lat much"
        );
        // North ⇒ more latitude.
        assert!(north.lat > datum().lat, "north should raise lat");
    }

    #[test]
    fn small_offset_matches_metres_per_degree_ballpark() {
        // 1 km north should change latitude by ~1000 / 111_320 deg ≈ 0.00898°.
        let gps = Gps::new(datum(), 0.0, 0.0, 0).unwrap();
        let g = gps.enu_to_geodetic(v(0.0, 1_000.0, 0.0));
        let dlat_deg = g.lat_deg() - datum().lat_deg();
        // Meridian length per degree near 47.6° is ~111.4 km; allow a band.
        assert!((0.0085..0.0093).contains(&dlat_deg), "Δlat = {dlat_deg}°");
    }

    #[test]
    fn roundtrip_geodetic_ecef_is_exact() {
        // ENU→ECEF→geodetic must recover a hand-built offset's geodetic point;
        // check the ECEF round-trip closes at the datum.
        let g = datum();
        let back = ecef_to_geodetic(geodetic_to_ecef(&g));
        assert!((back.lat - g.lat).abs() < 1e-12);
        assert!((back.lon - g.lon).abs() < 1e-12);
        assert!((back.alt - g.alt).abs() < 1e-6, "alt = {}", back.alt);
    }

    #[test]
    fn noise_is_deterministic_and_perturbs_the_offset() {
        let mut a = Gps::new(datum(), 2.0, 4.0, 7).unwrap();
        let mut b = Gps::new(datum(), 2.0, 4.0, 7).unwrap();
        let fa = a.sample(v(10.0, 20.0, 0.0));
        let fb = b.sample(v(10.0, 20.0, 0.0));
        assert_eq!(fa, fb, "same seed ⇒ identical fix");
        // The noisy ENU is near but not equal to the requested offset.
        assert!((fa.enu.x - 10.0).abs() < 10.0 && fa.enu.x != 10.0);
    }

    #[test]
    fn invalid_config_rejected() {
        // Latitude out of range.
        let bad = Geodetic {
            lat: 2.0, // > π/2
            lon: 0.0,
            alt: 0.0,
        };
        assert!(Gps::new(bad, 0.0, 0.0, 0).is_err());
        // Negative noise.
        assert!(Gps::new(datum(), -1.0, 0.0, 0).is_err());
        assert!(Gps::new(datum(), 0.0, f64::NAN, 0).is_err());
    }
}
