//! Ground-station **access / visibility** analysis: when can a station see a
//! satellite, and how high does it climb?
//!
//! Given a satellite **ephemeris** (a time-ordered series of ECI states) and a
//! ground [`GroundStation`] with a minimum-elevation mask, this finds the
//! **access windows** — the rise/set times and the peak-elevation geometry of
//! each pass — the core building block of contact scheduling, downlink
//! planning, sensor tasking and visibility studies (the bread-and-butter of an
//! STK-class space-domain tool).
//!
//! The geometry chain per sample is:
//!
//! 1. Rotate the satellite ECI position into the Earth-fixed frame
//!    ([`crate::frames::eci_to_ecef`]) at that sample's sidereal angle.
//! 2. Form the station→satellite vector in ECEF and rotate it into the
//!    station-local **topocentric SEZ** (South / East / Zenith) frame.
//! 3. The **elevation** is the angle of that vector above the local horizon;
//!    the satellite is *visible* when the elevation is at or above the mask.
//!
//! A contiguous run of visible samples is one access window; its boundaries are
//! refined to the mask-crossing by linear interpolation in the elevation, and
//! the peak elevation is found by a parabolic fit to the local maximum.
//!
//! # Honest scope
//!
//! Accuracy is inherited from the inputs: the ECI states come from the caller's
//! propagator (this crate's two-body / J2 [`crate::orbit3d::propagate`]), and
//! the frame rotation is the mean-sidereal model of [`crate::frames`] (no
//! precession / nutation / polar motion). **Rise/set and peak times are
//! resolved only to the ephemeris sampling**: the linear/parabolic refinement
//! interpolates *between* samples but cannot recover sub-sample structure, so a
//! denser ephemeris gives sharper window edges. There is no atmospheric
//! refraction or terrain-horizon masking — the mask is a single elevation angle
//! over a smooth ellipsoid. This is a visibility-planning product, not an
//! operational pointing/scheduling system.

use nalgebra::Vector3;

use crate::error::AstroError;
use crate::frames::{eci_to_ecef, geodetic_to_ecef, gmst_after};
use crate::orbit3d::StateVector;

/// A ground station / observer fixed to the rotating Earth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GroundStation {
    /// Geodetic latitude (rad), in `[-π/2, π/2]`.
    pub latitude: f64,
    /// East longitude (rad).
    pub longitude: f64,
    /// Height above the WGS-84 ellipsoid (m).
    pub altitude: f64,
    /// Minimum elevation angle for visibility — the *mask* (rad). The
    /// satellite is in view only when its elevation is at or above this.
    /// Typically a few degrees (horizon clutter / multipath); `0` is the
    /// geometric horizon.
    pub min_elevation: f64,
}

impl GroundStation {
    /// Construct a station from degrees (latitude, longitude, mask) plus an
    /// altitude in metres — the convenient human-facing form.
    pub fn from_degrees(lat_deg: f64, lon_deg: f64, alt_m: f64, min_elev_deg: f64) -> Self {
        Self {
            latitude: lat_deg.to_radians(),
            longitude: lon_deg.to_radians(),
            altitude: alt_m,
            min_elevation: min_elev_deg.to_radians(),
        }
    }

    /// The station's fixed position in the Earth-fixed (ECEF) frame (m).
    pub fn ecef(&self) -> Vector3<f64> {
        geodetic_to_ecef(self.latitude, self.longitude, self.altitude)
    }
}

/// One timestamped sample of a satellite ephemeris in the ECI frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EphemerisPoint {
    /// Seconds since the ephemeris epoch.
    pub time: f64,
    /// Inertial (ECI) state at that time.
    pub state: StateVector,
}

/// The instantaneous **look angles** of a satellite from a station: where to
/// point and how far.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LookAngles {
    /// Elevation above the local horizon (rad): `+π/2` straight up (zenith),
    /// `0` on the horizon, negative below it.
    pub elevation: f64,
    /// Azimuth (rad), measured clockwise from local north, in `[0, 2π)`.
    pub azimuth: f64,
    /// Slant range — straight-line station-to-satellite distance (m).
    pub range: f64,
}

/// One **access window**: a single contiguous pass during which the satellite
/// stays at or above the station's elevation mask.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AccessWindow {
    /// Rise time — when the satellite crosses *up* through the mask (s since
    /// epoch). Refined by linear interpolation between the bracketing samples.
    pub rise_time: f64,
    /// Set time — when it crosses *back down* through the mask (s).
    pub set_time: f64,
    /// Time of peak (maximum) elevation within the window (s).
    pub peak_time: f64,
    /// The peak elevation reached (rad).
    pub peak_elevation: f64,
}

impl AccessWindow {
    /// Window duration `set − rise` (s).
    pub fn duration(&self) -> f64 {
        self.set_time - self.rise_time
    }
}

/// Compute the satellite's [`LookAngles`] from a station, given the satellite's
/// **ECI** position and the Greenwich sidereal angle `gmst_rad` at that instant.
///
/// Builds the station→satellite vector in ECEF, rotates it into the local SEZ
/// (South-East-Zenith) topocentric frame, and reads off elevation, azimuth and
/// slant range. The SEZ basis at geodetic latitude `φ`, longitude `λ` is the
/// standard Vallado rotation `Ry(90° − φ)·Rz(λ)` applied to the ECEF
/// difference vector.
pub fn look_angles(station: &GroundStation, sat_eci: Vector3<f64>, gmst_rad: f64) -> LookAngles {
    let sat_ecef = eci_to_ecef(sat_eci, gmst_rad);
    let rho_ecef = sat_ecef - station.ecef(); // station -> satellite, ECEF

    let (sin_lat, cos_lat) = station.latitude.sin_cos();
    let (sin_lon, cos_lon) = station.longitude.sin_cos();

    // ECEF -> SEZ (South, East, Zenith).
    let s = sin_lat * cos_lon * rho_ecef.x + sin_lat * sin_lon * rho_ecef.y - cos_lat * rho_ecef.z;
    let e = -sin_lon * rho_ecef.x + cos_lon * rho_ecef.y;
    let z = cos_lat * cos_lon * rho_ecef.x + cos_lat * sin_lon * rho_ecef.y + sin_lat * rho_ecef.z;

    let range = (s * s + e * e + z * z).sqrt();
    if range < 1e-9 {
        // Satellite coincident with the station: no defined look direction.
        return LookAngles {
            elevation: std::f64::consts::FRAC_PI_2,
            azimuth: 0.0,
            range: 0.0,
        };
    }
    let elevation = (z / range).clamp(-1.0, 1.0).asin();
    // Azimuth clockwise from north. In SEZ, north is the −S direction.
    let mut azimuth = e.atan2(-s);
    if azimuth < 0.0 {
        azimuth += std::f64::consts::TAU;
    }
    LookAngles {
        elevation,
        azimuth,
        range,
    }
}

/// The elevation (rad) of `sat_eci` seen from `station` at sidereal angle
/// `gmst_rad` — a thin wrapper over [`look_angles`] for the common case where
/// only the elevation (visibility test) is needed.
pub fn elevation_of(station: &GroundStation, sat_eci: Vector3<f64>, gmst_rad: f64) -> f64 {
    look_angles(station, sat_eci, gmst_rad).elevation
}

/// Find all **access windows** of a satellite [`EphemerisPoint`] series from a
/// ground `station`, given the Greenwich sidereal angle `gmst0` at the
/// ephemeris epoch (e.g. from [`crate::frames::gmst`] at the start time).
///
/// Each sample's elevation is evaluated (the Earth-fixed frame is advanced from
/// `gmst0` at the sidereal rate to each sample time). A maximal run of samples
/// with elevation `≥ station.min_elevation` is one window; its rise/set edges
/// are refined to the mask-crossing by linear interpolation in elevation, and
/// the peak is sharpened by a 3-point parabolic fit. A pass already above the
/// mask at the first sample (or still above it at the last) is reported with
/// that endpoint as the rise/set time.
///
/// # Errors
///
/// Returns [`AstroError::InvalidParameter`] if any sample time or the station's
/// geometry is non-finite. An empty or single-point ephemeris yields an empty
/// window list (no pass can be bracketed), never an error or panic.
pub fn access_windows(
    ephemeris: &[EphemerisPoint],
    station: &GroundStation,
    gmst0: f64,
) -> Result<Vec<AccessWindow>, AstroError> {
    if !gmst0.is_finite()
        || !station.latitude.is_finite()
        || !station.longitude.is_finite()
        || !station.altitude.is_finite()
        || !station.min_elevation.is_finite()
    {
        return Err(AstroError::InvalidParameter(
            "access_windows: non-finite station geometry or gmst0",
        ));
    }

    // Precompute (time, elevation) for every sample.
    let mut elev = Vec::with_capacity(ephemeris.len());
    for p in ephemeris {
        if !p.time.is_finite() {
            return Err(AstroError::InvalidParameter(
                "access_windows: non-finite ephemeris sample time",
            ));
        }
        let theta = gmst_after(gmst0, p.time);
        elev.push((p.time, elevation_of(station, p.state.position, theta)));
    }

    let mask = station.min_elevation;
    let mut windows = Vec::new();
    if elev.len() < 2 {
        return Ok(windows); // cannot bracket a pass from <2 samples
    }

    // Linear-interpolate the time where elevation crosses the mask between two
    // bracketing samples.
    let cross_time = |i: usize, j: usize| -> f64 {
        let (t0, e0) = elev[i];
        let (t1, e1) = elev[j];
        let de = e1 - e0;
        if de.abs() < 1e-15 {
            t0
        } else {
            t0 + (mask - e0) / de * (t1 - t0)
        }
    };

    let mut in_pass = elev[0].1 >= mask;
    let mut rise_time = if in_pass { elev[0].0 } else { f64::NAN };

    for k in 1..elev.len() {
        let above = elev[k].1 >= mask;
        if above && !in_pass {
            // Rising edge between k-1 and k.
            rise_time = cross_time(k - 1, k);
            in_pass = true;
        } else if !above && in_pass {
            // Setting edge between k-1 and k.
            let set_time = cross_time(k - 1, k);
            windows.push(finish_window(&elev, rise_time, set_time));
            in_pass = false;
        }
    }
    // A pass still open at the end of the ephemeris closes at the last sample.
    if in_pass {
        let set_time = elev[elev.len() - 1].0;
        windows.push(finish_window(&elev, rise_time, set_time));
    }

    Ok(windows)
}

/// Assemble an [`AccessWindow`] from its refined rise/set times and the sampled
/// elevation series, locating the peak by a 3-point parabolic fit around the
/// highest in-window sample.
fn finish_window(elev: &[(f64, f64)], rise_time: f64, set_time: f64) -> AccessWindow {
    // Find the index of the maximum elevation sample within [rise, set].
    let mut best_i = 0usize;
    let mut best_e = f64::NEG_INFINITY;
    let mut found = false;
    for (idx, &(t, e)) in elev.iter().enumerate() {
        if t >= rise_time && t <= set_time && e > best_e {
            best_e = e;
            best_i = idx;
            found = true;
        }
    }
    if !found {
        // Degenerate (sub-sample) window: peak at the midpoint.
        let mid = 0.5 * (rise_time + set_time);
        return AccessWindow {
            rise_time,
            set_time,
            peak_time: mid,
            peak_elevation: best_e.max(0.0),
        };
    }

    // Parabolic refinement of the peak using neighbours when available.
    let (mut peak_time, mut peak_elev) = (elev[best_i].0, elev[best_i].1);
    if best_i > 0 && best_i + 1 < elev.len() {
        let (t0, e0) = elev[best_i - 1];
        let (t1, e1) = elev[best_i];
        let (t2, e2) = elev[best_i + 1];
        // Fit e(t) = parabola through the three points; vertex if concave.
        let denom = (t0 - t1) * (t0 - t2) * (t1 - t2);
        if denom.abs() > 1e-12 {
            let a = (t2 * (e1 - e0) + t1 * (e0 - e2) + t0 * (e2 - e1)) / denom;
            let b = (t2 * t2 * (e0 - e1) + t1 * t1 * (e2 - e0) + t0 * t0 * (e1 - e2)) / denom;
            if a < 0.0 {
                let tv = -b / (2.0 * a);
                if tv >= t0 && tv <= t2 {
                    let c = e1 - a * t1 * t1 - b * t1;
                    peak_time = tv;
                    peak_elev = a * tv * tv + b * tv + c;
                }
            }
        }
    }

    AccessWindow {
        rise_time,
        set_time,
        peak_time,
        peak_elevation: peak_elev,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{MU_EARTH, R_EARTH};
    use crate::frames::{gmst, J2000};
    use crate::orbit3d::{self, ClassicalElements};

    /// Build an ephemeris by propagating a COE over `total` seconds in `n`
    /// steps (two-body).
    fn ephemeris_from_coe(coe: &ClassicalElements, total: f64, n: usize) -> Vec<EphemerisPoint> {
        let mut state = orbit3d::coe_to_rv(coe).unwrap();
        let dt = total / n as f64;
        let mut out = Vec::with_capacity(n + 1);
        out.push(EphemerisPoint { time: 0.0, state });
        for i in 1..=n {
            state = orbit3d::propagate(&state, dt, 1, false).unwrap();
            out.push(EphemerisPoint {
                time: i as f64 * dt,
                state,
            });
        }
        out
    }

    #[test]
    fn satellite_directly_overhead_gives_near_90_deg_elevation() {
        // A station on the equator/prime meridian; place a satellite on the
        // local zenith (straight up along +x ECEF) at GMST = 0 (so ECI == ECEF
        // on the x-axis). Elevation must be ~90°.
        let station = GroundStation::from_degrees(0.0, 0.0, 0.0, 0.0);
        let alt = 800_000.0;
        let sat = Vector3::new(R_EARTH + alt, 0.0, 0.0); // straight up
        let la = look_angles(&station, sat, 0.0);
        assert!(
            (la.elevation.to_degrees() - 90.0).abs() < 1e-6,
            "elevation {} deg",
            la.elevation.to_degrees()
        );
        // Slant range equals the altitude when straight overhead.
        assert!((la.range - alt).abs() < 1.0, "range {} m", la.range);
    }

    #[test]
    fn satellite_below_the_horizon_has_negative_elevation_and_no_access() {
        // A satellite on the OPPOSITE side of the Earth from the station is far
        // below the local horizon → negative elevation, no access window.
        let station = GroundStation::from_degrees(0.0, 0.0, 0.0, 0.0);
        let alt = 800_000.0;
        let sat_behind = Vector3::new(-(R_EARTH + alt), 0.0, 0.0); // antipodal direction
        let elev = elevation_of(&station, sat_behind, 0.0);
        assert!(
            elev < 0.0,
            "elevation {} deg should be below horizon",
            elev.to_degrees()
        );

        // A flat ephemeris pinned on the far side yields zero windows.
        let eph: Vec<_> = (0..=20)
            .map(|i| EphemerisPoint {
                time: i as f64 * 60.0,
                state: StateVector {
                    position: sat_behind,
                    velocity: Vector3::zeros(),
                },
            })
            .collect();
        let win = access_windows(&eph, &station, 0.0).unwrap();
        assert!(win.is_empty(), "expected no access, got {win:?}");
    }

    #[test]
    fn circular_leo_pass_has_a_finite_sensible_duration() {
        // A 500 km circular inclined orbit. To make the pass deterministic
        // (independent of where the ground track happens to fall at this epoch)
        // we place the station directly beneath a chosen ephemeris sample's
        // sub-satellite point — guaranteeing a near-overhead pass — then confirm
        // there is at least one access window of a physically sensible length (a
        // few minutes, well under the orbital period) peaking near the zenith.
        let a = R_EARTH + 500_000.0;
        let coe = ClassicalElements {
            semi_major_axis: a,
            eccentricity: 0.0,
            inclination: 60.0_f64.to_radians(),
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let period = std::f64::consts::TAU * (a * a * a / MU_EARTH).sqrt();
        let eph = ephemeris_from_coe(&coe, 1.5 * period, 1500); // ~6 s sampling
        let gmst0 = gmst(J2000);

        // Pick a mid-pass sample and put the station under its sub-point. The
        // sub-satellite latitude/longitude is the geocentric sub-point at that
        // sample's sidereal angle (same convention as groundtrack::subpoint and
        // our eci_to_ecef rotation).
        let mid = eph.len() / 4; // a quarter orbit in — well inside the span
        let theta_mid = gmst_after(gmst0, eph[mid].time);
        let (sub_lat, sub_lon) = crate::groundtrack::subpoint(eph[mid].state.position, theta_mid);
        let station =
            GroundStation::from_degrees(sub_lat.to_degrees(), sub_lon.to_degrees(), 0.0, 5.0);

        let windows = access_windows(&eph, &station, gmst0).unwrap();
        assert!(!windows.is_empty(), "expected at least one LEO pass");

        // At least one window should peak high (near overhead, since the station
        // sits under the ground track).
        let max_peak = windows
            .iter()
            .map(|w| w.peak_elevation)
            .fold(0.0_f64, f64::max);
        assert!(
            max_peak.to_degrees() > 60.0,
            "best peak elevation only {}° for an overhead pass",
            max_peak.to_degrees()
        );

        for w in &windows {
            assert!(w.duration() > 0.0, "non-positive duration {w:?}");
            // A LEO pass is minutes long, never approaching the orbital period.
            assert!(
                w.duration() < 0.5 * period,
                "pass duration {} s implausibly long (period {} s)",
                w.duration(),
                period
            );
            // A near-overhead LEO pass is on the order of ~10 minutes; sanity
            // that the window is at least a minute (not a one-sample blip).
            assert!(
                w.duration() > 60.0 || w.peak_elevation < 10.0_f64.to_radians(),
                "high pass {} s too short",
                w.duration()
            );
            // Peak elevation is within the pass and above the mask.
            assert!(w.peak_elevation >= station.min_elevation - 1e-6);
            assert!(w.peak_time >= w.rise_time - 1e-6 && w.peak_time <= w.set_time + 1e-6);
            // Peak elevation cannot exceed 90°.
            assert!(w.peak_elevation <= std::f64::consts::FRAC_PI_2 + 1e-9);
        }
    }

    #[test]
    fn higher_mask_never_lengthens_a_pass() {
        // Raising the elevation mask can only shorten (or drop) a pass, never
        // extend it — a monotonicity sanity check on the window edges.
        let a = R_EARTH + 600_000.0;
        let coe = ClassicalElements {
            semi_major_axis: a,
            eccentricity: 0.0,
            inclination: 70.0_f64.to_radians(),
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let period = std::f64::consts::TAU * (a * a * a / MU_EARTH).sqrt();
        let eph = ephemeris_from_coe(&coe, 2.0 * period, 2000);
        let gmst0 = gmst(J2000);

        let low = GroundStation::from_degrees(50.0, 0.0, 0.0, 0.0);
        let high = GroundStation::from_degrees(50.0, 0.0, 0.0, 20.0);
        let total_low: f64 = access_windows(&eph, &low, gmst0)
            .unwrap()
            .iter()
            .map(AccessWindow::duration)
            .sum();
        let total_high: f64 = access_windows(&eph, &high, gmst0)
            .unwrap()
            .iter()
            .map(AccessWindow::duration)
            .sum();
        assert!(
            total_high <= total_low + 1e-6,
            "20° mask total {total_high}s exceeded 0° mask total {total_low}s"
        );
    }

    #[test]
    fn empty_and_single_point_ephemeris_are_graceful() {
        let station = GroundStation::from_degrees(0.0, 0.0, 0.0, 0.0);
        assert!(access_windows(&[], &station, 0.0).unwrap().is_empty());
        let one = [EphemerisPoint {
            time: 0.0,
            state: StateVector {
                position: Vector3::new(R_EARTH + 500_000.0, 0.0, 0.0),
                velocity: Vector3::zeros(),
            },
        }];
        assert!(access_windows(&one, &station, 0.0).unwrap().is_empty());
    }

    #[test]
    fn non_finite_inputs_error_not_panic() {
        let station = GroundStation::from_degrees(0.0, 0.0, 0.0, 0.0);
        let eph = [
            EphemerisPoint {
                time: 0.0,
                state: StateVector {
                    position: Vector3::new(R_EARTH + 500_000.0, 0.0, 0.0),
                    velocity: Vector3::zeros(),
                },
            },
            EphemerisPoint {
                time: f64::NAN,
                state: StateVector {
                    position: Vector3::new(R_EARTH + 500_000.0, 0.0, 0.0),
                    velocity: Vector3::zeros(),
                },
            },
        ];
        assert!(access_windows(&eph, &station, 0.0).is_err());
        assert!(access_windows(&[], &station, f64::INFINITY).is_err());
    }

    #[test]
    fn azimuth_points_the_right_way() {
        // From the equator/prime-meridian station, a satellite displaced to the
        // local east (toward +y at the station) should read azimuth ≈ 90°, and
        // one displaced north (toward +z) azimuth ≈ 0°.
        let station = GroundStation::from_degrees(0.0, 0.0, 0.0, 0.0);
        let base = R_EARTH + 700_000.0;
        // East: same x, +y offset, on the zenith line tilted east.
        let east = Vector3::new(base, 1.0e6, 0.0);
        let la_e = look_angles(&station, east, 0.0);
        assert!(
            (la_e.azimuth.to_degrees() - 90.0).abs() < 5.0,
            "az east {}",
            la_e.azimuth.to_degrees()
        );
        // North: +z offset (toward the pole).
        let north = Vector3::new(base, 0.0, 1.0e6);
        let la_n = look_angles(&station, north, 0.0);
        assert!(
            la_n.azimuth.to_degrees() < 5.0 || la_n.azimuth.to_degrees() > 355.0,
            "az north {}",
            la_n.azimuth.to_degrees()
        );
    }
}
