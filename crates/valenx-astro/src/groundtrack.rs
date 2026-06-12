//! Sub-satellite ground tracks: where a spacecraft is over the rotating
//! Earth.
//!
//! Converts an inertial (ECI) position to the **sub-satellite point** —
//! the geocentric latitude and longitude directly beneath the
//! spacecraft — by rotating into the Earth-fixed (ECEF) frame through
//! the current Greenwich rotation angle, then sampling that over an
//! orbit to produce a ground track.
//!
//! Scope: spherical Earth (geocentric latitude, not geodetic), and the
//! Earth-fixed frame is a simple rotation about the pole at `ω⊕`. Good
//! enough for coverage / visibility planning; not a geodetic-grade
//! ground-station product.

use serde::{Deserialize, Serialize};

use crate::constants::OMEGA_EARTH;
use crate::error::AstroError;
use crate::orbit3d::{self, ClassicalElements, StateVector};

/// Absolute ceiling on the number of ground-track samples a single call
/// may request, bounding the output `Vec`'s memory. A million points is
/// already a second-by-second track of an ~11-day repeat ground-track.
pub const MAX_GROUNDTRACK_SAMPLES: usize = 1_000_000;

/// A point on the ground track.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GroundPoint {
    /// Seconds since the start of the track.
    pub time: f64,
    /// Geocentric latitude (rad), in `[-π/2, π/2]`.
    pub latitude: f64,
    /// Longitude (rad), in `(-π, π]`, east-positive.
    pub longitude: f64,
}

/// Sub-satellite latitude / longitude (rad) for an ECI `position`, given
/// the Greenwich rotation angle `theta` (rad) of the Earth-fixed frame.
///
/// A zero or non-finite `position` has no defined sub-satellite point —
/// `position.z / ‖position‖` would be `0/0` (NaN) and silently poison the
/// latitude. In that degenerate case this returns `(0.0, 0.0)` rather than
/// a NaN coordinate.
pub fn subpoint(position: nalgebra::Vector3<f64>, theta: f64) -> (f64, f64) {
    let r = position.norm();
    if !r.is_finite() || r < 1e-12 {
        return (0.0, 0.0);
    }
    let (s, c) = theta.sin_cos();
    // Rotate inertial -> earth-fixed by -theta about z.
    let x_ecef = c * position.x + s * position.y;
    let y_ecef = -s * position.x + c * position.y;
    let latitude = (position.z / r).clamp(-1.0, 1.0).asin();
    let longitude = y_ecef.atan2(x_ecef);
    (latitude, longitude)
}

/// Sample the ground track of an orbit over `total_time` seconds in
/// `samples` steps, starting from Greenwich angle `theta0` (rad).
///
/// The orbit is propagated as a two-body (Keplerian) arc; the Earth
/// turns underneath at `ω⊕`, so successive passes shift west.
///
/// # Errors
///
/// Returns [`AstroError::InvalidIntegration`] if `total_time` is not
/// finite, [`AstroError::OutOfRange`] if `samples` exceeds
/// [`MAX_GROUNDTRACK_SAMPLES`], and [`AstroError::NonPhysicalState`] if
/// `coe` is non-physical (e.g. `e ≥ 1`, so the orbit cannot be sampled
/// as a closed ground track).
pub fn ground_track(
    coe: &ClassicalElements,
    theta0: f64,
    total_time: f64,
    samples: usize,
) -> Result<Vec<GroundPoint>, AstroError> {
    if !total_time.is_finite() {
        return Err(AstroError::InvalidIntegration("total_time not finite"));
    }
    if samples > MAX_GROUNDTRACK_SAMPLES {
        return Err(AstroError::OutOfRange {
            what: "samples",
            value: samples as u64,
            max: MAX_GROUNDTRACK_SAMPLES as u64,
        });
    }

    let mut state: StateVector = orbit3d::coe_to_rv(coe)?;
    let mut out = Vec::with_capacity(samples + 1);
    let dt_sample = if samples > 0 {
        total_time / samples as f64
    } else {
        total_time
    };

    for i in 0..=samples {
        let t = i as f64 * dt_sample;
        let theta = theta0 + OMEGA_EARTH * t;
        let (latitude, longitude) = subpoint(state.position, theta);
        out.push(GroundPoint {
            time: t,
            latitude,
            longitude,
        });
        // Advance the orbit by one sample interval (1 s RK4 sub-steps).
        if i < samples {
            let steps = dt_sample.max(1.0).round() as u64;
            let step_dt = dt_sample / steps as f64;
            state = orbit3d::propagate(&state, step_dt, steps, false)?;
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::R_EARTH;
    use nalgebra::Vector3;

    fn circular(inclination_deg: f64) -> ClassicalElements {
        ClassicalElements {
            semi_major_axis: R_EARTH + 500_000.0,
            eccentricity: 0.0,
            inclination: inclination_deg.to_radians(),
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        }
    }

    #[test]
    fn degenerate_position_is_no_op_not_nan() {
        // A zero (or non-finite) ECI position has no defined sub-point:
        // `position.z / norm()` is 0/0 -> NaN latitude pre-fix. The guard
        // returns (0, 0) ("undefined sub-point, hold") instead of poisoning
        // the track with NaN.
        for bad in [
            Vector3::zeros(),
            Vector3::new(f64::NAN, 0.0, 0.0),
            Vector3::new(f64::INFINITY, 0.0, 1.0),
        ] {
            let (lat, lon) = subpoint(bad, 0.7);
            assert_eq!((lat, lon), (0.0, 0.0), "degenerate position {bad:?}");
        }
    }

    #[test]
    fn subpoint_basics() {
        // On the +x axis with zero Greenwich angle -> lat 0, lon 0.
        let (lat, lon) = subpoint(Vector3::new(R_EARTH + 500_000.0, 0.0, 0.0), 0.0);
        assert!(lat.abs() < 1e-12 && lon.abs() < 1e-12);
        // Straight over the pole -> latitude +90°.
        let (latp, _) = subpoint(Vector3::new(0.0, 0.0, R_EARTH + 500_000.0), 0.0);
        assert!((latp.to_degrees() - 90.0).abs() < 1e-9);
    }

    #[test]
    fn earth_rotation_drifts_longitude_west() {
        // A fixed inertial point seen at a later Greenwich angle has a
        // smaller (more westward) longitude.
        let pos = Vector3::new(R_EARTH + 500_000.0, 0.0, 0.0);
        let (_, lon0) = subpoint(pos, 0.0);
        let (_, lon1) = subpoint(pos, 0.1);
        assert!(lon1 < lon0, "lon {lon0} -> {lon1}");
    }

    #[test]
    fn equatorial_orbit_stays_on_equator() {
        let track = ground_track(&circular(0.0), 0.0, 5_000.0, 50).expect("valid track");
        for p in &track {
            assert!(p.latitude.abs().to_degrees() < 0.1, "lat {}", p.latitude);
        }
    }

    #[test]
    fn rejects_absurd_sample_count_without_oom() {
        // The H2 repro: samples = usize::MAX would allocate / loop
        // unbounded. It must return an Err immediately.
        let err = ground_track(&circular(0.0), 0.0, 5_000.0, usize::MAX);
        assert!(
            matches!(
                err,
                Err(AstroError::OutOfRange {
                    what: "samples",
                    ..
                })
            ),
            "usize::MAX samples must be rejected, got {err:?}"
        );
        // Non-finite total_time is also rejected.
        assert!(ground_track(&circular(0.0), 0.0, f64::NAN, 10).is_err());
        // The cap itself is accepted in principle (use a tiny count here
        // to keep the test fast; the boundary is exercised by the cap
        // value, not by allocating a million points).
        assert!(ground_track(&circular(0.0), 0.0, 5_000.0, 10).is_ok());
    }

    #[test]
    fn inclined_orbit_reaches_its_inclination_in_latitude() {
        // Max sub-satellite latitude over an orbit equals the inclination.
        let coe = circular(51.6);
        let period = coe.period().unwrap();
        let track = ground_track(&coe, 0.0, period, 360).expect("valid track");
        let max_lat = track
            .iter()
            .map(|p| p.latitude.abs())
            .fold(0.0_f64, f64::max);
        assert!(
            (max_lat.to_degrees() - 51.6).abs() < 0.5,
            "max lat {}",
            max_lat.to_degrees()
        );
    }
}
