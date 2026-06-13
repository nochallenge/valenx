//! Geocentric equatorial coordinates — right ascension and declination.
//!
//! Where [`crate::groundtrack`] projects a spacecraft onto the *rotating* Earth
//! (latitude / longitude), this module gives its direction on the **inertial
//! celestial sphere**: the geocentric **right ascension** `α` (the equatorial
//! analogue of longitude, measured eastward from the vernal equinox / +x axis)
//! and **declination** `δ` (the analogue of latitude, above/below the equator),
//! plus the slant **range**. These are the angles an Earth-centred observer
//! tracks a satellite or star by, and the rates at which they sweep.
//!
//! For an ECI position `(x, y, z)` with `r = ‖(x,y,z)‖`:
//!
//! ```text
//! α = atan2(y, x)   (wrapped to [0, 2π))
//! δ = asin(z / r)
//! ```

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Geocentric equatorial coordinates of a point: right ascension, declination
/// and range.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Equatorial {
    /// Right ascension `α` (rad), in `[0, 2π)`, measured eastward in the
    /// equatorial plane from the +x (vernal-equinox) axis.
    pub right_ascension: f64,
    /// Declination `δ` (rad), in `[−π/2, π/2]`, the angle above (+) or below
    /// (−) the equatorial plane.
    pub declination: f64,
    /// Range `r` (m) — the geocentric distance to the point.
    pub range: f64,
}

/// The **geocentric equatorial coordinates** [`Equatorial`] of an ECI position.
///
/// `α = atan2(y, x)` wrapped to `[0, 2π)`, `δ = asin(z/r)`, `range = ‖r‖`. A
/// point on the +x axis has `α = 0, δ = 0`; the north pole has `δ = +π/2`. A
/// zero or non-finite position has no defined direction and returns all-zero
/// (rather than a `NaN` coordinate).
pub fn geocentric_equatorial(position: Vector3<f64>) -> Equatorial {
    let r = position.norm();
    if !r.is_finite() || r < 1e-12 {
        return Equatorial {
            right_ascension: 0.0,
            declination: 0.0,
            range: 0.0,
        };
    }
    let mut ra = position.y.atan2(position.x);
    if ra < 0.0 {
        ra += std::f64::consts::TAU; // wrap to [0, 2π)
    }
    let declination = (position.z / r).clamp(-1.0, 1.0).asin();
    Equatorial {
        right_ascension: ra,
        declination,
        range: r,
    }
}

/// The **angular rates** `(dα/dt, dδ/dt)` (rad/s) of the geocentric right
/// ascension and declination, from an ECI position and velocity.
///
/// Differentiating the [`geocentric_equatorial`] angles:
///
/// ```text
/// dα/dt = (x·vy − y·vx) / (x² + y²)
/// dδ/dt = (vz − ṙ·z/r) / ρ      with ṙ = (r·v)/r,  ρ = √(x² + y²)
/// ```
///
/// For a circular **equatorial** orbit `dα/dt` equals the mean motion
/// `n = √(μ/r³)` and `dδ/dt = 0`; at an ascending node an inclined orbit has
/// `dδ/dt > 0` (climbing north). Returns `(0, 0)` at or near the pole
/// (`ρ ≈ 0`, where right ascension is undefined) and for a degenerate position,
/// rather than a `NaN` rate.
pub fn geocentric_angular_rates(position: Vector3<f64>, velocity: Vector3<f64>) -> (f64, f64) {
    let (x, y, z) = (position.x, position.y, position.z);
    let (vx, vy, vz) = (velocity.x, velocity.y, velocity.z);
    let rho2 = x * x + y * y;
    let rho = rho2.sqrt();
    let r = position.norm();
    if !r.is_finite() || r < 1e-12 || rho < 1e-12 {
        return (0.0, 0.0);
    }
    let ra_dot = (x * vy - y * vx) / rho2;
    let r_dot = (x * vx + y * vy + z * vz) / r;
    let dec_dot = (vz - r_dot * z / r) / rho;
    (ra_dot, dec_dot)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{MU_EARTH, R_EARTH};
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, SQRT_2, TAU};

    #[test]
    fn equatorial_axis_points() {
        let r = 7.0e6;
        // +x → RA 0, Dec 0, range r.
        let e = geocentric_equatorial(Vector3::new(r, 0.0, 0.0));
        assert!(e.right_ascension.abs() < 1e-12 && e.declination.abs() < 1e-12);
        assert!((e.range - r).abs() < 1e-9);
        // +y → RA 90°.
        let e = geocentric_equatorial(Vector3::new(0.0, r, 0.0));
        assert!((e.right_ascension - FRAC_PI_2).abs() < 1e-12);
        // −x → RA 180°.
        assert!(
            (geocentric_equatorial(Vector3::new(-r, 0.0, 0.0)).right_ascension - PI).abs() < 1e-12
        );
        // North pole → Dec +90°.
        assert!(
            (geocentric_equatorial(Vector3::new(0.0, 0.0, r)).declination - FRAC_PI_2).abs()
                < 1e-12
        );
        // −y → RA wraps to 270° (in [0, 2π), not −90°).
        assert!(
            (geocentric_equatorial(Vector3::new(0.0, -r, 0.0)).right_ascension - 3.0 * FRAC_PI_2)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn equatorial_45_degree_worked_point() {
        // (1, 1, √2): ρ = √2, r = 2 → δ = asin(√2/2) = 45°, α = atan2(1,1) = 45°.
        let e = geocentric_equatorial(Vector3::new(1.0, 1.0, SQRT_2));
        assert!(
            (e.right_ascension - FRAC_PI_4).abs() < 1e-12,
            "RA {}",
            e.right_ascension
        );
        assert!(
            (e.declination - FRAC_PI_4).abs() < 1e-12,
            "Dec {}",
            e.declination
        );
        assert!((e.range - 2.0).abs() < 1e-12);
    }

    #[test]
    fn circular_equatorial_orbit_ra_advances_at_mean_motion() {
        // GROUND TRUTH: for a circular EQUATORIAL orbit the geocentric right
        // ascension advances at exactly the mean motion n = √(μ/r³), with the
        // declination held at 0.
        let r = R_EARTH + 500_000.0;
        let v = (MU_EARTH / r).sqrt(); // circular speed
        let n = (MU_EARTH / (r * r * r)).sqrt(); // mean motion
        let (ra_dot, dec_dot) =
            geocentric_angular_rates(Vector3::new(r, 0.0, 0.0), Vector3::new(0.0, v, 0.0));
        assert!((ra_dot - n).abs() / n < 1e-12, "RA-dot {ra_dot} != n {n}");
        assert!(dec_dot.abs() < 1e-12, "equatorial Dec-dot {dec_dot} != 0");
    }

    #[test]
    fn inclined_orbit_climbs_in_declination_at_ascending_node() {
        // At the ascending node (on +x, crossing the equator northbound) a
        // prograde inclined circular orbit has a +z velocity → declination
        // rising at exactly v·sin(i)/r (ṙ = 0 there).
        let r = R_EARTH + 500_000.0;
        let v = (MU_EARTH / r).sqrt();
        let inc = 45.0_f64.to_radians();
        let (ra_dot, dec_dot) = geocentric_angular_rates(
            Vector3::new(r, 0.0, 0.0),
            Vector3::new(0.0, v * inc.cos(), v * inc.sin()),
        );
        assert!(ra_dot > 0.0, "prograde RA increases: {ra_dot}");
        let expected = v * inc.sin() / r;
        assert!(
            (dec_dot - expected).abs() / expected < 1e-12,
            "Dec-dot {dec_dot} != v·sin(i)/r {expected}"
        );
    }

    #[test]
    fn equatorial_guards_degenerate_input() {
        // Zero position → all-zero (no defined direction), no NaN.
        let e = geocentric_equatorial(Vector3::zeros());
        assert_eq!((e.right_ascension, e.declination, e.range), (0.0, 0.0, 0.0));
        assert_eq!(
            geocentric_equatorial(Vector3::new(f64::NAN, 0.0, 0.0)).range,
            0.0
        );
        // At the pole the right ascension (and its rate) are undefined → 0.
        let (ra_dot, dec_dot) =
            geocentric_angular_rates(Vector3::new(0.0, 0.0, 7.0e6), Vector3::new(1.0, 0.0, 0.0));
        assert_eq!((ra_dot, dec_dot), (0.0, 0.0));
        // RA always lands in [0, 2π).
        for &(x, y) in &[(1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0), (1.0, -1.0)] {
            let ra = geocentric_equatorial(Vector3::new(x, y, 0.3)).right_ascension;
            assert!((0.0..TAU).contains(&ra), "RA {ra} out of [0, 2π)");
        }
    }
}
