//! Solar geometry: Sun direction, orbit-plane beta angle, and eclipse fraction.
//!
//! These tie an orbit to the Sun — the missing piece for power and thermal
//! budgeting. The chain is:
//!
//! 1. [`sun_direction_eci`] — the geocentric unit vector to the Sun in the
//!    Earth-centred inertial (equatorial) frame, from the low-precision
//!    *Astronomical Almanac* series (good to ~0.01°).
//! 2. [`orbit_normal_eci`] — the unit normal of an orbit plane from its
//!    inclination and RAAN.
//! 3. [`beta_angle`] — the angle of the Sun above the orbit plane.
//! 4. [`eclipse_fraction`] — the fraction of a circular orbit spent in Earth's
//!    shadow, from the standard cylindrical-shadow model (Vallado).
//!
//! [`solar_geometry`] runs the whole chain for a [`ClassicalElements`] orbit at
//! a Julian date.
//!
//! # Honest scope
//!
//! Research/educational. The Sun series is the low-precision analytic model
//! (no planetary perturbations, no nutation) and the shadow model is a
//! cylinder (no penumbra, circular-orbit assumption — eccentric orbits use the
//! semi-major-axis altitude). Good for mission-geometry budgeting and teaching,
//! not for operational eclipse timing.

use nalgebra::Vector3;

use crate::constants::R_EARTH;
use crate::orbit3d::ClassicalElements;

/// Julian Date of the J2000.0 epoch (2000-01-01 12:00 TT).
pub const J2000: f64 = 2_451_545.0;

/// The geocentric unit vector pointing at the Sun in the ECI (mean-equator,
/// equatorial) frame, for Julian Date `julian_date`.
///
/// Uses the low-precision *Astronomical Almanac* solar series (mean longitude +
/// equation-of-centre, projected through the mean obliquity). The returned
/// vector is exactly unit length; its `z` component is `sin(declination)`.
pub fn sun_direction_eci(julian_date: f64) -> Vector3<f64> {
    let n = julian_date - J2000; // days since J2000.0
    let mean_longitude = (280.460 + 0.985_647_4 * n).to_radians();
    let mean_anomaly = (357.528 + 0.985_600_3 * n).to_radians();
    // Ecliptic longitude (equation of centre); ecliptic latitude is ~0.
    let ecliptic_longitude = mean_longitude
        + (1.915_f64).to_radians() * mean_anomaly.sin()
        + (0.020_f64).to_radians() * (2.0 * mean_anomaly).sin();
    let obliquity = (23.439 - 0.000_000_4 * n).to_radians();

    let (sin_lambda, cos_lambda) = ecliptic_longitude.sin_cos();
    let (sin_eps, cos_eps) = obliquity.sin_cos();
    // Rotate the ecliptic-plane unit vector (cosλ, sinλ, 0) into the equatorial
    // frame about the x-axis by the obliquity. The result is already unit.
    Vector3::new(cos_lambda, cos_eps * sin_lambda, sin_eps * sin_lambda)
}

/// The unit normal of an orbit plane with `inclination` and `raan` (both in
/// radians), in the ECI frame: `(sin i sin Ω, −sin i cos Ω, cos i)`.
pub fn orbit_normal_eci(inclination: f64, raan: f64) -> Vector3<f64> {
    let (sin_i, cos_i) = inclination.sin_cos();
    let (sin_o, cos_o) = raan.sin_cos();
    Vector3::new(sin_i * sin_o, -sin_i * cos_o, cos_i)
}

/// The beta angle (radians, in `[−π/2, π/2]`): the angle of the Sun above the
/// orbit plane. `±π/2` means the Sun is along the orbit normal (the orbit is
/// edge-on to the Sun, maximally sunlit); `0` means the Sun lies in the orbit
/// plane. Inputs need not be pre-normalised.
pub fn beta_angle(orbit_normal: &Vector3<f64>, sun_direction: &Vector3<f64>) -> f64 {
    let n = orbit_normal.normalize();
    let s = sun_direction.normalize();
    n.dot(&s).clamp(-1.0, 1.0).asin()
}

/// The fraction of a circular orbit at `altitude` (metres above Earth's mean
/// radius) spent in Earth's shadow, given the `beta` angle (radians), using the
/// cylindrical-shadow model.
///
/// Returns `0.0` when `|beta|` is at or beyond the critical angle
/// `arcsin(R⊕ / (R⊕ + altitude))` (the orbit is fully sunlit). At `beta = 0`
/// the fraction equals that critical angle over `π`.
pub fn eclipse_fraction(beta: f64, altitude: f64) -> f64 {
    let r = R_EARTH + altitude;
    if r <= R_EARTH {
        // At or below the surface there is no orbit; treat as always shadowed
        // on the night side is undefined — return 0 rather than NaN.
        return 0.0;
    }
    let beta_star = (R_EARTH / r).asin();
    if beta.abs() >= beta_star {
        return 0.0;
    }
    let arg = ((r * r - R_EARTH * R_EARTH).sqrt() / (r * beta.cos())).clamp(-1.0, 1.0);
    arg.acos() / std::f64::consts::PI
}

/// The solar geometry of an orbit at a given Julian date: beta angle, eclipse
/// fraction, and sunlit fraction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolarGeometry {
    /// Beta angle (radians).
    pub beta: f64,
    /// Fraction of the orbit in Earth's shadow, in `[0, 1]`.
    pub eclipse_fraction: f64,
    /// Fraction of the orbit in sunlight, in `[0, 1]` (`1 − eclipse_fraction`).
    pub sunlit_fraction: f64,
}

/// Compute the [`SolarGeometry`] of `elements` at `julian_date`.
///
/// The eclipse model is circular, so the altitude used is
/// `semi_major_axis − R⊕`; for an eccentric orbit this is the mean-altitude
/// approximation.
pub fn solar_geometry(elements: &ClassicalElements, julian_date: f64) -> SolarGeometry {
    let sun = sun_direction_eci(julian_date);
    let normal = orbit_normal_eci(elements.inclination, elements.raan);
    let beta = beta_angle(&normal, &sun);
    let altitude = elements.semi_major_axis - R_EARTH;
    let eclipse = eclipse_fraction(beta, altitude);
    SolarGeometry {
        beta,
        eclipse_fraction: eclipse,
        sunlit_fraction: 1.0 - eclipse,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn sun_direction_is_unit_length() {
        for &jd in &[J2000, J2000 + 100.0, J2000 + 200.0, J2000 + 365.25] {
            let v = sun_direction_eci(jd);
            assert!((v.norm() - 1.0).abs() < 1e-12, "‖sun‖ = {}", v.norm());
        }
    }

    #[test]
    fn summer_solstice_declination_is_the_obliquity() {
        // 2000-06-21 12:00 ≈ JD 2451717.0; the Sun reaches its max declination.
        let dec = sun_direction_eci(2_451_717.0).z.asin().to_degrees();
        assert!(
            (dec - 23.44).abs() < 0.2,
            "solstice declination {dec:.3}° should be ~+23.44°"
        );
    }

    #[test]
    fn vernal_equinox_declination_is_near_zero() {
        // 2000-03-20 12:00 ≈ JD 2451624.0; the Sun crosses the equator.
        let dec = sun_direction_eci(2_451_624.0).z.asin().to_degrees();
        assert!(
            dec.abs() < 1.0,
            "equinox declination {dec:.3}° should be ~0°"
        );
    }

    #[test]
    fn declination_envelope_matches_obliquity() {
        // Over a year the Sun's declination peaks at the obliquity (~23.44°).
        let mut max_abs_dec = 0.0_f64;
        for day in 0..=366 {
            let dec = sun_direction_eci(J2000 + f64::from(day))
                .z
                .asin()
                .to_degrees();
            max_abs_dec = max_abs_dec.max(dec.abs());
        }
        assert!(
            (max_abs_dec - 23.44).abs() < 0.1,
            "declination envelope {max_abs_dec:.3}° should be the obliquity ~23.44°"
        );
    }

    #[test]
    fn orbit_normal_special_cases() {
        // Equatorial prograde orbit: normal is +Z.
        let eq = orbit_normal_eci(0.0, 0.0);
        assert!((eq - Vector3::new(0.0, 0.0, 1.0)).norm() < 1e-12);
        // Polar orbit, RAAN 0: normal is −Y.
        let polar = orbit_normal_eci(FRAC_PI_2, 0.0);
        assert!((polar - Vector3::new(0.0, -1.0, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn beta_is_90_when_sun_along_normal_and_0_when_in_plane() {
        let normal = Vector3::new(0.0, 0.0, 1.0);
        assert!((beta_angle(&normal, &Vector3::new(0.0, 0.0, 1.0)) - FRAC_PI_2).abs() < 1e-12);
        assert!(beta_angle(&normal, &Vector3::new(1.0, 0.0, 0.0)).abs() < 1e-12);
    }

    #[test]
    fn fully_sunlit_when_beta_exceeds_critical() {
        // h = 400 km -> critical beta ≈ 70.2°; at 80° there is no eclipse.
        let altitude = 400_000.0;
        assert_eq!(eclipse_fraction(80.0_f64.to_radians(), altitude), 0.0);
    }

    #[test]
    fn beta_zero_fraction_equals_critical_angle_over_pi() {
        // The closed-form check: at beta = 0 the shadow fraction is exactly
        // arcsin(R/r)/π.
        let altitude = 400_000.0;
        let r = R_EARTH + altitude;
        let expected = (R_EARTH / r).asin() / PI;
        let got = eclipse_fraction(0.0, altitude);
        assert!(
            (got - expected).abs() < 1e-9,
            "f(β=0) = {got:.6} should equal arcsin(R/r)/π = {expected:.6}"
        );
        // Sanity: a 400 km orbit spends ~39% of the period in shadow at β = 0.
        assert!(
            (0.36..0.42).contains(&got),
            "LEO β=0 eclipse fraction {got:.3}"
        );
    }

    #[test]
    fn solar_geometry_partitions_the_orbit() {
        let elements = ClassicalElements {
            semi_major_axis: R_EARTH + 700_000.0,
            eccentricity: 0.0,
            inclination: 98.0_f64.to_radians(), // sun-sync-ish
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let g = solar_geometry(&elements, J2000);
        assert!((g.eclipse_fraction + g.sunlit_fraction - 1.0).abs() < 1e-12);
        assert!((0.0..=1.0).contains(&g.eclipse_fraction));
        assert!(g.beta.abs() <= FRAC_PI_2 + 1e-12);
    }
}
