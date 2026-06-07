//! Keplerian orbital elements from a planar inertial state vector.
//!
//! Given a position and velocity in the Earth-centred inertial frame
//! (and `μ`), recover the conic the vehicle is on: semi-major axis,
//! eccentricity, apoapsis / periapsis radii and altitudes, specific
//! orbital energy and (for a bound orbit) the period.

use nalgebra::Vector2;
use serde::{Deserialize, Serialize};

use crate::constants::{MU_EARTH, R_EARTH};
use crate::error::AstroError;

/// The conic section a state vector is on.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct OrbitElements {
    /// Semi-major axis (m). Negative for a hyperbolic trajectory.
    pub semi_major_axis: f64,
    /// Eccentricity (dimensionless).
    pub eccentricity: f64,
    /// Apoapsis radius from Earth's centre (m); infinite if unbound.
    pub apoapsis_radius: f64,
    /// Periapsis radius from Earth's centre (m).
    pub periapsis_radius: f64,
    /// Apoapsis altitude above the equatorial radius (m); infinite if
    /// unbound.
    pub apoapsis_altitude: f64,
    /// Periapsis altitude above the equatorial radius (m). Negative
    /// means the orbit re-enters the atmosphere / impacts.
    pub periapsis_altitude: f64,
    /// Specific orbital energy `v²/2 − μ/r` (J/kg).
    pub specific_energy: f64,
    /// Orbital period (s) for a bound orbit, else `None`.
    pub period: Option<f64>,
    /// True when the orbit is gravitationally bound (energy < 0).
    pub is_bound: bool,
}

/// Compute the orbital elements for a planar state, using Earth's `μ`.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] for input that would
/// otherwise yield a silent `NaN`/`Inf` orbit: a zero or non-finite
/// position, a non-finite velocity, or the parabolic energy singularity
/// (specific energy ≈ 0, where the semi-major axis blows up).
pub fn elements(position: Vector2<f64>, velocity: Vector2<f64>) -> Result<OrbitElements, AstroError> {
    elements_with_mu(position, velocity, MU_EARTH)
}

/// Compute orbital elements for an arbitrary central body `μ`.
///
/// # Errors
///
/// As [`elements`], plus rejects a non-finite or non-positive `mu`.
pub fn elements_with_mu(
    position: Vector2<f64>,
    velocity: Vector2<f64>,
    mu: f64,
) -> Result<OrbitElements, AstroError> {
    if !mu.is_finite() || mu <= 0.0 {
        return Err(AstroError::NonPhysicalState("mu must be finite and > 0"));
    }
    let r = position.norm();
    if !r.is_finite() || r <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "position magnitude must be finite and > 0",
        ));
    }
    if !velocity.x.is_finite() || !velocity.y.is_finite() {
        return Err(AstroError::NonPhysicalState("velocity must be finite"));
    }
    let specific_energy = velocity.norm_squared() / 2.0 - mu / r;
    // The semi-major axis is a = -mu/(2·ε); ε ≈ 0 (parabolic) makes it
    // blow up to ±Inf. Reject that singular case rather than emit Inf.
    if specific_energy.abs() < f64::EPSILON {
        return Err(AstroError::NonPhysicalState(
            "parabolic energy singularity (specific energy ≈ 0)",
        ));
    }
    Ok(elements_with_mu_unchecked(position, velocity, mu))
}

/// Element-construction core without input validation. Internal use only,
/// for callers that pass a state already known to be finite, non-zero and
/// non-parabolic (e.g. the integrated ascent state). Mirrors the public
/// [`elements_with_mu`] math exactly.
pub(crate) fn elements_with_mu_unchecked(
    position: Vector2<f64>,
    velocity: Vector2<f64>,
    mu: f64,
) -> OrbitElements {
    let r = position.norm();
    let v2 = velocity.norm_squared();
    let specific_energy = v2 / 2.0 - mu / r;

    // Specific angular momentum (scalar in the plane): r × v.
    let h = position.x * velocity.y - position.y * velocity.x;

    // e² = 1 + 2 E h² / μ²  (guard tiny negatives from round-off).
    let e2 = 1.0 + 2.0 * specific_energy * h * h / (mu * mu);
    let eccentricity = e2.max(0.0).sqrt();

    let is_bound = specific_energy < 0.0;
    let semi_major_axis = -mu / (2.0 * specific_energy);

    let (apoapsis_radius, periapsis_radius, period) = if is_bound {
        let ra = semi_major_axis * (1.0 + eccentricity);
        let rp = semi_major_axis * (1.0 - eccentricity);
        let t = 2.0 * std::f64::consts::PI * (semi_major_axis.powi(3) / mu).sqrt();
        (ra, rp, Some(t))
    } else {
        // Parabolic / hyperbolic: periapsis is still finite.
        let rp = if eccentricity > 1.0 {
            semi_major_axis * (1.0 - eccentricity)
        } else {
            // Parabola: periapsis = h²/(2μ).
            h * h / (2.0 * mu)
        };
        (f64::INFINITY, rp, None)
    };

    OrbitElements {
        semi_major_axis,
        eccentricity,
        apoapsis_radius,
        periapsis_radius,
        apoapsis_altitude: apoapsis_radius - R_EARTH,
        periapsis_altitude: periapsis_radius - R_EARTH,
        specific_energy,
        period,
        is_bound,
    }
}

/// Circular orbital speed (m/s) at a given radius from Earth's centre.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `radius` is non-finite or
/// non-positive, which would otherwise yield a `NaN`/`Inf` speed.
pub fn circular_speed(radius: f64) -> Result<f64, AstroError> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "circular_speed radius must be finite and > 0",
        ));
    }
    Ok((MU_EARTH / radius).sqrt())
}

/// Escape speed (m/s) at a given radius from Earth's centre — the minimum speed
/// for an unbound (parabolic) trajectory, `v_esc = √(2·μ/r) = √2 · v_circ`. At this
/// speed the specific orbital energy is exactly zero, so the body just reaches
/// infinity with no kinetic energy to spare.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `radius` is non-finite or
/// non-positive, which would otherwise yield a `NaN`/`Inf` speed.
pub fn escape_speed(radius: f64) -> Result<f64, AstroError> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "escape_speed radius must be finite and > 0",
        ));
    }
    Ok((2.0 * MU_EARTH / radius).sqrt())
}

/// The **Keplerian orbital period** `T = 2π·√(a³/μ)` (s) of an orbit with semi-major
/// axis `semi_major_axis` `a` (m) about Earth — Kepler's third law, the time to
/// complete one revolution. For a circular orbit (`a = r`) it is equivalently
/// `2π·r / v_circ`; it grows with the 3/2 power of the semi-major axis, so a 400 km
/// LEO takes ≈ 92 min while a geostationary orbit takes a sidereal day.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `semi_major_axis` is non-finite or
/// non-positive.
pub fn orbital_period(semi_major_axis: f64) -> Result<f64, AstroError> {
    if !semi_major_axis.is_finite() || semi_major_axis <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "orbital_period semi_major_axis must be finite and > 0",
        ));
    }
    Ok(2.0 * std::f64::consts::PI * (semi_major_axis.powi(3) / MU_EARTH).sqrt())
}

/// The **orbital speed** `v = √(μ·(2/r − 1/a))` (m/s) at radius `radius` `r` (m) on an
/// orbit of semi-major axis `semi_major_axis` `a` (m) about Earth — the **vis-viva
/// equation**, expressing conservation of specific orbital energy `v²/2 − μ/r =
/// −μ/(2a)`. It is the general speed of which [`circular_speed`] (`r = a`) and
/// [`escape_speed`] (the `a → ∞` parabolic limit) are special cases: faster at
/// periapsis, slower at apoapsis, with the product `v·r` constant on the line of
/// apsides (the conserved angular momentum).
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `radius` or `semi_major_axis` is
/// non-finite or non-positive, or if `radius > 2a` (beyond the apoapsis the orbit
/// does not reach, where `v²` would be negative).
pub fn orbital_speed(radius: f64, semi_major_axis: f64) -> Result<f64, AstroError> {
    if !radius.is_finite()
        || !semi_major_axis.is_finite()
        || radius <= 0.0
        || semi_major_axis <= 0.0
    {
        return Err(AstroError::NonPhysicalState(
            "orbital_speed radius and semi_major_axis must be finite and > 0",
        ));
    }
    let v_sq = MU_EARTH * (2.0 / radius - 1.0 / semi_major_axis);
    if v_sq < 0.0 {
        return Err(AstroError::NonPhysicalState(
            "orbital_speed radius lies beyond the orbit (radius > 2·semi_major_axis)",
        ));
    }
    Ok(v_sq.sqrt())
}

/// The **specific orbital energy** `ε = −μ / (2·a)` (J/kg) of an orbit with semi-major
/// axis `semi_major_axis` `a` (m) about Earth — the conserved total mechanical energy
/// per unit mass (kinetic plus gravitational potential), the constant of the vis-viva
/// equation `v²/2 − μ/r = −μ/(2a)`. It is negative for a bound ellipse and rises toward
/// `0` as the orbit grows (`a → ∞`, the parabolic escape limit), so a higher (less
/// negative) energy means a larger orbit.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `semi_major_axis` is non-finite or
/// non-positive (a bound orbit requires `a > 0`).
pub fn specific_orbital_energy(semi_major_axis: f64) -> Result<f64, AstroError> {
    if !semi_major_axis.is_finite() || semi_major_axis <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "specific_orbital_energy semi_major_axis must be finite and > 0",
        ));
    }
    Ok(-MU_EARTH / (2.0 * semi_major_axis))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn specific_orbital_energy_matches_vis_viva() {
        // Threads orbital_speed via the vis-viva energy ε = v²/2 − μ/r, for any reachable
        // r on the orbit (all of a, 1.3a, 0.8a satisfy r < 2a).
        for &a in &[R_EARTH + 400_000.0, 1.0e7, 2.5e7] {
            let eps = specific_orbital_energy(a).unwrap();
            for &r in &[a, 1.3 * a, 0.8 * a] {
                let from_speed = orbital_speed(r, a).unwrap().powi(2) / 2.0 - MU_EARTH / r;
                assert!((eps - from_speed).abs() <= 1e-12 * eps.abs(), "ε = v²/2 − μ/r");
            }
        }

        // Worked: a 400 km circular LEO has ε ≈ −2.943e7 J/kg.
        let leo = specific_orbital_energy(R_EARTH + 400_000.0).unwrap();
        assert!((leo - (-2.943e7)).abs() / 2.943e7 < 1e-3, "LEO energy ≈ −29.4 MJ/kg, got {leo}");

        // Larger orbit → higher (less negative) energy; bound orbits have ε < 0.
        assert!(specific_orbital_energy(2.0e7).unwrap() > specific_orbital_energy(1.0e7).unwrap());
        assert!(specific_orbital_energy(1.0e7).unwrap() < 0.0, "bound orbit ε < 0");

        // Non-physical semi-major axis → error.
        assert!(specific_orbital_energy(0.0).is_err());
        assert!(specific_orbital_energy(-1.0).is_err());
        assert!(specific_orbital_energy(f64::NAN).is_err());
    }

    #[test]
    fn orbital_speed_matches_vis_viva_and_threads_circular_escape() {
        // Circular orbit (r = a): vis-viva reduces to the circular speed.
        for &r in &[R_EARTH, R_EARTH + 400_000.0, 4.2e7] {
            let v = orbital_speed(r, r).unwrap();
            assert!((v - circular_speed(r).unwrap()).abs() <= 1e-12 * v, "orbital_speed(a,a) = v_circ");
        }

        // Parabolic limit (a → ∞): vis-viva tends to the escape speed.
        let r = R_EARTH + 1.0e6;
        let v_far = orbital_speed(r, 1.0e15).unwrap();
        assert!(
            (v_far - escape_speed(r).unwrap()).abs() / escape_speed(r).unwrap() < 1e-6,
            "orbital_speed(r, ∞) → v_esc"
        );

        // Angular momentum: v·r is equal at periapsis and apoapsis of an ellipse.
        let a = 1.5 * R_EARTH;
        let e = 0.3;
        let v_peri = orbital_speed(a * (1.0 - e), a).unwrap();
        let v_apo = orbital_speed(a * (1.0 + e), a).unwrap();
        assert!(
            (v_peri * a * (1.0 - e) - v_apo * a * (1.0 + e)).abs() <= 1e-9 * (v_peri * a * (1.0 - e)),
            "v·r conserved across the apsides"
        );
        assert!(v_peri > v_apo, "faster at periapsis");

        // Worked: a 400 km circular LEO orbits at ≈ 7.67 km/s.
        let leo = orbital_speed(R_EARTH + 400_000.0, R_EARTH + 400_000.0).unwrap();
        assert!((leo - 7670.0).abs() < 20.0, "LEO speed ≈ 7.67 km/s, got {leo}");

        // Non-physical: r beyond apoapsis (r > 2a), non-positive, or NaN → Err.
        assert!(orbital_speed(3.0 * a, a).is_err(), "r > 2a → Err");
        assert!(orbital_speed(0.0, a).is_err());
        assert!(orbital_speed(r, -1.0).is_err());
        assert!(orbital_speed(f64::NAN, a).is_err());
    }

    #[test]
    fn orbital_period_matches_kepler_third_law() {
        use std::f64::consts::PI;
        // Threads circular_speed: for a circular orbit T = 2π·r / v_circ = 2π√(r³/μ).
        for &r in &[R_EARTH, R_EARTH + 400_000.0, 4.2e7] {
            let t = orbital_period(r).unwrap();
            let from_speed = 2.0 * PI * r / circular_speed(r).unwrap();
            assert!((t - from_speed).abs() / t < 1e-12, "T = 2πr/v_circ");
        }
        // Kepler's third law: T ∝ a^(3/2), so quadrupling a multiplies T by 8.
        let a = R_EARTH + 1.0e6;
        let ta = orbital_period(a).unwrap();
        let t4a = orbital_period(4.0 * a).unwrap();
        assert!((t4a - 8.0 * ta).abs() / ta < 1e-12, "T(4a) = 8·T(a)");
        // A 400 km LEO orbits in ≈ 92 min (≈ 5544 s).
        let leo = orbital_period(R_EARTH + 400_000.0).unwrap();
        assert!((leo - 5544.0).abs() < 60.0, "LEO period ≈ 92 min, got {leo} s");
        // Non-physical semi-major axis → error.
        assert!(orbital_period(0.0).is_err());
        assert!(orbital_period(-1.0).is_err());
        assert!(orbital_period(f64::NAN).is_err());
    }

    #[test]
    fn escape_speed_is_root_two_times_circular() {
        // v_esc = √2 · v_circ (threads circular_speed) and v_esc² = 2μ/r (zero energy).
        for &r in &[R_EARTH, R_EARTH + 400_000.0, 4.2e7] {
            let v_esc = escape_speed(r).unwrap();
            let v_circ = circular_speed(r).unwrap();
            assert!((v_esc - 2.0_f64.sqrt() * v_circ).abs() / v_esc < 1e-12, "v_esc = √2·v_circ");
            assert!(
                (v_esc * v_esc - 2.0 * MU_EARTH / r).abs() / (v_esc * v_esc) < 1e-12,
                "v_esc² = 2μ/r"
            );
        }
        // Earth-surface escape speed is the textbook ≈ 11.2 km/s.
        let surface = escape_speed(R_EARTH).unwrap();
        assert!((11_000.0..11_400.0).contains(&surface), "surface escape ≈ 11.2 km/s, got {surface}");
        // A higher orbit escapes more slowly.
        assert!(
            escape_speed(R_EARTH + 1.0e6).unwrap() < escape_speed(R_EARTH).unwrap(),
            "escape speed decreases with radius"
        );
        // Non-physical radius → error.
        assert!(escape_speed(0.0).is_err());
        assert!(escape_speed(-1.0).is_err());
        assert!(escape_speed(f64::NAN).is_err());
    }

    #[test]
    fn circular_orbit_has_zero_eccentricity() {
        let radius = R_EARTH + 400_000.0; // 400 km LEO
        let v = circular_speed(radius).expect("valid radius");
        let pos = Vector2::new(radius, 0.0);
        let vel = Vector2::new(0.0, v); // perpendicular -> circular
        let o = elements(pos, vel).expect("valid circular state");
        assert!(o.eccentricity < 1e-9, "e = {}", o.eccentricity);
        assert!(o.is_bound);
        assert!((o.apoapsis_altitude - 400_000.0).abs() < 1.0);
        assert!((o.periapsis_altitude - 400_000.0).abs() < 1.0);
        // Period of a 400 km orbit ≈ 5554 s.
        let t = o.period.unwrap();
        assert!((t - 5_554.0).abs() < 20.0, "period {t}");
    }

    #[test]
    fn elliptical_orbit_apo_peri_ordered() {
        let rp = R_EARTH + 200_000.0;
        // Speed a bit above circular at periapsis -> raises apoapsis.
        let v = circular_speed(rp).expect("valid radius") * 1.1;
        let pos = Vector2::new(rp, 0.0);
        let vel = Vector2::new(0.0, v);
        let o = elements(pos, vel).expect("valid elliptical state");
        assert!(o.is_bound);
        assert!(o.eccentricity > 0.0 && o.eccentricity < 1.0);
        assert!(o.apoapsis_radius > o.periapsis_radius);
        assert!((o.periapsis_radius - rp).abs() < 1.0);
    }

    #[test]
    fn escape_velocity_is_unbound() {
        let r = R_EARTH + 300_000.0;
        let v_esc = (2.0 * MU_EARTH / r).sqrt();
        let pos = Vector2::new(r, 0.0);
        let vel = Vector2::new(0.0, v_esc * 1.01);
        let o = elements(pos, vel).expect("valid hyperbolic state");
        assert!(!o.is_bound);
        assert!(o.apoapsis_radius.is_infinite());
        assert!(o.eccentricity > 1.0);
    }

    #[test]
    fn suborbital_has_negative_periapsis_altitude() {
        // Straight up at well below orbital speed -> periapsis is
        // deep inside the Earth (negative altitude): it comes back down.
        let r = R_EARTH + 100_000.0;
        let pos = Vector2::new(r, 0.0);
        let vel = Vector2::new(0.0, circular_speed(r).expect("valid radius") * 0.5);
        let o = elements(pos, vel).expect("valid suborbital state");
        assert!(o.periapsis_altitude < 0.0);
    }

    #[test]
    fn zero_position_is_rejected_not_nan() {
        // r = 0 -> mu/r is Inf -> specific_energy is -Inf and the whole
        // orbit used to come out NaN/Inf silently. Must be a clean Err.
        let r = elements(Vector2::zeros(), Vector2::new(0.0, 7_800.0));
        assert!(
            matches!(r, Err(AstroError::NonPhysicalState(_))),
            "zero position must be rejected, got {r:?}"
        );
    }

    #[test]
    fn non_finite_input_is_rejected() {
        let bad_pos = elements(Vector2::new(f64::NAN, 0.0), Vector2::new(0.0, 7_800.0));
        assert!(matches!(bad_pos, Err(AstroError::NonPhysicalState(_))));
        let r = R_EARTH + 400_000.0;
        let bad_vel = elements(Vector2::new(r, 0.0), Vector2::new(f64::INFINITY, 0.0));
        assert!(matches!(bad_vel, Err(AstroError::NonPhysicalState(_))));
    }

    #[test]
    fn parabolic_energy_singularity_is_rejected_not_inf() {
        // Exactly-escape speed -> specific energy == 0 -> a = -mu/0 = Inf
        // semi-major axis silently. Must be rejected instead.
        let r = R_EARTH + 300_000.0;
        let v_esc = (2.0 * MU_EARTH / r).sqrt();
        let pos = Vector2::new(r, 0.0);
        let vel = Vector2::new(0.0, v_esc); // energy exactly 0
        let res = elements(pos, vel);
        assert!(
            matches!(res, Err(AstroError::NonPhysicalState(_))),
            "parabolic singularity must be rejected, got {res:?}"
        );
    }

    #[test]
    fn circular_speed_rejects_non_positive_radius() {
        assert!(matches!(
            circular_speed(0.0),
            Err(AstroError::NonPhysicalState(_))
        ));
        assert!(matches!(
            circular_speed(-1.0),
            Err(AstroError::NonPhysicalState(_))
        ));
        assert!(matches!(
            circular_speed(f64::NAN),
            Err(AstroError::NonPhysicalState(_))
        ));
        // A normal radius still works.
        assert!(circular_speed(R_EARTH + 400_000.0).is_ok());
    }
}
