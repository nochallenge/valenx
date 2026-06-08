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

/// The **specific angular momentum of a circular orbit** `h = √(μ·r)` (m²/s) at radius
/// `radius` `r` (m) about Earth — the conserved angular momentum per unit mass, equal to
/// the circular speed times the radius (`h = v_circ·r`). It is the constant in Kepler's
/// second law: a circular orbit sweeps area at the steady rate `dA/dt = h/2`, so over one
/// [`orbital_period`] it covers the full disc `πr²`, giving `h·T = 2π·r²`.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `radius` is non-finite or non-positive.
pub fn circular_angular_momentum(radius: f64) -> Result<f64, AstroError> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "circular_angular_momentum radius must be finite and > 0",
        ));
    }
    Ok((MU_EARTH * radius).sqrt())
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

/// The **specific orbital energy from a state** `ε = ½·v² − μ/r` (J/kg) — the
/// instantaneous total energy per unit mass (kinetic plus gravitational potential) of a
/// body moving at speed `speed` `v` (m/s) at radius `radius` `r` (m) about Earth. This is
/// the vis-viva energy that fixes the orbit from a state vector: it is conserved along the
/// trajectory and equal to [`specific_orbital_energy`]`(a)` `= −μ/(2a)`. Its sign
/// classifies the orbit — `ε < 0` bound (elliptical), `ε = 0` parabolic (exactly escape
/// speed), `ε > 0` unbound (hyperbolic).
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `speed` is non-finite, or `radius` is
/// non-finite or non-positive.
pub fn specific_orbital_energy_from_state(speed: f64, radius: f64) -> Result<f64, AstroError> {
    if !speed.is_finite() || !radius.is_finite() || radius <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "specific_orbital_energy_from_state requires finite speed and radius > 0",
        ));
    }
    Ok(0.5 * speed * speed - MU_EARTH / radius)
}

/// The **semi-major axis from the specific orbital energy** `a = −μ / (2·ε)` (m) — the
/// inverse of [`specific_orbital_energy`]. Composed with
/// [`specific_orbital_energy_from_state`] it is the vis-viva *orbit determination* step:
/// a state vector's speed and radius give the energy `ε`, which fixes the orbit's size
/// `a`. Only a bound (elliptical) orbit has a finite semi-major axis — its energy is
/// negative — so `ε ≥ 0` (a parabolic or hyperbolic trajectory) is rejected.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `specific_energy` is non-finite or
/// non-negative (not a bound orbit).
pub fn semi_major_axis_from_energy(specific_energy: f64) -> Result<f64, AstroError> {
    if !specific_energy.is_finite() || specific_energy >= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "semi_major_axis_from_energy requires a bound orbit (specific energy < 0)",
        ));
    }
    Ok(-MU_EARTH / (2.0 * specific_energy))
}

/// The **hyperbolic excess speed** `v∞ = √(2·ε)` (m/s) — the residual speed an unbound
/// orbit keeps infinitely far from Earth, where the gravitational well has been fully
/// climbed and only the surplus specific energy `ε` remains as kinetic energy. Its square
/// is the characteristic energy `C3 = v∞²` that sizes an interplanetary departure. It is
/// the escape-regime complement of [`semi_major_axis_from_energy`]: composed with
/// [`specific_orbital_energy_from_state`], a state vector's energy maps to the orbit's
/// semi-major axis when bound (`ε < 0`) and to this excess speed when unbound (`ε > 0`).
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `specific_energy` is non-finite or
/// non-positive (a bound or parabolic orbit never reaches infinity with speed to spare).
pub fn hyperbolic_excess_speed(specific_energy: f64) -> Result<f64, AstroError> {
    if !specific_energy.is_finite() || specific_energy <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "hyperbolic_excess_speed requires an unbound orbit (specific energy > 0)",
        ));
    }
    Ok((2.0 * specific_energy).sqrt())
}

/// The **semi-major axis from the orbital period** `a = (μ·T² / (4π²))^(1/3)` (m) —
/// Kepler's third law inverted: given an observed period `period` `T` (s) it recovers
/// the size of the orbit about Earth. It is the inverse of [`orbital_period`] and the
/// relation that fixes mission altitudes — a sidereal-day period gives the
/// geostationary radius (≈ 42 164 km), the GPS 11 h 58 m period its ≈ 26 560 km orbit.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `period` is non-finite or non-positive.
pub fn semi_major_axis_from_period(period: f64) -> Result<f64, AstroError> {
    if !period.is_finite() || period <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "semi_major_axis_from_period period must be finite and > 0",
        ));
    }
    Ok((MU_EARTH * period * period / (4.0 * std::f64::consts::PI * std::f64::consts::PI)).cbrt())
}

/// The **circular-orbit radius from the circular speed** `r = μ / v²` (m) — the orbit
/// radius about Earth at which a body circles at speed `speed` `v` (m/s); the inverse of
/// [`circular_speed`]. Given a measured circular orbital speed it recovers the altitude,
/// and a faster orbit is a lower (smaller) one.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `speed` is non-finite or non-positive.
pub fn orbital_radius_from_circular_speed(speed: f64) -> Result<f64, AstroError> {
    if !speed.is_finite() || speed <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "orbital_radius_from_circular_speed speed must be finite and > 0",
        ));
    }
    Ok(MU_EARTH / (speed * speed))
}

/// The **orbit radius from the escape speed** `r = 2μ / v²` (m) — the radius about Earth
/// at which the local escape speed is `speed` `v` (m/s); the inverse of [`escape_speed`].
/// Since the escape speed is `√2` times the circular speed at the same radius, for a
/// given speed this radius is exactly twice [`orbital_radius_from_circular_speed`].
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `speed` is non-finite or non-positive.
pub fn orbital_radius_from_escape_speed(speed: f64) -> Result<f64, AstroError> {
    if !speed.is_finite() || speed <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "orbital_radius_from_escape_speed speed must be finite and > 0",
        ));
    }
    Ok(2.0 * MU_EARTH / (speed * speed))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orbital_radius_from_escape_speed_inverts_escape_speed() {
        // Round-trip: recover r from the escape speed it implies (the exact inverse of
        // √(2μ/r)).
        for &r in &[R_EARTH, R_EARTH + 400_000.0, 4.2164e7] {
            let recovered = orbital_radius_from_escape_speed(escape_speed(r).unwrap()).unwrap();
            assert!((recovered - r).abs() <= 1e-12 * r, "r = 2μ/v² inverts v = √(2μ/r)");
        }

        // Worked: Earth's surface escape speed ≈ 11.18 km/s → back to ≈ Earth's radius
        // (a rough textbook anchor at 1%; the exact relation is the round-trip).
        let r = orbital_radius_from_escape_speed(11_180.0).unwrap();
        assert!((r - R_EARTH).abs() / R_EARTH < 1e-2, "11.18 km/s → ≈ Earth's surface");

        // Threads orbital_radius_from_circular_speed: for the same speed the escape radius
        // is exactly twice the circular radius (escape speed is √2 × circular speed).
        for &v in &[3000.0, 7672.6, 11_180.0] {
            assert!(
                (orbital_radius_from_escape_speed(v).unwrap()
                    - 2.0 * orbital_radius_from_circular_speed(v).unwrap())
                .abs()
                    <= 1e-9 * orbital_radius_from_escape_speed(v).unwrap(),
                "escape radius = 2 × circular radius at the same speed"
            );
        }

        // Monotonic decreasing in speed; Err on non-physical input.
        assert!(
            orbital_radius_from_escape_speed(12000.0).unwrap()
                < orbital_radius_from_escape_speed(8000.0).unwrap()
        );
        assert!(orbital_radius_from_escape_speed(0.0).is_err());
        assert!(orbital_radius_from_escape_speed(-1.0).is_err());
        assert!(orbital_radius_from_escape_speed(f64::NAN).is_err());
    }

    #[test]
    fn orbital_radius_from_circular_speed_inverts_circular_speed() {
        // Round-trip: recover r from the circular speed it implies (the exact inverse
        // of √(μ/r)).
        for &r in &[R_EARTH, R_EARTH + 400_000.0, 4.2164e7] {
            let recovered = orbital_radius_from_circular_speed(circular_speed(r).unwrap()).unwrap();
            assert!((recovered - r).abs() <= 1e-12 * r, "r = μ/v² inverts v = √(μ/r)");
        }

        // Worked: a ≈ 7.67 km/s LEO orbital speed gives roughly a 400 km altitude (a
        // rough textbook anchor within 1%; the exact relationship is the round-trip).
        let r = orbital_radius_from_circular_speed(7672.6).unwrap();
        assert!(
            (r - (R_EARTH + 400_000.0)).abs() / (R_EARTH + 400_000.0) < 1e-2,
            "7.67 km/s → ≈ 400 km LEO, got r = {r}"
        );

        // A faster circular speed means a smaller orbit.
        assert!(
            orbital_radius_from_circular_speed(8000.0).unwrap()
                < orbital_radius_from_circular_speed(3000.0).unwrap(),
            "faster → smaller orbit"
        );

        // Non-physical speed → error.
        assert!(orbital_radius_from_circular_speed(0.0).is_err());
        assert!(orbital_radius_from_circular_speed(-1.0).is_err());
        assert!(orbital_radius_from_circular_speed(f64::NAN).is_err());
    }

    #[test]
    fn semi_major_axis_from_period_inverts_kepler_third_law() {
        // Round-trip: recover a from the period it implies (the exact inverse of
        // orbital_period).
        for &a in &[R_EARTH + 400_000.0, 1.5e7, 4.2164e7] {
            let recovered = semi_major_axis_from_period(orbital_period(a).unwrap()).unwrap();
            assert!((recovered - a).abs() <= 1e-12 * a, "a = (μT²/4π²)^⅓ inverts T = 2π√(a³/μ)");
        }

        // Worked: a sidereal-day period gives the geostationary radius ≈ 42 164 km.
        let geo = semi_major_axis_from_period(86164.0).unwrap();
        assert!((geo - 4.2164e7).abs() / 4.2164e7 < 1e-3, "GEO radius ≈ 42164 km, got {geo}");

        // Monotonic increasing in the period.
        assert!(
            semi_major_axis_from_period(7200.0).unwrap()
                < semi_major_axis_from_period(86164.0).unwrap(),
            "longer period → larger orbit"
        );

        // Non-physical period → error.
        assert!(semi_major_axis_from_period(0.0).is_err());
        assert!(semi_major_axis_from_period(-1.0).is_err());
        assert!(semi_major_axis_from_period(f64::NAN).is_err());
    }

    #[test]
    fn hyperbolic_excess_speed_is_the_residual_speed_of_an_unbound_orbit() {
        let r = R_EARTH + 400_000.0;

        // Threads specific_orbital_energy_from_state via the vis-viva limit
        // v∞ = √(v² − 2μ/r): a clearly-hyperbolic launch (v = 1.2·escape).
        let v = 1.2 * escape_speed(r).unwrap();
        let eps = specific_orbital_energy_from_state(v, r).unwrap();
        assert!(eps > 0.0, "1.2·escape is unbound");
        let vinf = hyperbolic_excess_speed(eps).unwrap();
        assert!(
            (vinf - (v * v - 2.0 * MU_EARTH / r).sqrt()).abs() <= 1e-9 * vinf,
            "v∞ = √(v² − 2μ/r)"
        );

        // Round-trips the energy: v∞² = 2ε.
        assert!((vinf.powi(2) - 2.0 * eps).abs() <= 1e-9 * (2.0 * eps), "v∞² = 2ε");

        // Worked: ε = 2e7 → v∞ = √(4e7) ≈ 6324.6 m/s.
        assert!(
            (hyperbolic_excess_speed(2.0e7).unwrap() - (4.0e7_f64).sqrt()).abs()
                <= 1e-9 * (4.0e7_f64).sqrt(),
            "√(2·2e7) = √4e7"
        );

        // Monotonic increasing in ε.
        assert!(hyperbolic_excess_speed(1.0e8).unwrap() > hyperbolic_excess_speed(2.0e7).unwrap());

        // Err for a bound/parabolic orbit (ε ≤ 0) — the complement of
        // semi_major_axis_from_energy.
        assert!(hyperbolic_excess_speed(-1.0e7).is_err());
        assert!(hyperbolic_excess_speed(0.0).is_err());
        assert!(hyperbolic_excess_speed(f64::NAN).is_err());
    }

    #[test]
    fn semi_major_axis_from_energy_inverts_the_orbit_energy() {
        // Round-trip inverting specific_orbital_energy: a → ε → a.
        for &a in &[7_000_000.0_f64, 4.2164e7, R_EARTH + 800_000.0] {
            let recovered = semi_major_axis_from_energy(specific_orbital_energy(a).unwrap()).unwrap();
            assert!((recovered - a).abs() <= 1e-9 * a, "a = −μ/2ε inverts ε = −μ/2a");
        }

        // Vis-viva orbit determination: recover a from a STATE (speed + radius), threading
        // specific_orbital_energy_from_state + orbital_speed.
        for &(r, a) in &[(7_000_000.0_f64, 7_000_000.0_f64), (6_800_000.0, 8_000_000.0)] {
            let energy = specific_orbital_energy_from_state(orbital_speed(r, a).unwrap(), r).unwrap();
            let recovered = semi_major_axis_from_energy(energy).unwrap();
            assert!((recovered - a).abs() <= 1e-9 * a, "state → ε → a");
        }

        // Worked: ε = −μ/(2·7e6) → a = 7e6.
        let a = semi_major_axis_from_energy(-MU_EARTH / (2.0 * 7_000_000.0)).unwrap();
        assert!((a - 7_000_000.0).abs() <= 1e-9 * 7_000_000.0, "a = −μ/2ε");

        // Err on a non-bound orbit (parabolic ε = 0, hyperbolic ε > 0) or NaN.
        assert!(semi_major_axis_from_energy(0.0).is_err());
        assert!(semi_major_axis_from_energy(1.0e6).is_err());
        assert!(semi_major_axis_from_energy(f64::NAN).is_err());
    }

    #[test]
    fn specific_orbital_energy_from_state_is_the_vis_viva_energy() {
        // Energy conservation: ½v² − μ/r at the vis-viva speed equals −μ/2a (the orbit's
        // own energy), threading orbital_speed + specific_orbital_energy.
        for &(r, a) in &[
            (7_000_000.0_f64, 7_000_000.0_f64),
            (6_800_000.0, 8_000_000.0),
            (4.2164e7, 4.2164e7),
        ] {
            let from_state =
                specific_orbital_energy_from_state(orbital_speed(r, a).unwrap(), r).unwrap();
            let from_a = specific_orbital_energy(a).unwrap();
            assert!(
                (from_state - from_a).abs() <= 1e-9 * from_a.abs(),
                "½v² − μ/r = −μ/2a (vis-viva invariant)"
            );
        }

        // At escape speed the total energy is exactly zero (parabolic).
        let e_esc =
            specific_orbital_energy_from_state(escape_speed(R_EARTH).unwrap(), R_EARTH).unwrap();
        assert!(e_esc.abs() < 1.0, "escape speed → ε = 0, got {e_esc}");

        // Sign classifies the orbit: bound at circular speed, unbound above escape.
        assert!(
            specific_orbital_energy_from_state(circular_speed(R_EARTH).unwrap(), R_EARTH).unwrap()
                < 0.0,
            "circular orbit is bound (ε < 0)"
        );
        assert!(
            specific_orbital_energy_from_state(escape_speed(R_EARTH).unwrap() * 1.5, R_EARTH)
                .unwrap()
                > 0.0,
            "hyperbolic is unbound (ε > 0)"
        );

        // Err on non-physical radius or speed.
        assert!(specific_orbital_energy_from_state(7000.0, 0.0).is_err());
        assert!(specific_orbital_energy_from_state(7000.0, -1.0).is_err());
        assert!(specific_orbital_energy_from_state(7000.0, f64::NAN).is_err());
        assert!(specific_orbital_energy_from_state(f64::NAN, 7.0e6).is_err());
    }

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
    fn circular_angular_momentum_threads_speed_and_keplers_area_law() {
        for &r in &[R_EARTH, R_EARTH + 400_000.0, 4.2164e7] {
            // Threads circular_speed: h = v_circ·r.
            let h = circular_angular_momentum(r).unwrap();
            assert!((h - circular_speed(r).unwrap() * r).abs() <= 1e-9 * h, "h = v_circ·r");

            // Threads orbital_period via Kepler's 2nd law: h·T = 2π·r² (area πr² swept per
            // period at dA/dt = h/2).
            let area_law = circular_angular_momentum(r).unwrap() * orbital_period(r).unwrap();
            assert!(
                (area_law - 2.0 * std::f64::consts::PI * r * r).abs()
                    <= 1e-9 * (2.0 * std::f64::consts::PI * r * r),
                "h·T = 2π·r²"
            );
        }

        // Monotonic increasing in radius.
        assert!(
            circular_angular_momentum(4.2164e7).unwrap()
                > circular_angular_momentum(R_EARTH).unwrap()
        );

        // Err on non-physical radius.
        assert!(circular_angular_momentum(0.0).is_err());
        assert!(circular_angular_momentum(-1.0).is_err());
        assert!(circular_angular_momentum(f64::NAN).is_err());
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
