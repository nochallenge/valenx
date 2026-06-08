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

/// The **specific gravitational potential** `Φ = −μ / r` (J/kg) at radius `radius` `r` (m)
/// from Earth's centre — the gravitational potential energy per unit mass, the depth of the
/// gravity well. It is negative and climbs toward `0` as `r → ∞`. It is the potential half
/// of the energy decomposition: the [`specific_orbital_energy_from_state`] is just the
/// kinetic term plus this, `ε = ½v² + Φ`; and escape is climbing out of it, so the
/// [`escape_speed`] satisfies `v_esc = √(−2Φ)` (zero total energy at the rim).
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `radius` is non-finite or non-positive.
pub fn gravitational_potential(radius: f64) -> Result<f64, AstroError> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "gravitational_potential radius must be finite and > 0",
        ));
    }
    Ok(-MU_EARTH / radius)
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

/// The **radius from vis-viva** `r = 2 / (v²/μ + 1/a)` (m) — the unique radius at which a
/// body on a bound orbit of semi-major axis `semi_major_axis` `a` (m) moves at speed `speed`
/// `v` (m/s). It inverts [`orbital_speed`] (vis-viva `v = √(μ(2/r − 1/a))`, which is strictly
/// decreasing in `r`, so the inverse is single-valued), and is the general case the circular
/// ([`orbital_radius_from_circular_speed`]) and escape ([`orbital_radius_from_escape_speed`])
/// inverses specialise.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `speed` is non-finite or negative, or
/// `semi_major_axis` is non-finite or non-positive.
pub fn orbital_radius_from_speed(speed: f64, semi_major_axis: f64) -> Result<f64, AstroError> {
    if !speed.is_finite()
        || speed < 0.0
        || !semi_major_axis.is_finite()
        || semi_major_axis <= 0.0
    {
        return Err(AstroError::NonPhysicalState(
            "orbital_radius_from_speed requires finite speed >= 0 and semi_major_axis > 0",
        ));
    }
    let inv_r = speed * speed / (2.0 * MU_EARTH) + 1.0 / (2.0 * semi_major_axis);
    Ok(1.0 / inv_r)
}

/// The **orbit eccentricity from its apsis radii** `e = (r_a − r_p) / (r_a + r_p)`
/// (dimensionless) — the shape of the conic specified by its apoapsis radius
/// `apoapsis_radius` `r_a` and periapsis radius `periapsis_radius` `r_p` (m). Together with
/// the semi-major axis `a = (r_a + r_p)/2` it fully fixes the orbit (`r_a = a(1+e)`,
/// `r_p = a(1−e)`); it is `0` for a circle (`r_a = r_p`) and tends to `1` as the orbit
/// becomes radial (`r_p → 0`). This recovers the eccentricity that [`elements`] computes from
/// a state vector.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if either radius is non-finite, the periapsis is
/// non-positive, or the apoapsis is smaller than the periapsis.
pub fn eccentricity_from_apsides(
    apoapsis_radius: f64,
    periapsis_radius: f64,
) -> Result<f64, AstroError> {
    if !apoapsis_radius.is_finite()
        || !periapsis_radius.is_finite()
        || periapsis_radius <= 0.0
        || apoapsis_radius < periapsis_radius
    {
        return Err(AstroError::NonPhysicalState(
            "eccentricity_from_apsides requires finite radii with apoapsis >= periapsis > 0",
        ));
    }
    Ok((apoapsis_radius - periapsis_radius) / (apoapsis_radius + periapsis_radius))
}

/// The **semi-major axis of a conic from its apsis radii** `a = (r_a + r_p) / 2`
/// (m) — the orbit size specified by its apoapsis radius `apoapsis_radius` `r_a`
/// and periapsis radius `periapsis_radius` `r_p` (m), the arithmetic mean of the
/// two turning-point distances. Together with the shape
/// [`eccentricity_from_apsides`] (which takes the same `(r_a, r_p)`) it fully fixes
/// the orbit — `r_a = a(1+e)`, `r_p = a(1−e)` — so `a` then feeds the period
/// [`orbital_period`], the energy [`specific_orbital_energy`], and the vis-viva
/// speed [`orbital_speed`]. For a circle (`r_a = r_p`) it is just the common radius.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if either radius is non-finite, the periapsis is
/// non-positive, or the apoapsis is smaller than the periapsis.
pub fn semi_major_axis_from_apsides(
    apoapsis_radius: f64,
    periapsis_radius: f64,
) -> Result<f64, AstroError> {
    if !apoapsis_radius.is_finite()
        || !periapsis_radius.is_finite()
        || periapsis_radius <= 0.0
        || apoapsis_radius < periapsis_radius
    {
        return Err(AstroError::NonPhysicalState(
            "semi_major_axis_from_apsides requires finite radii with apoapsis >= periapsis > 0",
        ));
    }
    Ok((apoapsis_radius + periapsis_radius) / 2.0)
}

/// The **semi-latus rectum of a conic from its orbital elements** `p = a·(1 − e²)`
/// (m) — the conic parameter that sets the orbit equation `r = p / (1 + e·cosθ)` and
/// the specific angular momentum `h = √(μ·p)`. `semi_major_axis` `a` (m) is the orbit
/// size and `eccentricity` `e` (`0 ≤ e < 1` for a bound ellipse) its shape. It is the
/// orbital radius at the ends of the latus rectum (true anomaly `±90°`) and lies between
/// the apsides (`r_p ≤ p ≤ r_a`); for a circle (`e = 0`) it is just `a`, and it equals
/// the harmonic mean of the apsis radii, `2·r_a·r_p / (r_a + r_p)`.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `a` is non-finite or non-positive, or `e`
/// is non-finite or outside the bound-ellipse range `0 ≤ e < 1`.
pub fn semi_latus_rectum_from_elements(
    semi_major_axis: f64,
    eccentricity: f64,
) -> Result<f64, AstroError> {
    if !semi_major_axis.is_finite()
        || semi_major_axis <= 0.0
        || !eccentricity.is_finite()
        || !(0.0..1.0).contains(&eccentricity)
    {
        return Err(AstroError::NonPhysicalState(
            "semi_latus_rectum_from_elements requires semi_major_axis > 0 and 0 <= eccentricity < 1",
        ));
    }
    Ok(semi_major_axis * (1.0 - eccentricity * eccentricity))
}

/// The **specific orbital angular momentum from its orbital elements**
/// `h = √(μ·a·(1 − e²)) = √(μ·p)` (m²/s) — the conserved angular momentum per unit mass of
/// the two-body orbit (twice the areal velocity, Kepler's second law), where `μ` is the
/// Earth gravitational parameter, `semi_major_axis` `a` (m) is the orbit size, and
/// `eccentricity` `e` (`0 ≤ e < 1` for a bound ellipse) its shape. Equivalently
/// `h = √(μ·p)` with `p` the [`semi_latus_rectum_from_elements`]; at the apsides it is
/// `h = r·v` (position and velocity are perpendicular there), and for a circle (`e = 0`)
/// it reduces to the [`circular_angular_momentum`] `√(μ·a)`.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if `a` is non-finite or non-positive, or `e`
/// is non-finite or outside the bound-ellipse range `0 ≤ e < 1`.
pub fn specific_angular_momentum_from_elements(
    semi_major_axis: f64,
    eccentricity: f64,
) -> Result<f64, AstroError> {
    if !semi_major_axis.is_finite()
        || semi_major_axis <= 0.0
        || !eccentricity.is_finite()
        || !(0.0..1.0).contains(&eccentricity)
    {
        return Err(AstroError::NonPhysicalState(
            "specific_angular_momentum_from_elements requires semi_major_axis > 0 and 0 <= eccentricity < 1",
        ));
    }
    Ok((MU_EARTH * semi_major_axis * (1.0 - eccentricity * eccentricity)).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eccentricity_from_apsides_specifies_the_conic() {
        // Worked: r_a=8e6, r_p=7e6 → e = 1e6/1.5e7 ≈ 0.066667.
        let e = eccentricity_from_apsides(8.0e6, 7.0e6).unwrap();
        assert!((e - 1.0e6 / 1.5e7).abs() <= 1e-9 * e, "e = (r_a−r_p)/(r_a+r_p)");

        // Circular: r_a == r_p → e == 0.
        assert_eq!(eccentricity_from_apsides(7.0e6, 7.0e6).unwrap(), 0.0, "circular → e = 0");

        // Reconstruction: a and e recover BOTH apsides (r_a = a(1+e), r_p = a(1−e)).
        for &(r_a, r_p) in &[(8.0e6, 7.0e6), (4.2e7, 7.0e6), (1.0e7, 1.0e7)] {
            let a = (r_a + r_p) / 2.0;
            let e = eccentricity_from_apsides(r_a, r_p).unwrap();
            assert!((a * (1.0 + e) - r_a).abs() <= 1e-9 * r_a, "r_a = a(1+e)");
            assert!((a * (1.0 - e) - r_p).abs() <= 1e-9 * r_p, "r_p = a(1−e)");
            assert!((0.0..1.0).contains(&e), "0 ≤ e < 1 for a bound ellipse");
        }

        // Near-radial: r_p → 0 gives e → 1.
        assert!(eccentricity_from_apsides(1.0e7, 1.0).unwrap() > 0.999, "near-radial → e ≈ 1");

        // Err on non-physical input.
        assert!(eccentricity_from_apsides(7.0e6, 0.0).is_err());
        assert!(eccentricity_from_apsides(7.0e6, 8.0e6).is_err()); // apoapsis < periapsis
        assert!(eccentricity_from_apsides(f64::NAN, 7.0e6).is_err());
    }

    #[test]
    fn semi_major_axis_from_apsides_averages_the_apsis_radii() {
        // (a) WORKED: a = (8e6 + 7e6)/2 = 7.5e6.
        let a = semi_major_axis_from_apsides(8.0e6, 7.0e6).unwrap();
        assert!((a - 7.5e6).abs() <= 1e-9 * a, "a = (r_a+r_p)/2 = 7.5e6");

        // (b) CIRCULAR: r_a == r_p → a == that radius (non-zero float → tolerance).
        let ac = semi_major_axis_from_apsides(7.0e6, 7.0e6).unwrap();
        assert!((ac - 7.0e6).abs() <= 1e-9 * 7.0e6, "circular → a = r");

        // (c) ROUND-TRIP threading #372 eccentricity_from_apsides (non-tautological):
        // a and e recover BOTH input apsides (r_a = a(1+e), r_p = a(1−e)).
        for &(r_a, r_p) in &[(8.0e6_f64, 7.0e6_f64), (4.2e7, 7.0e6), (1.0e7, 1.0e7)] {
            let a = semi_major_axis_from_apsides(r_a, r_p).unwrap();
            let e = eccentricity_from_apsides(r_a, r_p).unwrap();
            assert!((a * (1.0 + e) - r_a).abs() <= 1e-9 * r_a, "r_a = a(1+e)");
            assert!((a * (1.0 - e) - r_p).abs() <= 1e-9 * r_p, "r_p = a(1−e)");
        }

        // (d) BOUND: periapsis ≤ a ≤ apoapsis, and (e) the orbit is faster at
        // periapsis than apoapsis for the same a (vis-viva, threads orbital_speed).
        let (r_a, r_p) = (4.2e7_f64, 7.0e6_f64);
        let a = semi_major_axis_from_apsides(r_a, r_p).unwrap();
        assert!(r_p <= a && a <= r_a, "r_p ≤ a ≤ r_a");
        assert!(
            orbital_speed(r_p, a).unwrap() > orbital_speed(r_a, a).unwrap(),
            "periapsis speed > apoapsis speed"
        );

        // (f) Err on non-physical input (mirrors eccentricity_from_apsides).
        assert!(semi_major_axis_from_apsides(7.0e6, 0.0).is_err());
        assert!(semi_major_axis_from_apsides(7.0e6, 8.0e6).is_err()); // apoapsis < periapsis
        assert!(semi_major_axis_from_apsides(f64::NAN, 7.0e6).is_err());
    }

    #[test]
    fn semi_latus_rectum_from_elements_is_the_conic_parameter() {
        // (a) WORKED: a = 7e6, e = 0.1 → p = 7e6·(1 − 0.01) = 6.93e6.
        let p = semi_latus_rectum_from_elements(7.0e6, 0.1).unwrap();
        assert!((p - 6.93e6).abs() <= 1e-9 * p, "p = a(1−e²) = 6.93e6");

        // (b) CIRCLE: e = 0 → p = a.
        let pc = semi_latus_rectum_from_elements(7.0e6, 0.0).unwrap();
        assert!((pc - 7.0e6).abs() <= 1e-9 * 7.0e6, "circle → p = a");

        // (c) HARMONIC-MEAN cross-check threading #372 + #378 (non-tautological): for
        // apsides (r_a, r_p), p = a(1−e²) equals 2·r_a·r_p/(r_a+r_p), and (d) it lies
        // between the apsides.
        for &(r_a, r_p) in &[(8.0e6_f64, 7.0e6_f64), (4.2e7, 7.0e6), (1.0e7, 9.0e6)] {
            let a = semi_major_axis_from_apsides(r_a, r_p).unwrap();
            let e = eccentricity_from_apsides(r_a, r_p).unwrap();
            let p = semi_latus_rectum_from_elements(a, e).unwrap();
            let harmonic = 2.0 * r_a * r_p / (r_a + r_p);
            assert!((p - harmonic).abs() <= 1e-9 * harmonic, "p = 2·r_a·r_p/(r_a+r_p)");
            assert!(r_p <= p && p <= r_a, "r_p ≤ p ≤ r_a");
        }

        // (e) Err on non-physical input.
        assert!(semi_latus_rectum_from_elements(-1.0, 0.1).is_err()); // a ≤ 0
        assert!(semi_latus_rectum_from_elements(7.0e6, 1.0).is_err()); // e ≥ 1
        assert!(semi_latus_rectum_from_elements(7.0e6, -0.1).is_err()); // e < 0
        assert!(semi_latus_rectum_from_elements(f64::NAN, 0.1).is_err());
    }

    #[test]
    fn specific_angular_momentum_from_elements_is_sqrt_mu_p() {
        // (a) WORKED: h = √(μ·a·(1−e²)). a = 7e6, e = 0.1.
        let h = specific_angular_momentum_from_elements(7.0e6, 0.1).unwrap();
        let expected = (MU_EARTH * 7.0e6 * (1.0 - 0.01)).sqrt();
        assert!((h - expected).abs() <= 1e-9 * h, "h = √(μ·a(1−e²))");

        // (b) CROSS-CHECK threading semi_latus_rectum_from_elements (#384): h² = μ·p.
        for &(a, e) in &[(7.0e6_f64, 0.1_f64), (4.2e7, 0.0), (1.0e7, 0.3)] {
            let h = specific_angular_momentum_from_elements(a, e).unwrap();
            let p = semi_latus_rectum_from_elements(a, e).unwrap();
            assert!((h * h - MU_EARTH * p).abs() <= 1e-9 * (MU_EARTH * p), "h² = μ·p");
        }

        // (c) CIRCULAR cross-check threading circular_angular_momentum (#348): at e = 0,
        // h = √(μ·a) = circular_angular_momentum(a).
        for &a in &[7.0e6_f64, 4.2e7, 2.0e7] {
            assert!(
                (specific_angular_momentum_from_elements(a, 0.0).unwrap()
                    - circular_angular_momentum(a).unwrap())
                .abs()
                    <= 1e-9 * circular_angular_momentum(a).unwrap(),
                "h(a, 0) = circular_angular_momentum(a)"
            );
        }

        // (d) PERIAPSIS cross-check threading orbital_speed: h = r_p · v_p (position and
        // velocity are perpendicular at the apsis), with r_p = a(1−e).
        let (a, e) = (1.0e7_f64, 0.3_f64);
        let r_p = a * (1.0 - e);
        let h = specific_angular_momentum_from_elements(a, e).unwrap();
        assert!(
            (h - r_p * orbital_speed(r_p, a).unwrap()).abs() <= 1e-9 * h,
            "h = r_p · v_p at periapsis"
        );

        // (e) Err on non-physical input.
        assert!(specific_angular_momentum_from_elements(-1.0, 0.1).is_err());
        assert!(specific_angular_momentum_from_elements(7.0e6, 1.0).is_err());
        assert!(specific_angular_momentum_from_elements(f64::NAN, 0.1).is_err());
    }

    #[test]
    fn orbital_radius_from_speed_inverts_vis_viva() {
        // Round-trips orbital_speed for bound (r, a) pairs with r ≤ 2a.
        let r0 = R_EARTH + 400_000.0;
        for &(r, a) in &[(r0, r0), (r0, 1.2 * r0), (1.0e7, 2.0e7), (1.0e7, 1.0e7)] {
            let v = orbital_speed(r, a).unwrap();
            assert!(
                (orbital_radius_from_speed(v, a).unwrap() - r).abs() <= 1e-9 * r,
                "round-trip r = {r}"
            );
        }

        // Reduces to the circular case (a = r): the circular speed recovers r.
        let r = 7.0e6;
        assert!(
            (orbital_radius_from_speed(circular_speed(r).unwrap(), r).unwrap() - r).abs()
                <= 1e-9 * r,
            "circular: a = r"
        );

        // Monotonic: faster speed ⇒ smaller radius (vis-viva is decreasing in r).
        let a = 1.5e7;
        let v_slow = orbital_speed(1.2e7, a).unwrap();
        let v_fast = orbital_speed(8.0e6, a).unwrap();
        assert!(
            orbital_radius_from_speed(v_fast, a).unwrap()
                < orbital_radius_from_speed(v_slow, a).unwrap(),
            "faster → smaller radius"
        );

        // Err on non-physical input.
        assert!(orbital_radius_from_speed(7000.0, 0.0).is_err());
        assert!(orbital_radius_from_speed(7000.0, -1.0).is_err());
        assert!(orbital_radius_from_speed(-1.0, 1.0e7).is_err());
        assert!(orbital_radius_from_speed(f64::NAN, 1.0e7).is_err());
    }

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
    fn gravitational_potential_is_the_well_depth() {
        for &r in &[R_EARTH, R_EARTH + 400_000.0, 4.2164e7] {
            let phi = gravitational_potential(r).unwrap();

            // Threads specific_orbital_energy_from_state: ε = ½v² + Φ.
            let v = 5000.0;
            assert!(
                (specific_orbital_energy_from_state(v, r).unwrap() - (0.5 * v * v + phi)).abs()
                    <= 1e-9 * (0.5 * v * v + phi).abs(),
                "ε = ½v² + Φ"
            );

            // Threads escape_speed: v_esc = √(−2Φ).
            assert!(
                (escape_speed(r).unwrap() - (-2.0 * phi).sqrt()).abs()
                    <= 1e-9 * escape_speed(r).unwrap(),
                "v_esc = √(−2Φ)"
            );

            // Worked closed form Φ = −μ/r, and it is negative (a potential well).
            assert!((phi - (-MU_EARTH / r)).abs() <= 1e-9 * phi.abs(), "Φ = −μ/r");
            assert!(phi < 0.0, "potential well is negative");
        }

        // Monotonic: deeper (more negative) closer in, climbing toward 0 with r.
        assert!(
            gravitational_potential(R_EARTH).unwrap() < gravitational_potential(4.2164e7).unwrap(),
            "deeper well closer to Earth"
        );

        // Err on non-physical radius.
        assert!(gravitational_potential(0.0).is_err());
        assert!(gravitational_potential(-1.0).is_err());
        assert!(gravitational_potential(f64::NAN).is_err());
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
