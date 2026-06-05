//! Full 3-D orbital mechanics: classical orbital elements and a
//! two-body + J2 propagator.
//!
//! Where [`crate::orbit`] handles the *planar* ascent state, this module
//! works in the full three-dimensional Earth-centred inertial (ECI)
//! frame and carries the complete set of **classical orbital elements**
//! (COE): semi-major axis, eccentricity, inclination, right-ascension of
//! the ascending node (RAAN), argument of periapsis and true anomaly.
//!
//! It provides exact state ↔ element conversions (round-trip stable) and
//! an RK4 propagator with optional **J2 oblateness** — the dominant
//! perturbation in low Earth orbit, which makes the node regress and the
//! line of apsides rotate. The J2 secular rates are also given in closed
//! form so the propagator can be validated against them.
//!
//! Scope: this is the point-mass orbital layer (Phase 1 of the
//! launch-vehicle roadmap). It is not a full force model — no drag,
//! third-body, SRP, or higher-order geopotential yet — but J2 alone
//! captures the first-order LEO secular behaviour to good accuracy.

use std::f64::consts::TAU;

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::constants::{EARTH_ORBITAL_RATE, J2_EARTH, MU_EARTH, R_EARTH};
use crate::error::AstroError;
use crate::sim::check_step_count;

/// A 3-D inertial state vector (position + velocity).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StateVector {
    /// Position in the ECI frame (m).
    pub position: Vector3<f64>,
    /// Velocity in the ECI frame (m/s).
    pub velocity: Vector3<f64>,
}

/// Classical (Keplerian) orbital elements.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ClassicalElements {
    /// Semi-major axis (m).
    pub semi_major_axis: f64,
    /// Eccentricity (dimensionless).
    pub eccentricity: f64,
    /// Inclination (rad), in `[0, π]`.
    pub inclination: f64,
    /// Right ascension of the ascending node (rad), in `[0, 2π)`.
    pub raan: f64,
    /// Argument of periapsis (rad), in `[0, 2π)`.
    pub arg_periapsis: f64,
    /// True anomaly (rad), in `[0, 2π)`.
    pub true_anomaly: f64,
}

impl ClassicalElements {
    /// Apoapsis radius from Earth's centre (m).
    pub fn apoapsis_radius(&self) -> f64 {
        self.semi_major_axis * (1.0 + self.eccentricity)
    }

    /// Periapsis radius from Earth's centre (m).
    pub fn periapsis_radius(&self) -> f64 {
        self.semi_major_axis * (1.0 - self.eccentricity)
    }

    /// Orbital period (s) for a bound orbit (`a > 0`), else `None`.
    pub fn period(&self) -> Option<f64> {
        if self.semi_major_axis > 0.0 {
            Some(TAU * (self.semi_major_axis.powi(3) / MU_EARTH).sqrt())
        } else {
            None
        }
    }
}

/// Convert an inertial state vector to classical orbital elements.
///
/// Uses the standard angular-momentum / eccentricity / node-vector
/// construction. Circular and equatorial degeneracies fall back to the
/// usual conventions (undefined angles set to 0) so the result is always
/// finite.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] for input that would
/// otherwise yield silent `NaN`/`Inf` elements: a zero or non-finite
/// position (`r_mag` drives `μ/r`), a degenerate angular momentum
/// (`h_mag = 0`, i.e. rectilinear/zero motion, which makes the
/// inclination `NaN`), or the parabolic energy singularity.
pub fn rv_to_coe(state: &StateVector) -> Result<ClassicalElements, AstroError> {
    rv_to_coe_mu(state, MU_EARTH)
}

/// As [`rv_to_coe`] for an arbitrary central-body `μ`.
///
/// # Errors
///
/// As [`rv_to_coe`], plus rejects a non-finite or non-positive `mu`.
pub fn rv_to_coe_mu(state: &StateVector, mu: f64) -> Result<ClassicalElements, AstroError> {
    if !mu.is_finite() || mu <= 0.0 {
        return Err(AstroError::NonPhysicalState("mu must be finite and > 0"));
    }
    let r = state.position;
    let v = state.velocity;
    if !r.x.is_finite() || !r.y.is_finite() || !r.z.is_finite() {
        return Err(AstroError::NonPhysicalState("position must be finite"));
    }
    if !v.x.is_finite() || !v.y.is_finite() || !v.z.is_finite() {
        return Err(AstroError::NonPhysicalState("velocity must be finite"));
    }
    let r_mag = r.norm();
    if r_mag <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "position magnitude must be > 0",
        ));
    }
    let h_mag = r.cross(&v).norm();
    if h_mag <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "degenerate orbit: zero angular momentum (rectilinear motion)",
        ));
    }
    // a = -μ/(2ε); ε ≈ 0 (parabolic) blows the semi-major axis up.
    let energy = v.norm_squared() / 2.0 - mu / r_mag;
    if energy.abs() < f64::EPSILON {
        return Err(AstroError::NonPhysicalState(
            "parabolic energy singularity (specific energy ≈ 0)",
        ));
    }
    Ok(rv_to_coe_mu_unchecked(state, mu))
}

/// Element-recovery core without input validation. Internal use only,
/// for callers that pass a state already known to be finite with a
/// non-degenerate, non-parabolic orbit (e.g. an insertion or propagated
/// state). Mirrors the public [`rv_to_coe_mu`] math exactly.
pub(crate) fn rv_to_coe_mu_unchecked(state: &StateVector, mu: f64) -> ClassicalElements {
    let r = state.position;
    let v = state.velocity;
    let r_mag = r.norm();
    let v_mag = v.norm();

    let h = r.cross(&v);
    let h_mag = h.norm();

    // Node vector n = k × h.
    let k = Vector3::new(0.0, 0.0, 1.0);
    let node = k.cross(&h);
    let node_mag = node.norm();

    // Eccentricity vector.
    let e_vec = ((v_mag * v_mag - mu / r_mag) * r - r.dot(&v) * v) / mu;
    let ecc = e_vec.norm();

    let energy = v_mag * v_mag / 2.0 - mu / r_mag;
    // For non-parabolic orbits a = -μ/(2ε).
    let semi_major_axis = -mu / (2.0 * energy);

    let inclination = (h.z / h_mag).clamp(-1.0, 1.0).acos();

    // RAAN.
    let raan = if node_mag > 1e-12 {
        let mut o = (node.x / node_mag).clamp(-1.0, 1.0).acos();
        if node.y < 0.0 {
            o = TAU - o;
        }
        o
    } else {
        0.0 // equatorial: node undefined
    };

    // Argument of periapsis.
    let arg_periapsis = if node_mag > 1e-12 && ecc > 1e-12 {
        let mut w = (node.dot(&e_vec) / (node_mag * ecc))
            .clamp(-1.0, 1.0)
            .acos();
        if e_vec.z < 0.0 {
            w = TAU - w;
        }
        w
    } else {
        0.0 // circular or equatorial
    };

    // True anomaly.
    let true_anomaly = if ecc > 1e-12 {
        let mut nu = (e_vec.dot(&r) / (ecc * r_mag)).clamp(-1.0, 1.0).acos();
        if r.dot(&v) < 0.0 {
            nu = TAU - nu;
        }
        nu
    } else if node_mag > 1e-12 {
        // Circular inclined: argument of latitude from the node.
        let mut u = (node.dot(&r) / (node_mag * r_mag)).clamp(-1.0, 1.0).acos();
        if r.z < 0.0 {
            u = TAU - u;
        }
        u
    } else {
        // Circular equatorial: true longitude.
        let mut l = (r.x / r_mag).clamp(-1.0, 1.0).acos();
        if r.y < 0.0 {
            l = TAU - l;
        }
        l
    };

    ClassicalElements {
        semi_major_axis,
        eccentricity: ecc,
        inclination,
        raan,
        arg_periapsis,
        true_anomaly,
    }
}

/// Convert classical orbital elements to an inertial state vector.
///
/// # Errors
///
/// Returns [`AstroError::NonPhysicalState`] if the semi-latus rectum
/// `p = a(1 − e²)` is non-finite or non-positive (e.g. `e ≥ 1` with
/// `a > 0`, or non-finite elements), which would otherwise make the
/// perifocal speed `√(μ/p)` and the radius `NaN`/`Inf`.
pub fn coe_to_rv(coe: &ClassicalElements) -> Result<StateVector, AstroError> {
    coe_to_rv_mu(coe, MU_EARTH)
}

/// As [`coe_to_rv`] for an arbitrary central-body `μ`.
///
/// # Errors
///
/// As [`coe_to_rv`], plus rejects a non-finite or non-positive `mu`.
pub fn coe_to_rv_mu(coe: &ClassicalElements, mu: f64) -> Result<StateVector, AstroError> {
    if !mu.is_finite() || mu <= 0.0 {
        return Err(AstroError::NonPhysicalState("mu must be finite and > 0"));
    }
    let p = coe.semi_major_axis * (1.0 - coe.eccentricity * coe.eccentricity);
    if !p.is_finite() || p <= 0.0 {
        return Err(AstroError::NonPhysicalState(
            "semi-latus rectum p = a(1 - e²) must be finite and > 0",
        ));
    }
    Ok(coe_to_rv_mu_unchecked(coe, mu))
}

/// State-construction core without input validation. Internal use only,
/// for callers that pass elements already known to have a finite,
/// positive semi-latus rectum. Mirrors the public [`coe_to_rv_mu`] math
/// exactly.
pub(crate) fn coe_to_rv_mu_unchecked(coe: &ClassicalElements, mu: f64) -> StateVector {
    let p = coe.semi_major_axis * (1.0 - coe.eccentricity * coe.eccentricity);
    let (snu, cnu) = coe.true_anomaly.sin_cos();
    let r = p / (1.0 + coe.eccentricity * cnu);

    // Perifocal (PQW) frame.
    let r_pqw = Vector3::new(r * cnu, r * snu, 0.0);
    let sqrt_mu_p = (mu / p).sqrt();
    let v_pqw = Vector3::new(-sqrt_mu_p * snu, sqrt_mu_p * (coe.eccentricity + cnu), 0.0);

    // Rotate PQW -> ECI via Rz(Ω) Rx(i) Rz(ω).
    let (so, co) = coe.raan.sin_cos();
    let (si, ci) = coe.inclination.sin_cos();
    let (sw, cw) = coe.arg_periapsis.sin_cos();

    // Combined rotation matrix rows.
    let r11 = co * cw - so * sw * ci;
    let r12 = -co * sw - so * cw * ci;
    let r21 = so * cw + co * sw * ci;
    let r22 = -so * sw + co * cw * ci;
    let r31 = sw * si;
    let r32 = cw * si;

    let rotate = |v: Vector3<f64>| {
        Vector3::new(
            r11 * v.x + r12 * v.y,
            r21 * v.x + r22 * v.y,
            r31 * v.x + r32 * v.y,
        )
    };

    StateVector {
        position: rotate(r_pqw),
        velocity: rotate(v_pqw),
    }
}

/// Two-body (point-mass) gravitational acceleration in 3-D (m/s²).
pub fn two_body_accel(position: Vector3<f64>) -> Vector3<f64> {
    let r = position.norm();
    if r < 1.0 {
        return Vector3::zeros();
    }
    -MU_EARTH / (r * r * r) * position
}

/// J2 oblateness perturbing acceleration in the ECI frame (m/s²).
pub fn j2_accel(position: Vector3<f64>) -> Vector3<f64> {
    let r = position.norm();
    if r < 1.0 {
        return Vector3::zeros();
    }
    let (x, y, z) = (position.x, position.y, position.z);
    let factor = -1.5 * J2_EARTH * MU_EARTH * R_EARTH * R_EARTH / r.powi(5);
    let zr2 = 5.0 * z * z / (r * r);
    Vector3::new(
        factor * x * (1.0 - zr2),
        factor * y * (1.0 - zr2),
        factor * z * (3.0 - zr2),
    )
}

/// Propagate a 3-D state forward by `steps` RK4 steps of size `dt`,
/// optionally including the J2 perturbation.
///
/// # Errors
///
/// Returns [`AstroError::OutOfRange`] if `steps` exceeds
/// [`crate::sim::MAX_SIM_STEPS`].
pub fn propagate(
    state: &StateVector,
    dt: f64,
    steps: u64,
    with_j2: bool,
) -> Result<StateVector, AstroError> {
    check_step_count(steps)?;
    let accel = |pos: Vector3<f64>| {
        if with_j2 {
            two_body_accel(pos) + j2_accel(pos)
        } else {
            two_body_accel(pos)
        }
    };

    let mut s = *state;
    for _ in 0..steps {
        // RK4 on (r, v) with v̇ = a(r), ṙ = v.
        let k1r = s.velocity;
        let k1v = accel(s.position);
        let k2r = s.velocity + 0.5 * dt * k1v;
        let k2v = accel(s.position + 0.5 * dt * k1r);
        let k3r = s.velocity + 0.5 * dt * k2v;
        let k3v = accel(s.position + 0.5 * dt * k2r);
        let k4r = s.velocity + dt * k3v;
        let k4v = accel(s.position + dt * k3r);
        s.position += dt / 6.0 * (k1r + 2.0 * k2r + 2.0 * k3r + k4r);
        s.velocity += dt / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
    }
    Ok(s)
}

/// Mean motion `n = √(μ/a³)` and semi-latus rectum `p = a(1−e²)` for a
/// COE, returned only when both are well-defined for a **bound** orbit
/// (`a > 0` finite, `p > 0` finite). The secular J2 rates are physically
/// defined only for closed orbits, and a hand-built non-elliptic / non-
/// physical element set would otherwise drive `√(μ/a³)` or `(R⊕/p)²` to
/// `NaN`/`Inf`.
fn bound_n_and_p(coe: &ClassicalElements) -> Option<(f64, f64)> {
    let a = coe.semi_major_axis;
    if !a.is_finite() || a <= 0.0 {
        return None;
    }
    let p = a * (1.0 - coe.eccentricity * coe.eccentricity);
    if !p.is_finite() || p <= 0.0 {
        return None;
    }
    Some(((MU_EARTH / a.powi(3)).sqrt(), p))
}

/// Secular J2 rate of change of the RAAN (rad/s) — the nodal regression.
///
/// `dΩ/dt = -1.5 · n · J2 · (R⊕/p)² · cos i`, where `n = √(μ/a³)` and
/// `p = a(1−e²)`. Negative for prograde orbits (the node drifts west).
///
/// Returns `0.0` for a non-elliptic or non-physical element set
/// (`a ≤ 0` or `p ≤ 0`), for which the secular rate is undefined — rather
/// than the silent `NaN`/`Inf` the raw `√(μ/a³)` / `(R⊕/p)²` would give.
pub fn j2_raan_rate(coe: &ClassicalElements) -> f64 {
    let Some((n, p)) = bound_n_and_p(coe) else {
        return 0.0;
    };
    -1.5 * n * J2_EARTH * (R_EARTH / p).powi(2) * coe.inclination.cos()
}

/// Secular J2 rate of change of the argument of periapsis (rad/s).
///
/// `dω/dt = 1.5 · n · J2 · (R⊕/p)² · (2 − 2.5 sin²i)`.
///
/// Returns `0.0` for a non-elliptic or non-physical element set
/// (`a ≤ 0` or `p ≤ 0`), for which the secular rate is undefined — rather
/// than the silent `NaN`/`Inf` the raw `√(μ/a³)` / `(R⊕/p)²` would give.
pub fn j2_arg_periapsis_rate(coe: &ClassicalElements) -> f64 {
    let Some((n, p)) = bound_n_and_p(coe) else {
        return 0.0;
    };
    let si = coe.inclination.sin();
    1.5 * n * J2_EARTH * (R_EARTH / p).powi(2) * (2.0 - 2.5 * si * si)
}

/// The **sun-synchronous inclination** (rad) for an orbit of semi-major axis
/// `semi_major_axis` (m) and `eccentricity` — the inclination at which the J2
/// nodal regression [`j2_raan_rate`] exactly matches Earth's mean orbital rate
/// about the Sun ([`crate::constants::EARTH_ORBITAL_RATE`]), so the orbit plane
/// holds a fixed angle to the Sun and the ground track repeats at the same
/// local solar time.
///
/// Inverts the secular nodal-rate condition
/// `Ω̇ = −1.5·n·J2·(R⊕/p)²·cos i = Ω̇_sun` for `cos i`, giving
/// `i = arccos( −Ω̇_sun / [1.5·n·J2·(R⊕/p)²] )` with `n = √(μ/a³)` and
/// `p = a(1−e²)`. Because the required drift is prograde (eastward) the cosine
/// is negative, so the inclination is always **retrograde** (`> 90°`) — the
/// familiar ≈ 98° of a low-Earth sun-sync orbit.
///
/// Returns `None` when the elements are not a bound orbit (`a` not positive
/// finite, `e ∉ [0, 1)`), or when the geometry cannot reach the solar rate
/// (`|cos i| > 1`, i.e. an orbit too high for J2 to precess fast enough), where
/// no real inclination exists.
pub fn sun_synchronous_inclination(semi_major_axis: f64, eccentricity: f64) -> Option<f64> {
    let a = semi_major_axis;
    if !a.is_finite() || a <= 0.0 || !(0.0..1.0).contains(&eccentricity) {
        return None;
    }
    let p = a * (1.0 - eccentricity * eccentricity);
    if !p.is_finite() || p <= 0.0 {
        return None;
    }
    let n = (MU_EARTH / a.powi(3)).sqrt();
    let denom = 1.5 * n * J2_EARTH * (R_EARTH / p).powi(2);
    let cos_i = -EARTH_ORBITAL_RATE / denom;
    if !cos_i.is_finite() || cos_i.abs() > 1.0 {
        return None;
    }
    Some(cos_i.acos())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Smallest signed angular difference `a − b` wrapped to `(−π, π]`.
    fn angle_diff(a: f64, b: f64) -> f64 {
        let mut d = a - b;
        while d > PI {
            d -= TAU;
        }
        while d <= -PI {
            d += TAU;
        }
        d
    }

    fn iss_like() -> ClassicalElements {
        ClassicalElements {
            semi_major_axis: R_EARTH + 420_000.0,
            eccentricity: 0.001,
            inclination: 51.6_f64.to_radians(),
            raan: 0.3,
            arg_periapsis: 0.7,
            true_anomaly: 1.0,
        }
    }

    #[test]
    fn coe_rv_round_trip() {
        let coe = iss_like();
        let rv = coe_to_rv(&coe).expect("valid coe");
        let back = rv_to_coe(&rv).expect("valid state");
        assert!((back.semi_major_axis - coe.semi_major_axis).abs() < 1.0);
        assert!((back.eccentricity - coe.eccentricity).abs() < 1e-9);
        assert!(angle_diff(back.inclination, coe.inclination).abs() < 1e-9);
        assert!(angle_diff(back.raan, coe.raan).abs() < 1e-9);
        assert!(angle_diff(back.arg_periapsis, coe.arg_periapsis).abs() < 1e-7);
        assert!(angle_diff(back.true_anomaly, coe.true_anomaly).abs() < 1e-7);
    }

    #[test]
    fn inclination_recovered_from_state() {
        // A 60° inclined circular orbit must report i = 60°.
        let coe = ClassicalElements {
            semi_major_axis: R_EARTH + 500_000.0,
            eccentricity: 0.0,
            inclination: 60.0_f64.to_radians(),
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let rv = coe_to_rv(&coe).expect("valid coe");
        let back = rv_to_coe(&rv).expect("valid state");
        assert!((back.inclination.to_degrees() - 60.0).abs() < 1e-6);
    }

    #[test]
    fn two_body_propagation_conserves_energy_in_3d() {
        let coe = iss_like();
        let s0 = coe_to_rv(&coe).expect("valid coe");
        let e0 = s0.velocity.norm_squared() / 2.0 - MU_EARTH / s0.position.norm();
        let period = coe.period().unwrap();
        let dt = 1.0;
        let steps = (period / dt).round() as u64;
        let s1 = propagate(&s0, dt, steps, false).expect("valid step count");
        let e1 = s1.velocity.norm_squared() / 2.0 - MU_EARTH / s1.position.norm();
        assert!((e1 - e0).abs() / e0.abs() < 1e-6, "energy {e0} -> {e1}");
        // Returns near the start after one period.
        assert!((s1.position - s0.position).norm() < 5_000.0);
    }

    #[test]
    fn j2_nodal_regression_matches_analytic_rate() {
        // Propagate several orbits with J2 and compare the measured RAAN
        // drift to the closed-form secular rate.
        let coe = iss_like();
        let s0 = coe_to_rv(&coe).expect("valid coe");
        let period = coe.period().unwrap();
        let n_orbits = 5.0;
        let dt = 1.0;
        let steps = (n_orbits * period / dt).round() as u64;
        let s1 = propagate(&s0, dt, steps, true).expect("valid step count");
        let coe1 = rv_to_coe(&s1).expect("valid state");

        let elapsed = steps as f64 * dt;
        let measured_rate = angle_diff(coe1.raan, coe.raan) / elapsed;
        let analytic = j2_raan_rate(&coe);

        // Both must be negative (westward regression for a prograde orbit).
        assert!(
            analytic < 0.0 && measured_rate < 0.0,
            "rates {analytic} {measured_rate}"
        );
        let rel = (measured_rate - analytic).abs() / analytic.abs();
        assert!(rel < 0.05, "J2 RAAN rate off by {:.1}%", rel * 100.0);
    }

    #[test]
    fn j2_leaves_inclination_and_sma_secularly_unchanged() {
        // J2 has no secular effect on a or i — only periodic wobble.
        let coe = iss_like();
        let s0 = coe_to_rv(&coe).expect("valid coe");
        let period = coe.period().unwrap();
        let dt = 1.0;
        let steps = (3.0 * period / dt).round() as u64;
        let s1 = propagate(&s0, dt, steps, true).expect("valid step count");
        let coe1 = rv_to_coe(&s1).expect("valid state");
        assert!((coe1.inclination - coe.inclination).abs() < 1e-4);
        assert!((coe1.semi_major_axis - coe.semi_major_axis).abs() < 2_000.0);
    }

    #[test]
    fn propagate_rejects_absurd_step_count() {
        // u64::MAX steps would hang; expect a clean Err.
        let s0 = coe_to_rv(&iss_like()).expect("valid coe");
        let r = propagate(&s0, 1.0, u64::MAX, false);
        assert!(
            matches!(r, Err(AstroError::OutOfRange { what: "steps", .. })),
            "u64::MAX steps must be rejected, got {r:?}"
        );
        assert!(propagate(&s0, 1.0, 10, true).is_ok());
    }

    #[test]
    fn sun_synchronous_inclination_gives_expected_regression() {
        // A ~98° inclination LEO is sun-synchronous: the node should
        // regress eastward (positive) at roughly 360°/year ≈ 1.99e-7 rad/s.
        let coe = ClassicalElements {
            semi_major_axis: R_EARTH + 700_000.0,
            eccentricity: 0.001,
            inclination: 98.0_f64.to_radians(),
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let rate = j2_raan_rate(&coe);
        // Retrograde-inclination orbit -> positive nodal drift.
        assert!(rate > 0.0, "rate {rate}");
        let deg_per_day = rate.to_degrees() * 86_400.0;
        // Sun-sync target is ~+0.986°/day; this geometry lands in range.
        assert!((deg_per_day - 0.986).abs() < 0.5, "{deg_per_day} deg/day");
    }

    #[test]
    fn sun_synchronous_inclination_inverts_the_nodal_rate() {
        // The inverse of the regression test above: solve for the inclination
        // that makes the J2 nodal rate equal Earth's solar rate. A ~700 km
        // circular LEO is the textbook sun-synchronous case at i ≈ 98.2°.
        let a = R_EARTH + 700_000.0;
        let i = sun_synchronous_inclination(a, 0.0).expect("700 km is sun-syncable");
        let deg = i.to_degrees();
        assert!((deg - 98.2).abs() < 0.5, "sun-sync inclination ≈ 98.2°, got {deg}");
        // Sun-sync orbits are retrograde (cos i < 0 ⇒ i > 90°).
        assert!(deg > 90.0, "must be retrograde, got {deg}");

        // Round-trip: feeding that inclination back into j2_raan_rate reproduces
        // Earth's solar rate exactly (the defining sun-synchronous condition).
        let coe = ClassicalElements {
            semi_major_axis: a,
            eccentricity: 0.0,
            inclination: i,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        assert!(
            (j2_raan_rate(&coe) - EARTH_ORBITAL_RATE).abs() < 1e-15,
            "nodal rate must match the solar rate"
        );

        // A higher orbit needs a more retrograde (larger) inclination.
        let i_high =
            sun_synchronous_inclination(R_EARTH + 1_500_000.0, 0.0).expect("1500 km syncable");
        assert!(i_high > i, "higher orbit → larger inclination");

        // Too high for J2 to precess fast enough → no real inclination exists.
        assert!(sun_synchronous_inclination(R_EARTH + 10_000_000.0, 0.0).is_none());
        // Non-physical inputs → None (never a NaN).
        assert!(sun_synchronous_inclination(-1.0, 0.0).is_none());
        assert!(sun_synchronous_inclination(f64::NAN, 0.0).is_none());
        assert!(sun_synchronous_inclination(a, 1.5).is_none()); // unbound eccentricity
    }

    #[test]
    fn rv_to_coe_rejects_zero_position() {
        // r_mag = 0 -> mu/r_mag is Inf -> energy/e_vec come out NaN/Inf.
        let s = StateVector {
            position: Vector3::zeros(),
            velocity: Vector3::new(0.0, 7_800.0, 0.0),
        };
        let r = rv_to_coe(&s);
        assert!(
            matches!(r, Err(AstroError::NonPhysicalState(_))),
            "zero position must be rejected, got {r:?}"
        );
    }

    #[test]
    fn rv_to_coe_rejects_rectilinear_zero_angular_momentum() {
        // r ∥ v -> h = r × v = 0 -> inclination h.z/h_mag = 0/0 = NaN.
        let s = StateVector {
            position: Vector3::new(R_EARTH + 500_000.0, 0.0, 0.0),
            velocity: Vector3::new(100.0, 0.0, 0.0), // parallel to r
        };
        let r = rv_to_coe(&s);
        assert!(
            matches!(r, Err(AstroError::NonPhysicalState(_))),
            "rectilinear motion must be rejected, got {r:?}"
        );
    }

    #[test]
    fn rv_to_coe_rejects_non_finite_state() {
        let s = StateVector {
            position: Vector3::new(f64::NAN, 0.0, 0.0),
            velocity: Vector3::new(0.0, 7_800.0, 0.0),
        };
        assert!(matches!(
            rv_to_coe(&s),
            Err(AstroError::NonPhysicalState(_))
        ));
    }

    #[test]
    fn coe_to_rv_rejects_non_positive_semi_latus_rectum() {
        // e >= 1 with a > 0 -> p = a(1 - e²) <= 0 -> √(μ/p) is NaN/Inf.
        let parabolic = ClassicalElements {
            semi_major_axis: R_EARTH + 500_000.0,
            eccentricity: 1.0, // p = 0
            inclination: 0.5,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let r = coe_to_rv(&parabolic);
        assert!(
            matches!(r, Err(AstroError::NonPhysicalState(_))),
            "p = 0 must be rejected, got {r:?}"
        );
        let hyperbolic = ClassicalElements {
            eccentricity: 1.5, // p < 0 with a > 0
            ..parabolic
        };
        assert!(matches!(
            coe_to_rv(&hyperbolic),
            Err(AstroError::NonPhysicalState(_))
        ));
    }

    #[test]
    fn j2_rates_are_zero_no_op_for_non_physical_elements_not_nan() {
        // The secular J2 rates use n = √(μ/a³) and (R⊕/p)². A hand-built
        // element set with a non-positive semi-major axis or a non-positive
        // semi-latus rectum (e >= 1) drove those to NaN/Inf silently. They
        // now return 0 ("no secular drift for a non-elliptic orbit") — the
        // same a>0/p>0 well-definedness `period()`/`coe_to_rv` already use.
        let zero_a = ClassicalElements {
            semi_major_axis: 0.0,
            eccentricity: 0.1,
            inclination: 0.9,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        assert_eq!(j2_raan_rate(&zero_a), 0.0);
        assert_eq!(j2_arg_periapsis_rate(&zero_a), 0.0);

        let negative_a = ClassicalElements {
            semi_major_axis: -(R_EARTH + 500_000.0),
            ..zero_a
        };
        assert!(j2_raan_rate(&negative_a).is_finite());
        assert!(j2_arg_periapsis_rate(&negative_a).is_finite());
        assert_eq!(j2_raan_rate(&negative_a), 0.0);

        // e >= 1 with a > 0 -> p <= 0 -> (R/p)² is Inf.
        let para = ClassicalElements {
            semi_major_axis: R_EARTH + 500_000.0,
            eccentricity: 1.0,
            ..zero_a
        };
        assert_eq!(j2_raan_rate(&para), 0.0);
        assert_eq!(j2_arg_periapsis_rate(&para), 0.0);

        // A valid bound orbit still gets a real, non-zero rate.
        let leo = ClassicalElements {
            semi_major_axis: R_EARTH + 500_000.0,
            eccentricity: 0.001,
            inclination: 0.9,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        assert!(j2_raan_rate(&leo) != 0.0 && j2_raan_rate(&leo).is_finite());
    }
}
