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

    /// Orbital radius `r = a(1−e²)/(1 + e·cos ν)` (m) at true anomaly `nu` (rad) —
    /// the **polar equation of the conic**, the foundational relation behind the
    /// orbit's shape. Its `ν = 0` and `ν = π` values are exactly the periapsis and
    /// apoapsis radii; `ν = ±π/2` gives the semi-latus rectum `p = a(1−e²)`; a
    /// circular orbit (`e = 0`) returns `a` at every angle. (For an open orbit
    /// `e ≥ 1` it diverges to `±∞` at the asymptote `cos ν = −1/e`, as it should.)
    pub fn radius_at_true_anomaly(&self, nu: f64) -> f64 {
        let p = self.semi_major_axis * (1.0 - self.eccentricity * self.eccentricity);
        p / (1.0 + self.eccentricity * nu.cos())
    }

    /// The **outbound true anomaly** `ν = arccos((p/r − 1)/e)` (rad, in `[0, π]`)
    /// at which the orbit reaches radius `radius` `r` (m) — the inverse of
    /// [`radius_at_true_anomaly`](Self::radius_at_true_anomaly), which maps the
    /// other way (`ν → r`). It answers "where in the orbit does it pass through
    /// this radius?", the geometric basis of altitude-crossing and event timing
    /// (feed the result through the `ν → E → M → time` chain for the *when*). The
    /// inbound pass is the mirror `2π − ν`, since `r(ν) = r(−ν)`.
    ///
    /// `r = periapsis` gives `ν = 0`, `r = apoapsis` gives `ν = π`, and the
    /// semi-latus rectum `r = p = a(1−e²)` gives `ν = π/2`. Returns `None` when
    /// the orbit never reaches that radius (outside `[periapsis, apoapsis]`), and
    /// for input where the inverse is undefined — a circular orbit (`e = 0`, no
    /// apsides), an open orbit (`e ≥ 1`), or a non-finite `r`. The `arccos`
    /// argument is clamped to `[−1, 1]` so the apsidal boundaries are exact
    /// despite floating-point round-off.
    pub fn true_anomaly_at_radius(&self, radius: f64) -> Option<f64> {
        let e = self.eccentricity;
        if !radius.is_finite() || e <= 0.0 || e >= 1.0 {
            return None; // bound, non-circular orbits only (0 < e < 1)
        }
        let r_peri = self.semi_major_axis * (1.0 - e);
        let r_apo = self.semi_major_axis * (1.0 + e);
        if radius < r_peri - 1e-6 || radius > r_apo + 1e-6 {
            return None; // the orbit never reaches this radius
        }
        let p = self.semi_major_axis * (1.0 - e * e);
        let cos_nu = ((p / radius - 1.0) / e).clamp(-1.0, 1.0);
        Some(cos_nu.acos())
    }

    /// The orbital **velocity components** `(v_r, v_θ)` (m/s) at true anomaly
    /// `nu` (rad), in the rotating polar frame: the *radial* component
    /// `v_r = (μ/h)·e·sin ν` along the outward radius, and the *transverse*
    /// component `v_θ = (μ/h)·(1 + e·cos ν)` perpendicular to it (the direction
    /// of orbital motion), with `μ/h = √(μ/p)` and `p = a(1−e²)`.
    ///
    /// This is the velocity companion to
    /// [`radius_at_true_anomaly`](Self::radius_at_true_anomaly): together `r(ν)`
    /// and `(v_r, v_θ)(ν)` give the full in-plane state at any point of the
    /// orbit. The radial part vanishes at the apsides (`ν = 0, π`), where the
    /// motion is purely transverse and `v_θ` hits its extremes (fastest at
    /// periapsis, slowest at apoapsis); `v_r` is positive climbing out toward
    /// apoapsis and negative falling back in. The speed `√(v_r² + v_θ²)`
    /// reproduces the vis-viva law `√(μ(2/r − 1/a))`, and `v_r/v_θ` is the
    /// tangent of the flight-path angle. Uses Earth's `μ`; intended for closed
    /// orbits (`e < 1`).
    pub fn velocity_components_at_true_anomaly(&self, nu: f64) -> (f64, f64) {
        let e = self.eccentricity;
        let p = self.semi_major_axis * (1.0 - e * e);
        let mu_over_h = (MU_EARTH / p).sqrt(); // μ/h = √(μ/p)
        let v_r = mu_over_h * e * nu.sin();
        let v_theta = mu_over_h * (1.0 + e * nu.cos());
        (v_r, v_theta)
    }

    /// The orbital **specific angular momentum** `h = √(μ·a(1−e²))` (m²/s) —
    /// the angular momentum *per unit mass*, `h = |r × v|`, and the orbit's
    /// defining conserved quantity. It is constant everywhere along the path
    /// (Kepler's second law: the radius vector sweeps equal areas in equal
    /// times, at the areal rate `h/2`) and is the `μ/h` factor that scales the
    /// [`velocity_components_at_true_anomaly`](Self::velocity_components_at_true_anomaly).
    ///
    /// Equivalently `h = r·v_θ` at *every* true anomaly — orbital radius times
    /// transverse speed — so the large radius at apoapsis exactly offsets the
    /// small transverse speed there (and the reverse at periapsis), the product
    /// held fixed. In terms of the semi-latus rectum `p = a(1−e²)` it is simply
    /// `h = √(μ·p)`; a circular orbit (`e = 0`) gives `h = √(μ·a)`.
    ///
    /// Uses Earth's `μ`. Returns `None` for an orbit that is not bound and
    /// closed (`a ≤ 0`, or `e ≥ 1`), where the radicand `μ·a(1−e²)` is not a
    /// positive real and this closed form does not apply.
    pub fn specific_angular_momentum(&self) -> Option<f64> {
        let p = self.semi_major_axis * (1.0 - self.eccentricity * self.eccentricity);
        if self.semi_major_axis > 0.0 && p > 0.0 {
            Some((MU_EARTH * p).sqrt())
        } else {
            None
        }
    }

    /// The orbital **specific energy** `ε = −μ/(2a)` (J/kg) — the total orbital
    /// energy *per unit mass* (kinetic plus gravitational potential), the orbit's
    /// other conserved invariant alongside the
    /// [`specific_angular_momentum`](Self::specific_angular_momentum). It depends
    /// only on the semi-major axis, so a more energetic orbit is simply a *larger*
    /// one.
    ///
    /// Its sign classifies the conic: `ε < 0` for a **bound** ellipse (`a > 0`),
    /// `ε = 0` for the parabolic escape limit (`a → ∞`), and `ε > 0` for a
    /// **hyperbolic** flyby (`a < 0`) — so unlike the specific angular momentum it
    /// is meaningful for *every* orbit and is returned as a plain value, not an
    /// `Option`. It ties speed to radius through the **vis-viva** relation
    /// `½v² − μ/r = ε`, i.e. `v = √(μ(2/r − 1/a))`: the kinetic and potential terms
    /// trade off along the path while their sum stays fixed. Uses Earth's `μ`; the
    /// degenerate `a = 0` (a point orbit) gives `±∞`.
    pub fn specific_orbital_energy(&self) -> f64 {
        -MU_EARTH / (2.0 * self.semi_major_axis)
    }

    /// Solve **Kepler's equation** `M = E − e·sin E` for the eccentric anomaly
    /// `E` (rad) given the mean anomaly `mean_anomaly` `M` (rad), by
    /// Newton–Raphson iteration.
    ///
    /// This is the keystone of propagating an orbit *in time*: the mean anomaly
    /// `M = n·(t − t_p)` advances uniformly with time, but the geometry needs
    /// the eccentric anomaly `E`, and the link `M = E − e·sin E` is
    /// transcendental — it has no closed form and must be inverted numerically.
    /// With `E` in hand,
    /// [`true_anomaly_from_eccentric`](Self::true_anomaly_from_eccentric) gives
    /// `ν` and [`radius_at_true_anomaly`](Self::radius_at_true_anomaly) gives
    /// `r`, completing the `time → M → E → ν → r` chain. Newton's step
    /// `E ← E − (E − e·sin E − M)/(1 − e·cos E)` converges quadratically from the
    /// seed `E₀ = M + e·sin M`; the derivative `1 − e·cos E ≥ 1 − e > 0` never
    /// vanishes for a closed orbit, so the iteration is unconditionally stable.
    ///
    /// `M = 0` and `M = π` are fixed points (`E = 0`, `E = π`, where `sin E = 0`);
    /// a circular orbit (`e = 0`) collapses Kepler's equation to `E = M`.
    /// Defined for closed orbits (`0 ≤ e < 1`); a hyperbolic eccentricity
    /// (`e ≥ 1`), which needs the *hyperbolic* Kepler equation
    /// `M = e·sinh F − F`, or any non-finite `M`, yields `NaN`.
    pub fn eccentric_anomaly_from_mean(&self, mean_anomaly: f64) -> f64 {
        let e = self.eccentricity;
        if !(0.0..1.0).contains(&e) || !mean_anomaly.is_finite() {
            return f64::NAN;
        }
        let m = mean_anomaly;
        let mut ecc = m + e * m.sin(); // standard initial seed
        for _ in 0..64 {
            let delta = (ecc - e * ecc.sin() - m) / (1.0 - e * ecc.cos());
            ecc -= delta;
            if delta.abs() < 1e-14 {
                break;
            }
        }
        ecc
    }

    /// True anomaly `ν` (rad) from the **eccentric anomaly** `E` (rad) via the
    /// half-angle relation `ν = 2·atan2(√(1+e)·sin(E/2), √(1−e)·cos(E/2))`.
    ///
    /// This is the geometric half of propagating a Kepler orbit in time: once
    /// the time-driven *mean* anomaly `M = E − e·sin E` has been inverted for
    /// `E` (Kepler's equation), this turns `E` into the true anomaly `ν` that
    /// [`radius_at_true_anomaly`](Self::radius_at_true_anomaly) needs —
    /// completing the `time → M → E → ν → r` position-from-time chain (the
    /// reverse of the forward `ν → time` map). The half-angle `atan2` form is
    /// preferred over `cos ν = (cos E − e)/(1 − e·cos E)`
    /// because it resolves the quadrant directly, with no sign ambiguity past
    /// apoapsis.
    ///
    /// `E = 0` and `E = π` are fixed points (`ν = 0`, `ν = π`); a circular orbit
    /// (`e = 0`) collapses the map to the identity `ν = E`; and for a canonical
    /// `E ∈ [0, 2π)` the result lies in `[0, 2π)` and increases monotonically
    /// with `E`. (Defined for closed orbits, `e < 1`; the eccentric anomaly has
    /// no meaning for an open orbit, so `e ≥ 1` yields `NaN`.)
    pub fn true_anomaly_from_eccentric(&self, eccentric_anomaly: f64) -> f64 {
        let e = self.eccentricity;
        let half = eccentric_anomaly / 2.0;
        let y = (1.0 + e).sqrt() * half.sin();
        let x = (1.0 - e).sqrt() * half.cos();
        2.0 * y.atan2(x)
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
    fn radius_at_true_anomaly_traces_the_conic() {
        use std::f64::consts::PI;
        let coe = ClassicalElements {
            semi_major_axis: 7.0e6,
            eccentricity: 0.2,
            inclination: 0.0,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let (a, e) = (coe.semi_major_axis, coe.eccentricity);
        // ν=0 → periapsis a(1−e); ν=π → apoapsis a(1+e).
        assert!((coe.radius_at_true_anomaly(0.0) - coe.periapsis_radius()).abs() < 1e-6);
        assert!((coe.radius_at_true_anomaly(PI) - coe.apoapsis_radius()).abs() < 1e-6);
        assert!((coe.radius_at_true_anomaly(0.0) - a * (1.0 - e)).abs() < 1e-6, "perigee");
        assert!((coe.radius_at_true_anomaly(PI) - a * (1.0 + e)).abs() < 1e-6, "apogee");
        // ν=±π/2 → the semi-latus rectum p = a(1−e²).
        let p = a * (1.0 - e * e);
        assert!((coe.radius_at_true_anomaly(PI / 2.0) - p).abs() < 1e-6, "semi-latus rectum");
        // Symmetric about the line of apsides: r(ν) = r(−ν).
        assert!((coe.radius_at_true_anomaly(1.0) - coe.radius_at_true_anomaly(-1.0)).abs() < 1e-9);
        // A circular orbit has constant radius a at every true anomaly.
        let circ = ClassicalElements { eccentricity: 0.0, ..coe };
        for nu in [0.0, 1.0, PI, 2.5] {
            assert!((circ.radius_at_true_anomaly(nu) - a).abs() < 1e-6, "circular r at {nu}");
        }
    }

    #[test]
    fn true_anomaly_at_radius_inverts_the_orbit_equation() {
        use std::f64::consts::PI;
        let coe = ClassicalElements {
            semi_major_axis: 7.0e6,
            eccentricity: 0.2,
            inclination: 0.0,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let (a, e) = (coe.semi_major_axis, coe.eccentricity);
        // Apsides and the semi-latus rectum are the fixed points of the inverse.
        // (The apsidal checks use 1e-6: acos is ill-conditioned at cos = ±1, so a
        // ~1e-16 round-off in the argument shows up as a ~1e-8 angle error.)
        assert!(
            coe.true_anomaly_at_radius(coe.periapsis_radius()).unwrap().abs() < 1e-6,
            "perigee → ν=0"
        );
        assert!(
            (coe.true_anomaly_at_radius(coe.apoapsis_radius()).unwrap() - PI).abs() < 1e-6,
            "apogee → ν=π"
        );
        let p = a * (1.0 - e * e);
        assert!(
            (coe.true_anomaly_at_radius(p).unwrap() - PI / 2.0).abs() < 1e-9,
            "semi-latus rectum → ν=π/2"
        );
        // Round-trips with radius_at_true_anomaly (#144) across the outbound half.
        for nu in [0.3_f64, 1.0, 2.0, PI] {
            let r = coe.radius_at_true_anomaly(nu);
            let nu_back = coe.true_anomaly_at_radius(r).expect("r is reachable");
            // 1e-6: the ν=π round-trip touches the acos cos=−1 boundary.
            assert!((nu_back - nu).abs() < 1e-6, "round trip at ν={nu}: got {nu_back}");
        }
        // Radii the orbit never reaches → None.
        assert!(
            coe.true_anomaly_at_radius(coe.periapsis_radius() * 0.5).is_none(),
            "below perigee"
        );
        assert!(
            coe.true_anomaly_at_radius(coe.apoapsis_radius() * 2.0).is_none(),
            "above apogee"
        );
        // Undefined for a circular, hyperbolic, or non-finite case → None.
        assert!(
            ClassicalElements { eccentricity: 0.0, ..coe }.true_anomaly_at_radius(a).is_none(),
            "circular"
        );
        assert!(
            ClassicalElements { eccentricity: 1.5, ..coe }.true_anomaly_at_radius(a).is_none(),
            "hyperbolic"
        );
        assert!(coe.true_anomaly_at_radius(f64::NAN).is_none(), "non-finite r");
    }

    #[test]
    fn true_anomaly_from_eccentric_inverts_keplers_geometry() {
        use std::f64::consts::PI;
        let coe = ClassicalElements {
            semi_major_axis: 7.0e6,
            eccentricity: 0.5,
            inclination: 0.0,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let e = coe.eccentricity;
        // Apsides are fixed points: E=0 → ν=0, E=π → ν=π.
        assert!(coe.true_anomaly_from_eccentric(0.0).abs() < 1e-12, "E=0 → ν=0");
        assert!((coe.true_anomaly_from_eccentric(PI) - PI).abs() < 1e-12, "E=π → ν=π");
        // Worked point e=0.5, E=π/2 → ν=2π/3, cross-checked against the standard
        // cos ν = (cos E − e)/(1 − e·cos E) = −0.5.
        let nu = coe.true_anomaly_from_eccentric(PI / 2.0);
        assert!((nu - 2.0 * PI / 3.0).abs() < 1e-12, "e=0.5, E=π/2 → ν=2π/3, got {nu}");
        assert!((nu.cos() + 0.5).abs() < 1e-12, "cos ν = −0.5 cross-check");
        // Agrees with the closed cos-form at arbitrary E, stays in [0,2π), and
        // shares E's quadrant (same sign of sine).
        for ea in [0.3_f64, 1.0, 2.0, 3.0, 4.5, 6.0] {
            let nu = coe.true_anomaly_from_eccentric(ea);
            let cos_nu = (ea.cos() - e) / (1.0 - e * ea.cos());
            assert!((nu.cos() - cos_nu).abs() < 1e-9, "cos ν at E={ea}");
            assert!((0.0..TAU).contains(&nu), "ν in [0,2π) at E={ea}: {nu}");
            assert!(nu.sin() * ea.sin() >= 0.0, "ν shares E's quadrant at E={ea}");
        }
        // A circular orbit collapses the map to ν = E.
        let circ = ClassicalElements { eccentricity: 0.0, ..coe };
        for ea in [0.0_f64, 0.5, PI / 2.0, PI, 4.0, 5.5] {
            assert!((circ.true_anomaly_from_eccentric(ea) - ea).abs() < 1e-12, "circular ν=E at {ea}");
        }
    }

    #[test]
    fn eccentric_anomaly_from_mean_solves_keplers_equation() {
        use std::f64::consts::PI;
        let coe = ClassicalElements {
            semi_major_axis: 7.0e6,
            eccentricity: 0.3,
            inclination: 0.0,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let e = coe.eccentricity;
        // The returned E satisfies Kepler's equation M = E − e·sin E to tight tol.
        for m in [0.1_f64, 1.0, 2.5, 4.0, 5.9] {
            let ea = coe.eccentric_anomaly_from_mean(m);
            assert!((ea - e * ea.sin() - m).abs() < 1e-12, "Kepler residual at M={m}");
        }
        // Apsides are fixed points: M=0 → E=0, M=π → E=π (sin E = 0 there).
        assert!(coe.eccentric_anomaly_from_mean(0.0).abs() < 1e-12, "M=0 → E=0");
        assert!((coe.eccentric_anomaly_from_mean(PI) - PI).abs() < 1e-12, "M=π → E=π");
        // A circular orbit (e=0) collapses Kepler's equation to E = M.
        let circ = ClassicalElements { eccentricity: 0.0, ..coe };
        for m in [0.0_f64, 0.7, PI, 4.2] {
            assert!((circ.eccentric_anomaly_from_mean(m) - m).abs() < 1e-12, "circular E=M at {m}");
        }
        // Round-trips with the forward map M = E − e·sin E even at high eccentricity.
        let ecc = ClassicalElements { eccentricity: 0.9, ..coe };
        for e_true in [0.2_f64, 1.5, 3.0, 5.0] {
            let m = e_true - 0.9 * e_true.sin();
            assert!(
                (ecc.eccentric_anomaly_from_mean(m) - e_true).abs() < 1e-10,
                "round trip E={e_true}"
            );
        }
        // Closes the position-from-time chain: M → E → ν → r is finite, r > 0,
        // and the radius lies between perigee and apogee.
        let ea = coe.eccentric_anomaly_from_mean(1.0);
        let nu = coe.true_anomaly_from_eccentric(ea);
        let r = coe.radius_at_true_anomaly(nu);
        assert!(ea.is_finite() && nu.is_finite() && r.is_finite() && r > 0.0, "chain finite r>0");
        assert!(
            r >= coe.periapsis_radius() - 1.0 && r <= coe.apoapsis_radius() + 1.0,
            "r in [r_p, r_a]: {r}"
        );
        // Out of domain: a hyperbolic eccentricity and a non-finite M → NaN.
        let hyp = ClassicalElements { eccentricity: 1.5, ..coe };
        assert!(hyp.eccentric_anomaly_from_mean(1.0).is_nan(), "e≥1 → NaN");
        assert!(coe.eccentric_anomaly_from_mean(f64::NAN).is_nan(), "NaN M → NaN");
    }

    #[test]
    fn velocity_components_at_true_anomaly_decompose_the_orbital_speed() {
        use std::f64::consts::PI;
        let coe = ClassicalElements {
            semi_major_axis: 8.0e6,
            eccentricity: 0.25,
            inclination: 0.0,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let (a, e) = (coe.semi_major_axis, coe.eccentricity);
        let mu = MU_EARTH;
        // At perigee (ν=0): purely transverse, v_θ = √(μ/a·(1+e)/(1−e)) (max speed).
        let (vr0, vt0) = coe.velocity_components_at_true_anomaly(0.0);
        assert!(vr0.abs() < 1e-6, "v_r = 0 at perigee, got {vr0}");
        assert!((vt0 - (mu / a * (1.0 + e) / (1.0 - e)).sqrt()).abs() < 1e-3, "perigee speed");
        // At apogee (ν=π): purely transverse, v_θ = √(μ/a·(1−e)/(1+e)) (min speed).
        let (vrp, vtp) = coe.velocity_components_at_true_anomaly(PI);
        assert!(vrp.abs() < 1e-6, "v_r = 0 at apogee, got {vrp}");
        assert!((vtp - (mu / a * (1.0 - e) / (1.0 + e)).sqrt()).abs() < 1e-3, "apogee speed");
        assert!(vt0 > vtp, "faster at perigee than apogee");
        // The speed √(v_r²+v_θ²) reproduces vis-viva μ(2/r − 1/a) at the matching r.
        for nu in [0.3_f64, 1.0, 2.0, PI, 4.0, 5.5] {
            let (vr, vt) = coe.velocity_components_at_true_anomaly(nu);
            let speed_sq = vr * vr + vt * vt;
            let r = coe.radius_at_true_anomaly(nu);
            let vis_viva_sq = mu * (2.0 / r - 1.0 / a);
            assert!((speed_sq - vis_viva_sq).abs() / vis_viva_sq < 1e-12, "vis-viva at ν={nu}");
        }
        // Radial velocity is positive climbing out (0<ν<π) and negative falling in.
        assert!(coe.velocity_components_at_true_anomaly(1.0).0 > 0.0, "outbound v_r > 0");
        assert!(coe.velocity_components_at_true_anomaly(4.0).0 < 0.0, "inbound v_r < 0");
        // A circular orbit: no radial velocity, constant transverse speed √(μ/a).
        let circ = ClassicalElements { eccentricity: 0.0, ..coe };
        let v_circ = (mu / a).sqrt();
        for nu in [0.0_f64, 1.0, PI, 4.2] {
            let (vr, vt) = circ.velocity_components_at_true_anomaly(nu);
            assert!(vr.abs() < 1e-9, "circular v_r = 0 at {nu}");
            assert!((vt - v_circ).abs() < 1e-3, "circular v_θ = √(μ/a) at {nu}");
        }
    }

    #[test]
    fn specific_angular_momentum_is_the_conserved_h() {
        use std::f64::consts::PI;
        let coe = ClassicalElements {
            semi_major_axis: 8.0e6,
            eccentricity: 0.25,
            inclination: 0.0,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let (a, e) = (coe.semi_major_axis, coe.eccentricity);
        let h = coe.specific_angular_momentum().expect("bound orbit has h");
        // Closed form h = √(μ·a(1−e²)) = √(μ·p) (relative tol: h ≈ 5.5e10 m²/s).
        let expected = (MU_EARTH * a * (1.0 - e * e)).sqrt();
        assert!((h - expected).abs() / h < 1e-9, "closed form, got {h}");
        // Conserved: h = r·v_θ at *every* true anomaly — cross-checks #162's
        // velocity components and the conic radius (Kepler's second law).
        for nu in [0.0_f64, 0.7, PI / 2.0, 2.0, PI, 4.5, 5.7] {
            let r = coe.radius_at_true_anomaly(nu);
            let (_, v_theta) = coe.velocity_components_at_true_anomaly(nu);
            assert!((r * v_theta - h).abs() / h < 1e-12, "h = r·v_θ at ν={nu}");
        }
        // Apsidal form: h = r_peri·v_peri = r_apo·v_apo.
        let (_, vt_peri) = coe.velocity_components_at_true_anomaly(0.0);
        let (_, vt_apo) = coe.velocity_components_at_true_anomaly(PI);
        assert!((coe.periapsis_radius() * vt_peri - h).abs() / h < 1e-12, "h at periapsis");
        assert!((coe.apoapsis_radius() * vt_apo - h).abs() / h < 1e-12, "h at apoapsis");
        // A circular orbit (e=0): h = √(μ·a).
        let circ = ClassicalElements { eccentricity: 0.0, ..coe };
        let h_circ = circ.specific_angular_momentum().expect("circular orbit has h");
        assert!((h_circ - (MU_EARTH * a).sqrt()).abs() / h_circ < 1e-9, "circular h = √(μ·a)");
        // Same a, higher eccentricity carries less angular momentum (smaller p).
        let ecc = ClassicalElements { eccentricity: 0.6, ..coe };
        assert!(ecc.specific_angular_momentum().unwrap() < h, "more eccentric → smaller h");
        // Not bound/closed → None.
        let para = ClassicalElements { eccentricity: 1.0, ..coe };
        assert!(para.specific_angular_momentum().is_none(), "parabolic (e=1) → None");
        let hyp = ClassicalElements { eccentricity: 1.5, ..coe };
        assert!(hyp.specific_angular_momentum().is_none(), "hyperbolic (e>1) → None");
        let neg_a = ClassicalElements { semi_major_axis: -8.0e6, ..coe };
        assert!(neg_a.specific_angular_momentum().is_none(), "a ≤ 0 → None");
    }

    #[test]
    fn specific_orbital_energy_is_minus_mu_over_2a() {
        use std::f64::consts::PI;
        let coe = ClassicalElements {
            semi_major_axis: 8.0e6,
            eccentricity: 0.25,
            inclination: 0.0,
            raan: 0.0,
            arg_periapsis: 0.0,
            true_anomaly: 0.0,
        };
        let a = coe.semi_major_axis;
        let eps = coe.specific_orbital_energy();
        // Closed form ε = −μ/(2a) (relative tol; ε ≈ −2.49e7 J/kg).
        assert!(
            (eps - (-MU_EARTH / (2.0 * a))).abs() / eps.abs() < 1e-12,
            "closed form, got {eps}"
        );
        // A bound ellipse has negative specific energy.
        assert!(eps < 0.0, "bound orbit ε < 0");
        // Vis-viva cross-check: ½·v(ν)² − μ/r(ν) = ε at every true anomaly — ties
        // the energy to #162's velocity components and the conic radius.
        for nu in [0.0_f64, 0.7, PI / 2.0, 2.0, PI, 4.5, 5.7] {
            let r = coe.radius_at_true_anomaly(nu);
            let (v_r, v_theta) = coe.velocity_components_at_true_anomaly(nu);
            let speed_sq = v_r * v_r + v_theta * v_theta;
            let e_at_nu = 0.5 * speed_sq - MU_EARTH / r;
            assert!((e_at_nu - eps).abs() / eps.abs() < 1e-9, "vis-viva at ν={nu}");
        }
        // A larger orbit is more energetic (ε rises toward 0 as a grows).
        let bigger = ClassicalElements { semi_major_axis: 2.0e7, ..coe };
        assert!(
            bigger.specific_orbital_energy() > eps,
            "larger a → higher (less negative) energy"
        );
        // A hyperbolic orbit (a < 0) has positive specific energy (unbound).
        let hyper = ClassicalElements { semi_major_axis: -8.0e6, ..coe };
        assert!(hyper.specific_orbital_energy() > 0.0, "hyperbolic ε > 0");
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
