//! Quasi-static Archimedes buoyancy on a floating rigid body.
//!
//! Given a [wave field](crate::wave::OceanWaveField) and a rigid body, this
//! module computes the **submerged volume**, the **buoyant force**
//! `F_b = ρ_water · g · V_submerged` acting upward through the **centre of
//! buoyancy** (the centroid of the submerged volume), a simple **linear +
//! quadratic drag**, and the resulting **net force and torque**. A small
//! heave/pitch/roll integrator advances the body state under those forces.
//!
//! ## ⚠ Honesty / scope — quasi-static buoyancy, NOT seakeeping CFD
//!
//! This is a **Gerstner-wave + quasi-static (hydrostatic) buoyancy** model, the
//! kind used in games and first-cut engineering. It is **NOT** a seakeeping
//! RANS/CFD solver. Specifically:
//!
//! * Buoyancy is computed from the **instantaneous wave height under each sample
//!   point** (Froude–Krylov / hydrostatic pressure only). There is **no
//!   diffraction** (the body does not scatter the incident wave) and **no
//!   radiation** (the body's own motion does not generate outgoing waves).
//! * **Added mass and wave-radiation damping are not modelled** beyond the lumped
//!   linear+quadratic drag coefficients you supply — there is no
//!   frequency-dependent added-mass/`B(ω)` from a boundary-element/strip-theory
//!   solve.
//! * The waterplane restoring stiffness uses the **still-waterplane area**; the
//!   pitch/roll metacentric behaviour follows from the sample geometry, not a
//!   computed `GM`/`BM` from the hull offsets.
//!
//! What *is* pinned in the tests are the **exact hydrostatic checks** Archimedes'
//! principle must obey: a fully-submerged body feels `ρ g V`; a freely floating
//! body settles displacing its own weight; the small-amplitude heave oscillates
//! at the analytic natural frequency `ω_n = sqrt(ρ g A_wp / m)`; and the net
//! force is zero at equilibrium.

use crate::error::OceanError;
use crate::wave::OceanWaveField;
use nalgebra::Vector3;

/// A rigid body described by **sample volume points**.
///
/// The body is discretised into `points`, each carrying a small volume `volume`
/// and a body-frame position `position` (relative to the body's reference point,
/// typically its centre of mass). The submerged volume at a given pose and wave
/// state is the sum of the volumes of the points that lie below the water
/// surface; the centre of buoyancy is their volume-weighted centroid.
///
/// This representation handles arbitrary shapes and partial / asymmetric
/// immersion (the basis for pitch and roll restoring), at the cost of being a
/// quadrature: accuracy improves with more, smaller sample volumes.
#[derive(Debug, Clone)]
pub struct SampleBody {
    /// Body-frame sample positions (m), relative to the reference point.
    positions: Vec<Vector3<f64>>,
    /// Per-sample volumes (m³); each `> 0`. `positions` and `volumes` are
    /// parallel and equal-length.
    volumes: Vec<f64>,
    /// Total mass of the body (kg), `> 0`.
    mass: f64,
}

impl SampleBody {
    /// Build a sample-point body.
    ///
    /// * `samples` — `(position, volume)` pairs; every `volume` must be `> 0`
    ///   and finite, every `position` finite. At least one sample is required.
    /// * `mass` — total body mass (kg), `> 0` and finite.
    ///
    /// # Errors
    ///
    /// [`OceanError::InvalidConfig`] for an empty sample set, a non-positive
    /// `volume`/`mass`; [`OceanError::NonFinite`] for any non-finite input.
    pub fn new(samples: &[(Vector3<f64>, f64)], mass: f64) -> Result<Self, OceanError> {
        if !mass.is_finite() {
            return Err(OceanError::NonFinite(format!("mass = {mass}")));
        }
        if mass <= 0.0 {
            return Err(OceanError::InvalidConfig(format!(
                "mass must be > 0, got {mass}"
            )));
        }
        if samples.is_empty() {
            return Err(OceanError::InvalidConfig(
                "a SampleBody needs at least one sample point".to_string(),
            ));
        }
        let mut positions = Vec::with_capacity(samples.len());
        let mut volumes = Vec::with_capacity(samples.len());
        for (i, (p, v)) in samples.iter().enumerate() {
            if !p.iter().all(|c| c.is_finite()) {
                return Err(OceanError::NonFinite(format!("sample[{i}].position")));
            }
            if !v.is_finite() {
                return Err(OceanError::NonFinite(format!("sample[{i}].volume = {v}")));
            }
            if *v <= 0.0 {
                return Err(OceanError::InvalidConfig(format!(
                    "sample[{i}].volume must be > 0, got {v}"
                )));
            }
            positions.push(*p);
            volumes.push(*v);
        }
        Ok(Self {
            positions,
            volumes,
            mass,
        })
    }

    /// Total mass (kg).
    #[must_use]
    pub fn mass(&self) -> f64 {
        self.mass
    }

    /// Total (fully-immersed) volume `Σ vᵢ` (m³).
    #[must_use]
    pub fn total_volume(&self) -> f64 {
        self.volumes.iter().sum()
    }

    /// Number of sample points.
    #[must_use]
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Always `false` (construction requires at least one sample); provided to
    /// satisfy the `len`/`is_empty` clippy pairing.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Submerged volume and centre of buoyancy at a pose.
    ///
    /// The body reference point is at world position `ref_pos`; `roll`/`pitch`
    /// are small rotations (rad) about the world `+x`/`+z` axes applied to the
    /// body-frame sample positions (a small-angle pose suitable for the
    /// quasi-static model — heave dominates, pitch/roll are perturbations).
    /// Each rotated sample's world height is compared to the wave height there;
    /// fully-below samples count their whole volume.
    ///
    /// Returns `(v_submerged, centre_of_buoyancy_world)`. When nothing is
    /// submerged the centre of buoyancy is reported as `ref_pos` (the value is
    /// irrelevant because the volume — and hence the force — is zero).
    #[must_use]
    pub fn submerged(
        &self,
        ref_pos: &Vector3<f64>,
        roll: f64,
        pitch: f64,
        field: &OceanWaveField,
        t: f64,
    ) -> (f64, Vector3<f64>) {
        let (sr, cr) = roll.sin_cos();
        let (sp, cp) = pitch.sin_cos();
        let mut v_sub = 0.0;
        let mut moment = Vector3::zeros();
        for (p, &vol) in self.positions.iter().zip(&self.volumes) {
            // Small-angle: rotate about x (roll) then about z (pitch).
            // Rx: y,z mix; Rz: x,y mix. Applied in that order.
            let y1 = p.y * cr - p.z * sr;
            let z1 = p.y * sr + p.z * cr;
            let x2 = p.x * cp - y1 * sp;
            let y2 = p.x * sp + y1 * cp;
            let world = ref_pos + Vector3::new(x2, y2, z1);
            let surface = field.height_at(world.x, world.z, t);
            if world.y < surface {
                v_sub += vol;
                moment += vol * world;
            }
        }
        if v_sub <= 0.0 {
            (0.0, *ref_pos)
        } else {
            (v_sub, moment / v_sub)
        }
    }
}

/// A rigid body described by an **analytic prismatic hull**: a constant
/// waterplane (cross-sectional) area extruded to a fixed height.
///
/// This is the simplest body for which the buoyancy/heave physics has a closed
/// form, so it backs the analytic equilibrium and heave-frequency benchmarks.
/// The submerged volume as a function of the immersion depth `s` (how far the
/// keel is below the local water surface) is
///
/// ```text
/// V_sub(s) = A_wp · clamp(s, 0, height)
/// ```
///
/// i.e. linear in `s` until fully immersed, then capped at the total volume
/// `A_wp · height`. The restoring stiffness in heave is therefore exactly
/// `ρ g A_wp`, giving the heave natural frequency `ω_n = sqrt(ρ g A_wp / m)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HullBody {
    /// Still-waterplane (cross-sectional) area `A_wp` (m²), `> 0`.
    waterplane_area: f64,
    /// Vertical extent of the prism (m), `> 0`. The total volume is
    /// `waterplane_area * height`.
    height: f64,
    /// Total mass (kg), `> 0`.
    mass: f64,
}

impl HullBody {
    /// Build a prismatic hull.
    ///
    /// # Errors
    ///
    /// [`OceanError::InvalidConfig`] if `waterplane_area`, `height`, or `mass`
    /// is non-positive; [`OceanError::NonFinite`] for any non-finite input.
    pub fn new(waterplane_area: f64, height: f64, mass: f64) -> Result<Self, OceanError> {
        for (name, v) in [
            ("waterplane_area", waterplane_area),
            ("height", height),
            ("mass", mass),
        ] {
            if !v.is_finite() {
                return Err(OceanError::NonFinite(format!("{name} = {v}")));
            }
            if v <= 0.0 {
                return Err(OceanError::InvalidConfig(format!(
                    "{name} must be > 0, got {v}"
                )));
            }
        }
        Ok(Self {
            waterplane_area,
            height,
            mass,
        })
    }

    /// Waterplane area `A_wp` (m²).
    #[must_use]
    pub fn waterplane_area(&self) -> f64 {
        self.waterplane_area
    }

    /// Prism height (m).
    #[must_use]
    pub fn height(&self) -> f64 {
        self.height
    }

    /// Total volume `A_wp · height` (m³).
    #[must_use]
    pub fn total_volume(&self) -> f64 {
        self.waterplane_area * self.height
    }

    /// Total mass (kg).
    #[must_use]
    pub fn mass(&self) -> f64 {
        self.mass
    }

    /// Submerged volume for an immersion depth `s` (keel depth below the local
    /// surface, m): `A_wp · clamp(s, 0, height)`.
    #[must_use]
    pub fn submerged_volume(&self, immersion: f64) -> f64 {
        self.waterplane_area * immersion.clamp(0.0, self.height)
    }

    /// The **equilibrium draft** `T = m / (ρ A_wp)` (m): the immersion at which
    /// the displaced weight equals the body weight.
    ///
    /// `ρ A_wp > 0` by construction, so the divide is safe. The result is **not**
    /// clamped to `height`; a draft exceeding `height` means the body is too
    /// dense to float (it would sink), which the caller can detect by comparing
    /// to [`HullBody::height`].
    ///
    /// # Errors
    ///
    /// [`OceanError::InvalidConfig`] if `water_density` is non-positive;
    /// [`OceanError::NonFinite`] if it is not finite.
    pub fn equilibrium_draft(&self, water_density: f64) -> Result<f64, OceanError> {
        check_density(water_density)?;
        Ok(self.mass / (water_density * self.waterplane_area))
    }

    /// The analytic **heave natural frequency** `ω_n = sqrt(ρ g A_wp / m)`
    /// (rad/s) for small oscillations about the floating waterline.
    ///
    /// Valid only while the body is partially immersed (the linear regime); once
    /// fully submerged or fully emerged the restoring force saturates and this no
    /// longer applies. `m > 0` by construction.
    ///
    /// # Errors
    ///
    /// [`OceanError::InvalidConfig`] if `water_density` or `gravity` is
    /// non-positive; [`OceanError::NonFinite`] if either is not finite.
    pub fn heave_natural_frequency(
        &self,
        water_density: f64,
        gravity: f64,
    ) -> Result<f64, OceanError> {
        check_density(water_density)?;
        check_gravity(gravity)?;
        Ok((water_density * gravity * self.waterplane_area / self.mass).sqrt())
    }
}

/// Linear + quadratic translational drag coefficients (per axis the same).
///
/// The drag force is `F_d = −(c_lin · v + c_quad · |v| · v)`, the standard
/// low-order hydrodynamic damping form: a linear (skin-friction / radiation-
/// proxy) term dominant at low speed and a quadratic (form-drag) term dominant
/// at high speed. Both coefficients must be `≥ 0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Drag {
    /// Linear coefficient `c_lin ≥ 0` (N·s/m).
    pub linear: f64,
    /// Quadratic coefficient `c_quad ≥ 0` (N·s²/m²).
    pub quadratic: f64,
}

impl Drag {
    /// Build a drag model.
    ///
    /// # Errors
    ///
    /// [`OceanError::InvalidConfig`] if either coefficient is negative;
    /// [`OceanError::NonFinite`] for a non-finite input.
    pub fn new(linear: f64, quadratic: f64) -> Result<Self, OceanError> {
        for (name, v) in [("linear", linear), ("quadratic", quadratic)] {
            if !v.is_finite() {
                return Err(OceanError::NonFinite(format!("{name} = {v}")));
            }
            if v < 0.0 {
                return Err(OceanError::InvalidConfig(format!(
                    "{name} drag must be >= 0, got {v}"
                )));
            }
        }
        Ok(Self { linear, quadratic })
    }

    /// Zero drag.
    #[must_use]
    pub fn none() -> Self {
        Self {
            linear: 0.0,
            quadratic: 0.0,
        }
    }

    /// Drag force for a velocity `v`: `−(c_lin v + c_quad |v| v)`.
    #[must_use]
    pub fn force(&self, velocity: &Vector3<f64>) -> Vector3<f64> {
        let speed = velocity.norm();
        -(self.linear * velocity + self.quadratic * speed * velocity)
    }
}

/// The dynamic state of a floating body: world position and velocity of the
/// reference point plus small roll/pitch angles and their rates.
///
/// `position.y` is the heave coordinate. `roll` is rotation about world `+x`,
/// `pitch` about world `+z` (matching [`SampleBody::submerged`]). Yaw is not
/// modelled (a quasi-static buoyancy model has no yaw restoring).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyState {
    /// World position of the reference point (m).
    pub position: Vector3<f64>,
    /// World linear velocity of the reference point (m/s).
    pub velocity: Vector3<f64>,
    /// Roll angle about world `+x` (rad).
    pub roll: f64,
    /// Pitch angle about world `+z` (rad).
    pub pitch: f64,
    /// Roll rate (rad/s).
    pub roll_rate: f64,
    /// Pitch rate (rad/s).
    pub pitch_rate: f64,
}

impl BodyState {
    /// A body at rest at `position` with zero attitude and rates.
    #[must_use]
    pub fn at_rest(position: Vector3<f64>) -> Self {
        Self {
            position,
            velocity: Vector3::zeros(),
            roll: 0.0,
            pitch: 0.0,
            roll_rate: 0.0,
            pitch_rate: 0.0,
        }
    }
}

/// The net loads on the body at an instant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Loads {
    /// Net force on the reference point (N): buoyancy + gravity + drag.
    pub force: Vector3<f64>,
    /// Net torque about the reference point (N·m), from the buoyant force acting
    /// at the (offset) centre of buoyancy.
    pub torque: Vector3<f64>,
    /// Submerged volume at this instant (m³).
    pub submerged_volume: f64,
    /// Centre of buoyancy in world coordinates (m).
    pub center_of_buoyancy: Vector3<f64>,
}

/// The full configuration of a buoyancy simulation: the body, the environment
/// constants, and the drag model.
#[derive(Debug, Clone)]
pub struct BuoyancySim {
    body: SampleBody,
    water_density: f64,
    gravity: f64,
    drag: Drag,
    /// Rotational damping coefficient applied to roll/pitch rates (N·m·s/rad),
    /// `≥ 0`. Lumps rotational radiation/viscous damping.
    angular_drag: f64,
    /// Moment of inertia used for the roll/pitch integrator (kg·m²), `> 0`.
    /// A single scalar (isotropic in roll/pitch) keeps the quasi-static model
    /// simple; supply the body's representative `I` about a horizontal axis.
    inertia: f64,
}

impl BuoyancySim {
    /// Assemble a simulation.
    ///
    /// # Errors
    ///
    /// [`OceanError::InvalidConfig`] if `water_density`, `gravity`, or `inertia`
    /// is non-positive, or `angular_drag` is negative; [`OceanError::NonFinite`]
    /// for any non-finite environment constant.
    pub fn new(
        body: SampleBody,
        water_density: f64,
        gravity: f64,
        drag: Drag,
        angular_drag: f64,
        inertia: f64,
    ) -> Result<Self, OceanError> {
        check_density(water_density)?;
        check_gravity(gravity)?;
        if !inertia.is_finite() {
            return Err(OceanError::NonFinite(format!("inertia = {inertia}")));
        }
        if inertia <= 0.0 {
            return Err(OceanError::InvalidConfig(format!(
                "inertia must be > 0, got {inertia}"
            )));
        }
        if !angular_drag.is_finite() {
            return Err(OceanError::NonFinite(format!(
                "angular_drag = {angular_drag}"
            )));
        }
        if angular_drag < 0.0 {
            return Err(OceanError::InvalidConfig(format!(
                "angular_drag must be >= 0, got {angular_drag}"
            )));
        }
        Ok(Self {
            body,
            water_density,
            gravity,
            drag,
            angular_drag,
            inertia,
        })
    }

    /// The body.
    #[must_use]
    pub fn body(&self) -> &SampleBody {
        &self.body
    }

    /// Water density (kg/m³).
    #[must_use]
    pub fn water_density(&self) -> f64 {
        self.water_density
    }

    /// Gravity (m/s²).
    #[must_use]
    pub fn gravity(&self) -> f64 {
        self.gravity
    }

    /// Compute the instantaneous [`Loads`] on the body for a given state, wave
    /// field, and time.
    ///
    /// Buoyancy: `F_b = ρ g V_sub` upward at the centre of buoyancy. Gravity:
    /// `−m g` at the reference point. Drag: the [`Drag`] model on the linear
    /// velocity. Torque: `(cob − ref) × F_b` (the buoyant force's moment about
    /// the reference point) — this is what rights the body in pitch and roll.
    #[must_use]
    pub fn loads(&self, state: &BodyState, field: &OceanWaveField, t: f64) -> Loads {
        let (v_sub, cob) = self
            .body
            .submerged(&state.position, state.roll, state.pitch, field, t);
        let buoyant = Vector3::new(0.0, self.water_density * self.gravity * v_sub, 0.0);
        let weight = Vector3::new(0.0, -self.body.mass * self.gravity, 0.0);
        let drag = self.drag.force(&state.velocity);
        let force = buoyant + weight + drag;
        // Torque of the buoyant force about the reference point.
        let arm = cob - state.position;
        let torque = arm.cross(&buoyant);
        Loads {
            force,
            torque,
            submerged_volume: v_sub,
            center_of_buoyancy: cob,
        }
    }

    /// Advance the body one semi-implicit (symplectic) Euler step of `dt`
    /// seconds under the current loads.
    ///
    /// Translational: `v += (F/m) dt; x += v dt`. Rotational (roll about world
    /// `+x` from the torque's `x` component, pitch about world `+z` from the
    /// torque's `z` component), with the lumped scalar inertia and angular drag:
    /// `ω += (τ − c_ang ω)/I · dt; θ += ω dt`. Symplectic Euler is first-order
    /// (an `O(dt)` bias) and energy-stable for the oscillator, the same choice
    /// the SPH fluids solver makes.
    ///
    /// # Errors
    ///
    /// [`OceanError::InvalidConfig`] if `dt <= 0`; [`OceanError::NonFinite`] if
    /// `dt` is not finite or the resulting state would be non-finite (e.g. a
    /// caller fed a non-finite wave field / state).
    pub fn step(
        &self,
        state: &BodyState,
        field: &OceanWaveField,
        t: f64,
        dt: f64,
    ) -> Result<BodyState, OceanError> {
        if !dt.is_finite() {
            return Err(OceanError::NonFinite(format!("dt = {dt}")));
        }
        if dt <= 0.0 {
            return Err(OceanError::InvalidConfig(format!(
                "dt must be > 0, got {dt}"
            )));
        }
        let loads = self.loads(state, field, t);
        let mut next = *state;

        // Linear (semi-implicit Euler): m > 0 by construction.
        let accel = loads.force / self.body.mass;
        next.velocity += accel * dt;
        next.position += next.velocity * dt;

        // Rotational about x (roll) and z (pitch); I > 0 by construction.
        let roll_torque = loads.torque.x - self.angular_drag * state.roll_rate;
        let pitch_torque = loads.torque.z - self.angular_drag * state.pitch_rate;
        next.roll_rate += roll_torque / self.inertia * dt;
        next.pitch_rate += pitch_torque / self.inertia * dt;
        next.roll += next.roll_rate * dt;
        next.pitch += next.pitch_rate * dt;

        // Fail loud if numerics blew up (NaN wave field, etc.).
        let finite = next.position.iter().all(|v| v.is_finite())
            && next.velocity.iter().all(|v| v.is_finite())
            && next.roll.is_finite()
            && next.pitch.is_finite()
            && next.roll_rate.is_finite()
            && next.pitch_rate.is_finite();
        if !finite {
            return Err(OceanError::NonFinite(
                "integrated body state became non-finite".to_string(),
            ));
        }
        Ok(next)
    }
}

/// Buoyant force magnitude `ρ g V` (N) for a fully-known submerged volume.
///
/// A free function for the simplest Archimedes check: a body of submerged volume
/// `v_submerged` in fluid of density `ρ` under gravity `g` feels `ρ g V` upward.
///
/// # Errors
///
/// [`OceanError::InvalidConfig`] if `density`, `gravity`, or `v_submerged` is
/// negative (a zero submerged volume is allowed → zero force);
/// [`OceanError::NonFinite`] for any non-finite input.
pub fn archimedes_force(density: f64, gravity: f64, v_submerged: f64) -> Result<f64, OceanError> {
    check_density(density)?;
    check_gravity(gravity)?;
    if !v_submerged.is_finite() {
        return Err(OceanError::NonFinite(format!(
            "v_submerged = {v_submerged}"
        )));
    }
    if v_submerged < 0.0 {
        return Err(OceanError::InvalidConfig(format!(
            "v_submerged must be >= 0, got {v_submerged}"
        )));
    }
    Ok(density * gravity * v_submerged)
}

/// Validate a density argument (`> 0`, finite).
fn check_density(density: f64) -> Result<(), OceanError> {
    if !density.is_finite() {
        return Err(OceanError::NonFinite(format!("density = {density}")));
    }
    if density <= 0.0 {
        return Err(OceanError::InvalidConfig(format!(
            "density must be > 0, got {density}"
        )));
    }
    Ok(())
}

/// Validate a gravity argument (`> 0`, finite).
fn check_gravity(gravity: f64) -> Result<(), OceanError> {
    if !gravity.is_finite() {
        return Err(OceanError::NonFinite(format!("gravity = {gravity}")));
    }
    if gravity <= 0.0 {
        return Err(OceanError::InvalidConfig(format!(
            "gravity must be > 0, got {gravity}"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::wave::{OceanWaveField, SEAWATER_DENSITY, STANDARD_GRAVITY};

    const G: f64 = STANDARD_GRAVITY;
    const RHO: f64 = SEAWATER_DENSITY;

    /// A unit cube of `n×n×n` sample points spanning `[-half, half]³`, each
    /// carrying an equal share of the total volume `side³`.
    fn cube_body(side: f64, n: usize, mass: f64) -> SampleBody {
        let half = side / 2.0;
        let total = side * side * side;
        let per = total / (n * n * n) as f64;
        let mut samples = Vec::new();
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let f = |idx: usize| -half + side * (idx as f64 + 0.5) / n as f64;
                    samples.push((Vector3::new(f(i), f(j), f(k)), per));
                }
            }
        }
        SampleBody::new(&samples, mass).unwrap()
    }

    // ---- BENCHMARK-PIN (2a): fully-submerged body feels rho*g*V exactly ----

    #[test]
    fn fully_submerged_force_is_rho_g_v() {
        // archimedes_force is the exact statement.
        let v = 2.5;
        let f = archimedes_force(RHO, G, v).unwrap();
        assert!((f - RHO * G * v).abs() < 1e-9);

        // And the sample body, pushed deep below a flat sea, submerges its whole
        // volume so the buoyant force equals rho*g*V_total.
        let side = 1.0;
        let body = cube_body(side, 4, 800.0);
        let sim = BuoyancySim::new(body, RHO, G, Drag::none(), 0.0, 1.0).unwrap();
        let field = OceanWaveField::flat(0.0).unwrap();
        // Reference point 10 m down: every sample is well below the surface.
        let state = BodyState::at_rest(Vector3::new(0.0, -10.0, 0.0));
        let loads = sim.loads(&state, &field, 0.0);
        let v_total = sim.body().total_volume();
        assert!((loads.submerged_volume - v_total).abs() < 1e-12);
        assert!((loads.force.y - (RHO * G * v_total - sim.body().mass() * G)).abs() < 1e-9);
    }

    #[test]
    fn archimedes_rejects_bad_args() {
        assert!(archimedes_force(0.0, G, 1.0).is_err());
        assert!(archimedes_force(RHO, 0.0, 1.0).is_err());
        assert!(archimedes_force(RHO, G, -1.0).is_err());
        assert!(archimedes_force(f64::NAN, G, 1.0).is_err());
        assert_eq!(archimedes_force(RHO, G, 0.0).unwrap(), 0.0);
    }

    // ---- BENCHMARK-PIN (2b): floating equilibrium displaces own weight ----

    #[test]
    fn hull_equilibrium_displaces_own_weight() {
        let hull = HullBody::new(2.0, 3.0, 1500.0).unwrap();
        let draft = hull.equilibrium_draft(RHO).unwrap();
        // V_sub * rho * g == m * g  ⟺  rho * V_sub == m.
        let v_sub = hull.submerged_volume(draft);
        assert!((RHO * v_sub - hull.mass()).abs() < 1e-9);
        // Buoyant force balances weight.
        let fb = archimedes_force(RHO, G, v_sub).unwrap();
        assert!((fb - hull.mass() * G).abs() < 1e-6);
        // And the draft formula matches m/(rho A).
        assert!((draft - hull.mass() / (RHO * hull.waterplane_area())).abs() < 1e-12);
    }

    // ---- BENCHMARK-PIN (3): heave natural frequency ----

    #[test]
    fn hull_heave_natural_frequency_formula() {
        let hull = HullBody::new(2.0, 5.0, 1500.0).unwrap();
        let wn = hull.heave_natural_frequency(RHO, G).unwrap();
        let expected = (RHO * G * hull.waterplane_area() / hull.mass()).sqrt();
        assert!((wn - expected).abs() < 1e-12);
    }

    #[test]
    fn released_near_equilibrium_oscillates_at_heave_frequency() {
        // A prismatic sample body floating in a FLAT sea, displaced in heave and
        // released, must oscillate at omega_n = sqrt(rho g A / m).
        //
        // A finite set of sample points makes the submerged volume a *staircase*
        // function of depth (a whole layer enters the water at once), so the
        // restoring force only approximates the smooth `rho g A * ds` ramp when
        // the layer thickness is « the heave amplitude. We therefore (a) sample
        // the waterplane as a SINGLE column (exact for pure heave in a flat sea,
        // where every column is at the same height), and (b) make the layers very
        // thin relative to the heave amplitude so the staircase ≈ the ramp.
        let area = 1.0_f64; // 1 m^2 waterplane (one column carries it all)
        let height = 4.0;
        let mass = 500.0; // draft = m/(rho A) ≈ 0.488 m  (rho≈1025)
        let nz = 1000usize; // layer thickness = 4 mm
        let per = (area * height) / nz as f64;
        let mut samples = Vec::new();
        for i in 0..nz {
            let y = -height / 2.0 + height * (i as f64 + 0.5) / nz as f64;
            // single central column: full waterplane area, fine in y
            samples.push((Vector3::new(0.0, y, 0.0), per));
        }
        let body = SampleBody::new(&samples, mass).unwrap();
        let sim = BuoyancySim::new(body, RHO, G, Drag::none(), 0.0, 1.0).unwrap();
        let field = OceanWaveField::flat(0.0).unwrap();

        let omega_n = (RHO * G * area / mass).sqrt();
        let draft = mass / (RHO * area);
        // Equilibrium: keel at y_eq - height/2 sits 'draft' below surface(0).
        // Reference (center) equilibrium height:
        let y_eq = -draft + height / 2.0;

        // Release from rest displaced upward by 0.2 m — ~50 layers of swing, so
        // the volume staircase closely tracks the linear ramp, yet still partly
        // immersed throughout (draft ≈ 0.49 m).
        let amp = 0.2;
        let mut state = BodyState::at_rest(Vector3::new(0.0, y_eq + amp, 0.0));
        let dt = 5.0e-4;
        let steps = 40_000; // 20 s, omega*dt ≈ 2e-3 (accurate, stable)
                            // Track zero-up crossings of (y - y_eq) to measure the period.
        let mut prev = state.position.y - y_eq;
        let mut t = 0.0;
        let mut first_cross = None;
        let mut last_cross = None;
        let mut crossings = 0u32;
        for _ in 0..steps {
            state = sim.step(&state, &field, t, dt).unwrap();
            t += dt;
            let cur = state.position.y - y_eq;
            if prev <= 0.0 && cur > 0.0 {
                // upward zero crossing
                // linear interp for sub-step crossing time
                let frac = -prev / (cur - prev);
                let tc = t - dt + frac * dt;
                if first_cross.is_none() {
                    first_cross = Some(tc);
                }
                last_cross = Some(tc);
                crossings += 1;
            }
            prev = cur;
        }
        let (f0, fl) = (first_cross.unwrap(), last_cross.unwrap());
        assert!(crossings >= 3, "too few oscillations: {crossings}");
        let measured_period = (fl - f0) / (crossings - 1) as f64;
        let measured_omega = std::f64::consts::TAU / measured_period;
        // The buoyancy restoring is EXACTLY linear here (V_sub = A*s while the
        // body stays partially immersed), so the motion is true SHM; the only
        // departures from omega_n are the ~4 mm volume-quadrature staircase and
        // the symplectic-Euler period error (~(omega*dt)^2/24, negligible).
        // Measured omega lands within ~1.5% of the analytic value.
        let rel = (measured_omega - omega_n).abs() / omega_n;
        assert!(
            rel < 0.015,
            "heave omega {measured_omega} vs analytic {omega_n} (rel {rel})"
        );
    }

    // ---- BENCHMARK-PIN (4): zero net force at equilibrium ----

    #[test]
    fn zero_net_force_at_equilibrium() {
        // A prismatic sample body floating at its equilibrium draft in a FLAT
        // sea: buoyancy exactly cancels weight, net force ~ 0, and a step does
        // not move it (within the integrator's O(dt) bias, which is zero here
        // because velocity starts at zero and net force is ~0).
        let area = 1.0_f64;
        let side_xz = area.sqrt();
        let height = 4.0;
        let mass = 500.0;
        let nz = 200usize; // fine in y so the draft resolves well
        let per = (area * height) / nz as f64;
        let mut samples = Vec::new();
        for i in 0..nz {
            let y = -height / 2.0 + height * (i as f64 + 0.5) / nz as f64;
            let h = side_xz / 2.0;
            for (sx, sz) in [(-h, -h), (h, -h), (-h, h), (h, h)] {
                samples.push((Vector3::new(sx, y, sz), per / 4.0));
            }
        }
        let body = SampleBody::new(&samples, mass).unwrap();
        let sim = BuoyancySim::new(body, RHO, G, Drag::none(), 0.0, 1.0).unwrap();
        let field = OceanWaveField::flat(0.0).unwrap();
        let draft = mass / (RHO * area);
        let y_eq = -draft + height / 2.0;
        let state = BodyState::at_rest(Vector3::new(0.0, y_eq, 0.0));
        let loads = sim.loads(&state, &field, 0.0);
        // The y-quadrature resolves the draft to one layer thickness; net force
        // is within one layer's buoyancy of zero.
        let layer_force = RHO * G * (area * height / nz as f64);
        assert!(
            loads.force.y.abs() < layer_force,
            "net force {} exceeds one-layer tolerance {}",
            loads.force.y,
            layer_force
        );
        // Horizontal force is exactly zero (symmetry, flat sea).
        assert!(loads.force.x.abs() < 1e-9);
        assert!(loads.force.z.abs() < 1e-9);
    }

    #[test]
    fn drag_opposes_motion_and_is_zero_at_rest() {
        let drag = Drag::new(10.0, 2.0).unwrap();
        assert_eq!(drag.force(&Vector3::zeros()), Vector3::zeros());
        let v = Vector3::new(3.0, 0.0, 0.0);
        let f = drag.force(&v);
        // Opposes motion.
        assert!(f.x < 0.0);
        // Magnitude = c_lin*|v| + c_quad*|v|^2.
        let expected = 10.0 * 3.0 + 2.0 * 9.0;
        assert!((f.norm() - expected).abs() < 1e-9);
    }

    // ---- fail-loud config ----

    #[test]
    fn rejects_bad_config() {
        assert!(SampleBody::new(&[], 1.0).is_err());
        assert!(SampleBody::new(&[(Vector3::zeros(), 1.0)], 0.0).is_err());
        assert!(SampleBody::new(&[(Vector3::zeros(), -1.0)], 1.0).is_err());
        assert!(SampleBody::new(&[(Vector3::new(f64::NAN, 0.0, 0.0), 1.0)], 1.0).is_err());
        assert!(HullBody::new(0.0, 1.0, 1.0).is_err());
        assert!(HullBody::new(1.0, 0.0, 1.0).is_err());
        assert!(HullBody::new(1.0, 1.0, 0.0).is_err());
        assert!(HullBody::new(f64::INFINITY, 1.0, 1.0).is_err());
        assert!(Drag::new(-1.0, 0.0).is_err());
        assert!(Drag::new(0.0, -1.0).is_err());
        let body = cube_body(1.0, 2, 1.0);
        assert!(BuoyancySim::new(body.clone(), 0.0, G, Drag::none(), 0.0, 1.0).is_err());
        assert!(BuoyancySim::new(body.clone(), RHO, 0.0, Drag::none(), 0.0, 1.0).is_err());
        assert!(BuoyancySim::new(body.clone(), RHO, G, Drag::none(), -1.0, 1.0).is_err());
        assert!(BuoyancySim::new(body.clone(), RHO, G, Drag::none(), 0.0, 0.0).is_err());
        let sim = BuoyancySim::new(body, RHO, G, Drag::none(), 0.0, 1.0).unwrap();
        let field = OceanWaveField::flat(0.0).unwrap();
        let st = BodyState::at_rest(Vector3::zeros());
        assert!(sim.step(&st, &field, 0.0, 0.0).is_err());
        assert!(sim.step(&st, &field, 0.0, f64::NAN).is_err());
    }

    #[test]
    fn hull_too_dense_to_float_has_draft_exceeding_height() {
        // A body denser than water: equilibrium draft exceeds its height.
        let hull = HullBody::new(1.0, 1.0, 5000.0).unwrap(); // rho_body=5000 > 1025
        let draft = hull.equilibrium_draft(RHO).unwrap();
        assert!(draft > hull.height(), "draft {draft} should exceed height");
        // submerged_volume caps at total volume (it cannot displace more).
        assert!((hull.submerged_volume(draft) - hull.total_volume()).abs() < 1e-12);
    }
}
