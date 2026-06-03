//! Rigid-body rotational dynamics: the attitude core of a 6-DOF model.
//!
//! Integrates Euler's equations for a rigid body about its principal
//! axes, carrying the orientation as a unit quaternion. This is the
//! rotational half of a full 6-DOF flight model (the translational half
//! is the existing trajectory sim); coupling the two with thrust-vector
//! control, aerodynamic moments and a control law is the next step.
//!
//! What makes this brick *honest* is that rotational dynamics has exact
//! conservation laws to validate against: with no applied torque the
//! rotational kinetic energy and the inertial angular-momentum vector
//! are conserved, and an axisymmetric body precesses at the closed-form
//! body-cone rate. The tests pin all three.
//!
//! Scope: principal-axis (diagonal) inertia, ideal rigid body, torques
//! supplied externally (constant per propagation call). No structural
//! flexibility, fuel slosh, or a control loop yet.

use nalgebra::{Quaternion, UnitQuaternion, Vector3, Vector4};
use serde::{Deserialize, Serialize};

use crate::error::AstroError;
use crate::sim::check_step_count;

/// Principal moments of inertia (kg·m²) about the body x/y/z axes.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Inertia {
    /// Moment about the body x-axis.
    pub ix: f64,
    /// Moment about the body y-axis.
    pub iy: f64,
    /// Moment about the body z-axis.
    pub iz: f64,
}

impl Inertia {
    /// A new principal-axis inertia tensor.
    pub fn new(ix: f64, iy: f64, iz: f64) -> Self {
        Self { ix, iy, iz }
    }

    /// Validate the inertia tensor: every principal moment must be
    /// finite and strictly positive.
    ///
    /// # Errors
    ///
    /// Returns [`AstroError::InvalidPropulsion`] (field `ix`/`iy`/`iz`)
    /// for a zero, negative, or non-finite moment — Euler's equations
    /// divide by each moment, so such a value would produce a NaN/Inf
    /// angular acceleration.
    pub fn validate(&self) -> Result<(), AstroError> {
        for (field, value) in [("ix", self.ix), ("iy", self.iy), ("iz", self.iz)] {
            if !value.is_finite() || value <= 0.0 {
                return Err(AstroError::InvalidPropulsion {
                    index: usize::MAX,
                    field,
                    value,
                });
            }
        }
        Ok(())
    }
}

/// Rigid-body rotational state: body-frame angular velocity (rad/s) and
/// the orientation that rotates body-frame vectors into the inertial
/// frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AttitudeState {
    /// Angular velocity in the body frame (rad/s).
    pub angular_velocity: Vector3<f64>,
    /// Body→inertial orientation.
    pub attitude: UnitQuaternion<f64>,
}

impl AttitudeState {
    /// State from an angular velocity and an identity orientation.
    pub fn spinning(angular_velocity: Vector3<f64>) -> Self {
        Self {
            angular_velocity,
            attitude: UnitQuaternion::identity(),
        }
    }
}

/// Body-frame angular momentum `H = I·ω` (kg·m²/s).
pub fn angular_momentum_body(inertia: &Inertia, omega: Vector3<f64>) -> Vector3<f64> {
    Vector3::new(
        inertia.ix * omega.x,
        inertia.iy * omega.y,
        inertia.iz * omega.z,
    )
}

/// Inertial-frame angular momentum (kg·m²/s): the body momentum rotated
/// by the current attitude.
pub fn angular_momentum_inertial(inertia: &Inertia, state: &AttitudeState) -> Vector3<f64> {
    state
        .attitude
        .transform_vector(&angular_momentum_body(inertia, state.angular_velocity))
}

/// Rotational kinetic energy `½ ωᵀ I ω` (J).
pub fn kinetic_energy(inertia: &Inertia, omega: Vector3<f64>) -> f64 {
    0.5 * (inertia.ix * omega.x * omega.x
        + inertia.iy * omega.y * omega.y
        + inertia.iz * omega.z * omega.z)
}

/// Euler's equations: body-frame angular acceleration under an external
/// body-frame `torque` (N·m).
///
/// `I·ω̇ = (I·ω)×ω + M`, i.e. componentwise
/// `ω̇ₓ = [(Iy−Iz)ωyωz + Mₓ]/Iₓ`, and cyclic.
///
/// # Preconditions
///
/// Each principal moment must be finite and positive (the components are
/// divided by `ix`/`iy`/`iz`); pass an inertia that satisfies
/// [`Inertia::validate`]. The propagators validate the inertia before
/// calling this in their integration loop.
pub fn angular_acceleration(
    inertia: &Inertia,
    omega: Vector3<f64>,
    torque: Vector3<f64>,
) -> Vector3<f64> {
    let (ix, iy, iz) = (inertia.ix, inertia.iy, inertia.iz);
    Vector3::new(
        ((iy - iz) * omega.y * omega.z + torque.x) / ix,
        ((iz - ix) * omega.z * omega.x + torque.y) / iy,
        ((ix - iy) * omega.x * omega.y + torque.z) / iz,
    )
}

/// Quaternion kinematic derivative `q̇ = ½ q ⊗ (0, ω)`, returned as raw
/// coordinates `[i, j, k, w]`.
fn quat_derivative(q: &Quaternion<f64>, omega: Vector3<f64>) -> Vector4<f64> {
    let omega_pure = Quaternion::from_parts(0.0, omega);
    (q * omega_pure).coords * 0.5
}

/// Propagate the rotational state forward by `steps` RK4 steps of size
/// `dt` (s) under a constant body-frame `torque` (N·m). The attitude
/// quaternion is renormalised after each step.
///
/// # Errors
///
/// Returns [`AstroError::OutOfRange`] if `steps` exceeds
/// [`crate::sim::MAX_SIM_STEPS`], or [`AstroError::InvalidPropulsion`]
/// if the inertia is non-physical (see [`Inertia::validate`]) — a zero
/// or negative moment would otherwise integrate to NaN.
pub fn propagate(
    state: &AttitudeState,
    inertia: &Inertia,
    torque: Vector3<f64>,
    dt: f64,
    steps: u64,
) -> Result<AttitudeState, AstroError> {
    check_step_count(steps)?;
    inertia.validate()?;
    Ok(propagate_unchecked(state, inertia, torque, dt, steps))
}

/// Step core without the step-count cap. Internal use only, for callers
/// that pass a known-bounded `steps` (e.g. the per-step 6-DOF loop).
pub(crate) fn propagate_unchecked(
    state: &AttitudeState,
    inertia: &Inertia,
    torque: Vector3<f64>,
    dt: f64,
    steps: u64,
) -> AttitudeState {
    let mut omega = state.angular_velocity;
    let mut q = *state.attitude.quaternion();

    // Derivative of the combined (ω, q-coords) state.
    let deriv = |omega: Vector3<f64>, qc: Vector4<f64>| -> (Vector3<f64>, Vector4<f64>) {
        let q = Quaternion::from(qc);
        (
            angular_acceleration(inertia, omega, torque),
            quat_derivative(&q, omega),
        )
    };

    for _ in 0..steps {
        let qc = q.coords;
        let (k1o, k1q) = deriv(omega, qc);
        let (k2o, k2q) = deriv(omega + 0.5 * dt * k1o, qc + 0.5 * dt * k1q);
        let (k3o, k3q) = deriv(omega + 0.5 * dt * k2o, qc + 0.5 * dt * k2q);
        let (k4o, k4q) = deriv(omega + dt * k3o, qc + dt * k3q);

        omega += dt / 6.0 * (k1o + 2.0 * k2o + 2.0 * k3o + k4o);
        let new_qc = qc + dt / 6.0 * (k1q + 2.0 * k2q + 2.0 * k3q + k4q);
        q = Quaternion::from(new_qc);
    }

    AttitudeState {
        angular_velocity: omega,
        attitude: UnitQuaternion::from_quaternion(q),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rel(a: f64, b: f64) -> f64 {
        (a - b).abs() / b.abs().max(1e-30)
    }

    #[test]
    fn torque_free_conserves_kinetic_energy() {
        // A triaxial body tumbling with no torque keeps its rotational
        // kinetic energy.
        let inertia = Inertia::new(1.0, 2.0, 3.0);
        let s0 = AttitudeState::spinning(Vector3::new(1.0, 0.5, 0.3));
        let e0 = kinetic_energy(&inertia, s0.angular_velocity);
        let s1 = propagate(&s0, &inertia, Vector3::zeros(), 1e-3, 20_000).expect("valid steps"); // 20 s
        let e1 = kinetic_energy(&inertia, s1.angular_velocity);
        assert!(rel(e1, e0) < 1e-6, "energy {e0} -> {e1}");
    }

    #[test]
    fn torque_free_conserves_inertial_angular_momentum() {
        // The inertial angular-momentum *vector* is fixed without torque,
        // even though the body-frame components shuffle as it tumbles.
        let inertia = Inertia::new(1.0, 2.0, 3.0);
        let s0 = AttitudeState::spinning(Vector3::new(0.8, 0.6, 0.2));
        let h0 = angular_momentum_inertial(&inertia, &s0);
        let s1 = propagate(&s0, &inertia, Vector3::zeros(), 1e-3, 20_000).expect("valid steps");
        let h1 = angular_momentum_inertial(&inertia, &s1);
        assert!((h1 - h0).norm() / h0.norm() < 1e-5, "H {h0:?} -> {h1:?}");
        // The body-frame components really did change (genuine tumble).
        assert!((s1.angular_velocity - s0.angular_velocity).norm() > 1e-3);
    }

    #[test]
    fn axisymmetric_body_precesses_at_the_body_cone_rate() {
        // Symmetric body (Ix = Iy = It, Iz = Ia): the transverse angular
        // velocity rotates in the body frame at λ = ωz·(Ia − It)/It.
        let it = 1.0;
        let ia = 2.0;
        let inertia = Inertia::new(it, it, ia);
        let wz = 1.0;
        let s0 = AttitudeState::spinning(Vector3::new(0.1, 0.0, wz));
        let lambda = wz * (ia - it) / it; // 1.0 rad/s

        let t = 1.0;
        let dt = 1e-4;
        let s1 =
            propagate(&s0, &inertia, Vector3::zeros(), dt, (t / dt) as u64).expect("valid steps");

        // ωz is exactly constant for a symmetric body.
        assert!((s1.angular_velocity.z - wz).abs() < 1e-9);
        // Transverse component magnitude is preserved...
        let perp0 = (s0.angular_velocity.x.powi(2) + s0.angular_velocity.y.powi(2)).sqrt();
        let perp1 = (s1.angular_velocity.x.powi(2) + s1.angular_velocity.y.powi(2)).sqrt();
        assert!(rel(perp1, perp0) < 1e-6, "|ω⊥| {perp0} -> {perp1}");
        // ...and it has rotated by ≈ λ·t.
        let ang0 = s0.angular_velocity.y.atan2(s0.angular_velocity.x);
        let ang1 = s1.angular_velocity.y.atan2(s1.angular_velocity.x);
        let swept = (ang1 - ang0).rem_euclid(std::f64::consts::TAU);
        let expected = (lambda * t).rem_euclid(std::f64::consts::TAU);
        assert!(
            (swept - expected).abs() < 0.02,
            "swept {swept} vs {expected}"
        );
    }

    #[test]
    fn constant_torque_spins_up_linearly() {
        // Non-spinning body, constant torque about z: ωz = (Mz/Iz)·t.
        let inertia = Inertia::new(2.0, 2.0, 2.0);
        let s0 = AttitudeState::spinning(Vector3::zeros());
        let torque = Vector3::new(0.0, 0.0, 4.0);
        let t = 3.0;
        let dt = 1e-3;
        let s1 = propagate(&s0, &inertia, torque, dt, (t / dt) as u64).expect("valid steps");
        let expected = 4.0 / 2.0 * t; // 6 rad/s
        assert!(
            (s1.angular_velocity.z - expected).abs() < 1e-6,
            "ωz {}",
            s1.angular_velocity.z
        );
        assert!(s1.angular_velocity.x.abs() < 1e-9 && s1.angular_velocity.y.abs() < 1e-9);
    }

    #[test]
    fn propagate_rejects_nonphysical_inertia_instead_of_nan() {
        // Zero/negative/non-finite inertia divides to NaN in Euler's
        // equations; propagate must reject it up front.
        assert!(Inertia::new(1.0, 2.0, 3.0).validate().is_ok());
        assert!(Inertia::new(0.0, 2.0, 3.0).validate().is_err());
        assert!(Inertia::new(1.0, -2.0, 3.0).validate().is_err());
        assert!(Inertia::new(1.0, 2.0, f64::NAN).validate().is_err());

        let s0 = AttitudeState::spinning(Vector3::new(1.0, 0.5, 0.3));
        let bad = Inertia::new(0.0, 2.0, 3.0);
        let r = propagate(&s0, &bad, Vector3::zeros(), 1e-3, 10);
        assert!(
            matches!(r, Err(AstroError::InvalidPropulsion { .. })),
            "zero inertia must be rejected, got {r:?}"
        );
        // A finite spin really would have gone NaN with ix = 0.
        let naive = angular_acceleration(&bad, s0.angular_velocity, Vector3::zeros());
        assert!(naive.x.is_nan() || naive.x.is_infinite());
    }

    #[test]
    fn propagate_rejects_absurd_step_count() {
        // u64::MAX steps would hang; expect a clean Err.
        let inertia = Inertia::new(1.0, 2.0, 3.0);
        let s0 = AttitudeState::spinning(Vector3::new(1.0, 0.5, 0.3));
        let r = propagate(&s0, &inertia, Vector3::zeros(), 1e-3, u64::MAX);
        assert!(
            matches!(r, Err(AstroError::OutOfRange { what: "steps", .. })),
            "u64::MAX steps must be rejected, got {r:?}"
        );
        assert!(propagate(&s0, &inertia, Vector3::zeros(), 1e-3, 10).is_ok());
    }

    #[test]
    fn spin_about_a_principal_axis_stays_pure() {
        // A pure spin about a principal axis is an equilibrium: the
        // transverse rates stay zero.
        let inertia = Inertia::new(1.0, 2.0, 3.0);
        let s0 = AttitudeState::spinning(Vector3::new(0.0, 0.0, 5.0));
        let s1 = propagate(&s0, &inertia, Vector3::zeros(), 1e-3, 10_000).expect("valid steps");
        assert!(s1.angular_velocity.x.abs() < 1e-9);
        assert!(s1.angular_velocity.y.abs() < 1e-9);
        assert!((s1.angular_velocity.z - 5.0).abs() < 1e-9);
        // The attitude quaternion stays normalised through integration.
        assert!((s1.attitude.quaternion().norm() - 1.0).abs() < 1e-9);
    }
}
