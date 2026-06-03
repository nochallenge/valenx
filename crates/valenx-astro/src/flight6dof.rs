//! Coupled 6-DOF flight demonstrator: translational motion + rigid-body
//! attitude under a closed-loop pointing controller.
//!
//! This closes the loop the earlier bricks left open: the [`rigidbody`]
//! rotational core is now *driven by a control law* and *coupled to
//! translation* — thrust acts along the vehicle's body axis, and a
//! proportional-derivative (PD) controller torques the body to point
//! that axis where commanded. So the thrust direction is an *outcome of
//! the control loop*, not a free input.
//!
//! [`rigidbody`]: crate::rigidbody
//!
//! What makes this honest to validate: PD attitude control of a rigid
//! body is a textbook asymptotically-stable system, so the oracle is
//! control-theoretic behaviour — the closed loop drives the pointing
//! error to zero, settles without limit-cycling, and holds attitude
//! against a steady disturbance with bounded error. The tests pin that.
//!
//! Scope — a **v1 demonstrator, not flight-certified GNC**: rigid (no
//! flex/slosh), a single PD pointing law (no guidance/navigation filter,
//! no actuator limits or lag, roll about the thrust axis left free),
//! constant gravity, thrust along body +x. A production flight GNC stack
//! (guidance + estimation + control with actuator models, validated
//! against flight data) is the research-grade item in RFC 0010.

use nalgebra::{UnitQuaternion, Vector3};
use serde::{Deserialize, Serialize};

use crate::error::AstroError;
use crate::rigidbody::{self, AttitudeState, Inertia};
use crate::sim::check_step_count;

/// Full rigid-body state: translation + rotation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct State6dof {
    /// Inertial position (m).
    pub position: Vector3<f64>,
    /// Inertial velocity (m/s).
    pub velocity: Vector3<f64>,
    /// Body→inertial orientation.
    pub attitude: UnitQuaternion<f64>,
    /// Body-frame angular velocity (rad/s).
    pub angular_velocity: Vector3<f64>,
}

impl State6dof {
    /// At rest at the origin with a given orientation.
    pub fn at_rest(attitude: UnitQuaternion<f64>) -> Self {
        Self {
            position: Vector3::zeros(),
            velocity: Vector3::zeros(),
            attitude,
            angular_velocity: Vector3::zeros(),
        }
    }

    /// The current inertial direction the body +x ("nose") axis points.
    pub fn pointing(&self) -> Vector3<f64> {
        self.attitude.transform_vector(&Vector3::x())
    }
}

/// PD pointing-controller gains.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ControlGains {
    /// Proportional gain on the pointing error (N·m per unit sin-error).
    pub kp: f64,
    /// Derivative gain on the body angular rate (N·m·s).
    pub kd: f64,
}

impl ControlGains {
    /// Critically-damped-ish gains for a body of the given mean inertia.
    pub fn for_inertia(mean_inertia: f64, bandwidth: f64) -> Self {
        let kp = mean_inertia * bandwidth * bandwidth;
        let kd = 2.0 * mean_inertia * bandwidth;
        Self { kp, kd }
    }
}

/// A commanded pointing target whose norm is below this (or non-finite)
/// carries no usable direction: `target.normalize()` would yield a NaN
/// unit vector. Such a command is treated as "no valid command, hold".
const MIN_TARGET_NORM: f64 = 1e-12;

/// Angle (rad) between the body +x axis and a target inertial direction.
///
/// A near-zero or non-finite `target` carries no usable direction; rather
/// than normalize it into a NaN unit vector (which would poison the
/// reported error and any [`pointing_control_torque`] driven from it),
/// this returns `0.0` — "no valid command, nothing to correct".
pub fn pointing_error(attitude: &UnitQuaternion<f64>, target: Vector3<f64>) -> f64 {
    let norm = target.norm();
    if !norm.is_finite() || norm < MIN_TARGET_NORM {
        return 0.0;
    }
    let nose = attitude.transform_vector(&Vector3::x());
    let t = target / norm;
    nose.dot(&t).clamp(-1.0, 1.0).acos()
}

/// Body-frame control torque (N·m) for a PD law that drives the body +x
/// axis toward the inertial `target` direction.
///
/// `τ = Kp · (n̂ × t̂)|_body − Kd · ω`. The cross product is the rotation
/// axis (scaled by the sine of the error); damping on the body rate
/// suppresses overshoot. At a 180° error the cross product vanishes (an
/// unstable equilibrium) — start away from anti-parallel.
///
/// A near-zero or non-finite `target` carries no usable direction;
/// normalizing it would yield a NaN axis and a NaN torque that silently
/// corrupts the attitude integration. In that degenerate case the
/// controller issues **zero torque** ("no valid command, hold") rather
/// than acting on a NaN command.
pub fn pointing_control_torque(
    attitude: &UnitQuaternion<f64>,
    angular_velocity: Vector3<f64>,
    target: Vector3<f64>,
    gains: ControlGains,
) -> Vector3<f64> {
    let norm = target.norm();
    if !norm.is_finite() || norm < MIN_TARGET_NORM {
        return Vector3::zeros();
    }
    let nose = attitude.transform_vector(&Vector3::x());
    let axis_inertial = nose.cross(&(target / norm));
    let axis_body = attitude.inverse_transform_vector(&axis_inertial);
    gains.kp * axis_body - gains.kd * angular_velocity
}

/// One coupled step (s). Thrust acts along the body +x axis; the
/// controller torque (plus any `disturbance` body torque) drives the
/// attitude; `gravity` is a constant inertial acceleration (m/s²).
///
/// # Preconditions
///
/// `mass` must be finite and positive (thrust divides by it) and
/// `inertia` must satisfy [`Inertia::validate`]. [`simulate_pointing`]
/// validates both before stepping; call it for a checked entry point.
#[allow(clippy::too_many_arguments)]
pub fn step(
    state: &State6dof,
    inertia: &Inertia,
    thrust: f64,
    mass: f64,
    gravity: Vector3<f64>,
    target: Vector3<f64>,
    gains: ControlGains,
    disturbance: Vector3<f64>,
    dt: f64,
) -> State6dof {
    // Rotational update under control + disturbance torque.
    let torque = pointing_control_torque(&state.attitude, state.angular_velocity, target, gains)
        + disturbance;
    // A single step is trivially within the step-count cap, so use the
    // unchecked core to keep `step` infallible.
    let att = rigidbody::propagate_unchecked(
        &AttitudeState {
            angular_velocity: state.angular_velocity,
            attitude: state.attitude,
        },
        inertia,
        torque,
        dt,
        1,
    );

    // Translational update: thrust along the (step-start) body axis plus
    // gravity, integrated as constant acceleration over the step.
    let thrust_dir = state.pointing();
    let accel = gravity + thrust_dir * (thrust / mass);
    let position = state.position + state.velocity * dt + 0.5 * accel * dt * dt;
    let velocity = state.velocity + accel * dt;

    State6dof {
        position,
        velocity,
        attitude: att.attitude,
        angular_velocity: att.angular_velocity,
    }
}

/// Run the closed-loop demonstrator for `steps` steps holding a fixed
/// pointing `target`. Returns the final state and the **settled** peak
/// pointing error (rad) — the worst error over the last 20 % of the run,
/// so the initial slew is excluded but any steady offset or limit cycle
/// is captured.
///
/// # Errors
///
/// Returns [`AstroError::OutOfRange`] if `steps` exceeds
/// [`crate::sim::MAX_SIM_STEPS`], or [`AstroError::InvalidMass`] /
/// [`AstroError::InvalidPropulsion`] if `mass` is non-positive /
/// non-finite or the inertia is non-physical — either would otherwise
/// drive the translational or rotational update to Inf/NaN.
#[allow(clippy::too_many_arguments)]
pub fn simulate_pointing(
    initial: &State6dof,
    inertia: &Inertia,
    thrust: f64,
    mass: f64,
    gravity: Vector3<f64>,
    target: Vector3<f64>,
    gains: ControlGains,
    disturbance: Vector3<f64>,
    dt: f64,
    steps: u64,
) -> Result<(State6dof, f64), AstroError> {
    check_step_count(steps)?;
    inertia.validate()?;
    if !mass.is_finite() || mass <= 0.0 {
        return Err(AstroError::InvalidMass {
            index: usize::MAX,
            field: "mass",
            value: mass,
        });
    }
    let mut s = *initial;
    let settle_after = steps * 4 / 5; // last 20 % of the run
    let mut settled_peak = 0.0_f64;
    for i in 0..steps {
        s = step(
            &s,
            inertia,
            thrust,
            mass,
            gravity,
            target,
            gains,
            disturbance,
            dt,
        );
        if i >= settle_after {
            settled_peak = settled_peak.max(pointing_error(&s.attitude, target));
        }
    }
    Ok((s, settled_peak))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn body() -> Inertia {
        Inertia::new(1.0, 1.0, 1.0)
    }

    #[test]
    fn controller_slews_to_commanded_pointing() {
        // Nose starts at +x; command it to +y. The loop should slew over
        // and settle with a small steady pointing error.
        let gains = ControlGains::for_inertia(1.0, 3.0);
        let s0 = State6dof::at_rest(UnitQuaternion::identity());
        let target = Vector3::new(0.0, 1.0, 0.0);
        let (sf, settled_peak) = simulate_pointing(
            &s0,
            &body(),
            0.0,
            1_000.0,
            Vector3::zeros(),
            target,
            gains,
            Vector3::zeros(),
            1e-3,
            8_000, // 8 s
        )
        .expect("valid steps");
        // Settled to within ~1° and no residual tumble.
        assert!(
            settled_peak.to_degrees() < 1.0,
            "settled error {} deg",
            settled_peak.to_degrees()
        );
        assert!(
            sf.angular_velocity.norm() < 1e-2,
            "residual ω {}",
            sf.angular_velocity.norm()
        );
        // Final pointing really is ~+y.
        assert!(pointing_error(&sf.attitude, target).to_degrees() < 1.0);
    }

    #[test]
    fn thrust_follows_the_commanded_pointing() {
        // With thrust on and no gravity, once the nose tracks +y the
        // vehicle accelerates predominantly along +y.
        let gains = ControlGains::for_inertia(1.0, 3.0);
        let s0 = State6dof::at_rest(UnitQuaternion::identity());
        let target = Vector3::new(0.0, 1.0, 0.0);
        let (sf, _) = simulate_pointing(
            &s0,
            &body(),
            5_000.0, // thrust (N)
            1_000.0, // mass (kg)
            Vector3::zeros(),
            target,
            gains,
            Vector3::zeros(),
            1e-3,
            8_000,
        )
        .expect("valid steps");
        // Velocity is mostly +y (thrust steered there), with small x/z.
        assert!(sf.velocity.y > 0.0);
        let off_axis = (sf.velocity.x.powi(2) + sf.velocity.z.powi(2)).sqrt();
        assert!(
            off_axis / sf.velocity.y < 0.2,
            "off-axis frac {}",
            off_axis / sf.velocity.y
        );
    }

    #[test]
    fn rejects_a_steady_disturbance_with_bounded_error() {
        // A constant disturbance torque produces a bounded steady-state
        // pointing offset, not a divergence.
        let gains = ControlGains::for_inertia(1.0, 4.0);
        let s0 = State6dof::at_rest(UnitQuaternion::identity());
        let target = Vector3::new(1.0, 0.0, 0.0); // hold the initial pointing
        let disturbance = Vector3::new(0.0, 0.0, 2.0); // N·m about body z
        let (sf, settled_peak) = simulate_pointing(
            &s0,
            &body(),
            0.0,
            1_000.0,
            Vector3::zeros(),
            target,
            gains,
            disturbance,
            1e-3,
            10_000,
        )
        .expect("valid steps");
        // Bounded, small steady-state error (does not run away).
        assert!(
            settled_peak.to_degrees() < 10.0,
            "steady error {} deg",
            settled_peak.to_degrees()
        );
        assert!(sf.attitude.quaternion().norm().is_finite());
    }

    #[test]
    fn simulate_pointing_rejects_zero_mass_and_bad_inertia() {
        // thrust/mass with mass = 0 -> Inf; a zero inertia -> NaN. Both
        // must be rejected instead of integrating to a non-finite state.
        let gains = ControlGains::for_inertia(1.0, 3.0);
        let s0 = State6dof::at_rest(UnitQuaternion::identity());
        let zero_mass = simulate_pointing(
            &s0,
            &body(),
            5_000.0,
            0.0, // zero mass
            Vector3::zeros(),
            Vector3::x(),
            gains,
            Vector3::zeros(),
            1e-3,
            100,
        );
        assert!(
            matches!(zero_mass, Err(AstroError::InvalidMass { .. })),
            "zero mass must be rejected, got {zero_mass:?}"
        );
        let bad_inertia = simulate_pointing(
            &s0,
            &Inertia::new(0.0, 1.0, 1.0),
            0.0,
            1_000.0,
            Vector3::zeros(),
            Vector3::x(),
            gains,
            Vector3::zeros(),
            1e-3,
            100,
        );
        assert!(
            matches!(bad_inertia, Err(AstroError::InvalidPropulsion { .. })),
            "zero inertia must be rejected, got {bad_inertia:?}"
        );
    }

    #[test]
    fn simulate_pointing_rejects_absurd_step_count() {
        // u64::MAX steps would hang; expect a clean Err.
        let gains = ControlGains::for_inertia(1.0, 3.0);
        let s0 = State6dof::at_rest(UnitQuaternion::identity());
        let r = simulate_pointing(
            &s0,
            &body(),
            0.0,
            1_000.0,
            Vector3::zeros(),
            Vector3::x(),
            gains,
            Vector3::zeros(),
            1e-3,
            u64::MAX,
        );
        assert!(
            matches!(r, Err(AstroError::OutOfRange { what: "steps", .. })),
            "u64::MAX steps must be rejected, got {r:?}"
        );
    }

    // A zero (or non-finite) commanded target normalizes to a NaN unit
    // vector; pre-fix this poisoned the pointing error and the control
    // torque, silently corrupting the attitude integration. The guard now
    // makes a degenerate command a zero-error / zero-torque no-op ("hold").
    #[test]
    fn degenerate_target_is_no_op_not_nan() {
        let att = UnitQuaternion::identity();
        let gains = ControlGains::for_inertia(1.0, 3.0);
        for bad in [
            Vector3::zeros(),
            Vector3::new(f64::NAN, 0.0, 0.0),
            Vector3::new(f64::INFINITY, 0.0, 0.0),
            Vector3::new(1e-20, 0.0, 0.0), // below MIN_TARGET_NORM
        ] {
            let err = pointing_error(&att, bad);
            assert_eq!(err, 0.0, "degenerate target {bad:?} must give 0 error");
            let tau = pointing_control_torque(&att, Vector3::zeros(), bad, gains);
            assert_eq!(
                tau,
                Vector3::zeros(),
                "degenerate target {bad:?} must give zero torque"
            );
        }
    }

    // Valid (non-degenerate) commands are byte-identical to the original
    // `target.normalize()` formulation: the guard only diverts the
    // degenerate case.
    #[test]
    fn valid_target_pointing_unchanged() {
        let att = UnitQuaternion::from_euler_angles(0.3, -0.4, 0.5);
        let gains = ControlGains::for_inertia(2.0, 4.0);
        let omega = Vector3::new(0.1, -0.2, 0.05);
        let targets: [Vector3<f64>; 3] = [
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(3.0, -2.0, 1.0), // un-normalized, still valid
            Vector3::new(-1.0, 0.0, 0.0),
        ];
        for t in targets {
            // Recompute the original (pre-guard) expressions explicitly.
            let nose = att.transform_vector(&Vector3::x());
            let expect_err = nose.dot(&t.normalize()).clamp(-1.0, 1.0).acos();
            let axis_inertial = nose.cross(&t.normalize());
            let axis_body = att.inverse_transform_vector(&axis_inertial);
            let expect_tau = gains.kp * axis_body - gains.kd * omega;

            assert_eq!(pointing_error(&att, t), expect_err, "error for {t:?}");
            assert_eq!(
                pointing_control_torque(&att, omega, t, gains),
                expect_tau,
                "torque for {t:?}"
            );
        }
    }

    #[test]
    fn gravity_only_is_free_fall() {
        // No thrust: the body just falls under gravity (sanity on the
        // translational coupling).
        let gains = ControlGains::for_inertia(1.0, 3.0);
        let s0 = State6dof::at_rest(UnitQuaternion::identity());
        let g = Vector3::new(0.0, 0.0, -9.81);
        let (sf, _) = simulate_pointing(
            &s0,
            &body(),
            0.0,
            1_000.0,
            g,
            Vector3::x(),
            gains,
            Vector3::zeros(),
            1e-3,
            1_000, // 1 s
        )
        .expect("valid steps");
        // v ≈ g·t = -9.81 m/s after 1 s; z ≈ -½g t² = -4.905 m.
        assert!((sf.velocity.z + 9.81).abs() < 0.05, "vz {}", sf.velocity.z);
        assert!((sf.position.z + 4.905).abs() < 0.05, "z {}", sf.position.z);
    }
}
