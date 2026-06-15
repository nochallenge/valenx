//! Cam-follower rise motion laws (simple-harmonic and cycloidal).
//!
//! The entry point is [`RiseProfile`], a validated `(lift, beta, law)`
//! triple. Evaluate it at a cam angle with [`RiseProfile::at`] to get a
//! [`FollowerState`] carrying displacement, velocity, acceleration, and
//! jerk. See the crate-level docs for the closed-form expressions and a
//! note on their units (derivatives are with respect to the cam angle).

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::error::CamError;

/// The motion law used to interpolate the follower across the rise.
///
/// Both laws map a normalised position `x = theta / beta` in `[0, 1]`
/// onto a normalised displacement in `[0, 1]`; they differ in how the
/// higher derivatives behave at the ends of the interval.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Serialize, Deserialize)]
pub enum MotionLaw {
    /// Simple-harmonic motion (SHM).
    ///
    /// `s = (L/2)(1 - cos(pi x))`. Velocity is zero at both ends, but the
    /// acceleration is non-zero there, so joining this rise to a dwell
    /// introduces a step (infinite jerk) in acceleration.
    SimpleHarmonic,

    /// Cycloidal motion.
    ///
    /// `s = L(x - sin(2 pi x)/(2 pi))`. Both velocity and acceleration are
    /// zero at each end, making it the smoother law and the standard
    /// choice for high-speed cams; only a finite jerk discontinuity
    /// remains at the boundaries.
    Cycloidal,
}

impl MotionLaw {
    /// A stable, lower-case identifier for this law.
    pub fn name(&self) -> &'static str {
        match self {
            MotionLaw::SimpleHarmonic => "simple-harmonic",
            MotionLaw::Cycloidal => "cycloidal",
        }
    }
}

/// The kinematic state of the follower at one cam angle.
///
/// Every field is expressed with the cam angle `theta` as the independent
/// variable, so [`FollowerState::velocity`] has units of `length / radian`,
/// [`FollowerState::acceleration`] `length / radian^2`, and
/// [`FollowerState::jerk`] `length / radian^3`. Multiply by the relevant
/// power of the cam speed `omega` (rad/s) to obtain a time basis.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FollowerState {
    /// Cam angle at which this state was evaluated, in radians, measured
    /// from the start of the rise.
    pub theta: f64,
    /// Follower displacement `s(theta)`, in the same length unit as the
    /// configured lift. Ranges from `0` at the start to `lift` at the end.
    pub displacement: f64,
    /// First derivative `ds/dtheta`, in `length / radian`.
    pub velocity: f64,
    /// Second derivative `d2s/dtheta2`, in `length / radian^2`.
    pub acceleration: f64,
    /// Third derivative `d3s/dtheta3`, in `length / radian^3`.
    pub jerk: f64,
}

/// A validated cam-follower rise: a `lift`, a rise angle `beta`, and the
/// [`MotionLaw`] that connects them.
///
/// Construct one with [`RiseProfile::new`], which rejects non-finite,
/// negative, or zero inputs, then sample it with [`RiseProfile::at`].
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RiseProfile {
    lift: f64,
    beta: f64,
    law: MotionLaw,
}

impl RiseProfile {
    /// Build a rise profile from a `lift`, a rise angle `beta` (radians),
    /// and a [`MotionLaw`].
    ///
    /// # Errors
    ///
    /// Returns [`CamError::NotFinite`] if either `lift` or `beta` is `NaN`
    /// or infinite, [`CamError::BadParameter`] if `lift` is negative, and
    /// [`CamError::BadParameter`] if `beta` is not strictly positive
    /// (a zero-width rise has no defined kinematics).
    ///
    /// A zero `lift` is accepted and yields an everywhere-zero (flat)
    /// profile, which is a useful degenerate case for composition.
    pub fn new(lift: f64, beta: f64, law: MotionLaw) -> Result<Self, CamError> {
        if !lift.is_finite() {
            return Err(CamError::not_finite("lift", lift));
        }
        if !beta.is_finite() {
            return Err(CamError::not_finite("beta", beta));
        }
        if lift < 0.0 {
            return Err(CamError::negative("lift", lift));
        }
        if beta <= 0.0 {
            return Err(CamError::non_positive("beta", beta));
        }
        Ok(Self { lift, beta, law })
    }

    /// The configured lift (total rise displacement).
    pub fn lift(&self) -> f64 {
        self.lift
    }

    /// The configured rise angle `beta`, in radians.
    pub fn beta(&self) -> f64 {
        self.beta
    }

    /// The configured [`MotionLaw`].
    pub fn law(&self) -> MotionLaw {
        self.law
    }

    /// Evaluate the follower state at cam angle `theta` (radians, measured
    /// from the start of the rise).
    ///
    /// `theta` is **not** clamped to `[0, beta]`: the closed-form
    /// expressions are evaluated as given, so sampling outside the rise
    /// interval extrapolates the underlying trigonometric functions. Pass
    /// `theta` in `[0, beta]` to stay within the physical rise.
    pub fn at(&self, theta: f64) -> FollowerState {
        let x = theta / self.beta;
        let l = self.lift;
        let b = self.beta;
        match self.law {
            MotionLaw::SimpleHarmonic => {
                let (sin, cos) = (PI * x).sin_cos();
                FollowerState {
                    theta,
                    displacement: 0.5 * l * (1.0 - cos),
                    velocity: (PI * l / (2.0 * b)) * sin,
                    acceleration: (PI * PI * l / (2.0 * b * b)) * cos,
                    jerk: -(PI * PI * PI * l / (2.0 * b * b * b)) * sin,
                }
            }
            MotionLaw::Cycloidal => {
                let (sin, cos) = (2.0 * PI * x).sin_cos();
                FollowerState {
                    theta,
                    displacement: l * (x - sin / (2.0 * PI)),
                    velocity: (l / b) * (1.0 - cos),
                    acceleration: (2.0 * PI * l / (b * b)) * sin,
                    jerk: (4.0 * PI * PI * l / (b * b * b)) * cos,
                }
            }
        }
    }

    /// Sample the rise at `n` equally spaced cam angles spanning the full
    /// interval `[0, beta]` inclusive (endpoints included).
    ///
    /// # Errors
    ///
    /// Returns [`CamError::BadParameter`] if `n < 2`, since at least the
    /// two endpoints are required to span the interval.
    pub fn sample(&self, n: usize) -> Result<Vec<FollowerState>, CamError> {
        if n < 2 {
            return Err(CamError::BadParameter {
                name: "n",
                reason: "need at least 2 samples to span the interval",
                value: n as f64,
            });
        }
        let last = (n - 1) as f64;
        Ok((0..n)
            .map(|i| {
                let theta = self.beta * (i as f64) / last;
                self.at(theta)
            })
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point ground-truth comparisons.
    const EPS: f64 = 1e-9;

    /// Assert two floats agree to within [`EPS`]; never compares with `==`.
    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() < EPS,
            "expected {a} ~= {b}, diff {}",
            (a - b).abs()
        );
    }

    // ---- constructor validation -----------------------------------------

    #[test]
    fn rejects_non_positive_beta() {
        for beta in [0.0, -1.0, -0.001] {
            let err = RiseProfile::new(10.0, beta, MotionLaw::Cycloidal).unwrap_err();
            assert_eq!(err.code(), "camdynamics.bad_parameter");
        }
    }

    #[test]
    fn rejects_negative_lift() {
        let err = RiseProfile::new(-1.0, 1.0, MotionLaw::SimpleHarmonic).unwrap_err();
        assert_eq!(err.code(), "camdynamics.bad_parameter");
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert_eq!(
            RiseProfile::new(f64::NAN, 1.0, MotionLaw::Cycloidal)
                .unwrap_err()
                .code(),
            "camdynamics.not_finite"
        );
        assert_eq!(
            RiseProfile::new(1.0, f64::INFINITY, MotionLaw::Cycloidal)
                .unwrap_err()
                .code(),
            "camdynamics.not_finite"
        );
    }

    #[test]
    fn accepts_zero_lift_as_flat_profile() {
        let p = RiseProfile::new(0.0, 1.5, MotionLaw::Cycloidal).unwrap();
        for &th in &[0.0, 0.3, 0.75, 1.2, 1.5] {
            let s = p.at(th);
            close(s.displacement, 0.0);
            close(s.velocity, 0.0);
            close(s.acceleration, 0.0);
            close(s.jerk, 0.0);
        }
    }

    #[test]
    fn accessors_round_trip() {
        let p = RiseProfile::new(12.0, 2.0, MotionLaw::SimpleHarmonic).unwrap();
        close(p.lift(), 12.0);
        close(p.beta(), 2.0);
        assert_eq!(p.law(), MotionLaw::SimpleHarmonic);
        assert_eq!(p.law().name(), "simple-harmonic");
    }

    // ---- SHM ground truth -----------------------------------------------

    #[test]
    fn shm_endpoint_displacements() {
        // s(0) = 0 and s(beta) = lift, for any lift/beta.
        let lift = 25.0;
        let beta = 1.3;
        let p = RiseProfile::new(lift, beta, MotionLaw::SimpleHarmonic).unwrap();
        close(p.at(0.0).displacement, 0.0);
        close(p.at(beta).displacement, lift);
    }

    #[test]
    fn shm_velocity_zero_at_ends() {
        let beta = 0.9;
        let p = RiseProfile::new(40.0, beta, MotionLaw::SimpleHarmonic).unwrap();
        close(p.at(0.0).velocity, 0.0);
        close(p.at(beta).velocity, 0.0);
    }

    #[test]
    fn shm_midpoint_is_half_lift_and_peak_velocity() {
        // At x = 1/2: s = L/2, v = pi L / (2 beta) (the maximum), a = 0.
        let lift = 30.0;
        let beta = 1.1;
        let p = RiseProfile::new(lift, beta, MotionLaw::SimpleHarmonic).unwrap();
        let mid = p.at(beta / 2.0);
        close(mid.displacement, lift / 2.0);
        close(mid.velocity, PI * lift / (2.0 * beta));
        close(mid.acceleration, 0.0);
    }

    #[test]
    fn shm_acceleration_at_ends_is_nonzero_and_signed() {
        // a(0) = +pi^2 L / (2 beta^2), a(beta) = -pi^2 L / (2 beta^2).
        let lift = 18.0;
        let beta = 0.8;
        let p = RiseProfile::new(lift, beta, MotionLaw::SimpleHarmonic).unwrap();
        let amax = PI * PI * lift / (2.0 * beta * beta);
        close(p.at(0.0).acceleration, amax);
        close(p.at(beta).acceleration, -amax);
        // Non-zero end acceleration is the defining SHM "rough" feature.
        assert!(p.at(0.0).acceleration.abs() > 1.0);
    }

    #[test]
    fn shm_jerk_closed_form() {
        // j(x) = -(pi^3 L / (2 beta^3)) sin(pi x); peak magnitude at x=1/2.
        let lift = 10.0;
        let beta = 1.0;
        let p = RiseProfile::new(lift, beta, MotionLaw::SimpleHarmonic).unwrap();
        // At the ends sin(pi x) = 0 so jerk vanishes there.
        close(p.at(0.0).jerk, 0.0);
        close(p.at(beta).jerk, 0.0);
        // At x = 1/2, sin = 1 so jerk = -(pi^3 L)/(2 beta^3).
        let jmid = -(PI * PI * PI * lift) / (2.0 * beta * beta * beta);
        close(p.at(beta / 2.0).jerk, jmid);
    }

    // ---- cycloidal ground truth -----------------------------------------

    #[test]
    fn cycloidal_endpoint_displacements() {
        let lift = 22.0;
        let beta = 1.7;
        let p = RiseProfile::new(lift, beta, MotionLaw::Cycloidal).unwrap();
        close(p.at(0.0).displacement, 0.0);
        close(p.at(beta).displacement, lift);
    }

    #[test]
    fn cycloidal_velocity_zero_at_ends() {
        let beta = 1.25;
        let p = RiseProfile::new(50.0, beta, MotionLaw::Cycloidal).unwrap();
        close(p.at(0.0).velocity, 0.0);
        close(p.at(beta).velocity, 0.0);
    }

    #[test]
    fn cycloidal_acceleration_zero_at_ends() {
        // The defining smoothness property: a(0) = a(beta) = 0.
        let beta = 1.05;
        let p = RiseProfile::new(33.0, beta, MotionLaw::Cycloidal).unwrap();
        close(p.at(0.0).acceleration, 0.0);
        close(p.at(beta).acceleration, 0.0);
    }

    #[test]
    fn cycloidal_midpoint_displacement_and_velocity() {
        // At x = 1/2: sin(2 pi x) = sin(pi) = 0, so s = L/2 exactly.
        // v = (L/beta)(1 - cos(pi)) = 2 L / beta (the peak velocity).
        let lift = 16.0;
        let beta = 0.95;
        let p = RiseProfile::new(lift, beta, MotionLaw::Cycloidal).unwrap();
        let mid = p.at(beta / 2.0);
        close(mid.displacement, lift / 2.0);
        close(mid.velocity, 2.0 * lift / beta);
    }

    #[test]
    fn cycloidal_peak_acceleration_at_quarter() {
        // a(x) = (2 pi L / beta^2) sin(2 pi x); at x = 1/4, sin(pi/2)=1.
        let lift = 20.0;
        let beta = 1.0;
        let p = RiseProfile::new(lift, beta, MotionLaw::Cycloidal).unwrap();
        let amax = 2.0 * PI * lift / (beta * beta);
        close(p.at(beta / 4.0).acceleration, amax);
        close(p.at(3.0 * beta / 4.0).acceleration, -amax);
    }

    #[test]
    fn cycloidal_jerk_closed_form_at_ends() {
        // j(x) = (4 pi^2 L / beta^3) cos(2 pi x); at the ends cos = 1, so
        // jerk is at its (finite, equal) positive peak there.
        let lift = 14.0;
        let beta = 1.0;
        let p = RiseProfile::new(lift, beta, MotionLaw::Cycloidal).unwrap();
        let jmax = 4.0 * PI * PI * lift / (beta * beta * beta);
        close(p.at(0.0).jerk, jmax);
        close(p.at(beta).jerk, jmax);
        // Mid-rise cos(pi) = -1 -> negative peak.
        close(p.at(beta / 2.0).jerk, -jmax);
    }

    // ---- scaling laws ----------------------------------------------------

    #[test]
    fn acceleration_scales_linearly_with_lift() {
        // Doubling the lift doubles the acceleration everywhere (both laws).
        let beta = 1.2;
        for law in [MotionLaw::SimpleHarmonic, MotionLaw::Cycloidal] {
            let p1 = RiseProfile::new(10.0, beta, law).unwrap();
            let p2 = RiseProfile::new(20.0, beta, law).unwrap();
            for &th in &[0.1, 0.4, 0.7, 1.0] {
                let theta = th * beta;
                close(p2.at(theta).acceleration, 2.0 * p1.at(theta).acceleration);
            }
        }
    }

    #[test]
    fn acceleration_scales_as_inverse_beta_squared() {
        // Halving beta (same normalised position x) multiplies peak
        // acceleration by 4 = 1/(1/2)^2 for both laws.
        let lift = 12.0;
        // SHM peak |a| is at x=0; cycloidal peak |a| is at x=1/4.
        let cases = [
            (MotionLaw::SimpleHarmonic, 0.0_f64),
            (MotionLaw::Cycloidal, 0.25_f64),
        ];
        for (law, x) in cases {
            let big = RiseProfile::new(lift, 2.0, law).unwrap();
            let small = RiseProfile::new(lift, 1.0, law).unwrap();
            let a_big = big.at(2.0 * x).acceleration.abs();
            let a_small = small.at(1.0 * x).acceleration.abs();
            close(a_small, 4.0 * a_big);
        }
    }

    #[test]
    fn cycloidal_is_smoother_than_shm_at_ends() {
        // Concrete statement of "cycloidal is smoother": at the rise ends
        // SHM has finite non-zero acceleration while cycloidal is zero.
        let lift = 20.0;
        let beta = 1.0;
        let shm = RiseProfile::new(lift, beta, MotionLaw::SimpleHarmonic).unwrap();
        let cyc = RiseProfile::new(lift, beta, MotionLaw::Cycloidal).unwrap();
        for &theta in &[0.0, beta] {
            assert!(shm.at(theta).acceleration.abs() > 1.0);
            close(cyc.at(theta).acceleration, 0.0);
        }
    }

    // ---- numerical derivative cross-checks -------------------------------

    #[test]
    fn velocity_is_numerical_derivative_of_displacement() {
        // Finite-difference ds/dtheta should match the analytic velocity.
        let lift = 17.0;
        let beta = 1.4;
        let h = 1e-6;
        for law in [MotionLaw::SimpleHarmonic, MotionLaw::Cycloidal] {
            let p = RiseProfile::new(lift, beta, law).unwrap();
            for &x in &[0.15, 0.35, 0.55, 0.85] {
                let theta = x * beta;
                let num = (p.at(theta + h).displacement - p.at(theta - h).displacement) / (2.0 * h);
                let ana = p.at(theta).velocity;
                assert!(
                    (num - ana).abs() < 1e-4,
                    "{} d/dtheta mismatch at x={x}: num {num} vs ana {ana}",
                    law.name()
                );
            }
        }
    }

    #[test]
    fn acceleration_is_numerical_derivative_of_velocity() {
        let lift = 9.0;
        let beta = 1.1;
        let h = 1e-5;
        for law in [MotionLaw::SimpleHarmonic, MotionLaw::Cycloidal] {
            let p = RiseProfile::new(lift, beta, law).unwrap();
            for &x in &[0.2, 0.5, 0.8] {
                let theta = x * beta;
                let num = (p.at(theta + h).velocity - p.at(theta - h).velocity) / (2.0 * h);
                let ana = p.at(theta).acceleration;
                assert!(
                    (num - ana).abs() < 1e-3,
                    "{} a mismatch at x={x}: num {num} vs ana {ana}",
                    law.name()
                );
            }
        }
    }

    // ---- sampling --------------------------------------------------------

    #[test]
    fn sample_spans_endpoints_inclusive() {
        let lift = 8.0;
        let beta = 2.0;
        let p = RiseProfile::new(lift, beta, MotionLaw::Cycloidal).unwrap();
        let s = p.sample(11).unwrap();
        assert_eq!(s.len(), 11);
        close(s.first().unwrap().theta, 0.0);
        close(s.last().unwrap().theta, beta);
        close(s.first().unwrap().displacement, 0.0);
        close(s.last().unwrap().displacement, lift);
    }

    #[test]
    fn sample_rejects_too_few_points() {
        let p = RiseProfile::new(1.0, 1.0, MotionLaw::Cycloidal).unwrap();
        assert!(p.sample(1).is_err());
        assert!(p.sample(0).is_err());
        assert!(p.sample(2).is_ok());
    }

    #[test]
    fn cycloidal_displacement_monotonic_nondecreasing_over_rise() {
        // Physical sanity: the rise never reverses (v >= 0 throughout).
        let p = RiseProfile::new(10.0, 1.0, MotionLaw::Cycloidal).unwrap();
        let s = p.sample(101).unwrap();
        for w in s.windows(2) {
            assert!(
                w[1].displacement >= w[0].displacement - EPS,
                "non-monotonic: {} then {}",
                w[0].displacement,
                w[1].displacement
            );
            assert!(w[0].velocity >= -EPS, "negative velocity {}", w[0].velocity);
        }
    }
}
