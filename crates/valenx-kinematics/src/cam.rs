//! Cam-follower displacement laws for the *rise* segment.
//!
//! Two textbook motion programs map a cam rotation angle `θ ∈ [0, β]`
//! (where `β` is the rise duration, in radians of cam rotation) to a
//! follower displacement `s ∈ [0, L]`, where `L` is the total lift.
//! Both start at rest at the bottom dwell and finish at the top dwell.
//!
//! ## Simple-harmonic motion (SHM)
//!
//! ```text
//! s(θ)  = (L/2) · [1 − cos(π θ / β)]
//! s'(θ) = (π L) / (2β) · sin(π θ / β)              (displacement / cam-rad)
//! s''(θ)= (π² L) / (2β²) · cos(π θ / β)
//! ```
//!
//! SHM has *zero velocity* at both ends (`sin 0 = sin π = 0`) but a
//! *finite, non-zero acceleration* there (`cos 0 = 1`, `cos π = −1`),
//! so its acceleration is discontinuous against the neighbouring
//! dwells — fine at low speed, harsh at high speed.
//!
//! ## Cycloidal motion
//!
//! ```text
//! s(θ)  = L · [ θ/β − (1/2π) · sin(2π θ / β) ]
//! s'(θ) = (L/β) · [1 − cos(2π θ / β)]              (displacement / cam-rad)
//! s''(θ)= (2π L / β²) · sin(2π θ / β)
//! ```
//!
//! Cycloidal motion has *both* zero velocity *and* zero acceleration
//! at each end (`1 − cos 0 = 0`, `1 − cos 2π = 0`; `sin 0 = sin 2π =
//! 0`), so it blends smoothly into the dwells — the standard choice
//! for high-speed cams.
//!
//! The derivatives returned here are with respect to the cam angle
//! `θ`. Multiply by the (constant) cam angular velocity `ω` to obtain
//! time derivatives: `ds/dt = ω · s'(θ)`, `d²s/dt² = ω² · s''(θ)`.

use serde::{Deserialize, Serialize};

use crate::error::{require_non_negative, require_positive, KinematicsError};

/// Which follower motion program to evaluate.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum CamMotion {
    /// Simple-harmonic rise (zero end velocity, non-zero end accel).
    SimpleHarmonic,
    /// Cycloidal rise (zero end velocity *and* zero end acceleration).
    Cycloidal,
}

/// A cam *rise* segment: lift `L` achieved over a cam rotation of
/// `β` radians, under one of the [`CamMotion`] programs. Validated on
/// construction — `β` must be finite and strictly positive, `L` finite
/// and non-negative.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CamRise {
    /// Motion program for the segment.
    pub motion: CamMotion,
    /// Total lift `L` (displacement units) reached at `θ = β`.
    pub lift: f64,
    /// Rise duration `β` (radians of cam rotation).
    pub beta: f64,
}

/// Follower state at one cam angle: displacement and the first two
/// derivatives with respect to cam angle.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FollowerState {
    /// Displacement `s` (same units as the lift).
    pub displacement: f64,
    /// First derivative `ds/dθ` (displacement per cam-radian).
    pub velocity: f64,
    /// Second derivative `d²s/dθ²` (displacement per cam-radian²).
    pub acceleration: f64,
}

impl CamRise {
    /// Construct and validate a rise segment.
    ///
    /// # Errors
    /// [`KinematicsError::BadParameter`] if `beta` is not finite and
    /// `> 0`, or `lift` is not finite and `>= 0`.
    pub fn new(motion: CamMotion, lift: f64, beta: f64) -> Result<Self, KinematicsError> {
        require_non_negative("lift", lift)?;
        require_positive("beta", beta)?;
        Ok(Self { motion, lift, beta })
    }

    /// Evaluate displacement and the two cam-angle derivatives at cam
    /// angle `theta`.
    ///
    /// `theta` is clamped to `[0, β]`: angles before the rise read the
    /// bottom dwell (`s = 0`) and angles past it read the top dwell
    /// (`s = L`), matching how a rise segment is stitched between
    /// dwells in a full cam program.
    pub fn evaluate(&self, theta: f64) -> FollowerState {
        let t = theta.clamp(0.0, self.beta);
        match self.motion {
            CamMotion::SimpleHarmonic => self.eval_shm(t),
            CamMotion::Cycloidal => self.eval_cycloidal(t),
        }
    }

    /// Convenience: just the displacement `s(θ)`.
    pub fn rise(&self, theta: f64) -> f64 {
        self.evaluate(theta).displacement
    }

    /// Convenience: just the cam-angle velocity `ds/dθ`.
    pub fn velocity(&self, theta: f64) -> f64 {
        self.evaluate(theta).velocity
    }

    fn eval_shm(&self, t: f64) -> FollowerState {
        use std::f64::consts::PI;
        let l = self.lift;
        let b = self.beta;
        let phase = PI * t / b;
        FollowerState {
            displacement: 0.5 * l * (1.0 - phase.cos()),
            velocity: (PI * l) / (2.0 * b) * phase.sin(),
            acceleration: (PI * PI * l) / (2.0 * b * b) * phase.cos(),
        }
    }

    fn eval_cycloidal(&self, t: f64) -> FollowerState {
        use std::f64::consts::TAU;
        let l = self.lift;
        let b = self.beta;
        let ratio = t / b;
        let phase = TAU * ratio;
        FollowerState {
            displacement: l * (ratio - phase.sin() / TAU),
            velocity: (l / b) * (1.0 - phase.cos()),
            acceleration: (TAU * l / (b * b)) * phase.sin(),
            // note: 2π = TAU, so (2π L / β²) = TAU * l / b².
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_bad_parameters() {
        assert!(CamRise::new(CamMotion::SimpleHarmonic, 1.0, 0.0).is_err());
        assert!(CamRise::new(CamMotion::SimpleHarmonic, -1.0, 1.0).is_err());
        assert!(CamRise::new(CamMotion::Cycloidal, f64::NAN, 1.0).is_err());
        assert!(CamRise::new(CamMotion::Cycloidal, 0.0, 1.0).is_ok());
    }

    // --- SHM: rise(0)=0, rise(beta)=lift, monotonic ---

    #[test]
    fn shm_endpoints_hit_zero_and_lift() {
        let lift = 12.5;
        let beta = PI; // half a cam turn
        let cam = CamRise::new(CamMotion::SimpleHarmonic, lift, beta).unwrap();
        assert!(cam.rise(0.0).abs() < EPS, "rise(0) = {}", cam.rise(0.0));
        assert!(
            (cam.rise(beta) - lift).abs() < EPS,
            "rise(beta) = {}, want {lift}",
            cam.rise(beta)
        );
    }

    #[test]
    fn shm_midpoint_is_half_lift() {
        // s(β/2) = (L/2)(1 − cos(π/2)) = L/2.
        let lift = 10.0;
        let beta = 2.0;
        let cam = CamRise::new(CamMotion::SimpleHarmonic, lift, beta).unwrap();
        assert!((cam.rise(beta / 2.0) - lift / 2.0).abs() < EPS);
    }

    #[test]
    fn shm_is_monotonic_non_decreasing() {
        let cam = CamRise::new(CamMotion::SimpleHarmonic, 7.3, 1.4).unwrap();
        let n = 500;
        let mut prev = f64::NEG_INFINITY;
        for i in 0..=n {
            let theta = cam.beta * (i as f64) / (n as f64);
            let s = cam.rise(theta);
            assert!(
                s >= prev - 1e-12,
                "SHM not monotonic at theta = {theta}: {s} < {prev}"
            );
            prev = s;
        }
    }

    #[test]
    fn shm_has_zero_velocity_at_endpoints() {
        let cam = CamRise::new(CamMotion::SimpleHarmonic, 5.0, 1.1).unwrap();
        assert!(cam.velocity(0.0).abs() < EPS);
        assert!(cam.velocity(cam.beta).abs() < EPS);
    }

    #[test]
    fn shm_has_nonzero_acceleration_at_endpoints() {
        // The defining contrast with cycloidal: SHM end-accel != 0.
        let cam = CamRise::new(CamMotion::SimpleHarmonic, 5.0, 1.1).unwrap();
        let a0 = cam.evaluate(0.0).acceleration;
        let ab = cam.evaluate(cam.beta).acceleration;
        assert!(
            a0.abs() > 1e-3,
            "SHM start accel should be nonzero, got {a0}"
        );
        assert!(ab.abs() > 1e-3, "SHM end accel should be nonzero, got {ab}");
        // And they are equal/opposite: +π²L/2β² at start, −same at end.
        assert!((a0 + ab).abs() < EPS, "SHM end accels should be opposite");
    }

    // --- Cycloidal: rise(0)=0, rise(beta)=lift, monotonic, zero end velocity ---

    #[test]
    fn cycloidal_endpoints_hit_zero_and_lift() {
        let lift = 8.0;
        let beta = 1.5;
        let cam = CamRise::new(CamMotion::Cycloidal, lift, beta).unwrap();
        assert!(cam.rise(0.0).abs() < EPS);
        assert!((cam.rise(beta) - lift).abs() < EPS);
    }

    #[test]
    fn cycloidal_midpoint_is_half_lift() {
        // s(β/2) = L[1/2 − sin(π)/2π] = L/2 (sin π = 0).
        let lift = 9.0;
        let beta = 2.2;
        let cam = CamRise::new(CamMotion::Cycloidal, lift, beta).unwrap();
        assert!((cam.rise(beta / 2.0) - lift / 2.0).abs() < EPS);
    }

    #[test]
    fn cycloidal_is_monotonic_non_decreasing() {
        let cam = CamRise::new(CamMotion::Cycloidal, 6.0, 1.7).unwrap();
        let n = 500;
        let mut prev = f64::NEG_INFINITY;
        for i in 0..=n {
            let theta = cam.beta * (i as f64) / (n as f64);
            let s = cam.rise(theta);
            assert!(
                s >= prev - 1e-12,
                "cycloidal not monotonic at theta = {theta}: {s} < {prev}"
            );
            prev = s;
        }
    }

    #[test]
    fn cycloidal_has_zero_velocity_at_endpoints() {
        // The headline cycloidal property.
        let cam = CamRise::new(CamMotion::Cycloidal, 6.0, 1.7).unwrap();
        let v0 = cam.velocity(0.0);
        let vb = cam.velocity(cam.beta);
        assert!(v0.abs() < EPS, "cycloidal start velocity = {v0}");
        assert!(vb.abs() < EPS, "cycloidal end velocity = {vb}");
    }

    #[test]
    fn cycloidal_has_zero_acceleration_at_endpoints() {
        // Cycloidal additionally has zero end acceleration (its
        // advantage over SHM): a(0) = a(β) = (2πL/β²)·sin(0 or 2π) = 0.
        let cam = CamRise::new(CamMotion::Cycloidal, 6.0, 1.7).unwrap();
        assert!(cam.evaluate(0.0).acceleration.abs() < EPS);
        assert!(cam.evaluate(cam.beta).acceleration.abs() < EPS);
    }

    #[test]
    fn cycloidal_peak_velocity_matches_closed_form() {
        // Max velocity occurs at β/2 where cos(2π·½)=cos π=−1, giving
        // ds/dθ = (L/β)(1 − (−1)) = 2L/β.
        let lift = 6.0;
        let beta = 1.7;
        let cam = CamRise::new(CamMotion::Cycloidal, lift, beta).unwrap();
        let v_mid = cam.velocity(beta / 2.0);
        assert!((v_mid - 2.0 * lift / beta).abs() < EPS, "v_mid = {v_mid}");
    }

    // --- shared behaviour ---

    #[test]
    fn evaluate_clamps_outside_rise_to_dwells() {
        let lift = 4.0;
        let cam = CamRise::new(CamMotion::Cycloidal, lift, 1.0).unwrap();
        // Before the rise -> bottom dwell.
        assert!(cam.rise(-0.5).abs() < EPS);
        // After the rise -> top dwell.
        assert!((cam.rise(2.0) - lift).abs() < EPS);
    }

    #[test]
    fn zero_lift_is_flat_for_both_programs() {
        for motion in [CamMotion::SimpleHarmonic, CamMotion::Cycloidal] {
            let cam = CamRise::new(motion, 0.0, 1.0).unwrap();
            for i in 0..=10 {
                let theta = (i as f64) / 10.0;
                assert!(cam.rise(theta).abs() < EPS, "{motion:?} not flat");
            }
        }
    }
}
