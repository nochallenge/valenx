//! Ascent steering: a vertical-rise + pitch-kick + gravity-turn program.
//!
//! This is the classic open-loop launch profile. The vehicle holds
//! a vertical attitude through the dense low atmosphere, executes a
//! small "pitch kick" toward the downrange direction to start the
//! turn, then thrusts along its inertial velocity vector — a gravity
//! turn, where gravity itself does the pitching and the thrust never
//! fights the airflow. It is self-stabilising and needs only three
//! numbers, which is the right altitude for a v1 trajectory tool.

use nalgebra::Vector2;
use serde::{Deserialize, Serialize};

use crate::constants::OMEGA_EARTH;
use crate::error::AstroError;

/// Open-loop pitch program for a planar gravity-turn ascent.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GuidanceProgram {
    /// How long to hold a vertical attitude after liftoff (s).
    pub vertical_rise_time: f64,
    /// Pitch-over angle off vertical to kick the turn (degrees).
    pub pitch_kick_deg: f64,
    /// How long to hold the kicked attitude before handing over to the
    /// pure gravity turn (s).
    pub kick_duration: f64,
}

impl Default for GuidanceProgram {
    fn default() -> Self {
        Self {
            vertical_rise_time: 10.0,
            pitch_kick_deg: 6.0,
            kick_duration: 8.0,
        }
    }
}

impl GuidanceProgram {
    /// Unit thrust direction in the inertial plane at time `t` for a
    /// vehicle at `position` moving at `velocity`.
    ///
    /// The local "up" is the radial unit vector and the downrange
    /// ("east") direction is that rotated +90°, so a positive pitch
    /// kick tilts the thrust prograde for an eastward launch.
    pub fn thrust_direction(
        &self,
        position: Vector2<f64>,
        velocity: Vector2<f64>,
        t: f64,
    ) -> Vector2<f64> {
        let r_hat = safe_normalize(position, Vector2::new(1.0, 0.0));
        let east_hat = Vector2::new(-r_hat.y, r_hat.x);

        if t < self.vertical_rise_time {
            return r_hat;
        }
        if t < self.vertical_rise_time + self.kick_duration {
            let ang = self.pitch_kick_deg.to_radians();
            return safe_normalize(ang.cos() * r_hat + ang.sin() * east_hat, r_hat);
        }
        // Gravity turn: thrust prograde along the velocity *relative to
        // the rotating surface/air*, not the inertial velocity. The
        // Earth's eastward rotation (≈465 m/s at the equator) otherwise
        // swamps the small vertical velocity just after liftoff and the
        // turn would pitch almost flat at once. Relative velocity starts
        // vertical and tilts as the vehicle builds downrange speed over
        // the ground — exactly how a gravity turn is flown.
        let surface_vel = Vector2::new(-OMEGA_EARTH * position.y, OMEGA_EARTH * position.x);
        let relative = velocity - surface_vel;
        safe_normalize(relative, r_hat)
    }

    /// Validate the program parameters.
    pub fn validate(&self) -> Result<(), AstroError> {
        if !self.vertical_rise_time.is_finite() || self.vertical_rise_time < 0.0 {
            return Err(AstroError::InvalidGuidance("vertical_rise_time < 0"));
        }
        if !self.kick_duration.is_finite() || self.kick_duration < 0.0 {
            return Err(AstroError::InvalidGuidance("kick_duration < 0"));
        }
        if !self.pitch_kick_deg.is_finite()
            || self.pitch_kick_deg < 0.0
            || self.pitch_kick_deg >= 90.0
        {
            return Err(AstroError::InvalidGuidance("pitch_kick_deg out of [0,90)"));
        }
        Ok(())
    }
}

/// Normalize a vector, returning `fallback` if its norm is too small to
/// give a meaningful direction.
fn safe_normalize(v: Vector2<f64>, fallback: Vector2<f64>) -> Vector2<f64> {
    let n = v.norm();
    if n > 1e-12 {
        v / n
    } else {
        fallback
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vertical_during_rise_phase() {
        let g = GuidanceProgram::default();
        let pos = Vector2::new(6_378_137.0, 0.0);
        let vel = Vector2::new(0.0, 100.0);
        let dir = g.thrust_direction(pos, vel, 1.0);
        // Should point radially "up" (+x here).
        assert!((dir.x - 1.0).abs() < 1e-9);
        assert!(dir.y.abs() < 1e-9);
    }

    #[test]
    fn kick_tilts_prograde() {
        let g = GuidanceProgram {
            vertical_rise_time: 10.0,
            pitch_kick_deg: 10.0,
            kick_duration: 8.0,
        };
        let pos = Vector2::new(6_378_137.0, 0.0);
        let vel = Vector2::new(50.0, 0.0);
        let dir = g.thrust_direction(pos, vel, 12.0);
        // 10° off vertical: x = cos10°, y = sin10° (east).
        assert!((dir.x - 10f64.to_radians().cos()).abs() < 1e-9);
        assert!((dir.y - 10f64.to_radians().sin()).abs() < 1e-9);
        assert!((dir.norm() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn gravity_turn_follows_surface_relative_velocity() {
        let g = GuidanceProgram::default();
        let pos = Vector2::new(6_378_137.0, 0.0);
        // Surface (air) velocity at this point is ω×r = (0, ωR). Add a
        // (300, 400) surface-relative velocity on top; the gravity turn
        // must follow that relative vector -> direction (0.6, 0.8).
        let surface = Vector2::new(-OMEGA_EARTH * pos.y, OMEGA_EARTH * pos.x);
        let vel = surface + Vector2::new(300.0, 400.0);
        let dir = g.thrust_direction(pos, vel, 100.0);
        assert!((dir.x - 0.6).abs() < 1e-9, "dir {dir:?}");
        assert!((dir.y - 0.8).abs() < 1e-9, "dir {dir:?}");
    }

    #[test]
    fn always_unit_length() {
        let g = GuidanceProgram::default();
        let pos = Vector2::new(0.0, 6_378_137.0);
        for t in [0.0, 5.0, 11.0, 50.0, 500.0] {
            let vel = Vector2::new(t, 2.0 * t);
            let dir = g.thrust_direction(pos, vel, t);
            assert!((dir.norm() - 1.0).abs() < 1e-9, "t={t}");
        }
    }
}
