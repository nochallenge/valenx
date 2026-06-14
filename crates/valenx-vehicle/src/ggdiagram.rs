//! **g–g acceleration envelope** (friction-circle diagram) — the combined
//! longitudinal-vs-lateral grip budget a vehicle-dynamics engineer reads to see
//! how much braking/accelerating is left while cornering.
//!
//! A tyre can only deliver so much horizontal force; spend it sideways
//! (cornering) and less is left fore-and-aft (braking / accelerating). The
//! limit is the **friction circle**: `a_lon² + a_lat² ≤ a_max²`, with the
//! radius `a_max` the peak grip (speed-dependent through downforce). This
//! samples the boundary at a given speed: for each lateral g, the available
//! braking g (friction-limited) and acceleration g (also capped by the
//! power-limited straight-line acceleration).
//!
//! Reuses [`Car`]'s force model:
//! [`Car::max_lateral_g`] sets the friction-circle radius and
//! [`Car::acceleration_at`] the straight-line acceleration cap.
//!
//! Validated: the radius equals `max_lateral_g`; at zero lateral the full
//! braking grip is available; at maximum lateral no fore-aft g remains;
//! acceleration never exceeds braking (it is power-capped); and every boundary
//! point lies on (braking) or inside (acceleration) the friction circle.
//!
//! Honest scope: an isotropic friction-circle model on top of `Car`'s
//! point-mass forces — research / preliminary-design grade. It does not model
//! a measured (anisotropic) tyre `Fx–Fy` curve, load sensitivity per tyre, or
//! transient slip; a step toward, not an equal of, an Adams/Car tyre model.

use crate::Car;

/// Standard gravity (m/s²) for converting accelerations to g.
const G: f64 = 9.806_65;

/// One point on the g–g boundary at a fixed lateral demand.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GgPoint {
    /// Lateral acceleration (g).
    pub lateral_g: f64,
    /// Maximum simultaneous *acceleration* (forward) available (g) — the
    /// friction-circle value, capped by the power-limited straight-line value.
    pub accel_g: f64,
    /// Maximum simultaneous *braking* deceleration available (g, magnitude) —
    /// the friction-circle value (grip-limited).
    pub brake_g: f64,
}

/// The sampled g–g envelope at a speed: the friction-circle radius and the
/// boundary points from zero lateral up to the grip limit.
#[derive(Clone, Debug)]
pub struct GgEnvelope {
    /// Speed the envelope was evaluated at (m/s).
    pub speed: f64,
    /// Friction-circle radius — peak combined horizontal grip (g).
    pub friction_g: f64,
    /// Boundary points, lateral g ascending from 0 to `friction_g`.
    pub points: Vec<GgPoint>,
}

impl Car {
    /// Sample the g–g (friction-circle) envelope at `speed` (m/s) with `samples`
    /// lateral steps (`samples.max(2)`), from straight-line up to the grip
    /// limit. Acceleration is also capped by the power-limited straight-line
    /// acceleration; braking is grip-limited.
    pub fn gg_envelope(&self, speed: f64, samples: usize) -> GgEnvelope {
        let a_max = self.max_lateral_g(speed); // friction-circle radius (g)
        let straight = (self.acceleration_at(speed) / G).max(0.0); // power-limited accel (g)
        let n = samples.max(2);
        let mut points = Vec::with_capacity(n);
        for i in 0..n {
            let lat = a_max * i as f64 / (n - 1) as f64;
            let avail = (a_max * a_max - lat * lat).max(0.0).sqrt();
            points.push(GgPoint {
                lateral_g: lat,
                accel_g: avail.min(straight),
                brake_g: avail,
            });
        }
        GgEnvelope {
            speed,
            friction_g: a_max,
            points,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{hypercar, sports_car};

    #[test]
    fn radius_equals_max_lateral_g() {
        let car = sports_car();
        let v = 30.0;
        let env = car.gg_envelope(v, 16);
        assert!(
            (env.friction_g - car.max_lateral_g(v)).abs() < 1e-12,
            "radius {} vs max_lateral_g {}",
            env.friction_g,
            car.max_lateral_g(v)
        );
        assert!(env.friction_g > 0.5, "a sports car should pull > 0.5 g");
    }

    #[test]
    fn straight_line_and_cornering_extremes_are_correct() {
        let car = sports_car();
        let v = 30.0;
        let env = car.gg_envelope(v, 32);
        let first = env.points.first().unwrap();
        let last = env.points.last().unwrap();
        // At zero lateral, the full braking grip is available.
        assert!(first.lateral_g.abs() < 1e-9);
        assert!(
            (first.brake_g - env.friction_g).abs() < 1e-9,
            "straight-line braking should equal the friction radius"
        );
        // At maximum lateral, no fore-aft acceleration remains.
        assert!(
            (last.lateral_g - env.friction_g).abs() < 1e-9,
            "last point reaches the grip limit"
        );
        assert!(
            last.accel_g.abs() < 1e-6 && last.brake_g.abs() < 1e-6,
            "no long. g at the cornering limit: accel {}, brake {}",
            last.accel_g,
            last.brake_g
        );
    }

    #[test]
    fn every_point_respects_the_friction_circle() {
        let car = hypercar();
        let env = car.gg_envelope(40.0, 64);
        for p in &env.points {
            // Braking rides the circle exactly; acceleration is at or inside it.
            let brake_r = (p.lateral_g * p.lateral_g + p.brake_g * p.brake_g).sqrt();
            assert!(
                (brake_r - env.friction_g).abs() < 1e-6,
                "braking point off the circle: r {brake_r} vs {}",
                env.friction_g
            );
            let accel_r = (p.lateral_g * p.lateral_g + p.accel_g * p.accel_g).sqrt();
            assert!(
                accel_r <= env.friction_g + 1e-6,
                "accel point outside the circle: r {accel_r} vs {}",
                env.friction_g
            );
            // Acceleration is power-capped, so it never exceeds the grip-limited
            // braking value at the same lateral demand.
            assert!(
                p.accel_g <= p.brake_g + 1e-9,
                "accel {} > brake {}",
                p.accel_g,
                p.brake_g
            );
        }
    }

    #[test]
    fn available_longitudinal_g_decreases_with_lateral_demand() {
        let env = sports_car().gg_envelope(30.0, 24);
        for w in env.points.windows(2) {
            assert!(
                w[1].brake_g <= w[0].brake_g + 1e-9,
                "braking grip must fall as lateral rises: {} then {}",
                w[0].brake_g,
                w[1].brake_g
            );
        }
    }
}
