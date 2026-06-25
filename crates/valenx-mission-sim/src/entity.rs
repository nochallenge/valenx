//! Entities and their **movers** — pure analytic kinematics.
//!
//! An [`Entity`] is a point with a position, a side/team label, a liveness flag,
//! a sensor range, and an engagement range. How it *moves* is decided by its
//! [`Mover`]:
//!
//! * [`Mover::Static`] — never moves.
//! * [`Mover::ConstantVelocity`] — straight-line at a fixed velocity:
//!   `p(t) = p0 + v·(t − t0)`, evaluated **exactly** (the benchmark pin).
//! * [`Mover::Waypoint`] — follows an ordered list of waypoints at a constant
//!   ground speed, holding the final waypoint once the route is exhausted. This
//!   is piecewise-constant-velocity, so each leg is again `p0 + v·Δt` exactly.
//!
//! These are deliberately *analytic* movers, not a dynamics integrator: there is
//! no mass, force, drag, or control loop. A constructive mission sim needs to
//! know *where each entity is at time `t`*, and these closed forms answer that
//! exactly and reproducibly. (For force-based dynamics use `valenx-vehicle` /
//! `valenx-mbd`; for a kinematic sensor step-harness, `valenx-sensors`.) The
//! roadmap's `flight6dof` lives in `valenx-astro` as a rocket-pointing demo, not
//! a general mover, so it is intentionally not used here.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::{require_finite, require_non_negative, MissionError};

/// Which side an entity belongs to. Two sides are enough for Lanchester-style
/// force-on-force; the label is opaque to the framework.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Side {
    /// The "blue" side (conventionally friendly / own force).
    Blue,
    /// The "red" side (conventionally the opposing force).
    Red,
}

impl Side {
    /// The other side.
    #[must_use]
    pub fn opposite(self) -> Side {
        match self {
            Side::Blue => Side::Red,
            Side::Red => Side::Blue,
        }
    }
}

/// How an entity moves through space — analytic, closed-form kinematics.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Mover {
    /// Does not move; position is always its start position.
    Static,
    /// Straight-line constant velocity (m/s) from the entity's start position.
    ConstantVelocity(Vector3<f64>),
    /// Follow an ordered route of waypoints at a constant ground speed (m/s),
    /// holding the last waypoint after arrival.
    Waypoint {
        /// Ordered waypoints to visit (m). Must be non-empty.
        route: Vec<Vector3<f64>>,
        /// Constant ground speed along the route (m/s), finite and positive.
        speed_m_s: f64,
    },
}

/// A simulation entity: a moving (or static) point with a side, liveness, and
/// sensor / engagement ranges.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Entity {
    /// Start position at `t = 0` (m).
    pub start_position: Vector3<f64>,
    /// Which side the entity is on.
    pub side: Side,
    /// How the entity moves.
    pub mover: Mover,
    /// Detection range of this entity's sensor (m). `0` ⇒ carries no sensor.
    pub sensor_range_m: f64,
    /// Range within which this entity can *engage* an opponent (m). `0` ⇒
    /// cannot engage. Engagement outcomes are abstract & probabilistic — see
    /// [`crate::engagement`].
    pub engagement_range_m: f64,
    /// Per-engagement probability-of-kill *input* (`Pk`) this entity applies to
    /// an opponent within engagement range. An INPUT parameter in `[0, 1]`, not
    /// computed from any lethality model.
    pub pk: f64,
    /// Whether the entity is still alive (engagements may flip this to `false`).
    pub alive: bool,
}

impl Entity {
    /// Build a validated entity.
    ///
    /// # Errors
    ///
    /// - [`MissionError::NotFinite`] if any start-position component is not
    ///   finite, or (for a [`Mover::ConstantVelocity`]) any velocity component
    ///   is not finite, or (for a [`Mover::Waypoint`]) any waypoint component is
    ///   not finite;
    /// - [`MissionError::Negative`] if `sensor_range_m` or `engagement_range_m`
    ///   is negative / not finite;
    /// - [`MissionError::OutOfProbabilityRange`] if `pk` is outside `[0, 1]`;
    /// - [`MissionError::NonPositive`] if a waypoint mover's speed is not finite
    ///   and positive;
    /// - [`MissionError::EmptyRoute`] if a waypoint mover has no waypoints.
    pub fn new(
        start_position: Vector3<f64>,
        side: Side,
        mover: Mover,
        sensor_range_m: f64,
        engagement_range_m: f64,
        pk: f64,
    ) -> Result<Self, MissionError> {
        finite_vec("entity start_position", start_position)?;
        validate_mover(&mover)?;
        require_non_negative("sensor_range_m", sensor_range_m)?;
        require_non_negative("engagement_range_m", engagement_range_m)?;
        crate::error::require_probability("pk", pk)?;
        Ok(Self {
            start_position,
            side,
            mover,
            sensor_range_m,
            engagement_range_m,
            pk,
            alive: true,
        })
    }

    /// Position at simulated time `t` (s), evaluated analytically.
    ///
    /// * [`Mover::Static`] ⇒ the start position.
    /// * [`Mover::ConstantVelocity`] ⇒ `p0 + v·t` exactly.
    /// * [`Mover::Waypoint`] ⇒ the point reached after travelling `speed·t`
    ///   metres along the polyline, clamped to the final waypoint.
    ///
    /// A negative `t` is treated as `t = 0` (the sim never evaluates the past).
    #[must_use]
    pub fn position_at(&self, t: f64) -> Vector3<f64> {
        let t = t.max(0.0);
        match &self.mover {
            Mover::Static => self.start_position,
            Mover::ConstantVelocity(v) => self.start_position + v * t,
            Mover::Waypoint { route, speed_m_s } => {
                position_along_route(self.start_position, route, *speed_m_s, t)
            }
        }
    }

    /// Velocity at time `t` (s). Constant for the constant-velocity mover; the
    /// current-leg velocity for a waypoint mover (zero once the route is done);
    /// zero for a static entity.
    #[must_use]
    pub fn velocity_at(&self, t: f64) -> Vector3<f64> {
        let t = t.max(0.0);
        match &self.mover {
            Mover::Static => Vector3::zeros(),
            Mover::ConstantVelocity(v) => *v,
            Mover::Waypoint { route, speed_m_s } => {
                velocity_along_route(self.start_position, route, *speed_m_s, t)
            }
        }
    }
}

/// The point reached after travelling `speed·t` metres along the polyline
/// `start → route[0] → route[1] → …`, clamped to the last waypoint.
fn position_along_route(
    start: Vector3<f64>,
    route: &[Vector3<f64>],
    speed: f64,
    t: f64,
) -> Vector3<f64> {
    let mut remaining = speed * t; // metres of travel budget
    let mut from = start;
    for &to in route {
        let leg = to - from;
        let leg_len = leg.norm();
        if leg_len <= f64::EPSILON {
            // Degenerate zero-length leg: skip without consuming budget.
            from = to;
            continue;
        }
        if remaining <= leg_len {
            // Stop part-way along this leg.
            return from + leg * (remaining / leg_len);
        }
        remaining -= leg_len;
        from = to;
    }
    // Route exhausted: hold the final waypoint.
    from
}

/// The current-leg velocity vector (direction × speed) at travel time `t`, or
/// zero once the route is exhausted.
fn velocity_along_route(
    start: Vector3<f64>,
    route: &[Vector3<f64>],
    speed: f64,
    t: f64,
) -> Vector3<f64> {
    let mut remaining = speed * t;
    let mut from = start;
    for &to in route {
        let leg = to - from;
        let leg_len = leg.norm();
        if leg_len <= f64::EPSILON {
            from = to;
            continue;
        }
        if remaining < leg_len {
            return leg * (speed / leg_len);
        }
        remaining -= leg_len;
        from = to;
    }
    Vector3::zeros()
}

/// Validate a mover's numeric fields.
fn validate_mover(mover: &Mover) -> Result<(), MissionError> {
    match mover {
        Mover::Static => Ok(()),
        Mover::ConstantVelocity(v) => finite_vec("entity velocity", *v),
        Mover::Waypoint { route, speed_m_s } => {
            if route.is_empty() {
                return Err(MissionError::EmptyRoute);
            }
            for wp in route {
                finite_vec("waypoint", *wp)?;
            }
            crate::error::require_positive("waypoint speed_m_s", *speed_m_s)?;
            Ok(())
        }
    }
}

/// Validate that every component of a vector is finite.
fn finite_vec(name: &'static str, v: Vector3<f64>) -> Result<(), MissionError> {
    require_finite(name, v.x)?;
    require_finite(name, v.y)?;
    require_finite(name, v.z)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-9 * b.abs().max(1.0)
    }
    fn close_vec(a: Vector3<f64>, b: Vector3<f64>) -> bool {
        close(a.x, b.x) && close(a.y, b.y) && close(a.z, b.z)
    }

    // ---- BENCHMARK PIN: constant-velocity position is exactly p0 + v*t -------
    #[test]
    fn constant_velocity_position_is_p0_plus_vt() {
        let p0 = Vector3::new(10.0, -5.0, 2.0);
        let v = Vector3::new(3.0, 0.5, -1.0);
        let e = Entity::new(p0, Side::Blue, Mover::ConstantVelocity(v), 0.0, 0.0, 0.0).unwrap();
        for &t in &[0.0, 1.0, 2.5, 100.0, 1234.5] {
            assert!(close_vec(e.position_at(t), p0 + v * t), "t={t}");
            assert!(close_vec(e.velocity_at(t), v));
        }
    }

    #[test]
    fn static_entity_never_moves() {
        let p0 = Vector3::new(1.0, 2.0, 3.0);
        let e = Entity::new(p0, Side::Red, Mover::Static, 0.0, 0.0, 0.0).unwrap();
        assert!(close_vec(e.position_at(0.0), p0));
        assert!(close_vec(e.position_at(99.0), p0));
        assert!(close_vec(e.velocity_at(99.0), Vector3::zeros()));
    }

    #[test]
    fn waypoint_mover_walks_legs_and_holds_the_end() {
        // Start at origin; go to (100,0,0) then (100,100,0) at 10 m/s.
        let route = vec![
            Vector3::new(100.0, 0.0, 0.0),
            Vector3::new(100.0, 100.0, 0.0),
        ];
        let e = Entity::new(
            Vector3::zeros(),
            Side::Blue,
            Mover::Waypoint {
                route,
                speed_m_s: 10.0,
            },
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        // At t=0, at the start.
        assert!(close_vec(e.position_at(0.0), Vector3::zeros()));
        // 5 s -> 50 m along the first leg.
        assert!(close_vec(e.position_at(5.0), Vector3::new(50.0, 0.0, 0.0)));
        // 10 s -> exactly the first waypoint.
        assert!(close_vec(
            e.position_at(10.0),
            Vector3::new(100.0, 0.0, 0.0)
        ));
        // 15 s -> 50 m up the second leg.
        assert!(close_vec(
            e.position_at(15.0),
            Vector3::new(100.0, 50.0, 0.0)
        ));
        // 20 s -> the final waypoint; the velocity on the 2nd leg is +y*10.
        assert!(close_vec(
            e.position_at(20.0),
            Vector3::new(100.0, 100.0, 0.0)
        ));
        assert!(close_vec(e.velocity_at(15.0), Vector3::new(0.0, 10.0, 0.0)));
        // Past the end: holds, velocity zero.
        assert!(close_vec(
            e.position_at(1000.0),
            Vector3::new(100.0, 100.0, 0.0)
        ));
        assert!(close_vec(e.velocity_at(1000.0), Vector3::zeros()));
    }

    #[test]
    fn negative_time_is_clamped_to_zero() {
        let p0 = Vector3::new(4.0, 4.0, 4.0);
        let v = Vector3::new(1.0, 0.0, 0.0);
        let e = Entity::new(p0, Side::Blue, Mover::ConstantVelocity(v), 0.0, 0.0, 0.0).unwrap();
        assert!(close_vec(e.position_at(-10.0), p0));
    }

    #[test]
    fn rejects_degenerate_entities() {
        // Empty route.
        assert!(matches!(
            Entity::new(
                Vector3::zeros(),
                Side::Blue,
                Mover::Waypoint {
                    route: vec![],
                    speed_m_s: 5.0
                },
                0.0,
                0.0,
                0.0,
            ),
            Err(MissionError::EmptyRoute)
        ));
        // Non-positive waypoint speed.
        assert!(Entity::new(
            Vector3::zeros(),
            Side::Blue,
            Mover::Waypoint {
                route: vec![Vector3::new(1.0, 0.0, 0.0)],
                speed_m_s: 0.0
            },
            0.0,
            0.0,
            0.0,
        )
        .is_err());
        // Pk out of range.
        assert!(Entity::new(Vector3::zeros(), Side::Blue, Mover::Static, 0.0, 0.0, 1.5).is_err());
        // Negative sensor range.
        assert!(Entity::new(Vector3::zeros(), Side::Blue, Mover::Static, -1.0, 0.0, 0.0).is_err());
        // Non-finite velocity.
        assert!(Entity::new(
            Vector3::zeros(),
            Side::Blue,
            Mover::ConstantVelocity(Vector3::new(f64::NAN, 0.0, 0.0)),
            0.0,
            0.0,
            0.0,
        )
        .is_err());
    }

    #[test]
    fn side_opposite_round_trips() {
        assert_eq!(Side::Blue.opposite(), Side::Red);
        assert_eq!(Side::Red.opposite(), Side::Blue);
        assert_eq!(Side::Blue.opposite().opposite(), Side::Blue);
    }
}
