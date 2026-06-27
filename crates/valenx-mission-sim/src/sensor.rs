//! Range-based detection between two entities.
//!
//! Detection here is the same **pure geometry** the rest of valenx uses for
//! conflict / closest-point-of-approach analysis: an observer with a sensor of
//! range `R` detects a target the instant the target crosses inside `R`. The
//! heavy lifting is **reused from `valenx-uas`** — [`valenx_uas::detection_timeline`]
//! solves exactly "when does a constant-velocity track first enter range `R` of
//! a fixed sensor, when does it leave, and what is the closest approach" — so we
//! do not reimplement that quadratic.
//!
//! ## Two moving entities → one relative track
//!
//! `valenx-uas`'s detection geometry is written for a *fixed* sensor and a
//! moving track. Both of our entities may move, so we transform into the
//! observer's frame at the evaluation epoch: the relative position is
//! `p_target − p_observer` and the relative velocity is `v_target − v_observer`,
//! and the observer becomes a fixed sensor at the origin. Because both movers
//! are (piecewise) constant-velocity, the relative motion is constant-velocity
//! over the current leg, so the closed-form crossing time is exact — the
//! benchmark pin checks the geometric first-detect time against a hand
//! computation.
//!
//! This is **defensive sensing only**: a detection is "observer could first see
//! target at time `t`". There is no classification, tracking filter, fusion, or
//! any cueing of a weapon — entirely consistent with the dual-use boundary.

use nalgebra::Vector3;
use valenx_uas::{detection_timeline, ThreatTrack};

use crate::entity::Entity;
use crate::error::MissionError;

/// The geometric detection of one entity by another's sensor, expressed as a
/// time relative to the epoch at which it was computed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Detection {
    /// Whether the target is ever within sensor range at any `Δt >= 0` from the
    /// evaluation epoch.
    pub detected: bool,
    /// Time *after the evaluation epoch* (s) at which the observer first detects
    /// the target, if it does. `0` if already in range at the epoch.
    pub first_detect_dt_s: Option<f64>,
    /// Closest-point-of-approach distance over the (relative) straight-line
    /// motion (m).
    pub min_range_m: f64,
}

/// Compute when `observer` first detects `target`, evaluated at simulated time
/// `now` and reported as a Δt forward from `now`.
///
/// The observer's sensor range is [`Entity::sensor_range_m`]; an observer with a
/// zero sensor range (or that is not alive) never detects. Both entities are
/// frozen to their straight-line motion at `now` (their current constant-velocity
/// leg), transformed into the observer's frame, and handed to
/// [`valenx_uas::detection_timeline`].
///
/// # Errors
///
/// [`MissionError::Uas`] if the reused detection geometry rejects the derived
/// inputs (it validates finiteness and a positive range internally).
pub fn detect(observer: &Entity, target: &Entity, now: f64) -> Result<Detection, MissionError> {
    // A dead observer or one with no sensor never detects.
    if !observer.alive || observer.sensor_range_m <= 0.0 {
        let rel = target.position_at(now) - observer.position_at(now);
        return Ok(Detection {
            detected: false,
            first_detect_dt_s: None,
            min_range_m: rel.norm(),
        });
    }

    // Relative track in the observer's frame at `now`: the observer is a fixed
    // sensor at the origin, the target carries the *relative* velocity.
    let rel_pos = target.position_at(now) - observer.position_at(now);
    let rel_vel = target.velocity_at(now) - observer.velocity_at(now);
    let track = ThreatTrack::new(rel_pos, rel_vel)?;

    let tl = detection_timeline(&track, Vector3::zeros(), observer.sensor_range_m)?;
    Ok(Detection {
        detected: tl.detected,
        first_detect_dt_s: tl.first_detect_s,
        min_range_m: tl.min_range_m,
    })
}

/// The instantaneous slant range between two entities at simulated time `now`
/// (m) — a cheap helper used by the engagement check.
#[must_use]
pub fn range_between(a: &Entity, b: &Entity, now: f64) -> f64 {
    (a.position_at(now) - b.position_at(now)).norm()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entity::{Mover, Side};

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-6 * b.abs().max(1.0)
    }

    fn observer_with_range(r: f64) -> Entity {
        Entity::new(Vector3::zeros(), Side::Blue, Mover::Static, r, 0.0, 0.0).unwrap()
    }

    // ---- BENCHMARK PIN: first detection at the geometric range-crossing time -
    #[test]
    fn inbound_target_first_detected_at_range_crossing() {
        // Static observer at origin, sensor range 500 m. Target at (1000,0,0)
        // inbound at -100 m/s. It crosses x=500 at Δt = (1000-500)/100 = 5 s.
        let observer = observer_with_range(500.0);
        let target = Entity::new(
            Vector3::new(1000.0, 0.0, 0.0),
            Side::Red,
            Mover::ConstantVelocity(Vector3::new(-100.0, 0.0, 0.0)),
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        let d = detect(&observer, &target, 0.0).unwrap();
        assert!(d.detected);
        assert!(close(d.first_detect_dt_s.unwrap(), 5.0));
        assert!(d.min_range_m < 1e-6); // passes straight through the origin
    }

    #[test]
    fn relative_motion_accounts_for_a_moving_observer() {
        // Observer moves +x at 40 m/s; target moves -x at 60 m/s, 1000 m ahead.
        // Closing speed = 100 m/s, so with range 500 the crossing is at the same
        // Δt = 5 s as the static case — proving the relative-velocity transform.
        let observer = Entity::new(
            Vector3::zeros(),
            Side::Blue,
            Mover::ConstantVelocity(Vector3::new(40.0, 0.0, 0.0)),
            500.0,
            0.0,
            0.0,
        )
        .unwrap();
        let target = Entity::new(
            Vector3::new(1000.0, 0.0, 0.0),
            Side::Red,
            Mover::ConstantVelocity(Vector3::new(-60.0, 0.0, 0.0)),
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        let d = detect(&observer, &target, 0.0).unwrap();
        assert!(d.detected);
        assert!(close(d.first_detect_dt_s.unwrap(), 5.0));
    }

    #[test]
    fn target_already_in_range_detects_immediately() {
        let observer = observer_with_range(500.0);
        let target = Entity::new(
            Vector3::new(100.0, 0.0, 0.0),
            Side::Red,
            Mover::Static,
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        let d = detect(&observer, &target, 0.0).unwrap();
        assert!(d.detected);
        assert_eq!(d.first_detect_dt_s.unwrap(), 0.0);
    }

    #[test]
    fn target_that_stays_outside_is_never_detected() {
        // Passes no closer than 800 m to a 500 m sensor.
        let observer = observer_with_range(500.0);
        let target = Entity::new(
            Vector3::new(-1000.0, 800.0, 0.0),
            Side::Red,
            Mover::ConstantVelocity(Vector3::new(50.0, 0.0, 0.0)),
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        let d = detect(&observer, &target, 0.0).unwrap();
        assert!(!d.detected);
        assert!(d.first_detect_dt_s.is_none());
        assert!(close(d.min_range_m, 800.0));
    }

    #[test]
    fn sensorless_or_dead_observer_never_detects() {
        let mut blind = observer_with_range(0.0);
        let target = Entity::new(
            Vector3::new(10.0, 0.0, 0.0),
            Side::Red,
            Mover::Static,
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        assert!(!detect(&blind, &target, 0.0).unwrap().detected);

        // A dead observer with a sensor also never detects.
        blind = observer_with_range(500.0);
        blind.alive = false;
        assert!(!detect(&blind, &target, 0.0).unwrap().detected);
    }

    #[test]
    fn range_between_is_euclidean() {
        let a = observer_with_range(0.0);
        let b = Entity::new(
            Vector3::new(3.0, 4.0, 0.0),
            Side::Red,
            Mover::Static,
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        assert!(close(range_between(&a, &b, 0.0), 5.0));
    }
}
