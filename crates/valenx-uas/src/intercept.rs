//! Defensive counter-UAS **geometry & timeline** — pure kinematics.
//!
//! ## Dual-use boundary (read this)
//!
//! This module is **defensive and geometric only**. It answers two questions
//! about an inbound track, and *nothing else*:
//!
//! 1. **Detection timeline** — given a sensor of range `R` at a fixed site and
//!    a track flying a straight, constant-velocity path, *when* (and where)
//!    does the track first cross inside `R`, how long does it stay inside, and
//!    when does it leave? ([`detection_timeline`])
//! 2. **Intercept geometry** — given that track and an interceptor that can
//!    fly at a maximum speed `s` from a starting point, *what is the earliest
//!    time it could be co-located with the track, and at what point?*
//!    ([`time_to_intercept`])
//!
//! That is the **time and geometry of interception** — the information a
//! defender needs to understand whether and where a threat can be reached. It
//! is the constant-speed pursuit problem solved exactly. **There is no weapon,
//! no targeting solution, no guidance command, no lethality, and no notion of
//! "engagement" beyond two points coinciding in space-time.** What a defender
//! does at the rendezvous is entirely out of scope of this crate.
//!
//! This mirrors the dual-use posture of the civilian and air-traffic
//! conflict-detection / closest-point-of-approach math used for collision
//! avoidance: the same kinematics tell a drone how to *avoid* another aircraft
//! and tell a defender when an inbound track is reachable.
//!
//! ## The constant-speed pursuit quadratic
//!
//! Let the target start at `p_t` with constant velocity `v_t`, and the
//! interceptor start at `p_i` with a constant *speed* `s` (its heading is the
//! unknown we solve for — a constant-bearing / proportional-navigation closing
//! course in the limit). At intercept time `t` the target is at `p_t + v_t·t`;
//! the interceptor, flying straight at speed `s`, can reach any point within
//! `s·t` of `p_i`. Intercept is feasible exactly when the interceptor can be
//! co-located with the target:
//!
//! `‖(p_t − p_i) + v_t·t‖ = s·t`.
//!
//! Squaring gives the quadratic `a·t² + b·t + c = 0` with
//!
//! - `a = ‖v_t‖² − s²`
//! - `b = 2·(d · v_t)`   where `d = p_t − p_i`
//! - `c = ‖d‖²`
//!
//! The earliest non-negative root is the minimum time-to-intercept. Two
//! special cases the [`time_to_intercept`] solver guards exactly:
//!
//! - **Head-on, constant closing.** When the target flies straight at the
//!   interceptor, the closing speed is `s + ‖v_t‖` and `t = ‖d‖ / (s + ‖v_t‖)`
//!   — i.e. `range / v_rel`. (The quadratic reproduces this.)
//! - **Interceptor no faster than the target** (`a >= 0`). Then a stern chase
//!   may be impossible: if the target is opening, there is **no** non-negative
//!   root and the function returns `None` (no solution) rather than a bogus
//!   time.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::error::{require_finite, require_positive, UasError};

/// A threat track: a position and a *constant* velocity (m, m/s).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ThreatTrack {
    /// Current position (m).
    pub position: Vector3<f64>,
    /// Constant velocity (m/s).
    pub velocity: Vector3<f64>,
}

impl ThreatTrack {
    /// Build a validated threat track (all six components must be finite).
    ///
    /// # Errors
    ///
    /// [`UasError::NotFinite`] if any position or velocity component is not
    /// finite.
    pub fn new(position: Vector3<f64>, velocity: Vector3<f64>) -> Result<Self, UasError> {
        finite_vec("threat position", position)?;
        finite_vec("threat velocity", velocity)?;
        Ok(Self { position, velocity })
    }

    /// Position at time `t` (s): `p + v·t`.
    pub fn position_at(&self, t: f64) -> Vector3<f64> {
        self.position + self.velocity * t
    }

    /// Speed `‖v‖` (m/s).
    pub fn speed(&self) -> f64 {
        self.velocity.norm()
    }
}

/// An interceptor: a starting position and a maximum speed (m, m/s).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Interceptor {
    /// Starting position (m).
    pub position: Vector3<f64>,
    /// Maximum speed `s` (m/s), finite and positive.
    pub max_speed_m_s: f64,
}

impl Interceptor {
    /// Build a validated interceptor.
    ///
    /// # Errors
    ///
    /// [`UasError::NotFinite`] if a position component is not finite;
    /// [`UasError::NonPositive`] if `max_speed_m_s` is not finite and positive.
    pub fn new(position: Vector3<f64>, max_speed_m_s: f64) -> Result<Self, UasError> {
        finite_vec("interceptor position", position)?;
        let max_speed_m_s = require_positive("interceptor max_speed_m_s", max_speed_m_s)?;
        Ok(Self {
            position,
            max_speed_m_s,
        })
    }
}

/// A feasible intercept: the earliest time and the point of co-location.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct InterceptSolution {
    /// Earliest time-to-intercept (s), `>= 0`.
    pub time_s: f64,
    /// Intercept point (m) — where target and interceptor coincide.
    pub point: Vector3<f64>,
    /// Straight-line distance the interceptor must cover, `s·t` (m).
    pub interceptor_path_len_m: f64,
    /// Closing speed along the initial line of sight `v_rel` (m/s), reported
    /// for the head-on / constant-bearing intuition. Positive = closing.
    pub initial_closing_speed_m_s: f64,
}

/// Earliest constant-speed intercept of `track` by `interceptor`, or `None` if
/// no non-negative-time solution exists (e.g. the interceptor is no faster than
/// a target that is opening the range).
///
/// Solves the pursuit quadratic `a·t² + b·t + c = 0` (see the module docs) and
/// returns the **smallest non-negative root**. Every divide and the discriminant
/// `sqrt` are guarded:
///
/// - If the interceptor already starts *on* the target (`c = 0`), `t = 0`.
/// - If `a ≈ 0` (interceptor speed equals target speed), the equation is
///   *linear* `b·t + c = 0`; solved directly, with `b ≈ 0` (and `c > 0`)
///   meaning never-closing → `None`.
/// - If the discriminant is negative, there is no real intercept → `None`.
/// - Among real roots, the smallest `>= 0` (within a tiny tolerance) is taken;
///   if both are negative → `None`.
///
/// This is pure kinematics: it tells you *when and where* an intercept is
/// geometrically possible, not how to command a vehicle or what to do there.
pub fn time_to_intercept(
    track: &ThreatTrack,
    interceptor: &Interceptor,
) -> Result<Option<InterceptSolution>, UasError> {
    // Validate (the constructors already did, but accept bare structs too).
    finite_vec("threat position", track.position)?;
    finite_vec("threat velocity", track.velocity)?;
    finite_vec("interceptor position", interceptor.position)?;
    let s = require_positive("interceptor max_speed_m_s", interceptor.max_speed_m_s)?;

    let d = track.position - interceptor.position; // p_t - p_i
    let vt = track.velocity;

    let a = vt.norm_squared() - s * s;
    let b = 2.0 * d.dot(&vt);
    let c = d.norm_squared();

    // Initial closing speed along the line of sight (for reporting): the
    // target's velocity component TOWARD the interceptor is -(d̂ · v_t); a
    // positive value means the range is shrinking from the target's motion.
    let range0 = c.sqrt();
    let initial_closing_speed = if range0 > f64::EPSILON {
        -(d.dot(&vt)) / range0
    } else {
        0.0
    };

    // Already co-located.
    if c <= 0.0 {
        return Ok(Some(InterceptSolution {
            time_s: 0.0,
            point: track.position,
            interceptor_path_len_m: 0.0,
            initial_closing_speed_m_s: initial_closing_speed,
        }));
    }

    // Pick the earliest non-negative root, guarding the linear (a≈0) case and
    // the discriminant.
    const TOL: f64 = 1e-9;
    let t = if a.abs() < TOL {
        // Linear: b t + c = 0.  c > 0 here.
        if b.abs() < TOL {
            // No t dependence and c > 0: never closes.
            None
        } else {
            let t = -c / b;
            if t >= -TOL {
                Some(t.max(0.0))
            } else {
                None
            }
        }
    } else {
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            None
        } else {
            let sqrt_disc = disc.sqrt();
            let t1 = (-b - sqrt_disc) / (2.0 * a);
            let t2 = (-b + sqrt_disc) / (2.0 * a);
            smallest_nonneg(t1, t2, TOL)
        }
    };

    match t {
        None => Ok(None),
        Some(t) => {
            let point = track.position_at(t);
            Ok(Some(InterceptSolution {
                time_s: t,
                point,
                interceptor_path_len_m: s * t,
                initial_closing_speed_m_s: initial_closing_speed,
            }))
        }
    }
}

/// Pick the smallest root that is `>= -tol`, clamped to `0`, or `None` if both
/// are negative.
fn smallest_nonneg(t1: f64, t2: f64, tol: f64) -> Option<f64> {
    let (lo, hi) = if t1 <= t2 { (t1, t2) } else { (t2, t1) };
    if lo >= -tol {
        Some(lo.max(0.0))
    } else if hi >= -tol {
        Some(hi.max(0.0))
    } else {
        None
    }
}

/// When a fixed sensor of range `R` first detects a constant-velocity inbound
/// track, how long it stays inside, and when it leaves.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DetectionTimeline {
    /// Whether the track is ever within sensor range at any `t >= 0`.
    pub detected: bool,
    /// Time of first detection (s), `>= 0`, if `detected`. `0` if it starts
    /// inside the range.
    pub first_detect_s: Option<f64>,
    /// Time the track leaves the range again (s), if it does. `None` if it is
    /// still inbound/inside as `t → ∞` (e.g. heading straight in and stopping
    /// is not modelled — a constant-velocity track that ends inside stays).
    pub last_in_range_s: Option<f64>,
    /// Closest-point-of-approach distance to the sensor (m).
    pub min_range_m: f64,
    /// Time of closest approach (s), `>= 0`.
    pub time_of_closest_approach_s: f64,
}

/// Compute the [`DetectionTimeline`] for a track against a sensor of range `R`
/// sited at `sensor_position`.
///
/// The track flies `p(t) = p0 + v·t`. Detection is the interval of `t >= 0`
/// where `‖p(t) − sensor‖ <= R`. The squared range `‖p(t)−c‖² = R²` is a
/// quadratic in `t`; its real roots (clamped to `t >= 0`) bound the in-range
/// interval, and the unconstrained minimum of the squared range gives the
/// closest point of approach. All divides and the `sqrt` are guarded; a
/// stationary track (`v = 0`) is handled (it is either always in range or
/// never).
///
/// Pure geometry: this is *when the sensor could first see the track*, nothing
/// about classification, tracking, or response.
///
/// # Errors
///
/// [`UasError::NotFinite`] if a sensor component is non-finite;
/// [`UasError::NonPositive`] if `sensor_range_m` is not finite and positive.
pub fn detection_timeline(
    track: &ThreatTrack,
    sensor_position: Vector3<f64>,
    sensor_range_m: f64,
) -> Result<DetectionTimeline, UasError> {
    finite_vec("threat position", track.position)?;
    finite_vec("threat velocity", track.velocity)?;
    finite_vec("sensor position", sensor_position)?;
    let r = require_positive("sensor_range_m", sensor_range_m)?;

    let d0 = track.position - sensor_position; // relative position at t=0
    let v = track.velocity;

    // Squared range f(t) = ‖d0 + v t‖² = (v·v) t² + 2 (d0·v) t + (d0·d0).
    let a = v.norm_squared();
    let b = 2.0 * d0.dot(&v);
    let c = d0.norm_squared();
    let r2 = r * r;

    // Time of closest approach: minimum of f(t). If the track is moving
    // (a > 0), t* = -b/(2a), clamped to t >= 0; if stationary, t* = 0.
    const TOL: f64 = 1e-12;
    let tca = if a > TOL {
        (-b / (2.0 * a)).max(0.0)
    } else {
        0.0
    };
    let f_tca = (a * tca + b) * tca + c; // f(tca), Horner; >= 0
    let min_range = f_tca.max(0.0).sqrt();

    // In-range interval: solve f(t) = R²  ->  a t² + b t + (c - R²) = 0.
    let cc = c - r2;
    let starts_inside = c <= r2;

    // Find the [t_enter, t_leave] interval intersected with [0, ∞).
    let (detected, first_detect, last_in_range) = if a <= TOL {
        // Stationary track: constant range sqrt(c). In range for all t or none.
        if starts_inside {
            (true, Some(0.0), None) // always in range
        } else {
            (false, None, None)
        }
    } else {
        let disc = b * b - 4.0 * a * cc;
        if disc < 0.0 {
            // f(t) never reaches R²: either always inside (cc<0 impossible here
            // since min f = f_tca and disc<0 means f_tca > R²) => never in range.
            (false, None, None)
        } else {
            let sqrt_disc = disc.sqrt();
            let mut t_enter = (-b - sqrt_disc) / (2.0 * a);
            let mut t_leave = (-b + sqrt_disc) / (2.0 * a);
            if t_enter > t_leave {
                std::mem::swap(&mut t_enter, &mut t_leave);
            }
            // Intersect [t_enter, t_leave] with [0, ∞).
            if t_leave < 0.0 {
                // The whole in-range window is in the past.
                (false, None, None)
            } else {
                let first = t_enter.max(0.0);
                // first == 0 if it started inside; this matches starts_inside.
                let _ = starts_inside;
                (true, Some(first), Some(t_leave))
            }
        }
    };

    Ok(DetectionTimeline {
        detected,
        first_detect_s: first_detect,
        last_in_range_s: last_in_range,
        min_range_m: min_range,
        time_of_closest_approach_s: tca,
    })
}

/// Validate that every component of a vector is finite.
fn finite_vec(name: &'static str, v: Vector3<f64>) -> Result<(), UasError> {
    require_finite(name, v.x)?;
    require_finite(name, v.y)?;
    require_finite(name, v.z)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-6 * b.abs().max(1.0)
    }

    // ---- BENCHMARK PIN: head-on constant closing, t = range / v_rel -------
    #[test]
    fn head_on_intercept_is_range_over_closing_speed() {
        // Target 1000 m down +x, flying straight back toward the interceptor
        // at 20 m/s (-x). Interceptor at origin, 80 m/s. They close head-on at
        // 80 + 20 = 100 m/s, so t = 1000 / 100 = 10 s.
        let track = ThreatTrack::new(
            Vector3::new(1000.0, 0.0, 0.0),
            Vector3::new(-20.0, 0.0, 0.0),
        )
        .unwrap();
        let icpt = Interceptor::new(Vector3::zeros(), 80.0).unwrap();
        let sol = time_to_intercept(&track, &icpt).unwrap().unwrap();
        assert!(
            close(sol.time_s, 10.0),
            "head-on TTI should be range/v_rel = 10 s, got {}",
            sol.time_s
        );
        // The interceptor flies straight out +x to the meeting point.
        assert!(close(sol.point.x, 800.0)); // target moved 200 m back
        assert!(close(sol.interceptor_path_len_m, 800.0)); // 80 m/s * 10 s
        assert!(close(sol.initial_closing_speed_m_s, 20.0)); // target closing component
    }

    // ---- BENCHMARK PIN: perpendicular crossing yields the real quadratic root
    #[test]
    fn perpendicular_crossing_uses_the_quadratic_root() {
        // Target at (0, 1000, 0) crossing in +x at 30 m/s. Interceptor at
        // origin, 50 m/s. d = (0,1000,0), v_t=(30,0,0):
        //   a = 30^2 - 50^2 = 900 - 2500 = -1600
        //   b = 2*(d·v_t) = 0
        //   c = 1000^2 = 1e6
        //   disc = 0 - 4*(-1600)*1e6 = 6.4e9 ; sqrt = 80000
        //   roots = (0 ∓ 80000) / (2*-1600) = -25 and +25  -> t = 25 s.
        let track =
            ThreatTrack::new(Vector3::new(0.0, 1000.0, 0.0), Vector3::new(30.0, 0.0, 0.0)).unwrap();
        let icpt = Interceptor::new(Vector3::zeros(), 50.0).unwrap();
        let sol = time_to_intercept(&track, &icpt).unwrap().unwrap();
        assert!(
            close(sol.time_s, 25.0),
            "perpendicular crossing TTI should be the quadratic root 25 s, got {}",
            sol.time_s
        );
        // Verify the geometry closes: interceptor path length == s*t and the
        // meeting point is exactly reachable in that time.
        let meet = track.position_at(sol.time_s);
        let needed = (meet - icpt.position).norm();
        assert!(
            close(needed, icpt.max_speed_m_s * sol.time_s),
            "meeting point distance {needed} != s*t {}",
            icpt.max_speed_m_s * sol.time_s
        );
    }

    // ---- BENCHMARK PIN: target faster, opening -> no solution (None) ------
    #[test]
    fn opening_target_faster_than_interceptor_has_no_solution() {
        // Target ahead at (100,0,0) flying further away in +x at 60 m/s;
        // interceptor only 40 m/s. a = 60^2-40^2 > 0, target opening -> None.
        let track =
            ThreatTrack::new(Vector3::new(100.0, 0.0, 0.0), Vector3::new(60.0, 0.0, 0.0)).unwrap();
        let icpt = Interceptor::new(Vector3::zeros(), 40.0).unwrap();
        assert!(time_to_intercept(&track, &icpt).unwrap().is_none());
    }

    #[test]
    fn equal_speed_stern_chase_opening_is_none_but_closing_solves() {
        // Equal speeds (a≈0). Opening target ahead -> linear b t + c = 0 with
        // b>0, c>0 -> negative t -> None.
        let opening =
            ThreatTrack::new(Vector3::new(100.0, 0.0, 0.0), Vector3::new(50.0, 0.0, 0.0)).unwrap();
        let icpt = Interceptor::new(Vector3::zeros(), 50.0).unwrap();
        assert!(time_to_intercept(&opening, &icpt).unwrap().is_none());

        // Equal speeds but the target is closing head-on: linear solve gives a
        // finite time = range / (2 s) here (closing speed 2*50 along the line).
        let closing = ThreatTrack::new(
            Vector3::new(1000.0, 0.0, 0.0),
            Vector3::new(-50.0, 0.0, 0.0),
        )
        .unwrap();
        let sol = time_to_intercept(&closing, &icpt).unwrap().unwrap();
        assert!(close(sol.time_s, 10.0)); // 1000 / (50+50)
    }

    #[test]
    fn already_colocated_is_zero_time() {
        let track = ThreatTrack::new(Vector3::zeros(), Vector3::new(10.0, 0.0, 0.0)).unwrap();
        let icpt = Interceptor::new(Vector3::zeros(), 30.0).unwrap();
        let sol = time_to_intercept(&track, &icpt).unwrap().unwrap();
        assert_eq!(sol.time_s, 0.0);
    }

    // ---- Detection timeline ----------------------------------------------
    #[test]
    fn detection_timeline_inbound_track() {
        // Track at (1000,0,0) inbound at -100 m/s along x; sensor at origin,
        // range 500 m. It crosses x=500 at t=5 (enter), x=-500 at t=15 (leave),
        // closest approach is 0 at t=10.
        let track = ThreatTrack::new(
            Vector3::new(1000.0, 0.0, 0.0),
            Vector3::new(-100.0, 0.0, 0.0),
        )
        .unwrap();
        let tl = detection_timeline(&track, Vector3::zeros(), 500.0).unwrap();
        assert!(tl.detected);
        assert!(close(tl.first_detect_s.unwrap(), 5.0));
        assert!(close(tl.last_in_range_s.unwrap(), 15.0));
        assert!(tl.min_range_m < 1e-6); // passes through the sensor
        assert!(close(tl.time_of_closest_approach_s, 10.0));
    }

    #[test]
    fn detection_timeline_offset_track_min_range() {
        // Track on a line y=300 crossing in +x; sensor at origin range 500.
        // CPA is at x=0 (t such that x=0), min range 300 < 500 -> detected.
        let track = ThreatTrack::new(
            Vector3::new(-1000.0, 300.0, 0.0),
            Vector3::new(50.0, 0.0, 0.0),
        )
        .unwrap();
        let tl = detection_timeline(&track, Vector3::zeros(), 500.0).unwrap();
        assert!(tl.detected);
        assert!(close(tl.min_range_m, 300.0));
        // CPA when x = 0: from x=-1000 at 50 m/s -> t = 20 s.
        assert!(close(tl.time_of_closest_approach_s, 20.0));
    }

    #[test]
    fn detection_timeline_track_that_never_enters() {
        // Track passes no closer than 800 m to a 500 m sensor -> never detected.
        let track = ThreatTrack::new(
            Vector3::new(-1000.0, 800.0, 0.0),
            Vector3::new(50.0, 0.0, 0.0),
        )
        .unwrap();
        let tl = detection_timeline(&track, Vector3::zeros(), 500.0).unwrap();
        assert!(!tl.detected);
        assert!(tl.first_detect_s.is_none());
        assert!(close(tl.min_range_m, 800.0));
    }

    #[test]
    fn detection_timeline_starts_inside() {
        // Track starts 100 m from the sensor (inside 500 m) flying out; first
        // detection is t=0, and it leaves at some positive time.
        let track =
            ThreatTrack::new(Vector3::new(100.0, 0.0, 0.0), Vector3::new(50.0, 0.0, 0.0)).unwrap();
        let tl = detection_timeline(&track, Vector3::zeros(), 500.0).unwrap();
        assert!(tl.detected);
        assert_eq!(tl.first_detect_s.unwrap(), 0.0);
        // Leaves when x = 500: t = (500-100)/50 = 8 s.
        assert!(close(tl.last_in_range_s.unwrap(), 8.0));
    }

    #[test]
    fn detection_timeline_stationary_track() {
        // Stationary track inside range -> always in range (no leave time).
        let inside = ThreatTrack::new(Vector3::new(100.0, 0.0, 0.0), Vector3::zeros()).unwrap();
        let tl = detection_timeline(&inside, Vector3::zeros(), 500.0).unwrap();
        assert!(tl.detected);
        assert!(tl.last_in_range_s.is_none());
        // Stationary track outside range -> never.
        let outside = ThreatTrack::new(Vector3::new(900.0, 0.0, 0.0), Vector3::zeros()).unwrap();
        let tl = detection_timeline(&outside, Vector3::zeros(), 500.0).unwrap();
        assert!(!tl.detected);
    }

    // ---- Fail-loud: degenerate inputs -> Err/None, never a panic ----------
    #[test]
    fn rejects_degenerate_inputs() {
        let good_track =
            ThreatTrack::new(Vector3::new(100.0, 0.0, 0.0), Vector3::new(-10.0, 0.0, 0.0)).unwrap();

        // Zero / negative interceptor speed.
        assert!(Interceptor::new(Vector3::zeros(), 0.0).is_err());
        assert!(Interceptor::new(Vector3::zeros(), -5.0).is_err());
        // Non-finite components.
        assert!(ThreatTrack::new(Vector3::new(f64::NAN, 0.0, 0.0), Vector3::zeros()).is_err());
        assert!(Interceptor::new(Vector3::new(0.0, f64::INFINITY, 0.0), 10.0).is_err());

        // Zero sensor range.
        assert!(detection_timeline(&good_track, Vector3::zeros(), 0.0).is_err());
        assert!(detection_timeline(&good_track, Vector3::zeros(), -100.0).is_err());
        assert!(detection_timeline(&good_track, Vector3::new(f64::NAN, 0.0, 0.0), 100.0).is_err());
    }

    #[test]
    fn solutions_round_trip_through_serde() {
        let track = ThreatTrack::new(
            Vector3::new(1000.0, 0.0, 0.0),
            Vector3::new(-20.0, 0.0, 0.0),
        )
        .unwrap();
        let icpt = Interceptor::new(Vector3::zeros(), 80.0).unwrap();
        let sol = time_to_intercept(&track, &icpt).unwrap().unwrap();
        let json = serde_json::to_string(&sol).unwrap();
        let back: InterceptSolution = serde_json::from_str(&json).unwrap();
        assert!(close(back.time_s, sol.time_s));

        let tl = detection_timeline(&track, Vector3::zeros(), 500.0).unwrap();
        let json = serde_json::to_string(&tl).unwrap();
        let back: DetectionTimeline = serde_json::from_str(&json).unwrap();
        assert_eq!(back.detected, tl.detected);
    }
}
