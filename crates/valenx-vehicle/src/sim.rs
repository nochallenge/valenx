//! **Longitudinal performance trace** — a full speed/distance/acceleration
//! time-history of a launch, and the standard drag-strip metrics read off it.
//!
//! [`Car`] already answers single questions (top speed, time to a
//! target speed, braking distance). This records the *whole* launch as a series
//! of [`TracePoint`]s — speed, distance and instantaneous acceleration at each
//! instant — so you can plot the acceleration curve and extract every metric
//! from one pass: **0–100 km/h**, the **quarter-mile** elapsed time + trap
//! speed, and the **terminal speed**.
//!
//! The integrator is the same explicit forward-Euler step on
//! [`Car::acceleration_at`] that
//! [`Car::accelerate_to`] uses, so the two agree;
//! this just keeps the history and derives more from it. Works on any `Car`,
//! including the chassis of an [`EvCar`](crate::EvCar).
//!
//! Validated: the trace's 0–100 km/h time matches the independent
//! `accelerate_to`; the terminal speed matches `top_speed`; the quarter-mile is
//! self-consistent; speed and distance are monotonic; and a hypercar reaches
//! 100 km/h sooner than a sports car.
//!
//! Honest scope: a point-mass longitudinal model (the same forces `Car` uses —
//! power/traction-limited tractive force with weight transfer, aero drag,
//! rolling resistance) — research / preliminary-design grade. No gearshift
//! dynamics, tyre slip transients, or launch/clutch modelling; a step toward,
//! not an equal of, a full vehicle-dynamics suite (Adams/Car).

use crate::Car;

/// 100 km/h in m/s.
const KMH_100: f64 = 100.0 / 3.6;
/// A quarter mile in metres.
const QUARTER_MILE_M: f64 = 402.336;

/// One sample of a longitudinal launch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TracePoint {
    /// Elapsed time (s).
    pub time: f64,
    /// Speed (m/s).
    pub speed: f64,
    /// Distance covered (m).
    pub distance: f64,
    /// Instantaneous longitudinal acceleration (m/s²).
    pub accel: f64,
}

/// A recorded launch: the time-history plus drag-strip metric extraction.
#[derive(Clone, Debug)]
pub struct SprintTrace {
    /// The samples, in time order (first is `t = 0, v = 0`).
    pub points: Vec<TracePoint>,
}

impl SprintTrace {
    /// Time (s) to first reach speed `target` (m/s), linearly interpolated
    /// between bracketing samples. `None` if the launch never reaches it.
    pub fn time_to_speed(&self, target: f64) -> Option<f64> {
        for w in self.points.windows(2) {
            if w[1].speed >= target {
                let (v0, v1) = (w[0].speed, w[1].speed);
                let frac = if (v1 - v0).abs() > 1e-12 {
                    (target - v0) / (v1 - v0)
                } else {
                    0.0
                };
                return Some(w[0].time + frac * (w[1].time - w[0].time));
            }
        }
        None
    }

    /// Time (s) to first cover distance `target` (m), linearly interpolated.
    pub fn time_to_distance(&self, target: f64) -> Option<f64> {
        for w in self.points.windows(2) {
            if w[1].distance >= target {
                let (d0, d1) = (w[0].distance, w[1].distance);
                let frac = if (d1 - d0).abs() > 1e-12 {
                    (target - d0) / (d1 - d0)
                } else {
                    0.0
                };
                return Some(w[0].time + frac * (w[1].time - w[0].time));
            }
        }
        None
    }

    /// Speed (m/s) when distance `target` (m) is first reached, interpolated.
    pub fn speed_at_distance(&self, target: f64) -> Option<f64> {
        for w in self.points.windows(2) {
            if w[1].distance >= target {
                let (d0, d1) = (w[0].distance, w[1].distance);
                let frac = if (d1 - d0).abs() > 1e-12 {
                    (target - d0) / (d1 - d0)
                } else {
                    0.0
                };
                return Some(w[0].speed + frac * (w[1].speed - w[0].speed));
            }
        }
        None
    }

    /// 0–100 km/h time (s), if reached.
    pub fn zero_to_100_kmh(&self) -> Option<f64> {
        self.time_to_speed(KMH_100)
    }

    /// Quarter-mile `(elapsed time s, trap speed m/s)`, if reached.
    pub fn quarter_mile(&self) -> Option<(f64, f64)> {
        Some((
            self.time_to_distance(QUARTER_MILE_M)?,
            self.speed_at_distance(QUARTER_MILE_M)?,
        ))
    }

    /// The final (≈ terminal) speed of the trace (m/s).
    pub fn terminal_speed(&self) -> f64 {
        self.points.last().map(|p| p.speed).unwrap_or(0.0)
    }
}

impl Car {
    /// Integrate a standing-start launch, recording a [`SprintTrace`]. Steps by
    /// `dt` (s) until `max_time` (s), the car nears its top speed, or the
    /// available acceleration falls to ~0. Same forward-Euler integrator as
    /// [`Car::accelerate_to`].
    pub fn longitudinal_trace(&self, dt: f64, max_time: f64) -> SprintTrace {
        let mut points = Vec::new();
        let (mut v, mut t, mut d) = (0.0_f64, 0.0_f64, 0.0_f64);
        if !dt.is_finite() || dt <= 0.0 {
            points.push(TracePoint {
                time: 0.0,
                speed: 0.0,
                distance: 0.0,
                accel: self.acceleration_at(0.0),
            });
            return SprintTrace { points };
        }
        let top = self.top_speed() * 0.999;
        loop {
            let a = self.acceleration_at(v);
            points.push(TracePoint {
                time: t,
                speed: v,
                distance: d,
                accel: a,
            });
            if t >= max_time || v >= top || a <= 1e-4 {
                break;
            }
            v = (v + a * dt).min(top);
            d += v * dt;
            t += dt;
        }
        SprintTrace { points }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{hypercar, sports_car};

    #[test]
    fn trace_zero_to_100_matches_accelerate_to() {
        // The trace and the single-target accelerate_to use the same integrator,
        // so their 0–100 km/h times must agree.
        let car = sports_car();
        let trace = car.longitudinal_trace(1.0e-3, 60.0);
        let t_trace = trace.zero_to_100_kmh().expect("reaches 100 km/h");
        let t_direct = car.accelerate_to(KMH_100).time;
        assert!(
            (t_trace - t_direct).abs() < 0.05,
            "trace 0-100 {t_trace} vs accelerate_to {t_direct}"
        );
        assert!(t_trace > 0.5 && t_trace < 30.0, "0-100 sanity: {t_trace}s");
    }

    #[test]
    fn terminal_speed_matches_top_speed() {
        let car = sports_car();
        let trace = car.longitudinal_trace(1.0e-3, 120.0);
        let term = trace.terminal_speed();
        let top = car.top_speed();
        assert!(
            (term - top).abs() < 0.05 * top,
            "terminal {term} vs top_speed {top}"
        );
    }

    #[test]
    fn quarter_mile_is_self_consistent() {
        let car = sports_car();
        let trace = car.longitudinal_trace(1.0e-3, 120.0);
        let (et, trap) = trace.quarter_mile().expect("reaches the quarter mile");
        assert!(et > 0.0 && trap > 0.0, "ET {et}, trap {trap}");
        // ET equals the time-to-distance at 402.336 m by construction.
        assert!((trace.time_to_distance(QUARTER_MILE_M).unwrap() - et).abs() < 1e-9);
        // Trap speed cannot exceed the terminal speed.
        assert!(
            trap <= trace.terminal_speed() + 1e-6,
            "trap {trap} > terminal"
        );
        // A quick sanity envelope on a sports car's quarter mile.
        assert!((8.0..20.0).contains(&et), "1/4-mile ET {et}s out of range");
    }

    #[test]
    fn speed_and_distance_are_monotonic() {
        let trace = sports_car().longitudinal_trace(1.0e-3, 60.0);
        for w in trace.points.windows(2) {
            assert!(w[1].speed >= w[0].speed - 1e-9, "speed must not decrease");
            assert!(
                w[1].distance >= w[0].distance - 1e-9,
                "distance must not decrease"
            );
        }
    }

    #[test]
    fn hypercar_reaches_100_sooner_than_a_sports_car() {
        let hc = hypercar()
            .longitudinal_trace(1.0e-3, 60.0)
            .zero_to_100_kmh()
            .expect("hypercar reaches 100");
        let sc = sports_car()
            .longitudinal_trace(1.0e-3, 60.0)
            .zero_to_100_kmh()
            .expect("sports car reaches 100");
        assert!(
            hc < sc,
            "hypercar {hc}s should be quicker than sports car {sc}s"
        );
    }

    #[test]
    fn zero_dt_returns_a_single_initial_sample() {
        let trace = sports_car().longitudinal_trace(0.0, 10.0);
        assert_eq!(trace.points.len(), 1);
        assert_eq!(trace.points[0].time, 0.0);
        assert_eq!(trace.points[0].speed, 0.0);
    }
}
