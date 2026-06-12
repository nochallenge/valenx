//! Time-course simulation driver — feature 13.
//!
//! [`TimeCourse`] wraps an integrator behind a COPASI-style "task"
//! interface: it integrates a [`Model`] from `t0` to `t_end`, samples
//! the trajectory onto a **uniform output grid** (so plots and tables
//! line up regardless of the adaptive step the integrator chose), and
//! applies **discrete events** — instantaneous state changes scheduled
//! at fixed times.
//!
//! Events here are *time-triggered* (the common case: "add 10 units of
//! inducer at t = 100"). Each [`Event`] names a species index, an
//! amount and an [`EventOp`] (set / add / scale). Triggering on a
//! state condition (a species crossing a threshold) is a documented
//! v1 omission — it needs root-finding inside the step, which the
//! fixed-grid driver does not do.
//!
//! The driver picks the integrator from an [`Integrator`] enum so a
//! caller can swap RK4 / RK45 / BDF without rewriting the task.

use crate::error::{Result, SysbioError};
use crate::model::Model;
use crate::ode::integrate::{integrate_rk4, Bdf, Rk45, Trajectory};
use crate::ode::OdeSystem;

/// How a discrete event modifies a species amount.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventOp {
    /// Overwrite the amount with the event value.
    Set,
    /// Add the event value to the amount.
    Add,
    /// Multiply the amount by the event value.
    Scale,
}

/// A time-triggered discrete event.
#[derive(Debug, Clone, PartialEq)]
pub struct Event {
    /// Simulation time at which the event fires.
    pub time: f64,
    /// Species index the event acts on.
    pub species: usize,
    /// The operand value (see [`EventOp`]).
    pub value: f64,
    /// The modification operation.
    pub op: EventOp,
}

impl Event {
    /// Apply this event in place to a state vector.
    fn apply(&self, state: &mut [f64]) {
        if let Some(slot) = state.get_mut(self.species) {
            match self.op {
                EventOp::Set => *slot = self.value,
                EventOp::Add => *slot += self.value,
                EventOp::Scale => *slot *= self.value,
            }
        }
    }
}

/// Which integrator the [`TimeCourse`] driver should use.
#[derive(Debug, Clone)]
pub enum Integrator {
    /// Fixed-step RK4 with the given step size.
    Rk4 {
        /// Internal integration step.
        dt: f64,
    },
    /// Adaptive Dormand-Prince RK45.
    Rk45(Rk45),
    /// Implicit BDF for stiff systems.
    Bdf(Bdf),
}

impl Default for Integrator {
    fn default() -> Self {
        Integrator::Rk45(Rk45::default())
    }
}

/// A configured time-course task.
#[derive(Debug, Clone)]
pub struct TimeCourse {
    /// Start time.
    pub t0: f64,
    /// End time.
    pub t_end: f64,
    /// Number of *intervals* on the uniform output grid (so
    /// `n_points + 1` samples including both endpoints).
    pub n_points: usize,
    /// Integrator to use.
    pub integrator: Integrator,
    /// Discrete events, applied in time order.
    pub events: Vec<Event>,
}

impl TimeCourse {
    /// A default task: `[0, t_end]`, 100 output intervals, adaptive
    /// RK45, no events.
    pub fn new(t_end: f64) -> Self {
        TimeCourse {
            t0: 0.0,
            t_end,
            n_points: 100,
            integrator: Integrator::default(),
            events: Vec::new(),
        }
    }

    /// Builder: add a discrete event.
    pub fn with_event(mut self, ev: Event) -> Self {
        self.events.push(ev);
        self
    }

    /// Run the time course on `model`, returning a trajectory sampled
    /// on the uniform output grid.
    ///
    /// Internally the `[t0, t_end]` span is split at every event time;
    /// each sub-interval is integrated independently and the event is
    /// applied to the carried-over state at the boundary. Finally the
    /// concatenated trajectory is resampled (linear interpolation)
    /// onto the requested uniform grid.
    pub fn run(&self, model: &Model) -> Result<Trajectory> {
        model.validate()?;
        if self.t_end <= self.t0 {
            return Err(SysbioError::invalid("t_end", "t_end must exceed t0"));
        }
        if self.n_points == 0 {
            return Err(SysbioError::invalid(
                "n_points",
                "need at least one interval",
            ));
        }
        let sys = OdeSystem::from_model(model);

        // Collect in-range event times, sorted & de-duplicated.
        let mut boundaries: Vec<f64> = self
            .events
            .iter()
            .map(|e| e.time)
            .filter(|&t| t > self.t0 && t < self.t_end)
            .collect();
        boundaries.sort_by(|a, b| a.partial_cmp(b).unwrap());
        boundaries.dedup();

        // Integrate piecewise. `raw` accumulates a fine trajectory.
        let mut raw = Trajectory {
            times: vec![self.t0],
            states: vec![model.initial_state()],
        };
        let mut seg_start = self.t0;
        let mut state = model.initial_state();
        let segment_ends: Vec<f64> = boundaries
            .iter()
            .copied()
            .chain(std::iter::once(self.t_end))
            .collect();

        for &seg_end in &segment_ends {
            if seg_end <= seg_start {
                continue;
            }
            let sub = self.integrate_segment(&sys, &state, seg_start, seg_end)?;
            // Append, skipping the duplicated first sample.
            for k in 1..sub.len() {
                raw.times.push(sub.times[k]);
                raw.states.push(sub.states[k].clone());
            }
            state = sub.final_state().unwrap().to_vec();
            // Apply every event scheduled exactly at this boundary.
            for ev in &self.events {
                if (ev.time - seg_end).abs() < 1e-12 {
                    ev.apply(&mut state);
                }
            }
            // Record the post-event state so the output grid sees the
            // jump.
            if (seg_end - self.t_end).abs() > 1e-12 {
                raw.times.push(seg_end);
                raw.states.push(state.clone());
            }
            seg_start = seg_end;
        }

        Ok(self.resample(&raw))
    }

    /// Integrate a single event-free sub-interval with the chosen
    /// integrator.
    fn integrate_segment(&self, sys: &OdeSystem, y0: &[f64], a: f64, b: f64) -> Result<Trajectory> {
        match &self.integrator {
            Integrator::Rk4 { dt } => integrate_rk4(sys, y0, a, b, *dt, 1),
            Integrator::Rk45(r) => r.integrate(sys, y0, a, b),
            Integrator::Bdf(bdf) => bdf.integrate(sys, y0, a, b),
        }
    }

    /// Linearly resample a fine trajectory onto the uniform output
    /// grid of `n_points + 1` samples.
    fn resample(&self, raw: &Trajectory) -> Trajectory {
        let dim = raw.states.first().map(|s| s.len()).unwrap_or(0);
        let mut out = Trajectory::default();
        let span = self.t_end - self.t0;
        for k in 0..=self.n_points {
            let t = self.t0 + span * (k as f64) / (self.n_points as f64);
            out.times.push(t);
            out.states.push(interp(raw, t, dim));
        }
        out
    }
}

/// Linear interpolation of a fine trajectory at time `t`.
fn interp(raw: &Trajectory, t: f64, dim: usize) -> Vec<f64> {
    if raw.is_empty() {
        return vec![0.0; dim];
    }
    if t <= raw.times[0] {
        return raw.states[0].clone();
    }
    let last = raw.len() - 1;
    if t >= raw.times[last] {
        return raw.states[last].clone();
    }
    // Binary search for the bracketing pair.
    let mut lo = 0;
    let mut hi = last;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if raw.times[mid] <= t {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let (t0, t1) = (raw.times[lo], raw.times[hi]);
    let w = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
    (0..dim)
        .map(|i| raw.states[lo][i] * (1.0 - w) + raw.states[hi][i] * w)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RateLaw, Reaction, Species};

    fn decay_model(k: f64, a0: f64) -> Model {
        let mut m = Model::new("decay");
        let a = m.add_species(Species::new("A", a0));
        m.add_reaction(Reaction {
            id: "d".into(),
            reactants: vec![(a, 1.0)],
            products: vec![],
            rate_law: RateLaw::MassAction {
                k,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        m
    }

    #[test]
    fn uniform_grid_has_requested_point_count() {
        let m = decay_model(1.0, 10.0);
        let tc = TimeCourse {
            n_points: 20,
            ..TimeCourse::new(5.0)
        };
        let traj = tc.run(&m).unwrap();
        assert_eq!(traj.len(), 21);
        assert!((traj.times[0] - 0.0).abs() < 1e-12);
        assert!((traj.times[20] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn decay_time_course_is_monotone_decreasing() {
        let m = decay_model(0.8, 16.0);
        let traj = TimeCourse::new(6.0).run(&m).unwrap();
        let a = traj.series(0);
        for w in a.windows(2) {
            assert!(w[1] <= w[0] + 1e-9, "not decreasing: {w:?}");
        }
        // Analytic endpoint.
        let expect = 16.0 * (-0.8_f64 * 6.0).exp();
        assert!((a[a.len() - 1] - expect).abs() < 1e-3);
    }

    #[test]
    fn add_event_injects_material() {
        // Pure decay; at t=2 add 100 units. The post-event sample
        // must jump well above the pre-event decayed value.
        let m = decay_model(1.0, 10.0);
        let tc = TimeCourse {
            n_points: 200,
            ..TimeCourse::new(4.0)
        }
        .with_event(Event {
            time: 2.0,
            species: 0,
            value: 100.0,
            op: EventOp::Add,
        });
        let traj = tc.run(&m).unwrap();
        let a = traj.series(0);
        let peak = a.iter().cloned().fold(0.0_f64, f64::max);
        assert!(peak > 100.0, "event did not inject material: peak {peak}");
    }

    #[test]
    fn set_event_overwrites_state() {
        let m = decay_model(1.0, 10.0);
        let tc = TimeCourse {
            n_points: 100,
            integrator: Integrator::Rk4 { dt: 0.01 },
            ..TimeCourse::new(4.0)
        }
        .with_event(Event {
            time: 2.0,
            species: 0,
            value: 50.0,
            op: EventOp::Set,
        });
        let traj = tc.run(&m).unwrap();
        // Just after t=2 the value should be ~50 (then decays).
        let idx = traj.times.iter().position(|&t| t >= 2.0).unwrap();
        assert!(traj.states[idx][0] > 40.0);
    }

    #[test]
    fn rejects_zero_points() {
        let m = decay_model(1.0, 1.0);
        let tc = TimeCourse {
            n_points: 0,
            ..TimeCourse::new(1.0)
        };
        assert!(tc.run(&m).is_err());
    }
}
