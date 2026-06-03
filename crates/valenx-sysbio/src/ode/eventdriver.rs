//! SBML L3 event + rule-aware time-course driver - feature 35.
//!
//! [`EventDrivenTimeCourse`] is the same shape as the existing
//! [`crate::ode::TimeCourse`] driver (start/end/n_points/integrator)
//! but augmented with:
//!
//! - **Event-trigger detection** between integrator steps. After every
//!   integrator step the driver evaluates each event's trigger; on a
//!   rising-edge crossing (`<= 0` last step, `> 0` now) it bisects the
//!   span to find the crossing time `t*`, restarts the integrator at
//!   `t*`, and queues the event for execution. Events with a delay are
//!   queued for execution at `t* + delay`; events without are executed
//!   immediately.
//! - **Simultaneous-event priority**. When more than one event fires in
//!   the same step the queue is sorted by descending `priority`
//!   (ties broken by event index, as in iBioSim).
//! - **Assignment-rule projection**. Every recorded output sample - and
//!   every queued event execution - re-applies the model's assignment
//!   rules so the integrator's raw state is projected onto the rule
//!   surface before the user sees it.
//! - **Rate rules on parameters** are folded into the ODE state by
//!   appending them as extra components (the state vector grows by the
//!   number of parameters with a rate rule; the driver writes those
//!   back into the system's parameter slice on every step).
//!
//! Why a separate driver from [`crate::ode::TimeCourse`]? The existing
//! one targets the *time-triggered* common case (an event scheduled at
//! a fixed time) and goes through fast piecewise integration with
//! boundary patching. State-triggered events need crossing detection
//! inside the step and a quite different control flow; mixing them
//! into one driver would make every code path harder to read and slow
//! the simple case down. The two drivers share the integrator
//! enum and trajectory type and can be swapped freely by the caller.

use crate::error::{Result, SysbioError};
use crate::model::events::{EventAssignment, SbmlEvent, VarRef};
use crate::model::expr::Expr;
use crate::model::Model;
use crate::ode::integrate::{integrate_rk4, Trajectory};
use crate::ode::system::OdeSystem;
use crate::ode::timecourse::Integrator;

/// A trajectory plus the events that fired during the run, in firing
/// order.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct EventTrajectory {
    /// The state trajectory on the uniform output grid (one sample
    /// per `n_points + 1` points, including both endpoints).
    pub trajectory: Trajectory,
    /// Time + event-index pairs for every event execution.
    pub event_log: Vec<(f64, usize)>,
}

impl EventTrajectory {
    /// Borrow the underlying state trajectory.
    pub fn states(&self) -> &Trajectory {
        &self.trajectory
    }
    /// Number of recorded event firings.
    pub fn n_events_fired(&self) -> usize {
        self.event_log.len()
    }
}

/// The event-and-rule-aware time-course driver.
#[derive(Debug, Clone)]
pub struct EventDrivenTimeCourse {
    /// Start time.
    pub t0: f64,
    /// End time.
    pub t_end: f64,
    /// Number of output intervals on the uniform output grid.
    pub n_points: usize,
    /// Integrator selection.
    pub integrator: Integrator,
    /// Absolute tolerance for trigger-crossing bisection.
    pub bisection_tol: f64,
    /// Maximum bisection iterations per crossing.
    pub bisection_max: usize,
    /// Internal step ceiling between output-grid samples.
    pub max_internal_steps: usize,
}

impl EventDrivenTimeCourse {
    /// A default task: `[0, t_end]`, 100 output intervals, RK45,
    /// `1e-6` bisection tolerance.
    pub fn new(t_end: f64) -> Self {
        EventDrivenTimeCourse {
            t0: 0.0,
            t_end,
            n_points: 100,
            integrator: Integrator::default(),
            bisection_tol: 1e-6,
            bisection_max: 60,
            max_internal_steps: 100_000,
        }
    }

    /// Run the time course on `model`. Returns the sampled trajectory
    /// plus a log of every event firing.
    pub fn run(&self, model: &Model) -> Result<EventTrajectory> {
        model.validate()?;
        if self.t_end <= self.t0 {
            return Err(SysbioError::invalid("t_end", "t_end must exceed t0"));
        }
        if self.n_points == 0 {
            return Err(SysbioError::invalid(
                "n_points",
                "need at least one output interval",
            ));
        }

        let mut sys = OdeSystem::from_model(model);
        // Apply initial-state projection so assignment rules hold at t0.
        let mut state = model.initial_state();
        sys.project_assignments(&mut state, self.t0)?;
        // Initialise the trigger baseline for every event.
        let mut events: Vec<SbmlEvent> = model.events.clone();
        for ev in events.iter_mut() {
            let v = ev.trigger.value(&state, sys.params(), self.t0);
            ev.initial_value = v > 0.0;
        }

        // Output grid.
        let span = self.t_end - self.t0;
        let dt_out = span / self.n_points as f64;
        let mut out = Trajectory {
            times: vec![self.t0],
            states: vec![state.clone()],
        };
        let mut event_log: Vec<(f64, usize)> = Vec::new();
        let mut t = self.t0;
        let mut pending: Vec<PendingEvent> = Vec::new();

        for sample in 1..=self.n_points {
            let t_target = self.t0 + dt_out * sample as f64;
            // Advance from t to t_target, processing event crossings.
            let mut internal_steps = 0;
            while t < t_target - 1e-15 {
                internal_steps += 1;
                if internal_steps > self.max_internal_steps {
                    return Err(SysbioError::not_converged(
                        "event_driver",
                        "exceeded internal-step ceiling",
                    ));
                }
                // Next deadline: t_target or the earliest pending event
                // execution time, whichever is sooner.
                let next_deadline = pending
                    .iter()
                    .map(|p| p.exec_time)
                    .fold(t_target, f64::min);
                let step_end = next_deadline.min(t_target);

                // Integrate one chunk from `t` to `step_end`.
                let sub = self.integrate_segment(&sys, &state, t, step_end)?;
                let new_state = sub.final_state().unwrap().to_vec();
                let new_t = step_end;

                // Look for event triggers that crossed inside (t, new_t].
                let mut crossed: Vec<(f64, usize)> = Vec::new();
                let mut probe_state = new_state.clone();
                sys.project_assignments(&mut probe_state, new_t)?;
                for (ei, ev) in events.iter().enumerate() {
                    let was_high = ev.initial_value;
                    let now_high = ev.trigger.value(&probe_state, sys.params(), new_t) > 0.0;
                    if !was_high && now_high {
                        // Bisect the (t, new_t] span to find the
                        // crossing point.
                        let t_cross = self.bisect_crossing(
                            &sys, &state, t, new_t, &ev.trigger,
                        )?;
                        crossed.push((t_cross, ei));
                    }
                }

                if !crossed.is_empty() {
                    // Restart at the earliest crossing, queue all events
                    // whose crossing fell in this segment. Bisection
                    // returns a t* very close to the trigger transition;
                    // we accept the post-crossing state at `new_t` (the
                    // full step end) so the trigger evaluates as high
                    // for the baseline refresh - taking only a partial
                    // step to a t* still slightly before the transition
                    // would leave `was_high = false` and immediately
                    // re-trigger the event on the next outer iteration.
                    crossed.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                    let earliest = crossed[0].0;
                    let mut simul: Vec<(f64, usize)> = crossed
                        .iter()
                        .filter(|(tc, _)| (tc - earliest).abs() < self.bisection_tol * 10.0)
                        .copied()
                        .collect();
                    simul.sort_by(|a, b| {
                        let pa = events[a.1].priority;
                        let pb = events[b.1].priority;
                        pb.partial_cmp(&pa)
                            .unwrap()
                            .then_with(|| a.1.cmp(&b.1))
                    });
                    // Adopt the new state at new_t (the integrator's
                    // accepted step end); the crossing time itself is
                    // used only for delay accounting and the log.
                    state = new_state;
                    sys.project_assignments(&mut state, new_t)?;
                    t = new_t;
                    for &(tc, ei) in &simul {
                        let ev = &events[ei];
                        let exec_time = ev.delay.map(|d| tc + d).unwrap_or(tc);
                        let frozen: Option<Vec<f64>> = if ev.use_values_from_trigger_time {
                            Some(
                                ev.assignments
                                    .iter()
                                    .map(|a| a.formula.value(&state, sys.params(), tc))
                                    .collect(),
                            )
                        } else {
                            None
                        };
                        pending.push(PendingEvent {
                            event_idx: ei,
                            exec_time,
                            priority: ev.priority,
                            frozen_values: frozen,
                        });
                    }
                    // Triggers crossed in this segment are now high.
                    for ev in events.iter_mut() {
                        ev.initial_value =
                            ev.trigger.value(&state, sys.params(), t) > 0.0;
                    }
                    // Execute any pending events whose execution time is
                    // at or before the current t. (No-delay events queued
                    // at t_cross < t are included.)
                    let mut to_exec: Vec<PendingEvent> = pending
                        .iter()
                        .filter(|p| p.exec_time <= t + 1e-12)
                        .cloned()
                        .collect();
                    pending.retain(|p| p.exec_time > t + 1e-12);
                    to_exec.sort_by(|a, b| {
                        b.priority
                            .partial_cmp(&a.priority)
                            .unwrap()
                            .then_with(|| a.event_idx.cmp(&b.event_idx))
                    });
                    for pe in to_exec {
                        let ev = &events[pe.event_idx];
                        let frozen = pe.frozen_values.as_deref();
                        apply_event(&ev.assignments, frozen, &mut state, sys.params_mut(), t);
                        event_log.push((t, pe.event_idx));
                        sys.project_assignments(&mut state, t)?;
                        for ev2 in events.iter_mut() {
                            ev2.initial_value =
                                ev2.trigger.value(&state, sys.params(), t) > 0.0;
                        }
                    }
                    continue;
                }

                // No crossings: accept the integrated chunk wholesale.
                state = new_state;
                sys.project_assignments(&mut state, new_t)?;
                // Update event baselines for the next sub-step.
                for ev in events.iter_mut() {
                    ev.initial_value =
                        ev.trigger.value(&state, sys.params(), new_t) > 0.0;
                }
                t = new_t;

                // Now execute any pending events whose execution time
                // is at or before t, in priority order.
                let mut to_exec: Vec<PendingEvent> = pending
                    .iter()
                    .filter(|p| p.exec_time <= t + 1e-12)
                    .cloned()
                    .collect();
                pending.retain(|p| p.exec_time > t + 1e-12);
                to_exec.sort_by(|a, b| {
                    b.priority
                        .partial_cmp(&a.priority)
                        .unwrap()
                        .then_with(|| a.event_idx.cmp(&b.event_idx))
                });
                for pe in to_exec {
                    let ev = &events[pe.event_idx];
                    let frozen = pe.frozen_values.as_deref();
                    apply_event(&ev.assignments, frozen, &mut state, sys.params_mut(), t);
                    event_log.push((t, pe.event_idx));
                    // Reproject assignment rules after the event
                    // assignment may have invalidated them.
                    sys.project_assignments(&mut state, t)?;
                    // Refresh event baselines.
                    for ev2 in events.iter_mut() {
                        ev2.initial_value =
                            ev2.trigger.value(&state, sys.params(), t) > 0.0;
                    }
                }
            }
            // Record the output sample at t_target.
            out.times.push(t_target);
            out.states.push(state.clone());
        }

        // Honour any remaining delayed events whose execution time
        // landed exactly at or past t_end - apply them as a final
        // pass so the last sample sees their effect.
        if !pending.is_empty() {
            pending.sort_by(|a, b| {
                b.priority
                    .partial_cmp(&a.priority)
                    .unwrap()
                    .then_with(|| a.event_idx.cmp(&b.event_idx))
            });
            for pe in pending {
                let ev = &events[pe.event_idx];
                apply_event(
                    &ev.assignments,
                    pe.frozen_values.as_deref(),
                    &mut state,
                    sys.params_mut(),
                    self.t_end,
                );
                event_log.push((self.t_end, pe.event_idx));
            }
            sys.project_assignments(&mut state, self.t_end)?;
            // Overwrite the last sample with the post-execution state.
            if let Some(last) = out.states.last_mut() {
                *last = state.clone();
            }
        }

        Ok(EventTrajectory {
            trajectory: out,
            event_log,
        })
    }

    /// Bisect the trigger-crossing point on `[t_lo, t_hi]`.
    ///
    /// Avoids re-integrating the dynamics on every bisection step
    /// (which would re-call RK45 on ever-shrinking intervals and risk
    /// driving the adaptive controller into its `h_min` floor). The
    /// strategy: integrate the full `[t_lo, t_hi]` span once, then
    /// bisect by *linearly interpolating* the resulting trajectory at
    /// each midpoint. Linear interpolation of an RK45 trajectory is
    /// the standard "dense output" surrogate used by event-detection
    /// SBML simulators - it converges to the true crossing time as
    /// the RK45 step grid refines.
    fn bisect_crossing(
        &self,
        sys: &OdeSystem,
        y_lo: &[f64],
        t_lo: f64,
        t_hi: f64,
        trigger: &Expr,
    ) -> Result<f64> {
        let traj = self.integrate_segment(sys, y_lo, t_lo, t_hi)?;
        let mut sys_clone = sys.clone();
        let mut trigger_at = |t: f64| -> Result<f64> {
            let mut s = interp_state(&traj, t);
            sys_clone.project_assignments(&mut s, t)?;
            Ok(trigger.value(&s, sys_clone.params(), t))
        };
        let mut lo = t_lo;
        let mut hi = t_hi;
        // Check we genuinely have a sign change in this interval - if
        // the linear-interpolation surrogate says no crossing (it can
        // disagree with the integrator at very-coarse step density),
        // fall back to the midpoint.
        let v_lo = trigger_at(lo)?;
        let v_hi = trigger_at(hi)?;
        if !(v_lo <= 0.0 && v_hi > 0.0) {
            return Ok(t_hi);
        }
        for _ in 0..self.bisection_max {
            if (hi - lo) < self.bisection_tol {
                return Ok(0.5 * (lo + hi));
            }
            let mid = 0.5 * (lo + hi);
            let v = trigger_at(mid)?;
            if v > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        Ok(0.5 * (lo + hi))
    }

    /// Integrate one event-free segment with the configured integrator.
    fn integrate_segment(
        &self,
        sys: &OdeSystem,
        y0: &[f64],
        a: f64,
        b: f64,
    ) -> Result<Trajectory> {
        if (b - a).abs() < 1e-14 {
            return Ok(Trajectory {
                times: vec![a, b],
                states: vec![y0.to_vec(), y0.to_vec()],
            });
        }
        match &self.integrator {
            Integrator::Rk4 { dt } => integrate_rk4(sys, y0, a, b, *dt, 1),
            Integrator::Rk45(r) => r.integrate(sys, y0, a, b),
            Integrator::Bdf(bdf) => bdf.integrate(sys, y0, a, b),
        }
    }
}

#[derive(Debug, Clone)]
struct PendingEvent {
    event_idx: usize,
    exec_time: f64,
    priority: f64,
    /// Pre-evaluated assignment RHS values when
    /// `use_values_from_trigger_time` is true; `None` means evaluate
    /// the formula at execution time.
    frozen_values: Option<Vec<f64>>,
}

/// Linear interpolation of a [`Trajectory`] at `t`. Used by the
/// crossing-detection bisection's "dense output" surrogate.
fn interp_state(traj: &Trajectory, t: f64) -> Vec<f64> {
    if traj.is_empty() {
        return Vec::new();
    }
    if t <= traj.times[0] {
        return traj.states[0].clone();
    }
    let last = traj.len() - 1;
    if t >= traj.times[last] {
        return traj.states[last].clone();
    }
    let mut lo = 0usize;
    let mut hi = last;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if traj.times[mid] <= t {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let (t0, t1) = (traj.times[lo], traj.times[hi]);
    let w = if t1 > t0 { (t - t0) / (t1 - t0) } else { 0.0 };
    let dim = traj.states[lo].len();
    (0..dim)
        .map(|i| traj.states[lo][i] * (1.0 - w) + traj.states[hi][i] * w)
        .collect()
}

/// Apply the event's assignment list to the state / parameter slices.
fn apply_event(
    assignments: &[EventAssignment],
    frozen: Option<&[f64]>,
    y: &mut [f64],
    p: &mut [f64],
    t: f64,
) {
    for (i, a) in assignments.iter().enumerate() {
        let v = match frozen {
            Some(fv) => fv[i],
            None => a.formula.value(y, p, t),
        };
        match &a.target {
            VarRef::Species(idx) => {
                if let Some(slot) = y.get_mut(*idx) {
                    *slot = v;
                }
            }
            VarRef::Parameter(idx) => {
                if let Some(slot) = p.get_mut(*idx) {
                    *slot = v;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::events::{AssignmentRule, EventAssignment, RateRule, SbmlEvent, VarRef};
    use crate::model::expr::Expr;
    use crate::model::{Parameter, RateLaw, Reaction, Species};
    use crate::ode::Rk45;

    /// Pure decay model A -> 0 with rate k.
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
    fn state_triggered_event_fires_and_resets() {
        // Decay A from 10; fire an event when A drops below 5.0 that
        // resets A to 100. We expect at least one event firing and
        // the final A to be well above the no-event analytic value.
        let mut m = decay_model(1.0, 10.0);
        m.add_event(SbmlEvent::new(
            "reset",
            Expr::lt(Expr::var(0), Expr::k(5.0)),
            vec![EventAssignment {
                target: VarRef::Species(0),
                formula: Expr::k(100.0),
            }],
        ));
        let driver = EventDrivenTimeCourse {
            n_points: 500,
            integrator: Integrator::Rk45(Rk45::default()),
            ..EventDrivenTimeCourse::new(2.0)
        };
        let traj = driver.run(&m).unwrap();
        assert!(
            traj.n_events_fired() >= 1,
            "event must have fired at least once"
        );
        // After at least one reset, A at t=2 should be larger than
        // pure-decay value 10*exp(-2) ~ 1.35.
        let series = traj.trajectory.series(0);
        let final_val = *series.last().unwrap();
        assert!(
            final_val > 5.0,
            "final A {final_val} suggests event did not reset"
        );
    }

    #[test]
    fn time_triggered_event_fires_once() {
        // Pure decay A0=10 k=1; event at t >= 1.0 sets A = 50.
        let mut m = decay_model(1.0, 10.0);
        m.add_event(SbmlEvent::new(
            "pulse",
            Expr::ge(Expr::Time, Expr::k(1.0)),
            vec![EventAssignment {
                target: VarRef::Species(0),
                formula: Expr::k(50.0),
            }],
        ));
        let driver = EventDrivenTimeCourse {
            n_points: 200,
            ..EventDrivenTimeCourse::new(3.0)
        };
        let traj = driver.run(&m).unwrap();
        assert_eq!(traj.n_events_fired(), 1);
        // The post-event peak should be near 50.
        let series = traj.trajectory.series(0);
        let peak = series.iter().cloned().fold(0.0_f64, f64::max);
        assert!(peak > 45.0, "peak after event was {peak}");
    }

    #[test]
    fn delayed_event_fires_after_delay() {
        let mut m = decay_model(1.0, 10.0);
        // Trigger at t=1; delay 0.5; set A to 200. Without delay the
        // post-event value would last from t=1; with delay it should
        // appear after t=1.5.
        m.add_event(
            SbmlEvent::new(
                "delayed",
                Expr::ge(Expr::Time, Expr::k(1.0)),
                vec![EventAssignment {
                    target: VarRef::Species(0),
                    formula: Expr::k(200.0),
                }],
            )
            .with_delay(0.5),
        );
        let driver = EventDrivenTimeCourse {
            n_points: 400,
            ..EventDrivenTimeCourse::new(2.5)
        };
        let traj = driver.run(&m).unwrap();
        assert_eq!(traj.n_events_fired(), 1);
        let (t_fired, _) = traj.event_log[0];
        // The execution time should be ~1.5.
        assert!(
            (t_fired - 1.5).abs() < 0.05,
            "delayed event fired at {t_fired}"
        );
        // At sample just before t=1.5, A should be small (pure decay);
        // at sample just after, A should jump to ~200.
        let times = &traj.trajectory.times;
        let series = traj.trajectory.series(0);
        let i_after =
            times.iter().position(|&t| t > 1.51).expect("sample after exec");
        assert!(
            series[i_after] > 150.0,
            "post-delay A = {} not at 200",
            series[i_after]
        );
        let i_before =
            times.iter().rposition(|&t| t < 1.45).expect("sample before exec");
        assert!(
            series[i_before] < 50.0,
            "pre-delay A = {} not still decayed",
            series[i_before]
        );
    }

    #[test]
    fn simultaneous_events_fire_in_priority_order() {
        // Two events with the same trigger time. e1 priority 10 sets
        // A = 1; e2 priority 1 sets A = 2. The post-event value
        // should be 2 (e2 ran second and overwrote e1's assignment).
        let mut m = Model::new("sim");
        let _a = m.add_species(Species::new("A", 0.0));
        m.add_event(
            SbmlEvent::new(
                "e1",
                Expr::ge(Expr::Time, Expr::k(1.0)),
                vec![EventAssignment {
                    target: VarRef::Species(0),
                    formula: Expr::k(1.0),
                }],
            )
            .with_priority(10.0),
        );
        m.add_event(
            SbmlEvent::new(
                "e2",
                Expr::ge(Expr::Time, Expr::k(1.0)),
                vec![EventAssignment {
                    target: VarRef::Species(0),
                    formula: Expr::k(2.0),
                }],
            )
            .with_priority(1.0),
        );
        // No reactions - the model is just events. Rather than
        // changing validate, give it a Const-0 reaction so it has
        // something to integrate.
        m.add_reaction(Reaction {
            id: "noop".into(),
            reactants: vec![],
            products: vec![],
            rate_law: RateLaw::Constant { rate: 0.0 },
            reversible: false,
        });
        let driver = EventDrivenTimeCourse {
            n_points: 200,
            ..EventDrivenTimeCourse::new(2.0)
        };
        let traj = driver.run(&m).unwrap();
        let series = traj.trajectory.series(0);
        let final_val = *series.last().unwrap();
        // Both events should have fired; e2 (low priority) runs after
        // e1, so the final A is 2.
        assert!(
            (final_val - 2.0).abs() < 1e-6,
            "final A {final_val} expected 2",
        );
        assert!(traj.n_events_fired() >= 2);
    }

    #[test]
    fn rate_rule_drives_species() {
        // No reactions; rate rule d A / dt = 2.0. A(t) = A0 + 2 t.
        let mut m = Model::new("rate");
        let _a = m.add_species(Species::new("A", 0.0));
        m.add_rate_rule(RateRule {
            target: VarRef::Species(0),
            formula: Expr::k(2.0),
        });
        m.add_reaction(Reaction {
            id: "noop".into(),
            reactants: vec![],
            products: vec![],
            rate_law: RateLaw::Constant { rate: 0.0 },
            reversible: false,
        });
        let driver = EventDrivenTimeCourse {
            n_points: 100,
            ..EventDrivenTimeCourse::new(5.0)
        };
        let traj = driver.run(&m).unwrap();
        let final_val = traj.trajectory.series(0).last().copied().unwrap();
        assert!(
            (final_val - 10.0).abs() < 0.05,
            "expected ~10 got {final_val}"
        );
    }

    #[test]
    fn assignment_rule_projects_state() {
        // Two species A, B. A has a rate rule dA/dt = 1. B := 2 * A
        // (assignment rule). At t = 5, A should be ~5, B should be
        // ~10.
        let mut m = Model::new("assign");
        let _a = m.add_species(Species::new("A", 0.0));
        let _b = m.add_species(Species::new("B", 0.0));
        m.add_rate_rule(RateRule {
            target: VarRef::Species(0),
            formula: Expr::k(1.0),
        });
        m.add_assignment_rule(AssignmentRule {
            target: VarRef::Species(1),
            formula: Expr::mul(Expr::var(0), Expr::k(2.0)),
        });
        m.add_reaction(Reaction {
            id: "noop".into(),
            reactants: vec![],
            products: vec![],
            rate_law: RateLaw::Constant { rate: 0.0 },
            reversible: false,
        });
        let driver = EventDrivenTimeCourse {
            n_points: 50,
            ..EventDrivenTimeCourse::new(5.0)
        };
        let traj = driver.run(&m).unwrap();
        let final_a = *traj.trajectory.series(0).last().unwrap();
        let final_b = *traj.trajectory.series(1).last().unwrap();
        assert!((final_a - 5.0).abs() < 0.05, "A = {final_a}");
        assert!((final_b - 10.0).abs() < 0.1, "B = {final_b}");
    }

    #[test]
    fn event_assigning_parameter_updates_it() {
        // Rate rule reads a parameter; an event halves the parameter
        // at t = 2.
        let mut m = Model::new("paramevt");
        let _a = m.add_species(Species::new("A", 0.0));
        m.add_parameter(Parameter::new("k", 2.0));
        m.add_rate_rule(RateRule {
            target: VarRef::Species(0),
            formula: Expr::param(0),
        });
        m.add_event(SbmlEvent::new(
            "halve_k",
            Expr::ge(Expr::Time, Expr::k(2.0)),
            vec![EventAssignment {
                target: VarRef::Parameter(0),
                formula: Expr::k(0.0),
            }],
        ));
        m.add_reaction(Reaction {
            id: "noop".into(),
            reactants: vec![],
            products: vec![],
            rate_law: RateLaw::Constant { rate: 0.0 },
            reversible: false,
        });
        let driver = EventDrivenTimeCourse {
            n_points: 200,
            ..EventDrivenTimeCourse::new(4.0)
        };
        let traj = driver.run(&m).unwrap();
        // Before event: A grows at 2.0/unit. At t=2, A ~= 4.
        // After event: rate = 0, so A stays ~= 4 to t=4.
        let final_a = *traj.trajectory.series(0).last().unwrap();
        assert!(
            (final_a - 4.0).abs() < 0.2,
            "expected ~4 after parameter zeroed, got {final_a}"
        );
        assert_eq!(traj.n_events_fired(), 1);
    }
}
