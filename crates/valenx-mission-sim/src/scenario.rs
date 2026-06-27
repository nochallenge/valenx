//! A **scenario** — wiring entities, movers, sensors, and engagements onto the
//! discrete-event scheduler, running it to a stop time, and returning a timeline
//! of events plus final state and outcome metrics.
//!
//! ## How the scenario advances time
//!
//! The scenario is driven by the general [`crate::Scheduler`]: it seeds a repeating
//! [`Event::Tick`] at a fixed cadence (`tick_dt`) up to the stop time, and at
//! each tick it evaluates the analytic geometry of every ordered entity pair:
//!
//! * **Detection** — for each live observer with a sensor and each live target
//!   of the *opposite* side, [`crate::sensor::detect`] (which reuses
//!   `valenx-uas`'s exact range-crossing geometry) decides whether the target is
//!   within sensor range *at this tick*. The first tick at which a given
//!   observer→target pair is in range is recorded as a [`Event::Detection`] and
//!   that pair's time-to-first-detection metric is captured.
//! * **Engagement** — for each live attacker with a positive engagement range
//!   and each live opposing target within that range, one **abstract**
//!   probability-of-kill draw ([`crate::engagement::resolve_pk`]) is resolved
//!   against the attacker's `pk` input. On a hit the target's `alive` flag is
//!   cleared and an [`Event::Engagement`] (with the outcome) is recorded.
//!
//! The tick cadence makes the whole run deterministic and fully general (the
//! same loop serves logistics, epidemiology, or traffic models with different
//! entity semantics). It is an explicit time discretisation: detection times are
//! reported to within one tick. The *per-pair geometry* inside a tick is exact
//! (closed-form), so refining `tick_dt` converges to the analytic crossing time;
//! the [`crate::sensor`] benchmark pins that underlying geometry exactly.
//!
//! Engagement draws all come from a single seeded [`crate::SplitMix64`], so a
//! given seed replays an identical timeline. The PRNG is **not** used for any
//! security purpose.
//!
//! ## Dual-use posture
//!
//! Everything here is infrastructure + analysis: scheduling, geometry, and
//! aggregate/probabilistic outcomes. Engagement is the abstract `Pk` draw only —
//! **no** targeting, fire-control, lethality, or kill-chain logic. The framework
//! is equally a logistics / epidemiology / traffic / policy-wargaming engine.

use crate::engagement::{resolve_pk, EngagementOutcome};
use crate::entity::{Entity, Side};
use crate::error::{require_positive, MissionError};
use crate::scenario_metrics::OutcomeMetrics;
use crate::scheduler::Scheduler;
use crate::sensor::{detect, range_between};
use crate::SplitMix64;

/// An event recorded on the scenario timeline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Event {
    /// A scheduler heartbeat that re-evaluates geometry. Carried internally; it
    /// is *not* surfaced on the returned timeline (only detections/engagements
    /// are), but it is the payload the scheduler pops.
    Tick,
    /// The simulation reached its stop time.
    Stop,
    /// `observer` (by entity index) first detected `target` (by entity index).
    Detection {
        /// Index of the observing entity.
        observer: usize,
        /// Index of the detected entity.
        target: usize,
    },
    /// `attacker` resolved an abstract engagement against `target`, with the
    /// probabilistic outcome.
    Engagement {
        /// Index of the attacking entity.
        attacker: usize,
        /// Index of the engaged entity.
        target: usize,
        /// The abstract Pk-draw outcome.
        outcome: EngagementOutcome,
    },
}

/// A timestamped timeline entry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimelineEntry {
    /// Simulated time of the event (s).
    pub time: f64,
    /// The event.
    pub event: Event,
}

/// The result of running a [`Scenario`].
#[derive(Debug, Clone, PartialEq)]
pub struct ScenarioResult {
    /// Ordered timeline of detections and engagements (ticks/stop excluded).
    pub timeline: Vec<TimelineEntry>,
    /// Final state of every entity (positions advanced to the stop time, with
    /// liveness reflecting any engagements).
    pub final_entities: Vec<Entity>,
    /// Aggregate outcome metrics.
    pub metrics: OutcomeMetrics,
}

/// A constructive-simulation scenario: a set of entities plus run parameters.
#[derive(Debug, Clone)]
pub struct Scenario {
    entities: Vec<Entity>,
    stop_time_s: f64,
    tick_dt_s: f64,
    seed: u64,
}

impl Scenario {
    /// Build a validated scenario.
    ///
    /// `stop_time_s` and `tick_dt_s` must be finite and positive (a zero/empty
    /// scenario is still legal — see [`Scenario::run`] — but if you ask it to
    /// *run* with a cadence it must be a real cadence). `seed` seeds the single
    /// engagement PRNG.
    ///
    /// # Errors
    ///
    /// [`MissionError::NonPositive`] if `stop_time_s` or `tick_dt_s` is not
    /// finite and positive.
    pub fn new(
        entities: Vec<Entity>,
        stop_time_s: f64,
        tick_dt_s: f64,
        seed: u64,
    ) -> Result<Self, MissionError> {
        require_positive("stop_time_s", stop_time_s)?;
        require_positive("tick_dt_s", tick_dt_s)?;
        Ok(Self {
            entities,
            stop_time_s,
            tick_dt_s,
            seed,
        })
    }

    /// The entities (read-only).
    #[must_use]
    pub fn entities(&self) -> &[Entity] {
        &self.entities
    }

    /// The engagement-PRNG seed this scenario runs with (read-only).
    ///
    /// Used by [`monte_carlo`](fn@crate::monte_carlo) to deterministically derive
    /// the per-run seeds of a reproducible Monte-Carlo ensemble.
    #[must_use]
    pub fn base_seed(&self) -> u64 {
        self.seed
    }

    /// A copy of this scenario with a different engagement-PRNG `seed`.
    ///
    /// Everything else (entities, stop time, tick) is identical, so the run is
    /// the **same** abstract scenario with a fresh stochastic stream — exactly
    /// what [`monte_carlo`](fn@crate::monte_carlo) needs to sample the engagement model.
    #[must_use]
    pub fn with_seed(&self, seed: u64) -> Self {
        Self {
            entities: self.entities.clone(),
            stop_time_s: self.stop_time_s,
            tick_dt_s: self.tick_dt_s,
            seed,
        }
    }

    /// Run the scenario to its stop time, returning the timeline + final state +
    /// metrics.
    ///
    /// An **empty** scenario (no entities) returns a clean empty result with no
    /// panic and zero survivors. Detection and engagement are evaluated every
    /// `tick_dt_s` up to and including `stop_time_s`.
    ///
    /// # Errors
    ///
    /// Propagates [`MissionError`] from the reused detection geometry or the Pk
    /// draw (both validate their inputs); in practice the validated entities and
    /// finite tick make these unreachable, but they are surfaced rather than
    /// unwrapped.
    pub fn run(&self) -> Result<ScenarioResult, MissionError> {
        let mut entities = self.entities.clone();
        let mut rng = SplitMix64::new(self.seed);
        let mut timeline: Vec<TimelineEntry> = Vec::new();

        // Track which ordered (observer, target) pairs have already been
        // detected, so each first-detection is recorded once.
        let n = entities.len();
        let mut detected = vec![false; n * n];
        let mut first_detect_time: Vec<Option<f64>> = vec![None; n * n];

        // Seed the scheduler with ticks up to the stop time, plus a final Stop.
        let mut sched: Scheduler<Event> = Scheduler::new();
        let mut t = 0.0;
        // A small epsilon guards the float accumulation so the last tick at/just
        // below the stop time is still enqueued.
        while t <= self.stop_time_s + self.tick_dt_s * 1e-9 {
            sched.schedule(t.min(self.stop_time_s), Event::Tick)?;
            if t >= self.stop_time_s {
                break;
            }
            t += self.tick_dt_s;
        }
        sched.schedule(self.stop_time_s, Event::Stop)?;

        // Drive the event loop.
        while let Some(ev) = sched.pop() {
            let now = ev.time;
            match ev.payload {
                Event::Stop => break,
                Event::Tick => {
                    self.evaluate_tick(
                        &mut entities,
                        &mut rng,
                        &mut timeline,
                        &mut detected,
                        &mut first_detect_time,
                        now,
                    )?;
                }
                // Detection / Engagement are only ever produced *into* the
                // timeline, never scheduled as queue payloads here.
                Event::Detection { .. } | Event::Engagement { .. } => {}
            }
        }

        // Advance every entity to the stop time for the reported final state: a
        // frozen snapshot at the stop position, carrying its end-of-run liveness.
        let final_entities: Vec<Entity> = entities
            .iter()
            .map(|e| {
                let mut snap = e.clone();
                snap.start_position = e.position_at(self.stop_time_s);
                snap.mover = crate::entity::Mover::Static;
                snap
            })
            .collect();

        let metrics = OutcomeMetrics::compute(&entities, &first_detect_time, n);

        Ok(ScenarioResult {
            timeline,
            final_entities,
            metrics,
        })
    }

    /// Evaluate detections and engagements for all live opposing pairs at `now`.
    #[allow(clippy::too_many_arguments)]
    fn evaluate_tick(
        &self,
        entities: &mut [Entity],
        rng: &mut SplitMix64,
        timeline: &mut Vec<TimelineEntry>,
        detected: &mut [bool],
        first_detect_time: &mut [Option<f64>],
        now: f64,
    ) -> Result<(), MissionError> {
        let n = entities.len();

        // --- Detection pass (read-only on positions; records first crossings).
        for obs in 0..n {
            for tgt in 0..n {
                if obs == tgt {
                    continue;
                }
                let idx = obs * n + tgt;
                if detected[idx] {
                    continue;
                }
                let (observer, target) = (&entities[obs], &entities[tgt]);
                if !observer.alive || !target.alive {
                    continue;
                }
                if observer.side == target.side {
                    continue; // sensors look at the opposing side here
                }
                if observer.sensor_range_m <= 0.0 {
                    continue;
                }
                // "Detected at this tick" = in range right now.
                let in_range = range_between(observer, target, now) <= observer.sensor_range_m;
                // Cross-check with the geometry helper (also validates inputs).
                let _ = detect(observer, target, now)?;
                if in_range {
                    detected[idx] = true;
                    first_detect_time[idx] = Some(now);
                    timeline.push(TimelineEntry {
                        time: now,
                        event: Event::Detection {
                            observer: obs,
                            target: tgt,
                        },
                    });
                }
            }
        }

        // --- Engagement pass (abstract Pk draws; may clear `alive`).
        // Deterministic order: attacker index, then target index.
        for atk in 0..n {
            if !entities[atk].alive || entities[atk].engagement_range_m <= 0.0 {
                continue;
            }
            for tgt in 0..n {
                if atk == tgt || !entities[tgt].alive {
                    continue;
                }
                if entities[atk].side == entities[tgt].side {
                    continue;
                }
                let within = range_between(&entities[atk], &entities[tgt], now)
                    <= entities[atk].engagement_range_m;
                if !within {
                    continue;
                }
                let outcome = resolve_pk(rng, entities[atk].pk)?;
                if outcome == EngagementOutcome::Hit {
                    entities[tgt].alive = false;
                }
                timeline.push(TimelineEntry {
                    time: now,
                    event: Event::Engagement {
                        attacker: atk,
                        target: tgt,
                        outcome,
                    },
                });
            }
        }
        Ok(())
    }
}

/// Count survivors on a given side among a slice of entities.
#[must_use]
pub fn survivors_on(entities: &[Entity], side: Side) -> usize {
    entities
        .iter()
        .filter(|e| e.alive && e.side == side)
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::entity::{Mover, Side};
    use nalgebra::Vector3;

    // ---- BENCHMARK PIN: empty scenario -> clean empty result, no panic ------
    #[test]
    fn empty_scenario_runs_clean() {
        let scn = Scenario::new(vec![], 10.0, 1.0, 0).unwrap();
        let res = scn.run().unwrap();
        assert!(res.timeline.is_empty());
        assert!(res.final_entities.is_empty());
        assert_eq!(res.metrics.survivors_blue, 0);
        assert_eq!(res.metrics.survivors_red, 0);
        assert_eq!(res.metrics.detection_count, 0);
        assert!(res.metrics.time_to_first_detection_s.is_none());
    }

    #[test]
    fn zero_or_negative_run_params_are_rejected() {
        assert!(Scenario::new(vec![], 0.0, 1.0, 0).is_err());
        assert!(Scenario::new(vec![], 10.0, 0.0, 0).is_err());
        assert!(Scenario::new(vec![], -1.0, 1.0, 0).is_err());
    }

    #[test]
    fn detection_recorded_near_the_geometric_crossing_time() {
        // Blue sensor at origin (range 500); red inbound from (1000,0,0) at
        // -100 m/s crosses x=500 at t=5 s. With a fine tick the recorded
        // first-detection time should be ~5 s.
        let blue =
            Entity::new(Vector3::zeros(), Side::Blue, Mover::Static, 500.0, 0.0, 0.0).unwrap();
        let red = Entity::new(
            Vector3::new(1000.0, 0.0, 0.0),
            Side::Red,
            Mover::ConstantVelocity(Vector3::new(-100.0, 0.0, 0.0)),
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        let scn = Scenario::new(vec![blue, red], 20.0, 0.01, 0).unwrap();
        let res = scn.run().unwrap();
        let ttfd = res.metrics.time_to_first_detection_s.unwrap();
        assert!((ttfd - 5.0).abs() <= 0.02, "ttfd = {ttfd}");
        assert!(res.metrics.detection_count >= 1);
        // Refining the tick converges toward the exact 5.0 (never overshooting:
        // detection is registered at the first tick that is in range, which is
        // at or after the true crossing).
        let scn2 = Scenario::new(scn.entities().to_vec(), 20.0, 0.001, 0).unwrap();
        let ttfd2 = scn2
            .run()
            .unwrap()
            .metrics
            .time_to_first_detection_s
            .unwrap();
        assert!(
            ttfd2 >= 5.0 - 1e-9,
            "tick-resolved detection is >= the true crossing"
        );
        assert!(
            (ttfd2 - 5.0).abs() <= (ttfd - 5.0).abs() + 1e-12,
            "finer tick is no worse"
        );
    }

    #[test]
    fn pk_one_engagement_kills_the_target() {
        // Co-located opposing pair; blue has engagement range 100 and Pk=1, so
        // it annihilates red on the first tick it is in range.
        let blue =
            Entity::new(Vector3::zeros(), Side::Blue, Mover::Static, 0.0, 100.0, 1.0).unwrap();
        let red = Entity::new(
            Vector3::new(10.0, 0.0, 0.0),
            Side::Red,
            Mover::Static,
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        let scn = Scenario::new(vec![blue, red], 5.0, 1.0, 0).unwrap();
        let res = scn.run().unwrap();
        assert_eq!(res.metrics.survivors_blue, 1);
        assert_eq!(res.metrics.survivors_red, 0, "Pk=1 must kill red");
        assert!(res.timeline.iter().any(|e| matches!(
            e.event,
            Event::Engagement {
                outcome: EngagementOutcome::Hit,
                ..
            }
        )));
    }

    #[test]
    fn pk_zero_engagement_never_kills() {
        let blue =
            Entity::new(Vector3::zeros(), Side::Blue, Mover::Static, 0.0, 100.0, 0.0).unwrap();
        let red = Entity::new(
            Vector3::new(10.0, 0.0, 0.0),
            Side::Red,
            Mover::Static,
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        let scn = Scenario::new(vec![blue, red], 5.0, 1.0, 0).unwrap();
        let res = scn.run().unwrap();
        assert_eq!(res.metrics.survivors_red, 1, "Pk=0 must never kill");
        assert!(res.timeline.iter().all(|e| !matches!(
            e.event,
            Event::Engagement {
                outcome: EngagementOutcome::Hit,
                ..
            }
        )));
    }

    #[test]
    fn out_of_range_pair_never_engages_or_detects() {
        // 5 km apart, sensor/engagement ranges far too small -> nothing happens.
        let blue = Entity::new(
            Vector3::zeros(),
            Side::Blue,
            Mover::Static,
            100.0,
            100.0,
            1.0,
        )
        .unwrap();
        let red = Entity::new(
            Vector3::new(5000.0, 0.0, 0.0),
            Side::Red,
            Mover::Static,
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        let scn = Scenario::new(vec![blue, red], 10.0, 1.0, 0).unwrap();
        let res = scn.run().unwrap();
        assert!(res.timeline.is_empty());
        assert_eq!(res.metrics.survivors_blue, 1);
        assert_eq!(res.metrics.survivors_red, 1);
    }

    #[test]
    fn run_is_deterministic_for_a_fixed_seed() {
        // A 50/50 engagement; two runs at the same seed must match exactly.
        let blue =
            Entity::new(Vector3::zeros(), Side::Blue, Mover::Static, 0.0, 100.0, 0.5).unwrap();
        let red = Entity::new(
            Vector3::new(10.0, 0.0, 0.0),
            Side::Red,
            Mover::Static,
            0.0,
            0.0,
            0.0,
        )
        .unwrap();
        let scn = Scenario::new(vec![blue, red], 50.0, 1.0, 0xC0FFEE).unwrap();
        let a = scn.run().unwrap();
        let b = scn.run().unwrap();
        assert_eq!(a.timeline, b.timeline);
        assert_eq!(a.metrics, b.metrics);
    }

    #[test]
    fn survivors_helper_counts_per_side() {
        let mut es = vec![
            Entity::new(Vector3::zeros(), Side::Blue, Mover::Static, 0.0, 0.0, 0.0).unwrap(),
            Entity::new(Vector3::zeros(), Side::Blue, Mover::Static, 0.0, 0.0, 0.0).unwrap(),
            Entity::new(Vector3::zeros(), Side::Red, Mover::Static, 0.0, 0.0, 0.0).unwrap(),
        ];
        es[1].alive = false;
        assert_eq!(survivors_on(&es, Side::Blue), 1);
        assert_eq!(survivors_on(&es, Side::Red), 1);
    }
}
