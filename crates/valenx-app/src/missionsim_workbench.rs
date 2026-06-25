//! The right-side **Mission-Simulation Workbench** panel — a native front-end
//! over the in-house `valenx-mission-sim` crate (a general discrete-event /
//! agent constructive-simulation framework).
//!
//! This workbench is the **M2 track** of the valenx modeling-&-simulation
//! roadmap (AFSIM-class constructive mission/engagement simulation). It is a
//! *general* discrete-event / agent simulation front-end: the same engine
//! serves logistics, epidemiology, traffic flow, and policy wargaming as it
//! does a defensive force-on-force analysis. The user describes a small, tunable
//! demo scenario — a number of blue vs. red entities, their speeds / start
//! positions / headings, a sensor range, an engagement range, a
//! probability-of-kill (`Pk`), a stop time and a tick step — runs it to the stop
//! time, and reads back a **timeline**, the **final state**, and **outcome
//! metrics** (survivors per side, detection count, time-to-first-detection).
//! Alongside it runs a **Lanchester** aggregate sub-mode: two aggregate forces
//! `A`, `B` attrit each other under the century-old square-law ODE
//! `dA/dt = −b·B`, `dB/dt = −a·A`, and the workbench plots `A(t)` / `B(t)`.
//!
//! ## Dual-use boundary (abstract engagement only)
//!
//! This is a **general constructive simulation**: scheduling, geometry, and
//! aggregate / probabilistic outcomes. Engagement is kept deliberately
//! **abstract** — probability-of-kill (`Pk`) is an **INPUT** parameter and a
//! "hit" is just an abstract state change (an entity's `alive` flag); the
//! aggregate mode is the operations-research Lanchester ODE. **There is no
//! lethality model, no targeting / fire-control, and no kill-chain logic**
//! anywhere in this workbench or the `valenx-mission-sim` crate it drives. It is
//! force-on-force *bookkeeping* at the level of probabilities, geometry, and
//! aggregate counts — exactly the dual-use posture of the academic / think-tank
//! constructive sims the framework also serves.
//!
//! Mirrors the other workbenches (`uas_workbench`, `uq_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_missionsim_workbench`], toggled from the View menu
//! and openable by the agent bridge under the workbench id `"missionsim"`
//! (aliases `"mission"` / `"wargame"`; see [`crate::project_tabs::TabKind`] /
//! [`crate::agent_commands`]).
//!
//! Three painter views are drawn: (a) a **plan view** of the entity tracks over
//! the run (blue / red, sensor-range rings, detection markers), (b) a
//! **force-strength-vs-time** plot of the Lanchester `A(t)` / `B(t)` curves, and
//! (c) a **metrics readout**.
//!
//! Honesty: `valenx-mission-sim` is **research / educational grade** — the
//! movers are analytic kinematics (no dynamics, control, or terrain); detection
//! is geometric range-only line-of-sight (no occlusion, clutter, propagation, or
//! tracking filter); engagement is the abstract `Pk` / Lanchester abstraction
//! above, not a fidelity combat model; and the scenario loop is an explicit
//! fixed-cadence discretisation, so detection / engagement *times* resolve to
//! within one tick. Nothing here is accredited (VV&A). Every error from
//! `valenx-mission-sim` surfaces verbatim — the workbench never invents a number,
//! and degenerate parameters (a zero / negative stop time or tick, a `Pk`
//! outside `[0, 1]`, a negative range, a non-finite coordinate) show an in-panel
//! error, NOT a panic.

use eframe::egui;
use nalgebra::Vector3;
use valenx_mission_sim::{
    lanchester_square_step, monte_carlo, survivors_on, Entity, Event, ForceState, Mover,
    OutcomeStats, Scenario, ScenarioResult, Side,
};

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Editable scenario inputs: the two demo forces (counts / speeds / start
/// positions / headings), the shared sensor & engagement ranges and `Pk`, the
/// run parameters (stop time, tick), and a separate Lanchester aggregate group.
///
/// The demo scenario is a *deliberately simple, tunable* generator: `blue_count`
/// blue entities are placed on a vertical line near the origin, each carrying a
/// sensor (range `sensor_range_m`) and an engagement capability (range
/// `engagement_range_m`, kill probability `pk`); `red_count` red entities are
/// placed on a vertical line `red_standoff_m` away and march toward the blue
/// line at `red_speed_m_s` along the heading `red_heading_deg`. This is enough to
/// exercise detection (red crossing a blue sensor ring) and abstract engagement
/// (a red entity entering a blue engagement ring), while staying easy to reason
/// about.
#[derive(Clone, Debug)]
pub struct MissionSimParams {
    // -- Blue force --
    /// Number of blue entities (>= 0).
    pub blue_count: u32,
    /// Blue start x (m) — the blue line is at this x, spread in y.
    pub blue_x_m: f64,
    /// Blue line y-spacing between adjacent entities (m).
    pub blue_spacing_m: f64,
    /// Blue speed (m/s) along `blue_heading_deg`. `0` ⇒ static.
    pub blue_speed_m_s: f64,
    /// Blue heading (degrees, 0 = +x, 90 = +y).
    pub blue_heading_deg: f64,

    // -- Red force --
    /// Number of red entities (>= 0).
    pub red_count: u32,
    /// Red stand-off x (m) — the red line is at this x, spread in y.
    pub red_standoff_m: f64,
    /// Red line y-spacing between adjacent entities (m).
    pub red_spacing_m: f64,
    /// Red speed (m/s) along `red_heading_deg`. `0` ⇒ static.
    pub red_speed_m_s: f64,
    /// Red heading (degrees, 0 = +x, 90 = +y). Defaults to 180 (inbound −x).
    pub red_heading_deg: f64,

    // -- Shared sensing / engagement --
    /// Sensor detection range (m) carried by every blue entity. `0` ⇒ blind.
    pub sensor_range_m: f64,
    /// Engagement range (m) carried by every blue entity. `0` ⇒ cannot engage.
    pub engagement_range_m: f64,
    /// Abstract probability-of-kill input `Pk` in `[0, 1]` applied by blue
    /// entities. An INPUT parameter, not a lethality model.
    pub pk: f64,

    // -- Run parameters --
    /// Simulated stop time (s), finite and positive.
    pub stop_time_s: f64,
    /// Tick step (s), finite and positive; detection / engagement times resolve
    /// to within one tick.
    pub tick_dt_s: f64,
    /// Seed for the abstract engagement PRNG (determinism).
    pub seed: u64,

    // -- Lanchester aggregate sub-mode --
    /// Initial aggregate strength of force A.
    pub lanchester_a0: f64,
    /// Initial aggregate strength of force B.
    pub lanchester_b0: f64,
    /// Attrition-rate coefficient `a` (A's effectiveness against B).
    pub lanchester_a_rate: f64,
    /// Attrition-rate coefficient `b` (B's effectiveness against A).
    pub lanchester_b_rate: f64,
    /// Number of integration / plot sub-steps for the Lanchester curve (>= 1).
    pub lanchester_steps: usize,

    // -- Monte-Carlo engagement analysis --
    /// Number of Monte-Carlo replications of the constructive scenario (>= 1).
    /// Each replication is the SAME abstract scenario with a fresh seeded
    /// engagement stream; the spread of outcomes is pure statistics over the
    /// existing `Pk` draw (no new lethality / targeting).
    pub mc_runs: usize,
    /// Base seed for the Monte-Carlo ensemble. The per-run seeds are derived
    /// deterministically from this, so the same base seed replays an identical
    /// ensemble (and identical statistics).
    pub mc_seed: u64,
}

impl Default for MissionSimParams {
    fn default() -> Self {
        Self {
            blue_count: 3,
            blue_x_m: 0.0,
            blue_spacing_m: 150.0,
            blue_speed_m_s: 0.0,
            blue_heading_deg: 0.0,
            red_count: 3,
            red_standoff_m: 1500.0,
            red_spacing_m: 150.0,
            red_speed_m_s: 60.0,
            red_heading_deg: 180.0, // inbound, −x
            sensor_range_m: 800.0,
            engagement_range_m: 250.0,
            pk: 0.5,
            stop_time_s: 40.0,
            tick_dt_s: 0.25,
            seed: 0xC0FFEE,
            lanchester_a0: 100.0,
            lanchester_b0: 80.0,
            lanchester_a_rate: 0.05,
            lanchester_b_rate: 0.03,
            lanchester_steps: 200,
            mc_runs: 500,
            mc_seed: 0xC0FFEE,
        }
    }
}

impl MissionSimParams {
    /// Build the validated demo entities for the agent scenario, fail-loud.
    ///
    /// Blue entities sit on a vertical line at `blue_x_m`, centred in y; each
    /// carries the shared sensor / engagement ranges and `Pk`. Red entities sit
    /// on a vertical line at `red_standoff_m`, centred in y, and move along their
    /// heading at `red_speed_m_s`. Every [`Entity`] is built through
    /// `valenx-mission-sim`'s fail-loud constructor, so a bad range / `Pk` /
    /// non-finite coordinate surfaces that crate's error verbatim.
    fn build_entities(&self) -> Result<Vec<Entity>, String> {
        let mut entities = Vec::new();

        let line_offset = |count: u32, spacing: f64, i: u32| -> f64 {
            // Centre the line in y: i = 0..count maps to a symmetric spread.
            (i as f64 - (count as f64 - 1.0) / 2.0) * spacing
        };
        let heading_vec = |speed: f64, deg: f64| -> Vector3<f64> {
            let r = deg.to_radians();
            Vector3::new(speed * r.cos(), speed * r.sin(), 0.0)
        };

        // Blue force: sensors + engagement, optionally moving.
        let blue_vel = heading_vec(self.blue_speed_m_s, self.blue_heading_deg);
        let blue_mover = if self.blue_speed_m_s == 0.0 {
            Mover::Static
        } else {
            Mover::ConstantVelocity(blue_vel)
        };
        for i in 0..self.blue_count {
            let y = line_offset(self.blue_count, self.blue_spacing_m, i);
            let e = Entity::new(
                Vector3::new(self.blue_x_m, y, 0.0),
                Side::Blue,
                blue_mover.clone(),
                self.sensor_range_m,
                self.engagement_range_m,
                self.pk,
            )
            .map_err(|err| err.to_string())?;
            entities.push(e);
        }

        // Red force: inbound, no sensor / engagement of its own (a clean,
        // one-sided demo — blue detects & engages red).
        let red_vel = heading_vec(self.red_speed_m_s, self.red_heading_deg);
        let red_mover = if self.red_speed_m_s == 0.0 {
            Mover::Static
        } else {
            Mover::ConstantVelocity(red_vel)
        };
        for i in 0..self.red_count {
            let y = line_offset(self.red_count, self.red_spacing_m, i);
            let e = Entity::new(
                Vector3::new(self.red_standoff_m, y, 0.0),
                Side::Red,
                red_mover.clone(),
                0.0,
                0.0,
                0.0,
            )
            .map_err(|err| err.to_string())?;
            entities.push(e);
        }

        Ok(entities)
    }
}

// ---------------------------------------------------------------------------
// Simulation result
// ---------------------------------------------------------------------------

/// A sampled track for the plan view: a side plus the polyline its entity sweeps
/// from `t = 0` to the stop time, plus whether it ended alive.
#[derive(Clone, Debug)]
pub struct TrackPath {
    /// Which side this entity is on.
    pub side: Side,
    /// Sampled `(x, y)` positions over the run.
    pub points: Vec<[f64; 2]>,
    /// Sensor range of this entity (m); `0` ⇒ no ring drawn.
    pub sensor_range_m: f64,
    /// Whether the entity is alive at the stop time.
    pub alive_at_end: bool,
}

/// A detection marker for the plan view: where (the observer's position at the
/// detection time) and when.
#[derive(Clone, Copy, Debug)]
pub struct DetectionMarker {
    /// Observer position `(x, y)` at the detection time (m).
    pub at_xy: [f64; 2],
    /// Detection time (s).
    pub time_s: f64,
}

/// One Lanchester sample point: time and both force strengths.
#[derive(Clone, Copy, Debug)]
pub struct LanchesterPoint {
    /// Time (s).
    pub t: f64,
    /// Force A strength.
    pub a: f64,
    /// Force B strength.
    pub b: f64,
}

/// Cached mission-sim output for the painter + readouts.
#[derive(Default, Clone)]
pub struct MissionSimResult {
    /// Sampled tracks for the plan view.
    pub tracks: Vec<TrackPath>,
    /// Detection markers for the plan view.
    pub detections: Vec<DetectionMarker>,
    /// Blue survivors at the stop time.
    pub survivors_blue: usize,
    /// Red survivors at the stop time.
    pub survivors_red: usize,
    /// Total distinct first-detection count.
    pub detection_count: usize,
    /// Time-to-first-detection (s), if anything was detected.
    pub time_to_first_detection_s: Option<f64>,
    /// Number of engagement events recorded on the timeline.
    pub engagement_count: usize,
    /// Lanchester `A(t)` / `B(t)` curve samples.
    pub lanchester: Vec<LanchesterPoint>,
    /// Monte-Carlo engagement statistics, present once a Monte-Carlo ensemble has
    /// been run (the single-run pipeline leaves this `None`).
    pub mc_stats: Option<OutcomeStats>,
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the mission-sim workbench.
#[derive(Default)]
pub struct MissionSimWorkbenchState {
    /// User-editable parameters.
    pub params: MissionSimParams,
    /// Last successful result (populated after a successful run).
    pub result: Option<MissionSimResult>,
    /// Status / error line shown below the controls.
    pub status: String,
}

impl MissionSimWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`
    /// so an agent can discover the name space. Order follows the form.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "blue count",
            "blue x (m)",
            "blue spacing (m)",
            "blue speed (m/s)",
            "blue heading (deg)",
            "red count",
            "red standoff (m)",
            "red spacing (m)",
            "red speed (m/s)",
            "red heading (deg)",
            "sensor range (m)",
            "engagement range (m)",
            "Pk (0..1)",
            "stop time (s)",
            "tick dt (s)",
            "seed",
            "force A0",
            "force B0",
            "rate a",
            "rate b",
            "Lanchester steps",
            "Monte-Carlo runs N",
            "Monte-Carlo seed",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Captions match exactly what the form draws. Fail-loud
    /// on an unknown caption / wrong type (the bridge posts a `warn` note); no
    /// field is written on error and nothing panics. The count fields
    /// (`blue count` / `red count`, `Lanchester steps`) and `seed` read
    /// [`AgentValue::as_i64`]; every other caption is an `f64` drag value.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        let p = &mut self.params;
        match name {
            // -- Blue force --
            "blue count" => p.blue_count = parse_count(value, "blue count", 0, 64)? as u32,
            "blue x (m)" => p.blue_x_m = value.as_f64()?,
            "blue spacing (m)" => p.blue_spacing_m = value.as_f64()?,
            "blue speed (m/s)" => p.blue_speed_m_s = value.as_f64()?,
            "blue heading (deg)" => p.blue_heading_deg = value.as_f64()?,
            // -- Red force --
            "red count" => p.red_count = parse_count(value, "red count", 0, 64)? as u32,
            "red standoff (m)" => p.red_standoff_m = value.as_f64()?,
            "red spacing (m)" => p.red_spacing_m = value.as_f64()?,
            "red speed (m/s)" => p.red_speed_m_s = value.as_f64()?,
            "red heading (deg)" => p.red_heading_deg = value.as_f64()?,
            // -- Sensing & engagement --
            "sensor range (m)" => p.sensor_range_m = value.as_f64()?,
            "engagement range (m)" => p.engagement_range_m = value.as_f64()?,
            "Pk (0..1)" => p.pk = value.as_f64()?,
            // -- Run parameters --
            "stop time (s)" => p.stop_time_s = value.as_f64()?,
            "tick dt (s)" => p.tick_dt_s = value.as_f64()?,
            "seed" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("seed must be >= 0, got {n}"));
                }
                p.seed = n as u64;
            }
            // -- Lanchester aggregate --
            "force A0" => p.lanchester_a0 = value.as_f64()?,
            "force B0" => p.lanchester_b0 = value.as_f64()?,
            "rate a" => p.lanchester_a_rate = value.as_f64()?,
            "rate b" => p.lanchester_b_rate = value.as_f64()?,
            "Lanchester steps" => {
                p.lanchester_steps = parse_count(value, "Lanchester steps", 1, 5000)? as usize;
            }
            // -- Monte-Carlo engagement analysis --
            "Monte-Carlo runs N" => {
                p.mc_runs = parse_count(value, "Monte-Carlo runs N", 1, 1_000_000)? as usize;
            }
            "Monte-Carlo seed" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("Monte-Carlo seed must be >= 0, got {n}"));
                }
                p.mc_seed = n as u64;
            }
            other => return Err(format!("unknown mission-sim control: {other:?}")),
        }
        Ok(())
    }

    /// The current computed-result text for the agent `ReadReadout` bridge (see
    /// [`crate::agent_commands`]). This workbench keeps its result as a structured
    /// [`MissionSimResult`] and renders a one-line `status` summary (a `✔ …` line
    /// on success, a `⚠ …` line on error) — that same `status` string is returned
    /// here, with the **Monte-Carlo summary** (P(blue prevails) and the mean blue
    /// survivors ± 95% CI) appended once a Monte-Carlo ensemble has been run.
    /// `None` when nothing has run yet. Read-only — lets an agent read the answer
    /// back after driving a run, closing the live-driving loop.
    pub fn agent_readout(&self) -> Option<String> {
        if self.status.is_empty() {
            return None;
        }
        let mut out = self.status.clone();
        if let Some(stats) = self.result.as_ref().and_then(|r| r.mc_stats.as_ref()) {
            let bs = &stats.blue_survivors;
            out.push_str(&format!(
                " \u{00B7} MC[{} runs]: P(blue prevails) {:.3} \u{00B7} mean blue survivors {:.2} \
                 (95% CI {:.2}\u{2013}{:.2})",
                stats.runs, stats.p_blue_prevails, bs.mean, bs.ci95_lo, bs.ci95_hi
            ));
        }
        Some(out)
    }

    /// Run the full mission-sim pipeline: build + validate the demo entities, run
    /// the [`Scenario`] to the stop time (timeline + final state + metrics),
    /// sample the entity tracks for the plan view, and integrate the Lanchester
    /// aggregate curve.
    ///
    /// Every failure is returned as an `Err(String)` — no panics, no invented
    /// numbers. Degenerate inputs (a non-positive stop time / tick, a `Pk`
    /// outside `[0, 1]`, a negative range, a non-finite coordinate, a negative
    /// Lanchester rate / strength) surface `valenx-mission-sim`'s own error
    /// verbatim.
    pub fn run(&self) -> Result<MissionSimResult, String> {
        let p = &self.params;

        // --- Scenario (entities -> run -> timeline / final state / metrics) ---
        let entities = p.build_entities()?;
        let scenario = Scenario::new(entities, p.stop_time_s, p.tick_dt_s, p.seed)
            .map_err(|e| e.to_string())?;
        let res: ScenarioResult = scenario.run().map_err(|e| e.to_string())?;

        let mut tracks = sample_tracks(scenario.entities(), p.stop_time_s);
        // Fix up end-of-run liveness from the final state (index-aligned with the
        // input entities the tracks were sampled from).
        for (tr, fe) in tracks.iter_mut().zip(res.final_entities.iter()) {
            tr.alive_at_end = fe.alive;
        }
        let detections = detection_markers(scenario.entities(), &res);
        let engagement_count = res
            .timeline
            .iter()
            .filter(|e| matches!(e.event, Event::Engagement { .. }))
            .count();

        // --- Lanchester aggregate curve ---
        let lanchester = self.run_lanchester()?;

        // Cross-check survivor counts against the helper (cheap, also documents
        // the API in use).
        let survivors_blue = survivors_on(&res.final_entities, Side::Blue);
        let survivors_red = survivors_on(&res.final_entities, Side::Red);
        debug_assert_eq!(survivors_blue, res.metrics.survivors_blue);
        debug_assert_eq!(survivors_red, res.metrics.survivors_red);

        Ok(MissionSimResult {
            tracks,
            detections,
            survivors_blue: res.metrics.survivors_blue,
            survivors_red: res.metrics.survivors_red,
            detection_count: res.metrics.detection_count,
            time_to_first_detection_s: res.metrics.time_to_first_detection_s,
            engagement_count,
            lanchester,
            mc_stats: None,
        })
    }

    /// Run the **Monte-Carlo engagement analysis**: the single-run pipeline (for
    /// the plan view + Lanchester curve), then `mc_runs` replications of the
    /// *same* abstract scenario via `valenx-mission-sim`'s
    /// [`monte_carlo`], aggregated into an [`OutcomeStats`] (which side prevails,
    /// survivor spreads, exchange ratio, and a blue-survivor histogram).
    ///
    /// This adds **no** new lethality / targeting — it is pure statistics over the
    /// existing abstract `Pk` draw. Every failure is returned as `Err(String)`
    /// (no panics): a non-positive `mc_runs`, or any degenerate scenario
    /// parameter, surfaces `valenx-mission-sim`'s own error verbatim.
    pub fn run_monte_carlo(&self) -> Result<MissionSimResult, String> {
        let p = &self.params;
        if p.mc_runs == 0 {
            return Err("Monte-Carlo runs must be >= 1".to_string());
        }

        // Base single run populates the plan view + Lanchester curve as usual.
        let mut base = self.run()?;

        // Build the same scenario the single run used, re-seeded for the
        // ensemble, and sample it mc_runs times.
        let entities = p.build_entities()?;
        let scenario = Scenario::new(entities, p.stop_time_s, p.tick_dt_s, p.mc_seed)
            .map_err(|e| e.to_string())?;
        let stats = monte_carlo(&scenario, p.mc_runs).map_err(|e| e.to_string())?;

        base.mc_stats = Some(stats);
        Ok(base)
    }

    /// Integrate the Lanchester square-law ODE in `lanchester_steps` equal RK4
    /// sub-steps over the same stop time, sampling `A(t)` / `B(t)` at each step.
    ///
    /// Each step goes through `valenx-mission-sim`'s fail-loud
    /// [`lanchester_square_step`], so a negative rate / strength surfaces that
    /// crate's error verbatim rather than producing a bogus curve.
    fn run_lanchester(&self) -> Result<Vec<LanchesterPoint>, String> {
        let p = &self.params;
        if p.lanchester_steps == 0 {
            return Err("Lanchester steps must be >= 1".to_string());
        }
        let n = p.lanchester_steps;
        let dt = p.stop_time_s / n as f64;
        let mut state = ForceState {
            a: p.lanchester_a0,
            b: p.lanchester_b0,
        };
        let mut out = Vec::with_capacity(n + 1);
        out.push(LanchesterPoint {
            t: 0.0,
            a: state.a,
            b: state.b,
        });
        for i in 0..n {
            state = lanchester_square_step(state, p.lanchester_a_rate, p.lanchester_b_rate, dt)
                .map_err(|e| e.to_string())?;
            out.push(LanchesterPoint {
                t: dt * (i as f64 + 1.0),
                a: state.a,
                b: state.b,
            });
        }
        Ok(out)
    }
}

/// Sample each entity's track from `t = 0` to `stop` at a fixed number of points
/// for the plan view, reading the analytic mover off each input entity. Every
/// track is flagged `alive_at_end = true` here; the caller
/// ([`MissionSimWorkbenchState::run`]) fixes that up from the run's final state.
fn sample_tracks(entities: &[Entity], stop: f64) -> Vec<TrackPath> {
    const SAMPLES: usize = 48;
    entities
        .iter()
        .map(|e| {
            let points: Vec<[f64; 2]> = (0..=SAMPLES)
                .map(|i| {
                    let t = stop * i as f64 / SAMPLES as f64;
                    let p = e.position_at(t);
                    [p.x, p.y]
                })
                .collect();
            TrackPath {
                side: e.side,
                points,
                sensor_range_m: e.sensor_range_m,
                // Liveness is fixed up from the run result by the caller; default
                // alive here (geometry sampling only).
                alive_at_end: true,
            }
        })
        .collect()
}

/// Build plan-view detection markers from the run's timeline: for each
/// [`Event::Detection`] entry, place a marker at the observer's position at the
/// detection time.
fn detection_markers(entities: &[Entity], res: &ScenarioResult) -> Vec<DetectionMarker> {
    res.timeline
        .iter()
        .filter_map(|entry| {
            if let Event::Detection { observer, .. } = entry.event {
                entities.get(observer).map(|obs| {
                    let p = obs.position_at(entry.time);
                    DetectionMarker {
                        at_xy: [p.x, p.y],
                        time_s: entry.time,
                    }
                })
            } else {
                None
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the mission-sim workbench. A no-op unless toggled on via
/// View → Mission simulation.
///
/// Mirrors [`crate::uas_workbench::draw_uas_workbench`].
pub fn draw_missionsim_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_missionsim_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_missionsim_workbench",
        "Mission simulation (constructive)",
        missionsim_workbench_body,
    );
    if close {
        app.show_missionsim_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn missionsim_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    // The 21 controls + Monte-Carlo group + the viewport viz overflow the panel,
    // so wrap the whole body in a vertical scroll area (auto_shrink([false,false])
    // so it fills the panel width/height and shows a scrollbar when tall).
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| missionsim_workbench_body_inner(app, ui));
}

fn missionsim_workbench_body_inner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "General discrete-event / agent constructive simulation \u{00B7} \
             valenx-mission-sim  [research / educational \u{2014} analytic movers, range-only \
             detection; engagement is ABSTRACT (Pk input / Lanchester)]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let mut do_run = false;
    let mut do_run_mc = false;

    {
        let s = &mut app.missionsim;
        let p = &mut s.params;

        // --- Blue force -----------------------------------------------------
        ui.label(egui::RichText::new("Blue force").strong());
        egui::Grid::new("missionsim_blue_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                count_row(
                    ui,
                    "blue count",
                    &mut p.blue_count,
                    "Number of blue entities (each carries the shared sensor / engagement / Pk).",
                );
                drag_row(
                    ui,
                    "blue x (m)",
                    &mut p.blue_x_m,
                    10.0,
                    "Blue line x position; entities are spread in y.",
                );
                drag_row(
                    ui,
                    "blue spacing (m)",
                    &mut p.blue_spacing_m,
                    10.0,
                    "y-spacing between adjacent blue entities.",
                );
                drag_row(
                    ui,
                    "blue speed (m/s)",
                    &mut p.blue_speed_m_s,
                    1.0,
                    "Blue ground speed along the blue heading. 0 = static.",
                );
                drag_row(
                    ui,
                    "blue heading (deg)",
                    &mut p.blue_heading_deg,
                    5.0,
                    "Blue heading in degrees (0 = +x, 90 = +y).",
                );
            });

        // --- Red force ------------------------------------------------------
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Red force").strong());
        egui::Grid::new("missionsim_red_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                count_row(
                    ui,
                    "red count",
                    &mut p.red_count,
                    "Number of red entities (inbound; no sensor / engagement of their own in this demo).",
                );
                drag_row(ui, "red standoff (m)", &mut p.red_standoff_m, 25.0, "Red line x stand-off distance; entities are spread in y.");
                drag_row(ui, "red spacing (m)", &mut p.red_spacing_m, 10.0, "y-spacing between adjacent red entities.");
                drag_row(ui, "red speed (m/s)", &mut p.red_speed_m_s, 1.0, "Red ground speed along the red heading. 0 = static.");
                drag_row(ui, "red heading (deg)", &mut p.red_heading_deg, 5.0, "Red heading in degrees (0 = +x, 90 = +y; 180 = inbound -x).");
            });

        // --- Sensing & engagement ------------------------------------------
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Sensing & engagement (abstract)").strong());
        ui.label(
            egui::RichText::new(
                "Detection is range-only geometry; engagement is one ABSTRACT Pk draw (an \
                 input probability). No lethality, no targeting, no kill chain is modeled.",
            )
            .weak()
            .small(),
        );
        egui::Grid::new("missionsim_sensing_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                drag_row(
                    ui,
                    "sensor range (m)",
                    &mut p.sensor_range_m,
                    25.0,
                    "Blue sensor detection range. 0 = blind (no detections).",
                );
                drag_row(
                    ui,
                    "engagement range (m)",
                    &mut p.engagement_range_m,
                    10.0,
                    "Blue engagement range. 0 = cannot engage.",
                );
                drag_row(
                    ui,
                    "Pk (0..1)",
                    &mut p.pk,
                    0.02,
                    "Abstract probability-of-kill INPUT in [0, 1]. Not a lethality model.",
                );
            });

        // --- Run parameters -------------------------------------------------
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Run parameters").strong());
        egui::Grid::new("missionsim_run_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                drag_row(ui, "stop time (s)", &mut p.stop_time_s, 1.0, "Simulated stop time. Must be finite and positive.");
                drag_row(ui, "tick dt (s)", &mut p.tick_dt_s, 0.05, "Tick step; detection / engagement times resolve to within one tick. Must be > 0.");
                let lbl = ui.label("seed");
                ui.add(egui::DragValue::new(&mut p.seed).speed(1.0))
                    .labelled_by(lbl.id)
                    .on_hover_text("Seed for the abstract engagement PRNG. Same seed -> identical run.");
                ui.end_row();
            });

        // --- Lanchester aggregate sub-mode ----------------------------------
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Lanchester aggregate (square law)").strong());
        ui.label(
            egui::RichText::new(
                "Operations-research ODE dA/dt = -b*B, dB/dt = -a*A over two aggregate forces. \
                 a / b are abstract effectiveness coefficients, not lethality models.",
            )
            .weak()
            .small(),
        );
        egui::Grid::new("missionsim_lanchester_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                drag_row(
                    ui,
                    "force A0",
                    &mut p.lanchester_a0,
                    1.0,
                    "Initial aggregate strength of force A. Must be >= 0.",
                );
                drag_row(
                    ui,
                    "force B0",
                    &mut p.lanchester_b0,
                    1.0,
                    "Initial aggregate strength of force B. Must be >= 0.",
                );
                drag_row(
                    ui,
                    "rate a",
                    &mut p.lanchester_a_rate,
                    0.005,
                    "Attrition coefficient a (A's effectiveness against B). Must be >= 0.",
                );
                drag_row(
                    ui,
                    "rate b",
                    &mut p.lanchester_b_rate,
                    0.005,
                    "Attrition coefficient b (B's effectiveness against A). Must be >= 0.",
                );
                let lbl = ui.label("Lanchester steps");
                ui.add(
                    egui::DragValue::new(&mut p.lanchester_steps)
                        .speed(1)
                        .range(1..=5000),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Number of RK4 integration / plot sub-steps. Must be >= 1.");
                ui.end_row();
            });

        // --- Monte-Carlo engagement analysis -------------------------------
        ui.add_space(4.0);
        ui.label(egui::RichText::new("Monte-Carlo engagement analysis").strong());
        ui.label(
            egui::RichText::new(
                "Replicate the SAME abstract scenario N times with fresh seeded engagement \
                 streams and aggregate the spread of outcomes. Pure statistics over the Pk \
                 input \u{2014} no new lethality, targeting, or kill chain.",
            )
            .weak()
            .small(),
        );
        egui::Grid::new("missionsim_montecarlo_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("Monte-Carlo runs N");
                ui.add(
                    egui::DragValue::new(&mut p.mc_runs)
                        .speed(10)
                        .range(1..=1_000_000),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Number of replications N of the scenario. More runs tighten the \
                     confidence intervals (Monte-Carlo error is O(1/sqrt(N))). Must be >= 1.",
                );
                ui.end_row();

                let lbl = ui.label("Monte-Carlo seed");
                ui.add(egui::DragValue::new(&mut p.mc_seed).speed(1.0))
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Base seed for the ensemble. The per-run seeds derive from it, so the \
                         same base seed replays identical statistics.",
                    );
                ui.end_row();
            });

        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui
                .button(egui::RichText::new("Run").strong())
                .on_hover_text(
                    "Run the constructive scenario (timeline + final state + metrics) and \
                     integrate the Lanchester aggregate curve.",
                )
                .clicked()
            {
                do_run = true;
            }
            if ui
                .button(egui::RichText::new("Run Monte-Carlo").strong())
                .on_hover_text(
                    "Run the single scenario (plan view + Lanchester) AND replicate it N times \
                     to estimate P(blue prevails), survivor spreads, and the outcome histogram.",
                )
                .clicked()
            {
                do_run_mc = true;
            }
        });
    }

    // --- Execute (outside borrow) -------------------------------------------
    if do_run {
        run_and_store(app);
    }
    if do_run_mc {
        run_monte_carlo_and_store(app);
    }

    // --- Status line ---------------------------------------------------------
    let s = &app.missionsim;
    if !s.status.is_empty() {
        ui.add_space(6.0);
        let color = if s.status.starts_with('\u{26A0}') {
            egui::Color32::from_rgb(220, 120, 60)
        } else {
            egui::Color32::from_rgb(90, 180, 110)
        };
        ui.label(egui::RichText::new(&s.status).color(color).strong());
    }

    // --- Visualisation -------------------------------------------------------
    ui.add_space(6.0);
    ui.separator();
    draw_missionsim_viz(s, ui);
}

/// A labelled `DragValue` row in a 2-column grid: a caption cell, the drag value
/// (`labelled_by` the caption so it carries an accessible name), then `end_row`.
/// Mirrors the uas workbench's per-row caption pattern.
fn drag_row(ui: &mut egui::Ui, caption: &str, value: &mut f64, speed: f64, hint: &str) {
    let lbl = ui.label(caption);
    ui.add(egui::DragValue::new(value).speed(speed))
        .labelled_by(lbl.id)
        .on_hover_text(hint);
    ui.end_row();
}

/// Like [`drag_row`] but for a `u32` count (clamped to a sane range).
fn count_row(ui: &mut egui::Ui, caption: &str, value: &mut u32, hint: &str) {
    let lbl = ui.label(caption);
    ui.add(egui::DragValue::new(value).speed(1).range(0..=64))
        .labelled_by(lbl.id)
        .on_hover_text(hint);
    ui.end_row();
}

/// Read an [`crate::agent_commands::AgentValue`] as an integer count for a named
/// control and validate it against `[lo, hi]` (inclusive), the same bounds the
/// matching `DragValue` enforces in the UI. Fail-loud so an out-of-range count
/// becomes a `warn` note rather than silently clamping. Shared by the count /
/// step captions in [`MissionSimWorkbenchState::agent_set`].
fn parse_count(
    value: &crate::agent_commands::AgentValue,
    caption: &str,
    lo: i64,
    hi: i64,
) -> Result<i64, String> {
    let n = value.as_i64()?;
    if !(lo..=hi).contains(&n) {
        return Err(format!("{caption} must be in {lo}..={hi}, got {n}"));
    }
    Ok(n)
}

/// Run the pipeline and fold the result (or error) into the workbench status.
/// Factored out so the Run button (and tests) can share it.
pub(crate) fn run_and_store(app: &mut ValenxApp) {
    let s = &mut app.missionsim;
    match s.run() {
        Ok(res) => {
            let ttfd = res
                .time_to_first_detection_s
                .map(|t| format!("{t:.2}s"))
                .unwrap_or_else(|| "none".to_string());
            s.status = format!(
                "\u{2714} survivors B/R {}/{} \u{00B7} {} detections (first {}) \u{00B7} {} engagements",
                res.survivors_blue, res.survivors_red, res.detection_count, ttfd, res.engagement_count
            );
            s.result = Some(res);
        }
        Err(e) => {
            s.status = format!("\u{26A0} {e}");
            s.result = None;
        }
    }
}

/// Run the Monte-Carlo pipeline and fold the result (or error) into the workbench
/// status. The status carries the single-run summary plus the Monte-Carlo
/// headline (P(blue prevails) and mean blue survivors ± 95% CI). Factored out so
/// the Run Monte-Carlo button (and tests) can share it.
pub(crate) fn run_monte_carlo_and_store(app: &mut ValenxApp) {
    let s = &mut app.missionsim;
    match s.run_monte_carlo() {
        Ok(res) => {
            let mc = res
                .mc_stats
                .as_ref()
                .map(|st| {
                    let bs = &st.blue_survivors;
                    format!(
                        " \u{00B7} MC[{} runs] P(blue) {:.3} \u{00B7} mean blue {:.2} (95% CI {:.2}\u{2013}{:.2})",
                        st.runs, st.p_blue_prevails, bs.mean, bs.ci95_lo, bs.ci95_hi
                    )
                })
                .unwrap_or_default();
            s.status = format!(
                "\u{2714} survivors B/R {}/{} \u{00B7} {} detections \u{00B7} {} engagements{}",
                res.survivors_blue,
                res.survivors_red,
                res.detection_count,
                res.engagement_count,
                mc
            );
            s.result = Some(res);
        }
        Err(e) => {
            s.status = format!("\u{26A0} {e}");
            s.result = None;
        }
    }
}

// ---------------------------------------------------------------------------
// 2-D visualisation (plan view + Lanchester plot + metrics readout)
// ---------------------------------------------------------------------------

fn draw_missionsim_viz(s: &MissionSimWorkbenchState, ui: &mut egui::Ui) {
    let Some(res) = &s.result else {
        ui.label(
            egui::RichText::new(
                "press \"Run\" to simulate the scenario (plan view + metrics) and the \
                 Lanchester aggregate curve",
            )
            .weak(),
        );
        return;
    };

    draw_plan_view(res, ui);
    ui.add_space(8.0);
    draw_lanchester_plot(res, ui);
    ui.add_space(8.0);
    draw_metrics_readout(res, ui);
    if let Some(stats) = &res.mc_stats {
        ui.add_space(8.0);
        ui.separator();
        draw_monte_carlo_panel(stats, ui);
    }
}

/// View (a): a top-down **plan view** — every entity's track over the run (blue /
/// red), each blue sensor-range ring (drawn at the entity's start position), and
/// a marker at each first-detection. Pure geometry; nothing about engagement
/// targeting.
fn draw_plan_view(res: &MissionSimResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Plan view (entity tracks)").strong());
    ui.label(
        egui::RichText::new(
            "blue = blue side \u{00B7} red = red side \u{00B7} cyan rings = blue sensor range \
             \u{00B7} yellow \u{00D7} = first detection (geometry only)",
        )
        .weak()
        .small(),
    );

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(460.0), 260.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(10, 18, 28));

    if res.tracks.is_empty() {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "no entities",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    // World extent over all track points, all sensor rings, and all detections.
    let mut min = [f64::INFINITY, f64::INFINITY];
    let mut max = [f64::NEG_INFINITY, f64::NEG_INFINITY];
    let mut expand = |x: f64, y: f64| {
        if x.is_finite() && y.is_finite() {
            min[0] = min[0].min(x);
            min[1] = min[1].min(y);
            max[0] = max[0].max(x);
            max[1] = max[1].max(y);
        }
    };
    for tr in &res.tracks {
        for &[x, y] in &tr.points {
            expand(x, y);
        }
        if tr.sensor_range_m > 0.0 {
            if let Some(&[x, y]) = tr.points.first() {
                expand(x + tr.sensor_range_m, y + tr.sensor_range_m);
                expand(x - tr.sensor_range_m, y - tr.sensor_range_m);
            }
        }
    }
    for d in &res.detections {
        expand(d.at_xy[0], d.at_xy[1]);
    }
    if !(min[0].is_finite() && min[1].is_finite() && max[0].is_finite() && max[1].is_finite()) {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "non-finite geometry",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    // Square the world so rings stay circular.
    let span_x = (max[0] - min[0]).max(1.0);
    let span_y = (max[1] - min[1]).max(1.0);
    let span = span_x.max(span_y) * 1.12;
    let cx = (min[0] + max[0]) / 2.0;
    let cy = (min[1] + max[1]) / 2.0;
    let world_lo = [cx - span / 2.0, cy - span / 2.0];

    let margin = 14.0_f32;
    let inner = rect.shrink(margin);
    let side = inner.width().min(inner.height());
    let plot = egui::Rect::from_center_size(inner.center(), egui::vec2(side, side));
    let to_px = |x: f64, y: f64| -> egui::Pos2 {
        let fx = ((x - world_lo[0]) / span).clamp(0.0, 1.0) as f32;
        // World +y up; screen +y down -> invert.
        let fy = ((y - world_lo[1]) / span).clamp(0.0, 1.0) as f32;
        egui::pos2(plot.left() + fx * side, plot.bottom() - fy * side)
    };

    let blue = egui::Color32::from_rgb(90, 150, 240);
    let blue_dead = egui::Color32::from_rgb(70, 90, 130);
    let red = egui::Color32::from_rgb(220, 90, 80);
    let red_dead = egui::Color32::from_rgb(130, 70, 65);

    // Sensor rings first (under the tracks).
    for tr in &res.tracks {
        if tr.sensor_range_m > 0.0 {
            if let Some(&[x, y]) = tr.points.first() {
                let c = to_px(x, y);
                let ring_r = (tr.sensor_range_m / span) as f32 * side;
                painter.circle_stroke(
                    c,
                    ring_r,
                    egui::Stroke::new(1.0, egui::Color32::from_rgb(60, 150, 170)),
                );
            }
        }
    }

    // Tracks + endpoint markers.
    for tr in &res.tracks {
        let (col, dead) = match tr.side {
            Side::Blue => (blue, blue_dead),
            Side::Red => (red, red_dead),
        };
        let line_col = if tr.alive_at_end { col } else { dead };
        let pts: Vec<egui::Pos2> = tr.points.iter().map(|&[x, y]| to_px(x, y)).collect();
        for w in pts.windows(2) {
            painter.line_segment([w[0], w[1]], egui::Stroke::new(1.3, line_col));
        }
        // Start (hollow) and end (filled) markers.
        if let Some(first) = pts.first() {
            painter.circle_stroke(*first, 3.0, egui::Stroke::new(1.0, line_col));
        }
        if let Some(last) = pts.last() {
            if tr.alive_at_end {
                painter.circle_filled(*last, 3.5, col);
            } else {
                // A small X for a killed entity.
                let d = 3.5;
                painter.line_segment(
                    [*last + egui::vec2(-d, -d), *last + egui::vec2(d, d)],
                    egui::Stroke::new(1.6, dead),
                );
                painter.line_segment(
                    [*last + egui::vec2(-d, d), *last + egui::vec2(d, -d)],
                    egui::Stroke::new(1.6, dead),
                );
            }
        }
    }

    // Detection markers (yellow X).
    for det in &res.detections {
        let c = to_px(det.at_xy[0], det.at_xy[1]);
        let d = 4.0;
        let y = egui::Color32::from_rgb(235, 215, 90);
        painter.line_segment(
            [c + egui::vec2(-d, -d), c + egui::vec2(d, d)],
            egui::Stroke::new(1.6, y),
        );
        painter.line_segment(
            [c + egui::vec2(-d, d), c + egui::vec2(d, -d)],
            egui::Stroke::new(1.6, y),
        );
    }
}

/// View (b): the **Lanchester** force-strength-vs-time plot — `A(t)` and `B(t)`
/// over the run. The classic square-law trajectories: the stronger force pulls
/// away while the weaker collapses.
fn draw_lanchester_plot(res: &MissionSimResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Lanchester force strength vs time").strong());
    ui.label(
        egui::RichText::new(
            "green = force A(t) \u{00B7} orange = force B(t) \u{00B7} x = time (s)",
        )
        .weak()
        .small(),
    );

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(460.0), 180.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(14, 22, 34));

    if res.lanchester.len() < 2 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "too few Lanchester points",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    let t_lo = 0.0;
    let mut t_hi = f64::NEG_INFINITY;
    let mut y_hi = f64::NEG_INFINITY;
    for lp in &res.lanchester {
        if lp.t.is_finite() {
            t_hi = t_hi.max(lp.t);
        }
        if lp.a.is_finite() {
            y_hi = y_hi.max(lp.a);
        }
        if lp.b.is_finite() {
            y_hi = y_hi.max(lp.b);
        }
    }
    if !(t_hi.is_finite() && y_hi.is_finite()) || t_hi <= t_lo {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "non-finite Lanchester curve",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }
    let y_lo = 0.0;
    if (y_hi - y_lo).abs() < 1e-12 {
        y_hi += 1.0;
    }

    let margin = 26.0_f32;
    let inner = rect.shrink(margin);
    let to_px = |t: f64, y: f64| -> egui::Pos2 {
        let fx = ((t - t_lo) / (t_hi - t_lo)).clamp(0.0, 1.0) as f32;
        let fy = ((y - y_lo) / (y_hi - y_lo)).clamp(0.0, 1.0) as f32;
        egui::pos2(
            inner.left() + fx * inner.width(),
            inner.bottom() - fy * inner.height(),
        )
    };

    // Axes.
    painter.line_segment(
        [
            egui::pos2(inner.left(), inner.bottom()),
            egui::pos2(inner.right(), inner.bottom()),
        ],
        egui::Stroke::new(1.0, egui::Color32::from_gray(70)),
    );
    painter.line_segment(
        [
            egui::pos2(inner.left(), inner.top()),
            egui::pos2(inner.left(), inner.bottom()),
        ],
        egui::Stroke::new(1.0, egui::Color32::from_gray(70)),
    );

    // A(t) and B(t) polylines.
    let a_pts: Vec<egui::Pos2> = res.lanchester.iter().map(|lp| to_px(lp.t, lp.a)).collect();
    let b_pts: Vec<egui::Pos2> = res.lanchester.iter().map(|lp| to_px(lp.t, lp.b)).collect();
    for w in a_pts.windows(2) {
        painter.line_segment(
            [w[0], w[1]],
            egui::Stroke::new(1.8, egui::Color32::from_rgb(110, 210, 130)),
        );
    }
    for w in b_pts.windows(2) {
        painter.line_segment(
            [w[0], w[1]],
            egui::Stroke::new(1.8, egui::Color32::from_rgb(235, 165, 80)),
        );
    }

    // Axis labels.
    painter.text(
        egui::pos2(inner.center().x, rect.bottom() - 2.0),
        egui::Align2::CENTER_BOTTOM,
        "time (s)",
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(150),
    );
    painter.text(
        egui::pos2(rect.left() + 2.0, inner.center().y),
        egui::Align2::LEFT_CENTER,
        "strength",
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(150),
    );
}

/// View (c): a metrics readout — survivors per side, the detection count and
/// time-to-first-detection, the engagement count, and the Lanchester end state.
fn draw_metrics_readout(res: &MissionSimResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Outcome metrics").strong());
    egui::Grid::new("missionsim_metrics_grid")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(egui::RichText::new(v).monospace());
                ui.end_row();
            };
            row(ui, "survivors (blue)", format!("{}", res.survivors_blue));
            row(ui, "survivors (red)", format!("{}", res.survivors_red));
            row(ui, "detections", format!("{}", res.detection_count));
            row(
                ui,
                "time to first detection",
                res.time_to_first_detection_s
                    .map(|t| format!("{t:.2} s"))
                    .unwrap_or_else(|| "\u{2014}".to_string()),
            );
            row(ui, "engagements", format!("{}", res.engagement_count));
            if let (Some(first), Some(last)) = (res.lanchester.first(), res.lanchester.last()) {
                row(
                    ui,
                    "Lanchester A (start \u{2192} end)",
                    format!("{:.1} \u{2192} {:.1}", first.a, last.a),
                );
                row(
                    ui,
                    "Lanchester B (start \u{2192} end)",
                    format!("{:.1} \u{2192} {:.1}", first.b, last.b),
                );
            }
        });
}

/// View (d): the **Monte-Carlo engagement-analysis** panel, drawn in the same
/// viewport with the same painter idiom as the plan view / Lanchester plot. It
/// shows three things over the `N`-run ensemble:
///
/// * a **P(blue prevails) gauge** — a horizontal bar filled to the win fraction;
/// * the blue-survivor **histogram** — one bar per integer survivor count, bar
///   height proportional to the number of runs landing in that bin; and
/// * a **mean ± 95% CI readout** for blue / red survivors and the exchange ratio,
///   plus the outcome probabilities.
///
/// Pure statistics over the abstract `Pk` draw — nothing here is a lethality,
/// targeting, or kill-chain model.
fn draw_monte_carlo_panel(stats: &OutcomeStats, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Monte-Carlo engagement analysis").strong());
    ui.label(
        egui::RichText::new(format!(
            "{} replications \u{00B7} green gauge = P(blue prevails) \u{00B7} bars = blue-survivor \
             histogram (abstract Pk statistics only)",
            stats.runs
        ))
        .weak()
        .small(),
    );

    // --- P(blue prevails) gauge ------------------------------------------
    let available = ui.available_size();
    let gauge_w = available.x.min(460.0);
    let (grect, _) = ui.allocate_exact_size(egui::vec2(gauge_w, 26.0), egui::Sense::hover());
    let gp = ui.painter_at(grect);
    gp.rect_filled(grect, 3.0, egui::Color32::from_rgb(20, 28, 40));
    let frac = stats.p_blue_prevails.clamp(0.0, 1.0) as f32;
    let fill =
        egui::Rect::from_min_size(grect.min, egui::vec2(grect.width() * frac, grect.height()));
    gp.rect_filled(fill, 3.0, egui::Color32::from_rgb(90, 150, 240));
    gp.text(
        grect.center(),
        egui::Align2::CENTER_CENTER,
        format!("P(blue prevails) = {:.3}", stats.p_blue_prevails),
        egui::FontId::monospace(12.0),
        egui::Color32::WHITE,
    );

    // --- Blue-survivor histogram -----------------------------------------
    ui.add_space(6.0);
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(460.0), 160.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(12, 20, 30));

    let hist = &stats.blue_survivor_histogram;
    let max_count = hist.iter().copied().max().unwrap_or(0);
    if hist.is_empty() || max_count == 0 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "no histogram data",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
    } else {
        let margin = 26.0_f32;
        let inner = rect.shrink(margin);
        // Baseline axis.
        painter.line_segment(
            [
                egui::pos2(inner.left(), inner.bottom()),
                egui::pos2(inner.right(), inner.bottom()),
            ],
            egui::Stroke::new(1.0, egui::Color32::from_gray(70)),
        );
        let n_bins = hist.len();
        let slot_w = inner.width() / n_bins as f32;
        let bar_w = (slot_w * 0.7).max(1.0);
        let bar_col = egui::Color32::from_rgb(90, 150, 240);
        for (k, &count) in hist.iter().enumerate() {
            let h = (count as f32 / max_count as f32) * inner.height();
            let cx = inner.left() + slot_w * (k as f32 + 0.5);
            let bar = egui::Rect::from_min_max(
                egui::pos2(cx - bar_w / 2.0, inner.bottom() - h),
                egui::pos2(cx + bar_w / 2.0, inner.bottom()),
            );
            painter.rect_filled(bar, 0.0, bar_col);
            // Bin label (the survivor count) under the axis.
            painter.text(
                egui::pos2(cx, inner.bottom() + 2.0),
                egui::Align2::CENTER_TOP,
                format!("{k}"),
                egui::FontId::monospace(10.0),
                egui::Color32::from_gray(150),
            );
            // Count above non-empty bars.
            if count > 0 {
                painter.text(
                    egui::pos2(cx, inner.bottom() - h - 2.0),
                    egui::Align2::CENTER_BOTTOM,
                    format!("{count}"),
                    egui::FontId::monospace(9.0),
                    egui::Color32::from_gray(180),
                );
            }
        }
        // Axis captions.
        painter.text(
            egui::pos2(inner.center().x, rect.bottom() - 1.0),
            egui::Align2::CENTER_BOTTOM,
            "blue survivors",
            egui::FontId::monospace(10.0),
            egui::Color32::from_gray(150),
        );
        painter.text(
            egui::pos2(rect.left() + 2.0, inner.top()),
            egui::Align2::LEFT_TOP,
            "runs",
            egui::FontId::monospace(10.0),
            egui::Color32::from_gray(150),
        );
    }

    // --- Mean ± CI readout ------------------------------------------------
    ui.add_space(4.0);
    egui::Grid::new("missionsim_mc_grid")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(egui::RichText::new(v).monospace());
                ui.end_row();
            };
            row(
                ui,
                "P(blue prevails)",
                format!("{:.3}", stats.p_blue_prevails),
            );
            row(
                ui,
                "P(red prevails)",
                format!("{:.3}", stats.p_red_prevails),
            );
            row(ui, "P(draw)", format!("{:.3}", stats.p_draw));
            let fmt = |s: &valenx_mission_sim::SummaryStat| {
                format!(
                    "{:.3} \u{00B1} {:.3} (95% CI {:.3}\u{2013}{:.3})",
                    s.mean, s.std, s.ci95_lo, s.ci95_hi
                )
            };
            row(
                ui,
                "blue survivors (mean\u{00B1}sd)",
                fmt(&stats.blue_survivors),
            );
            row(
                ui,
                "red survivors (mean\u{00B1}sd)",
                fmt(&stats.red_survivors),
            );
            row(
                ui,
                "exchange ratio (mean\u{00B1}sd)",
                fmt(&stats.exchange_ratio),
            );
        });
}

// ---------------------------------------------------------------------------
// Tests (unit + headless_ui_tests, mirroring uas_workbench)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_run_succeeds_and_is_populated() {
        let s = MissionSimWorkbenchState::default();
        let res = s.run().expect("default mission-sim run should succeed");
        // 3 blue + 3 red entities -> 6 tracks.
        assert_eq!(res.tracks.len(), 6, "one track per entity");
        // Inbound red crossing 800 m blue sensors should be detected.
        assert!(res.detection_count >= 1, "inbound red should be detected");
        assert!(
            res.time_to_first_detection_s.is_some(),
            "a detection implies a first-detection time"
        );
        // Lanchester curve has steps + 1 samples.
        assert_eq!(res.lanchester.len(), s.params.lanchester_steps + 1);
    }

    #[test]
    fn lanchester_stronger_side_pulls_ahead() {
        // A0 > B0 with equal-ish rates -> A ends above B.
        let s = MissionSimWorkbenchState::default();
        let res = s.run().expect("run");
        let last = res.lanchester.last().expect("end point");
        assert!(last.a >= last.b, "stronger force A should end >= B");
    }

    #[test]
    fn run_is_deterministic_for_a_fixed_seed() {
        // Same seed -> identical metrics (the engagement draws are seeded).
        let mut s = MissionSimWorkbenchState::default();
        // Bring red right onto the blue engagement ring so 50/50 draws happen.
        s.params.pk = 0.5;
        let a = s.run().expect("run a");
        let b = s.run().expect("run b");
        assert_eq!(a.survivors_blue, b.survivors_blue);
        assert_eq!(a.survivors_red, b.survivors_red);
        assert_eq!(a.detection_count, b.detection_count);
        assert_eq!(a.engagement_count, b.engagement_count);
        assert_eq!(a.time_to_first_detection_s, b.time_to_first_detection_s);
    }

    #[test]
    fn pk_one_engagement_removes_red() {
        // Red marches into a blue engagement ring with Pk=1 -> red survivors drop.
        let mut s = MissionSimWorkbenchState::default();
        s.params.pk = 1.0;
        s.params.engagement_range_m = 400.0;
        s.params.stop_time_s = 60.0;
        let res = s.run().expect("run");
        assert!(
            res.survivors_red < s.params.red_count as usize,
            "Pk=1 inside engagement range should remove some red"
        );
    }

    #[test]
    fn empty_forces_run_clean() {
        let mut s = MissionSimWorkbenchState::default();
        s.params.blue_count = 0;
        s.params.red_count = 0;
        let res = s.run().expect("empty scenario should still run");
        assert!(res.tracks.is_empty());
        assert_eq!(res.survivors_blue, 0);
        assert_eq!(res.survivors_red, 0);
        assert_eq!(res.detection_count, 0);
    }

    // ---- degenerate-param tests — must return Err, NOT panic ----

    #[test]
    fn zero_stop_time_returns_err() {
        let mut s = MissionSimWorkbenchState::default();
        s.params.stop_time_s = 0.0;
        assert!(
            s.run().is_err(),
            "zero stop time must return Err, not panic"
        );
    }

    #[test]
    fn zero_tick_returns_err() {
        let mut s = MissionSimWorkbenchState::default();
        s.params.tick_dt_s = 0.0;
        assert!(s.run().is_err(), "zero tick must return Err, not panic");
    }

    #[test]
    fn pk_out_of_range_returns_err() {
        let mut s = MissionSimWorkbenchState::default();
        s.params.pk = 1.5;
        assert!(s.run().is_err(), "Pk > 1 must return Err, not panic");
        s.params.pk = -0.1;
        assert!(s.run().is_err(), "Pk < 0 must return Err, not panic");
    }

    #[test]
    fn negative_sensor_range_returns_err() {
        let mut s = MissionSimWorkbenchState::default();
        s.params.sensor_range_m = -10.0;
        assert!(
            s.run().is_err(),
            "negative sensor range must return Err, not panic"
        );
    }

    #[test]
    fn negative_lanchester_rate_returns_err() {
        let mut s = MissionSimWorkbenchState::default();
        s.params.lanchester_a_rate = -0.1;
        assert!(
            s.run().is_err(),
            "negative Lanchester rate must return Err, not panic"
        );
    }

    #[test]
    fn zero_lanchester_steps_returns_err() {
        let mut s = MissionSimWorkbenchState::default();
        s.params.lanchester_steps = 0;
        assert!(
            s.run().is_err(),
            "zero Lanchester steps must return Err, not panic"
        );
    }

    // ---- Monte-Carlo engagement analysis ----

    #[test]
    fn monte_carlo_populates_stats_and_single_run_views() {
        // run_monte_carlo must populate the MC stats AND keep the single-run plan
        // view + Lanchester curve (so the viewport still draws them).
        let s = MissionSimWorkbenchState::default();
        let res = s.run_monte_carlo().expect("monte-carlo run should succeed");
        let stats = res.mc_stats.expect("MC stats present");
        assert_eq!(stats.runs, s.params.mc_runs);
        assert_eq!(res.tracks.len(), 6, "single-run plan view still populated");
        assert_eq!(res.lanchester.len(), s.params.lanchester_steps + 1);
        // Histogram bins sum to N (one bin per possible blue-survivor count).
        assert_eq!(
            stats.blue_survivor_histogram.iter().sum::<usize>(),
            s.params.mc_runs,
            "histogram bin counts must sum to N"
        );
        // 3 blue entities -> survivors in 0..=3 -> 4 bins.
        assert_eq!(stats.blue_survivor_histogram.len(), 4);
        // Probabilities partition unity.
        let total = stats.p_blue_prevails + stats.p_red_prevails + stats.p_draw;
        assert!(
            (total - 1.0).abs() < 1e-12,
            "outcome probabilities sum to 1"
        );
    }

    #[test]
    fn monte_carlo_is_reproducible_for_a_fixed_seed() {
        // Same mc_seed -> bit-identical statistics.
        let s = MissionSimWorkbenchState::default();
        let a = s.run_monte_carlo().expect("run a").mc_stats.unwrap();
        let b = s.run_monte_carlo().expect("run b").mc_stats.unwrap();
        assert_eq!(
            a.p_blue_prevails, b.p_blue_prevails,
            "fixed seed reproduces P(blue prevails)"
        );
        assert_eq!(
            a.blue_survivors.mean, b.blue_survivors.mean,
            "fixed seed reproduces mean blue survivors"
        );
        assert_eq!(a.blue_survivor_histogram, b.blue_survivor_histogram);
    }

    #[test]
    fn monte_carlo_pk_zero_blue_never_loses_an_entity() {
        // Pk = 0: blue never removes a red AND red has no engagement of its own,
        // so every run ends 3 blue vs 3 red -> every run is a draw and all 500
        // runs land in the "3 blue survivors" histogram bin.
        let mut s = MissionSimWorkbenchState::default();
        s.params.pk = 0.0;
        let stats = s.run_monte_carlo().expect("run").mc_stats.unwrap();
        assert_eq!(stats.p_draw, 1.0, "Pk=0 symmetric survival -> all draws");
        assert_eq!(stats.p_blue_prevails, 0.0);
        assert!((stats.blue_survivors.mean - 3.0).abs() < 1e-12);
        assert_eq!(
            stats.blue_survivors.std, 0.0,
            "deterministic -> zero spread"
        );
        assert_eq!(stats.blue_survivor_histogram[3], s.params.mc_runs);
    }

    #[test]
    fn monte_carlo_zero_runs_returns_err() {
        let mut s = MissionSimWorkbenchState::default();
        s.params.mc_runs = 0;
        assert!(
            s.run_monte_carlo().is_err(),
            "zero Monte-Carlo runs must return Err, not panic"
        );
    }

    #[test]
    fn monte_carlo_degenerate_scenario_returns_err() {
        // A degenerate single-run parameter must fail the whole MC loud, not panic.
        let mut s = MissionSimWorkbenchState::default();
        s.params.stop_time_s = 0.0;
        assert!(s.run_monte_carlo().is_err(), "zero stop time -> Err");
        s.params.stop_time_s = 40.0;
        s.params.pk = 2.0;
        assert!(s.run_monte_carlo().is_err(), "Pk > 1 -> Err");
    }

    #[test]
    fn agent_readout_includes_mc_summary_after_monte_carlo() {
        // After a Monte-Carlo run the agent readout must carry the MC headline so
        // an agent can read P(blue prevails) + mean blue survivors back.
        let mut s = MissionSimWorkbenchState::default();
        let res = s.run_monte_carlo().expect("run");
        s.status = "\u{2714} ok".to_string();
        s.result = Some(res);
        let readout = s.agent_readout().expect("readout present");
        assert!(
            readout.contains("P(blue prevails)"),
            "MC readout must report P(blue prevails): {readout}"
        );
        assert!(
            readout.contains("mean blue survivors"),
            "MC readout must report mean blue survivors +/- CI: {readout}"
        );
    }

    #[test]
    fn agent_set_drives_monte_carlo_controls() {
        use crate::agent_commands::AgentValue;
        let mut s = MissionSimWorkbenchState::default();
        s.agent_set("Monte-Carlo runs N", &AgentValue::Int(250))
            .expect("set runs");
        assert_eq!(s.params.mc_runs, 250);
        s.agent_set("Monte-Carlo seed", &AgentValue::Int(123))
            .expect("set seed");
        assert_eq!(s.params.mc_seed, 123);
        // Fail-loud on a bad value (no panic, field unchanged).
        assert!(s
            .agent_set("Monte-Carlo runs N", &AgentValue::Int(0))
            .is_err());
        assert_eq!(
            s.params.mc_runs, 250,
            "rejected set must not change the field"
        );
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_missionsim_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_missionsim_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_missionsim_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_missionsim_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_missionsim_workbench = true;
        let res = app.missionsim.run().expect("run should succeed");
        app.missionsim.result = Some(res);
        app.missionsim.status = "\u{2714} test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_monte_carlo_panel_without_panic() {
        // Exercise the Monte-Carlo viewport panel (gauge + histogram bars + CI
        // readout) end-to-end through the accesskit draw path.
        let mut app = ValenxApp::default();
        app.show_missionsim_workbench = true;
        app.missionsim.params.mc_runs = 64;
        run_monte_carlo_and_store(&mut app);
        assert!(
            app.missionsim
                .result
                .as_ref()
                .is_some_and(|r| r.mc_stats.is_some()),
            "Monte-Carlo run should populate mc_stats"
        );
        assert!(
            app.missionsim.status.contains("P(blue)"),
            "status line should carry the MC headline: {}",
            app.missionsim.status
        );
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_missionsim_workbench = true;
        // Trigger an error state (zero stop time is fail-loud in run()).
        app.missionsim.params.stop_time_s = 0.0;
        let result = app.missionsim.run();
        app.missionsim.status = match result {
            Err(e) => format!("\u{26A0} {e}"),
            Ok(_) => "\u{26A0} simulated error for testing".to_string(),
        };
        app.missionsim.result = None;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by` its
        // caption (egui clears a DragValue's own Name), so an AI / screen reader
        // can find the control by caption text. Each `labelled_by` target must
        // RESOLVE to a real named caption node, not a dangling id.
        let mut app = ValenxApp::default();
        app.show_missionsim_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // Many numeric controls across blue / red / sensing / run / Lanchester /
        // Monte-Carlo; a conservative lower bound that all are present and named.
        assert!(
            spin_buttons.len() >= 20,
            "expected many numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );
        assert!(
            spin_buttons.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every DragValue's labelled_by must point at a named caption node"
        );

        // Representative captions present as named accessibility nodes.
        for caption in [
            "blue count",
            "red count",
            "red standoff (m)",
            "red speed (m/s)",
            "sensor range (m)",
            "engagement range (m)",
            "Pk (0..1)",
            "stop time (s)",
            "tick dt (s)",
            "seed",
            "force A0",
            "force B0",
            "rate a",
            "rate b",
            "Lanchester steps",
            "Monte-Carlo runs N",
            "Monte-Carlo seed",
        ] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }

        // The Run button must be a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| {
                n.role() == Role::Button && n.name().is_some_and(|s| s.contains("Run"))
            }),
            "the Run button must be a named, invokable node"
        );
    }

    #[test]
    fn degenerate_params_show_error_not_panic() {
        // When the stop time is 0 or the Pk is out of range the workbench must
        // surface the error in-panel, not panic.
        let mut state = MissionSimWorkbenchState::default();
        state.params.stop_time_s = 0.0;
        assert!(
            state.run().is_err(),
            "zero stop time must produce Err, not panic"
        );
        state.params.stop_time_s = 40.0;
        state.params.pk = 2.0;
        assert!(state.run().is_err(), "Pk > 1 must produce Err, not panic");
    }

    #[test]
    fn determinism_same_seed_same_metrics() {
        // A 50/50 engagement scenario: two runs at the same seed must yield
        // bit-identical metrics (proves the seeded engagement draws replay).
        let s = MissionSimWorkbenchState::default();
        let a = s.run().expect("run a");
        let b = s.run().expect("run b");
        assert_eq!(a.survivors_blue, b.survivors_blue);
        assert_eq!(a.survivors_red, b.survivors_red);
        assert_eq!(a.detection_count, b.detection_count);
        assert_eq!(a.engagement_count, b.engagement_count);
        assert_eq!(a.time_to_first_detection_s, b.time_to_first_detection_s);
    }

    #[test]
    fn agent_bridge_missionsim_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge uses for
        //   `OpenWorkbench { id: "missionsim" }` (and the aliases):
        //   1. TabKind::from_id("missionsim") -> Some(TabKind::MissionSim)
        //   2. set_workbench_flag(app, "missionsim", true) -> flag set
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup + aliases + case/whitespace tolerance.
        assert_eq!(
            TabKind::from_id("missionsim"),
            Some(TabKind::MissionSim),
            "\"missionsim\" must resolve to TabKind::MissionSim"
        );
        assert_eq!(TabKind::from_id("MISSIONSIM"), Some(TabKind::MissionSim));
        assert_eq!(
            TabKind::from_id("  missionsim  "),
            Some(TabKind::MissionSim)
        );
        assert_eq!(TabKind::from_id("mission"), Some(TabKind::MissionSim));
        assert_eq!(TabKind::from_id("wargame"), Some(TabKind::MissionSim));

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_missionsim_workbench);
        set_workbench_flag(&mut app, "missionsim", true);
        assert!(
            app.show_missionsim_workbench,
            "set_workbench_flag(\"missionsim\", true) must set show_missionsim_workbench"
        );
        set_workbench_flag(&mut app, "missionsim", false);
        assert!(!app.show_missionsim_workbench);
    }
}
