//! [`Trace`] — the *output* of a V&V run, and [`run_scenario`] that produces it.
//!
//! A [`Trace`] is the full record of one scenario run: the initial
//! [`VehicleState`] at `t = 0` plus the synchronised [`SensorFrame`] stream the
//! `valenx-sensors` harness emitted, one per tick. Every [`crate::Requirement`]
//! is a predicate over a `Trace`. The trace deliberately keeps the harness's
//! ground-truth [`VehicleState`] at every step (it lives inside each
//! [`SensorFrame`]), since most safety requirements are checked against the true
//! pose, not against a noisy estimate.

use valenx_sensors::{Harness, Scene, SensorFrame, VehicleState};

use crate::error::VnvError;
use crate::scenario::Scenario;

/// The full record of one scenario run.
#[derive(Debug, Clone, PartialEq)]
pub struct Trace {
    /// The name of the scenario this trace came from.
    pub scenario: String,
    /// The initial vehicle state at `t = 0` (before any step).
    pub initial_state: VehicleState,
    /// The analytic world the run took place in (a clone of the scenario's
    /// [`Scene`]). Carried on the trace so clearance/collision requirements are
    /// self-contained — they measure the ground-truth pose against this scene
    /// without any external context. `Scene` is `Clone + PartialEq`, so this
    /// keeps [`Trace`] comparable and self-describing.
    pub scene: Scene,
    /// The synchronised sensor frames, one per tick, in time order. Each frame
    /// carries the ground-truth [`VehicleState`] at its time (`frame.state`).
    pub frames: Vec<SensorFrame>,
}

impl Trace {
    /// The number of stepped frames in the trace.
    #[must_use]
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Whether the trace has no stepped frames. A trace from a valid scenario
    /// always has at least one; requirements treat an empty trace as a mismatch.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// The final simulated time (s), or `0.0` for an empty trace.
    #[must_use]
    pub fn final_time(&self) -> f64 {
        self.frames.last().map_or(0.0, |f| f.time)
    }

    /// The ground-truth vehicle state at every tick, in order (does **not**
    /// include the `t = 0` initial state — use [`Trace::initial_state`] for
    /// that, or [`Trace::states_with_initial`] for the full path).
    #[must_use]
    pub fn states(&self) -> Vec<VehicleState> {
        self.frames.iter().map(|f| f.state).collect()
    }

    /// The full ground-truth path: the initial state followed by the state at
    /// every tick.
    #[must_use]
    pub fn states_with_initial(&self) -> Vec<VehicleState> {
        let mut v = Vec::with_capacity(self.frames.len() + 1);
        v.push(self.initial_state);
        v.extend(self.frames.iter().map(|f| f.state));
        v
    }

    /// Whether any frame in the trace carries a LiDAR scan (used by
    /// LiDAR-inspecting requirements to detect a trace/requirement mismatch).
    #[must_use]
    pub fn has_lidar(&self) -> bool {
        self.frames.iter().any(|f| f.lidar.is_some())
    }
}

/// Drive the `valenx-sensors` [`Harness`] through a [`Scenario`] step-by-step
/// and collect the result into a [`Trace`].
///
/// The scenario is validated first (fail loud on a bad command sequence,
/// non-finite state, or non-finite parameter). A fresh harness is built from the
/// scenario's initial state, cloned scene, and cloned sensor set, so the run is
/// independent and — because every sensor owns a seeded PRNG — reproducible:
/// running the same scenario twice yields byte-identical traces.
///
/// # Errors
/// - The scenario's own [`Scenario::validate`] error, or
/// - a [`VnvError::Harness`] if the harness rejects a step (it should not, given
///   validation passed, but the error is surfaced rather than swallowed).
pub fn run_scenario(scenario: &Scenario) -> Result<Trace, VnvError> {
    scenario.validate()?;

    let mut harness = Harness::new(scenario.initial_state, scenario.scene.clone());
    if let Some(lidar) = scenario.sensors.lidar.clone() {
        harness = harness.with_lidar(lidar);
    }
    if let Some(imu) = scenario.sensors.imu.clone() {
        harness = harness.with_imu(imu);
    }
    if let Some(gps) = scenario.sensors.gps.clone() {
        harness = harness.with_gps(gps);
    }

    let dt = scenario.commands.dt();
    let commands = scenario.commands.iter_commands();
    let mut frames = Vec::with_capacity(commands.len());
    for command in &commands {
        frames.push(harness.step(command, dt)?);
    }

    Ok(Trace {
        scenario: scenario.name.clone(),
        initial_state: scenario.initial_state,
        scene: scenario.scene.clone(),
        frames,
    })
}
