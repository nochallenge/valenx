//! [`Scenario`] / [`ScenarioSuite`] — the *inputs* to a V&V run.
//!
//! A [`Scenario`] is a fully specified, reproducible test case for the
//! `valenx-sensors` autonomy harness: an initial kinematic [`VehicleState`], a
//! set of sensor configurations to attach, an analytic [`Scene`] environment,
//! and a finite [`CommandSeq`] (the control/event sequence to play out, with a
//! fixed time step). Running it yields a [`crate::Trace`]. A [`ScenarioSuite`]
//! is a named, non-empty collection of scenarios — e.g. one grid of swept
//! parameters, or a curated regression set.
//!
//! Each scenario also carries an optional set of named scalar **parameters**
//! ([`Scenario::params`]). These are *labels*, not behaviour — they record where
//! in a parameter space the scenario sits (e.g. `obstacle_x = 8.0`,
//! `approach_speed = 3.0`) so [`crate::coverage`] and [`crate::sweep`] can report
//! exactly which cells of a grid were exercised.

use std::collections::BTreeMap;

use valenx_sensors::{Command, Gps, Imu, Lidar, Scene, VehicleState};

use crate::error::VnvError;

/// How a scenario's control inputs are played out over time.
///
/// Either a single [`Command`] held for `steps` ticks of length `dt`, or an
/// explicit list of per-tick commands (an event sequence). Both are *finite* —
/// a V&V run always terminates.
#[derive(Debug, Clone, PartialEq)]
pub enum CommandSeq {
    /// Hold one constant command for `steps` ticks of length `dt` seconds.
    Constant {
        /// The command applied every tick.
        command: Command,
        /// Time step per tick (s, finite and `> 0`).
        dt: f64,
        /// Number of ticks (`> 0`).
        steps: usize,
    },
    /// Play an explicit, finite list of per-tick commands, each over `dt`
    /// seconds. The trace will have exactly `commands.len()` stepped frames.
    Explicit {
        /// One command per tick, in order (non-empty).
        commands: Vec<Command>,
        /// Time step per tick (s, finite and `> 0`).
        dt: f64,
    },
}

impl CommandSeq {
    /// The number of ticks this sequence will produce.
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            CommandSeq::Constant { steps, .. } => *steps,
            CommandSeq::Explicit { commands, .. } => commands.len(),
        }
    }

    /// Whether the sequence is empty (zero ticks). A valid sequence never is;
    /// [`Scenario::validate`] rejects an empty one.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The fixed per-tick time step.
    #[must_use]
    pub fn dt(&self) -> f64 {
        match self {
            CommandSeq::Constant { dt, .. } | CommandSeq::Explicit { dt, .. } => *dt,
        }
    }

    /// The total simulated duration (`len · dt` seconds).
    #[must_use]
    pub fn duration(&self) -> f64 {
        self.len() as f64 * self.dt()
    }

    /// Validate that the sequence is finite, non-empty, with a finite positive
    /// `dt` and only finite command components.
    fn validate(&self) -> Result<(), VnvError> {
        let dt = self.dt();
        if !(dt.is_finite() && dt > 0.0) {
            return Err(VnvError::InvalidConfig(format!(
                "command-sequence dt must be finite and > 0, got {dt}"
            )));
        }
        if self.is_empty() {
            return Err(VnvError::InvalidConfig(
                "command sequence must have at least one tick".into(),
            ));
        }
        let check = |c: &Command| {
            c.accel_body.iter().all(|x| x.is_finite())
                && c.angular_rate_body.iter().all(|x| x.is_finite())
        };
        let ok = match self {
            CommandSeq::Constant { command, .. } => check(command),
            CommandSeq::Explicit { commands, .. } => commands.iter().all(check),
        };
        if !ok {
            return Err(VnvError::NonFinite(
                "command sequence has a non-finite component".into(),
            ));
        }
        Ok(())
    }

    /// Iterate the per-tick commands in order (cloning for the `Constant` case).
    pub(crate) fn iter_commands(&self) -> Vec<Command> {
        match self {
            CommandSeq::Constant { command, steps, .. } => vec![*command; *steps],
            CommandSeq::Explicit { commands, .. } => commands.clone(),
        }
    }
}

/// Which sensors a scenario attaches to the harness, by configuration.
///
/// Each is optional; a scenario can run with any subset (including none, for a
/// pure ground-truth kinematic check). Sensors are *moved into* the harness when
/// the scenario runs, so this struct is the reusable recipe and
/// [`crate::run_scenario`] builds fresh sensor instances per run from it being
/// cloned — keeping every run independent and reproducible.
#[derive(Debug, Clone, Default)]
pub struct SensorSet {
    /// An optional LiDAR (already built, with its seed baked in).
    pub lidar: Option<Lidar>,
    /// An optional IMU.
    pub imu: Option<Imu>,
    /// An optional GPS.
    pub gps: Option<Gps>,
}

impl SensorSet {
    /// An empty sensor set (ground-truth-only run).
    #[must_use]
    pub fn none() -> Self {
        Self::default()
    }

    /// Attach a LiDAR (builder style).
    #[must_use]
    pub fn with_lidar(mut self, lidar: Lidar) -> Self {
        self.lidar = Some(lidar);
        self
    }

    /// Attach an IMU (builder style).
    #[must_use]
    pub fn with_imu(mut self, imu: Imu) -> Self {
        self.imu = Some(imu);
        self
    }

    /// Attach a GPS (builder style).
    #[must_use]
    pub fn with_gps(mut self, gps: Gps) -> Self {
        self.gps = Some(gps);
        self
    }

    /// Whether a LiDAR is configured (used by requirement/coverage checks that
    /// need to know what a trace *could* contain).
    #[must_use]
    pub fn has_lidar(&self) -> bool {
        self.lidar.is_some()
    }
}

/// A single, fully specified, reproducible V&V test case.
#[derive(Debug, Clone)]
pub struct Scenario {
    /// A short human-readable name (for reports).
    pub name: String,
    /// The vehicle's initial kinematic state.
    pub initial_state: VehicleState,
    /// The sensors to attach for this run.
    pub sensors: SensorSet,
    /// The analytic world the run takes place in.
    pub scene: Scene,
    /// The finite control/event sequence to play out.
    pub commands: CommandSeq,
    /// Named scalar parameters tagging *where* this scenario sits in a swept
    /// parameter space (labels for coverage/sweep, not behaviour). A
    /// `BTreeMap` so the ordering — and any derived coverage key — is
    /// deterministic.
    pub params: BTreeMap<String, f64>,
}

impl Scenario {
    /// Build a minimal scenario (no sensors, no parameters) from a name, an
    /// initial state, a scene, and a command sequence. Use the `with_*` /
    /// `set_param` builders to add sensors and parameter tags.
    #[must_use]
    pub fn new(
        name: impl Into<String>,
        initial_state: VehicleState,
        scene: Scene,
        commands: CommandSeq,
    ) -> Self {
        Self {
            name: name.into(),
            initial_state,
            sensors: SensorSet::none(),
            scene,
            commands,
            params: BTreeMap::new(),
        }
    }

    /// Attach a sensor set (builder style).
    #[must_use]
    pub fn with_sensors(mut self, sensors: SensorSet) -> Self {
        self.sensors = sensors;
        self
    }

    /// Tag a named scalar parameter (builder style). Overwrites any existing
    /// value for `key`.
    #[must_use]
    pub fn with_param(mut self, key: impl Into<String>, value: f64) -> Self {
        self.params.insert(key.into(), value);
        self
    }

    /// Tag a named scalar parameter in place.
    pub fn set_param(&mut self, key: impl Into<String>, value: f64) {
        self.params.insert(key.into(), value);
    }

    /// Validate the scenario: a finite non-empty command sequence and finite
    /// initial-state / parameter values. (Sensor configs were already validated
    /// at their own construction.)
    ///
    /// # Errors
    /// [`VnvError::InvalidConfig`] / [`VnvError::NonFinite`] as appropriate.
    pub fn validate(&self) -> Result<(), VnvError> {
        self.commands.validate()?;
        let s = &self.initial_state;
        let finite_state = s.position.iter().all(|x| x.is_finite())
            && s.velocity.iter().all(|x| x.is_finite())
            && s.angular_rate.iter().all(|x| x.is_finite());
        if !finite_state {
            return Err(VnvError::NonFinite(format!(
                "scenario '{}' has a non-finite initial state",
                self.name
            )));
        }
        for (k, v) in &self.params {
            if !v.is_finite() {
                return Err(VnvError::NonFinite(format!(
                    "scenario '{}' parameter '{k}' is non-finite ({v})",
                    self.name
                )));
            }
        }
        Ok(())
    }
}

/// A named, non-empty collection of [`Scenario`]s evaluated together — e.g. a
/// swept parameter grid or a curated regression set.
#[derive(Debug, Clone)]
pub struct ScenarioSuite {
    /// A short human-readable name (for reports).
    pub name: String,
    /// The scenarios in the suite (non-empty for a valid suite).
    pub scenarios: Vec<Scenario>,
}

impl ScenarioSuite {
    /// Build a suite from a name and a list of scenarios.
    ///
    /// This does *not* validate; call [`ScenarioSuite::validate`] (or rely on
    /// [`crate::sweep::run_suite`], which validates first) before running.
    #[must_use]
    pub fn new(name: impl Into<String>, scenarios: Vec<Scenario>) -> Self {
        Self {
            name: name.into(),
            scenarios,
        }
    }

    /// Number of scenarios in the suite.
    #[must_use]
    pub fn len(&self) -> usize {
        self.scenarios.len()
    }

    /// Whether the suite is empty (invalid — a suite must have ≥ 1 scenario).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.scenarios.is_empty()
    }

    /// Validate the suite: non-empty, and every scenario valid.
    ///
    /// # Errors
    /// [`VnvError::InvalidConfig`] if the suite is empty; otherwise the first
    /// scenario validation error.
    pub fn validate(&self) -> Result<(), VnvError> {
        if self.is_empty() {
            return Err(VnvError::InvalidConfig(format!(
                "scenario suite '{}' is empty (need ≥ 1 scenario)",
                self.name
            )));
        }
        for s in &self.scenarios {
            s.validate()?;
        }
        Ok(())
    }
}
