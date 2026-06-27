//! # valenx-autonomy-vnv — a V&V framework for autonomous systems
//!
//! An **in-house verification & validation (V&V) methodology** layered on the
//! [`valenx_sensors`] autonomy harness. It gives you the moving parts of a
//! test-&-evaluation campaign for an autonomous vehicle, as plain, deterministic
//! Rust:
//!
//! * **[`Scenario`] / [`ScenarioSuite`]** — a reproducible test case (an initial
//!   kinematic [`VehicleState`], a [`SensorSet`] to attach, an analytic [`Scene`]
//!   environment, and a finite [`CommandSeq`]), and a named collection of them.
//! * **[`Requirement`] / [`RequirementSet`]** — a safety/performance property as
//!   a predicate over a run, returning **pass/fail *and* a signed margin**
//!   ([`MinClearance`](Requirement::MinClearance),
//!   [`DetectByTime`](Requirement::DetectByTime),
//!   [`StayInBounds`](Requirement::StayInBounds),
//!   [`NoCollision`](Requirement::NoCollision)).
//! * **[`run_scenario`] → [`Trace`]** — drive the harness step-by-step into a
//!   full run record (initial state + the per-tick [`SensorFrame`] stream + the
//!   scene), then **[`evaluate`] → [`VnvReport`]** scores a requirement set
//!   against it (per-requirement verdict + margin + an overall AND).
//! * **Coverage** — [`parameter_coverage`] reports the exact fraction of a
//!   discrete [`ParamGrid`]'s cells a suite exercised, and
//!   [`requirement_coverage`] reports which requirements were *exercised* and
//!   *triggered* (driven to a failure) across a batch of runs.
//! * **Sweep** — [`run_suite`] runs a whole suite, aggregating **pass-rate**,
//!   **worst-case (minimum) margins** per requirement, and coverage; with
//!   [`grid_suite_auto`] (full discrete grid) and [`monte_carlo_suite`] (seeded,
//!   reproducible continuous sampling) to *generate* the suite.
//!
//! ## Margin sign convention
//!
//! Every requirement reports a signed margin with one uniform meaning: **≥ 0 ⇒
//! satisfied (that much slack); < 0 ⇒ violated (by that much)**, and `pass` is
//! exactly `margin >= 0`. This is what makes the sweep's *minimum* margin
//! meaningful — the worst case is the smallest margin, whether a nearest miss or
//! the deepest violation.
//!
//! ## Determinism
//!
//! Like `valenx-sensors`, the framework takes **no `rand` dependency**: the
//! Monte-Carlo sweep draws from the in-house seeded [`SplitMix64`], and sensor
//! noise comes from each sensor's own seeded generator. The same scenario yields
//! a byte-identical [`Trace`]; the same seed yields a byte-identical MC suite and
//! pass-rate. Bad configuration **fails loud** ([`VnvError`]): an empty suite, a
//! non-finite parameter or threshold, an empty command sequence, a zero-length
//! grid axis, or a requirement that cannot apply to a trace (a `DetectByTime`
//! over a LiDAR-less run, or any requirement over an empty trace) is a
//! recoverable error, never a silent `NaN` or a bogus pass.
//!
//! ## Honesty / scope — read this
//!
//! This is a **V&V *methodology* / harness over the analytic, model-grade
//! `valenx_sensors` models** (see that crate's own scope caveats: the sensors
//! are graphics-grade, the vehicle is kinematic, the world is analytic
//! surfaces). What this crate *validates* is the **framework logic itself** —
//! that requirements evaluate correctly (a violating trace fails with the right
//! margin, a safe one passes), that coverage is counted correctly (100% when a
//! suite spans every cell, an exact fraction with a gap), and that a seeded
//! sweep yields a reproducible pass-rate. The benchmark-pinned tests in each
//! module assert exactly these properties against hand-built ground truth.
//!
//! It is **not** a certified safety case, a real-world autonomy assurance
//! argument, or evidence that any actual autonomous system is safe. Running a
//! real assurance campaign would substitute hardware-calibrated sensor models, a
//! dynamics simulator, a real autonomy stack under test, a traceable
//! requirements baseline, and far richer coverage criteria — all of which this
//! crate is the clean, reproducible *scaffolding* for, not a replacement for.
//! This is **defensive test-&-evaluation** tooling (M9 autonomy V&V): it scores
//! whether a simulated run met stated requirements; it cues no weapon and makes
//! no targeting decision.
//!
//! ## Example
//!
//! ```
//! use nalgebra::Vector3;
//! use valenx_sensors::{Command, Scene, Sphere, VehicleState};
//! use valenx_autonomy_vnv::{
//!     evaluate, run_scenario, CommandSeq, Requirement, RequirementSet, Scenario,
//! };
//!
//! // A sphere obstacle 10 m ahead; drive straight at it for 0.5 s but stop
//! // short, so the vehicle stays clear.
//! let mut scene = Scene::new();
//! scene.push_sphere(Sphere::new(Vector3::new(10.0, 0.0, 0.0), 1.0).unwrap());
//! let scenario = Scenario::new(
//!     "approach-and-hold",
//!     VehicleState::default(),
//!     scene,
//!     CommandSeq::Constant {
//!         command: Command { accel_body: Vector3::new(1.0, 0.0, 0.0), ..Default::default() },
//!         dt: 0.1,
//!         steps: 5,
//!     },
//! );
//!
//! let trace = run_scenario(&scenario).unwrap();
//! // Require ≥ 2 m clearance from the (sphere-surface) obstacle.
//! let reqs = RequirementSet::new(vec![Requirement::MinClearance { d: 2.0 }]);
//! let report = evaluate(&reqs, &trace).unwrap();
//! assert!(report.overall_pass, "vehicle stayed clear");
//! assert!(report.worst_margin().unwrap() > 0.0);
//! ```

#![forbid(unsafe_code)]

pub mod coverage;
pub mod report;
pub mod requirement;
pub mod scenario;
pub mod sweep;
pub mod trace;

mod error;

pub use error::VnvError;

pub use coverage::{
    parameter_coverage, requirement_coverage, ParamCoverage, ParamGrid, RequirementCoverage,
};
pub use report::{evaluate, VnvReport};
pub use requirement::{Aabb, Requirement, RequirementOutcome, RequirementSet};
pub use scenario::{CommandSeq, Scenario, ScenarioSuite, SensorSet};
pub use sweep::{
    grid_suite, grid_suite_auto, monte_carlo_suite, run_suite, SampleAxis, SweepResult,
};
pub use trace::{run_scenario, Trace};

// Re-export the upstream harness types a caller needs to build scenarios, so a
// downstream user can `use valenx_autonomy_vnv::*` without also reaching into
// `valenx_sensors` for the basics.
pub use valenx_sensors::{
    Command, Gps, Imu, Lidar, LidarConfig, Scene, SensorFrame, SplitMix64, VehicleState,
};
