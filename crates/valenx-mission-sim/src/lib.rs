//! # valenx-mission-sim — a general discrete-event / agent constructive-simulation framework
//!
//! The **M2 track** of the valenx defense / modeling-&-simulation roadmap
//! (AFSIM-class mission/engagement constructive simulation). It is a **general**
//! discrete-event / agent simulation framework: the same engine serves
//! logistics, epidemiology, traffic flow, and policy wargaming as it does a
//! defensive mission analysis. valenx builds the **infrastructure + analysis**,
//! not weapons.
//!
//! ## Dual-use boundary (hard gate — read this)
//!
//! Engagement outcomes are kept **abstract and probabilistic**:
//!
//! - **Probability-of-kill (`Pk`) is an INPUT parameter**, supplied by the
//!   scenario in `[0, 1]` — never computed from a weapon/target physics model.
//!   An "engagement" is one seeded Bernoulli draw against that input, and a
//!   "hit" is just an abstract state change (an entity's `alive` flag). `Pk = 1`
//!   always hits, `Pk = 0` never.
//! - **Aggregate attrition uses Lanchester's square law** — the century-old
//!   operations-research ODE `dA/dt = −b·B`, `dB/dt = −a·A` over two aggregate
//!   forces, with its conserved quantity `a·A² − b·B²` and closed-form solution.
//!
//! There is **no** detailed lethality, **no** targeting / fire-control, and
//! **no** kill-chain logic anywhere in the crate. It is force-on-force
//! *bookkeeping* at the level of probabilities, geometry, and aggregate counts —
//! exactly the dual-use posture of the academic / think-tank constructive sims
//! the framework also serves.
//!
//! ## What it reuses (no geometry is reimplemented)
//!
//! Sensor detection reuses [`valenx_uas::detection_timeline`] — valenx's exact
//! constant-velocity range-crossing geometry (itself the civilian
//! conflict / closest-point-of-approach math). Two moving entities are reduced
//! to a single relative track and handed to that solver, so the crate does not
//! re-derive the detection quadratic.
//!
//! ## The modules
//!
//! - [`scheduler`] — a deterministic min-heap discrete-event [`Scheduler`]:
//!   events ordered by simulated time (ties broken by insertion order), an
//!   event may schedule future events, **no wall clock** anywhere.
//! - [`entity`] — [`Entity`] state (position, side, liveness, sensor /
//!   engagement ranges, the `Pk` input) and analytic [`Mover`]s
//!   (constant-velocity and waypoint-follow), all closed-form.
//! - [`sensor`] — range-based [`detect`]ion built on the reused `valenx-uas`
//!   geometry.
//! - [`engagement`] — the abstract [`resolve_pk`] draw and the
//!   [`lanchester_run`] / [`lanchester_square_step`] attrition ODE.
//! - [`scenario`] — a [`Scenario`] that wires entities + movers + sensors +
//!   engagements onto the scheduler, runs to a stop time, and returns a
//!   [`ScenarioResult`] (timeline + final state + [`OutcomeMetrics`]).
//!
//! ## Determinism
//!
//! Everything is seeded and reproducible: the scheduler uses simulated time
//! only, and every stochastic engagement draw comes from one in-crate
//! [`SplitMix64`] (no `rand` dependency). The same seed replays a bit-for-bit
//! identical timeline on every run and machine. The PRNG is **not** used for any
//! security purpose.
//!
//! ## Honest scope
//!
//! Research / educational grade. The movers are analytic kinematics (no
//! dynamics, control, or terrain); detection is geometric line-of-sight by range
//! only (the analytic caveats of `valenx-uas` / `valenx-sensors` apply — no
//! occlusion, clutter, propagation, or tracking filter); engagement is the
//! abstract `Pk` / Lanchester abstraction above, not a fidelity combat model.
//! The scenario loop is an explicit fixed-cadence time discretisation, so
//! detection/engagement *times* are resolved to within one tick (the per-pair
//! geometry inside a tick is exact). Nothing here is accredited (VV&A) or a
//! substitute for a validated mission-level analysis tool.
//!
//! ## Example
//!
//! ```
//! use nalgebra::Vector3;
//! use valenx_mission_sim::{Entity, Mover, Scenario, Side};
//!
//! // A blue sensor at the origin (500 m range) and an inbound red track.
//! let blue = Entity::new(
//!     Vector3::zeros(), Side::Blue, Mover::Static,
//!     500.0, /* sensor_range */ 0.0 /* engage_range */, 0.0 /* pk */,
//! ).unwrap();
//! let red = Entity::new(
//!     Vector3::new(1000.0, 0.0, 0.0), Side::Red,
//!     Mover::ConstantVelocity(Vector3::new(-100.0, 0.0, 0.0)),
//!     0.0, 0.0, 0.0,
//! ).unwrap();
//!
//! let scenario = Scenario::new(vec![blue, red], 20.0 /* stop */, 0.01 /* tick */, 7).unwrap();
//! let result = scenario.run().unwrap();
//!
//! // Red crosses 500 m at t = 5 s; the framework records the first detection
//! // there (to within one tick).
//! let ttfd = result.metrics.time_to_first_detection_s.unwrap();
//! assert!((ttfd - 5.0).abs() <= 0.02);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod engagement;
pub mod entity;
pub mod error;
pub mod scenario;
pub mod scheduler;
pub mod sensor;

mod rng;
mod scenario_metrics;

pub use engagement::{
    lanchester_run, lanchester_square_step, resolve_pk, square_law_invariant, EngagementOutcome,
    ForceState,
};
pub use entity::{Entity, Mover, Side};
pub use error::MissionError;
pub use rng::SplitMix64;
pub use scenario::{survivors_on, Event, Scenario, ScenarioResult, TimelineEntry};
pub use scenario_metrics::OutcomeMetrics;
pub use scheduler::{ScheduledEvent, Scheduler};
pub use sensor::{detect, range_between, Detection};
