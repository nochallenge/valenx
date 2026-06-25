//! # valenx-uas — small-UAS design → performance → trade study (+ defensive counter-UAS geometry)
//!
//! The **M1 track** of the valenx defense / modeling-&-simulation roadmap: a
//! *fast, iterative* small-UAS design-and-simulate loop. It is the **same tool
//! a civilian drone designer uses** — multirotor and fixed-wing vehicle
//! assembly, integrated performance (hover power, endurance, range, payload),
//! and mission/trade studies (sweep design parameters, Pareto fronts) — plus a
//! **defensive** counter-UAS layer that is *pure intercept geometry and a
//! detection timeline*.
//!
//! ## Dual-use boundary (hard gate)
//!
//! This crate is **defensive / design only**. It builds the
//! civilian-drone-identical design tools and the defensive counter-UAS
//! *detect / track / intercept-**geometry*** — the time and geometry of
//! interception and the timeline on which a sensor first sees an inbound
//! track. It contains **no weapon employment, no targeting, no lethality**
//! model of any kind. The intercept module ([`intercept`]) answers only
//! "*when* and *where* could an interceptor reach this track, given its
//! speed?" and "*when* does a sensor of range `R` first detect it?" — pure
//! kinematics. What happens at the rendezvous is out of scope and not modeled.
//!
//! ## What it composes (no physics is reimplemented)
//!
//! The performance numbers come from valenx's existing, individually validated
//! in-house aerodynamics crates — this crate *integrates* them into one vehicle
//! and a trade-study / counter-UAS workflow:
//!
//! - [`valenx_drone::Multirotor`] — multirotor **hover** by actuator-disk
//!   (momentum) theory: ideal hover power `P = T^1.5 / sqrt(2·ρ·A)`, induced
//!   velocity, disk loading, thrust-to-weight, hover endurance.
//! - [`valenx_rotor::Rotor`] — blade-element-momentum-theory (BEMT) rotor
//!   solve, used as an **independent cross-check** of the momentum-theory hover
//!   power where a blade geometry is supplied
//!   ([`MultirotorUas::bemt_hover_power_w`]).
//! - [`valenx_fixedwing::Aircraft`] — fixed-wing point performance from a
//!   parabolic drag polar: the maximum lift-to-drag ratio `(L/D)max`, stall
//!   speed and drag polar that drive the electric-Breguet range.
//!
//! Battery energy is taken directly as installed watt-hours with a usable
//! fraction (depth-of-discharge / reserve). For series/parallel pack *sizing*
//! from a cell up to that watt-hour figure, see `valenx-batterypack`.
//!
//! ## The modules
//!
//! - [`vehicle`] — [`MultirotorUas`] and [`FixedWingUas`]: validated,
//!   fail-loud vehicle definitions and their integrated performance
//!   ([`MultirotorPerformance`], [`FixedWingPerformance`]).
//! - [`trade`] — a generic N-parameter [`trade::TradeStudy`]: sweep design
//!   points, evaluate each to a set of objectives, and extract the
//!   **Pareto (non-dominated) front**.
//! - [`intercept`] — defensive counter-UAS geometry: the constant-speed
//!   pursuit quadratic ([`intercept::time_to_intercept`]), the intercept
//!   point, and a sensor [`intercept::detection_timeline`].
//!
//! ## Honest scope
//!
//! Research / educational grade. The performance models inherit exactly the
//! caveats of the crates they compose (ideal hover momentum theory; parabolic
//! point performance; no motor/ESC efficiency curves beyond a lumped factor,
//! no voltage sag, no transient or gust loads, no controls). The electric
//! range uses the closed-form electric-Breguet relation (fixed-wing) and a
//! steady cruise-power balance (multirotor). The counter-UAS layer is *exact
//! constant-velocity kinematics*, not a tracking-filter or sensor-fusion model.
//! Nothing here is flight-grade or accredited; numbers are first estimates, not
//! a substitute for detailed analysis, test data, or certification.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod intercept;
pub mod trade;
pub mod vehicle;

pub use error::UasError;
pub use intercept::{
    detection_timeline, time_to_intercept, DetectionTimeline, InterceptSolution, Interceptor,
    ThreatTrack,
};
pub use trade::{DesignPoint, Objective, ParetoFront, TradeStudy};
pub use vehicle::{
    FixedWingPerformance, FixedWingUas, MultirotorPerformance, MultirotorUas, GRAVITY,
    SEA_LEVEL_AIR_DENSITY,
};
