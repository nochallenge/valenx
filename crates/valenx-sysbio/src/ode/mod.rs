//! Deterministic (ODE) simulation — features 8 through 13.
//!
//! The continuous-deterministic half of the simulation layer. A
//! reaction-network [`Model`](crate::model::Model) is turned into a
//! system of ordinary differential equations and integrated.
//!
//! - [`system`] — [`OdeSystem`]: assembles `dy/dt = S·v(y)` from a
//!   model and supplies a finite-difference Jacobian (feature 8).
//! - [`integrate`] — three integrators: fixed-step RK4, adaptive
//!   Dormand-Prince RK45, and an implicit BDF for stiff systems
//!   (features 9, 10, 11).
//! - [`steady`] — a damped-Newton steady-state solver (feature 12).
//! - [`timecourse`] — a COPASI-style time-course task with discrete
//!   events and uniform output sampling (feature 13).
//! - [`linalg`] — small dense linear-algebra helpers shared by the
//!   numerical layer.

pub mod eventdriver;
pub mod integrate;
pub mod linalg;
pub mod steady;
pub mod system;
pub mod timecourse;

pub use eventdriver::{EventDrivenTimeCourse, EventTrajectory};
pub use integrate::{integrate_rk4, rk4_step, Bdf, Rk45, Trajectory};
pub use steady::{steady_state, SteadyState};
pub use system::OdeSystem;
pub use timecourse::{Event, EventOp, Integrator, TimeCourse};
