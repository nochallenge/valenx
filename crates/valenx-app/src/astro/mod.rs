//! The Astro / Launch workbench's logic + layout sub-modules.
//!
//! Mirrors the CFD-side [`crate::aero`] split:
//!
//! - [`model`] — pure, non-UI form state, the `valenx-astro`
//!   [`Vehicle`](valenx_astro::Vehicle) / [`AscentConfig`](valenx_astro::AscentConfig)
//!   builders, unit conversions and result formatters (fully
//!   `#[test]`-coverable without an egui context).
//! - [`run`] — the synchronous Run action (the bounded RK4 ascent runs
//!   on click; no background thread).
//! - [`panels`] — egui layout for the two tabs (Ascent + Planners) and
//!   the four closed-form mission planners.
//!
//! The panel host itself lives in [`crate::astro_workbench`].

pub mod model;
pub mod panels;
pub mod run;
