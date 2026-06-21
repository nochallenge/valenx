//! # valenx-app-core
//!
//! The self-contained leaf modules of the Valenx desktop application
//! ([`valenx-app`]) — the helpers that have **no** reference to the
//! root `ValenxApp` struct and no `eframe::App` dispatch.
//!
//! This is stage A1 of the `valenx-app` crate split (see
//! `docs/refactor/2026-06-20-valenx-app-split.md`). Pulling these
//! pure helpers out of the ~115k-line `valenx-app` crate means editing
//! one of them recompiles this small leaf crate instead of the whole
//! app. The big `ValenxApp` struct and its `impl` files move in a
//! later slice.
//!
//! Everything here is intentionally `ValenxApp`-free:
//!
//! - [`theme`] — token-driven egui palette (`ThemeVariant`, `apply`).
//! - [`tooltips`] / [`panel_help`] — tooltip + contextual-help text.
//! - [`histograms`] — mesh-quality histogram rendering.
//! - [`menu_ui`] — the long-menu vertical-scroll wrapper.
//! - [`workbench_ui`] — the shared right-panel header chrome.
//! - [`solver_parse`] — the `solver` id + `[sweep.derived]` parsers.
//! - [`time_format`] — `valenx_fields::TimeKey` snapshot labels.

#![forbid(unsafe_code)]
#![allow(missing_docs)] // relaxed during pre-alpha; mirrors valenx-app

pub mod histograms;
pub mod menu_ui;
pub mod panel_help;
pub mod solver_parse;
pub mod theme;
pub mod time_format;
pub mod tooltips;
pub mod workbench_ui;
