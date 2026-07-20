//! # valenx-buoyancy
//!
//! Closed-form hydrostatics for buoyancy and the *initial* (small-angle)
//! transverse stability of floating bodies: Archimedes' principle, the
//! upright draft of a rectangular barge, and the metacentric height
//! `GM` that decides whether a floating body is stable.
//!
//! ## What
//!
//! Three topic modules, each a handful of validated pure functions:
//!
//! - [`archimedes`] — buoyant force `F = rho * g * V`, the displaced
//!   volume needed to support a weight, mean density `m / V`, and the
//!   sink / neutral / rise verdict for a fully submerged body.
//! - [`floating`] — the [`FloatingBox`] flotation model: waterplane
//!   area, displaced volume, upright draft, freeboard, and whether the
//!   box stays afloat.
//! - [`stability`] — the rectangular-waterplane second moment, the
//!   metacentric radius `BM = I_wp / V`, the metacentric height
//!   `GM = BM - BG`, and a [`StabilityVerdict`] from the sign of `GM`.
//!
//! ## Model
//!
//! The governing equations, all standard textbook hydrostatics:
//!
//! ```text
//! F_b = rho * g * V                 (buoyant force)
//! V   = m / rho                     (displaced volume to float a mass m)
//! T   = V / (L * B)                 (upright draft of a box, waterplane L*B)
//! I   = L * B^3 / 12                (rectangular waterplane second moment)
//! BM  = I_wp / V_disp               (metacentric radius)
//! GM  = BM - BG                     (metacentric height; stable iff GM > 0)
//! ```
//!
//! The rectangular-barge metacentric radius reduces to the closed form
//! `BM = B^2 / (12 * T)`, independent of length, and a fully submerged
//! body is neutrally buoyant exactly when its mean density equals the
//! fluid density. These two analytic results are the ground truth the
//! crate's tests check against.
//!
//! ## Honest scope
//!
//! Research / educational grade: standard textbook closed-form
//! hydrostatics, validated against analytic ground truth — NOT a
//! naval-architecture or production-certified engineering tool. It
//! assumes a single homogeneous incompressible fluid, prismatic /
//! rectangular waterplanes, rigid weightless-free-surface loading, and
//! small heel angles only. It does not model hull form factors,
//! free-surface or cargo-shift corrections, large-angle righting-arm
//! (`GZ`) curves, damaged stability, dynamic seakeeping, structural or
//! fatigue limits, safety factors, or any classification-society or
//! stability-code compliance. Do not use it to certify a real vessel.
//!
//! ## Example
//!
//! ```
//! use valenx_buoyancy::archimedes::{buoyant_force, STANDARD_GRAVITY, RHO_FRESH_WATER};
//! use valenx_buoyancy::floating::FloatingBox;
//! use valenx_buoyancy::stability::{barge_metacentric_height, stability_verdict, StabilityVerdict};
//!
//! // A 10 m x 4 m x 2 m pontoon massing 40 t in fresh water.
//! let pontoon = FloatingBox::new(10.0, 4.0, 2.0, 40_000.0)?;
//!
//! // It displaces 40 m^3 and sits 1 m deep, with 1 m of freeboard.
//! let draft = pontoon.draft(RHO_FRESH_WATER)?;
//! assert!((draft - 1.0).abs() < 1e-9);
//! assert!(pontoon.floats(RHO_FRESH_WATER)?);
//!
//! // The buoyant force balances its ~392 kN weight.
//! let f = buoyant_force(RHO_FRESH_WATER, draft * pontoon.waterplane_area(), STANDARD_GRAVITY)?;
//! assert!((f - 40_000.0 * STANDARD_GRAVITY).abs() < 1e-6);
//!
//! // With the centre of gravity 0.5 m above the centre of buoyancy,
//! // GM = B^2/(12 T) - BG = 16/12 - 0.5 = 0.833 m > 0: stable.
//! let gm = barge_metacentric_height(4.0, draft, 0.5)?;
//! assert!((gm - (16.0 / 12.0 - 0.5)).abs() < 1e-9);
//! assert_eq!(stability_verdict(gm, 1e-6)?, StabilityVerdict::Stable);
//! # Ok::<(), valenx_buoyancy::BuoyancyError>(())
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod archimedes;
pub mod error;
pub mod floating;
pub mod stability;

// --- Convenience re-exports of the most-used types --------------------

pub use archimedes::{SubmergedVerdict, RHO_FRESH_WATER, RHO_SEA_WATER, STANDARD_GRAVITY};
pub use error::{BuoyancyError, ErrorCategory, Result};
pub use floating::FloatingBox;
pub use stability::StabilityVerdict;
