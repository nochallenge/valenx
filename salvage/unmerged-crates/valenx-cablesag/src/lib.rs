//! # valenx-cablesag
//!
//! Closed-form sag and tension of a suspended cable between two level
//! supports, in both the parabolic and catenary idealisations.
//!
//! ## What
//!
//! Given a uniform load `w`, a horizontal span `L`, and a horizontal
//! tension component `H`, this crate computes the mid-span sag, the
//! deflected profile, the cable arc length, and the maximum tension —
//! and quantifies the gap between the two classic models:
//!
//! `parabola` — load uniform per unit *horizontal* length (a cable
//! carrying a heavy uniform deck, the suspension-bridge idealisation).
//! The exact shape is a parabola.
//!
//! `catenary` — load uniform per unit *arc* length (a bare cable or
//! chain hanging under its own weight). The exact shape is a `cosh`.
//!
//! `compare` — both sags side by side plus their ratio and the
//! leading-order series gap, so a caller can judge when the simpler
//! parabola suffices.
//!
//! ## Model
//!
//! Origin at the low point, supports at `x = +/- L/2`, `a = H / w`:
//!
//! Parabola: `y(x) = w x^2 / (2 H)`,
//! `d = w L^2 / (8 H)`,
//! `T = sqrt(H^2 + (w L / 2)^2)`.
//!
//! Catenary: `y(x) = a (cosh(x/a) - 1)`,
//! `d = a (cosh(L/2a) - 1)`,
//! `s = 2 a sinh(L/2a)`,
//! `T = H cosh(L/2a) = w (a + d)`.
//!
//! In every case the support tension is the maximum and satisfies
//! `T >= H`, with equality only in the `w -> 0` limit. The two models
//! coincide in the shallow-cable limit `L / 2a -> 0`, where the catenary
//! sag expands to `w L^2 / (8 H) (1 + (L/2a)^2 / 12 + ...)`.
//!
//! References: H. M. Irvine, *Cable Structures*, MIT Press, 1981, ch. 1.
//!
//! ## Honest scope
//!
//! Research/educational grade: standard textbook closed-form models,
//! validated against analytic ground truth (hand-worked `cosh`/`sinh`
//! values and known limiting cases) — NOT a clinical/medical or
//! production-certified engineering tool. It has no component
//! tolerances, safety factors, fatigue or thermal limits, elastic
//! stretch, wind/ice loading, support-level offsets, or building-code
//! compliance, and must not be used for structural certification.
//!
//! ## Example
//!
//! ```
//! use valenx_cablesag::{parabola, catenary, compare};
//!
//! // A light deck cable: w = 2 N/m, L = 100 m, H = 2500 N.
//! let d = parabola::sag(2.0, 100.0, 2500.0).unwrap();
//! assert!((d - 1.0).abs() < 1e-12); // exactly 1 m of sag
//!
//! let t = parabola::max_tension(2.0, 100.0, 2500.0).unwrap();
//! assert!(t >= 2500.0); // support tension never below H
//!
//! // The same cable hung under self-weight sags a touch more.
//! let cat = catenary::sag(2.0, 100.0, 2500.0).unwrap();
//! assert!(cat >= d);
//!
//! // How much more, as a ratio:
//! let cmp = compare::sag_comparison(2.0, 100.0, 2500.0).unwrap();
//! assert!(cmp.ratio >= 1.0);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;

pub mod catenary;
pub mod compare;
pub mod parabola;

pub use error::CableError;
