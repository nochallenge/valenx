//! # valenx-doseresponse
//!
//! Closed-form pharmacology dose-response models built on the Hill
//! equation, with validated inputs and analytic round-trip inverses.
//!
//! ## What
//!
//! Given a maximal effect `Emax`, a half-maximal concentration `EC50`,
//! and a Hill slope `n`, this crate computes the response to a drug
//! concentration `C`, the inverse (the concentration for a target
//! effect), the potency log-scale `pEC50`, single-site receptor
//! occupancy, and derived metrics such as the `ECx` concentrations, the
//! local slope of the curve, and the `EC90 / EC10` dynamic range. Every
//! public function validates its inputs and returns a [`Result`].
//!
//! The three modules are:
//!
//! 1. [`hill`] — the Hill equation [`HillCurve`], its fractional form,
//!    and the inverse.
//! 2. [`occupancy`](mod@occupancy) — the `n = 1` Langmuir occupancy
//!    isotherm and `pEC50` conversions.
//! 3. [`metrics`] — `ECx`, the analytic slope, and the dynamic range
//!    (methods on [`HillCurve`]).
//!
//! ## Model
//!
//! The governing equation is the Hill (Hill-Langmuir) equation
//!
//! ```text
//! E(C) = Emax * C^n / (EC50^n + C^n)
//! ```
//!
//! whose dimensionless fractional form is `f(C) = C^n / (EC50^n + C^n)`
//! in `[0, 1)`. It is monotonically increasing in `C`, equals `Emax / 2`
//! exactly at `C = EC50`, tends to `0` as `C -> 0`, and approaches
//! `Emax` as `C -> infinity`. Inverting gives
//! `C = EC50 * (E / (Emax - E))^(1/n)`. The potency log-scale is
//! `pEC50 = -log10(EC50)` (with `EC50` in molar), the steepest slope is
//! `E'(EC50) = n * Emax / (4 * EC50)`, and the 10-to-90-percent span is
//! `EC90 / EC10 = 81^(1/n)`.
//!
//! ## Honest scope
//!
//! This is research / educational grade: standard textbook closed-form
//! pharmacology models validated against analytic ground truth. It is
//! NOT a clinical, medical, or production-certified tool — it has no
//! dosing recommendations, pharmacokinetic absorption / distribution /
//! metabolism / excretion modelling, inter-patient variability, safety
//! margins, or regulatory-grade fitting. Do not use it for clinical
//! decisions.
//!
//! ## Example
//!
//! ```rust
//! use valenx_doseresponse::{HillCurve, p_ec50};
//!
//! // A drug with 100-unit max effect, EC50 = 5 (concentration units),
//! // and a Hill slope of 2.
//! let curve = HillCurve::new(100.0, 5.0, 2.0).unwrap();
//!
//! // At C = EC50 the response is exactly Emax / 2.
//! let half = curve.response(5.0).unwrap();
//! assert!((half - 50.0).abs() < 1e-9);
//!
//! // Invert: which concentration gives 80% of Emax?
//! let c80 = curve.ec_percent(80.0).unwrap();
//! assert!((curve.response(c80).unwrap() - 80.0).abs() < 1e-9);
//!
//! // Steepest slope is n*Emax/(4*EC50) = 2*100/(4*5) = 10.
//! assert!((curve.slope_at(5.0).unwrap() - 10.0).abs() < 1e-9);
//!
//! // EC50 of 1 nM (1e-9 M) is a pEC50 of 9.
//! assert!((p_ec50(1.0e-9).unwrap() - 9.0).abs() < 1e-9);
//! ```
#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod hill;
pub mod metrics;
pub mod occupancy;

pub use error::{DoseResponseError, ErrorCategory};
pub use hill::HillCurve;
pub use occupancy::{ec50_from_p, occupancy, p_ec50};

/// Crate-wide result alias: `Result<T, DoseResponseError>`.
pub type Result<T> = std::result::Result<T, DoseResponseError>;
