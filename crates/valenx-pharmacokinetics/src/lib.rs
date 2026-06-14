//! # valenx-pharmacokinetics
//!
//! Closed-form pharmacokinetic (PK) and pharmacodynamic (PD) models for
//! Valenx — the analytic textbook equations for drug disposition and
//! dose-response, with validated inputs.
//!
//! ## What
//!
//! Three small, independent models, each constructed from validated
//! parameters and queried by pure functions:
//!
//! - [`OneCompartmentIv`] — the one-compartment IV-bolus model. From a
//!   dose, a volume of distribution `V`, and a clearance `CL` it gives
//!   the elimination rate `k = CL/V`, the concentration-time curve
//!   `C(t) = dose/V·exp(-k·t)`, the half-life `ln2/k`, the total exposure
//!   `AUC = dose/CL`, the peak `Cmax = dose/V`, and the time at which the
//!   concentration falls to a chosen threshold.
//! - [`TwoCompartment`] — the biexponential disposition curve
//!   `C(t) = A·exp(-α·t) + B·exp(-β·t)`, built directly from its four
//!   macro constants, with intercept `C(0) = A + B`, terminal half-life
//!   `ln2/β`, and `AUC = A/α + B/β`.
//! - [`DoseResponse`] — the Emax / Hill sigmoidal dose-response
//!   `E(d) = Emax·d^n / (EC50^n + d^n)`.
//!
//! All quantities are unit-agnostic: pick a consistent set (e.g. mg, L,
//! L/h, h) and the outputs follow; the crate never converts units.
//!
//! ## Model
//!
//! The one-compartment IV bolus treats the body as a single well-stirred
//! volume cleared by a first-order process, so the amount obeys
//! `dA/dt = -k·A` and the concentration is a single decaying exponential.
//! The two-compartment model adds a peripheral tissue compartment, giving
//! the sum of a fast distribution exponential and a slow elimination
//! exponential; this crate parameterises that curve by the **macro
//! constants** `A, α, B, β` that one fits to plasma data, rather than
//! re-deriving the micro-rate constants. The Emax / Hill equation is the
//! standard saturable, sigmoidal receptor-occupancy form of the
//! concentration-effect relationship.
//!
//! Each model validates its physical parameters at construction time
//! (positive volume, clearance, EC50, Emax, and Hill coefficient;
//! non-negative dose) through [`PkError`], so the analytic expressions
//! stay well-defined.
//!
//! ## Honest scope
//!
//! Research/educational grade: these are the closed-form **textbook
//! PK/PD equations** for the simplest compartmental models, validated
//! against their own analytic identities (`C(0) = dose/V`,
//! `C(t½) = C(0)/2`, `AUC = dose/CL`, `E(EC50) = Emax/2`,
//! two-compartment `C(0) = A + B`). It is **NOT a clinical, medical, or
//! production dosing tool**: there is no absorption / first-pass model,
//! no non-linear (Michaelis-Menten) elimination, no inter-individual
//! variability or population PK, no Bayesian fitting, and no
//! dose-individualisation. Do not use it for patient care or any decision
//! affecting a real subject.
//!
//! ## Example
//!
//! ```rust
//! use valenx_pharmacokinetics::OneCompartmentIv;
//!
//! // dose 100 mg, V 10 L, CL 5 L/h.
//! let pk = OneCompartmentIv::new(100.0, 10.0, 5.0).unwrap();
//! assert!((pk.cmax() - 10.0).abs() < 1e-12); // dose / V
//! assert!((pk.auc() - 20.0).abs() < 1e-12); //  dose / CL
//!
//! // Concentration halves after one half-life.
//! let c0 = pk.concentration(0.0).unwrap();
//! let c_half = pk.concentration(pk.half_life()).unwrap();
//! assert!((c_half - c0 / 2.0).abs() < 1e-12);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod doseresponse;
pub mod error;
pub mod onecompartment;
pub mod twocompartment;

pub use doseresponse::DoseResponse;
pub use error::{PkError, Result};
pub use onecompartment::OneCompartmentIv;
pub use twocompartment::TwoCompartment;
