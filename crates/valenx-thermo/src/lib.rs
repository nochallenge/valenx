//! `valenx-thermo` — in-house industrial fluid thermodynamics: equations of
//! state and vapor–liquid phase behavior.
//!
//! # What this crate provides (V1)
//!
//! * [`Fluid`] — a pure component described by its critical constants
//!   (`Tc`, `Pc`) and Pitzer acentric factor `ω`, plus a small library of
//!   common fluids (CO₂, N₂, CH₄, H₂O, C₃H₈) with literature constants.
//! * [`Eos`] — a two-parameter cubic equation of state, either
//!   [`EosModel::PengRobinson`] or [`EosModel::Srk`], giving pressure
//!   `P(T, V)`, the compressibility factor `Z` (via the closed-form cubic
//!   root), and the fugacity coefficient `φ`.
//! * [`saturation_pressure`] — the pure-component vapor pressure `Psat(T)`,
//!   solved by equating liquid and vapor fugacities (Newton iteration on
//!   `ln P`, seeded by the Wilson correlation [`wilson_psat`]).
//!
//! Everything is implemented from first principles in pure Rust with no
//! numerical dependencies, so it builds anywhere the workspace does.
//!
//! # Why a cubic EoS (and not PC-SAFT / `feos`)
//!
//! The mature `feos`/`feos-core` crates implement PC-SAFT and Helmholtz
//! functionals, but at the time of writing they require `nalgebra 0.34` /
//! `ndarray 0.17` (plus `num-dual` and `quantity`), whereas the Valenx
//! workspace pins `nalgebra 0.33` / `ndarray 0.15`. Adopting `feos` would force
//! a second major version of `nalgebra` into the workspace, against the
//! single-pinned-version policy. A two-parameter cubic EoS covers the V1 scope
//! (EoS + phase behavior) exactly, with zero new dependencies and full control
//! over validation, so it is the in-house choice here. PC-SAFT can be layered
//! in later behind a feature flag once the workspace's linear-algebra pins move.
//!
//! # Example
//!
//! ```
//! use valenx_thermo::{saturation_pressure, Eos, EosModel, Fluid};
//!
//! let co2 = Fluid::carbon_dioxide();
//! let eos = Eos::new(co2, EosModel::PengRobinson);
//!
//! // Vapor compressibility of CO₂ at 350 K, 5 MPa.
//! let z = eos.z_roots(350.0, 5.0e6).unwrap().vapor;
//! assert!(z > 0.7 && z < 1.0);
//!
//! // Saturation pressure of CO₂ at 273.15 K (experiment ≈ 3.49 MPa).
//! let psat = saturation_pressure(&eos, 273.15).unwrap();
//! assert!((psat - 3.49e6).abs() / 3.49e6 < 0.05);
//! ```

#![forbid(unsafe_code)]

mod eos;
mod error;
mod fluid;
mod phase;

pub use eos::{Eos, EosModel, ZRoots};
pub use error::{Result, ThermoError, GAS_CONSTANT};
pub use fluid::Fluid;
pub use phase::{saturation_pressure, wilson_psat};
