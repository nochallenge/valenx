//! # valenx-fracture
//!
//! Linear-elastic fracture mechanics (LEFM) — the closed-form Mode-I
//! relations that connect an applied stress, a crack, and a material's
//! fracture toughness.
//!
//! ## What
//!
//! Four textbook relations, each a small validated function:
//!
//! - **Stress-intensity factor** [`stress_intensity`] —
//!   `K = Y·σ·√(πa)`. The amplitude of the crack-tip stress
//!   singularity; the master variable of LEFM.
//! - **Critical crack length** [`critical_crack_length`] —
//!   `a_c = (1/π)·(K_Ic / (Y·σ))²`. The flaw size at which `K`
//!   reaches the toughness under a given stress.
//! - **Fracture stress** [`fracture_stress`] —
//!   `σ_f = K_Ic / (Y·√(πa))`. The stress at which a crack of a given
//!   size goes unstable.
//! - **Irwin plastic-zone radius** [`plastic_zone_radius`] —
//!   `r_p = (1/2π)·(K / σ_y)²`. The crack-tip yielded region.
//!
//! A [`Material`] bundles the two properties these need: the Mode-I
//! fracture toughness `K_Ic` and the yield strength `σ_y`.
//!
//! ```
//! use valenx_fracture::{Material, critical_crack_length, fracture_stress, stress_intensity};
//!
//! // 7075-T6 aluminium: K_Ic = 24 MPa·√m, sigma_y = 470 MPa.
//! let alloy = Material::new(24.0, 470.0).unwrap();
//!
//! // Central through-crack (Y = 1) under 200 MPa of remote tension.
//! let a_c = critical_crack_length(&alloy, 1.0, 200.0).unwrap();
//!
//! // At the critical length, K has climbed exactly to the toughness.
//! let k = stress_intensity(1.0, 200.0, a_c).unwrap();
//! assert!((k - alloy.fracture_toughness).abs() < 1e-9);
//!
//! // ...which is the same statement as: the fracture stress at a_c is 200 MPa.
//! let s = fracture_stress(&alloy, 1.0, a_c).unwrap();
//! assert!((s - 200.0).abs() < 1e-6);
//! ```
//!
//! ## Model
//!
//! Pure **Mode-I, small-scale-yielding linear-elastic fracture
//! mechanics** (Griffith 1921 / Irwin 1957). The body is linear-elastic;
//! a crack of length `a` raises the local stress field to a `1/√r`
//! singularity whose strength is the stress-intensity factor `K`. Fast
//! fracture is the single criterion `K ≥ K_Ic`. The geometry / finite-size
//! correction enters through one dimensionless factor `Y` (`1.0` for a
//! central through-crack in an infinite plate, `≈1.12` for an edge crack),
//! supplied by the caller. The plastic-zone estimate is Irwin's first-order
//! plane-stress form. Units are not fixed by the types but must be mutually
//! consistent (the shipped examples use `MPa`, `m`, `MPa·√m`).
//!
//! Inputs are validated up front: a non-positive toughness, geometry
//! factor or denominator stress, or a negative / non-finite length, is
//! rejected with a [`FractureError`] rather than producing a silent
//! `NaN`/`Inf`.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are well-established **textbook
//! closed-form** models, exercised against analytic ground truth in the
//! unit tests (linear-in-`σ` and `√a` scaling of `K`, the
//! `a_c`-to-`σ_f` round trip, the `1/σ²` fall of the critical crack
//! length, the `K²` growth of the plastic zone). It is deliberately a
//! first-cut single-criterion calculator:
//!
//! - **Mode I only.** No mixed-mode (`K_II` / `K_III`) combination, no
//!   J-integral, no elastic-plastic / large-scale-yielding fracture.
//! - **The geometry factor `Y` is an input**, not derived from a handbook
//!   of `Y(a/W)` solutions — supply the right value for your geometry.
//! - **No fatigue-crack growth** (no Paris-law `da/dN`), no R-curve /
//!   stable-tearing, no residual-stress or environmental effects.
//!
//! It is **not** a clinical, medical, or production-engineering tool and
//! must not be used for structural-integrity certification or any
//! safety-of-life assessment.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod crack;
pub mod error;
pub mod material;

pub use crack::{critical_crack_length, fracture_stress, plastic_zone_radius, stress_intensity};
pub use error::{FractureError, Result};
pub use material::Material;
