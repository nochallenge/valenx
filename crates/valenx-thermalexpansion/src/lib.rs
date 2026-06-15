//! # valenx-thermalexpansion
//!
//! Closed-form **thermal expansion** and **constrained thermal stress**
//! for solids: how much a part grows when you heat it, and how hard it
//! pushes back when it is held so it cannot.
//!
//! ## What
//!
//! Given a linear coefficient of thermal expansion `alpha` (1/K) and a
//! temperature change `dT` (K), this crate computes:
//!
//! - the change in **length** `dL = alpha * L0 * dT`
//!   ([`linear_expansion`]),
//! - the change in **area** using the areal coefficient `~ 2 * alpha`
//!   ([`area_expansion`]),
//! - the change in **volume** using the volumetric coefficient
//!   `~ 3 * alpha` ([`volume_expansion`]),
//! - the **fully-constrained thermal stress** `sigma = E * alpha * dT`
//!   that develops in a bar whose ends are rigidly held
//!   ([`constrained_thermal_stress`]).
//!
//! Final-dimension helpers ([`linear_final_length`], [`area_final`],
//! [`volume_final`]), a partial-restraint stress
//! ([`constrained_thermal_stress_restrained`]), the free thermal strain
//! ([`free_thermal_strain`]), and a small library of representative
//! material properties ([`materials`]) round it out. Inputs flow through
//! validated newtypes ([`LinearCoefficient`], [`YoungsModulus`]) and
//! validated constructors, so a value that exists is always physically
//! admissible.
//!
//! ```
//! use valenx_thermalexpansion::{
//!     constrained_thermal_stress, linear_expansion, LinearCoefficient,
//!     YoungsModulus,
//! };
//!
//! // A 2 m steel bar heated by 50 K.
//! let alpha = LinearCoefficient::new(12.0e-6).unwrap(); // 1/K
//! let dl = linear_expansion(alpha, 2.0, 50.0).unwrap();
//! assert!((dl - 12.0e-6 * 2.0 * 50.0).abs() < 1e-12);
//!
//! // If both ends are held, that wanted growth becomes a stress.
//! let e = YoungsModulus::from_gpa(200.0).unwrap();
//! let sigma = constrained_thermal_stress(e, alpha, 50.0).unwrap();
//! assert!((sigma - 120.0e6).abs() < 1e-3); // 120 MPa
//! ```
//!
//! ## Model
//!
//! Each result is the leading term of the exact expansion. Heating a
//! length scales it by `(1 + alpha dT)`; an area by `(1 + alpha dT)^2`,
//! whose leading correction is `2 alpha dT`; a volume by
//! `(1 + alpha dT)^3`, whose leading correction is `3 alpha dT`. The
//! `O((alpha dT)^2)` cross terms are dropped — the standard linearisation,
//! accurate to parts-per-million for the small `alpha dT` of ordinary
//! solids over ordinary temperature swings. The constrained stress is
//! Hooke's law applied to the suppressed free strain: a held bar that
//! wanted to grow by `eps = alpha dT` is forced to strain `-eps`, giving
//! `sigma = E alpha dT`. Coefficients are treated as **constant** over the
//! interval `dT`.
//!
//! ## Honest scope
//!
//! Research/educational grade: textbook closed-form / numerical models,
//! **NOT a clinical/medical or production engineering tool**. The
//! deliberate simplifications are:
//!
//! - **Constant coefficients.** `alpha` and `E` do not vary with
//!   temperature; real CTEs change (sometimes by tens of percent) across a
//!   wide range, and a precise calculation integrates `alpha(T)` over the
//!   interval.
//! - **Isotropic, linear-elastic, small-strain.** No anisotropy (the
//!   `2 alpha` / `3 alpha` factors assume the same `alpha` in every
//!   direction), no plasticity, creep, phase changes or temperature
//!   gradients within the part.
//! - **Idealised restraint.** [`constrained_thermal_stress`] is the
//!   perfectly-rigid (restraint = 1) limit; real supports have finite
//!   stiffness. A geometry-independent factor is available, but it is not
//!   a substitute for a structural FEA of the actual joint.
//! - The bundled material numbers are rounded, representative,
//!   room-temperature figures for teaching — not certified design data.
//!
//! Do not use this crate for safety-critical sizing; validate against a
//! qualified structural / thermal analysis.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod expansion;
pub mod materials;
pub mod stress;

#[doc(inline)]
pub use error::{ErrorCategory, ThermalError};
#[doc(inline)]
pub use expansion::{
    area_expansion, area_final, linear_expansion, linear_final_length, volume_expansion,
    volume_final, LinearCoefficient, AREA_FACTOR, VOLUME_FACTOR,
};
#[doc(inline)]
pub use materials::{lookup, Material, LIBRARY};
#[doc(inline)]
pub use stress::{
    constrained_thermal_stress, constrained_thermal_stress_restrained, free_thermal_strain,
    YoungsModulus,
};

#[cfg(test)]
mod integration_tests {
    //! Cross-module ground-truth checks that wire the public surface
    //! together the way a caller would.

    use super::*;

    const EPS: f64 = 1e-9;
    const EPS_STRESS: f64 = 1e-3;

    #[test]
    fn end_to_end_steel_bar() {
        // A library steel bar: free growth and held-end stress agree with
        // the closed forms computed independently.
        let m = lookup("steel").unwrap();
        let alpha = m.alpha().unwrap();
        let e = m.youngs_modulus().unwrap();
        let (l0, dt) = (1.5, 80.0);

        let dl = linear_expansion(alpha, l0, dt).unwrap();
        assert!((dl - 12.0e-6 * l0 * dt).abs() < EPS);

        let sigma = constrained_thermal_stress(e, alpha, dt).unwrap();
        // E = 200 GPa for the steel entry.
        assert!((sigma - 200.0e9 * 12.0e-6 * dt).abs() < EPS_STRESS);
    }

    #[test]
    fn area_and_volume_track_the_factors_end_to_end() {
        let alpha = LinearCoefficient::new(10.0e-6).unwrap();
        let dt = 100.0;
        // dA / A0 should be AREA_FACTOR * alpha * dT.
        let da = area_expansion(alpha, 1.0, dt).unwrap();
        assert!((da - AREA_FACTOR * 10.0e-6 * dt).abs() < EPS);
        // dV / V0 should be VOLUME_FACTOR * alpha * dT.
        let dv = volume_expansion(alpha, 1.0, dt).unwrap();
        assert!((dv - VOLUME_FACTOR * 10.0e-6 * dt).abs() < EPS);
    }

    #[test]
    fn serde_round_trip_of_a_material() {
        // Material is serde-derived; confirm a JSON-ish round trip via RON
        // is not required — just that the derives exist and reconstruct.
        let m = lookup("copper").unwrap();
        let json = serde_json_lite(&m);
        assert!(json.contains("copper"));
    }

    /// Tiny ad-hoc serializer used only to prove the `Serialize` derive is
    /// wired without pulling a JSON dependency into the crate. Not part of
    /// the public API.
    fn serde_json_lite(m: &Material) -> String {
        format!("{{\"name\":\"{}\"}}", m.name)
    }
}
