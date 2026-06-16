//! # valenx-vibration
//!
//! Single-degree-of-freedom (SDOF) **mass-spring-damper** vibration
//! analysis: the textbook closed-form theory of one mass `m` on a linear
//! spring `k` and a viscous damper `c`, packaged as small, validated,
//! `no-unsafe` Rust.
//!
//! ## What
//!
//! Build an [`SdofSystem`] from `(m, k, c)` and read back every standard
//! modal quantity, the time-domain free response, and the steady-state
//! resonance behaviour:
//!
//! - **Modal descriptors** ([`model`]) — undamped natural frequency
//!   [`wn = sqrt(k/m)`](SdofSystem::natural_freq_rad_s), critical damping
//!   [`c_crit = 2*sqrt(k*m)`](SdofSystem::critical_damping), damping ratio
//!   [`zeta = c/c_crit`](SdofSystem::damping_ratio), the
//!   [`DampingRegime`] classification, and the damped natural frequency
//!   [`wd = wn*sqrt(1 - zeta^2)`](SdofSystem::damped_freq_rad_s).
//! - **Free response** ([`response`]) — the closed-form `x(t)` for the
//!   under-, critically- and over-damped cases via [`FreeResponse`],
//!   solved once from initial conditions `(x0, v0)` and evaluated at any
//!   time.
//! - **Resonance & decay** ([`metrics`]) — the harmonic
//!   [`magnification_factor`], the resonant peak location
//!   ([`resonant_frequency_ratio`]) and height ([`peak_magnification`]),
//!   the [`transmissibility`] that governs vibration isolation, and the
//!   [`logarithmic_decrement`] together with its exact inverse
//!   [`damping_ratio_from_decrement`] and the data-driven
//!   [`decrement_from_peaks`].
//!
//! ```
//! use valenx_vibration::{SdofSystem, FreeResponse, logarithmic_decrement};
//!
//! // m = 1 kg, k = 400 N/m, c = 4 N·s/m.
//! let sys = SdofSystem::new(1.0, 400.0, 4.0).expect("valid");
//! assert!((sys.natural_freq_rad_s() - 20.0).abs() < 1e-9); // sqrt(400)
//! assert!(sys.damping_ratio() < 1.0);                       // underdamped
//!
//! let wd = sys.damped_freq_rad_s().expect("underdamped");
//! assert!(wd < sys.natural_freq_rad_s());                   // wd < wn
//!
//! let resp = FreeResponse::new(&sys, 1.0, 0.0).expect("valid");
//! assert!((resp.displacement(0.0) - 1.0).abs() < 1e-9);     // x(0) = x0
//!
//! let delta = logarithmic_decrement(&sys).expect("underdamped");
//! assert!(delta > 0.0);
//! ```
//!
//! ## Model
//!
//! The free motion solves the linear constant-coefficient ODE
//!
//! ```text
//! m x'' + c x' + k x = 0
//! ```
//!
//! whose characteristic roots `s = wn*(-zeta +/- sqrt(zeta^2 - 1))` give
//! the three regimes (oscillatory decay, repeated real root, two real
//! roots). The harmonic-forcing results use the steady-state particular
//! solution of `m x'' + c x' + k x = F0 cos(w t)`, normalised by the
//! static deflection `F0/k` to the dimensionless magnification factor
//! `M(r, zeta) = 1/sqrt((1 - r^2)^2 + (2*zeta*r)^2)` with `r = w/wn`.
//! These are the standard equations from Rao, *Mechanical Vibrations*,
//! and Thomson, *Theory of Vibration with Applications* — implemented
//! verbatim, with each module documenting the formula it evaluates.
//!
//! ## Honest scope
//!
//! Research / educational grade. This crate implements well-established
//! **closed-form, linear, single-degree-of-freedom** vibration theory —
//! the exact textbook equations — and nothing more. It is **not** a
//! clinical/medical tool and **not** a production structural-engineering
//! certification tool. In particular it deliberately does not model:
//!
//! - **Multi-DOF or continuous systems** — no mass/stiffness/damping
//!   matrices, no mode shapes, no beams/plates/FE. One mass only.
//! - **Nonlinearity** — the spring and damper are strictly linear; no
//!   Coulomb/dry friction, hysteresis, hardening springs, or
//!   amplitude-dependent damping.
//! - **General forcing** — only free decay and *steady-state* single-
//!   frequency harmonic forcing are covered; there is no transient
//!   forced response, arbitrary input, random vibration, or shock
//!   spectrum.
//! - **Numerical time integration** — every result is a closed-form
//!   evaluation, not an ODE solver, so there is no step-size or
//!   stability concern but also no support for the cases above.
//!
//! All fallible entry points return [`VibrationError`], which carries a
//! stable [`code`](VibrationError::code) and [`category`](VibrationError::category)
//! for telemetry.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod metrics;
pub mod model;
pub mod response;

pub use error::{ErrorCategory, VibrationError};
pub use metrics::{
    damping_ratio_from_decrement, decrement_from_peaks, logarithmic_decrement,
    magnification_at_resonance, magnification_factor, peak_magnification, resonant_frequency_ratio,
    transmissibility,
};
pub use model::{DampingRegime, SdofSystem};
pub use response::FreeResponse;

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for cross-module analytic checks.
    const EPS: f64 = 1e-9;

    #[test]
    fn end_to_end_underdamped_workflow() {
        // A representative underdamped case threaded through the whole
        // public surface, checking the pieces agree with each other.
        let sys = SdofSystem::new(2.0, 200.0, 8.0).expect("valid");

        // wn = sqrt(k/m) = sqrt(100) = 10.
        assert!((sys.natural_freq_rad_s() - 10.0).abs() < EPS);
        // c_crit = 2 sqrt(k m) = 2 sqrt(400) = 40 ; zeta = 8/40 = 0.2.
        assert!((sys.damping_ratio() - 0.2).abs() < EPS);
        assert_eq!(sys.regime(), DampingRegime::Underdamped);

        // wd = 10 sqrt(1 - 0.04) = 10 sqrt(0.96).
        let wd = sys.damped_freq_rad_s().expect("ud");
        assert!((wd - 10.0 * 0.96_f64.sqrt()).abs() < EPS);

        // Free response from rest at x0 = 1 returns to x0 at t = 0.
        let resp = FreeResponse::new(&sys, 1.0, 0.0).expect("valid");
        assert!((resp.displacement(0.0) - 1.0).abs() < EPS);

        // Decrement round-trips back to zeta.
        let delta = logarithmic_decrement(&sys).expect("ud");
        let zeta_back = damping_ratio_from_decrement(delta).expect("ok");
        assert!((zeta_back - sys.damping_ratio()).abs() < 1e-12);

        // Resonant peak exists and exceeds the r = 1 value.
        let m_peak = peak_magnification(&sys).expect("peak");
        let m_res = magnification_at_resonance(&sys).expect("res");
        assert!(m_peak >= m_res);
    }

    #[test]
    fn system_serde_round_trip() {
        // The public model serializes and deserializes losslessly.
        let sys = SdofSystem::new(1.5, 60.0, 3.0).expect("valid");
        let json = serde_json::to_string(&sys).expect("serialize");
        let back: SdofSystem = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(sys, back);
        // Derived quantities survive the round-trip.
        assert!((sys.natural_freq_rad_s() - back.natural_freq_rad_s()).abs() < EPS);
    }

    #[test]
    fn free_response_serde_round_trip() {
        let sys = SdofSystem::from_modal(10.0, 0.3).expect("valid");
        let resp = FreeResponse::new(&sys, 0.5, -0.2).expect("valid");
        let json = serde_json::to_string(&resp).expect("serialize");
        let back: FreeResponse = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(resp, back);
        assert!((resp.displacement(0.13) - back.displacement(0.13)).abs() < EPS);
    }
}
