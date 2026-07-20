//! # valenx-fin
//!
//! Closed-form **extended-surface (straight fin) heat transfer** for a
//! constant cross-section fin with an insulated (adiabatic) tip.
//!
//! ## What
//!
//! Given the convection coefficient, fin geometry, and material conductivity,
//! this crate computes the four quantities that characterise a straight fin:
//! the fin parameter `m`, the dimensionless length `mL`, the fin efficiency
//! `eta`, the fin effectiveness `epsilon_f`, and the dissipated heat rate `Q`
//! (with its long-fin ceiling `Q_max`). Every input is validated and every
//! function returns a [`Result`].
//!
//! ## Model
//!
//! The governing equation is the 1-D constant-area fin equation
//! `d^2 theta / dx^2 - m^2 theta = 0` with `theta = T - T_inf`. Its closed-form
//! solution for an insulated tip yields:
//!
//! ```text
//! m       = sqrt( h * P / (k * A) )                 fin parameter   [1/m]
//! eta     = tanh(mL) / (mL)                          efficiency      [-]
//! eps_f   = sqrt( k * P / (h * A) ) * tanh(mL)       effectiveness   [-]
//! Q       = sqrt( h * P * k * A ) * theta_b * tanh(mL)   heat rate   [W]
//! Q_max   = sqrt( h * P * k * A ) * theta_b           long-fin Q     [W]
//! ```
//!
//! with `h` the convection coefficient (W/m^2 K), `P` the wetted perimeter
//! (m), `k` the conductivity (W/m K), `A` the cross-sectional area (m^2), `L`
//! the fin length (m), and `theta_b = T_base - T_inf` the base excess (K).
//! The analytic limits are exact and are used as test ground truth: as
//! `mL -> 0`, `eta -> 1`; as `mL` grows large, `tanh(mL) -> 1` so
//! `eta -> 1/(mL)` and `Q -> Q_max`.
//!
//! ## Honest scope
//!
//! This is research/educational grade: standard textbook closed-form models,
//! validated against analytic ground truth, NOT a clinical/medical or
//! production-certified engineering tool (no component tolerances, safety
//! factors, fatigue/thermal limits, or code compliance). It assumes 1-D
//! conduction, a uniform `h`, an insulated tip, constant properties, and no
//! radiation — outside those assumptions the numbers are illustrative only.
//!
//! ## Example
//!
//! ```
//! use valenx_fin::{dimensionless_length, efficiency, fin_parameter, heat_rate};
//!
//! // An aluminium pin fin: h = 100 W/m^2K, P = 0.1 m, k = 200 W/mK,
//! // A = 1e-4 m^2, L = 0.05 m, base 80 K above the fluid.
//! let m = fin_parameter(100.0, 0.1, 200.0, 1.0e-4).unwrap();
//! let ml = dimensionless_length(100.0, 0.1, 200.0, 1.0e-4, 0.05).unwrap();
//! let eta = efficiency(ml).unwrap();
//! let q = heat_rate(100.0, 0.1, 200.0, 1.0e-4, 80.0, ml).unwrap();
//!
//! assert!((m - 500.0_f64.sqrt()).abs() < 1e-9);   // sqrt(hP/kA)
//! assert!(eta > 0.0 && eta <= 1.0);
//! assert!(q > 0.0);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod efficiency;
pub mod error;
pub mod geometry;
pub mod heat_rate;

pub use efficiency::{effectiveness, efficiency};
pub use error::{FinError, Result};
pub use geometry::{dimensionless_length, fin_parameter};
pub use heat_rate::{heat_rate, max_heat_rate};
