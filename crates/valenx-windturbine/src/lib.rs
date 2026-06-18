//! # valenx-windturbine
//!
//! Closed-form **wind-turbine aerodynamics**: how much power is in the
//! wind, how much an ideal rotor can extract from it, how fast the
//! blade tips are moving relative to the wind, and the idealised
//! power-versus-wind-speed curve of a controlled machine.
//!
//! ## What
//!
//! Four small, independently-usable pieces:
//!
//! - [`power`] — actuator-disc power extraction. The available wind
//!   power [`power::available_power`] (`1/2 rho A v^3`), the
//!   [`power::BETZ_LIMIT`] (`16/27 ~ 0.593`), and the captured shaft
//!   power [`power::extracted_power`] (`1/2 rho A v^3 Cp`), with the
//!   power coefficient `Cp` validated against the Betz limit, plus the
//!   sizing inverse [`power::rotor_radius_for_power`]
//!   (`R = sqrt(2 P / (rho pi v^3 Cp))`) that sizes the rotor for a
//!   target shaft power, and the complementary operating-speed inverse
//!   [`power::wind_speed_for_power`] (`v = cbrt(2 P / (rho A Cp))`) — the
//!   rated wind speed a fixed disc needs for that power.
//! - [`tsr`] — the tip-speed ratio [`tsr::tip_speed_ratio`]
//!   (`lambda = omega R / v`), plus rad/s ↔ rev/min converters.
//! - [`curve`] — the idealised [`curve::PowerCurve`] with cut-in,
//!   rated, and cut-out break-points and its three operating
//!   [`curve::Region`]s.
//! - [`error`] — the [`error::WindTurbineError`] taxonomy returned by
//!   every fallible call.
//!
//! ## Model
//!
//! Air of density `rho` (kg/m³) flowing at free-stream speed `v` (m/s)
//! through a rotor disc of swept area `A = pi R^2` (m²) carries a
//! kinetic-energy flux
//!
//! ```text
//! P_avail = 1/2 * rho * A * v^3.
//! ```
//!
//! One-dimensional momentum theory (Betz, 1919) caps the extractable
//! fraction at `Cp_max = 16/27 ~ 0.5926`: slowing the air more than
//! that reduces the mass flow faster than it raises the per-unit-mass
//! energy drop, so the product peaks. The shaft power is
//!
//! ```text
//! P = 1/2 * rho * A * v^3 * Cp,   0 <= Cp <= 16/27,
//! ```
//!
//! and the tip-speed ratio `lambda = omega R / v` selects the operating
//! point on a turbine's `Cp(lambda)` curve. A real controller then maps
//! wind speed to electrical power as a piecewise curve: zero below
//! cut-in, a cube-law ramp up to rated, a flat plateau at rated power up
//! to cut-out, and zero above cut-out (see [`curve`]).
//!
//! ## Honest scope
//!
//! Research/educational grade. Every formula here is the **textbook
//! closed form** from one-dimensional actuator-disc momentum theory —
//! the same equations in any wind-energy primer — and the unit tests
//! check them against analytic, hand-computed values (the cube-in-speed
//! and linear-in-area scaling, `16/27 ~ 0.593`, extracted `<=` Betz
//! `<=` available, the three power-curve regions, and `lambda`). It is
//! **not** a blade-element-momentum (BEM) aerodynamics solver, does not
//! model airfoil lift/drag polars, wake/induction factors, tip and hub
//! losses, dynamic stall, tower shadow, turbulence, or structural
//! loads, and is **NOT a clinical/medical/production engineering tool**.
//! Do not use it to certify, site, or load-rate a real machine.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod curve;
pub mod error;
pub mod power;
pub mod tsr;

pub use curve::{PowerCurve, Region};
pub use error::{ErrorCategory, WindTurbineError};
pub use power::{
    available_power, betz_power, extracted_power, rotor_radius_for_power, swept_area, validate_cp,
    wind_speed_for_power, AIR_DENSITY_SEA_LEVEL, BETZ_LIMIT,
};
pub use tsr::{rad_per_s_to_rpm, rpm_to_rad_per_s, tip_speed, tip_speed_ratio};
