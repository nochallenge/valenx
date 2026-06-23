//! # valenx-rotor
//!
//! Propeller / rotor performance by **blade-element-momentum theory
//! (BEMT)** — the headline blade-element aerodynamics capability for the
//! Valenx aero stack (the sibling `valenx-drone` crate covers only lumped
//! hover momentum theory).
//!
//! ## What
//!
//! Discretise a blade into radial elements, give each a chord, twist and
//! airfoil polar, then for each element couple 2-D blade-element
//! aerodynamics with annular momentum theory and solve for the local
//! inflow angle. Integrating the element loads over the span yields the
//! rotor's thrust, torque, shaft power, propeller efficiency and hover
//! figure of merit.
//!
//! ```
//! use valenx_rotor::{Rotor, Polar};
//!
//! // A 2-blade, 0.15 m-radius propeller with a tapered, twisted blade.
//! let radii  = [0.03, 0.06, 0.09, 0.12, 0.15];
//! let chords = [0.025, 0.022, 0.018, 0.014, 0.008];
//! let twist  = [
//!     25_f64.to_radians(), 18_f64.to_radians(), 13_f64.to_radians(),
//!     10_f64.to_radians(), 8_f64.to_radians(),
//! ];
//! let rotor = Rotor::from_slices(2, 0.15, 0.02, &radii, &chords, &twist).unwrap();
//!
//! // Forward flight: 5000 rpm, 5 m/s axial inflow, sea-level air.
//! let perf = rotor.solve(5000.0, 5.0, 1.225).unwrap();
//! assert!(perf.thrust_n > 0.0 && perf.power_w > 0.0);
//!
//! // Hover (V = 0) reports a physical figure of merit in (0, 1].
//! let hover = rotor.solve(6000.0, 0.0, 1.225).unwrap();
//! assert!(hover.figure_of_merit > 0.0 && hover.figure_of_merit <= 1.0);
//! ```
//!
//! ## Model (clean-room from the standard equations)
//!
//! Per radial element at radius `r` (`hub <= r <= R`), chord `c(r)`, twist
//! `beta(r)` (rad), `n` blades, `Omega = rpm * 2*pi/60`, axial freestream
//! `V`, density `rho`:
//!
//! - solidity `sigma = n*c / (2*pi*r)`; angle of attack `alpha = beta - phi`
//!   with `phi` the unknown inflow angle;
//! - section polar `Cl(alpha)`, `Cd(alpha)` — either the built-in analytic
//!   thin-airfoil polar or a caller-supplied table ([`Polar`]);
//! - rotor-frame force coefficients (propeller convention)
//!   `Cn = Cl cos phi - Cd sin phi`, `Ct = Cl sin phi + Cd cos phi`;
//! - Prandtl tip loss `F_tip = (2/pi) acos(exp(-(n/2)(R-r)/(r |sin phi|)))`,
//!   hub loss analogously, total `F = F_tip F_hub` (acos argument clamped,
//!   `sin phi` and the exp argument guarded);
//! - the inflow angle is found by balancing the elemental thrust two ways
//!   — blade element `dT_be/dr = 0.5 rho W^2 (n c) Cn` against annular
//!   momentum `dT_mom/dr = 4 pi rho r F U_a v_i` (with `U_a = Omega r tan phi`
//!   the through-disk axial speed, `v_i = U_a - V`, `W^2 = U_a^2 + U_t^2`),
//!   plus a Glauert/Buhl high-thrust correction in the windmill-brake
//!   state (turbine-convention induction `a = 1 - U_a/V` above ~0.4). The
//!   residual `g(phi) = dT_be/dr - dT_mom/dr` is solved on
//!   `(eps, pi/2 - eps)` by scanning for the first physical sign change
//!   and refining with **bisection** (guaranteed convergence, capped
//!   iterations). This thrust-balance residual is well-posed in hover
//!   (`V = 0`), unlike the bare `tan phi = V(1-a)/(Omega r (1+a'))` form
//!   which there degenerates to `phi = 0`.
//!
//! Integrating `dT/dr` and `dQ/dr` (trapezoid) gives thrust `T`, torque
//! `Q`, power `P = Q Omega`, efficiency `eta = T V / P` (for `V > 0`) and
//! the hover figure of merit `FM = (T^1.5 / sqrt(2 rho A)) / P` (`A = pi R^2`).
//!
//! ## Validation
//!
//! The crate's tests check the physical correctness gates: the hover
//! figure of merit stays in `(0, 1]` (it cannot beat the actuator-disk
//! limit) — for the 2-blade 0.15 m example prop at 6000 rpm it computes
//! **FM ~ 0.65**, in the plausible 0.6-0.8 band of a real rotor; thrust
//! increases monotonically with rpm; the forward-flight efficiency is
//! physical (`0 < eta < 1`, ~0.54 for the example at 6000 rpm / 8 m/s);
//! zero / negative / non-finite inputs return `Err` (never a NaN or
//! panic); and a degenerate element at the hub does not panic. These are
//! self-consistency and sanity checks, **not** a calibration against
//! measured propeller data — see the honesty note below.
//!
//! ## Robustness & honesty
//!
//! Every constructor and entry point validates that inputs are finite and
//! in-domain (positive `R`, `c`, `rho`, `rpm`, `n`; `hub < R`) and returns
//! a [`RotorError`] otherwise — never a silent `NaN`/`inf`. Every divide,
//! the `acos` domain and the no-convergence case are guarded; a station
//! whose inflow root cannot be bracketed or converged returns
//! [`RotorError::NoConvergence`] (iterations are capped).
//!
//! This is **research/educational grade**: 1-D strip theory with the
//! standard engineering corrections (Prandtl loss, Glauert/Buhl
//! high-induction). It omits radial/3-D and vortex effects, unsteady and
//! dynamic-stall aerodynamics, compressibility, and Reynolds-number
//! variation of the polar; the analytic polar is a thin-airfoil
//! idealisation. The trends it produces (thrust vs rpm, the FM / efficiency
//! being physical) are reliable, but the **absolute magnitudes are not
//! calibrated** — with the idealised analytic polar (low profile drag, no
//! Reynolds effects) it tends to be optimistic. It is a reasonable first
//! estimate, **not** a substitute for a vortex-lattice / CFD analysis, and
//! quantitative use needs measured airfoil data ([`Polar::Table`]) and
//! validation against test data.
//!
//! The math is implemented clean-room from public-domain BEMT textbook
//! equations; the residual-and-loss formulation follows the standard
//! approach also used by open-source tools such as NREL/BYU CCBlade, but
//! the code here is an independent reimplementation from the equations.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod airfoil;
mod bemt;
mod error;

pub use airfoil::{Polar, PolarSample, DEFAULT_CD0, DEFAULT_CL_MAX, DEFAULT_K};
pub use bemt::{BladeStation, ElementResult, Rotor, RotorPerformance};
pub use error::{ErrorCategory, RotorError};

/// Nominal sea-level air density (kg/m^3, ISA 15 C). A convenient default
/// for the `air_density` argument to [`Rotor::solve`].
pub const SEA_LEVEL_AIR_DENSITY: f64 = 1.225;
