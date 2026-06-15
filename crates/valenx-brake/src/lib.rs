//! # valenx-brake
//!
//! Textbook closed-form models for the three quantities every brake
//! calculation starts from: the friction **torque** a disc/caliper
//! brake develops, the **tension ratio** a band (or drum) brake
//! multiplies across its wrap angle, and the **kinetic energy** a
//! moving mass dumps into the brakes when it stops.
//!
//! ## What
//!
//! Three independent, dependency-free calculators:
//!
//! - [`disc`] — a disc / caliper brake. The clamp force `F` is applied
//!   by `n_pads` friction faces at an effective radius `r_eff`, giving
//!   the retarding torque `T = mu * F * n_pads * r_eff`. Helpers invert
//!   the relation to recover the clamp force or the radius needed for a
//!   target torque.
//! - [`band`] — a band / drum (capstan) brake. The belt friction
//!   ("capstan" / Euler–Eytelwein) equation relates the tight-side and
//!   slack-side tensions across a wrap angle `theta`:
//!   `T1 / T2 = exp(mu * theta)`. The net braking force at the drum
//!   surface is `T1 - T2`, and the braking torque is `(T1 - T2) * r`.
//! - [`energy`] — the rigid-body kinetic energy a translating mass
//!   carries, `E = 0.5 * m * v^2`, which (ignoring losses) is the heat
//!   the brakes must absorb to bring it to rest, plus the stopping
//!   distance implied by a constant deceleration.
//!
//! ## Model
//!
//! Everything here is **Coulomb dry friction** with a single constant
//! coefficient `mu`, plus **Newtonian rigid-body energy**. Concretely:
//!
//! - **Disc torque.** Each of the `n_pads` pad faces presses on the
//!   rotor with normal (clamp) force `F`; the Coulomb friction force per
//!   face is `mu * F`, acting at the mean (effective) radius `r_eff`, so
//!   the per-face torque is `mu * F * r_eff` and the total is
//!   `T = mu * F * n_pads * r_eff`. A floating single-piston caliper
//!   squeezes the rotor from both sides, so a typical caliper has
//!   `n_pads = 2`.
//! - **Band tension ratio.** Integrating the differential friction of a
//!   flexible band wrapping a drum gives the classic capstan result
//!   `T1 = T2 * exp(mu * theta)`, with `theta` the contact (wrap) angle
//!   in **radians**. The drum-surface braking force is `T1 - T2` and the
//!   braking torque about the drum axis is `(T1 - T2) * r`.
//! - **Braking energy & distance.** A mass `m` moving at speed `v` has
//!   translational kinetic energy `E = 0.5 * m * v^2`. Under a constant
//!   deceleration `a > 0` it stops in distance `d = v^2 / (2 a)` and
//!   time `t = v / a` — the same `v^2 / (2 a)` that the work–energy
//!   theorem gives from `E = m * a * d`.
//!
//! All quantities are SI: metres, kilograms, seconds, newtons,
//! newton-metres, joules. Angles are radians.
//!
//! ```
//! use valenx_brake::{disc, band, energy};
//!
//! // A 2-pad caliper: 8 kN clamp, mu = 0.4, 0.12 m effective radius.
//! let t = disc::disc_torque(0.4, 8_000.0, 2, 0.12).unwrap();
//! assert!((t - 768.0).abs() < 1e-9); // 0.4*8000*2*0.12 = 768 N·m
//!
//! // A band brake wrapping 270° with mu = 0.3.
//! let theta = 270.0_f64.to_radians();
//! let ratio = band::tension_ratio(0.3, theta).unwrap();
//! assert!((ratio - (0.3 * theta).exp()).abs() < 1e-12);
//!
//! // Stopping a 1500 kg car from 30 m/s.
//! let e = energy::kinetic_energy(1500.0, 30.0).unwrap();
//! assert!((e - 675_000.0).abs() < 1e-6); // 0.5*1500*900 = 675 kJ
//! ```
//!
//! ## Honest scope
//!
//! This is **research/educational grade**. Every formula here is the
//! genuine textbook closed-form result and is checked against analytic
//! ground truth in the unit tests, but the model is deliberately the
//! idealised first-order one:
//!
//! - A **single constant friction coefficient** `mu`. Real pad/rotor
//!   and band/drum friction varies strongly with temperature, speed,
//!   contact pressure and wear (fade, glazing, the µ–v curve); none of
//!   that is modelled.
//! - **No thermal model.** [`energy`] reports the kinetic energy that
//!   must be absorbed but says nothing about rotor temperature rise,
//!   heat-capacity sizing, cooling, or fade onset.
//! - **No mechanics beyond the friction interface.** Hydraulic /
//!   pneumatic actuation, caliper compliance, self-servo / self-energising
//!   effects in drum shoes, brake balance, ABS/EBD, and dynamic weight
//!   transfer are all out of scope.
//! - **Rigid-body translational energy only** — rotational inertia of
//!   wheels/driveline, rolling resistance, aero drag, and grade are not
//!   included in the stopping-distance helper.
//!
//! It is a calculator for the canonical brake equations, **not** a
//! clinical/medical or production engineering sizing tool, and must not
//! be used to certify or build a real braking system.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod band;
pub mod disc;
pub mod energy;
pub mod error;

pub use error::{BrakeError, ErrorCategory};
