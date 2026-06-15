//! # valenx-hemodynamics
//!
//! Closed-form cardiovascular hemodynamics — the textbook relations that
//! connect vessel geometry, blood viscosity, flow, pressure and the
//! lumped behaviour of the whole circulation.
//!
//! ## What
//!
//! A small, dependency-light library of the standard first-order
//! hemodynamic formulas, each a pure function that validates its inputs
//! and returns a [`Result`]:
//!
//! - **Single-vessel flow** ([`flow`]) — the Hagen-Poiseuille
//!   volumetric flow [`flow::poiseuille_flow`], the matching hydraulic
//!   [`flow::vascular_resistance`] and its Ohm-analogue inverse
//!   [`flow::flow_from_resistance`], and the Poiseuille
//!   [`flow::wall_shear_stress`] on the vessel wall.
//! - **Whole-circulation lumped relations** ([`circulation`]) — cardiac
//!   output [`circulation::cardiac_output`], mean arterial pressure
//!   [`circulation::mean_arterial_pressure`] and the inverse
//!   [`circulation::systemic_vascular_resistance`].
//! - **Windkessel decay** ([`windkessel`]) — the 2-element (RC)
//!   diastolic pressure relaxation [`windkessel::windkessel_pressure`]
//!   and its [`windkessel::time_constant`].
//!
//! ## Model
//!
//! For a rigid cylindrical vessel of radius `r` and length `L`, a fluid
//! of dynamic viscosity `mu`, and a pressure drop `dP`:
//!
//! ```text
//! Q   = pi * r^4 * dP / (8 * mu * L)      (Hagen-Poiseuille flow)
//! R   = 8 * mu * L / (pi * r^4)           (vascular resistance, Q = dP / R)
//! tau_wall = 4 * mu * Q / (pi * r^3)      (wall shear stress)
//! ```
//!
//! For the lumped systemic circuit driven by the heart, with heart rate
//! `HR`, stroke volume `SV`, cardiac output `CO` and systemic vascular
//! resistance `SVR`:
//!
//! ```text
//! CO  = HR * SV                           (cardiac output)
//! MAP = CO * SVR                          (mean arterial pressure)
//! ```
//!
//! And the 2-element Windkessel diastolic decay, with the time constant
//! `tau = R * C` (resistance `R` in parallel with arterial compliance
//! `C`):
//!
//! ```text
//! P(t) = P0 * exp(-t / (R * C))           (pressure relaxation)
//! ```
//!
//! All formulas are unit-agnostic: feed any self-consistent unit system
//! (SI throughout, or clinical units such as mmHg / mL / min) and the
//! result carries the matching units. The doc comments on each function
//! spell out the SI choice.
//!
//! ```
//! use valenx_hemodynamics::{
//!     cardiac_output, mean_arterial_pressure, poiseuille_flow, windkessel_pressure,
//! };
//!
//! // A pressure drop drives Poiseuille flow through an arteriole.
//! let q = poiseuille_flow(2.0e-3, 1000.0, 3.5e-3, 0.05).expect("valid");
//! assert!(q > 0.0);
//!
//! // The heart's output sets the mean arterial pressure across the circuit.
//! let co = cardiac_output(1.2, 7.0e-5).expect("valid"); // 72 bpm, 70 mL
//! let map = mean_arterial_pressure(co, 1.2e8).expect("valid");
//! assert!(map > 0.0);
//!
//! // During diastole the arterial pressure relaxes toward zero.
//! let p = windkessel_pressure(80.0, 1.2, 1.2e8, 1.0e-8).expect("valid");
//! assert!(p < 80.0);
//! ```
//!
//! ## Honest scope
//!
//! These are **research/educational-grade** closed-form models —
//! well-established textbook formulas (Hagen-Poiseuille flow, the
//! electrical / Ohm analogue of the circulation, and the lumped
//! 2-element Windkessel), reproduced exactly. They are **NOT a
//! clinical, medical, or production engineering tool** and must not be
//! used to make any diagnostic or treatment decision.
//!
//! The idealisations are deliberate and load-bearing:
//!
//! - **Steady, fully-developed, laminar flow** of an **incompressible
//!   Newtonian fluid** in a **straight, rigid, circular** tube. Real
//!   blood is a non-Newtonian shear-thinning suspension; real flow is
//!   pulsatile and can be turbulent or entrance-region; real vessels
//!   taper, curve, branch and distend.
//! - The circulation relations lump the entire systemic tree into a
//!   single resistance and a single compliance and take the venous
//!   reference pressure as zero, so `MAP = CO * SVR` omits the central
//!   venous pressure offset and all regional / pulmonary detail.
//! - The 2-element Windkessel models **diastole only** (zero inflow) and
//!   has no characteristic (aortic) impedance term, unlike the 3- and
//!   4-element variants; it does not reproduce the systolic upstroke or
//!   wave-reflection features of a real pressure waveform.
//!
//! Within those stated assumptions every number is the genuine
//! closed-form answer, and each function checks its domain so an invalid
//! input yields a typed [`HemodynamicsError`] rather than a silent
//! `NaN` or infinity.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod circulation;
pub mod error;
pub mod flow;
pub mod windkessel;

pub use circulation::{cardiac_output, mean_arterial_pressure, systemic_vascular_resistance};
pub use error::{ErrorCategory, HemodynamicsError};
pub use flow::{flow_from_resistance, poiseuille_flow, vascular_resistance, wall_shear_stress};
pub use windkessel::{time_constant, windkessel_pressure};
