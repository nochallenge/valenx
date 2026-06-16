//! # valenx-bjt
//!
//! DC biasing of bipolar junction transistors (BJTs).
//!
//! ## What
//!
//! Closed-form analysis of the operating point (the quiescent or
//! "Q-point") of a single BJT in the two textbook bias topologies:
//!
//! - **Fixed-base bias** — one resistor `Rb` from the supply into the
//!   base, plus a collector resistor `Rc` (and an optional emitter
//!   resistor `Re`).
//! - **Voltage-divider bias** — the classic four-resistor network
//!   (`R1`, `R2`, `Rc`, `Re`) whose base side is reduced to a Thevenin
//!   source before solving the base loop.
//!
//! Both produce a [`bias::OperatingPoint`] holding the base, collector
//! and emitter currents, the emitter voltage drop, the
//! collector-emitter voltage, and the [`model::Region`] the device
//! sits in (active vs. saturation).
//!
//! Each topology also reports its DC **bias stability factor**
//! `S(ICO) = ∂Ic/∂Ico` ([`bias::FixedBias::stability_factor`],
//! [`bias::DividerBias::stability_factor`]) — the standard figure of
//! merit for how strongly leakage / temperature drift would move the
//! quiescent point, bounded by `1 <= S <= beta + 1`.
//!
//! ## Model
//!
//! The crate uses the standard large-signal DC model with two
//! simplifying assumptions, valid for hand analysis:
//!
//! - A **constant** base-emitter turn-on drop `VBE` (e.g. 0.7 V for
//!   silicon) whenever the device conducts.
//! - A **constant** forward current gain `beta` (the small-signal and
//!   DC betas are taken equal), so the terminal currents obey
//!
//! > `Ic = beta * Ib`,  `Ie = (beta + 1) * Ib`,  `Ie = Ic + Ib`.
//!
//! For voltage-divider bias the base divider `R1` / `R2` is replaced by
//! its Thevenin equivalent `Vth = Vcc * R2 / (R1 + R2)`,
//! `Rth = R1 || R2`, and Kirchhoff's voltage law around the base loop
//!
//! > `Vth = Ib*Rth + VBE + Ie*Re`
//!
//! is solved (with `Ie = (beta + 1) * Ib`) for the base current
//!
//! > `Ib = (Vth - VBE) / (Rth + (beta + 1) * Re)`.
//!
//! The emitter voltage is `VE = Ie * Re` and the collector-emitter
//! voltage is `Vce = Vcc - Ic*Rc - Ie*Re`. The transistor is reported
//! as [`model::Region::Active`] while `Vce > Vce_sat`; once the
//! active-region collector current would drive `Vce` down to or below
//! the saturation floor the device is reported
//! [`model::Region::Saturation`] and the collector current is clamped
//! to the saturation value `Ic_sat = (Vcc - Vce_sat) / (Rc + Re)`.
//!
//! Differentiating the same base-loop solution with respect to the
//! reverse-saturation current `Ico` gives the bias **stability factor**
//!
//! > `S(ICO) = (beta + 1)(RB + Re) / (RB + (beta + 1) Re)`
//!
//! (with `RB = Rb` for a fixed bias and `RB = Rth` for a divider), equal
//! to `(1 + beta) / (1 + beta Re / (Re + RB))`. It depends only on the
//! network and `beta`, collapses to `beta + 1` for a bare fixed bias
//! (`Re = 0`) and to `1` for an ideal emitter bias (`RB = 0`), and obeys
//! `1 <= S <= beta + 1`.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are the textbook constant-`VBE`,
//! constant-`beta` hand-analysis equations (Sedra-Smith, Razavi,
//! Boylestad) — they are **not** a clinical/medical/production
//! engineering tool. The model deliberately ignores the Early effect
//! (finite output resistance / base-width modulation), the exponential
//! Ebers-Moll / Gummel-Poon `Ic = Is * exp(Vbe/Vt)` law, temperature
//! and `beta` spread, leakage, high-level injection, and breakdown.
//! Numbers it returns are first-order design estimates for learning and
//! sizing, not a substitute for a SPICE simulation or measured device
//! data when qualifying real silicon.
//!
//! The reported stability factor `S(ICO)` is the first-order figure of
//! merit for exactly that ignored leakage / temperature / `beta` drift:
//! it quantifies how strongly such drift *would* move the Q-point in a
//! given network, but the crate still does not model the drift itself.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bias;
pub mod error;
pub mod model;

pub use bias::{DividerBias, FixedBias, OperatingPoint};
pub use error::{BjtError, ErrorCategory};
pub use model::{Region, Transistor};
