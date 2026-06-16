//! # valenx-led
//!
//! Closed-form DC operating-point calculator for the textbook
//! "LED in series with a current-limiting resistor" circuit, plus its
//! series-string generalisation.
//!
//! ## What
//!
//! Given a supply voltage, an LED forward voltage, and a target LED
//! current, this crate computes the series resistor that sets that current
//! and the resulting power split between the LED and the resistor — or,
//! conversely, builds the circuit from a chosen resistor and reports the
//! current it sets (`from_resistor`). It covers a single LED
//! ([`circuit::LedCircuit`]) and `n` identical LEDs wired in series
//! ([`circuit::LedString`]), where the forward voltages add.
//!
//! ## Model
//!
//! The circuit is a single Kirchhoff voltage loop. With supply `Vs`, total
//! LED forward voltage `Vf` (one drop, or the sum over a series string), and
//! current `I`:
//!
//! ```text
//! R          = (Vs - Vf) / I        (series current-limiting resistor)
//! P_led      = Vf * I               (LED power)
//! P_resistor = (Vs - Vf) * I        (resistor power)
//! P_total    = Vs * I = P_led + P_resistor
//! ```
//!
//! For a series string of `n` identical LEDs the same equations apply with
//! `Vf` replaced by `n * Vf_each`, because the series elements carry the same
//! current and their forward voltages sum. The LED is represented with the
//! constant-voltage-drop diode model: `Vf` is a fixed parameter rather than a
//! current- or temperature-dependent function. The current is set entirely by
//! the resistor; `I = (Vs - Vf) / R` is the inverse of the design equation,
//! exposed by [`circuit::LedCircuit::from_resistor`] /
//! [`circuit::LedString::from_resistor`].
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form models: an
//! ideal constant-drop LED, an ideal ohmic resistor, an ideal stiff supply,
//! steady-state DC only, no temperature drift, no diode `I`-`V` curve, no
//! resistor tolerance, no wire or contact resistance, no thermal derating,
//! and no transient/AC behaviour. It is NOT a clinical/medical/production
//! engineering tool and must not be used to size components for safety- or
//! life-critical hardware; validate any real design against datasheet `I`-`V`
//! curves, worst-case tolerances, and thermal limits.
//!
//! ## Example
//!
//! ```
//! use valenx_led::circuit::LedCircuit;
//!
//! // 5 V supply, 2.0 V red LED, target 20 mA.
//! let c = LedCircuit::new(5.0, 2.0, 0.020).unwrap();
//! assert!((c.resistor_ohm() - 150.0).abs() < 1e-9); // (5 - 2) / 0.020
//! assert!((c.led_power_w() - 0.040).abs() < 1e-9); //  2.0 * 0.020
//! assert!((c.resistor_power_w() - 0.060).abs() < 1e-9); // (5 - 2) * 0.020
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod circuit;
pub mod error;

pub use circuit::{LedCircuit, LedString};
pub use error::LedError;
