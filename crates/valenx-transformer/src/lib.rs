//! # valenx-transformer — ideal two-winding transformer relations
//!
//! Closed-form relations for an ideal (and simple efficiency-derated)
//! two-winding electrical transformer, built around the turns ratio
//! `a = Np / Ns`.
//!
//! ## What
//!
//! - [`ratio`] — the turns ratio [`TurnsRatio`] and the ideal voltage /
//!   current relations `a = Vp/Vs = Is/Ip`, plus step-up / step-down /
//!   isolation classification.
//! - [`power`] — apparent power `P = V * I`, ideal power conservation
//!   `Pin = Pout`, and an [`Efficiency`] newtype giving
//!   `eta = Pout/Pin` in `(0, 1]` with output / input / loss helpers.
//! - [`impedance`] — reflected (referred) load impedance
//!   `Zp = a^2 * Zs`, its inverse, and the impedance-matching ratio
//!   `a = sqrt(Zsource / Zload)`.
//! - [`error`] — the [`TransformerError`] taxonomy with stable
//!   [`code`](TransformerError::code) / [`category`](TransformerError::category)
//!   accessors.
//!
//! ## Model
//!
//! For `Np` primary turns and `Ns` secondary turns the ideal relations
//! are
//!
//! ```text
//! a   = Np / Ns
//! Vp / Vs = a          (voltage scales with the turns ratio)
//! Ip / Is = 1 / a      (current scales inversely)
//! Pin = Vp Ip = Vs Is = Pout      (lossless: apparent power conserved)
//! Zp  = a^2 Zs         (a load reflects by the square of the ratio)
//! ```
//!
//! A real transformer is modelled only to the extent of a single
//! user-supplied efficiency `eta`, giving `Pout = eta Pin` and
//! `Ploss = (1 - eta) Pin`. Every public function validates its inputs
//! and returns [`Result<_, TransformerError>`](Result); float results
//! are computed in closed form, never iterated.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are the textbook closed-form
//! ideal-transformer relations (with an efficiency knob), suitable for
//! teaching, sanity-checking, and back-of-the-envelope sizing. They are
//! NOT a clinical, medical, or production electrical-engineering tool.
//! In particular this crate does not model: magnetising and leakage
//! reactance, winding resistance, core hysteresis / eddy-current loss
//! curves, saturation, frequency dependence, phase / complex two-port
//! behaviour, thermal rise, or insulation and safety margins. Do not
//! size real hardware from it.
//!
//! ## Example
//!
//! ```
//! use valenx_transformer::{TurnsRatio, Efficiency};
//! use valenx_transformer::impedance::reflect_to_primary;
//!
//! // 240-turn primary, 24-turn secondary => a = 10 (step-down).
//! let xfmr = TurnsRatio::from_turns(240.0, 24.0).unwrap();
//! assert!((xfmr.ratio() - 10.0).abs() < 1e-12);
//! assert!(xfmr.is_step_down());
//!
//! // Step 230 V down to the secondary.
//! let vs = xfmr.secondary_voltage(230.0).unwrap();
//! assert!((vs - 23.0).abs() < 1e-12);
//!
//! // An 8-ohm load looks like 800 ohm from the primary.
//! let zp = reflect_to_primary(&xfmr, 8.0).unwrap();
//! assert!((zp - 800.0).abs() < 1e-9);
//!
//! // A 97%-efficient unit drawing 500 W in delivers 485 W out.
//! let eta = Efficiency::new(0.97).unwrap();
//! assert!((eta.output_power(500.0).unwrap() - 485.0).abs() < 1e-9);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod impedance;
pub mod power;
pub mod ratio;

pub use error::{ErrorCategory, TransformerError};
pub use power::{apparent_power, ideal_output_power, Efficiency};
pub use ratio::TurnsRatio;

/// Crate-wide result alias: `Result<T, TransformerError>`.
pub type Result<T> = std::result::Result<T, TransformerError>;
