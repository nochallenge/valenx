//! # valenx-acidbase
//!
//! Native aqueous **acid-base equilibrium** models for Valenx.
//!
//! ## What
//!
//! Closed-form, textbook acid-base chemistry over dilute aqueous
//! solutions at a single temperature:
//!
//! - [`mod@ph`] — the logarithmic [`ph::ph`] / [`ph::poh`] definitions,
//!   the hydronium / hydroxide conversions, and the water autoionization
//!   relation `pH + pOH = pKw`.
//! - [`strong`] — pH and pOH of strong (fully dissociating) monoprotic
//!   acids and bases.
//! - [`weak`] — the weak monoprotic acid / base equilibrium, both the
//!   classic `[H+] ~= sqrt(Ka*C)` approximation and the exact quadratic
//!   solution, plus fraction-dissociated.
//! - [`buffer`] — the Henderson-Hasselbalch equation and the Van Slyke
//!   buffer capacity of a weak acid / conjugate-base buffer.
//!
//! ## Model
//!
//! All quantities use molar concentration (`mol / L`). The free-proton
//! activity is approximated by its concentration, so
//! `pH = -log10([H+])` exactly (the activity-coefficient correction of
//! real, non-ideal solutions is **not** modelled). Equilibrium is the
//! standard mass-action treatment:
//!
//! - Water:    `Kw = [H+][OH-]`, with `pKw = -log10(Kw)`. The default
//!   `Kw = 1.0e-14` (`pKw = 14`) holds at 25 C; pass a different `Kw`
//!   for another temperature.
//! - Strong acid `HA -> H+ + A-`: `[H+] = C` (dissociation is complete).
//! - Weak acid `HA <-> H+ + A-`: `Ka = [H+][A-] / [HA]`. The
//!   small-`x` approximation drops the `-x` in the denominator to give
//!   `[H+] = sqrt(Ka*C)`; [`weak::ph_weak_acid_exact`] keeps it and
//!   solves the resulting quadratic.
//! - Buffer: the log form of the same equilibrium,
//!   `pH = pKa + log10([A-] / [HA])`.
//!
//! Dilute, ideal-solution assumptions: a single weak acid (no diprotic
//! or polyprotic ladder), no ionic-strength / activity correction, and
//! water self-ionization folded in only through `Kw`. The strong / weak
//! solvers assume the acid dominates over the `1.0e-7 M` of water protons
//! (the usual textbook regime, valid for `C` well above `1.0e-6 M`); the
//! validated constructors reject non-positive inputs but do not enforce
//! this dilution bound, so results in the extreme-dilution limit are
//! documented as out of scope rather than guarded.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are the **textbook closed-form
//! and numerical models** taught in introductory chemistry — useful for
//! teaching, sanity checks, and back-of-the-envelope estimates. This is
//! **NOT a clinical, medical, or production engineering tool**: it does
//! not model activity coefficients, ionic strength, temperature beyond a
//! caller-supplied `Kw`, polyprotic equilibria, complexation, or
//! precipitation, and must not be used to make medical, safety, or
//! regulatory decisions.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod buffer;
pub mod error;
pub mod ph;
pub mod strong;
pub mod weak;

pub use buffer::{buffer_capacity, henderson_hasselbalch, Buffer};
pub use error::{AcidBaseError, Result};
pub use ph::{h_from_ph, oh_from_poh, ph, ph_from_poh, poh, poh_from_ph, KW_25C, PKW_25C};
pub use strong::{ph_strong_acid, ph_strong_base};
pub use weak::{fraction_dissociated, ph_weak_acid, ph_weak_acid_exact, WeakAcid};
