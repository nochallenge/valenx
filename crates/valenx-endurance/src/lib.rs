//! # valenx-endurance — respiration & endurance physiology
//!
//! Closed-form respiratory and exercise-physiology equations for the
//! Valenx workspace: the oxygen-transport chain (the O2-hemoglobin
//! dissociation curve, arterial oxygen content, oxygen delivery), the
//! blood-lactate response to exercise intensity, and a linear
//! oxygen-uptake (VO2) model. Pure algorithms — no external processes,
//! no data files, no platform dependencies.
//!
//! ## What
//!
//! - **Oxygen transport** ([`oxygen`]) — hemoglobin saturation via the
//!   Hill equation ([`oxygen::HillCurve`] / [`oxygen::saturation`]),
//!   arterial oxygen content ([`oxygen::cao2`]), and oxygen delivery
//!   ([`oxygen::do2`]).
//! - **Lactate** ([`lactate`]) — blood lactate as a function of
//!   exercise intensity, flat near a resting baseline below the lactate
//!   threshold and rising exponentially above it
//!   ([`lactate::LactateModel`] / [`lactate::lactate_at`]).
//! - **VO2** ([`vo2`]) — steady-state oxygen uptake as a linear
//!   function of mechanical power ([`vo2::Vo2Model`] /
//!   [`vo2::vo2_at`]) and its expression relative to VO2max
//!   ([`vo2::percent_vo2max`]).
//!
//! ## Model
//!
//! Each module documents its governing equation in full; the headline
//! relations are:
//!
//! - Hill dissociation: `saturation = po2^n / (p50^n + po2^n)`, with
//!   defaults `p50 = 26.6` mmHg and `n = 2.7`.
//! - Arterial content: `cao2 = 1.34 * hb * sat + 0.003 * po2` (mL O2/dL).
//! - Oxygen delivery: `do2 = cardiac_output * cao2 * 10` (mL O2/min).
//! - Lactate: a resting baseline below the threshold (default `0.75` of
//!   VO2max) that turns up exponentially above it.
//! - VO2: `vo2 = resting + slope * power`, with a default slope of
//!   `10.8` mL O2/min per watt.
//!
//! Every fallible entry point validates its inputs and returns
//! [`Result<_, EnduranceError>`](error::EnduranceError); the error type
//! exposes stable [`code`](error::EnduranceError::code) and
//! [`category`](error::EnduranceError::category) accessors for telemetry.
//!
//! ## Honest scope
//!
//! Research and educational grade only. These are simplified textbook
//! equations with representative default constants — population-average
//! curve shapes, not a physiological model of any individual. They
//! ignore CO2, pH (the Bohr effect), 2,3-BPG, temperature, carbon
//! monoxide, altitude acclimatisation, and inter-subject variation, and
//! the lactate and VO2 shapes are illustrative rather than fitted. This
//! crate is **not** a clinical, medical, diagnostic, or production tool
//! and must not be used for diagnosis, treatment, dosing, or any
//! decision affecting a real person. Quantitative outputs are for
//! exploration and teaching only.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod lactate;
pub mod oxygen;
pub mod vo2;

pub use error::{EnduranceError, ErrorCategory};
pub use lactate::{lactate_at, LactateModel};
pub use oxygen::{cao2, do2, saturation, HillCurve};
pub use vo2::{percent_vo2max, vo2_at, Vo2Model};

#[cfg(test)]
mod tests {
    //! Cross-module integration checks that exercise the whole
    //! oxygen-transport chain end to end with physiological reference
    //! numbers.

    use super::*;

    #[test]
    fn resting_oxygen_chain_reference() {
        // A healthy adult at rest: arterial po2 ~95 mmHg, hb 15 g/dL,
        // cardiac output 5 L/min. Saturation should be high (>0.95),
        // arterial content roughly 19-21 mL/dL, and delivery near the
        // textbook ~1000 mL/min figure.
        let sat = saturation(95.0).unwrap();
        assert!(sat > 0.95, "arterial saturation {sat} unexpectedly low");

        let content = cao2(15.0, sat, 95.0).unwrap();
        assert!(
            (18.0..=21.0).contains(&content),
            "arterial content {content} outside physiological band"
        );

        let delivery = do2(5.0, content).unwrap();
        assert!(
            (900.0..=1100.0).contains(&delivery),
            "oxygen delivery {delivery} outside textbook band"
        );
    }

    #[test]
    fn intensity_drives_lactate_and_vo2_together() {
        // Above the lactate threshold both lactate and the fraction of
        // VO2max should be elevated relative to an easy effort.
        let vo2max = vo2_at(400.0).unwrap();

        let easy_vo2 = vo2_at(120.0).unwrap();
        let hard_vo2 = vo2_at(360.0).unwrap();
        let easy_pct = percent_vo2max(easy_vo2, vo2max).unwrap();
        let hard_pct = percent_vo2max(hard_vo2, vo2max).unwrap();
        assert!(
            hard_pct > easy_pct,
            "hard {hard_pct}% not above easy {easy_pct}%"
        );

        // Map those VO2max fractions (as fractions, not percent) onto the
        // lactate curve; the hard effort sits above threshold and so
        // shows more lactate.
        let easy_lac = lactate_at(easy_pct / 100.0).unwrap();
        let hard_lac = lactate_at(hard_pct / 100.0).unwrap();
        assert!(
            hard_lac > easy_lac,
            "hard-effort lactate {hard_lac} not above easy {easy_lac}"
        );
    }
}
