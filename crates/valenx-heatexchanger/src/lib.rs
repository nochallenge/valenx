//! # valenx-heatexchanger
//!
//! Two-stream heat-exchanger thermal analysis by the two standard
//! textbook methods: the **log-mean-temperature-difference (LMTD)**
//! method and the **effectiveness-NTU** method.
//!
//! ## What
//!
//! - [`FlowArrangement`] — counterflow vs. parallel-flow (co-current)
//!   stream layout, the one degree of freedom shared by both methods.
//! - LMTD method ([`lmtd`] module): [`TerminalTemperatures`] (the four
//!   terminal temperatures, validated), [`lmtd::lmtd`] for the log-mean
//!   temperature difference, and [`lmtd::duty`] for `Q = U * A * LMTD`.
//! - Effectiveness-NTU method ([`ntu`] module): [`NtuProblem`] (the two
//!   heat-capacity rates, `UA`, and the inlet temperatures, validated),
//!   the derived quantities `Cr`, `NTU`, `qmax`, the closed-form
//!   [`ntu::effectiveness`], the duty [`ntu::duty`] = `eps * qmax`, and
//!   a one-call [`ntu::solve`] that also returns the implied outlet
//!   temperatures from each stream's energy balance. The inverse
//!   *sizing* direction is [`ntu::ntu_from_effectiveness`]: the `NTU`
//!   (hence `UA = NTU * Cmin`) needed to reach a target effectiveness.
//!
//! ## Model
//!
//! Let the hot and cold streams have heat-capacity rates `Ch` and `Cc`
//! (each `C = m_dot * c_p`, units W/K), with `Cmin = min(Ch, Cc)`,
//! `Cmax = max(Ch, Cc)`, and ratio `Cr = Cmin / Cmax`.
//!
//! LMTD method — with terminal approaches `dT1`, `dT2` (which depend on
//! the arrangement):
//!
//! ```text
//! LMTD = (dT1 - dT2) / ln(dT1 / dT2)        (LMTD = dT1 when dT1 == dT2)
//! Q    = U * A * LMTD
//! ```
//!
//! Effectiveness-NTU method:
//!
//! ```text
//! NTU  = U * A / Cmin
//! qmax = Cmin * (Th_in - Tc_in)
//! Q    = eps * qmax
//! ```
//!
//! with the closed-form effectiveness relations
//!
//! ```text
//! counterflow:  eps = (1 - E) / (1 - Cr E),  E = exp(-NTU (1 - Cr)),  Cr < 1
//!               eps = NTU / (1 + NTU),                                 Cr = 1
//! parallel:     eps = (1 - exp(-NTU (1 + Cr))) / (1 + Cr)
//! ```
//!
//! These invert in closed form for sizing — given a target `eps`, the
//! required `NTU` is `ln((1 - eps Cr)/(1 - eps))/(1 - Cr)` (counterflow,
//! `Cr < 1`), `eps/(1 - eps)` (counterflow, `Cr = 1`), or
//! `-ln(1 - eps (1 + Cr))/(1 + Cr)` (parallel), defined only below each
//! arrangement's limiting effectiveness (`1` counterflow,
//! `1/(1 + Cr)` parallel).
//!
//! Both reduce to `eps = 1 - exp(-NTU)` at `Cr = 0` (a phase-changing
//! stream, `Cmax -> inf`). For identical terminal temperatures the
//! counterflow arrangement gives the larger LMTD, and for identical
//! `NTU` / `Cr` it gives the larger effectiveness — properties that the
//! crate's tests check against hand calculations.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are the closed-form,
//! well-established heat-transfer relations from standard undergraduate
//! textbooks (Incropera, *Fundamentals of Heat and Mass Transfer*;
//! Cengel, *Heat and Mass Transfer*), implemented as pure deterministic
//! arithmetic with no external processes. It is **not** a clinical /
//! medical or production engineering tool: it carries no fluid-property
//! tables, no fouling factors, no cross-flow / multi-pass correction
//! factors `F`, and no mechanical, pressure-drop, or vibration rating.
//! `U` and the stream capacities are taken as given inputs. Do not use
//! it to size or certify real equipment.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod arrangement;
pub mod error;
pub mod lmtd;
pub mod ntu;

pub use arrangement::FlowArrangement;
pub use error::{ErrorCategory, HeatExchangerError};
pub use lmtd::TerminalTemperatures;
pub use ntu::{NtuProblem, NtuResult};

#[cfg(test)]
mod cross_method_tests {
    //! Cross-checks that the LMTD and effectiveness-NTU methods agree on
    //! the duty of one fully-specified exchanger — the central
    //! self-consistency guarantee tying the two modules together.

    use super::*;

    #[test]
    fn lmtd_and_ntu_agree_on_duty() {
        // Pick capacities, UA and inlets, solve via NTU to get outlets,
        // then feed those outlets to the LMTD method. Both must report
        // the same Q.
        let arr = FlowArrangement::Counterflow;
        let problem = NtuProblem::new(2000.0, 1000.0, 1500.0, 100.0, 20.0).unwrap();
        let result = ntu::solve(&problem, arr).unwrap();

        let temps = TerminalTemperatures::new(
            problem.hot_in,
            result.hot_out,
            problem.cold_in,
            result.cold_out,
        )
        .unwrap();

        // UA / LMTD identity: Q = UA * LMTD, and here UA = 1500.
        let q_lmtd = lmtd::duty(
            problem.ua_w_per_k, // treat the whole UA as U with A = 1
            1.0,
            &temps,
            arr,
        )
        .unwrap();

        assert!(
            (q_lmtd - result.q_w).abs() < 1e-6,
            "LMTD duty {q_lmtd} vs NTU duty {}",
            result.q_w
        );
    }

    #[test]
    fn parallel_lmtd_and_ntu_agree_on_duty() {
        let arr = FlowArrangement::ParallelFlow;
        let problem = NtuProblem::new(1800.0, 1200.0, 1000.0, 90.0, 15.0).unwrap();
        let result = ntu::solve(&problem, arr).unwrap();

        let temps = TerminalTemperatures::new(
            problem.hot_in,
            result.hot_out,
            problem.cold_in,
            result.cold_out,
        )
        .unwrap();

        let q_lmtd = lmtd::duty(problem.ua_w_per_k, 1.0, &temps, arr).unwrap();
        assert!(
            (q_lmtd - result.q_w).abs() < 1e-6,
            "LMTD duty {q_lmtd} vs NTU duty {}",
            result.q_w
        );
    }
}
