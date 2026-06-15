//! # valenx-pidtuning
//!
//! Closed-form PID controller tuning by the Ziegler-Nichols ultimate-gain
//! (closed-loop) method.
//!
//! ## What
//!
//! Turn a single closed-loop experiment into controller settings. The
//! experiment yields an ultimate gain `Ku` (the proportional gain at
//! which the loop first sustains constant-amplitude oscillation) and an
//! ultimate period `Tu` (the period of that oscillation). From that pair
//! this crate emits P, PI, and PID settings in both the standard
//! time-constant form `(Kp, Ti, Td)` and the parallel independent-gain
//! form `(Kp, Ki, Kd)`.
//!
//! Pipeline:
//!
//! - [`UltimateMeasurement::new`] validates the `(Ku, Tu)` pair (finite,
//!   strictly positive).
//! - [`ZieglerNichols`] wraps the measurement and applies the rules.
//! - [`ZieglerNichols::p`] / [`ZieglerNichols::pi`] /
//!   [`ZieglerNichols::pid`] read off each table row as a [`Gains`].
//! - [`Gains::ki`] / [`Gains::kd`] convert to parallel form.
//!
//! ## Model
//!
//! The classic 1942 Ziegler-Nichols closed-loop table for the standard
//! controller `u = Kp ( e + (1/Ti) integral e dt + Td de/dt )`:
//!
//! | Controller | `Kp`       | `Ti`       | `Td`       |
//! | ---------- | ---------- | ---------- | ---------- |
//! | P          | `0.5 Ku`   | infinite   | `0`        |
//! | PI         | `0.45 Ku`  | `Tu / 1.2` | `0`        |
//! | PID        | `0.6 Ku`   | `Tu / 2`   | `Tu / 8`   |
//!
//! Parallel-form gains follow from `Ki = Kp / Ti` and `Kd = Kp Td`. A
//! P controller's infinite `Ti` makes `Ki` collapse to exactly zero; a
//! PI controller's `Td = 0` makes `Kd` zero.
//!
//! ## Honest scope
//!
//! Research/educational grade. This crate evaluates the textbook
//! closed-form Ziegler-Nichols constants and performs the algebraic
//! standard-to-parallel conversion: nothing here identifies a plant,
//! checks a stability margin, models actuator saturation or measurement
//! noise, or simulates a loop. The ultimate-gain rules are deliberately
//! aggressive (they aim for roughly quarter-amplitude decay) and the
//! resulting tunes are lightly damped. This is NOT a clinical/medical or
//! production engineering tuning tool; do not deploy its output to a real
//! control loop without proper plant modelling, robustness analysis, and
//! validation.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod tuning;
pub mod ultimate;

pub use error::{ErrorCategory, PidTuningError};
pub use tuning::{ControllerKind, Gains, ZieglerNichols};
pub use ultimate::UltimateMeasurement;

#[cfg(test)]
mod tests {
    use super::*;

    /// Tolerance for floating-point ground-truth comparisons.
    const EPS: f64 = 1e-12;

    /// Build a Ziegler-Nichols tuner from raw `(Ku, Tu)`, panicking on a
    /// validation failure (tests use only valid inputs).
    fn zn(ku: f64, tu: f64) -> ZieglerNichols {
        ZieglerNichols::new(UltimateMeasurement::new(ku, tu).expect("valid measurement"))
    }

    #[test]
    fn measurement_round_trips_inputs() {
        let m = UltimateMeasurement::new(2.0, 3.0).expect("valid");
        assert!((m.ultimate_gain() - 2.0).abs() < EPS);
        assert!((m.ultimate_period() - 3.0).abs() < EPS);
    }

    #[test]
    fn classic_pid_matches_ground_truth() {
        // Ku = 10, Tu = 4 s. Hand-computed table values:
        //   Kp = 0.6 * 10 = 6.0
        //   Ti = 4 / 2     = 2.0
        //   Td = 4 / 8     = 0.5
        let g = zn(10.0, 4.0).pid();
        assert_eq!(g.kind, ControllerKind::Pid);
        assert!((g.kp() - 6.0).abs() < EPS, "kp was {}", g.kp());
        assert!(
            (g.integral_time() - 2.0).abs() < EPS,
            "ti was {}",
            g.integral_time()
        );
        assert!(
            (g.derivative_time() - 0.5).abs() < EPS,
            "td was {}",
            g.derivative_time()
        );
    }

    #[test]
    fn pid_ti_is_tu_over_two_and_td_is_tu_over_eight() {
        // Independent of Ku, the time constants are pure functions of Tu.
        for &tu in &[0.5_f64, 1.0, 4.0, 7.3, 120.0] {
            let g = zn(3.0, tu).pid();
            assert!(
                (g.integral_time() - tu / 2.0).abs() < EPS,
                "Ti != Tu/2 for tu={tu}"
            );
            assert!(
                (g.derivative_time() - tu / 8.0).abs() < EPS,
                "Td != Tu/8 for tu={tu}"
            );
            // Cross-check the ratio Ti/Td = 4 exactly (Tu cancels).
            assert!(
                (g.integral_time() / g.derivative_time() - 4.0).abs() < EPS,
                "Ti/Td != 4 for tu={tu}"
            );
        }
    }

    #[test]
    fn p_only_is_half_ku() {
        // Kp = 0.5 * Ku across a range of ultimate gains.
        for &ku in &[0.1_f64, 1.0, 2.0, 10.0, 42.5] {
            let g = zn(ku, 5.0).p();
            assert_eq!(g.kind, ControllerKind::P);
            assert!((g.kp() - 0.5 * ku).abs() < EPS, "kp != 0.5*Ku for ku={ku}");
        }
    }

    #[test]
    fn p_only_has_no_integral_or_derivative_action() {
        let g = zn(8.0, 2.0).p();
        assert!(g.integral_time().is_infinite());
        assert!(g.integral_time() > 0.0, "reset time should be +inf");
        assert!((g.derivative_time() - 0.0).abs() < EPS);
        // Parallel form: infinite Ti => Ki collapses to exactly 0.
        assert!((g.ki() - 0.0).abs() < EPS, "ki should be 0, was {}", g.ki());
        assert!((g.kd() - 0.0).abs() < EPS, "kd should be 0, was {}", g.kd());
        assert!(!g.kind.has_integral());
        assert!(!g.kind.has_derivative());
    }

    #[test]
    fn pi_gains_match_ground_truth() {
        // Ku = 20, Tu = 6 s:
        //   Kp = 0.45 * 20 = 9.0
        //   Ti = 6 / 1.2   = 5.0
        //   Td = 0
        let g = zn(20.0, 6.0).pi();
        assert_eq!(g.kind, ControllerKind::Pi);
        assert!((g.kp() - 9.0).abs() < EPS, "kp was {}", g.kp());
        assert!(
            (g.integral_time() - 5.0).abs() < EPS,
            "ti was {}",
            g.integral_time()
        );
        assert!(
            (g.derivative_time() - 0.0).abs() < EPS,
            "td was {}",
            g.derivative_time()
        );
    }

    #[test]
    fn pi_kp_is_point_four_five_ku_and_ti_is_tu_over_1_2() {
        for &(ku, tu) in &[(1.0_f64, 1.0_f64), (5.0, 12.0), (33.0, 0.9)] {
            let g = zn(ku, tu).pi();
            assert!(
                (g.kp() - 0.45 * ku).abs() < EPS,
                "kp != 0.45*Ku for ku={ku}"
            );
            assert!(
                (g.integral_time() - tu / 1.2).abs() < EPS,
                "Ti != Tu/1.2 for tu={tu}"
            );
            assert!(g.derivative_time().abs() < EPS, "Td should be 0");
        }
    }

    #[test]
    fn pi_has_integral_but_no_derivative() {
        let g = zn(4.0, 3.0).pi();
        assert!(g.kind.has_integral());
        assert!(!g.kind.has_derivative());
        assert!(g.integral_time().is_finite());
        // Ki = Kp / Ti = 0.45*4 / (3/1.2) = 1.8 / 2.5 = 0.72.
        assert!((g.ki() - 0.72).abs() < EPS, "ki was {}", g.ki());
        assert!((g.kd() - 0.0).abs() < EPS);
    }

    #[test]
    fn higher_ku_gives_higher_kp_for_every_structure() {
        // Strict monotonicity in Ku at fixed Tu, across P / PI / PID.
        let low = zn(5.0, 4.0);
        let high = zn(15.0, 4.0);
        assert!(high.p().kp() > low.p().kp());
        assert!(high.pi().kp() > low.pi().kp());
        assert!(high.pid().kp() > low.pid().kp());
        // And the time constants do NOT depend on Ku.
        assert!((high.pid().integral_time() - low.pid().integral_time()).abs() < EPS);
        assert!((high.pid().derivative_time() - low.pid().derivative_time()).abs() < EPS);
    }

    #[test]
    fn longer_tu_gives_longer_time_constants() {
        // Strict monotonicity in Tu at fixed Ku for Ti and Td (PID).
        let fast = zn(10.0, 2.0).pid();
        let slow = zn(10.0, 8.0).pid();
        assert!(slow.integral_time() > fast.integral_time());
        assert!(slow.derivative_time() > fast.derivative_time());
        // Kp unchanged when only Tu changes.
        assert!((slow.kp() - fast.kp()).abs() < EPS);
    }

    #[test]
    fn parallel_form_pid_round_trips() {
        // Ki = Kp/Ti, Kd = Kp*Td. For Ku=10, Tu=4: Kp=6, Ti=2, Td=0.5.
        //   Ki = 6 / 2   = 3.0
        //   Kd = 6 * 0.5 = 3.0
        let g = zn(10.0, 4.0).pid();
        assert!((g.ki() - 3.0).abs() < EPS, "ki was {}", g.ki());
        assert!((g.kd() - 3.0).abs() < EPS, "kd was {}", g.kd());
        // Reconstruct standard form from parallel gains.
        let ti_back = g.kp() / g.ki();
        let td_back = g.kd() / g.kp();
        assert!((ti_back - g.integral_time()).abs() < EPS);
        assert!((td_back - g.derivative_time()).abs() < EPS);
    }

    #[test]
    fn z_n_ratios_kp_pid_over_kp_p_is_one_point_two() {
        // PID Kp / P Kp = 0.6Ku / 0.5Ku = 1.2 exactly (Ku cancels).
        let t = zn(7.0, 3.0);
        let ratio = t.pid().kp() / t.p().kp();
        assert!((ratio - 1.2).abs() < EPS, "ratio was {ratio}");
        // PI Kp / P Kp = 0.45 / 0.5 = 0.9.
        let ratio_pi = t.pi().kp() / t.p().kp();
        assert!((ratio_pi - 0.9).abs() < EPS, "ratio_pi was {ratio_pi}");
    }

    #[test]
    fn gains_dispatch_matches_dedicated_accessors() {
        let t = zn(9.0, 5.0);
        assert_eq!(t.gains(ControllerKind::P), t.p());
        assert_eq!(t.gains(ControllerKind::Pi), t.pi());
        assert_eq!(t.gains(ControllerKind::Pid), t.pid());
    }

    #[test]
    fn measurement_is_accessible_from_tuner() {
        let t = zn(2.5, 1.5);
        assert!((t.measurement().ultimate_gain() - 2.5).abs() < EPS);
        assert!((t.measurement().ultimate_period() - 1.5).abs() < EPS);
    }

    #[test]
    fn rejects_non_positive_gain() {
        for bad in [0.0, -1.0, -0.001] {
            let err = UltimateMeasurement::new(bad, 4.0).expect_err("should reject Ku");
            match err {
                PidTuningError::NonPositive { name, value } => {
                    assert_eq!(name, "Ku");
                    assert!((value - bad).abs() < EPS);
                }
            }
            assert_eq!(err.code(), "pidtuning.non-positive");
            assert_eq!(err.category(), ErrorCategory::Input);
        }
    }

    #[test]
    fn rejects_non_positive_period() {
        let err = UltimateMeasurement::new(4.0, 0.0).expect_err("should reject Tu");
        match err {
            PidTuningError::NonPositive { name, value } => {
                assert_eq!(name, "Tu");
                assert!(value.abs() < EPS);
            }
        }
    }

    #[test]
    fn rejects_non_finite_inputs() {
        assert!(UltimateMeasurement::new(f64::NAN, 1.0).is_err());
        assert!(UltimateMeasurement::new(f64::INFINITY, 1.0).is_err());
        assert!(UltimateMeasurement::new(1.0, f64::NAN).is_err());
        assert!(UltimateMeasurement::new(1.0, f64::INFINITY).is_err());
    }

    #[test]
    fn serde_round_trips_gains() {
        let g = zn(10.0, 4.0).pid();
        let json = serde_json::to_string(&g).expect("serialize");
        let back: Gains = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(g, back);
    }
}
