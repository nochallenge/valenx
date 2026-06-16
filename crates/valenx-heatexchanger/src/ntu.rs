//! Effectiveness-NTU method.
//!
//! The effectiveness-NTU (number of transfer units) method sizes or
//! rates a heat exchanger from the two stream heat-capacity rates
//! `C = m_dot * c_p` (W/K), the overall conductance `UA` (W/K) and the
//! two inlet temperatures — without needing the outlet temperatures up
//! front.
//!
//! Definitions (Incropera/Cengel):
//!
//! - `Cmin = min(Ch, Cc)`, `Cmax = max(Ch, Cc)`
//! - heat-capacity-rate ratio `Cr = Cmin / Cmax` in `[0, 1]`
//! - `NTU = UA / Cmin`
//! - maximum possible duty `qmax = Cmin * (Th_in - Tc_in)`
//! - effectiveness `eps = Q / qmax`, so `Q = eps * qmax`
//!
//! Closed-form effectiveness relations (the two modelled here):
//!
//! Counterflow:
//!
//! `eps = (1 - exp(-NTU (1 - Cr))) / (1 - Cr exp(-NTU (1 - Cr)))` for
//! `Cr < 1`, and the removable limit `eps = NTU / (1 + NTU)` for
//! `Cr = 1`.
//!
//! Parallel flow:
//!
//! `eps = (1 - exp(-NTU (1 + Cr))) / (1 + Cr)`.
//!
//! Both reduce, at `Cr = 0` (one stream condensing/evaporating, so
//! `Cmax -> inf`), to `eps = 1 - exp(-NTU)`.

use serde::{Deserialize, Serialize};

use crate::arrangement::FlowArrangement;
use crate::error::HeatExchangerError;

/// A validated effectiveness-NTU problem: the two heat-capacity rates,
/// the overall conductance `UA`, and the two inlet temperatures.
///
/// Construct with [`NtuProblem::new`]. The stored `c_hot` / `c_cold`
/// are the stream heat-capacity rates `C = m_dot * c_p` in W/K; a
/// condensing or evaporating stream (effectively infinite capacity) is
/// modelled by passing [`f64::INFINITY`] for that stream, which yields
/// `Cr = 0`.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NtuProblem {
    /// Hot-stream heat-capacity rate `Ch = m_dot * c_p` (W/K).
    pub c_hot: f64,
    /// Cold-stream heat-capacity rate `Cc = m_dot * c_p` (W/K).
    pub c_cold: f64,
    /// Overall conductance `UA = U * A` (W/K).
    pub ua_w_per_k: f64,
    /// Hot-stream inlet temperature.
    pub hot_in: f64,
    /// Cold-stream inlet temperature.
    pub cold_in: f64,
}

impl NtuProblem {
    /// Validate and build an [`NtuProblem`].
    ///
    /// # Errors
    ///
    /// Returns [`HeatExchangerError::BadParameter`] if either capacity
    /// rate is non-positive or NaN, if at least one of the two
    /// capacities is not finite (both infinite is rejected — `Cmin`
    /// would be infinite), or if `UA` is negative / non-finite.
    /// Returns [`HeatExchangerError::InconsistentTemperatures`] if the
    /// hot inlet does not strictly exceed the cold inlet.
    pub fn new(
        c_hot: f64,
        c_cold: f64,
        ua_w_per_k: f64,
        hot_in: f64,
        cold_in: f64,
    ) -> Result<Self, HeatExchangerError> {
        for (name, c) in [("c_hot", c_hot), ("c_cold", c_cold)] {
            if c.is_nan() || c <= 0.0 {
                return Err(HeatExchangerError::BadParameter {
                    name,
                    reason: format!("heat-capacity rate must be > 0, got {c}"),
                });
            }
        }
        if !c_hot.is_finite() && !c_cold.is_finite() {
            return Err(HeatExchangerError::BadParameter {
                name: "c_hot",
                reason: "at most one stream may have infinite capacity".to_string(),
            });
        }
        if !(ua_w_per_k.is_finite() && ua_w_per_k >= 0.0) {
            return Err(HeatExchangerError::BadParameter {
                name: "ua_w_per_k",
                reason: format!("must be finite and >= 0, got {ua_w_per_k}"),
            });
        }
        for (name, v) in [("hot_in", hot_in), ("cold_in", cold_in)] {
            if !v.is_finite() {
                return Err(HeatExchangerError::BadParameter {
                    name,
                    reason: format!("must be finite, got {v}"),
                });
            }
        }
        if hot_in <= cold_in {
            return Err(HeatExchangerError::InconsistentTemperatures(format!(
                "hot inlet ({hot_in}) must exceed cold inlet ({cold_in})"
            )));
        }
        Ok(Self {
            c_hot,
            c_cold,
            ua_w_per_k,
            hot_in,
            cold_in,
        })
    }

    /// Minimum heat-capacity rate `Cmin = min(Ch, Cc)` (W/K).
    pub fn c_min(&self) -> f64 {
        self.c_hot.min(self.c_cold)
    }

    /// Maximum heat-capacity rate `Cmax = max(Ch, Cc)` (W/K). May be
    /// [`f64::INFINITY`] for a phase-changing stream.
    pub fn c_max(&self) -> f64 {
        self.c_hot.max(self.c_cold)
    }

    /// Heat-capacity-rate ratio `Cr = Cmin / Cmax` in `[0, 1]`. Equals
    /// `0` when one stream has infinite capacity.
    pub fn capacity_ratio(&self) -> f64 {
        let cmax = self.c_max();
        if cmax.is_infinite() {
            0.0
        } else {
            self.c_min() / cmax
        }
    }

    /// Number of transfer units `NTU = UA / Cmin` (dimensionless).
    pub fn ntu(&self) -> f64 {
        self.ua_w_per_k / self.c_min()
    }

    /// Maximum possible duty `qmax = Cmin * (Th_in - Tc_in)` (W).
    pub fn q_max(&self) -> f64 {
        self.c_min() * (self.hot_in - self.cold_in)
    }
}

/// Effectiveness `eps` from `NTU`, the capacity ratio `cr`, and the
/// flow arrangement, evaluated directly from the closed-form relations.
///
/// This is the lower-level entry point; [`NtuProblem`] callers normally
/// use [`effectiveness`]. The result is clamped into `[0, 1]` to absorb
/// floating-point round-off at the extremes (it is analytically already
/// in that range).
///
/// # Errors
///
/// Returns [`HeatExchangerError::BadParameter`] if `ntu` is negative or
/// non-finite, or if `cr` is outside `[0, 1]`.
pub fn effectiveness_from_ntu(
    ntu: f64,
    cr: f64,
    arrangement: FlowArrangement,
) -> Result<f64, HeatExchangerError> {
    if !(ntu.is_finite() && ntu >= 0.0) {
        return Err(HeatExchangerError::BadParameter {
            name: "ntu",
            reason: format!("must be finite and >= 0, got {ntu}"),
        });
    }
    if !(cr.is_finite() && (0.0..=1.0).contains(&cr)) {
        return Err(HeatExchangerError::BadParameter {
            name: "cr",
            reason: format!("capacity ratio must be in [0, 1], got {cr}"),
        });
    }

    let eps = match arrangement {
        FlowArrangement::Counterflow => counterflow_effectiveness(ntu, cr),
        FlowArrangement::ParallelFlow => parallel_effectiveness(ntu, cr),
    };
    Ok(eps.clamp(0.0, 1.0))
}

/// Counterflow effectiveness with the `Cr -> 1` limit handled.
fn counterflow_effectiveness(ntu: f64, cr: f64) -> f64 {
    // Cr == 1 removable limit: eps = NTU / (1 + NTU).
    if (1.0 - cr).abs() <= 1e-12 {
        return ntu / (1.0 + ntu);
    }
    let e = (-ntu * (1.0 - cr)).exp();
    (1.0 - e) / (1.0 - cr * e)
}

/// Parallel-flow effectiveness.
fn parallel_effectiveness(ntu: f64, cr: f64) -> f64 {
    (1.0 - (-ntu * (1.0 + cr)).exp()) / (1.0 + cr)
}

/// The `NTU` required to reach a target effectiveness `eps` at capacity
/// ratio `cr` and flow arrangement — the exact inverse of
/// [`effectiveness_from_ntu`], i.e. the *sizing* direction of the
/// effectiveness-NTU method (the required conductance is then
/// `UA = NTU * Cmin`).
///
/// The closed-form inversions (Incropera, Table 11.5) are
///
/// ```text
/// counterflow, Cr < 1:  NTU = ln((1 - eps Cr) / (1 - eps)) / (1 - Cr)
/// counterflow, Cr = 1:  NTU = eps / (1 - eps)
/// parallel:             NTU = -ln(1 - eps (1 + Cr)) / (1 + Cr)
/// ```
///
/// (the counterflow `Cr < 1` form already covers `Cr = 0`, giving the
/// shared `NTU = -ln(1 - eps)`). A finite `NTU` only exists below the
/// arrangement's limiting effectiveness — `1` for counterflow and
/// `1 / (1 + Cr)` for parallel flow — so a target at or above that limit
/// is rejected rather than returning an infinity.
///
/// # Errors
///
/// Returns [`HeatExchangerError::BadParameter`] if `eps` is negative or
/// non-finite, if `cr` is outside `[0, 1]`, or if `eps` is at or above
/// the maximum effectiveness the arrangement can reach as `NTU -> inf`.
pub fn ntu_from_effectiveness(
    eps: f64,
    cr: f64,
    arrangement: FlowArrangement,
) -> Result<f64, HeatExchangerError> {
    if !(eps.is_finite() && eps >= 0.0) {
        return Err(HeatExchangerError::BadParameter {
            name: "eps",
            reason: format!("effectiveness must be finite and >= 0, got {eps}"),
        });
    }
    if !(cr.is_finite() && (0.0..=1.0).contains(&cr)) {
        return Err(HeatExchangerError::BadParameter {
            name: "cr",
            reason: format!("capacity ratio must be in [0, 1], got {cr}"),
        });
    }
    let eps_max = match arrangement {
        FlowArrangement::Counterflow => 1.0,
        FlowArrangement::ParallelFlow => 1.0 / (1.0 + cr),
    };
    if eps >= eps_max {
        return Err(HeatExchangerError::BadParameter {
            name: "eps",
            reason: format!(
                "effectiveness {eps} is unreachable for {arrangement:?} at Cr = {cr}; \
                 the limit as NTU -> inf is {eps_max}"
            ),
        });
    }
    let ntu = match arrangement {
        FlowArrangement::Counterflow => {
            if (1.0 - cr).abs() <= 1e-12 {
                eps / (1.0 - eps)
            } else {
                ((1.0 - eps * cr) / (1.0 - eps)).ln() / (1.0 - cr)
            }
        }
        FlowArrangement::ParallelFlow => -(1.0 - eps * (1.0 + cr)).ln() / (1.0 + cr),
    };
    Ok(ntu)
}

/// Effectiveness for a validated [`NtuProblem`].
///
/// # Errors
///
/// Propagates any error from [`effectiveness_from_ntu`] (which cannot
/// occur for a well-formed [`NtuProblem`], whose `NTU >= 0` and
/// `Cr in [0, 1]` are guaranteed by construction).
pub fn effectiveness(
    problem: &NtuProblem,
    arrangement: FlowArrangement,
) -> Result<f64, HeatExchangerError> {
    effectiveness_from_ntu(problem.ntu(), problem.capacity_ratio(), arrangement)
}

/// Actual heat duty `Q = eps * qmax` for the problem and arrangement.
///
/// # Errors
///
/// Propagates any error from [`effectiveness`].
pub fn duty(problem: &NtuProblem, arrangement: FlowArrangement) -> Result<f64, HeatExchangerError> {
    let eps = effectiveness(problem, arrangement)?;
    Ok(eps * problem.q_max())
}

/// Convenience bundle of the derived quantities for a solved problem.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct NtuResult {
    /// Capacity ratio `Cr` in `[0, 1]`.
    pub capacity_ratio: f64,
    /// Number of transfer units `NTU`.
    pub ntu: f64,
    /// Effectiveness `eps` in `[0, 1]`.
    pub effectiveness: f64,
    /// Maximum possible duty `qmax` (W).
    pub q_max_w: f64,
    /// Actual duty `Q = eps * qmax` (W).
    pub q_w: f64,
    /// Hot-stream outlet temperature `Th_out = Th_in - Q / Ch`.
    pub hot_out: f64,
    /// Cold-stream outlet temperature `Tc_out = Tc_in + Q / Cc`.
    pub cold_out: f64,
}

/// Solve an [`NtuProblem`] into an [`NtuResult`], including the implied
/// outlet temperatures from energy balance on each stream.
///
/// The outlet temperatures use `Q = Ch (Th_in - Th_out) = Cc (Tc_out -
/// Tc_in)`; for a phase-changing (infinite-capacity) stream the
/// corresponding outlet equals its inlet, as the limit `Q / inf -> 0`
/// gives.
///
/// # Errors
///
/// Propagates any error from [`effectiveness`].
pub fn solve(
    problem: &NtuProblem,
    arrangement: FlowArrangement,
) -> Result<NtuResult, HeatExchangerError> {
    let eps = effectiveness(problem, arrangement)?;
    let q_max = problem.q_max();
    let q = eps * q_max;
    // Q / C -> 0 for an infinite-capacity stream, so guard the division.
    let hot_drop = if problem.c_hot.is_infinite() {
        0.0
    } else {
        q / problem.c_hot
    };
    let cold_rise = if problem.c_cold.is_infinite() {
        0.0
    } else {
        q / problem.c_cold
    };
    Ok(NtuResult {
        capacity_ratio: problem.capacity_ratio(),
        ntu: problem.ntu(),
        effectiveness: eps,
        q_max_w: q_max,
        q_w: q,
        hot_out: problem.hot_in - hot_drop,
        cold_out: problem.cold_in + cold_rise,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derived_quantities_are_correct() {
        // Ch = 2000, Cc = 1000 -> Cmin = 1000, Cr = 0.5.
        // UA = 1500 -> NTU = 1.5. qmax = 1000 * (100 - 20) = 80_000 W.
        let p = NtuProblem::new(2000.0, 1000.0, 1500.0, 100.0, 20.0).unwrap();
        assert!((p.c_min() - 1000.0).abs() < 1e-9);
        assert!((p.c_max() - 2000.0).abs() < 1e-9);
        assert!((p.capacity_ratio() - 0.5).abs() < 1e-12);
        assert!((p.ntu() - 1.5).abs() < 1e-12);
        assert!((p.q_max() - 80_000.0).abs() < 1e-6);
    }

    #[test]
    fn effectiveness_is_within_unit_interval() {
        // Sweep NTU and Cr for both arrangements; eps must stay in
        // [0, 1] everywhere.
        for arr in [FlowArrangement::Counterflow, FlowArrangement::ParallelFlow] {
            let mut ntu = 0.0;
            while ntu <= 12.0 {
                let mut cr = 0.0;
                while cr <= 1.0 {
                    let eps = effectiveness_from_ntu(ntu, cr, arr).unwrap();
                    assert!(
                        (0.0..=1.0).contains(&eps),
                        "eps={eps} out of range at ntu={ntu}, cr={cr}, {arr:?}"
                    );
                    cr += 0.1;
                }
                ntu += 0.25;
            }
        }
    }

    #[test]
    fn cr_zero_reduces_to_one_minus_exp_neg_ntu() {
        // Validation requirement: Cr = 0 gives eps = 1 - exp(-NTU) for
        // BOTH arrangements.
        for arr in [FlowArrangement::Counterflow, FlowArrangement::ParallelFlow] {
            for &ntu in &[0.0, 0.5, 1.0, 2.0, 3.5, 8.0] {
                let eps = effectiveness_from_ntu(ntu, 0.0, arr).unwrap();
                let expected = 1.0 - (-ntu).exp();
                assert!(
                    (eps - expected).abs() < 1e-12,
                    "eps={eps} expected={expected} at ntu={ntu}, {arr:?}"
                );
            }
        }
    }

    #[test]
    fn infinite_capacity_stream_gives_cr_zero_branch() {
        // A condensing hot stream (Ch = inf) -> Cr = 0, so the duty
        // follows eps = 1 - exp(-NTU) with Cmin = Cc.
        let p = NtuProblem::new(f64::INFINITY, 1000.0, 2000.0, 150.0, 30.0).unwrap();
        assert!((p.capacity_ratio() - 0.0).abs() < 1e-15);
        assert!((p.c_min() - 1000.0).abs() < 1e-9);
        let ntu = p.ntu(); // 2000 / 1000 = 2.
        assert!((ntu - 2.0).abs() < 1e-12);
        let eps = effectiveness(&p, FlowArrangement::Counterflow).unwrap();
        assert!((eps - (1.0 - (-2.0_f64).exp())).abs() < 1e-12, "eps={eps}");
    }

    #[test]
    fn counterflow_eps_approaches_one_as_ntu_grows_when_cr_below_one() {
        // NTU -> inf gives eps -> 1 for Cr < 1 (counterflow).
        let eps = effectiveness_from_ntu(50.0, 0.6, FlowArrangement::Counterflow).unwrap();
        assert!(eps > 0.999_999, "eps={eps}");
        assert!(eps <= 1.0, "eps={eps}");
    }

    #[test]
    fn counterflow_cr_one_limit_is_finite() {
        // For Cr = 1 the naive formula is 0/0; the limit eps =
        // NTU/(1+NTU) must be used and stay < 1.
        let ntu = 4.0;
        let eps = effectiveness_from_ntu(ntu, 1.0, FlowArrangement::Counterflow).unwrap();
        let expected = ntu / (1.0 + ntu); // 0.8
        assert!((eps - expected).abs() < 1e-12, "eps={eps}");
        assert!(eps < 1.0, "eps={eps}");
    }

    #[test]
    fn counterflow_matches_textbook_value() {
        // Incropera worked-value class case: NTU = 1.5, Cr = 0.5.
        // e = exp(-1.5 * 0.5) = exp(-0.75) = 0.4723665...
        // eps = (1 - e) / (1 - 0.5 e)
        //     = (0.5276335) / (0.7638168) = 0.690783...
        let eps = effectiveness_from_ntu(1.5, 0.5, FlowArrangement::Counterflow).unwrap();
        let e = (-0.75_f64).exp();
        let expected = (1.0 - e) / (1.0 - 0.5 * e);
        assert!((eps - expected).abs() < 1e-12, "eps={eps}");
        assert!((eps - 0.690_783).abs() < 1e-5, "eps={eps}");
    }

    #[test]
    fn parallel_matches_closed_form() {
        // NTU = 1.5, Cr = 0.5 -> eps = (1 - exp(-1.5*1.5)) / 1.5
        //                            = (1 - exp(-2.25)) / 1.5
        //                            = (1 - 0.105399) / 1.5
        //                            = 0.596400...
        let eps = effectiveness_from_ntu(1.5, 0.5, FlowArrangement::ParallelFlow).unwrap();
        let expected = (1.0 - (-2.25_f64).exp()) / 1.5;
        assert!((eps - expected).abs() < 1e-12, "eps={eps}");
        assert!((eps - 0.596_400).abs() < 1e-5, "eps={eps}");
    }

    #[test]
    fn counterflow_effectiveness_at_least_parallel() {
        // Documented property: counterflow eps >= parallel eps for the
        // same NTU and Cr (strict when Cr > 0 and NTU > 0).
        for &ntu in &[0.5, 1.0, 2.0, 4.0] {
            for &cr in &[0.2, 0.5, 0.8, 1.0] {
                let ec = effectiveness_from_ntu(ntu, cr, FlowArrangement::Counterflow).unwrap();
                let ep = effectiveness_from_ntu(ntu, cr, FlowArrangement::ParallelFlow).unwrap();
                assert!(
                    ec >= ep - 1e-12,
                    "counterflow {ec} < parallel {ep} at ntu={ntu}, cr={cr}"
                );
            }
        }
    }

    #[test]
    fn duty_is_eps_times_qmax() {
        let p = NtuProblem::new(2000.0, 1000.0, 1500.0, 100.0, 20.0).unwrap();
        let eps = effectiveness(&p, FlowArrangement::Counterflow).unwrap();
        let q = duty(&p, FlowArrangement::Counterflow).unwrap();
        assert!((q - eps * p.q_max()).abs() < 1e-6, "q={q}");
    }

    #[test]
    fn solve_outlets_satisfy_energy_balance() {
        // Q must equal Ch*(Th_in - Th_out) and Cc*(Tc_out - Tc_in).
        let p = NtuProblem::new(2000.0, 1000.0, 1500.0, 100.0, 20.0).unwrap();
        let r = solve(&p, FlowArrangement::Counterflow).unwrap();
        let q_from_hot = p.c_hot * (p.hot_in - r.hot_out);
        let q_from_cold = p.c_cold * (r.cold_out - p.cold_in);
        assert!(
            (q_from_hot - r.q_w).abs() < 1e-6,
            "hot balance {q_from_hot}"
        );
        assert!(
            (q_from_cold - r.q_w).abs() < 1e-6,
            "cold balance {q_from_cold}"
        );
        // Each stream must stay within the inlet temperature span. Note
        // that in counterflow the cold OUTLET may legitimately exceed
        // the hot OUTLET (the defining advantage of counterflow), so the
        // physical bounds are per-stream, not a hot_out > cold_out
        // ordering.
        assert!(
            r.hot_out >= p.cold_in,
            "hot_out {} below cold_in",
            r.hot_out
        );
        assert!(r.hot_out <= p.hot_in, "hot_out {} above hot_in", r.hot_out);
        assert!(
            r.cold_out >= p.cold_in,
            "cold_out {} below cold_in",
            r.cold_out
        );
        assert!(
            r.cold_out <= p.hot_in,
            "cold_out {} above hot_in",
            r.cold_out
        );
        // Here the cold outlet does exceed the hot outlet — assert it so
        // the counterflow behaviour is documented as intended.
        assert!(
            r.cold_out > r.hot_out,
            "expected cold_out {} > hot_out {} for this counterflow case",
            r.cold_out,
            r.hot_out
        );
    }

    #[test]
    fn solve_infinite_hot_stream_keeps_hot_outlet_at_inlet() {
        // Condensing hot stream: its temperature does not drop.
        let p = NtuProblem::new(f64::INFINITY, 1000.0, 2000.0, 150.0, 30.0).unwrap();
        let r = solve(&p, FlowArrangement::Counterflow).unwrap();
        assert!((r.hot_out - 150.0).abs() < 1e-9, "hot_out={}", r.hot_out);
        // Cold stream picks up all the duty.
        let q_from_cold = p.c_cold * (r.cold_out - p.cold_in);
        assert!((q_from_cold - r.q_w).abs() < 1e-6);
    }

    #[test]
    fn rejects_nonpositive_capacity() {
        assert!(NtuProblem::new(0.0, 1000.0, 1500.0, 100.0, 20.0).is_err());
        assert!(NtuProblem::new(2000.0, -5.0, 1500.0, 100.0, 20.0).is_err());
    }

    #[test]
    fn rejects_both_infinite_capacities() {
        let err = NtuProblem::new(f64::INFINITY, f64::INFINITY, 1500.0, 100.0, 20.0).unwrap_err();
        assert_eq!(err.code(), "heatexchanger.bad_parameter");
    }

    #[test]
    fn rejects_hot_inlet_not_above_cold() {
        let err = NtuProblem::new(2000.0, 1000.0, 1500.0, 20.0, 50.0).unwrap_err();
        assert_eq!(err.code(), "heatexchanger.inconsistent_temperatures");
    }

    #[test]
    fn effectiveness_from_ntu_rejects_bad_cr() {
        assert!(effectiveness_from_ntu(1.0, 1.5, FlowArrangement::Counterflow).is_err());
        assert!(effectiveness_from_ntu(1.0, -0.1, FlowArrangement::Counterflow).is_err());
    }

    #[test]
    fn effectiveness_from_ntu_rejects_negative_ntu() {
        assert!(effectiveness_from_ntu(-1.0, 0.5, FlowArrangement::ParallelFlow).is_err());
    }

    #[test]
    fn zero_ntu_gives_zero_effectiveness() {
        for arr in [FlowArrangement::Counterflow, FlowArrangement::ParallelFlow] {
            for &cr in &[0.0, 0.5, 1.0] {
                let eps = effectiveness_from_ntu(0.0, cr, arr).unwrap();
                assert!(eps.abs() < 1e-12, "eps={eps} at cr={cr}, {arr:?}");
            }
        }
    }

    #[test]
    fn ntu_from_effectiveness_inverts_effectiveness_from_ntu() {
        // GOLD identity: forward then inverse recovers NTU exactly, swept
        // over both arrangements, a range of NTU, and Cr including the
        // 0 and 1 limits.
        for arr in [FlowArrangement::Counterflow, FlowArrangement::ParallelFlow] {
            for &cr in &[0.0, 0.3, 0.6, 1.0] {
                for &ntu in &[0.25, 0.5, 1.0, 2.0, 4.0, 7.0] {
                    let eps = effectiveness_from_ntu(ntu, cr, arr).unwrap();
                    let back = ntu_from_effectiveness(eps, cr, arr).unwrap();
                    assert!(
                        (back - ntu).abs() < 1e-9 * ntu.max(1.0),
                        "round-trip NTU {back} vs {ntu} at cr={cr}, {arr:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn ntu_from_effectiveness_matches_closed_forms() {
        // Counterflow NTU = 1.5, Cr = 0.5 -> eps = 0.690783..., inverts
        // back to 1.5.
        let arr = FlowArrangement::Counterflow;
        let eps = effectiveness_from_ntu(1.5, 0.5, arr).unwrap();
        assert!((ntu_from_effectiveness(eps, 0.5, arr).unwrap() - 1.5).abs() < 1e-9);

        // Parallel NTU = 1.5, Cr = 0.5 -> eps = 0.596400..., inverts to 1.5.
        let arr = FlowArrangement::ParallelFlow;
        let epsp = effectiveness_from_ntu(1.5, 0.5, arr).unwrap();
        assert!((ntu_from_effectiveness(epsp, 0.5, arr).unwrap() - 1.5).abs() < 1e-9);
    }

    #[test]
    fn cr_zero_inverse_is_minus_ln_one_minus_eps() {
        // At Cr = 0 both arrangements invert to NTU = -ln(1 - eps).
        for arr in [FlowArrangement::Counterflow, FlowArrangement::ParallelFlow] {
            for &eps in &[0.1, 0.5, 0.9, 0.99] {
                let ntu = ntu_from_effectiveness(eps, 0.0, arr).unwrap();
                let expected = -(1.0 - eps).ln();
                assert!(
                    (ntu - expected).abs() < 1e-12,
                    "ntu={ntu} at eps={eps}, {arr:?}"
                );
            }
        }
    }

    #[test]
    fn counterflow_cr_one_inverse_is_eps_over_one_minus_eps() {
        // Cr = 1 counterflow: eps = NTU/(1+NTU) inverts to NTU = eps/(1-eps).
        let arr = FlowArrangement::Counterflow;
        let ntu = ntu_from_effectiveness(0.8, 1.0, arr).unwrap();
        assert!((ntu - 4.0).abs() < 1e-12, "ntu={ntu}"); // 0.8 / 0.2
    }

    #[test]
    fn zero_effectiveness_needs_zero_ntu() {
        for arr in [FlowArrangement::Counterflow, FlowArrangement::ParallelFlow] {
            for &cr in &[0.0, 0.5, 1.0] {
                let ntu = ntu_from_effectiveness(0.0, cr, arr).unwrap();
                assert!(ntu.abs() < 1e-12, "ntu={ntu} at cr={cr}, {arr:?}");
            }
        }
    }

    #[test]
    fn inverse_matches_a_solved_problem() {
        // Solve a real problem for its effectiveness, then confirm the
        // inverse recovers the problem's own NTU.
        for arr in [FlowArrangement::Counterflow, FlowArrangement::ParallelFlow] {
            let p = NtuProblem::new(2000.0, 1000.0, 1500.0, 100.0, 20.0).unwrap();
            let r = solve(&p, arr).unwrap();
            let ntu = ntu_from_effectiveness(r.effectiveness, p.capacity_ratio(), arr).unwrap();
            assert!((ntu - p.ntu()).abs() < 1e-9, "ntu={ntu} vs {}", p.ntu());
        }
    }

    #[test]
    fn inverse_rejects_unreachable_effectiveness() {
        // Parallel limit is 1/(1+Cr): at Cr = 0.5 that is 0.6667, so 0.8
        // is unreachable.
        assert!(ntu_from_effectiveness(0.8, 0.5, FlowArrangement::ParallelFlow).is_err());
        // Counterflow can approach 1 but never reach it.
        assert!(ntu_from_effectiveness(1.0, 0.5, FlowArrangement::Counterflow).is_err());
        // Just below the parallel limit is fine.
        assert!(ntu_from_effectiveness(0.66, 0.5, FlowArrangement::ParallelFlow).is_ok());
    }

    #[test]
    fn inverse_rejects_bad_domain() {
        assert!(ntu_from_effectiveness(-0.1, 0.5, FlowArrangement::Counterflow).is_err());
        assert!(ntu_from_effectiveness(f64::NAN, 0.5, FlowArrangement::Counterflow).is_err());
        assert!(ntu_from_effectiveness(0.5, 1.5, FlowArrangement::Counterflow).is_err());
        assert!(ntu_from_effectiveness(0.5, -0.1, FlowArrangement::ParallelFlow).is_err());
    }
}
