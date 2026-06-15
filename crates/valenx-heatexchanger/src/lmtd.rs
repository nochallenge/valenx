//! Log-mean-temperature-difference (LMTD) method.
//!
//! For a two-stream exchanger with overall heat-transfer coefficient
//! `U` and area `A`, the duty is
//!
//! `Q = U * A * LMTD`
//!
//! where the log-mean temperature difference combines the two terminal
//! approaches `dT1` and `dT2`:
//!
//! `LMTD = (dT1 - dT2) / ln(dT1 / dT2)`
//!
//! In the limit `dT1 -> dT2` the expression is removable and equals the
//! common approach. The terminal approaches depend on the flow
//! arrangement (see [`crate::FlowArrangement`]).

use serde::{Deserialize, Serialize};

use crate::arrangement::FlowArrangement;
use crate::error::HeatExchangerError;

/// The four terminal temperatures of a two-stream heat exchanger, in a
/// single consistent unit (degrees Celsius or kelvin — only differences
/// enter the formulae, so either works).
///
/// Construct with [`TerminalTemperatures::new`], which validates that
/// the hot stream is hotter than the cold stream at inlet and that each
/// stream's outlet lies on the physically expected side of its inlet
/// (the hot stream cools, the cold stream warms).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TerminalTemperatures {
    /// Hot-stream inlet temperature.
    pub hot_in: f64,
    /// Hot-stream outlet temperature.
    pub hot_out: f64,
    /// Cold-stream inlet temperature.
    pub cold_in: f64,
    /// Cold-stream outlet temperature.
    pub cold_out: f64,
}

impl TerminalTemperatures {
    /// Validate and build a [`TerminalTemperatures`].
    ///
    /// # Errors
    ///
    /// Returns [`HeatExchangerError::InconsistentTemperatures`] if the
    /// hot inlet is not strictly above the cold inlet, if the hot
    /// stream does not cool (`hot_out > hot_in`), or if the cold stream
    /// does not warm (`cold_out < cold_in`).
    pub fn new(
        hot_in: f64,
        hot_out: f64,
        cold_in: f64,
        cold_out: f64,
    ) -> Result<Self, HeatExchangerError> {
        for (name, v) in [
            ("hot_in", hot_in),
            ("hot_out", hot_out),
            ("cold_in", cold_in),
            ("cold_out", cold_out),
        ] {
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
        if hot_out > hot_in {
            return Err(HeatExchangerError::InconsistentTemperatures(format!(
                "hot stream must cool: hot_out ({hot_out}) > hot_in ({hot_in})"
            )));
        }
        if cold_out < cold_in {
            return Err(HeatExchangerError::InconsistentTemperatures(format!(
                "cold stream must warm: cold_out ({cold_out}) < cold_in ({cold_in})"
            )));
        }
        Ok(Self {
            hot_in,
            hot_out,
            cold_in,
            cold_out,
        })
    }

    /// Terminal approaches `(dT1, dT2)` for the given arrangement.
    ///
    /// - Counterflow: `dT1 = hot_in - cold_out`, `dT2 = hot_out -
    ///   cold_in` — the hot inlet is paired with the cold outlet at one
    ///   end and the hot outlet with the cold inlet at the other.
    /// - Parallel flow: `dT1 = hot_in - cold_in`, `dT2 = hot_out -
    ///   cold_out` — both inlets meet at one end, both outlets at the
    ///   other.
    pub fn approaches(&self, arrangement: FlowArrangement) -> (f64, f64) {
        match arrangement {
            FlowArrangement::Counterflow => {
                (self.hot_in - self.cold_out, self.hot_out - self.cold_in)
            }
            FlowArrangement::ParallelFlow => {
                (self.hot_in - self.cold_in, self.hot_out - self.cold_out)
            }
        }
    }
}

/// Relative tolerance below which the two approaches are treated as
/// equal and the removable `dT1 == dT2` limit `LMTD = dT1` is used
/// instead of evaluating the logarithm.
const LMTD_EQUAL_REL_TOL: f64 = 1e-9;

/// Log-mean temperature difference for the given terminal temperatures
/// and flow arrangement.
///
/// Uses the removable limit `LMTD = dT1` when the two approaches are
/// equal to within a small relative tolerance.
///
/// # Errors
///
/// Returns [`HeatExchangerError::Degenerate`] if either approach is
/// non-positive (which would make the logarithm undefined). For valid
/// [`TerminalTemperatures`] this cannot happen for the counterflow case
/// driving a positive duty, but a degenerate parallel-flow case where
/// the streams have already crossed is rejected here.
pub fn lmtd(
    temps: &TerminalTemperatures,
    arrangement: FlowArrangement,
) -> Result<f64, HeatExchangerError> {
    let (dt1, dt2) = temps.approaches(arrangement);
    if dt1 <= 0.0 || dt2 <= 0.0 {
        return Err(HeatExchangerError::Degenerate(format!(
            "both terminal approaches must be positive, got dT1={dt1}, dT2={dt2}"
        )));
    }
    // Removable singularity at dT1 == dT2.
    if (dt1 - dt2).abs() <= LMTD_EQUAL_REL_TOL * dt1.max(dt2) {
        return Ok(0.5 * (dt1 + dt2));
    }
    Ok((dt1 - dt2) / (dt1 / dt2).ln())
}

/// Heat duty from the LMTD method: `Q = U * A * LMTD`.
///
/// `u_w_per_m2k` is the overall heat-transfer coefficient `U`
/// (W/m^2/K) and `area_m2` is the heat-transfer area `A` (m^2). The
/// result is in watts.
///
/// # Errors
///
/// Returns [`HeatExchangerError::BadParameter`] if `U` or `A` is
/// non-positive or non-finite, and propagates any error from [`lmtd`].
pub fn duty(
    u_w_per_m2k: f64,
    area_m2: f64,
    temps: &TerminalTemperatures,
    arrangement: FlowArrangement,
) -> Result<f64, HeatExchangerError> {
    if !(u_w_per_m2k.is_finite() && u_w_per_m2k > 0.0) {
        return Err(HeatExchangerError::BadParameter {
            name: "u_w_per_m2k",
            reason: format!("must be finite and > 0, got {u_w_per_m2k}"),
        });
    }
    if !(area_m2.is_finite() && area_m2 > 0.0) {
        return Err(HeatExchangerError::BadParameter {
            name: "area_m2",
            reason: format!("must be finite and > 0, got {area_m2}"),
        });
    }
    let lmtd_value = lmtd(temps, arrangement)?;
    Ok(u_w_per_m2k * area_m2 * lmtd_value)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> TerminalTemperatures {
        // Hot 120 -> 80, cold 20 -> 60. Symmetric: each stream changes
        // by 40 degrees.
        TerminalTemperatures::new(120.0, 80.0, 20.0, 60.0).unwrap()
    }

    #[test]
    fn rejects_hot_not_above_cold() {
        let err = TerminalTemperatures::new(20.0, 18.0, 25.0, 22.0).unwrap_err();
        assert_eq!(err.code(), "heatexchanger.inconsistent_temperatures");
    }

    #[test]
    fn rejects_hot_stream_warming() {
        assert!(TerminalTemperatures::new(100.0, 110.0, 20.0, 40.0).is_err());
    }

    #[test]
    fn rejects_cold_stream_cooling() {
        assert!(TerminalTemperatures::new(100.0, 80.0, 40.0, 30.0).is_err());
    }

    #[test]
    fn counterflow_approaches_are_correct() {
        let t = sample();
        let (dt1, dt2) = t.approaches(FlowArrangement::Counterflow);
        // dT1 = 120 - 60 = 60 ; dT2 = 80 - 20 = 60.
        assert!((dt1 - 60.0).abs() < 1e-12, "dt1 = {dt1}");
        assert!((dt2 - 60.0).abs() < 1e-12, "dt2 = {dt2}");
    }

    #[test]
    fn parallel_approaches_are_correct() {
        let t = sample();
        let (dt1, dt2) = t.approaches(FlowArrangement::ParallelFlow);
        // dT1 = 120 - 20 = 100 ; dT2 = 80 - 60 = 20.
        assert!((dt1 - 100.0).abs() < 1e-12, "dt1 = {dt1}");
        assert!((dt2 - 20.0).abs() < 1e-12, "dt2 = {dt2}");
    }

    #[test]
    fn equal_approaches_use_removable_limit() {
        // Counterflow on the symmetric sample gives dT1 == dT2 == 60,
        // so LMTD must equal 60 exactly (no 0/0).
        let t = sample();
        let l = lmtd(&t, FlowArrangement::Counterflow).unwrap();
        assert!((l - 60.0).abs() < 1e-9, "lmtd = {l}");
    }

    #[test]
    fn parallel_lmtd_matches_hand_calculation() {
        // dT1 = 100, dT2 = 20 -> LMTD = (100 - 20)/ln(100/20)
        //                              = 80 / ln(5) = 80 / 1.6094379...
        //                              = 49.70679...
        let t = sample();
        let l = lmtd(&t, FlowArrangement::ParallelFlow).unwrap();
        let expected = 80.0 / 5.0_f64.ln();
        assert!((l - expected).abs() < 1e-9, "lmtd = {l}");
        assert!((l - 49.706_790).abs() < 1e-5, "lmtd = {l}");
    }

    #[test]
    fn counterflow_lmtd_at_least_parallel_for_same_terminals() {
        // A documented property: for identical terminal temperatures the
        // counterflow LMTD is >= the parallel-flow LMTD. Use an
        // asymmetric case so the two differ strictly.
        let t = TerminalTemperatures::new(150.0, 90.0, 30.0, 70.0).unwrap();
        let lc = lmtd(&t, FlowArrangement::Counterflow).unwrap();
        let lp = lmtd(&t, FlowArrangement::ParallelFlow).unwrap();
        assert!(lc > lp, "counterflow {lc} should exceed parallel {lp}");
    }

    #[test]
    fn duty_is_product_of_u_a_and_lmtd() {
        let t = sample();
        let q = duty(500.0, 2.0, &t, FlowArrangement::Counterflow).unwrap();
        // LMTD = 60 -> Q = 500 * 2 * 60 = 60_000 W.
        assert!((q - 60_000.0).abs() < 1e-6, "q = {q}");
    }

    #[test]
    fn duty_rejects_nonpositive_area() {
        let t = sample();
        let err = duty(500.0, 0.0, &t, FlowArrangement::Counterflow).unwrap_err();
        assert_eq!(err.code(), "heatexchanger.bad_parameter");
    }

    #[test]
    fn duty_rejects_nonpositive_u() {
        let t = sample();
        assert!(duty(-1.0, 2.0, &t, FlowArrangement::Counterflow).is_err());
    }
}
