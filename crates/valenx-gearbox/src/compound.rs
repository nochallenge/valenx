//! A compound gear train: several [`GearStage`]s in series.
//!
//! ## Model
//!
//! When stages are chained so each stage's output shaft drives the next
//! stage's input, the overall ratio is the **product** of the stage
//! ratios, and the overall efficiency is the product of the stage
//! efficiencies:
//!
//! ```text
//! ratio_total      = ratio_1 * ratio_2 * ... * ratio_n
//! efficiency_total = eff_1   * eff_2   * ... * eff_n
//! ```
//!
//! Speed and torque transform exactly as for a single stage but using
//! these totals:
//!
//! ```text
//! output_speed  = input_speed  / ratio_total
//! output_torque = input_torque * ratio_total * efficiency_total
//! ```
//!
//! This is why multi-stage reducers reach large reductions from modest
//! per-stage ratios: three 4:1 stages give 4 * 4 * 4 = 64:1.

use serde::{Deserialize, Serialize};

use crate::error::GearboxError;
use crate::stage::GearStage;

/// An ordered series of [`GearStage`]s, input shaft first.
///
/// Build with [`CompoundTrain::new`]. The train must contain at least
/// one stage; an empty train has no defined ratio and is rejected with
/// [`GearboxError::EmptyTrain`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CompoundTrain {
    stages: Vec<GearStage>,
}

impl CompoundTrain {
    /// Build a compound train from stages ordered input-shaft first.
    ///
    /// # Errors
    ///
    /// [`GearboxError::EmptyTrain`] if `stages` is empty.
    pub fn new(stages: Vec<GearStage>) -> Result<Self, GearboxError> {
        if stages.is_empty() {
            return Err(GearboxError::EmptyTrain);
        }
        Ok(Self { stages })
    }

    /// The stages, in input-to-output order.
    pub fn stages(&self) -> &[GearStage] {
        &self.stages
    }

    /// Number of stages in the train.
    pub fn len(&self) -> usize {
        self.stages.len()
    }

    /// Always `false`: a [`CompoundTrain`] is guaranteed non-empty by
    /// construction. Provided to satisfy the `clippy::len_without_is_empty`
    /// lint and to read naturally at call sites.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Overall ratio: the product of every stage ratio.
    pub fn ratio(&self) -> f64 {
        self.stages.iter().map(GearStage::ratio).product()
    }

    /// Overall efficiency: the product of every stage efficiency.
    ///
    /// Always in `(0, 1]` since each factor is in `(0, 1]`.
    pub fn efficiency(&self) -> f64 {
        self.stages.iter().map(GearStage::efficiency).product()
    }

    /// Output speed for a given `input_speed`: `input_speed / ratio`.
    pub fn output_speed(&self, input_speed: f64) -> f64 {
        input_speed / self.ratio()
    }

    /// Output torque for a given `input_torque`:
    /// `input_torque * ratio * efficiency`.
    pub fn output_torque(&self, input_torque: f64) -> f64 {
        input_torque * self.ratio() * self.efficiency()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    fn stage(t_in: u32, t_out: u32) -> GearStage {
        GearStage::new(t_in, t_out).unwrap()
    }

    #[test]
    fn ratio_is_product_of_stages() {
        // Three 4:1 stages -> 64:1.
        let t = CompoundTrain::new(vec![stage(10, 40), stage(10, 40), stage(10, 40)]).unwrap();
        assert!((t.ratio() - 64.0).abs() < EPS);
        assert_eq!(t.len(), 3);
        assert!(!t.is_empty());
    }

    #[test]
    fn mixed_ratio_product() {
        // 2:1 then 3:1 then 1.5:1 -> 9:1.
        let t = CompoundTrain::new(vec![stage(10, 20), stage(10, 30), stage(20, 30)]).unwrap();
        assert!((t.ratio() - 9.0).abs() < EPS);
    }

    #[test]
    fn reduction_then_overdrive_cancels() {
        // 4:1 then 1/4 -> exactly 1:1 overall.
        let t = CompoundTrain::new(vec![stage(10, 40), stage(40, 10)]).unwrap();
        assert!((t.ratio() - 1.0).abs() < EPS);
    }

    #[test]
    fn single_stage_train_matches_stage() {
        let s = stage(13, 52);
        let t = CompoundTrain::new(vec![s]).unwrap();
        assert!((t.ratio() - s.ratio()).abs() < EPS);
        assert!((t.output_speed(900.0) - s.output_speed(900.0)).abs() < EPS);
    }

    #[test]
    fn efficiency_is_product_of_stages() {
        let t = CompoundTrain::new(vec![
            GearStage::with_efficiency(10, 40, 0.9).unwrap(),
            GearStage::with_efficiency(10, 40, 0.8).unwrap(),
        ])
        .unwrap();
        // 0.9 * 0.8 = 0.72.
        assert!((t.efficiency() - 0.72).abs() < EPS);
    }

    #[test]
    fn speed_and_torque_use_totals() {
        // 4:1 * 3:1 = 12:1, efficiency 0.95 * 0.95 = 0.9025.
        let t = CompoundTrain::new(vec![
            GearStage::with_efficiency(10, 40, 0.95).unwrap(),
            GearStage::with_efficiency(10, 30, 0.95).unwrap(),
        ])
        .unwrap();
        let ratio = 12.0;
        let eff = 0.95 * 0.95;
        assert!((t.output_speed(1200.0) - 1200.0 / ratio).abs() < 1e-9);
        assert!((t.output_torque(20.0) - 20.0 * ratio * eff).abs() < 1e-9);
    }

    #[test]
    fn lossless_compound_conserves_power() {
        let t = CompoundTrain::new(vec![stage(10, 40), stage(10, 30)]).unwrap();
        let in_speed = 1200.0;
        let in_torque = 15.0;
        let p_in = in_speed * in_torque;
        let p_out = t.output_speed(in_speed) * t.output_torque(in_torque);
        assert!((p_out - p_in).abs() < 1e-9);
    }

    #[test]
    fn empty_train_rejected() {
        let err = CompoundTrain::new(vec![]).unwrap_err();
        assert_eq!(err.code(), "gearbox.empty_train");
        assert!(matches!(err, GearboxError::EmptyTrain));
    }

    #[test]
    fn serde_round_trip() {
        let t = CompoundTrain::new(vec![
            GearStage::with_efficiency(11, 37, 0.97).unwrap(),
            stage(12, 48),
        ])
        .unwrap();
        let json = serde_json::to_string(&t).unwrap();
        let back: CompoundTrain = serde_json::from_str(&json).unwrap();
        assert_eq!(t, back);
    }
}
