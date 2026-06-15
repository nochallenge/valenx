//! A single meshing gear pair (one reduction stage).
//!
//! ## Model
//!
//! For an input pinion with `teeth_in` teeth driving an output gear
//! with `teeth_out` teeth, the **gear ratio** is
//!
//! ```text
//! ratio = teeth_out / teeth_in
//! ```
//!
//! A ratio above `1` is a *reduction* (output turns slower than the
//! input); below `1` is an *overdrive*. Speed and torque transform
//! inversely:
//!
//! ```text
//! output_speed  = input_speed  / ratio
//! output_torque = input_torque * ratio * efficiency
//! ```
//!
//! In the ideal lossless case (`efficiency = 1`) power is conserved:
//! `output_speed * output_torque == input_speed * input_torque`, since
//! the `1/ratio` on speed cancels the `ratio` on torque. A real stage
//! has `efficiency` in `(0, 1)`, which scales the delivered torque (and
//! hence the output power) down while leaving the kinematic speed
//! relationship exact.

use serde::{Deserialize, Serialize};

use crate::error::{check_efficiency, GearboxError};

/// One meshing gear pair: an input pinion driving an output gear,
/// with a mechanical efficiency.
///
/// Construct with [`GearStage::new`] (lossless) or
/// [`GearStage::with_efficiency`]. Both validate their inputs and never
/// panic.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct GearStage {
    teeth_in: u32,
    teeth_out: u32,
    efficiency: f64,
}

impl GearStage {
    /// Build an **ideal (lossless)** stage from the two tooth counts.
    ///
    /// Efficiency defaults to `1.0`. Both tooth counts must be at least
    /// `1`; a count of `0` yields [`GearboxError::ZeroTeeth`].
    pub fn new(teeth_in: u32, teeth_out: u32) -> Result<Self, GearboxError> {
        Self::with_efficiency(teeth_in, teeth_out, 1.0)
    }

    /// Build a stage with an explicit mechanical `efficiency` in the
    /// half-open interval `(0, 1]`.
    ///
    /// # Errors
    ///
    /// - [`GearboxError::ZeroTeeth`] if either tooth count is `0`.
    /// - [`GearboxError::EfficiencyOutOfRange`] if `efficiency` is not
    ///   in `(0, 1]`.
    /// - [`GearboxError::NotFinite`] if `efficiency` is `NaN`/`±∞`.
    pub fn with_efficiency(
        teeth_in: u32,
        teeth_out: u32,
        efficiency: f64,
    ) -> Result<Self, GearboxError> {
        if teeth_in == 0 {
            return Err(GearboxError::ZeroTeeth { name: "teeth_in" });
        }
        if teeth_out == 0 {
            return Err(GearboxError::ZeroTeeth { name: "teeth_out" });
        }
        let efficiency = check_efficiency(efficiency)?;
        Ok(Self {
            teeth_in,
            teeth_out,
            efficiency,
        })
    }

    /// Number of teeth on the input pinion.
    pub fn teeth_in(&self) -> u32 {
        self.teeth_in
    }

    /// Number of teeth on the output gear.
    pub fn teeth_out(&self) -> u32 {
        self.teeth_out
    }

    /// Mechanical efficiency in `(0, 1]`.
    pub fn efficiency(&self) -> f64 {
        self.efficiency
    }

    /// Gear ratio `teeth_out / teeth_in`.
    ///
    /// Greater than `1` for a reduction, less than `1` for an
    /// overdrive, exactly `1` for a 1:1 pass-through.
    pub fn ratio(&self) -> f64 {
        self.teeth_out as f64 / self.teeth_in as f64
    }

    /// Output speed for a given `input_speed` (any consistent angular
    /// unit: rpm, rad/s, ...). Equals `input_speed / ratio`.
    ///
    /// The result carries the same unit as the input; only the
    /// magnitude is scaled, so this is also valid for a signed
    /// (directional) speed.
    pub fn output_speed(&self, input_speed: f64) -> f64 {
        input_speed / self.ratio()
    }

    /// Output torque for a given `input_torque` (any consistent unit:
    /// N·m, lbf·ft, ...). Equals `input_torque * ratio * efficiency`.
    ///
    /// A reduction (`ratio > 1`) multiplies torque; efficiency then
    /// scales the result down to account for friction losses.
    pub fn output_torque(&self, input_torque: f64) -> f64 {
        input_torque * self.ratio() * self.efficiency
    }

    /// Fraction of input *power* delivered to the output.
    ///
    /// Because the kinematic speed ratio and the torque ratio are exact
    /// reciprocals of each other, the only power loss is the efficiency
    /// factor itself; this method simply returns [`Self::efficiency`].
    /// It exists so call sites can read intent (`power_efficiency`)
    /// without re-deriving the cancellation by hand.
    pub fn power_efficiency(&self) -> f64 {
        self.efficiency
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::ErrorCategory;

    const EPS: f64 = 1e-12;

    #[test]
    fn ratio_is_out_over_in() {
        // 40-tooth gear driven by a 10-tooth pinion -> 4:1 reduction.
        let s = GearStage::new(10, 40).unwrap();
        assert!((s.ratio() - 4.0).abs() < EPS);
    }

    #[test]
    fn overdrive_ratio_below_one() {
        // 12-tooth output, 36-tooth input -> 1/3 (overdrive).
        let s = GearStage::new(36, 12).unwrap();
        assert!((s.ratio() - (1.0 / 3.0)).abs() < EPS);
    }

    #[test]
    fn unity_ratio_passthrough() {
        let s = GearStage::new(20, 20).unwrap();
        assert!((s.ratio() - 1.0).abs() < EPS);
        assert!((s.output_speed(1500.0) - 1500.0).abs() < EPS);
        assert!((s.output_torque(50.0) - 50.0).abs() < EPS);
    }

    #[test]
    fn speed_and_torque_are_inverse() {
        // 4:1 reduction: speed quartered, torque quadrupled (ideal).
        let s = GearStage::new(10, 40).unwrap();
        let in_speed = 1000.0;
        let in_torque = 25.0;
        assert!((s.output_speed(in_speed) - 250.0).abs() < EPS);
        assert!((s.output_torque(in_torque) - 100.0).abs() < EPS);
        // Ideal power conservation: P_out == P_in.
        let p_in = in_speed * in_torque;
        let p_out = s.output_speed(in_speed) * s.output_torque(in_torque);
        assert!((p_out - p_in).abs() < 1e-9);
    }

    #[test]
    fn efficiency_reduces_output_torque() {
        let ideal = GearStage::new(10, 40).unwrap();
        let lossy = GearStage::with_efficiency(10, 40, 0.9).unwrap();
        // Same kinematic speed regardless of efficiency.
        assert!((lossy.output_speed(1000.0) - ideal.output_speed(1000.0)).abs() < EPS);
        // Torque scaled by 0.9: 25 * 4 * 0.9 = 90.
        assert!((lossy.output_torque(25.0) - 90.0).abs() < EPS);
        // Lossy torque strictly below ideal.
        assert!(lossy.output_torque(25.0) < ideal.output_torque(25.0));
    }

    #[test]
    fn efficiency_one_is_lossless() {
        let s = GearStage::with_efficiency(15, 45, 1.0).unwrap();
        // 3:1, lossless -> torque exactly tripled.
        assert!((s.output_torque(10.0) - 30.0).abs() < EPS);
        assert!((s.power_efficiency() - 1.0).abs() < EPS);
    }

    #[test]
    fn lossy_power_out_below_in() {
        let s = GearStage::with_efficiency(10, 30, 0.85).unwrap();
        let in_speed = 600.0;
        let in_torque = 12.0;
        let p_in = in_speed * in_torque;
        let p_out = s.output_speed(in_speed) * s.output_torque(in_torque);
        assert!((p_out - 0.85 * p_in).abs() < 1e-9);
        assert!(p_out < p_in);
    }

    #[test]
    fn zero_teeth_in_rejected() {
        let err = GearStage::new(0, 10).unwrap_err();
        assert_eq!(err.code(), "gearbox.zero_teeth");
        assert_eq!(err.category(), ErrorCategory::Input);
        assert!(matches!(
            err,
            GearboxError::ZeroTeeth { name } if name == "teeth_in"
        ));
    }

    #[test]
    fn zero_teeth_out_rejected() {
        let err = GearStage::new(10, 0).unwrap_err();
        assert!(matches!(
            err,
            GearboxError::ZeroTeeth { name } if name == "teeth_out"
        ));
    }

    #[test]
    fn efficiency_above_one_rejected() {
        let err = GearStage::with_efficiency(10, 40, 1.0001).unwrap_err();
        assert_eq!(err.code(), "gearbox.efficiency_out_of_range");
        assert_eq!(err.category(), ErrorCategory::Config);
    }

    #[test]
    fn efficiency_zero_rejected() {
        assert!(GearStage::with_efficiency(10, 40, 0.0).is_err());
    }

    #[test]
    fn efficiency_negative_rejected() {
        assert!(GearStage::with_efficiency(10, 40, -0.5).is_err());
    }

    #[test]
    fn efficiency_nan_rejected() {
        let err = GearStage::with_efficiency(10, 40, f64::NAN).unwrap_err();
        assert_eq!(err.code(), "gearbox.not_finite");
    }

    #[test]
    fn serde_round_trip() {
        let s = GearStage::with_efficiency(11, 37, 0.97).unwrap();
        let json = serde_json::to_string(&s).unwrap();
        let back: GearStage = serde_json::from_str(&json).unwrap();
        assert_eq!(s, back);
    }
}
