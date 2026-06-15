//! Sprocket-pair specification: the two sprockets and the roller-chain
//! pitch that connect a driver shaft to a driven shaft.
//!
//! A [`SprocketPair`] is the geometric heart of a single-stage chain
//! drive. The [`ratio`](SprocketPair::ratio) it exposes — the driven
//! tooth count over the driver tooth count — is the speed-reduction (or
//! step-up) factor every downstream calculation in [`crate::drive`]
//! builds on.

use crate::error::ChainDriveError;
use serde::{Deserialize, Serialize};

/// Practical minimum number of teeth on a roller-chain sprocket.
///
/// Below roughly this many teeth the chordal (polygonal) action of the
/// chain on the sprocket becomes severe and the joint is not a usable
/// power-transmission drive. Used by [`SprocketPair::new`] to reject
/// pathological inputs (a 1- or 2-tooth "sprocket") that would otherwise
/// yield a numerically valid but physically meaningless ratio.
pub const MIN_SPROCKET_TEETH: u32 = 7;

/// A driver sprocket, a driven sprocket, and the common chain pitch.
///
/// "Driver" is the input sprocket fixed to the powered shaft; "driven"
/// is the output sprocket. Both run the **same** roller chain, so they
/// share a single [`pitch_mm`](Self::pitch_mm) — the distance between
/// adjacent roller (pin) centres, measured along the chain.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SprocketPair {
    /// Tooth count of the **driver** (input) sprocket.
    pub driver_teeth: u32,
    /// Tooth count of the **driven** (output) sprocket.
    pub driven_teeth: u32,
    /// Chain pitch in millimetres — the roller-to-roller spacing shared
    /// by both sprockets. Standard ANSI 40 chain, for instance, is
    /// 12.7 mm (1/2 in).
    pub pitch_mm: f64,
}

impl SprocketPair {
    /// Build a validated sprocket pair.
    ///
    /// # Errors
    ///
    /// Returns [`ChainDriveError::TeethTooFew`] if either sprocket has
    /// fewer than [`MIN_SPROCKET_TEETH`] teeth, and
    /// [`ChainDriveError::BadParameter`] if `pitch_mm` is not a finite,
    /// strictly positive number.
    pub fn new(
        driver_teeth: u32,
        driven_teeth: u32,
        pitch_mm: f64,
    ) -> Result<Self, ChainDriveError> {
        if driver_teeth < MIN_SPROCKET_TEETH {
            return Err(ChainDriveError::TeethTooFew {
                name: "driver_teeth",
                got: driver_teeth,
                min: MIN_SPROCKET_TEETH,
            });
        }
        if driven_teeth < MIN_SPROCKET_TEETH {
            return Err(ChainDriveError::TeethTooFew {
                name: "driven_teeth",
                got: driven_teeth,
                min: MIN_SPROCKET_TEETH,
            });
        }
        if !pitch_mm.is_finite() || pitch_mm <= 0.0 {
            return Err(ChainDriveError::bad_parameter(
                "pitch_mm",
                format!("must be a finite value > 0, got {pitch_mm}"),
            ));
        }
        Ok(Self {
            driver_teeth,
            driven_teeth,
            pitch_mm,
        })
    }

    /// Speed-reduction ratio — driven teeth ÷ driver teeth (`N2 / N1`).
    ///
    /// This is the canonical chain / gear ratio. A value `> 1` is a
    /// **reduction** (output turns slower than input but with more
    /// torque); a value `< 1` is a **step-up** (overdrive). Because the
    /// same chain links engage both sprockets, the ratio depends only on
    /// the tooth counts, never on the pitch.
    ///
    /// The tooth counts are validated `>= MIN_SPROCKET_TEETH` at
    /// construction, so the denominator is always non-zero.
    pub fn ratio(&self) -> f64 {
        self.driven_teeth as f64 / self.driver_teeth as f64
    }

    /// Pitch (reference) diameter of a sprocket with `teeth` teeth, in
    /// millimetres.
    ///
    /// The chain rollers seat on a pitch polygon whose circumscribed
    /// circle has diameter `d = p / sin(π / z)`, where `p` is the chain
    /// pitch and `z` the tooth count. This is the diameter at which the
    /// chain-line tension acts, and hence the lever arm that converts
    /// chain tension to shaft torque.
    fn pitch_diameter_mm(&self, teeth: u32) -> f64 {
        self.pitch_mm / (std::f64::consts::PI / teeth as f64).sin()
    }

    /// Pitch diameter of the **driver** sprocket, in millimetres.
    pub fn driver_pitch_diameter_mm(&self) -> f64 {
        self.pitch_diameter_mm(self.driver_teeth)
    }

    /// Pitch diameter of the **driven** sprocket, in millimetres.
    pub fn driven_pitch_diameter_mm(&self) -> f64 {
        self.pitch_diameter_mm(self.driven_teeth)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn pair(n1: u32, n2: u32) -> SprocketPair {
        SprocketPair::new(n1, n2, 12.7).expect("valid pair")
    }

    #[test]
    fn ratio_is_driven_over_driver() {
        // 17-tooth driver, 34-tooth driven -> exactly 2:1 reduction.
        let p = pair(17, 34);
        assert!((p.ratio() - 2.0).abs() < EPS, "ratio = {}", p.ratio());
    }

    #[test]
    fn ratio_is_independent_of_pitch() {
        let a = SprocketPair::new(15, 45, 12.7).unwrap();
        let b = SprocketPair::new(15, 45, 25.4).unwrap();
        assert!((a.ratio() - b.ratio()).abs() < EPS);
        assert!((a.ratio() - 3.0).abs() < EPS);
    }

    #[test]
    fn step_up_ratio_is_below_one() {
        // Big driver, small driven -> overdrive.
        let p = pair(40, 20);
        assert!(p.ratio() < 1.0);
        assert!((p.ratio() - 0.5).abs() < EPS);
    }

    #[test]
    fn pitch_diameter_matches_closed_form() {
        // Independent reference value: d = p / sin(pi/z).
        // For z = 17, p = 12.7 mm: sin(pi/17) = 0.1837495...,
        // d = 12.7 / 0.1837495 = 69.1156... mm.
        let p = pair(17, 34);
        let expected = 12.7 / (std::f64::consts::PI / 17.0).sin();
        assert!((p.driver_pitch_diameter_mm() - expected).abs() < 1e-6);
        assert!((p.driver_pitch_diameter_mm() - 69.115_6).abs() < 1e-3);
    }

    #[test]
    fn larger_sprocket_has_larger_pitch_diameter() {
        let p = pair(17, 34);
        assert!(p.driven_pitch_diameter_mm() > p.driver_pitch_diameter_mm());
    }

    #[test]
    fn rejects_too_few_driver_teeth() {
        let err = SprocketPair::new(3, 20, 12.7).unwrap_err();
        assert_eq!(err.code(), "chaindrive.teeth_too_few");
    }

    #[test]
    fn rejects_too_few_driven_teeth() {
        let err = SprocketPair::new(20, 4, 12.7).unwrap_err();
        assert_eq!(err.code(), "chaindrive.teeth_too_few");
    }

    #[test]
    fn rejects_non_positive_pitch() {
        let err = SprocketPair::new(17, 34, 0.0).unwrap_err();
        assert_eq!(err.code(), "chaindrive.bad_parameter");
        let err = SprocketPair::new(17, 34, -1.0).unwrap_err();
        assert_eq!(err.code(), "chaindrive.bad_parameter");
        let err = SprocketPair::new(17, 34, f64::NAN).unwrap_err();
        assert_eq!(err.code(), "chaindrive.bad_parameter");
    }

    #[test]
    fn serde_round_trip() {
        let p = pair(17, 34);
        let json = serde_json::to_string(&p).expect("serialize");
        let back: SprocketPair = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(p, back);
    }
}
