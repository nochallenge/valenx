//! GD&T verification — check an actual measured value against an ASME
//! Y14.5 feature-control frame.
//!
//! ## v1 scope
//!
//! Real GD&T verification requires hundreds of surface sample points
//! and characteristic-specific evaluators (flatness wants a min-zone
//! plane fit, cylindricity wants a min-zone cylinder fit, position
//! wants the diameter of the tolerance cylinder around the datum
//! axis, etc.). v1 doesn't ship those evaluators — instead it
//! supports the most common Form/Profile/Runout characteristics which
//! reduce to a *single scalar* the user has already produced via
//! [`crate::Measurement`]:
//!
//! - Straightness / Flatness / Circularity / Cylindricity /
//!   ProfileLine / ProfileSurface / CircularRunout / TotalRunout —
//!   the actual value is interpreted as the **diameter / width of the
//!   tolerance zone the measured form fits inside**. Pass iff
//!   `actual <= tolerance_value`.
//!
//! All other characteristics (Perpendicularity, Angularity,
//! Parallelism, Position, Concentricity, Symmetry) return
//! [`crate::InspectError::GdtNotImplemented`]; they need a datum-frame
//! aware evaluator which is on the Phase 25.5 roadmap.

use serde::{Deserialize, Serialize};
use valenx_techdraw::gdt::{GdtSymbol, GeometricCharacteristic};

use crate::error::InspectError;
use crate::report::CheckResult;

/// Verification rule expressing how to interpret the GD&T
/// `tolerance_value` against an actual scalar.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum GdtRule {
    /// The actual form-zone width must be `<=` the spec value.
    ZoneWidthMax,
    /// Datum-frame evaluator required — not in v1.
    RequiresDatumFrame,
}

impl GdtRule {
    /// Pick the v1 rule for a characteristic.
    pub fn for_characteristic(c: GeometricCharacteristic) -> Self {
        use GeometricCharacteristic::*;
        match c {
            Straightness | Flatness | Circularity | Cylindricity | ProfileLine
            | ProfileSurface | CircularRunout | TotalRunout => Self::ZoneWidthMax,
            Perpendicularity | Angularity | Parallelism | Position | Concentricity | Symmetry => {
                Self::RequiresDatumFrame
            }
        }
    }
}

/// One actual-vs-frame verification job.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GdtCheck {
    /// The ASME Y14.5 feature control frame (from Phase 18 TechDraw).
    pub frame: GdtSymbol,
    /// Measured scalar — meaning depends on
    /// [`GdtRule::for_characteristic`].
    pub actual_value: f64,
}

impl GdtCheck {
    /// Build a check with a frame and a measured value.
    pub fn new(frame: GdtSymbol, actual_value: f64) -> Self {
        Self {
            frame,
            actual_value,
        }
    }

    /// Evaluate the check.
    ///
    /// Returns:
    /// - `Ok(CheckResult::Pass)` when the rule applies and actual is
    ///   within the spec.
    /// - `Ok(CheckResult::Fail)` when actual is out of spec.
    /// - `Err(InspectError::BadParameter)` when the tolerance value
    ///   string fails to parse as `f64` (after stripping a leading
    ///   diameter sign `⌀`).
    /// - `Err(InspectError::GdtNotImplemented)` when the characteristic
    ///   needs a datum-frame evaluator (Phase 25.5).
    pub fn verify(&self) -> Result<CheckResult, InspectError> {
        let rule = GdtRule::for_characteristic(self.frame.geometric_characteristic);
        match rule {
            GdtRule::ZoneWidthMax => {
                let spec = parse_tolerance_value(&self.frame.tolerance_value)?;
                if spec <= 0.0 {
                    return Err(InspectError::BadParameter {
                        name: "tolerance_value",
                        reason: format!("must be > 0 (got {spec})"),
                    });
                }
                if self.actual_value <= spec {
                    Ok(CheckResult::Pass)
                } else {
                    Ok(CheckResult::Fail)
                }
            }
            GdtRule::RequiresDatumFrame => Err(InspectError::GdtNotImplemented {
                characteristic: characteristic_label(self.frame.geometric_characteristic),
                reason: "needs datum-frame aware evaluator (Phase 25.5)".into(),
            }),
        }
    }
}

/// Parse a tolerance-value string like `"0.1"` or `"⌀0.05"` to a
/// positive `f64`. Errors are surfaced as `BadParameter`.
pub fn parse_tolerance_value(s: &str) -> Result<f64, InspectError> {
    let trimmed = s.trim().trim_start_matches('⌀').trim();
    trimmed
        .parse::<f64>()
        .map_err(|e| InspectError::BadParameter {
            name: "tolerance_value",
            reason: format!("could not parse {s:?}: {e}"),
        })
}

/// Stable static-string label for a characteristic — used by error
/// messages so the static-lifetime `characteristic` field on
/// [`InspectError::GdtNotImplemented`] is satisfied.
fn characteristic_label(c: GeometricCharacteristic) -> &'static str {
    use GeometricCharacteristic::*;
    match c {
        Straightness => "Straightness",
        Flatness => "Flatness",
        Circularity => "Circularity",
        Cylindricity => "Cylindricity",
        ProfileLine => "ProfileLine",
        ProfileSurface => "ProfileSurface",
        Perpendicularity => "Perpendicularity",
        Angularity => "Angularity",
        Parallelism => "Parallelism",
        Position => "Position",
        Concentricity => "Concentricity",
        Symmetry => "Symmetry",
        CircularRunout => "CircularRunout",
        TotalRunout => "TotalRunout",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_techdraw::gdt::GdtSymbol;

    #[test]
    fn flatness_within_spec_passes() {
        let frame = GdtSymbol::new(
            [0.0, 0.0],
            GeometricCharacteristic::Flatness,
            "0.1",
        );
        let c = GdtCheck::new(frame, 0.05);
        assert_eq!(c.verify().unwrap(), CheckResult::Pass);
    }

    #[test]
    fn flatness_out_of_spec_fails() {
        let frame = GdtSymbol::new(
            [0.0, 0.0],
            GeometricCharacteristic::Flatness,
            "0.05",
        );
        let c = GdtCheck::new(frame, 0.07);
        assert_eq!(c.verify().unwrap(), CheckResult::Fail);
    }

    #[test]
    fn position_returns_not_implemented() {
        let frame = GdtSymbol::new(
            [0.0, 0.0],
            GeometricCharacteristic::Position,
            "⌀0.1",
        );
        let c = GdtCheck::new(frame, 0.05);
        assert!(matches!(
            c.verify().unwrap_err(),
            InspectError::GdtNotImplemented { .. }
        ));
    }

    #[test]
    fn diameter_sign_strips() {
        assert_eq!(parse_tolerance_value("⌀0.1").unwrap(), 0.1);
        assert_eq!(parse_tolerance_value(" 0.05 ").unwrap(), 0.05);
    }

    #[test]
    fn bad_tolerance_string_errors() {
        assert!(parse_tolerance_value("not a number").is_err());
    }

    #[test]
    fn rule_dispatch_is_complete() {
        // Every characteristic gets *some* rule.
        use GeometricCharacteristic::*;
        let all = [
            Straightness,
            Flatness,
            Circularity,
            Cylindricity,
            ProfileLine,
            ProfileSurface,
            Perpendicularity,
            Angularity,
            Parallelism,
            Position,
            Concentricity,
            Symmetry,
            CircularRunout,
            TotalRunout,
        ];
        for c in all {
            let _ = GdtRule::for_characteristic(c);
        }
    }
}
