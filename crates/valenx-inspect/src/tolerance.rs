//! Dimensional tolerance bands.
//!
//! A [`Tolerance`] expresses a "nominal ± deviations" band. The two
//! deviations are signed (lower may be negative). `evaluate(value)`
//! returns one of [`crate::CheckResult`] variants.

use serde::{Deserialize, Serialize};

use crate::error::InspectError;
use crate::report::CheckResult;

/// Nominal-with-deviations tolerance band.
///
/// Semantics:
/// - `nominal` — the target value.
/// - `upper_dev` — signed upper deviation (typically ≥ 0).
/// - `lower_dev` — signed lower deviation (typically ≤ 0).
/// - `warn_fraction` — the inner fraction of the band that's still a
///   Pass; values outside this inner band but still inside the full
///   band yield [`CheckResult::Warning`]. Set to `1.0` to disable.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Tolerance {
    /// Target value.
    pub nominal: f64,
    /// Upper deviation (signed, usually positive).
    pub upper_dev: f64,
    /// Lower deviation (signed, usually negative).
    pub lower_dev: f64,
    /// Inner-band Pass fraction; 0.0..=1.0. Outside the inner band
    /// but inside the full band → `Warning`.
    pub warn_fraction: f64,
}

impl Tolerance {
    /// Symmetric band around `nominal` (no warning band).
    pub fn symmetric(nominal: f64, dev: f64) -> Self {
        Self {
            nominal,
            upper_dev: dev.abs(),
            lower_dev: -dev.abs(),
            warn_fraction: 1.0,
        }
    }

    /// Asymmetric band. Returns `Err(BadParameter)` when
    /// `upper_dev < lower_dev`.
    pub fn asymmetric(nominal: f64, upper_dev: f64, lower_dev: f64) -> Result<Self, InspectError> {
        if upper_dev < lower_dev {
            return Err(InspectError::BadParameter {
                name: "upper_dev",
                reason: "upper_dev must be >= lower_dev".into(),
            });
        }
        Ok(Self {
            nominal,
            upper_dev,
            lower_dev,
            warn_fraction: 1.0,
        })
    }

    /// Replace the warn-fraction (clamped to `0.0..=1.0`).
    pub fn with_warn_fraction(mut self, f: f64) -> Self {
        self.warn_fraction = f.clamp(0.0, 1.0);
        self
    }

    /// Cheap "is in full band?" predicate.
    pub fn accepts(&self, value: f64) -> bool {
        let dev = value - self.nominal;
        dev >= self.lower_dev && dev <= self.upper_dev
    }

    /// Full evaluation including the inner Pass band.
    pub fn evaluate(&self, value: f64) -> CheckResult {
        let dev = value - self.nominal;
        if dev < self.lower_dev || dev > self.upper_dev {
            return CheckResult::Fail;
        }
        let inner_upper = self.upper_dev * self.warn_fraction;
        let inner_lower = self.lower_dev * self.warn_fraction;
        if dev >= inner_lower && dev <= inner_upper {
            CheckResult::Pass
        } else {
            CheckResult::Warning
        }
    }

    /// Full width of the band — `upper_dev - lower_dev`.
    pub fn width(&self) -> f64 {
        self.upper_dev - self.lower_dev
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symmetric_pass_fail() {
        let t = Tolerance::symmetric(10.0, 0.1);
        assert_eq!(t.evaluate(10.0), CheckResult::Pass);
        assert_eq!(t.evaluate(10.05), CheckResult::Pass);
        assert_eq!(t.evaluate(10.11), CheckResult::Fail);
        assert_eq!(t.evaluate(9.89), CheckResult::Fail);
    }

    #[test]
    fn asymmetric_band() {
        let t = Tolerance::asymmetric(10.0, 0.2, -0.05).unwrap();
        assert!(t.accepts(10.15));
        assert!(!t.accepts(9.9));
    }

    #[test]
    fn warn_band_triggers_warning() {
        let t = Tolerance::symmetric(10.0, 0.1).with_warn_fraction(0.5);
        // Inner band is ±0.05; outer is ±0.10.
        assert_eq!(t.evaluate(10.0), CheckResult::Pass);
        assert_eq!(t.evaluate(10.04), CheckResult::Pass);
        assert_eq!(t.evaluate(10.07), CheckResult::Warning);
        assert_eq!(t.evaluate(10.15), CheckResult::Fail);
    }

    #[test]
    fn bad_band_errors() {
        assert!(Tolerance::asymmetric(10.0, -0.1, 0.1).is_err());
    }

    #[test]
    fn width_is_sum_of_devs() {
        let t = Tolerance::asymmetric(10.0, 0.2, -0.1).unwrap();
        assert!((t.width() - 0.3).abs() < 1e-12);
    }
}
