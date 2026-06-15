//! Gearbox error taxonomy.
//!
//! A single [`GearboxError`] enum covers every failure mode in the
//! crate. Validated constructors elsewhere in the crate return these
//! variants instead of panicking, so callers can recover from bad
//! tooth counts, out-of-range efficiencies, or degenerate planetary
//! geometry.

use thiserror::Error;

/// Errors raised by gear-train construction and analysis.
#[derive(Debug, Error)]
pub enum GearboxError {
    /// A tooth count was zero. Gear ratios divide by tooth counts, so a
    /// count of zero is undefined.
    #[error("tooth count `{name}` must be >= 1, got 0")]
    ZeroTeeth {
        /// Which tooth count was zero (e.g. `"teeth_in"`).
        name: &'static str,
    },

    /// An efficiency was outside the half-open interval `(0, 1]`.
    ///
    /// Efficiency is the fraction of input power delivered to the
    /// output. It must be positive (a stage that delivers nothing is
    /// not a stage) and at most `1` (a passive gear train cannot create
    /// power). The ideal lossless case is exactly `1.0`.
    #[error("efficiency must be in (0, 1], got {value}")]
    EfficiencyOutOfRange {
        /// The offending value.
        value: f64,
    },

    /// A numeric input was not finite (`NaN` or `±∞`).
    #[error("parameter `{name}` must be finite, got {value}")]
    NotFinite {
        /// Which parameter was non-finite.
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A compound train was built with no stages. The overall ratio of
    /// an empty train is undefined.
    #[error("compound train must contain at least one stage, got 0")]
    EmptyTrain,

    /// A planetary train had a non-positive sun or ring tooth count, or
    /// a ring not larger than the sun. A physical fixed-ring epicyclic
    /// set requires `ring > sun >= 1` so the planets fit between them.
    #[error("degenerate planetary geometry: {reason}")]
    DegeneratePlanetary {
        /// Human-readable explanation.
        reason: String,
    },
}

/// Coarse error category, useful for UI grouping and metrics.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Bad user-supplied input (tooth counts, train shape).
    Input,
    /// Out-of-range tunable knob (efficiency).
    Config,
    /// Algorithm / geometry domain violation.
    Algorithm,
}

impl GearboxError {
    /// Stable, dotted identifier for logs and telemetry.
    ///
    /// The returned string is part of the crate's public contract and
    /// will not change for a given variant across patch releases.
    pub fn code(&self) -> &'static str {
        match self {
            GearboxError::ZeroTeeth { .. } => "gearbox.zero_teeth",
            GearboxError::EfficiencyOutOfRange { .. } => "gearbox.efficiency_out_of_range",
            GearboxError::NotFinite { .. } => "gearbox.not_finite",
            GearboxError::EmptyTrain => "gearbox.empty_train",
            GearboxError::DegeneratePlanetary { .. } => "gearbox.degenerate_planetary",
        }
    }

    /// Coarse category for the variant.
    pub fn category(&self) -> ErrorCategory {
        match self {
            GearboxError::ZeroTeeth { .. }
            | GearboxError::EmptyTrain
            | GearboxError::NotFinite { .. } => ErrorCategory::Input,
            GearboxError::EfficiencyOutOfRange { .. } => ErrorCategory::Config,
            GearboxError::DegeneratePlanetary { .. } => ErrorCategory::Algorithm,
        }
    }
}

/// Validate that `value` lies in the efficiency interval `(0, 1]`.
///
/// Returns the value unchanged on success, or
/// [`GearboxError::EfficiencyOutOfRange`] / [`GearboxError::NotFinite`]
/// otherwise. Shared by every stage constructor so the rule is defined
/// in exactly one place.
pub(crate) fn check_efficiency(value: f64) -> Result<f64, GearboxError> {
    if !value.is_finite() {
        return Err(GearboxError::NotFinite {
            name: "efficiency",
            value,
        });
    }
    if value <= 0.0 || value > 1.0 {
        return Err(GearboxError::EfficiencyOutOfRange { value });
    }
    Ok(value)
}
