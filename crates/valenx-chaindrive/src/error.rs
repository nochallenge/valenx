//! Roller-chain drive error taxonomy.
//!
//! A single [`ChainDriveError`] enum covers every way the inputs to a
//! chain-drive calculation can be invalid. The validated constructors on
//! the spec / calculation types build values of this type rather than
//! panicking, so callers always get a typed, categorised failure.

use thiserror::Error;

/// Errors raised by roller-chain drive calculations.
#[derive(Debug, Error)]
pub enum ChainDriveError {
    /// A continuous parameter (pitch, speed, torque, centre distance)
    /// was non-finite or outside its allowed range.
    #[error("bad parameter `{name}`: {reason}")]
    BadParameter {
        /// Parameter name, e.g. `"pitch_mm"`.
        name: &'static str,
        /// Human-readable reason the value was rejected.
        reason: String,
    },

    /// A sprocket tooth count was below the practical minimum.
    ///
    /// Real roller-chain sprockets need enough teeth for the chain to
    /// wrap without excessive polygonal (chordal) action; fewer than a
    /// handful of teeth is not a usable drive, so the constructors reject
    /// it rather than returning a meaningless ratio.
    #[error("tooth count `{name}` too low: got {got}, need >= {min}")]
    TeethTooFew {
        /// Which sprocket, e.g. `"driver_teeth"`.
        name: &'static str,
        /// The supplied (too-small) count.
        got: u32,
        /// The enforced minimum.
        min: u32,
    },

    /// The two-sprocket geometry is degenerate: the centre distance is
    /// too small for the sprockets to coexist, so no valid chain wrap
    /// exists.
    #[error("degenerate geometry: {0}")]
    Degenerate(String),
}

impl ChainDriveError {
    /// Build a [`ChainDriveError::BadParameter`] from any displayable
    /// reason. Keeps the call sites in the calculation code terse.
    pub fn bad_parameter(name: &'static str, reason: impl Into<String>) -> Self {
        Self::BadParameter {
            name,
            reason: reason.into(),
        }
    }

    /// Build a [`ChainDriveError::Degenerate`] from any displayable
    /// message.
    pub fn degenerate(message: impl Into<String>) -> Self {
        Self::Degenerate(message.into())
    }

    /// Stable kebab-cased identifier, suitable for logs / telemetry that
    /// must not depend on the human-readable `Display` text.
    pub fn code(&self) -> &'static str {
        match self {
            ChainDriveError::BadParameter { .. } => "chaindrive.bad_parameter",
            ChainDriveError::TeethTooFew { .. } => "chaindrive.teeth_too_few",
            ChainDriveError::Degenerate(_) => "chaindrive.degenerate",
        }
    }

    /// Coarse category for grouping failures in a UI or report.
    pub fn category(&self) -> ErrorCategory {
        match self {
            ChainDriveError::BadParameter { .. } => ErrorCategory::Config,
            ChainDriveError::TeethTooFew { .. } => ErrorCategory::Input,
            ChainDriveError::Degenerate(_) => ErrorCategory::Algorithm,
        }
    }
}

/// Coarse error category.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Discrete user input (a tooth count).
    Input,
    /// A tunable continuous knob (pitch, speed, torque, centre distance).
    Config,
    /// An algorithm-domain failure (degenerate geometry).
    Algorithm,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_are_stable_and_distinct() {
        let bad = ChainDriveError::bad_parameter("pitch_mm", "must be > 0");
        let few = ChainDriveError::TeethTooFew {
            name: "driver_teeth",
            got: 2,
            min: 7,
        };
        let deg = ChainDriveError::degenerate("centre distance too small");

        assert_eq!(bad.code(), "chaindrive.bad_parameter");
        assert_eq!(few.code(), "chaindrive.teeth_too_few");
        assert_eq!(deg.code(), "chaindrive.degenerate");

        assert_eq!(bad.category(), ErrorCategory::Config);
        assert_eq!(few.category(), ErrorCategory::Input);
        assert_eq!(deg.category(), ErrorCategory::Algorithm);
    }

    #[test]
    fn display_includes_context() {
        let few = ChainDriveError::TeethTooFew {
            name: "driver_teeth",
            got: 3,
            min: 7,
        };
        let text = few.to_string();
        assert!(text.contains("driver_teeth"), "got: {text}");
        assert!(text.contains('3'), "got: {text}");
    }
}
