//! Error type for the riveted-joint calculator.

use thiserror::Error;

/// Shorthand for `Result<T, RivetError>`.
pub type Result<T> = core::result::Result<T, RivetError>;

/// Anything that can go wrong building a rivet, plate, or joint, or
/// evaluating its strength.
///
/// This enum is `#[non_exhaustive]`: new variants may be added in a
/// future release without it being a breaking change, so downstream
/// `match` arms must include a wildcard.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum RivetError {
    /// A geometric or material quantity that must be strictly positive
    /// and finite was zero, negative, or non-finite (`NaN`/`±∞`). These
    /// feed directly into areas, stresses, and ratios, so a bad value
    /// would otherwise produce a silent `NaN`/`Inf` answer.
    ///
    /// Carries the parameter name and the offending value.
    #[error("parameter `{name}` must be finite and positive, got {value}")]
    NotPositive {
        /// Which quantity was bad (e.g. `"diameter"`, `"thickness"`).
        name: &'static str,
        /// The offending value.
        value: f64,
    },

    /// A count that must be at least one (rivets per row / number of
    /// rows / shear planes) was zero.
    #[error("count `{name}` must be at least 1, got 0")]
    ZeroCount {
        /// Which count was zero (e.g. `"rivets_per_row"`).
        name: &'static str,
    },

    /// The net plate section is non-positive: the rivet holes in a row
    /// remove at least the full plate width, leaving no material to
    /// carry tension. Equivalent to `width <= rivets_per_row * diameter`.
    #[error(
        "net width is non-positive: width {width} m with {holes} hole(s) of \
         diameter {diameter} m removes {removed} m, leaving no net section"
    )]
    NetSectionNonPositive {
        /// Gross plate width, metres.
        width: f64,
        /// Number of holes across the critical row.
        holes: u32,
        /// Hole (rivet) diameter, metres.
        diameter: f64,
        /// Total material removed by the holes, metres.
        removed: f64,
    },
}

impl RivetError {
    /// Validate that `value` is finite and strictly positive, returning
    /// it on success or a [`RivetError::NotPositive`] carrying `name`.
    ///
    /// Used by every validated constructor in the crate so the rule is
    /// stated once.
    pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<f64> {
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(RivetError::NotPositive { name, value })
        }
    }

    /// Validate that `value` is at least one, returning it on success or
    /// a [`RivetError::ZeroCount`] carrying `name`.
    pub(crate) fn require_count(name: &'static str, value: u32) -> Result<u32> {
        if value >= 1 {
            Ok(value)
        } else {
            Err(RivetError::ZeroCount { name })
        }
    }
}
