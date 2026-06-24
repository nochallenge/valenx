//! The crate-wide error type.

use thiserror::Error;

/// Errors produced by the particle-fluid toolkit.
///
/// Every constructor and configuration validates its parameters and returns
/// [`Result<_, FluidError>`] rather than panicking or silently producing a
/// `NaN`. A bad configuration — a non-positive smoothing length `h`, a
/// non-positive particle mass or rest density, a non-finite particle position —
/// is therefore a recoverable error caught up front: *fail loud, fail early*.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum FluidError {
    /// A configuration parameter was out of its valid range — e.g. a
    /// non-positive smoothing length `h`, a non-positive particle mass or rest
    /// density, a non-positive sound speed, or a non-positive grid cell size.
    #[error("invalid configuration: {0}")]
    InvalidConfig(String),

    /// A supplied domain / boundary was degenerate — e.g. a box whose `min`
    /// is not strictly less than its `max` on some axis.
    #[error("invalid domain: {0}")]
    InvalidDomain(String),

    /// A non-finite value (`NaN` / `±∞`) reached an API that requires finite
    /// inputs — e.g. a particle position, velocity, or gravity vector.
    #[error("non-finite value: {0}")]
    NonFinite(String),
}
