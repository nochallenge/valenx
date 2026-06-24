//! The crate-wide error type.

use thiserror::Error;

/// Errors produced by the simulated-sensor toolkit.
///
/// Every sensor constructor validates its parameters and returns
/// [`Result<_, SensorError>`] rather than panicking or silently producing a
/// `NaN`, so a bad configuration (a non-positive focal length, a zero LiDAR beam
/// count, a negative noise standard deviation, …) is a recoverable error caught
/// at build time — *fail loud, fail early*.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum SensorError {
    /// A configuration parameter was out of its valid range — e.g. a
    /// non-positive focal length, a zero or negative sensor size, a beam count
    /// of zero, or `min_range >= max_range`.
    #[error("invalid configuration: {0}")]
    InvalidConfig(String),

    /// A noise parameter (a standard deviation or bias magnitude) was negative,
    /// or otherwise not usable as a `N(·, std)` spread.
    #[error("invalid noise parameter: {0}")]
    InvalidNoise(String),

    /// A supplied geometry was degenerate — e.g. a plane with a zero normal, a
    /// sphere with a non-positive radius, or a triangle whose vertices are
    /// collinear (zero area).
    #[error("degenerate geometry: {0}")]
    DegenerateGeometry(String),

    /// A non-finite value (`NaN` / `±∞`) reached an API that requires finite
    /// inputs — e.g. a camera point or a vehicle state component.
    #[error("non-finite value: {0}")]
    NonFinite(String),
}
