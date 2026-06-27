//! The crate-wide error type.

use thiserror::Error;

/// Errors produced by the ocean wave-field and buoyancy toolkit.
///
/// Every constructor and configuration validates its parameters and returns
/// [`Result<_, OceanError>`] rather than panicking or silently producing a
/// `NaN`. A bad configuration — a non-positive wavelength or amplitude, a
/// non-positive water density or gravity, a non-finite position — is therefore a
/// recoverable error caught up front: *fail loud, fail early*.
#[derive(Debug, Error, Clone, PartialEq)]
#[non_exhaustive]
pub enum OceanError {
    /// A configuration parameter was out of its valid range — e.g. a
    /// non-positive wavelength, amplitude, steepness, water density, gravity, or
    /// body mass, or a non-positive waterplane area or hull volume.
    #[error("invalid configuration: {0}")]
    InvalidConfig(String),

    /// A non-finite value (`NaN` / `±∞`) reached an API that requires finite
    /// inputs — e.g. a wave direction, a sample point, or a body-state
    /// component.
    #[error("non-finite value: {0}")]
    NonFinite(String),
}
