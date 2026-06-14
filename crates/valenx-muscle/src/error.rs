//! Error type for invalid muscle parameters.

use thiserror::Error;

/// Something was wrong with the parameters describing a [`crate::Muscle`].
#[derive(Debug, Clone, PartialEq, Error)]
pub enum MuscleError {
    /// A quantity that must be strictly positive (and finite) was not.
    #[error("{field} must be a positive finite number, got {value}")]
    NotPositive {
        /// Which parameter was invalid.
        field: &'static str,
        /// The offending value.
        value: f64,
    },
    /// The pennation angle must lie in `[0, 90)` degrees.
    #[error("pennation angle must be in [0, 90) degrees, got {0}")]
    BadPennation(f64),
}
