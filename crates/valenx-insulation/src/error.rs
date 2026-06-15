//! Error taxonomy for the building-insulation heat-loss models.
//!
//! Every fallible constructor in this crate funnels its rejection
//! through [`InsulationError`]. The variants distinguish a bad scalar
//! input (a non-positive thickness, conductivity, coefficient, or
//! area) from a structurally empty wall assembly, so a caller can
//! branch on the failure mode without string-matching the message.

use thiserror::Error;

/// Error raised while validating heat-loss model inputs.
///
/// The physics here is one-dimensional steady-state conduction, where
/// every quantity (thickness, conductivity, surface coefficient,
/// area, temperature difference) is a strictly positive real. A zero
/// or negative value is non-physical and is rejected at construction
/// time rather than silently producing a `NaN`, an infinity, or a
/// negative heat-loss rate downstream.
#[derive(Debug, Error)]
pub enum InsulationError {
    /// A scalar input was not strictly positive (or not finite).
    ///
    /// Carries the offending parameter name and the value that was
    /// supplied so the message points straight at the bad input.
    #[error("non-positive value for `{name}`: expected a finite value > 0, got {value}")]
    NonPositive {
        /// Name of the rejected parameter (e.g. `"thickness_m"`).
        name: &'static str,
        /// The value that failed the strictly-positive-and-finite test.
        value: f64,
    },

    /// A composite wall was assembled with no thermal layers and no
    /// surface films, so its total resistance would be zero and the
    /// derived U-value infinite.
    #[error("empty wall assembly: a wall must contain at least one layer or surface film")]
    EmptyAssembly,
}

/// Coarse classification of an [`InsulationError`], mirroring the
/// `ErrorCategory` convention used across the Valenx workspace.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// The caller supplied an invalid value or an empty assembly.
    Input,
    /// The requested computation is outside the model's valid domain.
    Algorithm,
}

impl InsulationError {
    /// Construct a [`InsulationError::NonPositive`] if `value` is not a
    /// finite, strictly positive real; otherwise return `Ok(value)`.
    ///
    /// This is the single validated gate used by every public
    /// constructor in the crate, which keeps the strictly-positive rule
    /// in exactly one place.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_insulation::error::InsulationError;
    ///
    /// assert!(InsulationError::require_positive("thickness_m", 0.1).is_ok());
    /// assert!(InsulationError::require_positive("thickness_m", 0.0).is_err());
    /// assert!(InsulationError::require_positive("thickness_m", -1.0).is_err());
    /// assert!(InsulationError::require_positive("thickness_m", f64::NAN).is_err());
    /// ```
    pub fn require_positive(name: &'static str, value: f64) -> Result<f64, InsulationError> {
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(InsulationError::NonPositive { name, value })
        }
    }

    /// Stable, kebab-cased identifier for this error, suitable for
    /// logs and machine matching.
    pub fn code(&self) -> &'static str {
        match self {
            InsulationError::NonPositive { .. } => "insulation.non-positive",
            InsulationError::EmptyAssembly => "insulation.empty-assembly",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            InsulationError::NonPositive { .. } => ErrorCategory::Input,
            InsulationError::EmptyAssembly => ErrorCategory::Input,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn require_positive_accepts_positive_finite() {
        let v = InsulationError::require_positive("x", 2.5).unwrap();
        assert!((v - 2.5).abs() < 1e-12);
    }

    #[test]
    fn require_positive_rejects_zero_negative_and_nonfinite() {
        assert!(InsulationError::require_positive("x", 0.0).is_err());
        assert!(InsulationError::require_positive("x", -0.001).is_err());
        assert!(InsulationError::require_positive("x", f64::NAN).is_err());
        assert!(InsulationError::require_positive("x", f64::INFINITY).is_err());
    }

    #[test]
    fn codes_and_categories_are_stable() {
        let np = InsulationError::NonPositive {
            name: "x",
            value: -1.0,
        };
        assert_eq!(np.code(), "insulation.non-positive");
        assert_eq!(np.category(), ErrorCategory::Input);

        let empty = InsulationError::EmptyAssembly;
        assert_eq!(empty.code(), "insulation.empty-assembly");
        assert_eq!(empty.category(), ErrorCategory::Input);
    }

    #[test]
    fn display_mentions_the_parameter_name() {
        let e = InsulationError::NonPositive {
            name: "conductivity_w_per_m_k",
            value: 0.0,
        };
        let msg = format!("{e}");
        assert!(msg.contains("conductivity_w_per_m_k"));
    }
}
