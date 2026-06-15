//! Error type for `valenx-diffusion`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, DiffusionError>`]. The variants are deliberately coarse:
//! a 1-D diffusion caller really only ever hits three failure classes —
//!
//! 1. a physical / geometric parameter is out of range
//!    ([`DiffusionError::BadParameter`]) — a non-positive diffusion
//!    coefficient, a zero-width cell, a domain with fewer than two
//!    nodes, a negative evaluation time;
//! 2. the requested explicit time step would make the forward-time
//!    centred-space scheme unstable
//!    ([`DiffusionError::Unstable`]) — `dt > dx^2 / (2 D)`;
//! 3. two arrays that must share a length do not
//!    ([`DiffusionError::DimensionMismatch`]) — an initial-condition
//!    vector whose length differs from the grid node count.
//!
//! Use [`DiffusionError::code`] for stable log / telemetry tagging and
//! [`DiffusionError::category`] to bucket failures into `input` /
//! `numerics` without matching every variant. The pattern mirrors
//! `valenx-cfd-native`'s `CfdError` and the other native numerical
//! crates.

use thiserror::Error;

/// Errors produced by `valenx-diffusion`.
#[derive(Debug, Clone, PartialEq, Error)]
pub enum DiffusionError {
    /// A physical or geometric parameter was out of its admissible
    /// range — a non-positive diffusion coefficient, a non-positive
    /// cell spacing, a grid with fewer than two nodes, a negative
    /// time, a non-finite value.
    #[error("invalid diffusion parameter `{name}`: {reason}")]
    BadParameter {
        /// Name of the offending parameter (e.g. `"D"`, `"dx"`, `"t"`).
        name: &'static str,
        /// Human-readable explanation.
        reason: String,
    },

    /// The requested explicit time step violates the forward-time
    /// centred-space stability limit `dt <= dx^2 / (2 D)`. The fields
    /// report the offered step and the largest stable step so the
    /// caller can clamp or sub-cycle.
    #[error(
        "explicit step unstable: dt = {dt} exceeds the FTCS limit dx^2/(2 D) = {dt_max} \
         (dx = {dx}, D = {d})"
    )]
    Unstable {
        /// The time step the caller asked for.
        dt: f64,
        /// The largest stable time step, `dx^2 / (2 D)`.
        dt_max: f64,
        /// The cell spacing in use.
        dx: f64,
        /// The diffusion coefficient in use.
        d: f64,
    },

    /// Two quantities that must agree on a length disagree — an
    /// initial-condition vector whose length differs from the grid's
    /// node count.
    #[error("dimension mismatch: {reason}")]
    DimensionMismatch {
        /// Human-readable reason, naming both sides and their sizes.
        reason: String,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on the error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// User-supplied input is wrong (bad parameter, size mismatch).
    Input,
    /// A numerical-stability constraint was violated.
    Numerics,
}

impl DiffusionError {
    /// A stable, dot-namespaced identifier suitable for log / telemetry
    /// tagging. Format: `"diffusion.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            DiffusionError::BadParameter { .. } => "diffusion.bad_parameter",
            DiffusionError::Unstable { .. } => "diffusion.unstable",
            DiffusionError::DimensionMismatch { .. } => "diffusion.dimension_mismatch",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            DiffusionError::BadParameter { .. } | DiffusionError::DimensionMismatch { .. } => {
                "input"
            }
            DiffusionError::Unstable { .. } => "numerics",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            DiffusionError::BadParameter { .. } | DiffusionError::DimensionMismatch { .. } => {
                ErrorCategory::Input
            }
            DiffusionError::Unstable { .. } => ErrorCategory::Numerics,
        }
    }

    /// Convenience constructor for [`DiffusionError::BadParameter`].
    pub fn bad_parameter(name: &'static str, reason: impl Into<String>) -> Self {
        DiffusionError::BadParameter {
            name,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`DiffusionError::DimensionMismatch`].
    pub fn dimension(reason: impl Into<String>) -> Self {
        DiffusionError::DimensionMismatch {
            reason: reason.into(),
        }
    }
}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, DiffusionError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_are_stable() {
        assert_eq!(
            DiffusionError::bad_parameter("D", "must be positive").code(),
            "diffusion.bad_parameter"
        );
        assert_eq!(
            DiffusionError::Unstable {
                dt: 1.0,
                dt_max: 0.5,
                dx: 1.0,
                d: 1.0,
            }
            .code(),
            "diffusion.unstable"
        );
        assert_eq!(
            DiffusionError::dimension("3 vs 4").code(),
            "diffusion.dimension_mismatch"
        );
    }

    #[test]
    fn categories_match_variants() {
        assert_eq!(
            DiffusionError::bad_parameter("dx", "zero").category(),
            "input"
        );
        assert_eq!(
            DiffusionError::bad_parameter("dx", "zero").category_enum(),
            ErrorCategory::Input
        );
        assert_eq!(
            DiffusionError::dimension("x").category_enum(),
            ErrorCategory::Input
        );
        let unstable = DiffusionError::Unstable {
            dt: 2.0,
            dt_max: 0.5,
            dx: 1.0,
            d: 1.0,
        };
        assert_eq!(unstable.category(), "numerics");
        assert_eq!(unstable.category_enum(), ErrorCategory::Numerics);
    }

    #[test]
    fn display_carries_context() {
        let e = DiffusionError::bad_parameter("D", "must be positive");
        let msg = e.to_string();
        assert!(msg.contains('D'), "got: {msg}");
        assert!(msg.contains("must be positive"), "got: {msg}");

        let e = DiffusionError::Unstable {
            dt: 2.0,
            dt_max: 0.5,
            dx: 1.0,
            d: 1.0,
        };
        let msg = e.to_string();
        assert!(msg.contains('2'), "got: {msg}");
        assert!(msg.contains("0.5"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(DiffusionError::dimension("x"));
        assert!(err.to_string().contains('x'));
    }
}
