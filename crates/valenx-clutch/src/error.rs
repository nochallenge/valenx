//! Error taxonomy for the friction-clutch models.
//!
//! Every fallible constructor in this crate returns a [`ClutchError`].
//! The variants are intentionally coarse: a clutch capacity calculation
//! has only a handful of physically-meaningful ways to be ill-posed
//! (non-positive geometry, an inverted radius pair, a negative
//! coefficient of friction, a fractional surface count, and so on), so
//! the taxonomy enumerates exactly those.

use thiserror::Error;

/// Errors raised while validating clutch inputs or evaluating the
/// torque / power models.
///
/// Returned by the validated constructors ([`crate::clutch::ClutchGeometry::new`],
/// [`crate::clutch::FrictionClutch::new`]) and by the operating-condition
/// helpers that reject non-physical inputs.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ClutchError {
    /// A scalar parameter was outside its physically-admissible range.
    ///
    /// `name` is the offending parameter (a `'static` identifier so it
    /// can be matched on), `value` is what was supplied, and `reason`
    /// describes the constraint that was violated.
    #[error("invalid parameter `{name}` = {value}: {reason}")]
    InvalidParameter {
        /// Identifier of the offending parameter.
        name: &'static str,
        /// The value that was supplied.
        value: f64,
        /// Human-readable description of the violated constraint.
        reason: &'static str,
    },

    /// The inner radius was not strictly smaller than the outer radius.
    ///
    /// The annular friction face is the region `ri <= r <= ro`; if
    /// `ri >= ro` that region is empty (or inverted) and the model has
    /// no meaning. Stored in millimetres, mirroring the constructor
    /// inputs.
    #[error(
        "inner radius ({inner_mm} mm) must be strictly less than outer radius ({outer_mm} mm)"
    )]
    InvertedRadii {
        /// Inner radius of the friction annulus, in millimetres.
        inner_mm: f64,
        /// Outer radius of the friction annulus, in millimetres.
        outer_mm: f64,
    },

    /// The number of friction surfaces in contact was not a positive
    /// integer.
    ///
    /// A single-plate clutch has two faces in contact (`N = 2`); a
    /// multi-plate clutch has more. `N` must be at least one and a whole
    /// number.
    #[error("friction-surface count must be a positive whole number, got {0}")]
    InvalidSurfaceCount(f64),
}

impl ClutchError {
    /// A short, stable, kebab-cased identifier for the error variant.
    ///
    /// Useful for logging, metrics labels, or mapping to UI messages
    /// without matching on the (translatable) `Display` text.
    ///
    /// ```
    /// use valenx_clutch::error::ClutchError;
    /// let err = ClutchError::InvalidSurfaceCount(0.0);
    /// assert_eq!(err.code(), "clutch.invalid-surface-count");
    /// ```
    pub fn code(&self) -> &'static str {
        match self {
            ClutchError::InvalidParameter { .. } => "clutch.invalid-parameter",
            ClutchError::InvertedRadii { .. } => "clutch.inverted-radii",
            ClutchError::InvalidSurfaceCount(_) => "clutch.invalid-surface-count",
        }
    }
}
