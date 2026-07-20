//! Error taxonomy for the epicyclic gear-train kinematics.
//!
//! Every constructor in this crate funnels its input validation through
//! the helpers here so that non-finite or out-of-domain values are
//! rejected at the boundary rather than silently producing `NaN`-laden
//! results downstream.

use thiserror::Error;

/// Errors raised while validating gear-train inputs or solving the
/// kinematic relations.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum EpicyclicError {
    /// A tooth count was not a strictly positive integer.
    ///
    /// Physical gears need at least one tooth; the kinematic ratios
    /// divide by tooth counts, so zero or negative values are rejected.
    #[error("non-positive tooth count `{name}` = {value} (must be >= 1)")]
    NonPositiveTeeth {
        /// Which tooth count failed (`"sun"`, `"planet"`, `"ring"`).
        name: &'static str,
        /// The offending value.
        value: i64,
    },

    /// The standard coaxial meshing constraint `ring = sun + 2*planet`
    /// is violated, so the proposed sun/planet/ring teeth cannot form a
    /// concentric simple planetary set.
    #[error(
        "meshing constraint violated: ring {ring} != sun {sun} + 2*planet {planet} \
         (expected ring = {expected})"
    )]
    MeshingConstraint {
        /// Sun tooth count.
        sun: u32,
        /// Planet tooth count.
        planet: u32,
        /// Ring (annulus) tooth count supplied.
        ring: u32,
        /// The value `sun + 2*planet` that `ring` should have equalled.
        expected: u32,
    },

    /// A supplied angular velocity (or other real input) was not finite.
    #[error("non-finite value for `{name}`")]
    NonFinite {
        /// Which quantity failed (e.g. `"omega_sun"`).
        name: &'static str,
    },

    /// The Willis relation became singular for the requested unknown,
    /// i.e. the coefficient that multiplies the unknown vanished, so no
    /// unique solution exists.
    #[error("singular Willis relation: {reason}")]
    Singular {
        /// Human-readable explanation of why the solve is degenerate.
        reason: &'static str,
    },
}

impl EpicyclicError {
    /// Stable, machine-readable identifier for this error, suitable for
    /// logging or matching in tests without depending on the `Display`
    /// wording.
    pub fn code(&self) -> &'static str {
        match self {
            EpicyclicError::NonPositiveTeeth { .. } => "epicyclic.non_positive_teeth",
            EpicyclicError::MeshingConstraint { .. } => "epicyclic.meshing_constraint",
            EpicyclicError::NonFinite { .. } => "epicyclic.non_finite",
            EpicyclicError::Singular { .. } => "epicyclic.singular",
        }
    }
}

/// Validate that an integer tooth count is strictly positive and return
/// it as a `u32`.
///
/// # Errors
///
/// Returns [`EpicyclicError::NonPositiveTeeth`] when `value < 1`.
pub fn require_positive_teeth(name: &'static str, value: i64) -> Result<u32, EpicyclicError> {
    if value < 1 {
        return Err(EpicyclicError::NonPositiveTeeth { name, value });
    }
    Ok(value as u32)
}

/// Validate that a real input is finite (neither infinite nor `NaN`).
///
/// # Errors
///
/// Returns [`EpicyclicError::NonFinite`] when `value` is not finite.
pub fn require_finite(name: &'static str, value: f64) -> Result<f64, EpicyclicError> {
    if !value.is_finite() {
        return Err(EpicyclicError::NonFinite { name });
    }
    Ok(value)
}
