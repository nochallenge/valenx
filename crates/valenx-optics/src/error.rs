//! Geometric-optics error taxonomy.
//!
//! Every fallible constructor and computation in this crate funnels its
//! failure modes through [`OpticsError`]. The variants are deliberately
//! coarse — a bad physical input (a non-positive refractive index, an
//! out-of-range angle), a singular geometry (a focal length of zero, an
//! object placed exactly at the front focal point so the image runs off
//! to infinity), or a physically impossible request (asking for the
//! critical angle of refraction when `n1 <= n2`, where total internal
//! reflection cannot occur).

use thiserror::Error;

/// Errors raised by the geometric-optics models.
#[derive(Debug, Error)]
pub enum OpticsError {
    /// A physical input fell outside its valid domain.
    ///
    /// Examples: a refractive index that is not strictly greater than
    /// zero (vacuum is the floor, `n = 1`, and real media are `> 1`),
    /// or an incidence angle outside `[0, 90]` degrees.
    #[error("invalid parameter `{name}` = {value}: {reason}")]
    InvalidParameter {
        /// Name of the offending parameter.
        name: &'static str,
        /// The numeric value that was supplied.
        value: f64,
        /// Human-readable explanation of why it was rejected.
        reason: &'static str,
    },

    /// Total internal reflection occurs, so no refracted ray exists.
    ///
    /// Raised by [`crate::refraction::refract_angle`] when the ray
    /// travels into a less-dense medium (`n1 > n2`) at an incidence
    /// angle at or beyond the critical angle. Snell's law would demand
    /// `sin t2 > 1`, which has no real solution: the interface acts as a
    /// perfect mirror.
    #[error(
        "total internal reflection: incidence {incidence_deg} deg >= critical {critical_deg} deg \
         (n1 = {n1} into n2 = {n2}); no refracted ray exists"
    )]
    TotalInternalReflection {
        /// Refractive index of the incident (denser) medium.
        n1: f64,
        /// Refractive index of the transmitting (rarer) medium.
        n2: f64,
        /// The incidence angle that triggered TIR, in degrees.
        incidence_deg: f64,
        /// The critical angle for this interface, in degrees.
        critical_deg: f64,
    },

    /// The requested image is at infinity (or the geometry is singular).
    ///
    /// Raised by the thin-lens solvers when the object sits exactly at
    /// the front focal point (`do = f`), so `1/di = 1/f - 1/do = 0` and
    /// the rays emerge collimated — the image distance diverges and the
    /// magnification is unbounded. Also covers a focal length of exactly
    /// zero, for which the thin-lens equation degenerates.
    #[error("singular image geometry: {reason}")]
    SingularGeometry {
        /// Human-readable explanation of the singularity.
        reason: &'static str,
    },
}

/// Coarse category for an [`OpticsError`], useful for UI grouping.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A user-supplied physical input was out of range.
    Input,
    /// A physically valid but special-case geometry (TIR, image at
    /// infinity) that has no ordinary numeric answer.
    Physical,
}

impl OpticsError {
    /// Stable, dot-separated identifier for this error.
    ///
    /// The string is part of the crate's API surface: callers may match
    /// on it for logging or localisation without depending on the
    /// `Display` text, which is free to change.
    pub fn code(&self) -> &'static str {
        match self {
            OpticsError::InvalidParameter { .. } => "optics.invalid_parameter",
            OpticsError::TotalInternalReflection { .. } => "optics.total_internal_reflection",
            OpticsError::SingularGeometry { .. } => "optics.singular_geometry",
        }
    }

    /// Coarse [`ErrorCategory`] for this error.
    pub fn category(&self) -> ErrorCategory {
        match self {
            OpticsError::InvalidParameter { .. } => ErrorCategory::Input,
            OpticsError::TotalInternalReflection { .. } | OpticsError::SingularGeometry { .. } => {
                ErrorCategory::Physical
            }
        }
    }
}

/// Validate that a refractive index is physically admissible.
///
/// Refractive indices in this crate must be finite and strictly greater
/// than zero. (Real transparent media have `n >= 1`, but values in
/// `(0, 1)` are accepted so the same checker serves exotic / effective
/// indices; the floor that actually matters for the models is `> 0` so
/// that ratios and `asin` arguments stay well defined.)
///
/// # Errors
///
/// Returns [`OpticsError::InvalidParameter`] if `n` is not finite or is
/// not strictly positive.
pub fn validate_index(name: &'static str, n: f64) -> Result<f64, OpticsError> {
    if !n.is_finite() {
        return Err(OpticsError::InvalidParameter {
            name,
            value: n,
            reason: "refractive index must be finite",
        });
    }
    if n <= 0.0 {
        return Err(OpticsError::InvalidParameter {
            name,
            value: n,
            reason: "refractive index must be strictly positive",
        });
    }
    Ok(n)
}

/// Validate that an angle (in degrees) lies in the physical range
/// `[0, 90]` for an angle measured from the surface normal.
///
/// # Errors
///
/// Returns [`OpticsError::InvalidParameter`] if `deg` is not finite or
/// falls outside the closed interval `[0, 90]`.
pub fn validate_angle_deg(name: &'static str, deg: f64) -> Result<f64, OpticsError> {
    if !deg.is_finite() {
        return Err(OpticsError::InvalidParameter {
            name,
            value: deg,
            reason: "angle must be finite",
        });
    }
    if !(0.0..=90.0).contains(&deg) {
        return Err(OpticsError::InvalidParameter {
            name,
            value: deg,
            reason: "angle from the normal must be within [0, 90] degrees",
        });
    }
    Ok(deg)
}
