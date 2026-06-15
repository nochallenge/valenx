//! Soil material properties: friction angle, cohesion, unit weight.
//!
//! [`SoilProperties`] bundles the three soil parameters the Terzaghi
//! equation consumes. Build it through [`SoilProperties::new`], which
//! validates each input against its physical domain.

use serde::{Deserialize, Serialize};

use crate::error::SoilBearingError;

/// Drained soil parameters for a bearing-capacity calculation.
///
/// Units are the caller's responsibility but must be self-consistent
/// with the rest of the calculation (see the crate-level
/// [documentation](crate) for the recommended kN / m / kPa system):
///
/// - `friction_angle_deg`: effective angle of internal friction, degrees.
/// - `cohesion`: effective cohesion `c` (a stress, e.g. kPa).
/// - `unit_weight`: effective unit weight `gamma` (force/volume, e.g. kN/m^3).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SoilProperties {
    friction_angle_deg: f64,
    cohesion: f64,
    unit_weight: f64,
}

impl SoilProperties {
    /// Build validated soil properties.
    ///
    /// # Parameters
    ///
    /// - `friction_angle_deg` must be finite and in `[0, 90)` degrees.
    ///   `90` (and above) is rejected because the bearing-capacity
    ///   factors diverge as `phi -> 90 deg`.
    /// - `cohesion` must be finite and non-negative.
    /// - `unit_weight` must be finite and non-negative.
    ///
    /// # Errors
    ///
    /// Returns [`SoilBearingError::NotFinite`] for any NaN/infinite
    /// input, or [`SoilBearingError::InvalidParameter`] for an
    /// out-of-domain value.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_soilbearing::SoilProperties;
    ///
    /// let sand = SoilProperties::new(32.0, 0.0, 18.5).unwrap();
    /// assert!((sand.friction_angle_deg() - 32.0).abs() < 1e-12);
    ///
    /// assert!(SoilProperties::new(90.0, 0.0, 18.0).is_err());
    /// assert!(SoilProperties::new(30.0, -1.0, 18.0).is_err());
    /// ```
    pub fn new(
        friction_angle_deg: f64,
        cohesion: f64,
        unit_weight: f64,
    ) -> Result<Self, SoilBearingError> {
        let friction_angle_deg = SoilBearingError::require_finite("phi_deg", friction_angle_deg)?;
        let cohesion = SoilBearingError::require_finite("cohesion", cohesion)?;
        let unit_weight = SoilBearingError::require_finite("unit_weight", unit_weight)?;

        if !(0.0..90.0).contains(&friction_angle_deg) {
            return Err(SoilBearingError::invalid(
                "phi_deg",
                friction_angle_deg,
                "friction angle must be in the half-open range [0, 90) degrees",
            ));
        }
        if cohesion < 0.0 {
            return Err(SoilBearingError::invalid(
                "cohesion",
                cohesion,
                "cohesion must be non-negative",
            ));
        }
        if unit_weight < 0.0 {
            return Err(SoilBearingError::invalid(
                "unit_weight",
                unit_weight,
                "unit weight must be non-negative",
            ));
        }

        Ok(SoilProperties {
            friction_angle_deg,
            cohesion,
            unit_weight,
        })
    }

    /// Effective angle of internal friction, in degrees.
    pub fn friction_angle_deg(&self) -> f64 {
        self.friction_angle_deg
    }

    /// Effective angle of internal friction, in radians.
    pub fn friction_angle_rad(&self) -> f64 {
        self.friction_angle_deg.to_radians()
    }

    /// Effective cohesion `c` (stress units).
    pub fn cohesion(&self) -> f64 {
        self.cohesion
    }

    /// Effective unit weight `gamma` (force-per-volume units).
    pub fn unit_weight(&self) -> f64 {
        self.unit_weight
    }

    /// Whether this soil is cohesionless (`c = 0`), e.g. clean sand.
    ///
    /// In a cohesionless soil the `c * Nc` term of the Terzaghi
    /// equation vanishes regardless of `Nc`.
    pub fn is_cohesionless(&self) -> bool {
        self.cohesion == 0.0
    }
}
