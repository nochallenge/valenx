//! Material fracture properties — Mode-I fracture toughness and yield
//! strength, with a validated constructor.

use crate::error::{FractureError, Result};
use serde::{Deserialize, Serialize};

/// A material described by the two properties the LEFM calculators need:
/// its plane-strain Mode-I **fracture toughness** `K_Ic` and its tensile
/// **yield strength** `σ_y`.
///
/// Units are not fixed by the type, but every field must share a
/// consistent set: if `K_Ic` is in `MPa·√m` then `σ_y` must be in `MPa`
/// and crack lengths in `m`. The shipped test values use that SI-ish set
/// (e.g. a 7075-T6 aluminium alloy: `K_Ic ≈ 24 MPa·√m`, `σ_y ≈ 470 MPa`).
///
/// Construct with [`Material::new`], which rejects non-positive or
/// non-finite inputs.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Material {
    /// Plane-strain Mode-I fracture toughness `K_Ic` (e.g. `MPa·√m`).
    /// Strictly positive.
    pub fracture_toughness: f64,
    /// Tensile yield strength `σ_y` (e.g. `MPa`). Strictly positive.
    pub yield_strength: f64,
}

impl Material {
    /// Build a [`Material`] from a fracture toughness and a yield strength.
    ///
    /// # Errors
    ///
    /// Returns [`FractureError::NonPositive`] if either argument is not a
    /// finite, strictly-positive number.
    ///
    /// ```
    /// use valenx_fracture::Material;
    /// let alu = Material::new(24.0, 470.0).unwrap();
    /// assert_eq!(alu.fracture_toughness, 24.0);
    /// ```
    pub fn new(fracture_toughness: f64, yield_strength: f64) -> Result<Self> {
        require_positive("fracture_toughness", fracture_toughness)?;
        require_positive("yield_strength", yield_strength)?;
        Ok(Self {
            fracture_toughness,
            yield_strength,
        })
    }
}

/// Validate that `value` is finite and strictly positive, returning a
/// [`FractureError::NonPositive`] tagged with `name` otherwise.
///
/// Shared by the material constructor and the crack-mechanics functions so
/// the "must be > 0" contract reads the same everywhere.
pub(crate) fn require_positive(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(FractureError::NonPositive { name, value })
    }
}

/// Validate that `value` is finite and non-negative, returning a
/// [`FractureError::InvalidLength`] tagged with `name` otherwise.
///
/// Used for crack lengths and plate dimensions, which may legitimately be
/// zero in some expressions but never negative or `NaN`.
pub(crate) fn require_non_negative_length(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(FractureError::InvalidLength { name, value })
    }
}

/// Validate that `value` is finite and non-negative, returning a
/// [`FractureError::InvalidStress`] tagged with `name` otherwise.
///
/// The Mode-I formulae assume an opening (non-negative) far-field tension.
pub(crate) fn require_non_negative_stress(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(FractureError::InvalidStress { name, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_material_round_trips_fields() {
        let m = Material::new(24.0, 470.0).expect("valid");
        assert!((m.fracture_toughness - 24.0).abs() < 1e-12);
        assert!((m.yield_strength - 470.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_non_positive_toughness() {
        let err = Material::new(0.0, 470.0).unwrap_err();
        assert_eq!(
            err,
            FractureError::NonPositive {
                name: "fracture_toughness",
                value: 0.0
            }
        );
    }

    #[test]
    fn rejects_negative_yield() {
        let err = Material::new(24.0, -1.0).unwrap_err();
        assert_eq!(err.code(), "fracture.non_positive");
    }

    #[test]
    fn rejects_non_finite() {
        assert!(Material::new(f64::NAN, 470.0).is_err());
        assert!(Material::new(24.0, f64::INFINITY).is_err());
    }
}
