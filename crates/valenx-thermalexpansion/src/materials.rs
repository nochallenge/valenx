//! A small library of representative material properties.
//!
//! Each [`Material`] bundles a *representative* near-room-temperature
//! linear coefficient of thermal expansion and Young's modulus so the
//! expansion and stress formulas can be driven from a name. The values are
//! rounded textbook figures for teaching and rough estimation, not
//! certified design data; real coefficients vary with temperature, alloy,
//! temper and processing.

use crate::error::ThermalError;
use crate::expansion::LinearCoefficient;
use crate::stress::YoungsModulus;
use serde::{Deserialize, Serialize};

/// A named material with representative thermal / elastic properties.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Material {
    /// Lower-case key used by [`lookup`].
    pub name: &'static str,
    /// Linear coefficient of thermal expansion, in 1/K.
    pub alpha_per_kelvin: f64,
    /// Young's modulus, in pascals (Pa).
    pub youngs_modulus_pa: f64,
}

impl Material {
    /// The validated [`LinearCoefficient`] for this material.
    ///
    /// # Errors
    ///
    /// Returns [`ThermalError::NonPositive`] / [`ThermalError::NonFinite`]
    /// only if the stored value were invalid; the built-in table is all
    /// positive and finite, so this never fails for library entries.
    pub fn alpha(&self) -> Result<LinearCoefficient, ThermalError> {
        LinearCoefficient::new(self.alpha_per_kelvin)
    }

    /// The validated [`YoungsModulus`] for this material.
    ///
    /// # Errors
    ///
    /// As [`Material::alpha`]: never fails for the built-in table.
    pub fn youngs_modulus(&self) -> Result<YoungsModulus, ThermalError> {
        YoungsModulus::new(self.youngs_modulus_pa)
    }
}

/// The built-in material table.
///
/// Keys are lower-case. Values are representative room-temperature
/// figures: CTE in 1/K and Young's modulus in Pa.
pub const LIBRARY: &[Material] = &[
    Material {
        name: "aluminium",
        alpha_per_kelvin: 23.1e-6,
        youngs_modulus_pa: 69.0e9,
    },
    Material {
        name: "steel",
        alpha_per_kelvin: 12.0e-6,
        youngs_modulus_pa: 200.0e9,
    },
    Material {
        name: "copper",
        alpha_per_kelvin: 16.5e-6,
        youngs_modulus_pa: 117.0e9,
    },
    Material {
        name: "invar",
        alpha_per_kelvin: 1.2e-6,
        youngs_modulus_pa: 141.0e9,
    },
    Material {
        name: "fused-silica",
        alpha_per_kelvin: 0.55e-6,
        youngs_modulus_pa: 73.0e9,
    },
];

/// Look up a material by (case-insensitive) name in the built-in
/// [`LIBRARY`].
///
/// # Errors
///
/// Returns [`ThermalError::UnknownMaterial`] if no entry matches `name`.
pub fn lookup(name: &str) -> Result<Material, ThermalError> {
    let key = name.trim().to_ascii_lowercase();
    LIBRARY
        .iter()
        .find(|m| m.name == key)
        .copied()
        .ok_or(ThermalError::UnknownMaterial { name: key })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expansion::linear_expansion;

    const EPS: f64 = 1e-9;

    #[test]
    fn lookup_is_case_insensitive_and_trims() {
        let a = lookup("Aluminium").unwrap();
        let b = lookup("  aluminium ").unwrap();
        assert_eq!(a, b);
        assert_eq!(a.name, "aluminium");
    }

    #[test]
    fn lookup_unknown_errors() {
        let err = lookup("unobtainium").unwrap_err();
        assert!(matches!(err, ThermalError::UnknownMaterial { .. }));
        assert_eq!(err.code(), "thermalexpansion.unknown-material");
    }

    #[test]
    fn all_library_entries_are_constructible() {
        // Every stored property must pass the newtype validators.
        for m in LIBRARY {
            let alpha = m.alpha().expect("library alpha valid");
            let e = m.youngs_modulus().expect("library E valid");
            assert!(alpha.per_kelvin() > 0.0);
            assert!(e.pascals() > 0.0);
        }
    }

    #[test]
    fn invar_expands_far_less_than_aluminium() {
        // Physical sanity: same bar, same heating, low-CTE Invar moves a
        // small fraction of what aluminium does (its whole reason to exist).
        let dt = 100.0;
        let l0 = 1.0;
        let al = linear_expansion(lookup("aluminium").unwrap().alpha().unwrap(), l0, dt).unwrap();
        let invar = linear_expansion(lookup("invar").unwrap().alpha().unwrap(), l0, dt).unwrap();
        assert!(invar > 0.0 && al > 0.0);
        assert!(invar < al / 10.0, "Invar should be <1/10 of aluminium");
    }

    #[test]
    fn library_alpha_drives_the_formula() {
        let steel = lookup("steel").unwrap();
        let dl = linear_expansion(steel.alpha().unwrap(), 10.0, 100.0).unwrap();
        let expected = 12.0e-6 * 10.0 * 100.0;
        assert!((dl - expected).abs() < EPS, "got {dl}, expected {expected}");
    }
}
