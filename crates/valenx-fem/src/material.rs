//! FEM material model + a hardcoded library of common engineering
//! materials.

use serde::{Deserialize, Serialize};

/// Elastoplastic properties — the extra constants the
/// [`crate::plasticity`] J2 solver needs on top of the elastic ones.
///
/// A material with `plasticity = None` is purely elastic; supplying a
/// [`PlasticProperties`] enables the von Mises radial-return update.
/// The hardening law is **linear isotropic**: the yield surface
/// expands uniformly as
///
/// ```text
///   σ_y(ε̄ᵖ) = σ_y0 + H · ε̄ᵖ
/// ```
///
/// with `ε̄ᵖ` the accumulated equivalent plastic strain and `H` the
/// (constant) plastic hardening modulus.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
pub struct PlasticProperties {
    /// Initial yield stress `σ_y0` in Pa — the von Mises stress at
    /// which the material first yields.
    pub yield_stress: f64,
    /// Linear isotropic hardening modulus `H` in Pa — the slope of the
    /// yield stress against accumulated plastic strain. `H = 0` is
    /// perfect (non-hardening) plasticity.
    pub hardening_modulus: f64,
}

impl PlasticProperties {
    /// The current yield stress at an accumulated equivalent plastic
    /// strain `eq_plastic_strain` — `σ_y0 + H·ε̄ᵖ`.
    #[inline]
    pub fn yield_stress_at(&self, eq_plastic_strain: f64) -> f64 {
        self.yield_stress + self.hardening_modulus * eq_plastic_strain.max(0.0)
    }
}

/// Isotropic linear-elastic material properties — with optional
/// [`PlasticProperties`] for the J2 plasticity solver.
///
/// The four core fields are the ones every FEA solver needs for linear
/// static analysis: Young's modulus, Poisson's ratio, density, and
/// thermal conductivity. The optional `plasticity` field carries the
/// yield stress + hardening modulus the [`crate::plasticity`] solver
/// reads; it defaults to `None` (a purely elastic material), so every
/// existing elastic / modal / thermal solve is unaffected.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct FemMaterial {
    /// User-facing material name. Used as the input-deck section
    /// label (`*MATERIAL, NAME=Steel` in CalculiX).
    pub name: String,
    /// Young's modulus E in Pa.
    pub youngs_modulus: f64,
    /// Poisson's ratio (dimensionless, typically 0.0..=0.5).
    pub poisson_ratio: f64,
    /// Density in kg/m^3.
    pub density: f64,
    /// Thermal conductivity in W/(m·K).
    pub thermal_conductivity: f64,
    /// Optional elastoplastic constants. `None` → a purely elastic
    /// material. `Some(_)` enables the J2 radial-return plasticity
    /// update in [`crate::plasticity`].
    #[serde(default)]
    pub plasticity: Option<PlasticProperties>,
}

impl Default for FemMaterial {
    fn default() -> Self {
        Self {
            name: "Steel_AISI_1045".into(),
            youngs_modulus: 205e9,
            poisson_ratio: 0.29,
            density: 7850.0,
            thermal_conductivity: 49.8,
            plasticity: None,
        }
    }
}

/// Hardcoded library of common engineering materials. The desktop
/// shell renders these in a combobox; the user can also author a
/// custom material that overrides the library entries.
pub fn material_library() -> Vec<FemMaterial> {
    vec![
        // Steel AISI 1045 — generic medium-carbon engineering steel.
        // Yield ≈ 530 MPa, a representative post-yield hardening slope.
        FemMaterial {
            plasticity: Some(PlasticProperties {
                yield_stress: 530e6,
                hardening_modulus: 2.0e9,
            }),
            ..FemMaterial::default()
        },
        // Aluminium 6061-T6 — common aerospace / structural alloy.
        FemMaterial {
            name: "Aluminium_6061_T6".into(),
            youngs_modulus: 68.9e9,
            poisson_ratio: 0.33,
            density: 2700.0,
            thermal_conductivity: 167.0,
            // Yield ≈ 276 MPa for T6 temper.
            plasticity: Some(PlasticProperties {
                yield_stress: 276e6,
                hardening_modulus: 0.7e9,
            }),
        },
        // Titanium Ti-6Al-4V — grade 5 titanium alloy.
        FemMaterial {
            name: "Titanium_Ti6Al4V".into(),
            youngs_modulus: 113.8e9,
            poisson_ratio: 0.342,
            density: 4429.0,
            thermal_conductivity: 6.7,
            // Yield ≈ 880 MPa.
            plasticity: Some(PlasticProperties {
                yield_stress: 880e6,
                hardening_modulus: 1.0e9,
            }),
        },
        // ABS plastic — generic injection-moulding plastic.
        FemMaterial {
            name: "ABS_Plastic".into(),
            youngs_modulus: 2.3e9,
            poisson_ratio: 0.35,
            density: 1050.0,
            thermal_conductivity: 0.17,
            plasticity: None,
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn library_includes_steel_aluminium_titanium_abs() {
        let lib = material_library();
        let names: Vec<&str> = lib.iter().map(|m| m.name.as_str()).collect();
        assert!(names.contains(&"Steel_AISI_1045"));
        assert!(names.contains(&"Aluminium_6061_T6"));
        assert!(names.contains(&"Titanium_Ti6Al4V"));
        assert!(names.contains(&"ABS_Plastic"));
    }

    #[test]
    fn default_material_is_steel() {
        let m = FemMaterial::default();
        assert_eq!(m.name, "Steel_AISI_1045");
        assert!((m.youngs_modulus - 205e9).abs() < 1e-3);
        assert!((m.poisson_ratio - 0.29).abs() < 1e-9);
    }

    #[test]
    fn material_serializes_to_json() {
        let m = FemMaterial::default();
        let j = serde_json::to_string(&m).unwrap();
        let m2: FemMaterial = serde_json::from_str(&j).unwrap();
        assert_eq!(m, m2);
    }
}
