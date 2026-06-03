//! PBR material.

use serde::{Deserialize, Serialize};

/// Stable material identifier — string for editor friendliness.
pub type MaterialId = String;

/// Disney-style PBR material approximation.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Material {
    /// Display name.
    pub name: String,
    /// Linear-RGB base / diffuse colour.
    pub diffuse_color: [f32; 3],
    /// Linear-RGB specular tint.
    pub specular_color: [f32; 3],
    /// Surface roughness in `[0, 1]`.
    pub roughness: f32,
    /// Metallic in `[0, 1]`.
    pub metallic: f32,
    /// Index of refraction (1.0 = air, 1.5 = glass, …).
    pub ior: f32,
    /// Linear-RGB emissive colour (W/m²·sr-equivalent).
    pub emissive: [f32; 3],
}

impl Default for Material {
    /// Sensible default: matte white plastic.
    fn default() -> Self {
        Self {
            name: "default".into(),
            diffuse_color: [0.8, 0.8, 0.8],
            specular_color: [0.04, 0.04, 0.04],
            roughness: 0.6,
            metallic: 0.0,
            ior: 1.5,
            emissive: [0.0, 0.0, 0.0],
        }
    }
}

impl Material {
    /// Build a named matte material.
    pub fn matte(name: impl Into<String>, rgb: [f32; 3]) -> Self {
        Self {
            name: name.into(),
            diffuse_color: rgb,
            ..Default::default()
        }
    }

    /// Build a chrome-like polished metal.
    pub fn polished_metal(name: impl Into<String>, tint: [f32; 3]) -> Self {
        Self {
            name: name.into(),
            diffuse_color: [0.0, 0.0, 0.0],
            specular_color: tint,
            roughness: 0.05,
            metallic: 1.0,
            ior: 1.5,
            emissive: [0.0, 0.0, 0.0],
        }
    }

    /// Build a glass-like dielectric.
    pub fn glass(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            diffuse_color: [0.95, 0.95, 0.95],
            specular_color: [0.04, 0.04, 0.04],
            roughness: 0.0,
            metallic: 0.0,
            ior: 1.52,
            emissive: [0.0, 0.0, 0.0],
        }
    }
}
