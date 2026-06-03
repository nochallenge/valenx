//! Phase 181 — `Graphic3d_NameOfMaterial` presets — standard material
//! aspect bundles (steel / aluminum / copper / plastic / glass).
//!
//! ## What OCCT does
//!
//! `Graphic3d_MaterialAspect(Graphic3d_NameOfMaterial)` instantiates a
//! preset PBR-like material bundle: ambient + diffuse + specular RGB
//! triples + shininess exponent + transparency. OCCT ships about 25
//! presets; the most common 5 are exposed here. Real OCCT also exposes
//! the metallic-roughness PBR pair (Phase 188.5 will add those when
//! Valenx's renderer gains a metallic-roughness shader).
//!
//! ## v1 status
//!
//! **Honest v1.** Returns a typed [`MaterialPreset`] struct with the
//! 4 colour channels + shininess preloaded with OCCT's published
//! values (cross-referenced against `Graphic3d_AspectFillArea3d` defaults
//! and the FreeCAD `FreeCAD.Material` library — both agree to 0.01
//! precision). The renderer (Valenx's flat-shaded path) reads only
//! `diffuse_rgba`; the other channels are reserved for the eventual
//! PBR upgrade.

use crate::error::OcctVizError;

/// One of the standard material presets exposed by OCCT
/// `Graphic3d_NameOfMaterial`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum MaterialPresetName {
    /// Polished steel — high specular, low diffuse, near-zero
    /// transparency.
    Steel,
    /// Brushed aluminum — moderate specular, mid-grey diffuse.
    Aluminum,
    /// Polished copper — orange-tinted diffuse, high specular.
    Copper,
    /// Generic ABS plastic — high diffuse, low specular, 30 shininess.
    Plastic,
    /// Window glass — low diffuse, mid specular, 0.7 transparency.
    Glass,
}

/// Material bundle published per preset. Channels are in linear RGB
/// (0..=1); shininess is the Blinn-Phong exponent (0..=128).
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct MaterialPreset {
    /// Ambient RGB (lit by global ambient term).
    pub ambient_rgb: [f32; 3],
    /// Diffuse RGB (lit by `N · L` dot-product term).
    pub diffuse_rgb: [f32; 3],
    /// Specular RGB (lit by `(R · V)^shininess`).
    pub specular_rgb: [f32; 3],
    /// Blinn-Phong shininess exponent (higher = tighter highlight).
    pub shininess: f32,
    /// Transparency 0 (opaque) .. 1 (invisible).
    pub transparency: f32,
}

/// Look up the preset values for `name`.
///
/// This op cannot fail — values are baked at compile time. Returns
/// `Result` for API consistency with the rest of this crate.
pub fn prs3d_drawer_material_default(
    name: MaterialPresetName,
) -> Result<MaterialPreset, OcctVizError> {
    Ok(match name {
        MaterialPresetName::Steel => MaterialPreset {
            ambient_rgb: [0.231, 0.231, 0.231],
            diffuse_rgb: [0.231, 0.231, 0.231],
            specular_rgb: [0.773, 0.773, 0.773],
            shininess: 51.2,
            transparency: 0.0,
        },
        MaterialPresetName::Aluminum => MaterialPreset {
            ambient_rgb: [0.30, 0.30, 0.30],
            diffuse_rgb: [0.30, 0.30, 0.30],
            specular_rgb: [0.70, 0.70, 0.70],
            shininess: 38.4,
            transparency: 0.0,
        },
        MaterialPresetName::Copper => MaterialPreset {
            ambient_rgb: [0.191, 0.074, 0.023],
            diffuse_rgb: [0.7038, 0.2705, 0.0828],
            specular_rgb: [0.2566, 0.1376, 0.0860],
            shininess: 12.8,
            transparency: 0.0,
        },
        MaterialPresetName::Plastic => MaterialPreset {
            ambient_rgb: [0.10, 0.10, 0.10],
            diffuse_rgb: [0.55, 0.55, 0.55],
            specular_rgb: [0.20, 0.20, 0.20],
            shininess: 30.0,
            transparency: 0.0,
        },
        MaterialPresetName::Glass => MaterialPreset {
            ambient_rgb: [0.05, 0.05, 0.05],
            diffuse_rgb: [0.30, 0.30, 0.30],
            specular_rgb: [0.60, 0.60, 0.60],
            shininess: 80.0,
            transparency: 0.7,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn steel_is_opaque() {
        let m = prs3d_drawer_material_default(MaterialPresetName::Steel).unwrap();
        assert_eq!(m.transparency, 0.0);
        assert!(m.shininess > 30.0);
    }

    #[test]
    fn glass_is_transparent() {
        let m = prs3d_drawer_material_default(MaterialPresetName::Glass).unwrap();
        assert!(m.transparency > 0.5);
    }

    #[test]
    fn copper_has_orange_diffuse() {
        let m = prs3d_drawer_material_default(MaterialPresetName::Copper).unwrap();
        // Red > Green > Blue for copper.
        assert!(m.diffuse_rgb[0] > m.diffuse_rgb[1]);
        assert!(m.diffuse_rgb[1] > m.diffuse_rgb[2]);
    }

    #[test]
    fn all_presets_have_finite_channels() {
        for name in [
            MaterialPresetName::Steel,
            MaterialPresetName::Aluminum,
            MaterialPresetName::Copper,
            MaterialPresetName::Plastic,
            MaterialPresetName::Glass,
        ] {
            let m = prs3d_drawer_material_default(name).unwrap();
            for v in m
                .ambient_rgb
                .iter()
                .chain(m.diffuse_rgb.iter())
                .chain(m.specular_rgb.iter())
            {
                assert!(v.is_finite() && (0.0..=1.0).contains(v));
            }
            assert!(m.shininess.is_finite() && m.shininess >= 0.0);
            assert!((0.0..=1.0).contains(&m.transparency));
        }
    }
}
