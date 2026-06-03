//! Render engines.

use serde::{Deserialize, Serialize};

/// Which renderer to dispatch to.
#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub enum RenderEngine {
    /// LuxCoreRender / LuxRender — PBRT-style physically-based path
    /// tracer.
    LuxRender,
    /// Cycles — Blender's path tracer (standalone build).
    Cycles,
    /// POV-Ray — Persistence of Vision SDL-based renderer.
    PovRay,
    /// In-app egui+wgpu viewport screenshot (no subprocess).
    #[default]
    Native,
}

impl RenderEngine {
    /// UI dropdown label.
    pub fn label(self) -> &'static str {
        match self {
            Self::LuxRender => "LuxRender",
            Self::Cycles => "Cycles",
            Self::PovRay => "POV-Ray",
            Self::Native => "Native (viewport screenshot)",
        }
    }

    /// Canonical file extension for the engine's scene file.
    pub fn scene_file_extension(self) -> &'static str {
        match self {
            Self::LuxRender => "lxs",
            Self::Cycles => "xml",
            Self::PovRay => "pov",
            Self::Native => "png",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_engine_is_native() {
        assert_eq!(RenderEngine::default(), RenderEngine::Native);
    }

    #[test]
    fn every_variant_has_a_label_and_extension() {
        // One assertion per variant so all four `match` arms of both
        // accessors are exercised.
        for (engine, label, ext) in [
            (RenderEngine::LuxRender, "LuxRender", "lxs"),
            (RenderEngine::Cycles, "Cycles", "xml"),
            (RenderEngine::PovRay, "POV-Ray", "pov"),
            (RenderEngine::Native, "Native (viewport screenshot)", "png"),
        ] {
            assert_eq!(engine.label(), label);
            assert_eq!(engine.scene_file_extension(), ext);
        }
    }

    #[test]
    fn engine_round_trips_through_ron() {
        for engine in [
            RenderEngine::LuxRender,
            RenderEngine::Cycles,
            RenderEngine::PovRay,
            RenderEngine::Native,
        ] {
            let ron = ron::to_string(&engine).unwrap();
            let back: RenderEngine = ron::from_str(&ron).unwrap();
            assert_eq!(engine, back);
        }
    }

    #[test]
    fn engine_is_copy_and_comparable() {
        let a = RenderEngine::Cycles;
        let b = a; // Copy
        assert_eq!(a, b);
        assert_ne!(RenderEngine::Cycles, RenderEngine::PovRay);
        assert!(format!("{a:?}").contains("Cycles"));
    }
}
