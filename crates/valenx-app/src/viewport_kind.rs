//! [`ViewportKind`] — the palette of viewport implementations the
//! central panel can host.
//!
//! # Extension guide
//!
//! To add a new domain viewport:
//! 1. Add a variant here.
//! 2. Create a `viewport_<name>` module under `crates/valenx-app/src/`.
//! 3. Wire the dispatch in the `CentralPanel` closure in `update.rs`.
//! 4. Add a `for_<workbench>` constructor here if the workbench has a
//!    natural default viewport.

/// Which viewport implementation to render in the central panel.
///
/// The active kind is stored as `crate::ValenxApp::active_viewport` and
/// the user can change it via **View → Central viewport**. Each workbench
/// that has a preferred viewport kind sets it when it is first enabled.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ViewportKind {
    /// Standard Blender-style 3D orbit / pan / zoom render.
    ///
    /// Backed by the `wgpu` GPU renderer (if available) or the egui
    /// painter software fallback. Used by CAD, mesh, FEM, aero, and
    /// astro workbenches.
    #[default]
    Viewport3D,

    /// 2D annotated-sequence viewport for DNA / RNA work.
    ///
    /// Shows a linear feature-track map and / or a circular plasmid
    /// diagram. Driven by [`crate::viewport_2d`]; data is drawn from
    /// `valenx-bioseq`.
    Viewport2dDna,
}

impl ViewportKind {
    /// Human-readable label for menus and the status bar.
    pub fn label(self) -> &'static str {
        match self {
            ViewportKind::Viewport3D => "3D Viewport",
            ViewportKind::Viewport2dDna => "2D DNA / Plasmid Viewport",
        }
    }

    /// The viewport kind the Genetics Workbench prefers. Called when
    /// the user first enables the workbench so the central panel
    /// immediately shows a DNA-relevant view.
    pub fn for_genetics() -> Self {
        ViewportKind::Viewport2dDna
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_3d() {
        assert_eq!(ViewportKind::default(), ViewportKind::Viewport3D);
    }

    #[test]
    fn for_genetics_is_2d_dna() {
        assert_eq!(ViewportKind::for_genetics(), ViewportKind::Viewport2dDna);
    }

    #[test]
    fn labels_are_non_empty() {
        for kind in [ViewportKind::Viewport3D, ViewportKind::Viewport2dDna] {
            assert!(!kind.label().is_empty(), "{kind:?} has empty label");
        }
    }
}
