//! Phase 187 — `Prs3d_Drawer::SetHiddenLineAspect` — toggle hidden-
//! line visualization (HLR — Hidden Line Removal).
//!
//! ## What OCCT does
//!
//! `Prs3d_Drawer::SetTypeOfHLR(Prs3d_TypeOfHLR_PolyAlgo)` enables a
//! second wire-rendering pass that draws normally-hidden edges (the
//! ones behind opaque faces) as dashed lines on top of the regular
//! visible-edge layer. The HLR algorithm runs on the CPU
//! (`HLRBRep_Algo`) and emits two edge layers: `Visible` and `Hidden`,
//! each rendered with their own [`crate::prs3d_drawer_line_style()`]
//! aspect.
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 187.5). The crate already ships a
//! real hidden-line-removal pass — `valenx_techdraw::hlr::classify_edges`
//! (the Phase 5 HLR pipeline) — which projects a solid through a
//! camera matrix and partitions its edges into `visible` and `hidden`
//! sets via a depth-grid z-buffer. This function wraps that pass into
//! the OCCT `SetHiddenLineAspect` display semantics:
//!
//! - `enabled == true` — both edge layers are returned, so the
//!   renderer can draw the hidden edges (typically dashed) on top of
//!   the visible ones.
//! - `enabled == false` — only the `visible` layer is populated;
//!   `hidden` comes back empty (the OCCT default, hidden lines off).
//!
//! Edge segments are 2D drawing-plane millimetre coordinates, ready
//! to feed a wire-overlay pass.

use nalgebra::Matrix4;
use valenx_cad::Solid;

use crate::error::OcctVizError;

/// The two edge layers produced by the hidden-line display pass.
#[derive(Clone, Debug, Default)]
pub struct HiddenLineDisplay {
    /// Edges in front of every opaque face — always drawn.
    pub visible: Vec<[(f64, f64); 2]>,
    /// Edges occluded by an opaque face. Populated only when the
    /// caller passed `enabled = true`; the renderer draws these
    /// dashed.
    pub hidden: Vec<[(f64, f64); 2]>,
}

/// Compute the visible/hidden edge layers for `solid` viewed through
/// `camera`, honouring the `enabled` toggle.
///
/// `camera` is a world→clip 4x4 matrix (the product of the view and
/// projection matrices).
///
/// # Errors
///
/// - [`OcctVizError::Render`] when the underlying HLR pass fails
///   (empty solid, tessellation failure).
///
/// # Example
///
/// ```
/// use nalgebra::Matrix4;
/// use valenx_occt_viz::prs3d_drawer_hidden_line_display::prs3d_drawer_hidden_line_display;
/// let cube = valenx_cad::box_solid(10.0, 10.0, 10.0).unwrap();
/// // A perspective camera looking at the cube from the front.
/// let cam = Matrix4::new_perspective(1.0, 0.9, 0.1, 1000.0)
///     * Matrix4::new_translation(&nalgebra::Vector3::new(0.0, 0.0, -40.0));
/// let layers = prs3d_drawer_hidden_line_display(&cube, &cam, true).unwrap();
/// // With hidden lines enabled both layers can carry edges.
/// assert!(!layers.visible.is_empty());
/// ```
pub fn prs3d_drawer_hidden_line_display(
    solid: &Solid,
    camera: &Matrix4<f64>,
    enabled: bool,
) -> Result<HiddenLineDisplay, OcctVizError> {
    let (visible, hidden) = valenx_techdraw::hlr::classify_edges(solid, camera)
        .map_err(|e| OcctVizError::render(format!("hidden-line removal: {e}")))?;

    Ok(HiddenLineDisplay {
        visible,
        // OCCT's SetHiddenLineAspect default is "off" — when disabled
        // the hidden layer is simply not surfaced for drawing.
        hidden: if enabled { hidden } else { Vec::new() },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enabled_returns_both_layers() {
        // A cube viewed through a perspective-ish camera has both
        // front-facing and occluded edges.
        let cube = valenx_cad::box_solid(10.0, 10.0, 10.0).unwrap();
        let cam = Matrix4::new_perspective(1.0, 0.9, 0.1, 1000.0)
            * Matrix4::new_translation(&nalgebra::Vector3::new(0.0, 0.0, -40.0));
        let layers = prs3d_drawer_hidden_line_display(&cube, &cam, true).unwrap();
        assert!(!layers.visible.is_empty(), "a cube has visible edges");
        // A cube always has back edges occluded by the front faces.
        assert!(!layers.hidden.is_empty(), "a cube has hidden edges");
    }

    #[test]
    fn disabled_suppresses_the_hidden_layer() {
        let cube = valenx_cad::box_solid(10.0, 10.0, 10.0).unwrap();
        let cam = Matrix4::new_perspective(1.0, 0.9, 0.1, 1000.0)
            * Matrix4::new_translation(&nalgebra::Vector3::new(0.0, 0.0, -40.0));
        let layers = prs3d_drawer_hidden_line_display(&cube, &cam, false).unwrap();
        assert!(
            layers.hidden.is_empty(),
            "hidden layer must be empty when display is off"
        );
        // Visible edges are still produced regardless of the toggle.
        assert!(!layers.visible.is_empty());
    }
}
