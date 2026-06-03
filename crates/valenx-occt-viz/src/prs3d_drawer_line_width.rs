//! Phase 185 — `Prs3d_Drawer::SetLineAspect` width — per-object edge
//! line width in pixels.
//!
//! ## What OCCT does
//!
//! `Prs3d_Drawer::SetLineAspect(Graphic3d_AspectLine3d::SetWidth)`
//! sets the per-object edge stroke width in pixels (capped at 1..=10
//! on most hardware via `GL_MAX_LINE_WIDTH`). Used for "make the
//! selected edge 3px so it stands out" workflows. Real OCCT applies
//! the same width to all subshapes — for per-edge widths you need
//! [`crate::prs3d_drawer_line_style()`] paired with
//! [`crate::prs3d_drawer_edge_color()`] (the {colour, style, width}
//! triplet is the unit of override).
//!
//! ## v1 status
//!
//! **Honest v1.** Validates `0.5 ≤ width ≤ 10.0` (range matches
//! `wgpu::Limits::max_line_width` floor across vendors) and returns
//! the value. egui's `Stroke::width` reads f32 directly so the egui-
//! paint path is wired; the wgpu wireframe pass uses `wgpu::PolygonMode::Line`
//! with the device's configured default width (per-object override
//! requires changing `wgpu::RenderPipelineDescriptor::primitive.polygon_mode`
//! which the wgpu pipeline doesn't yet expose dynamically; Phase 188.5
//! adds the per-object pipeline variant).

use crate::error::OcctVizError;

/// Minimum permitted line width (matches OCCT's
/// `Graphic3d_AspectLine3d::SetWidth` floor).
pub const MIN_LINE_WIDTH: f32 = 0.5;
/// Maximum permitted line width (matches the most-restrictive
/// `GL_MAX_LINE_WIDTH` across Valenx's supported hardware tier).
pub const MAX_LINE_WIDTH: f32 = 10.0;

/// Validate and return the requested line width.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `width` is not finite or outside
///   `[MIN_LINE_WIDTH, MAX_LINE_WIDTH]`.
pub fn prs3d_drawer_line_width(width: f32) -> Result<f32, OcctVizError> {
    if !width.is_finite() {
        return Err(OcctVizError::bad_input("width", "must be finite"));
    }
    if !(MIN_LINE_WIDTH..=MAX_LINE_WIDTH).contains(&width) {
        return Err(OcctVizError::bad_input(
            "width",
            format!("must be in [{MIN_LINE_WIDTH}, {MAX_LINE_WIDTH}] (got {width})"),
        ));
    }
    Ok(width)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_below_floor() {
        let err = prs3d_drawer_line_width(0.1).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_above_cap() {
        let err = prs3d_drawer_line_width(20.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn accepts_typical_widths() {
        assert_eq!(prs3d_drawer_line_width(1.0).unwrap(), 1.0);
        assert_eq!(prs3d_drawer_line_width(2.0).unwrap(), 2.0);
        assert_eq!(prs3d_drawer_line_width(5.0).unwrap(), 5.0);
    }
}
