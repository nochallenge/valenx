//! Phase 183 — `Prs3d_Drawer::SetWireAspect` per-edge colour
//! override.
//!
//! ## What OCCT does
//!
//! `Prs3d_Drawer::SetWireAspect(Graphic3d_AspectLine3d)` overrides the
//! line colour for one or more `TopoDS_Edge` owners. Common use:
//! highlight "fixed" edges in a sketch as black, "computed" edges as
//! blue, "construction" edges as dashed grey. Internally OCCT keys the
//! override on edge index.
//!
//! ## v1 status
//!
//! **Honest v1.** Same validation pattern as
//! [`crate::prs3d_drawer_face_color()`]: builds an
//! [`EdgeColorRgba`] from validated RGBA. Storage is the caller's
//! responsibility; the wgpu wireframe pass per-edge upload arrives in
//! Phase 200.5 alongside face colouring and picking.

use crate::error::OcctVizError;

/// RGBA colour bundled for edge overrides. Channels are in linear RGB
/// (0..=1).
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct EdgeColorRgba {
    /// Red 0..=1.
    pub r: f32,
    /// Green 0..=1.
    pub g: f32,
    /// Blue 0..=1.
    pub b: f32,
    /// Alpha 0..=1.
    pub a: f32,
}

/// Build an [`EdgeColorRgba`] from RGBA components.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] when any channel is out of `[0, 1]`
///   or non-finite.
pub fn prs3d_drawer_edge_color(
    r: f32,
    g: f32,
    b: f32,
    a: f32,
) -> Result<EdgeColorRgba, OcctVizError> {
    for (name, v) in [("rgba.r", r), ("rgba.g", g), ("rgba.b", b), ("rgba.a", a)] {
        if !v.is_finite() {
            return Err(OcctVizError::bad_input(name, "must be finite"));
        }
        if !(0.0..=1.0).contains(&v) {
            return Err(OcctVizError::bad_input(
                name,
                format!("must be in [0, 1] (got {v})"),
            ));
        }
    }
    Ok(EdgeColorRgba { r, g, b, a })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_out_of_range_blue() {
        let err = prs3d_drawer_edge_color(0.5, 0.5, 2.0, 1.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn round_trips() {
        let c = prs3d_drawer_edge_color(0.1, 0.2, 0.3, 0.4).unwrap();
        assert_eq!(c.r, 0.1);
        assert_eq!(c.a, 0.4);
    }
}
