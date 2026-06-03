//! Phase 182 — `Prs3d_Drawer::SetShadingAspect` per-face colour
//! override.
//!
//! ## What OCCT does
//!
//! `Prs3d_Drawer::SetShadingAspect(Graphic3d_AspectFillArea3d)`
//! overrides the diffuse fill colour for one or more `TopoDS_Face`
//! owners independent of the parent shape's material preset. Used for
//! "highlight this face yellow" / "colour-code by region" workflows.
//! Internally OCCT stores the override in a `TColStd_DataMapOfInteger3d`
//! keyed by face index.
//!
//! ## v1 status
//!
//! **Honest v1.** Validates the RGBA values (each in [0, 1], finite)
//! and returns the validated override packed into a [`FaceColorRgba`].
//! Storage is the caller's responsibility — typically the Mesh Toolbox
//! holds a `HashMap<(parent_id, face_index), FaceColorRgba>` that
//! `valenx_app::viewport` consults during shading. The actual
//! per-face triangle-colour upload to the wgpu vertex buffer is Phase
//! 200.5 (paired with the picking-pass face tagging since the same
//! tri-range table is needed for both).

use crate::error::OcctVizError;

/// RGBA colour bundled for face overrides. Channels are in linear RGB
/// (0..=1).
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct FaceColorRgba {
    /// Red 0..=1.
    pub r: f32,
    /// Green 0..=1.
    pub g: f32,
    /// Blue 0..=1.
    pub b: f32,
    /// Alpha 0..=1 (1 = opaque).
    pub a: f32,
}

/// Build a [`FaceColorRgba`] from RGBA components, validating each
/// channel is finite and in `[0, 1]`.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] when any channel is out of range or
///   non-finite.
pub fn prs3d_drawer_face_color(
    r: f32,
    g: f32,
    b: f32,
    a: f32,
) -> Result<FaceColorRgba, OcctVizError> {
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
    Ok(FaceColorRgba { r, g, b, a })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_negative() {
        let err = prs3d_drawer_face_color(-0.1, 0.5, 0.5, 1.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_above_one() {
        let err = prs3d_drawer_face_color(0.5, 1.5, 0.5, 1.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_nan_alpha() {
        let err = prs3d_drawer_face_color(0.5, 0.5, 0.5, f32::NAN).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn round_trips_valid_color() {
        let c = prs3d_drawer_face_color(0.2, 0.4, 0.6, 0.8).unwrap();
        assert_eq!(c.r, 0.2);
        assert_eq!(c.g, 0.4);
        assert_eq!(c.b, 0.6);
        assert_eq!(c.a, 0.8);
    }

    #[test]
    fn accepts_boundary_values() {
        let zero = prs3d_drawer_face_color(0.0, 0.0, 0.0, 0.0).unwrap();
        let one = prs3d_drawer_face_color(1.0, 1.0, 1.0, 1.0).unwrap();
        assert_eq!(zero.r, 0.0);
        assert_eq!(one.r, 1.0);
    }
}
