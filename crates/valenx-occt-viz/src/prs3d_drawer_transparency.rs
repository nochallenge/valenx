//! Phase 184 — `Prs3d_Drawer::SetTransparency` per-object
//! transparency.
//!
//! ## What OCCT does
//!
//! `Prs3d_Drawer::SetTransparency(Standard_Real)` overrides the alpha
//! channel for the entire object (0 = opaque, 1 = invisible). Distinct
//! from [`crate::prs3d_drawer_face_color()`] which sets *per-face* RGB+A
//! — this one is a single scalar applied to all faces in the object.
//! Transparent objects are sorted back-to-front before draw to make
//! standard alpha blending look right.
//!
//! ## v1 status
//!
//! **Honest v1.** Validates `transparency ∈ [0, 1]` and returns the
//! validated value. The renderer reads this once per object per
//! frame; the back-to-front sort is already implemented in
//! `valenx_app::viewport`'s painter's-algorithm depth sort for the
//! egui-paint path. The wgpu path will gain alpha sorting in Phase
//! 188.5 (deferred until the wgpu pipeline gains a transparency
//! flag — currently every solid is opaque on the wgpu path).

use crate::error::OcctVizError;

/// Set the transparency for an object.
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `transparency` is not finite or
///   outside `[0, 1]`.
pub fn prs3d_drawer_transparency(transparency: f32) -> Result<f32, OcctVizError> {
    if !transparency.is_finite() {
        return Err(OcctVizError::bad_input("transparency", "must be finite"));
    }
    if !(0.0..=1.0).contains(&transparency) {
        return Err(OcctVizError::bad_input(
            "transparency",
            format!("must be in [0, 1] (got {transparency})"),
        ));
    }
    Ok(transparency)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_negative() {
        let err = prs3d_drawer_transparency(-0.1).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_above_one() {
        let err = prs3d_drawer_transparency(1.5).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_nan() {
        let err = prs3d_drawer_transparency(f32::NAN).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn accepts_boundary_values() {
        assert_eq!(prs3d_drawer_transparency(0.0).unwrap(), 0.0);
        assert_eq!(prs3d_drawer_transparency(1.0).unwrap(), 1.0);
        assert_eq!(prs3d_drawer_transparency(0.5).unwrap(), 0.5);
    }
}
