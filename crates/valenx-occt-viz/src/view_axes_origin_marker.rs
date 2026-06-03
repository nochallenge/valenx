//! Phase 199 — `Prs3d_Drawer::SetDatumAspect` — XYZ origin axes
//! always visible at world origin.
//!
//! ## What OCCT does
//!
//! `Prs3d_Drawer::SetDatumAspect(Prs3d_DatumAspect)` enables a small
//! colour-coded trihedron at the world origin (red=X arrow,
//! green=Y arrow, blue=Z arrow with text labels). The trihedron stays
//! at fixed *screen size* regardless of camera distance — OCCT
//! computes a per-frame scale factor to keep it visually constant
//! at ~50 px on screen.
//!
//! ## v1 status
//!
//! **Honest v1.** Returns the validated [`OriginMarkerConfig`] that
//! the caller (typically the View menu's "Show World Axes" toggle)
//! stores in app state. The actual rendering happens in the
//! viewport's egui painter pass — three line segments + three
//! `text_at` calls. The per-frame screen-constant scaling is the
//! caller's responsibility (use `camera.distance * pixels_per_unit`
//! at the trihedron's world position to compute the line length).

use crate::error::OcctVizError;

/// Configuration for the world-origin trihedron marker.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct OriginMarkerConfig {
    /// Whether the marker is visible.
    pub visible: bool,
    /// Per-axis arm length in pixels (caller multiplies by per-frame
    /// scale to get world units). Default 50.
    pub arm_pixels: f32,
    /// Whether to draw the "X" / "Y" / "Z" letter labels at the
    /// arm tips.
    pub draw_labels: bool,
}

impl Default for OriginMarkerConfig {
    fn default() -> Self {
        Self {
            visible: true,
            arm_pixels: 50.0,
            draw_labels: true,
        }
    }
}

/// Build a validated [`OriginMarkerConfig`].
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `arm_pixels` is non-finite or
///   outside `[5.0, 500.0]`.
pub fn view_axes_origin_marker(
    visible: bool,
    arm_pixels: f32,
    draw_labels: bool,
) -> Result<OriginMarkerConfig, OcctVizError> {
    if !arm_pixels.is_finite() {
        return Err(OcctVizError::bad_input("arm_pixels", "must be finite"));
    }
    if !(5.0..=500.0).contains(&arm_pixels) {
        return Err(OcctVizError::bad_input(
            "arm_pixels",
            format!("must be in [5, 500] (got {arm_pixels})"),
        ));
    }
    Ok(OriginMarkerConfig {
        visible,
        arm_pixels,
        draw_labels,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_visible_with_labels() {
        let d = OriginMarkerConfig::default();
        assert!(d.visible);
        assert!(d.draw_labels);
        assert_eq!(d.arm_pixels, 50.0);
    }

    #[test]
    fn rejects_too_small() {
        let err = view_axes_origin_marker(true, 1.0, true).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_too_large() {
        let err = view_axes_origin_marker(true, 1000.0, true).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn round_trips_valid_config() {
        let c = view_axes_origin_marker(false, 100.0, false).unwrap();
        assert!(!c.visible);
        assert!(!c.draw_labels);
        assert_eq!(c.arm_pixels, 100.0);
    }
}
