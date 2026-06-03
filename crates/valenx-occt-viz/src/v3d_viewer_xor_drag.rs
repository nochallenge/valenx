//! Phase 170 — `V3d_View::DrawXORDragRect()` — overlay rubber-band
//! rectangle for box-selection feedback.
//!
//! ## What OCCT does
//!
//! `V3d_View` exposes `DrawXORDragRect(x1, y1, x2, y2)` which draws a
//! 2D screen-space rectangle in XOR blending mode (so the next call to
//! the same coordinates erases it — the trick that pre-shader hardware
//! used to animate a drag-rectangle without re-rendering the scene).
//! Modern OCCT layers an `AIS_RubberBand` on top of the same primitive.
//!
//! ## v1 status
//!
//! **Honest v1.** Returns a normalized [`DragRect`] that the caller
//! (the viewport's `painter.rect_stroke` call in
//! `valenx_app::viewport`) can paint as a 1-pixel dashed-yellow
//! outline. egui's painter naturally blends overlays on top of the
//! already-rendered viewport image, so the XOR-erase trick OCCT used
//! is unnecessary — egui re-paints every frame anyway. This op only
//! validates and normalizes the coordinates (swaps min/max so the
//! rectangle has positive extent).

use crate::error::OcctVizError;

/// Screen-pixel rectangle with normalized min/max coordinates
/// (`min.x <= max.x`, `min.y <= max.y`).
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct DragRect {
    /// Top-left corner in screen pixels.
    pub min: [f32; 2],
    /// Bottom-right corner in screen pixels.
    pub max: [f32; 2],
}

impl DragRect {
    /// Width (`max.x - min.x`).
    pub fn width(&self) -> f32 {
        self.max[0] - self.min[0]
    }
    /// Height (`max.y - min.y`).
    pub fn height(&self) -> f32 {
        self.max[1] - self.min[1]
    }
}

/// Build a normalized [`DragRect`] from two arbitrary corners.
///
/// The two corners may be in any order — the result always has
/// `min <= max` on both axes. Zero-extent rectangles are valid
/// (corner-snap during a click-without-drag).
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if any coordinate is non-finite.
pub fn v3d_viewer_xor_drag(p1: [f32; 2], p2: [f32; 2]) -> Result<DragRect, OcctVizError> {
    for (i, v) in p1.iter().chain(p2.iter()).enumerate() {
        if !v.is_finite() {
            return Err(OcctVizError::bad_input(
                "corner",
                format!("component {i} is non-finite"),
            ));
        }
    }
    Ok(DragRect {
        min: [p1[0].min(p2[0]), p1[1].min(p2[1])],
        max: [p1[0].max(p2[0]), p1[1].max(p2[1])],
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalizes_reverse_order() {
        let r = v3d_viewer_xor_drag([5.0, 10.0], [0.0, 0.0]).unwrap();
        assert_eq!(r.min, [0.0, 0.0]);
        assert_eq!(r.max, [5.0, 10.0]);
        assert_eq!(r.width(), 5.0);
        assert_eq!(r.height(), 10.0);
    }

    #[test]
    fn rejects_nan() {
        let err = v3d_viewer_xor_drag([0.0, f32::NAN], [1.0, 1.0]).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn zero_extent_is_valid() {
        let r = v3d_viewer_xor_drag([5.0, 5.0], [5.0, 5.0]).unwrap();
        assert_eq!(r.width(), 0.0);
        assert_eq!(r.height(), 0.0);
    }
}
