//! Detected-feature point type.

/// A detected feature point in image coordinates.
///
/// Coordinates use the same convention as [`crate::GrayImage`]: `(0, 0)`
/// is the top-left corner, `x` to the right, `y` downward. FAST keypoints
/// land on integer pixel centres, but the fields are `f32` so a later
/// sub-pixel refinement stage can store fractional positions without a
/// type change.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Keypoint {
    /// Sub-pixel column coordinate (pixels from the left edge).
    pub x: f32,
    /// Sub-pixel row coordinate (pixels from the top edge).
    pub y: f32,
    /// Detector response strength used for ranking / non-maximum
    /// suppression. For FAST this is the corner score `V` (the largest
    /// brightness threshold at which the contiguous-arc test still
    /// passes); larger is a stronger corner.
    pub score: f32,
    /// Dominant orientation in radians, measured from the +x axis and
    /// growing toward +y, in `(-π, π]`. Computed from the intensity
    /// centroid of a circular patch; `0.0` before orientation assignment.
    pub angle: f32,
}

impl Keypoint {
    /// Construct a keypoint with an unset (`0.0`) orientation. The angle
    /// is filled in later by the descriptor stage.
    #[inline]
    #[must_use]
    pub fn new(x: f32, y: f32, score: f32) -> Self {
        Self {
            x,
            y,
            score,
            angle: 0.0,
        }
    }
}
