//! Phase 200 — `Aspect_GridDrawer` rectangular XY ground-plane grid.
//!
//! **This is the final phase of the 200-phase FreeCAD-parity
//! roadmap.** Together with the other 39 modules in this crate +
//! the 160 phases that preceded it across Rounds 1-5, it
//! completes Valenx's intended OCCT/FreeCAD feature parity surface.
//!
//! ## What OCCT does
//!
//! `V3d_Viewer::ActivateGrid(Aspect_GridType_Rectangular,
//! Aspect_GridDrawMode_Lines)` draws a 2D grid in the world XY plane
//! (Z=0) as a guide for spatial reference. Configurable per-cell
//! spacing, line colour, and either lines-only or lines+points
//! rendering. Some viewports use the grid as a sketch-snap target;
//! Valenx's sketcher does that separately (`valenx_sketch`) so the
//! grid here is purely visual.
//!
//! ## v1 status
//!
//! **Honest v1.** Returns a validated [`GridConfig`] the caller
//! stores in app state. The viewport's egui painter draws the grid as
//! a series of horizontal + vertical line segments projected from
//! world space via the current camera. The actual line drawing is
//! caller responsibility; this op just validates the cell spacing
//! and extent.

use crate::error::OcctVizError;

/// Configuration for the XY ground-plane grid.
#[derive(Copy, Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GridConfig {
    /// Whether the grid is visible.
    pub visible: bool,
    /// World-units between adjacent grid lines. Must be > 0.
    pub cell_spacing: f32,
    /// Grid extends from `-half_extent` to `+half_extent` on both X
    /// and Y. Must be > 0.
    pub half_extent: f32,
    /// Whether to draw bold major lines every `major_every` cells
    /// (0 = no major lines).
    pub major_every: u32,
}

impl Default for GridConfig {
    fn default() -> Self {
        Self {
            visible: true,
            cell_spacing: 10.0,
            half_extent: 100.0,
            major_every: 10,
        }
    }
}

/// Build a validated [`GridConfig`].
///
/// # Errors
///
/// - [`OcctVizError::BadInput`] if `cell_spacing` or `half_extent`
///   is ≤ 0 or non-finite, or `half_extent < cell_spacing` (not
///   even one cell fits).
pub fn view_grid_floor_xy(
    visible: bool,
    cell_spacing: f32,
    half_extent: f32,
    major_every: u32,
) -> Result<GridConfig, OcctVizError> {
    if !cell_spacing.is_finite() || cell_spacing <= 0.0 {
        return Err(OcctVizError::bad_input(
            "cell_spacing",
            format!("must be > 0 (got {cell_spacing})"),
        ));
    }
    if !half_extent.is_finite() || half_extent <= 0.0 {
        return Err(OcctVizError::bad_input(
            "half_extent",
            format!("must be > 0 (got {half_extent})"),
        ));
    }
    if half_extent < cell_spacing {
        return Err(OcctVizError::bad_input(
            "half_extent",
            format!("must be >= cell_spacing ({cell_spacing}); got {half_extent}"),
        ));
    }
    Ok(GridConfig {
        visible,
        cell_spacing,
        half_extent,
        major_every,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_sensible() {
        let g = GridConfig::default();
        assert!(g.visible);
        assert!(g.cell_spacing > 0.0);
        assert!(g.half_extent > g.cell_spacing);
        assert_eq!(g.major_every, 10);
    }

    #[test]
    fn rejects_zero_spacing() {
        let err = view_grid_floor_xy(true, 0.0, 100.0, 10).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_negative_extent() {
        let err = view_grid_floor_xy(true, 10.0, -100.0, 10).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_extent_smaller_than_spacing() {
        let err = view_grid_floor_xy(true, 100.0, 50.0, 10).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn rejects_nan_spacing() {
        let err = view_grid_floor_xy(true, f32::NAN, 100.0, 10).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn round_trips_valid_config() {
        let g = view_grid_floor_xy(true, 5.0, 50.0, 5).unwrap();
        assert!(g.visible);
        assert_eq!(g.cell_spacing, 5.0);
        assert_eq!(g.half_extent, 50.0);
        assert_eq!(g.major_every, 5);
    }

    #[test]
    fn major_every_zero_is_valid() {
        // 0 = no major lines, all cells equal.
        let g = view_grid_floor_xy(true, 10.0, 100.0, 0).unwrap();
        assert_eq!(g.major_every, 0);
    }
}
