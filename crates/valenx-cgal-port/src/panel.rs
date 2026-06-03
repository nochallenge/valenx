//! UI panel state.

use crate::voronoi::VoronoiCell;

/// Workbench-panel state.
#[derive(Default)]
pub struct CgalPanelState {
    /// Last input 2D point set.
    pub points_2d: Vec<[f64; 2]>,
    /// Last Delaunay triangulation.
    pub last_triangulation: Vec<(usize, usize, usize)>,
    /// Last convex hull (2D).
    pub last_convex_hull: Vec<[f64; 2]>,
    /// Last alpha-shape boundary.
    pub last_alpha_shape: Vec<[f64; 2]>,
    /// Last Voronoi cells.
    pub last_voronoi: Vec<VoronoiCell>,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
}

impl CgalPanelState {
    /// Empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record status.
    pub fn set_status(&mut self, msg: impl Into<String>) {
        self.last_status = Some(msg.into());
        self.last_error = None;
    }

    /// Record error.
    pub fn set_error(&mut self, msg: impl Into<String>) {
        self.last_error = Some(msg.into());
        self.last_status = None;
    }
}
