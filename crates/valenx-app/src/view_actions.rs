//! Visual / camera-related methods on [`ValenxApp`]. Split out of
//! `lib.rs` as part of the structural refactor to keep each `impl
//! ValenxApp` block focused on a single concern.

use valenx_viz::OrbitCamera;

use crate::mesh_loader::mesh_bounding_box;
use crate::viewport::ShadingMode;
use crate::ValenxApp;

impl ValenxApp {
    /// Mutable borrow of the orbit camera — exposed so panel code can
    /// drive yaw/pitch/dolly without going through every action.
    pub fn camera_mut(&mut self) -> &mut OrbitCamera {
        &mut self.camera
    }

    /// Flip the viewport's shading mode between [`ShadingMode::Shaded`]
    /// and [`ShadingMode::Wireframe`].
    pub fn toggle_shading(&mut self) {
        self.shading = match self.shading {
            ShadingMode::Shaded => ShadingMode::Wireframe,
            ShadingMode::Wireframe => ShadingMode::Shaded,
        };
    }

    /// Collapse / expand the bottom Residuals / Log dock.
    ///
    /// When collapsed the panel renders only its header strip (the tab
    /// selectors + the named toggle button) and skips the content body,
    /// so it shrinks to a thin bar. Backs the header's AI-drivable
    /// "Collapse panel" / "Expand panel" button.
    pub fn toggle_bottom_panel(&mut self) {
        self.bottom_panel_collapsed = !self.bottom_panel_collapsed;
    }

    /// Collapse / expand the left Browser panel.
    ///
    /// When collapsed the panel shrinks to a thin vertical bar holding
    /// only the named "Expand panel" button; the heavy browser body
    /// (open-tabs list, navigator, Cases / Geometry / Mesh / Results) is
    /// skipped. Mirrors [`Self::toggle_bottom_panel`]. Backs the
    /// AI-drivable "Collapse panel" / "Expand panel" button.
    pub fn toggle_browser_panel(&mut self) {
        self.browser_collapsed = !self.browser_collapsed;
    }

    /// Reframe the camera around the loaded STL's bounding box.
    /// No-op when nothing is loaded.
    pub fn frame_current_stl(&mut self) {
        if let Some(stl) = &self.stl {
            if let Some((min, max)) = stl.mesh.bounding_box() {
                self.camera.frame_bounds(min, max);
            }
        }
    }

    /// Frame the camera around the loaded mesh's bounding box.
    pub fn frame_current_mesh(&mut self) {
        if let Some(loaded) = &self.mesh {
            if let Some((min, max)) = mesh_bounding_box(&loaded.mesh) {
                self.camera.frame_bounds(min, max);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::ValenxApp;

    #[test]
    fn toggle_shading_round_trips() {
        let mut app = ValenxApp::default();
        let original = app.shading;
        app.toggle_shading();
        assert_ne!(app.shading, original);
        app.toggle_shading();
        assert_eq!(app.shading, original);
    }

    #[test]
    fn toggle_bottom_panel_flips_collapsed_flag() {
        let mut app = ValenxApp::default();
        // Defaults to expanded (false) via `#[derive(Default)]`.
        assert!(!app.bottom_panel_collapsed);
        app.toggle_bottom_panel();
        assert!(app.bottom_panel_collapsed, "first toggle collapses");
        app.toggle_bottom_panel();
        assert!(!app.bottom_panel_collapsed, "second toggle re-expands");
    }

    #[test]
    fn toggle_browser_panel_flips_collapsed_flag() {
        let mut app = ValenxApp::default();
        // Defaults to expanded (false) via `#[derive(Default)]`.
        assert!(!app.browser_collapsed);
        app.toggle_browser_panel();
        assert!(app.browser_collapsed, "first toggle collapses");
        app.toggle_browser_panel();
        assert!(!app.browser_collapsed, "second toggle re-expands");
    }
}
