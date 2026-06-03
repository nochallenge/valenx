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
}
