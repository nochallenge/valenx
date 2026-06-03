//! Orbit camera — the turntable model the ViewCube drives.
//!
//! Designed for the common CAD pattern: hold middle-mouse to orbit,
//! scroll to zoom, shift+middle to pan. The camera stores a
//! `target` point (center of orbit), `distance` from target, and two
//! angles (`azimuth_deg`, `elevation_deg`). The eye position is
//! derived; this lets the ViewCube snap to canonical views without
//! losing context about what's being looked at.
//!
//! Math is in a right-handed world (X right, Y up, Z towards viewer),
//! consistent with the Fusion 360 convention Valenx mirrors.

use nalgebra::{Matrix4, Point3, Vector3};

/// One of the canonical camera orientations exposed by the ViewCube.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ViewDirection {
    Front,
    Back,
    Top,
    Bottom,
    Left,
    Right,
    Iso,
}

impl ViewDirection {
    /// `(azimuth_deg, elevation_deg)` for this canonical view.
    pub fn angles(self) -> (f32, f32) {
        match self {
            ViewDirection::Front => (0.0, 0.0),
            ViewDirection::Back => (180.0, 0.0),
            ViewDirection::Right => (90.0, 0.0),
            ViewDirection::Left => (-90.0, 0.0),
            ViewDirection::Top => (0.0, 90.0),
            ViewDirection::Bottom => (0.0, -90.0),
            ViewDirection::Iso => (45.0, 35.264), // arctan(1/sqrt(2))
        }
    }
}

/// Projection mode for the camera — perspective (foreshortening) or
/// orthographic (parallel). OCCT's `V3d_View::SetType` toggles the
/// same distinction; CAD users switch to orthographic for
/// dimension-true drafting views.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum ProjectionMode {
    /// Standard perspective frustum — distant geometry shrinks.
    #[default]
    Perspective,
    /// Parallel / orthographic frustum — no foreshortening.
    Orthographic,
}

/// The orbit camera state.
#[derive(Clone, Debug)]
pub struct OrbitCamera {
    pub target: Point3<f32>,
    pub distance: f32,
    pub azimuth_deg: f32,
    pub elevation_deg: f32,
    pub fov_y_deg: f32,
    pub near: f32,
    pub far: f32,
    /// Projection mode used by [`OrbitCamera::projection_matrix`].
    pub projection_mode: ProjectionMode,
}

impl Default for OrbitCamera {
    fn default() -> Self {
        Self {
            target: Point3::origin(),
            distance: 13.0,
            azimuth_deg: 45.0,
            elevation_deg: 25.0,
            fov_y_deg: 35.0,
            near: 0.01,
            far: 10_000.0,
            projection_mode: ProjectionMode::Perspective,
        }
    }
}

impl OrbitCamera {
    /// Current eye position, derived from target + angles + distance.
    pub fn eye(&self) -> Point3<f32> {
        let az = self.azimuth_deg.to_radians();
        let el = self.elevation_deg.to_radians();
        let r = self.distance;
        // Y-up spherical to cartesian:
        //   x = r cos(el) sin(az)
        //   y = r sin(el)
        //   z = r cos(el) cos(az)
        let x = r * el.cos() * az.sin();
        let y = r * el.sin();
        let z = r * el.cos() * az.cos();
        self.target + Vector3::new(x, y, z)
    }

    /// Zoom by a fractional amount (0.1 = "zoom in by 10%"). Positive
    /// zooms in (reduces distance), negative zooms out. Clamped so
    /// the camera can't invert through the target.
    pub fn zoom(&mut self, frac: f32) {
        self.distance = (self.distance * (1.0 - frac)).max(1e-4);
    }

    /// Orbit by a screen-pixel delta. Horizontal movement changes
    /// azimuth, vertical changes elevation (clamped).
    pub fn orbit(&mut self, dx_deg: f32, dy_deg: f32) {
        self.azimuth_deg = (self.azimuth_deg + dx_deg) % 360.0;
        self.elevation_deg = (self.elevation_deg + dy_deg).clamp(-89.9, 89.9);
    }

    /// Snap to a canonical ViewCube direction.
    pub fn set_view(&mut self, dir: ViewDirection) {
        let (az, el) = dir.angles();
        self.azimuth_deg = az;
        self.elevation_deg = el;
    }

    /// Frame the camera around an axis-aligned bounding box so the
    /// whole thing fits in view. Useful after loading an STL.
    pub fn frame_bounds(&mut self, min: [f32; 3], max: [f32; 3]) {
        let center = [
            (min[0] + max[0]) * 0.5,
            (min[1] + max[1]) * 0.5,
            (min[2] + max[2]) * 0.5,
        ];
        let size = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
        let diag = (size[0] * size[0] + size[1] * size[1] + size[2] * size[2]).sqrt();
        self.target = Point3::new(center[0], center[1], center[2]);
        // Geometry-by-FoV distance with a small margin.
        let half = (self.fov_y_deg * 0.5).to_radians().tan();
        self.distance = (diag * 0.5 / half).max(1e-3) * 1.15;
    }

    /// Right-handed view matrix (world → view).
    pub fn view_matrix(&self) -> Matrix4<f32> {
        Matrix4::look_at_rh(&self.eye(), &self.target, &Vector3::y())
    }

    /// Projection matrix (view → clip) for the given aspect ratio
    /// (width / height).
    ///
    /// Branches on [`OrbitCamera::projection_mode`]:
    ///
    /// - [`ProjectionMode::Perspective`] — a standard perspective
    ///   frustum from `fov_y_deg` + `near` / `far`.
    /// - [`ProjectionMode::Orthographic`] — a parallel frustum whose
    ///   half-height equals `distance * tan(fov_y / 2)`. That is the
    ///   visible footprint perspective showed *at the target plane*,
    ///   so toggling projection keeps the model the same on-screen
    ///   size (the "swap projection without zooming" UX).
    pub fn projection_matrix(&self, aspect: f32) -> Matrix4<f32> {
        let aspect = aspect.max(1e-6);
        match self.projection_mode {
            ProjectionMode::Perspective => Matrix4::new_perspective(
                aspect,
                self.fov_y_deg.to_radians(),
                self.near,
                self.far,
            ),
            ProjectionMode::Orthographic => {
                // Half-height of the perspective footprint at the
                // target plane (distance away from the eye).
                let half_h = (self.distance * (self.fov_y_deg * 0.5).to_radians().tan())
                    .max(1e-6);
                let half_w = half_h * aspect;
                Matrix4::new_orthographic(
                    -half_w,
                    half_w,
                    -half_h,
                    half_h,
                    self.near,
                    self.far,
                )
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_angles_round_trip() {
        let mut cam = OrbitCamera::default();
        cam.set_view(ViewDirection::Front);
        assert_eq!(cam.azimuth_deg, 0.0);
        assert_eq!(cam.elevation_deg, 0.0);
        cam.set_view(ViewDirection::Top);
        assert_eq!(cam.elevation_deg, 90.0);
    }

    #[test]
    fn orbit_clamps_elevation() {
        let mut cam = OrbitCamera::default();
        cam.orbit(0.0, 500.0);
        assert!(cam.elevation_deg <= 89.9);
        cam.orbit(0.0, -1000.0);
        assert!(cam.elevation_deg >= -89.9);
    }

    #[test]
    fn zoom_never_inverts() {
        let mut cam = OrbitCamera::default();
        cam.zoom(2.0); // would take distance negative without clamp
        assert!(cam.distance > 0.0);
    }

    #[test]
    fn frame_bounds_centers_target() {
        let mut cam = OrbitCamera::default();
        cam.frame_bounds([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]);
        assert!((cam.target.x - 5.0).abs() < 1e-5);
        assert!((cam.target.y - 5.0).abs() < 1e-5);
        assert!((cam.target.z - 5.0).abs() < 1e-5);
        assert!(cam.distance > 0.0);
    }

    #[test]
    fn view_matrix_is_finite() {
        let cam = OrbitCamera::default();
        let m = cam.view_matrix();
        assert!(m.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn default_projection_mode_is_perspective() {
        assert_eq!(OrbitCamera::default().projection_mode, ProjectionMode::Perspective);
    }

    #[test]
    fn orthographic_projection_matrix_is_finite_and_distinct() {
        let persp = OrbitCamera {
            projection_mode: ProjectionMode::Perspective,
            ..Default::default()
        };
        let ortho = OrbitCamera {
            projection_mode: ProjectionMode::Orthographic,
            ..Default::default()
        };

        let mp = persp.projection_matrix(1.5);
        let mo = ortho.projection_matrix(1.5);
        assert!(mp.iter().all(|v| v.is_finite()));
        assert!(mo.iter().all(|v| v.is_finite()));
        // An orthographic matrix has its bottom-right element == 1 and
        // m[(3,2)] == 0; a perspective matrix has m[(3,3)] == 0 and
        // m[(3,2)] == -1. This distinguishes the two projections.
        assert!((mo[(3, 3)] - 1.0).abs() < 1e-6, "ortho m[3,3] should be 1");
        assert!(mo[(3, 2)].abs() < 1e-6, "ortho m[3,2] should be 0");
        assert!(mp[(3, 3)].abs() < 1e-6, "perspective m[3,3] should be 0");
    }
}
