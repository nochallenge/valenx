//! 3D → 2D screen projection helper.
//!
//! Used by the interim wireframe viewport that draws through egui's
//! 2D painter before the `wgpu` render pass lands. The math is
//! standard: multiply a world-space point by view × projection,
//! perform the perspective divide, remap from NDC (-1..1) to screen
//! pixel space. Triangles whose vertices all fall behind the near
//! plane are culled.
//!
//! Keeping this module self-contained (it only depends on `nalgebra`
//! and the already-public `OrbitCamera`) means the viewport that
//! uses it doesn't need to know anything about the rasteriser behind
//! it — once `wgpu` lands, this module keeps serving the tests and
//! fallback rendering paths.

use nalgebra::{Matrix4, Vector4};

use crate::camera::OrbitCamera;

/// A 2D point in the viewport's pixel coordinate system. `(0, 0)` is
/// the top-left of the viewport rect.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScreenPoint {
    pub x: f32,
    pub y: f32,
    /// Depth in normalized device coordinates (0 = near, 1 = far).
    /// Useful for depth-sorted painters.
    pub depth: f32,
}

/// Result of projecting a single world-space point. `None` means the
/// point fell behind the near plane and should be clipped.
pub fn project_point(
    camera: &OrbitCamera,
    viewport_width: f32,
    viewport_height: f32,
    world: [f32; 3],
) -> Option<ScreenPoint> {
    if viewport_width <= 0.0 || viewport_height <= 0.0 {
        return None;
    }
    let aspect = viewport_width / viewport_height;
    let mvp: Matrix4<f32> = camera.projection_matrix(aspect) * camera.view_matrix();
    let clip = mvp * Vector4::new(world[0], world[1], world[2], 1.0);
    if clip.w <= 1e-6 {
        // Behind the camera.
        return None;
    }
    // Perspective divide → NDC in [-1, 1]^3.
    let ndc_x = clip.x / clip.w;
    let ndc_y = clip.y / clip.w;
    let ndc_z = clip.z / clip.w;
    // NDC → screen. Y-flip because egui's y grows downward.
    let sx = (ndc_x * 0.5 + 0.5) * viewport_width;
    let sy = (1.0 - (ndc_y * 0.5 + 0.5)) * viewport_height;
    Some(ScreenPoint {
        x: sx,
        y: sy,
        depth: ndc_z * 0.5 + 0.5,
    })
}

/// Project all three vertices of a triangle. Returns `None` if any
/// vertex is behind the near plane — a coarse cull that's good
/// enough for a wireframe preview; a full implementation clips
/// against the near plane.
pub fn project_triangle(
    camera: &OrbitCamera,
    viewport_width: f32,
    viewport_height: f32,
    vertices: &[[f32; 3]; 3],
) -> Option<[ScreenPoint; 3]> {
    let a = project_point(camera, viewport_width, viewport_height, vertices[0])?;
    let b = project_point(camera, viewport_width, viewport_height, vertices[1])?;
    let c = project_point(camera, viewport_width, viewport_height, vertices[2])?;
    Some([a, b, c])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::camera::ViewDirection;

    #[test]
    fn origin_maps_near_center_from_front() {
        let mut cam = OrbitCamera::default();
        cam.set_view(ViewDirection::Front);
        cam.target = nalgebra::Point3::origin();
        cam.distance = 10.0;
        let p = project_point(&cam, 800.0, 600.0, [0.0, 0.0, 0.0]).expect("projects");
        // Origin should map to roughly the middle of the viewport.
        assert!((p.x - 400.0).abs() < 2.0, "got {}", p.x);
        assert!((p.y - 300.0).abs() < 2.0, "got {}", p.y);
        assert!(p.depth > 0.0 && p.depth < 1.0);
    }

    #[test]
    fn behind_camera_culls() {
        // Front view looks down +Z towards origin from z = +10.
        // A point at z = +20 is *behind* the camera (camera eye is
        // at +10 looking at origin, so the far side is z < 0).
        // Ensure we get None for points outside the view frustum's
        // near side.
        let mut cam = OrbitCamera::default();
        cam.set_view(ViewDirection::Front);
        cam.target = nalgebra::Point3::origin();
        cam.distance = 10.0;
        // Place a point at z = +30 — which is behind the camera eye
        // at (0, 0, +10) looking at origin (looks along -Z).
        let culled = project_point(&cam, 800.0, 600.0, [0.0, 0.0, 30.0]);
        assert!(culled.is_none(), "expected cull, got {culled:?}");
    }

    #[test]
    fn zero_viewport_returns_none() {
        let cam = OrbitCamera::default();
        assert!(project_point(&cam, 0.0, 600.0, [0.0, 0.0, 0.0]).is_none());
        assert!(project_point(&cam, 800.0, 0.0, [0.0, 0.0, 0.0]).is_none());
    }
}
