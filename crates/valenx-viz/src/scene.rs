//! Pure helpers for the viewport scene furniture (grid, axes, gizmo).
//!
//! No egui / wgpu here — just math, so it's unit-testable and shared
//! between the GPU grid (adaptive spacing) and the egui HUD overlay
//! (ray-pick for cursor coordinates, gizmo→view mapping).

use nalgebra::{Matrix4, Point3, Vector3, Vector4};

use crate::{OrbitCamera, ViewDirection};

/// Pick a "nice" minor grid spacing (1 / 2 / 5 × 10ⁿ) for a camera at
/// `distance`, targeting roughly ten minor cells across the view.
pub fn nice_grid_spacing(distance: f32) -> f32 {
    let raw = (distance.max(1e-4) / 10.0) as f64;
    let base = 10f64.powf(raw.log10().floor());
    let m = raw / base;
    let nice = if m < 1.5 {
        1.0
    } else if m < 3.5 {
        2.0
    } else if m < 7.5 {
        5.0
    } else {
        10.0
    };
    (nice * base) as f32
}

/// Adaptive LOD parameters for the GPU ground grid.
///
/// Returns `(minor_a, blend_t)`:
/// - `minor_a`: the current minor grid spacing (same as `nice_grid_spacing(distance)`)
/// - `blend_t`: 0.0 → minor lines fully opaque; 1.0 → minor lines fully faded
///
/// As the camera zooms out within each "nice" band, `blend_t` rises from 0 to 1,
/// smoothly fading the minor lines before the spacing snaps to the next coarser
/// band. The shader implements a three-level crossfade:
/// - minor_a lines fade out (α = (1-t) × 0.45)
/// - major_a = minor_a × 10 transitions from major brightness (0.85) down to minor
///   brightness (0.45) — becoming the "new minor" of the next level
/// - major_b = minor_a × 100 fades in (α = t × 0.85) — becoming the new major
///
/// This gives continuous, pop-free transitions at every zoom level.
pub fn grid_lod_params(distance: f32) -> (f32, f32) {
    let minor_a = nice_grid_spacing(distance);
    let raw = distance.max(1e-4) as f64 / 10.0;
    let base = 10f64.powf(raw.log10().floor());
    let m = raw / base;
    let blend_t: f64 = if m < 1.5 {
        (m - 1.0) / 0.5
    } else if m < 3.5 {
        (m - 1.5) / 2.0
    } else if m < 7.5 {
        (m - 3.5) / 4.0
    } else {
        (m - 7.5) / 2.5
    };
    (minor_a, blend_t.clamp(0.0, 1.0) as f32)
}

/// A world-space ray.
#[derive(Copy, Clone, Debug)]
pub struct Ray {
    pub origin: Point3<f32>,
    pub dir: Vector3<f32>,
}

/// Unproject a screen point (pixels, origin top-left) into a world-space
/// ray using the camera's inverse view-projection.
pub fn ray_from_screen(cam: &OrbitCamera, w: f32, h: f32, screen: [f32; 2]) -> Ray {
    let aspect = (w / h.max(1.0)).max(1e-6);
    let vp = cam.projection_matrix(aspect) * cam.view_matrix();
    let inv = vp.try_inverse().unwrap_or_else(Matrix4::identity);
    let ndc_x = 2.0 * screen[0] / w.max(1.0) - 1.0;
    let ndc_y = 1.0 - 2.0 * screen[1] / h.max(1.0);
    let unproj = |z: f32| {
        let p = inv * Vector4::new(ndc_x, ndc_y, z, 1.0);
        let iw = if p.w.abs() < 1e-9 { 1.0 } else { p.w };
        Point3::new(p.x / iw, p.y / iw, p.z / iw)
    };
    let near = unproj(0.0);
    let far = unproj(1.0);
    Ray {
        origin: near,
        dir: (far - near).normalize(),
    }
}

/// Intersect a ray with the `y = 0` ground plane. `None` if the ray is
/// parallel to the plane or points away from it.
pub fn intersect_ground_y0(r: &Ray) -> Option<Point3<f32>> {
    if r.dir.y.abs() < 1e-6 {
        return None;
    }
    let t = -r.origin.y / r.dir.y;
    if t < 0.0 {
        return None;
    }
    Some(r.origin + r.dir * t)
}

/// Snap a ground-plane point to the nearest grid intersection at
/// `spacing` — X and Z independently, Y preserved. Returns the point
/// unchanged when `spacing <= 0`. This is the 3D-viewport analogue of
/// the Draft workbench's 2D `grid_snap`: it locks a free cursor onto
/// the same lattice the GPU grid draws, so a click lands exactly on a
/// grid node (Fusion-style snap).
pub fn snap_ground_point(p: Point3<f32>, spacing: f32) -> Point3<f32> {
    if spacing <= 0.0 {
        return p;
    }
    Point3::new(
        (p.x / spacing).round() * spacing,
        p.y,
        (p.z / spacing).round() * spacing,
    )
}

/// The six clickable faces of the corner orientation gizmo.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum GizmoFace {
    Front,
    Back,
    Left,
    Right,
    Top,
    Bottom,
}

/// Map a clicked gizmo face to the canonical camera view it snaps to.
pub fn gizmo_view_for_face(f: GizmoFace) -> ViewDirection {
    match f {
        GizmoFace::Front => ViewDirection::Front,
        GizmoFace::Back => ViewDirection::Back,
        GizmoFace::Left => ViewDirection::Left,
        GizmoFace::Right => ViewDirection::Right,
        GizmoFace::Top => ViewDirection::Top,
        GizmoFace::Bottom => ViewDirection::Bottom,
    }
}

/// Screen-space 2D direction `[x_right, y_down]` of each world axis
/// (X, Y, Z order) under `cam`, paired with the view-space depth (z;
/// `> 0` means the axis tip points toward the viewer). Drives the corner
/// orientation gizmo; orthographic, so perspective is ignored.
pub fn gizmo_axis_screen_dirs(cam: &OrbitCamera) -> [([f32; 2], f32); 3] {
    let view = cam.view_matrix();
    let axes = [
        Vector3::new(1.0_f32, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ];
    let mut out = [([0.0_f32, 0.0], 0.0_f32); 3];
    for (i, a) in axes.iter().enumerate() {
        let v = view.transform_vector(a);
        out[i] = ([v.x, -v.y], v.z);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Point3, Vector3};

    #[test]
    fn spacing_is_nice_and_scales_with_distance() {
        assert_eq!(nice_grid_spacing(10.0), 1.0);
        assert_eq!(nice_grid_spacing(100.0), 10.0);
        assert!((nice_grid_spacing(0.5) - 0.05).abs() < 1e-6);
        assert!(nice_grid_spacing(1000.0) >= nice_grid_spacing(100.0));
    }

    #[test]
    fn grid_lod_minor_matches_nice_spacing() {
        for d in [0.5_f32, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0] {
            let (minor, _) = grid_lod_params(d);
            assert_eq!(
                minor,
                nice_grid_spacing(d),
                "minor mismatch at distance={d}"
            );
        }
    }

    #[test]
    fn grid_lod_blend_t_is_zero_at_band_start() {
        // At distance=10 (raw=1.0, exactly the start of the 1.0 band), blend_t=0.
        let (_, t) = grid_lod_params(10.0);
        assert!(t < 0.05, "blend_t at band start should be near 0, got {t}");
    }

    #[test]
    fn grid_lod_blend_t_rises_within_band() {
        // As distance increases within a band the blend_t should increase.
        let (_, t_lo) = grid_lod_params(10.0);
        let (_, t_hi) = grid_lod_params(12.0);
        assert!(t_hi > t_lo, "blend_t should rise as distance grows, lo={t_lo}, hi={t_hi}");
    }

    #[test]
    fn grid_lod_blend_t_in_0_1() {
        for d in [0.1_f32, 1.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0] {
            let (_, t) = grid_lod_params(d);
            assert!((0.0..=1.0).contains(&t), "blend_t out of [0,1] at d={d}: {t}");
        }
    }

    #[test]
    fn grid_lod_continuity_across_band_boundary() {
        // blend_t should reset to ~0 after a band transition, and the
        // minor spacing should have snapped. The visual grid is continuous
        // because the old minor (fully faded at t→1) equals the old major
        // that the next band starts showing at t=0.
        let (m_before, t_before) = grid_lod_params(14.9); // just before band transition at 15.0
        let (m_after, t_after) = grid_lod_params(15.1); // just after
        assert!(
            t_before > 0.9,
            "blend_t should be near 1.0 just before a band boundary, got {t_before}"
        );
        assert!(
            t_after < 0.1,
            "blend_t should be near 0.0 just after a band boundary, got {t_after}"
        );
        assert!(
            m_after > m_before,
            "minor spacing should increase at band boundary: {m_before} → {m_after}"
        );
    }

    #[test]
    fn center_ray_hits_ground_near_target() {
        let cam = OrbitCamera {
            target: Point3::origin(),
            ..Default::default()
        };
        let r = ray_from_screen(&cam, 100.0, 100.0, [50.0, 50.0]);
        let hit = intersect_ground_y0(&r).expect("center ray should hit ground");
        assert!(hit.y.abs() < 1e-3, "hit.y = {}", hit.y);
    }

    #[test]
    fn parallel_ray_misses_ground() {
        let r = Ray {
            origin: Point3::new(0.0, 5.0, 0.0),
            dir: Vector3::new(1.0, 0.0, 0.0),
        };
        assert!(intersect_ground_y0(&r).is_none());
    }

    #[test]
    fn gizmo_faces_map_to_view_directions() {
        assert_eq!(gizmo_view_for_face(GizmoFace::Top), ViewDirection::Top);
        assert_eq!(gizmo_view_for_face(GizmoFace::Front), ViewDirection::Front);
        assert_eq!(gizmo_view_for_face(GizmoFace::Right), ViewDirection::Right);
    }

    #[test]
    fn gizmo_axis_dirs_match_front_view_convention() {
        let mut cam = OrbitCamera::default();
        cam.set_view(ViewDirection::Front); // az 0, el 0 → looking down -Z
        let d = gizmo_axis_screen_dirs(&cam);
        // +X points right on screen.
        assert!(d[0].0[0] > 0.9 && d[0].0[1].abs() < 0.1, "X dir {:?}", d[0].0);
        // +Y points up on screen (screen y is down → negative).
        assert!(d[1].0[1] < -0.9 && d[1].0[0].abs() < 0.1, "Y dir {:?}", d[1].0);
        // +Z points toward the viewer (positive view-space depth).
        assert!(d[2].1 > 0.5, "Z depth {}", d[2].1);
    }

    #[test]
    fn snap_ground_point_rounds_x_and_z_keeps_y() {
        let p = Point3::new(1.2, 0.0, 3.7);
        let s = snap_ground_point(p, 1.0);
        assert_eq!(s, Point3::new(1.0, 0.0, 4.0));
        // Half-unit grid.
        let s2 = snap_ground_point(Point3::new(0.49, 0.0, 0.51), 0.5);
        assert_eq!(s2, Point3::new(0.5, 0.0, 0.5));
        // Y is preserved, not snapped.
        let s3 = snap_ground_point(Point3::new(2.4, 0.123, -2.4), 1.0);
        assert_eq!(s3, Point3::new(2.0, 0.123, -2.0));
        // spacing <= 0 → unchanged.
        assert_eq!(snap_ground_point(p, 0.0), p);
        assert_eq!(snap_ground_point(p, -1.0), p);
    }
}
