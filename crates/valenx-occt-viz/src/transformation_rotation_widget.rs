//! Phase 193 — `AIS_Manipulator` rotation rings — interactive XYZ
//! rotation manipulator (drag-ring widget).
//!
//! ## What OCCT does
//!
//! Same family as [`crate::transformation_local_axis_widget()`] but
//! with three colour-coded circles instead of arrows. Hovering a ring
//! shows its rotation axis; dragging projects the cursor onto the
//! ring's tangent and rotates the target about the corresponding
//! axis. Snap-to-15° is a common option (held with Shift in OCCT;
//! Valenx will gain the same modifier).
//!
//! ## v1 status — real ring geometry + drag-state machine
//!
//! Both renderer-independent halves are shipped:
//!
//! - **Geometry.** [`RotationGizmo::ring_geometry`] builds the actual
//!   circle polyline for any axis, ready to wrap in a
//!   [`crate::ais_interactive_context::Pickable::Polyline`].
//! - **Drag math.** The [`RotationGizmo`] drag-state machine
//!   intersects the cursor ray with the ring's plane (normal = axis),
//!   measures the polar angle of the hit relative to a reference
//!   direction, and reports the signed angle swept since drag start.
//!   Optional 15° snapping mirrors OCCT's Shift modifier.
//!
//! The entry point [`transformation_rotation_widget`] constructs an
//! active gizmo at a world pose.

use nalgebra::Vector3;

use crate::ais_interactive_context::Ray;
use crate::error::OcctVizError;
use crate::transformation_local_axis_widget::GizmoAxis;

/// Rotation manipulator — three axis rings at a world origin, with a
/// single-axis drag-state machine that reports a swept angle.
#[derive(Clone, Debug)]
pub struct RotationGizmo {
    /// World-space origin the rings are centred on.
    pub origin: Vector3<f64>,
    /// Ring radius.
    pub radius: f64,
    /// Snap the reported angle to 15° increments when `true`.
    pub snap_15deg: bool,
    /// Active drag, if a drag is in progress.
    drag: Option<RotDragState>,
}

/// State captured at `begin_drag`.
#[derive(Copy, Clone, Debug)]
struct RotDragState {
    axis: GizmoAxis,
    /// Polar angle (radians) of the cursor at the drag start.
    start_angle: f64,
}

impl RotationGizmo {
    /// A gizmo at `origin` with the default ring radius and no snap.
    pub fn new(origin: Vector3<f64>) -> Self {
        Self {
            origin,
            radius: 1.0,
            snap_15deg: false,
            drag: None,
        }
    }

    /// Ring polyline for one axis — a closed circle of `segments`
    /// points in the plane perpendicular to the axis. Returned as a
    /// point list ready to wrap in a
    /// [`crate::ais_interactive_context::Pickable::Polyline`] (the
    /// last point repeats the first so the loop closes).
    pub fn ring_geometry(&self, axis: GizmoAxis) -> Vec<[f64; 3]> {
        const SEGMENTS: usize = 48;
        let (u, v) = ring_basis(axis);
        let mut pts: Vec<[f64; 3]> = Vec::with_capacity(SEGMENTS + 1);
        for k in 0..=SEGMENTS {
            let a = std::f64::consts::TAU * k as f64 / SEGMENTS as f64;
            let p = self.origin + (u * a.cos() + v * a.sin()) * self.radius;
            pts.push([p.x, p.y, p.z]);
        }
        pts
    }

    /// Begin a drag on `axis`. `ray` is the cursor ray at click. If
    /// the ray is parallel to the ring plane (a grazing view) the
    /// drag still starts but the first angle is taken as 0.
    pub fn begin_drag(&mut self, axis: GizmoAxis, ray: &Ray) {
        let start_angle = self.cursor_angle(axis, ray).unwrap_or(0.0);
        self.drag = Some(RotDragState { axis, start_angle });
    }

    /// Update the active drag with the current cursor `ray`. Returns
    /// the signed rotation angle (radians) swept since the drag start,
    /// or `None` if no drag is active. With [`snap_15deg`](Self::snap_15deg)
    /// set the angle is snapped to the nearest 15°.
    pub fn update_drag(&mut self, ray: &Ray) -> Option<f64> {
        let drag = self.drag?;
        let now = self.cursor_angle(drag.axis, ray)?;
        let mut delta = wrap_angle(now - drag.start_angle);
        if self.snap_15deg {
            let step = std::f64::consts::FRAC_PI_2 / 6.0; // 15°
            delta = (delta / step).round() * step;
        }
        Some(delta)
    }

    /// End the active drag, returning the final swept angle, or `None`
    /// if no drag was active.
    pub fn end_drag(&mut self, last_ray: &Ray) -> Option<f64> {
        let result = self.update_drag(last_ray);
        self.drag = None;
        result
    }

    /// Whether a drag is currently in progress.
    pub fn is_dragging(&self) -> bool {
        self.drag.is_some()
    }

    /// Polar angle of the cursor ray's intersection with the ring's
    /// plane, measured in the ring's `(u, v)` basis. `None` when the
    /// ray is parallel to the plane.
    fn cursor_angle(&self, axis: GizmoAxis, ray: &Ray) -> Option<f64> {
        let normal = axis.direction();
        // Intersect the ray with the plane (normal, through origin).
        let denom = ray.direction.dot(&normal);
        if denom.abs() < 1e-9 {
            return None; // ray parallel to the ring plane
        }
        let t = (self.origin - ray.origin).dot(&normal) / denom;
        let hit = ray.origin + ray.direction * t;
        let rel = hit - self.origin;
        let (u, v) = ring_basis(axis);
        let x = rel.dot(&u);
        let y = rel.dot(&v);
        Some(y.atan2(x))
    }
}

/// Orthonormal basis `(u, v)` spanning the plane perpendicular to the
/// axis. Consistent so `ring_geometry` and `cursor_angle` agree.
fn ring_basis(axis: GizmoAxis) -> (Vector3<f64>, Vector3<f64>) {
    let n = axis.direction();
    let helper = if n.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let u = n.cross(&helper).normalize();
    let v = n.cross(&u);
    (u, v)
}

/// Wrap an angle into `(-π, π]`.
fn wrap_angle(mut a: f64) -> f64 {
    use std::f64::consts::{PI, TAU};
    while a > PI {
        a -= TAU;
    }
    while a <= -PI {
        a += TAU;
    }
    a
}

/// Construct an active rotation manipulator at the world pose
/// `(target_x, target_y, target_z)` and return it.
///
/// # Errors
///
/// [`OcctVizError::BadInput`] if any coordinate is non-finite.
pub fn transformation_rotation_widget(
    target_x: f32,
    target_y: f32,
    target_z: f32,
) -> Result<RotationGizmo, OcctVizError> {
    for (n, v) in [
        ("target_x", target_x),
        ("target_y", target_y),
        ("target_z", target_z),
    ] {
        if !v.is_finite() {
            return Err(OcctVizError::bad_input(n, "must be finite"));
        }
    }
    Ok(RotationGizmo::new(Vector3::new(
        target_x as f64,
        target_y as f64,
        target_z as f64,
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn rejects_inf_z() {
        let err = transformation_rotation_widget(0.0, 0.0, f32::INFINITY).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn valid_input_builds_a_gizmo() {
        let g = transformation_rotation_widget(1.0, 2.0, 3.0).unwrap();
        assert!((g.origin - Vector3::new(1.0, 2.0, 3.0)).norm() < 1e-9);
        assert!(!g.is_dragging());
    }

    #[test]
    fn ring_geometry_is_a_closed_circle() {
        let g = RotationGizmo::new(Vector3::zeros());
        let ring = g.ring_geometry(GizmoAxis::Z);
        // 48 segments → 49 points (loop closed).
        assert_eq!(ring.len(), 49);
        // First and last coincide.
        let first = ring[0];
        let last = ring[48];
        assert!(
            (first[0] - last[0]).abs() < 1e-9 && (first[1] - last[1]).abs() < 1e-9,
            "ring loop should close"
        );
        // The Z-axis ring lies in the z=0 plane and has radius 1.
        for p in &ring {
            assert!(p[2].abs() < 1e-9, "Z ring off the z=0 plane: {}", p[2]);
            let r = (p[0] * p[0] + p[1] * p[1]).sqrt();
            assert!((r - 1.0).abs() < 1e-9, "ring radius drift: {r}");
        }
    }

    #[test]
    fn quarter_turn_drag_about_z_reports_90deg() {
        let mut g = RotationGizmo::new(Vector3::zeros());
        // Drag the Z ring. Rays come straight down -Z so they hit the
        // z=0 plane directly at their (x, y).
        let ray_at = |x: f64, y: f64| Ray {
            origin: Vector3::new(x, y, 5.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        // Start: cursor on +X side of the ring.
        g.begin_drag(GizmoAxis::Z, &ray_at(1.0, 0.0));
        // End: cursor on +Y side — a 90° sweep.
        let delta = g.update_drag(&ray_at(0.0, 1.0)).unwrap();
        assert!(
            (delta.abs() - FRAC_PI_2).abs() < 1e-6,
            "expected ~90°, got {delta} rad"
        );
    }

    #[test]
    fn snap_rounds_to_15_degrees() {
        let mut g = RotationGizmo::new(Vector3::zeros());
        g.snap_15deg = true;
        let ray_at = |x: f64, y: f64| Ray {
            origin: Vector3::new(x, y, 5.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        g.begin_drag(GizmoAxis::Z, &ray_at(1.0, 0.0));
        // Cursor swept to ~20° — should snap to 15°.
        let a = 20f64.to_radians();
        let delta = g.update_drag(&ray_at(a.cos(), a.sin())).unwrap();
        let snapped_to = 15f64.to_radians();
        assert!(
            (delta.abs() - snapped_to).abs() < 1e-6,
            "expected snap to 15°, got {} deg",
            delta.to_degrees()
        );
    }

    #[test]
    fn end_drag_clears_state() {
        let mut g = RotationGizmo::new(Vector3::zeros());
        let ray_at = |x: f64, y: f64| Ray {
            origin: Vector3::new(x, y, 5.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        g.begin_drag(GizmoAxis::Z, &ray_at(1.0, 0.0));
        let _ = g.end_drag(&ray_at(-1.0, 0.0));
        assert!(!g.is_dragging());
    }

    #[test]
    fn update_with_no_drag_returns_none() {
        let mut g = RotationGizmo::new(Vector3::zeros());
        let ray = Ray {
            origin: Vector3::new(1.0, 0.0, 5.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        assert!(g.update_drag(&ray).is_none());
    }

    #[test]
    fn wrap_angle_keeps_range() {
        assert!((wrap_angle(3.0 * PI) - PI).abs() < 1e-9);
        assert!((wrap_angle(-3.0 * PI) - PI).abs() < 1e-9);
        assert!((wrap_angle(0.5) - 0.5).abs() < 1e-12);
    }
}
