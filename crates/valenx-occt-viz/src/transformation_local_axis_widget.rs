//! Phase 192 — `AIS_Manipulator` translation arrows — interactive XYZ
//! axis manipulator (drag-arrow widget).
//!
//! ## What OCCT does
//!
//! `AIS_Manipulator` attaches three colour-coded arrows (red=X,
//! green=Y, blue=Z) to a target shape's origin. Hovering an arrow
//! highlights it; dragging projects the cursor onto the arrow's
//! axis and translates the target along it. Visually it's a tri-axis
//! "gizmo" — same UX as Blender/Maya/3ds Max object-mode translate.
//!
//! ## v1 status — real gizmo geometry + constraint-aware drag
//!
//! Both halves of the manipulator that don't need a live GPU are
//! shipped here:
//!
//! - **Geometry.** [`TranslationGizmo::axis_geometry`] builds the
//!   actual arrow mesh for any axis — a cylindrical shaft plus a
//!   conical head, faceted into triangles. [`TranslationGizmo::plane_geometry`]
//!   builds the small corner quad for a plane handle. The renderer
//!   uploads them to an overlay pass; the picker ray-casts them (a
//!   [`crate::ais_interactive_context::Pickable::Mesh`]).
//! - **Constraint-aware drag math.** The [`TranslationGizmo`]
//!   drag-state machine constrains motion two ways:
//!   - **Axis** — [`begin_drag`](TranslationGizmo::begin_drag)
//!     projects the cursor ray onto the picked axis line (the
//!     closest-point-of-two-lines solve), so motion is locked to one
//!     axis. This is OCCT's arrow-handle drag.
//!   - **Plane** — [`begin_plane_drag`](TranslationGizmo::begin_plane_drag)
//!     intersects the cursor ray with the picked coordinate plane, so
//!     motion is locked to two axes. This is OCCT's plane-handle
//!     drag.
//!
//!   Either drag mode optionally **snaps to a grid increment**
//!   ([`snap_increment`](TranslationGizmo::snap_increment)) — the
//!   translation analogue of the rotation gizmo's 15° snap — so the
//!   gizmo can place geometry at exact coordinates, the constraint
//!   that makes a manipulator useful for precise CAD work.
//!
//! The thin entry point [`transformation_local_axis_widget`]
//! constructs an active gizmo at a world pose.
//!
//! ### Honest scope
//!
//! "Constraint-aware" here means the **motion constraints** of the
//! manipulator itself — axis lock, plane lock, grid snap. It does
//! *not* mean assembly-constraint propagation (dragging one part and
//! having mates / joints re-solve the rest of the assembly): that is
//! a constraint-solver problem on the assembly graph and stays a
//! genuine Tier-3 item.

use nalgebra::Vector3;

use crate::ais_interactive_context::Ray;
use crate::error::OcctVizError;

/// The three manipulator axes.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, serde::Serialize, serde::Deserialize)]
pub enum GizmoAxis {
    /// World +X (red arrow).
    X,
    /// World +Y (green arrow).
    Y,
    /// World +Z (blue arrow).
    Z,
}

impl GizmoAxis {
    /// Unit direction of the axis.
    pub fn direction(self) -> Vector3<f64> {
        match self {
            GizmoAxis::X => Vector3::new(1.0, 0.0, 0.0),
            GizmoAxis::Y => Vector3::new(0.0, 1.0, 0.0),
            GizmoAxis::Z => Vector3::new(0.0, 0.0, 1.0),
        }
    }

    /// Conventional gizmo colour (RGB, 0..1).
    pub fn color(self) -> [f32; 3] {
        match self {
            GizmoAxis::X => [0.9, 0.2, 0.2],
            GizmoAxis::Y => [0.2, 0.8, 0.2],
            GizmoAxis::Z => [0.25, 0.45, 0.95],
        }
    }
}

/// The three manipulator planes — the plane-constrained drag handles.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, serde::Serialize, serde::Deserialize)]
pub enum GizmoPlane {
    /// The world XY plane (normal +Z) — motion locked to X and Y.
    Xy,
    /// The world YZ plane (normal +X) — motion locked to Y and Z.
    Yz,
    /// The world ZX plane (normal +Y) — motion locked to Z and X.
    Zx,
}

impl GizmoPlane {
    /// Unit normal of the plane.
    pub fn normal(self) -> Vector3<f64> {
        match self {
            GizmoPlane::Xy => Vector3::new(0.0, 0.0, 1.0),
            GizmoPlane::Yz => Vector3::new(1.0, 0.0, 0.0),
            GizmoPlane::Zx => Vector3::new(0.0, 1.0, 0.0),
        }
    }

    /// The two in-plane axis directions `(a, b)` spanning the plane.
    pub fn axes(self) -> (Vector3<f64>, Vector3<f64>) {
        match self {
            GizmoPlane::Xy => (GizmoAxis::X.direction(), GizmoAxis::Y.direction()),
            GizmoPlane::Yz => (GizmoAxis::Y.direction(), GizmoAxis::Z.direction()),
            GizmoPlane::Zx => (GizmoAxis::Z.direction(), GizmoAxis::X.direction()),
        }
    }
}

/// Translation manipulator — three axis arrows + three plane handles
/// at a world origin, with a constraint-aware drag-state machine.
#[derive(Clone, Debug)]
pub struct TranslationGizmo {
    /// World-space origin the arrows emanate from.
    pub origin: Vector3<f64>,
    /// Overall length of each arrow (shaft + head).
    pub size: f64,
    /// When `Some(step)`, the reported translation is snapped to a
    /// grid of spacing `step` (per axis). `None` is free motion. The
    /// translation analogue of the rotation gizmo's 15° snap.
    pub snap_increment: Option<f64>,
    /// Active drag, if a drag is in progress.
    drag: Option<DragState>,
}

/// What the active drag is constrained to — one axis, or one plane.
#[derive(Copy, Clone, Debug)]
enum DragConstraint {
    /// Locked to a single axis. Carries the cursor parameter `t`
    /// along the axis line at the drag start.
    Axis { axis: GizmoAxis, start_t: f64 },
    /// Locked to a coordinate plane. Carries the cursor's in-plane
    /// hit point at the drag start.
    Plane {
        plane: GizmoPlane,
        start_hit: Vector3<f64>,
    },
}

/// State captured at the start of a drag so updates can report a
/// delta relative to the drag origin.
#[derive(Copy, Clone, Debug)]
struct DragState {
    constraint: DragConstraint,
    /// Gizmo origin at the drag start.
    start_origin: Vector3<f64>,
}

impl TranslationGizmo {
    /// A gizmo at `origin` with the default arrow size and free
    /// (un-snapped) motion.
    pub fn new(origin: Vector3<f64>) -> Self {
        Self {
            origin,
            size: 1.0,
            snap_increment: None,
            drag: None,
        }
    }

    /// Arrow mesh for one axis — a faceted shaft cylinder plus a
    /// conical head — returned as a flat triangle list (three points
    /// per triangle) ready to wrap in a
    /// [`crate::ais_interactive_context::Pickable::Mesh`].
    pub fn axis_geometry(&self, axis: GizmoAxis) -> Vec<[f64; 3]> {
        const SEGMENTS: usize = 12;
        let dir = axis.direction();
        // Build an orthonormal frame (dir, u, v).
        let helper = if dir.x.abs() < 0.9 {
            Vector3::new(1.0, 0.0, 0.0)
        } else {
            Vector3::new(0.0, 1.0, 0.0)
        };
        let u = dir.cross(&helper).normalize();
        let v = dir.cross(&u);

        let shaft_len = self.size * 0.8;
        let shaft_r = self.size * 0.03;
        let head_len = self.size * 0.2;
        let head_r = self.size * 0.08;
        let o = self.origin;

        let ring = |center: Vector3<f64>, radius: f64| -> Vec<Vector3<f64>> {
            (0..SEGMENTS)
                .map(|k| {
                    let a = std::f64::consts::TAU * k as f64 / SEGMENTS as f64;
                    center + (u * a.cos() + v * a.sin()) * radius
                })
                .collect()
        };

        let base_ring = ring(o, shaft_r);
        let top_ring = ring(o + dir * shaft_len, shaft_r);
        let head_base = ring(o + dir * shaft_len, head_r);
        let tip = o + dir * (shaft_len + head_len);

        let mut tris: Vec<[f64; 3]> = Vec::new();
        let push = |tris: &mut Vec<[f64; 3]>, a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>| {
            tris.push([a.x, a.y, a.z]);
            tris.push([b.x, b.y, b.z]);
            tris.push([c.x, c.y, c.z]);
        };
        // Shaft side wall.
        for k in 0..SEGMENTS {
            let n = (k + 1) % SEGMENTS;
            push(&mut tris, base_ring[k], base_ring[n], top_ring[n]);
            push(&mut tris, base_ring[k], top_ring[n], top_ring[k]);
        }
        // Cone side wall.
        for k in 0..SEGMENTS {
            let n = (k + 1) % SEGMENTS;
            push(&mut tris, head_base[k], head_base[n], tip);
        }
        tris
    }

    /// Plane handle geometry for one coordinate plane — a small square
    /// quad set back from the origin in the plane, returned as a flat
    /// triangle list (two triangles) ready to wrap in a
    /// [`crate::ais_interactive_context::Pickable::Mesh`]. This is the
    /// clickable handle for a plane-constrained drag.
    pub fn plane_geometry(&self, plane: GizmoPlane) -> Vec<[f64; 3]> {
        let (a, b) = plane.axes();
        // A quad inset from the origin: corner at 0.25·size, side
        // 0.25·size — the conventional small offset square.
        let near = self.size * 0.25;
        let far = self.size * 0.5;
        let o = self.origin;
        let c00 = o + a * near + b * near;
        let c10 = o + a * far + b * near;
        let c11 = o + a * far + b * far;
        let c01 = o + a * near + b * far;
        let p = |v: Vector3<f64>| [v.x, v.y, v.z];
        vec![
            p(c00), p(c10), p(c11), // triangle 1
            p(c00), p(c11), p(c01), // triangle 2
        ]
    }

    /// Begin an axis-constrained drag on `axis`. `ray` is the cursor
    /// ray at the moment of click. Subsequent
    /// [`update_drag`](Self::update_drag) calls report the translation
    /// delta relative to this start, locked to the axis.
    pub fn begin_drag(&mut self, axis: GizmoAxis, ray: &Ray) {
        let start_t = closest_param_on_axis(self.origin, axis.direction(), ray);
        self.drag = Some(DragState {
            constraint: DragConstraint::Axis { axis, start_t },
            start_origin: self.origin,
        });
    }

    /// Begin a plane-constrained drag on `plane`. `ray` is the cursor
    /// ray at click. Subsequent [`update_drag`](Self::update_drag)
    /// calls report a translation delta that stays *in* the plane
    /// (the third axis is locked).
    ///
    /// If the cursor ray grazes the plane (nearly parallel) the drag
    /// still begins but the start hit falls back to the gizmo origin.
    pub fn begin_plane_drag(&mut self, plane: GizmoPlane, ray: &Ray) {
        let start_hit =
            ray_plane_hit(self.origin, plane.normal(), ray).unwrap_or(self.origin);
        self.drag = Some(DragState {
            constraint: DragConstraint::Plane { plane, start_hit },
            start_origin: self.origin,
        });
    }

    /// Update the active drag with the current cursor `ray`, moving
    /// the gizmo origin under the active constraint. Returns the
    /// *total* translation vector from the drag start, or `None` if no
    /// drag is active.
    ///
    /// When [`snap_increment`](Self::snap_increment) is set the delta
    /// is rounded to the grid before being applied, so the gizmo lands
    /// on exact coordinates.
    pub fn update_drag(&mut self, ray: &Ray) -> Option<Vector3<f64>> {
        let drag = self.drag?;
        let raw_delta = match drag.constraint {
            DragConstraint::Axis { axis, start_t } => {
                let dir = axis.direction();
                let now_t = closest_param_on_axis(drag.start_origin, dir, ray);
                dir * (now_t - start_t)
            }
            DragConstraint::Plane { plane, start_hit } => {
                // Intersect the ray with the plane through the *drag
                // start origin* — the constraint plane does not move
                // with the gizmo during the drag.
                let hit = ray_plane_hit(drag.start_origin, plane.normal(), ray)
                    .unwrap_or(start_hit);
                hit - start_hit
            }
        };
        let delta = self.snap_delta(raw_delta);
        self.origin = drag.start_origin + delta;
        Some(delta)
    }

    /// End the active drag. Returns the final translation, or `None`
    /// if no drag was active.
    pub fn end_drag(&mut self) -> Option<Vector3<f64>> {
        let drag = self.drag.take()?;
        Some(self.origin - drag.start_origin)
    }

    /// Whether a drag is currently in progress.
    pub fn is_dragging(&self) -> bool {
        self.drag.is_some()
    }

    /// Round a translation delta to the grid set by
    /// [`snap_increment`](Self::snap_increment), per component. With
    /// no increment set the delta passes through unchanged.
    fn snap_delta(&self, delta: Vector3<f64>) -> Vector3<f64> {
        match self.snap_increment {
            Some(step) if step > 0.0 => Vector3::new(
                (delta.x / step).round() * step,
                (delta.y / step).round() * step,
                (delta.z / step).round() * step,
            ),
            _ => delta,
        }
    }
}

/// Closest point of the cursor ray to the axis line through
/// `(origin, dir)`, returned as the signed parameter `t` along the
/// axis (so the world point is `origin + dir * t`).
fn closest_param_on_axis(origin: Vector3<f64>, dir: Vector3<f64>, ray: &Ray) -> f64 {
    // Two-line closest-point solve: line A = axis, line B = ray.
    let d1 = dir;
    let d2 = ray.direction;
    let r = origin - ray.origin;
    let a = d1.dot(&d1);
    let b = d1.dot(&d2);
    let c = d2.dot(&d2);
    let d = d1.dot(&r);
    let e = d2.dot(&r);
    let denom = a * c - b * b;
    if denom.abs() < 1e-12 {
        // Ray parallel to the axis — project the ray origin instead.
        return -d / a.max(1e-12);
    }
    // t on the axis line.
    (b * e - c * d) / denom
}

/// World hit point of the cursor `ray` against the plane through
/// `plane_origin` with unit `plane_normal`. `None` when the ray is
/// (nearly) parallel to the plane.
fn ray_plane_hit(
    plane_origin: Vector3<f64>,
    plane_normal: Vector3<f64>,
    ray: &Ray,
) -> Option<Vector3<f64>> {
    let denom = ray.direction.dot(&plane_normal);
    if denom.abs() < 1e-9 {
        return None; // grazing — ray parallel to the plane
    }
    let t = (plane_origin - ray.origin).dot(&plane_normal) / denom;
    Some(ray.origin + ray.direction * t)
}

/// Construct an active translation manipulator at the world pose
/// `(target_x, target_y, target_z)` and return it.
///
/// This is the entry point OCCT spells `AIS_Manipulator::SetPart` +
/// `Attach`; the returned [`TranslationGizmo`] carries the geometry
/// and the drag-state machine.
///
/// # Errors
///
/// [`OcctVizError::BadInput`] if any coordinate is non-finite.
pub fn transformation_local_axis_widget(
    target_x: f32,
    target_y: f32,
    target_z: f32,
) -> Result<TranslationGizmo, OcctVizError> {
    for (n, v) in [
        ("target_x", target_x),
        ("target_y", target_y),
        ("target_z", target_z),
    ] {
        if !v.is_finite() {
            return Err(OcctVizError::bad_input(n, "must be finite"));
        }
    }
    Ok(TranslationGizmo::new(Vector3::new(
        target_x as f64,
        target_y as f64,
        target_z as f64,
    )))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_nan_y() {
        let err = transformation_local_axis_widget(0.0, f32::NAN, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_viz.bad_input");
    }

    #[test]
    fn valid_input_builds_a_gizmo() {
        let g = transformation_local_axis_widget(10.0, 20.0, 30.0).unwrap();
        assert!((g.origin - Vector3::new(10.0, 20.0, 30.0)).norm() < 1e-9);
        assert!(!g.is_dragging());
    }

    #[test]
    fn axis_geometry_produces_triangles() {
        let g = TranslationGizmo::new(Vector3::zeros());
        let tris = g.axis_geometry(GizmoAxis::X);
        // 12-segment shaft (2 tris each) + 12-segment cone (1 each) =
        // 36 triangles = 108 points.
        assert_eq!(tris.len(), 108);
        // The arrow points along +X — the furthest vertex is near the
        // tip at x ≈ size.
        let max_x = tris.iter().map(|p| p[0]).fold(f64::MIN, f64::max);
        assert!(max_x > 0.9, "tip near x=1, got {max_x}");
    }

    #[test]
    fn drag_along_x_translates_along_x() {
        let mut g = TranslationGizmo::new(Vector3::zeros());
        // A ray that crosses the X axis at x=0 (pointing -Z from above
        // the origin).
        let ray0 = Ray {
            origin: Vector3::new(0.0, 0.0, 5.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        g.begin_drag(GizmoAxis::X, &ray0);
        assert!(g.is_dragging());
        // A ray that crosses the X axis at x=3.
        let ray1 = Ray {
            origin: Vector3::new(3.0, 0.0, 5.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        let delta = g.update_drag(&ray1).unwrap();
        assert!((delta - Vector3::new(3.0, 0.0, 0.0)).norm() < 1e-6, "delta={delta:?}");
        // The origin moved along X only.
        assert!((g.origin - Vector3::new(3.0, 0.0, 0.0)).norm() < 1e-6);
    }

    #[test]
    fn drag_locks_to_the_axis_ignoring_perpendicular_motion() {
        let mut g = TranslationGizmo::new(Vector3::zeros());
        let ray0 = Ray {
            origin: Vector3::new(0.0, 0.0, 5.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        g.begin_drag(GizmoAxis::X, &ray0);
        // Cursor moved diagonally — crosses X axis at x=2 but also has
        // Y offset; the gizmo must only move along X.
        let ray1 = Ray {
            origin: Vector3::new(2.0, 9.0, 5.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        let delta = g.update_drag(&ray1).unwrap();
        assert!(delta.y.abs() < 1e-6, "Y must stay locked, got {}", delta.y);
        assert!(delta.z.abs() < 1e-6, "Z must stay locked, got {}", delta.z);
        assert!((delta.x - 2.0).abs() < 1e-6);
    }

    #[test]
    fn end_drag_returns_total_and_clears_state() {
        let mut g = TranslationGizmo::new(Vector3::zeros());
        let ray0 = Ray {
            origin: Vector3::new(0.0, 0.0, 5.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        g.begin_drag(GizmoAxis::Z, &ray0);
        // For Z-axis drag we need a ray not parallel to Z.
        let ray1 = Ray {
            origin: Vector3::new(5.0, 0.0, 4.0),
            direction: Vector3::new(-1.0, 0.0, 0.0),
        };
        g.update_drag(&ray1);
        let total = g.end_drag().unwrap();
        assert!(!g.is_dragging());
        // Movement was locked to Z.
        assert!(total.x.abs() < 1e-6 && total.y.abs() < 1e-6);
    }

    #[test]
    fn update_with_no_drag_returns_none() {
        let mut g = TranslationGizmo::new(Vector3::zeros());
        let ray = Ray {
            origin: Vector3::new(0.0, 0.0, 5.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        assert!(g.update_drag(&ray).is_none());
        assert!(g.end_drag().is_none());
    }

    #[test]
    fn plane_geometry_lies_in_the_plane() {
        let g = TranslationGizmo::new(Vector3::zeros());
        // The XY plane handle: every triangle vertex has z = 0.
        let xy = g.plane_geometry(GizmoPlane::Xy);
        assert_eq!(xy.len(), 6, "plane handle is a 2-triangle quad");
        for v in &xy {
            assert!(v[2].abs() < 1e-9, "XY handle vertex off the plane: {v:?}");
        }
        // The YZ plane handle: every vertex has x = 0.
        let yz = g.plane_geometry(GizmoPlane::Yz);
        for v in &yz {
            assert!(v[0].abs() < 1e-9, "YZ handle vertex off the plane: {v:?}");
        }
    }

    #[test]
    fn plane_drag_moves_within_the_plane_only() {
        // Drag the XY plane handle. The cursor rays come straight down
        // -Z so they hit the z=0 plane at their (x, y). The gizmo must
        // translate in X and Y but never in Z.
        let mut g = TranslationGizmo::new(Vector3::zeros());
        let ray_at = |x: f64, y: f64| Ray {
            origin: Vector3::new(x, y, 8.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        g.begin_plane_drag(GizmoPlane::Xy, &ray_at(0.0, 0.0));
        assert!(g.is_dragging());
        // Cursor moved to (3, -2) in the plane.
        let delta = g.update_drag(&ray_at(3.0, -2.0)).unwrap();
        assert!((delta.x - 3.0).abs() < 1e-6, "X delta {}", delta.x);
        assert!((delta.y + 2.0).abs() < 1e-6, "Y delta {}", delta.y);
        assert!(delta.z.abs() < 1e-6, "Z must stay locked in a plane drag");
    }

    #[test]
    fn plane_drag_on_yz_translates_y_and_z() {
        // Drag the YZ plane handle with rays travelling -X, so they hit
        // the x=0 plane at their (y, z).
        let mut g = TranslationGizmo::new(Vector3::zeros());
        let ray_at = |y: f64, z: f64| Ray {
            origin: Vector3::new(8.0, y, z),
            direction: Vector3::new(-1.0, 0.0, 0.0),
        };
        g.begin_plane_drag(GizmoPlane::Yz, &ray_at(0.0, 0.0));
        let delta = g.update_drag(&ray_at(5.0, 4.0)).unwrap();
        assert!(delta.x.abs() < 1e-6, "X must stay locked");
        assert!((delta.y - 5.0).abs() < 1e-6);
        assert!((delta.z - 4.0).abs() < 1e-6);
    }

    #[test]
    fn snap_increment_quantises_an_axis_drag() {
        let mut g = TranslationGizmo::new(Vector3::zeros());
        g.snap_increment = Some(1.0); // 1-unit grid
        let ray_at = |x: f64| Ray {
            origin: Vector3::new(x, 0.0, 5.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        g.begin_drag(GizmoAxis::X, &ray_at(0.0));
        // Cursor at x = 2.7 → should snap to the grid value 3.0.
        let delta = g.update_drag(&ray_at(2.7)).unwrap();
        assert!(
            (delta.x - 3.0).abs() < 1e-9,
            "2.7 should snap to 3.0, got {}",
            delta.x
        );
        // Cursor at x = 2.2 → should snap down to 2.0.
        let delta = g.update_drag(&ray_at(2.2)).unwrap();
        assert!((delta.x - 2.0).abs() < 1e-9, "2.2 should snap to 2.0");
    }

    #[test]
    fn snap_increment_quantises_a_plane_drag() {
        let mut g = TranslationGizmo::new(Vector3::zeros());
        g.snap_increment = Some(0.5);
        let ray_at = |x: f64, y: f64| Ray {
            origin: Vector3::new(x, y, 8.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        g.begin_plane_drag(GizmoPlane::Xy, &ray_at(0.0, 0.0));
        // Cursor at (1.3, 0.9) → snaps to (1.5, 1.0) on the 0.5 grid.
        let delta = g.update_drag(&ray_at(1.3, 0.9)).unwrap();
        assert!((delta.x - 1.5).abs() < 1e-9, "X snap, got {}", delta.x);
        assert!((delta.y - 1.0).abs() < 1e-9, "Y snap, got {}", delta.y);
    }

    #[test]
    fn no_snap_increment_leaves_motion_continuous() {
        // With snap_increment None the delta is exact (not quantised).
        let mut g = TranslationGizmo::new(Vector3::zeros());
        assert!(g.snap_increment.is_none());
        let ray_at = |x: f64| Ray {
            origin: Vector3::new(x, 0.0, 5.0),
            direction: Vector3::new(0.0, 0.0, -1.0),
        };
        g.begin_drag(GizmoAxis::X, &ray_at(0.0));
        let delta = g.update_drag(&ray_at(2.7)).unwrap();
        assert!((delta.x - 2.7).abs() < 1e-6, "free motion should be exact");
    }

    #[test]
    fn gizmo_plane_axes_are_orthonormal_and_span_the_plane() {
        for plane in [GizmoPlane::Xy, GizmoPlane::Yz, GizmoPlane::Zx] {
            let (a, b) = plane.axes();
            let n = plane.normal();
            assert!((a.norm() - 1.0).abs() < 1e-12);
            assert!((b.norm() - 1.0).abs() < 1e-12);
            assert!(a.dot(&b).abs() < 1e-12, "in-plane axes must be ⟂");
            // a × b should equal the plane normal (right-handed).
            assert!((a.cross(&b) - n).norm() < 1e-12);
        }
    }
}
