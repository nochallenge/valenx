//! [`View`] and [`ViewKind`] — a single projected drawing view.
//!
//! A view is one window onto a 3D solid: it picks a camera (front,
//! top, iso, custom), projects the solid through that camera onto the
//! drawing plane, classifies edges as visible or hidden, and places
//! the result at a `(x, y)` position on the sheet at a chosen scale.
//!
//! Heavy lifting (edge extraction, projection math, HLR) lives in the
//! sibling [`crate::projection`] and [`crate::hlr`] modules — this
//! file owns the plain-data view definition + [`View::generate`]
//! which is the convenience wrapper that calls both.

use nalgebra::{Matrix4, Vector3};
use serde::{Deserialize, Serialize};

use crate::error::TechDrawError;
use crate::{hlr, projection};

/// Camera orientation for a view. The six orthographic cardinals plus
/// an isometric and an explicit free-form custom camera.
///
/// Each variant maps to a 4×4 view matrix that places the camera in
/// world space looking at the origin; the orthographic projection on
/// the drawing plane is then the simple `(x, y)` of the
/// view-space-transformed point. Z is preserved for HLR depth tests.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum ViewKind {
    /// Looking along -Y (the +Y face points at the viewer). World
    /// `(x, z)` maps to drawing `(x, y)`.
    Front,
    /// Looking down -Z. World `(x, y)` maps to drawing `(x, y)`.
    Top,
    /// Looking along -X. World `(y, z)` maps to drawing `(x, y)`.
    Right,
    /// Looking along +Y (opposite of [`ViewKind::Front`]).
    Back,
    /// Looking up +Z (opposite of [`ViewKind::Top`]).
    Bottom,
    /// Looking along +X (opposite of [`ViewKind::Right`]).
    Left,
    /// SE / NE / front-up isometric projection (eye at unit cube
    /// corner `(1, 1, 1)` looking at origin).
    Isometric,
    /// Free-form camera defined by eye + target + up vectors.
    Custom {
        /// Camera position in world space.
        eye: Vector3<f64>,
        /// World-space point the camera looks at.
        target: Vector3<f64>,
        /// World-space "up" vector for the camera frame.
        up: Vector3<f64>,
    },
}

impl ViewKind {
    /// Short human label for UI display.
    pub fn label(&self) -> &'static str {
        match self {
            ViewKind::Front => "Front",
            ViewKind::Top => "Top",
            ViewKind::Right => "Right",
            ViewKind::Back => "Back",
            ViewKind::Bottom => "Bottom",
            ViewKind::Left => "Left",
            ViewKind::Isometric => "Isometric",
            ViewKind::Custom { .. } => "Custom",
        }
    }

    /// 4×4 view matrix transforming world-space points into the
    /// camera's coordinate frame. The projection on the drawing plane
    /// is the resulting `(x, y)` (z is the depth used by the HLR pass).
    ///
    /// Built with the standard right-handed look-at convention: -Z is
    /// "forward", +Y is "up", +X is "right" in camera space.
    pub fn camera_matrix(&self) -> Matrix4<f64> {
        let (eye, target, up) = self.frame_vectors();
        look_at_rh(eye, target, up)
    }

    /// The (eye, target, up) used by `camera_matrix`. Pulled out so
    /// [`crate::projection`] can read the eye for backface culling
    /// without rebuilding the matrix.
    pub fn frame_vectors(&self) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>) {
        // Distance is symbolic here — any non-zero magnitude works
        // since orthographic projection ignores depth. 100 keeps the
        // numbers in a clean range.
        const D: f64 = 100.0;
        let zero = Vector3::zeros();
        match *self {
            ViewKind::Front => (
                Vector3::new(0.0, -D, 0.0),
                zero,
                Vector3::new(0.0, 0.0, 1.0),
            ),
            ViewKind::Back => (Vector3::new(0.0, D, 0.0), zero, Vector3::new(0.0, 0.0, 1.0)),
            ViewKind::Top => (Vector3::new(0.0, 0.0, D), zero, Vector3::new(0.0, 1.0, 0.0)),
            ViewKind::Bottom => (
                Vector3::new(0.0, 0.0, -D),
                zero,
                Vector3::new(0.0, -1.0, 0.0),
            ),
            ViewKind::Right => (Vector3::new(D, 0.0, 0.0), zero, Vector3::new(0.0, 0.0, 1.0)),
            ViewKind::Left => (
                Vector3::new(-D, 0.0, 0.0),
                zero,
                Vector3::new(0.0, 0.0, 1.0),
            ),
            ViewKind::Isometric => (Vector3::new(D, D, D), zero, Vector3::new(0.0, 0.0, 1.0)),
            ViewKind::Custom { eye, target, up } => (eye, target, up),
        }
    }
}

/// Right-handed look-at matrix. Identical math to
/// `nalgebra::Matrix4::look_at_rh` but written out so we don't pull a
/// dependency on `nalgebra::Isometry3` just to compute one matrix.
fn look_at_rh(eye: Vector3<f64>, target: Vector3<f64>, up: Vector3<f64>) -> Matrix4<f64> {
    let f = (target - eye).normalize(); // forward (toward target)
    let s = f.cross(&up).normalize(); // right
    let u = s.cross(&f); // recomputed up (orthonormal)
                         // Column-major: each column is one basis vector + translation row.
    let mut m = Matrix4::identity();
    m[(0, 0)] = s.x;
    m[(0, 1)] = s.y;
    m[(0, 2)] = s.z;
    m[(0, 3)] = -s.dot(&eye);
    m[(1, 0)] = u.x;
    m[(1, 1)] = u.y;
    m[(1, 2)] = u.z;
    m[(1, 3)] = -u.dot(&eye);
    m[(2, 0)] = -f.x;
    m[(2, 1)] = -f.y;
    m[(2, 2)] = -f.z;
    m[(2, 3)] = f.dot(&eye);
    m
}

/// A single view on a drawing sheet.
///
/// `visible_edges` and `hidden_edges` are populated by
/// [`View::generate`]; before the first generation pass both are
/// empty. The 2D coordinates inside the edge lists are in millimeters
/// **in the view's local frame** (origin at the view's projected
/// centroid). The export pipeline adds `position` + scale at render
/// time to place the view on the sheet.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct View {
    /// Stable index in the parent drawing. Assigned by
    /// [`crate::Drawing::add_view`] — callers shouldn't set it
    /// directly.
    pub id: usize,
    /// Camera orientation (Front / Top / Iso / Custom).
    pub kind: ViewKind,
    /// Scale factor applied at render time (`scale = 0.5` halves
    /// every edge length on the sheet). Must be > 0.
    pub scale: f64,
    /// Placement on the sheet in millimeters — the view's local
    /// origin lands here.
    pub position: [f64; 2],
    /// Index into the parent drawing's source-solid registry. The UI
    /// layer keeps a parallel `Vec<Solid>` and resolves the index at
    /// render time; the drawing itself doesn't own the solid because
    /// the same Solid is usually already owned by the feature tree.
    pub source_solid_index: usize,
    /// Visible (non-occluded) edge segments in local mm.
    pub visible_edges: Vec<[(f64, f64); 2]>,
    /// Hidden (occluded) edge segments in local mm.
    pub hidden_edges: Vec<[(f64, f64); 2]>,
}

impl View {
    /// Create an empty view at `position` on the sheet at the given
    /// scale. `source_solid_index` defaults to 0 — callers that work
    /// with multiple solids set it explicitly.
    pub fn new(kind: ViewKind, scale: f64, position: [f64; 2]) -> Self {
        Self {
            id: 0,
            kind,
            scale,
            position,
            source_solid_index: 0,
            visible_edges: Vec::new(),
            hidden_edges: Vec::new(),
        }
    }

    /// Populate `visible_edges` + `hidden_edges` by projecting `solid`
    /// through this view's camera and running the approximate HLR
    /// pass. Idempotent: calling twice with the same solid yields the
    /// same edge lists (overwrites the previous contents).
    ///
    /// Returns [`TechDrawError::BadViewParameter`] when `self.scale`
    /// is non-positive.
    pub fn generate(&mut self, solid: &valenx_cad::Solid) -> Result<(), TechDrawError> {
        if !self.scale.is_finite() || self.scale <= 0.0 {
            return Err(TechDrawError::BadViewParameter {
                name: "scale",
                reason: format!("must be > 0, got {}", self.scale),
            });
        }
        let cam = self.kind.camera_matrix();
        let (visible, hidden) = hlr::classify_edges(solid, &cam)?;
        self.visible_edges = visible;
        self.hidden_edges = hidden;
        Ok(())
    }

    /// Axis-aligned bounding box of every projected edge endpoint
    /// (both visible and hidden). Returns `None` when the view has no
    /// edges yet — useful for [`Self::auto_dimension_bbox`].
    pub fn bbox(&self) -> Option<([f64; 2], [f64; 2])> {
        let mut min = [f64::INFINITY; 2];
        let mut max = [f64::NEG_INFINITY; 2];
        let walk = |arr: &Vec<[(f64, f64); 2]>, min: &mut [f64; 2], max: &mut [f64; 2]| {
            for seg in arr {
                for p in seg {
                    if p.0 < min[0] {
                        min[0] = p.0;
                    }
                    if p.1 < min[1] {
                        min[1] = p.1;
                    }
                    if p.0 > max[0] {
                        max[0] = p.0;
                    }
                    if p.1 > max[1] {
                        max[1] = p.1;
                    }
                }
            }
        };
        walk(&self.visible_edges, &mut min, &mut max);
        walk(&self.hidden_edges, &mut min, &mut max);
        if min[0].is_finite() && max[0].is_finite() {
            Some((min, max))
        } else {
            None
        }
    }

    /// Project a world-space point through this view's camera. The
    /// returned `(x, y)` is in local-view millimeters (same frame as
    /// `visible_edges` / `hidden_edges`).
    pub fn project(&self, world: Vector3<f64>) -> [f64; 2] {
        let cam = self.kind.camera_matrix();
        projection::project_point(world, &cam)
    }

    /// Generate two simple linear dimensions covering the bounding
    /// box's width and height. Returns an empty list when the view
    /// has no projected edges yet.
    pub fn auto_dimension_bbox(&self) -> Vec<crate::Dimension> {
        let Some((min, max)) = self.bbox() else {
            return Vec::new();
        };
        let width = max[0] - min[0];
        let height = max[1] - min[1];
        vec![
            crate::Dimension::Linear {
                from: [min[0], min[1]],
                to: [max[0], min[1]],
                offset: -10.0,
                value: width,
            },
            crate::Dimension::Linear {
                from: [min[0], min[1]],
                to: [min[0], max[1]],
                offset: -10.0,
                value: height,
            },
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tiny ad-hoc tolerance helper so tests don't pull `approx` as a
    /// new workspace dependency.
    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn view_kind_labels_are_distinct() {
        let labels = [
            ViewKind::Front.label(),
            ViewKind::Top.label(),
            ViewKind::Right.label(),
            ViewKind::Back.label(),
            ViewKind::Bottom.label(),
            ViewKind::Left.label(),
            ViewKind::Isometric.label(),
        ];
        let unique: std::collections::HashSet<_> = labels.iter().collect();
        assert_eq!(unique.len(), labels.len(), "labels should be unique");
    }

    #[test]
    fn camera_matrix_front_view_projects_xz_to_xy() {
        // From the front (eye on -Y looking at origin, up = +Z), a
        // world point at (1, 0, 2) should project to drawing
        // coordinates (1, 2). The matrix transforms world → camera
        // space; we read the (x, y) of the resulting vector.
        let cam = ViewKind::Front.camera_matrix();
        let p = nalgebra::Vector4::new(1.0, 0.0, 2.0, 1.0);
        let v = cam * p;
        assert!(close(v.x, 1.0, 1e-9), "x got {}", v.x);
        assert!(close(v.y, 2.0, 1e-9), "y got {}", v.y);
    }

    #[test]
    fn camera_matrix_top_view_projects_xy_to_xy() {
        // Top view: eye on +Z looking down at origin, up = +Y.
        // World (1, 2, 0) → drawing (1, 2).
        let cam = ViewKind::Top.camera_matrix();
        let p = nalgebra::Vector4::new(1.0, 2.0, 0.0, 1.0);
        let v = cam * p;
        assert!(close(v.x, 1.0, 1e-9), "x got {}", v.x);
        assert!(close(v.y, 2.0, 1e-9), "y got {}", v.y);
    }

    #[test]
    fn view_new_starts_with_no_edges() {
        let v = View::new(ViewKind::Front, 1.0, [10.0, 20.0]);
        assert_eq!(v.scale, 1.0);
        assert_eq!(v.position, [10.0, 20.0]);
        assert!(v.visible_edges.is_empty());
        assert!(v.hidden_edges.is_empty());
    }

    #[test]
    fn generate_rejects_bad_scale() {
        use valenx_cad::primitives::box_solid;
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let mut v = View::new(ViewKind::Front, -1.0, [0.0, 0.0]);
        let e = v.generate(&cube).unwrap_err();
        match e {
            TechDrawError::BadViewParameter { name, .. } => assert_eq!(name, "scale"),
            other => panic!("wrong error: {other:?}"),
        }
    }
}
