//! Phase 171 — `AIS_InteractiveContext` — central registry of
//! selectable objects with hover + click states.
//!
//! ## What OCCT does
//!
//! `AIS_InteractiveContext` is the lookup table the viewer consults
//! every frame to decide which objects are visible, hover-highlighted,
//! selected, or owner-attached (sub-shape selection). It owns the
//! presentation list (`PrsMgr_PresentableObject`), the selection
//! manager (`SelectMgr_SelectionManager`), and the highlight states.
//! Typical lifecycle: `Display(obj)` → user hovers → `MoveTo(x, y)`
//! updates `DynamicHighlight` → click → `Select()` promotes the
//! highlight to a selection.
//!
//! ## v1 status — real registry + CPU picking infrastructure
//!
//! The registry is keyed by integer IDs (allocated monotonically)
//! with four per-object states (`Hidden` / `Visible` / `Hovered` /
//! `Selected`).
//!
//! Phase 171-180 adds the **picking substrate** that the
//! `ais_select_*` family is built on. OCCT's default selection path
//! (`SelectionMode 0`) is a CPU ray-cast against each object's
//! registered selection primitives — *not* a GPU pass. This module
//! ships exactly that:
//!
//! - Each object may carry a [`Pickable`] geometry payload — a
//!   triangle mesh (face owners), a polyline (edge owners), or a
//!   point (vertex owners). Register it with [`InteractiveContext::display_geometry`].
//! - The context holds the current camera as a [`PickView`]
//!   (view-projection matrix + viewport). The renderer refreshes it
//!   every frame via [`InteractiveContext::set_view`].
//! - [`InteractiveContext::ray_at`] unprojects a screen pixel into a
//!   world-space ray; the `ais_select_*` ops ray-cast that against
//!   the registered geometry.
//!
//! Geometry-less [`InteractiveContext::display`] is retained for
//! callers (and the highlight ops) that only need the state machine.

use std::collections::HashMap;

use nalgebra::{Matrix4, Vector3, Vector4};

use crate::error::OcctVizError;

/// Per-object visibility / highlight state mirror of OCCT's
/// `AIS_DisplayStatus` + dynamic/selected highlight flags.
#[derive(
    Copy, Clone, Debug, Eq, PartialEq, Hash, Default, serde::Serialize, serde::Deserialize,
)]
#[serde(rename_all = "lowercase")]
pub enum ObjectState {
    /// Object is hidden from the viewer.
    Hidden,
    /// Object is visible but not highlighted.
    #[default]
    Visible,
    /// Cursor is hovering over the object (yellow outline in v1).
    Hovered,
    /// Object is in the active selection set (orange outline in v1).
    Selected,
}

/// Geometry payload a registered object can carry so the `ais_select_*`
/// pickers can ray-cast against it. Mirrors OCCT's three owner kinds
/// (face / edge / vertex).
///
/// The `Tagged*` variants additionally label each primitive with a
/// **sub-element index** — this is how OCCT's subshape selection
/// (`StdSelect_BRepOwner_for_Faces` / `_Edges` / `_Vertices`) works:
/// the picking pass finds the nearest triangle / segment / point and
/// reports which face / edge / vertex it belongs to. Use them when
/// registering a solid whose individual faces / edges should be
/// pickable.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub enum Pickable {
    /// Triangle-mesh geometry — flat `[a, b, c, a, b, c, …]` vertex
    /// list (one entry of three points per triangle). Picked by
    /// ray-triangle intersection (face / solid owners).
    Mesh {
        /// Triangle vertices, three per triangle.
        triangles: Vec<[f64; 3]>,
    },
    /// Polyline geometry — an ordered point list. Picked by
    /// ray-to-segment proximity (edge owners).
    Polyline {
        /// Ordered polyline vertices.
        points: Vec<[f64; 3]>,
    },
    /// A single point — picked by ray-to-point proximity (vertex
    /// owners).
    Point {
        /// World-space position.
        position: [f64; 3],
    },
    /// A triangle mesh whose triangles are labelled with a per-face
    /// index: `triangles[3i..3i+3]` belongs to face `face_tags[i]`.
    /// Picked by `ais_select_face` to return the hit face's index.
    TaggedMesh {
        /// Triangle vertices, three per triangle.
        triangles: Vec<[f64; 3]>,
        /// `face_tags[i]` is the sub-face index of triangle `i`.
        face_tags: Vec<usize>,
    },
    /// A collection of polyline edges, each labelled with an edge
    /// index. `edges[k]` is `(edge_index, ordered_points)`. Picked by
    /// `ais_select_edge` to return the hit edge's index.
    TaggedEdges {
        /// Each `(edge_index, polyline_points)`.
        edges: Vec<(usize, Vec<[f64; 3]>)>,
    },
    /// A collection of vertices, each labelled with a vertex index.
    /// `vertices[k]` is `(vertex_index, position)`. Picked by
    /// `ais_select_vertex` to return the hit vertex's index.
    TaggedVertices {
        /// Each `(vertex_index, position)`.
        vertices: Vec<(usize, [f64; 3])>,
    },
}

impl Pickable {
    /// Centroid of the geometry — used as the winding-test point for
    /// lasso/polygon selection.
    pub fn centroid(&self) -> Vector3<f64> {
        let avg = |pts: &[[f64; 3]]| -> Vector3<f64> {
            if pts.is_empty() {
                return Vector3::zeros();
            }
            let mut c = Vector3::zeros();
            for p in pts {
                c += Vector3::new(p[0], p[1], p[2]);
            }
            c / pts.len() as f64
        };
        match self {
            Pickable::Mesh { triangles } | Pickable::TaggedMesh { triangles, .. } => avg(triangles),
            Pickable::Polyline { points } => avg(points),
            Pickable::Point { position } => Vector3::new(position[0], position[1], position[2]),
            Pickable::TaggedEdges { edges } => {
                let all: Vec<[f64; 3]> = edges
                    .iter()
                    .flat_map(|(_, pts)| pts.iter().copied())
                    .collect();
                avg(&all)
            }
            Pickable::TaggedVertices { vertices } => {
                let all: Vec<[f64; 3]> = vertices.iter().map(|(_, p)| *p).collect();
                avg(&all)
            }
        }
    }
}

/// A world-space ray — origin + unit direction.
#[derive(Copy, Clone, Debug)]
pub struct Ray {
    /// Ray origin (camera eye for a perspective pick).
    pub origin: Vector3<f64>,
    /// Unit ray direction.
    pub direction: Vector3<f64>,
}

/// The current camera, captured for picking: the combined
/// view-projection matrix and the pixel viewport size. Refreshed each
/// frame by the renderer.
#[derive(Copy, Clone, Debug)]
pub struct PickView {
    /// World → clip-space matrix (`projection * view`).
    pub view_proj: Matrix4<f64>,
    /// Viewport width in pixels.
    pub width: f32,
    /// Viewport height in pixels.
    pub height: f32,
}

/// Result of a successful pick.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PickHit {
    /// ID of the picked object.
    pub id: usize,
    /// Ray parameter `t` at the hit point (distance along the ray).
    pub distance: f64,
}

/// In-memory registry of objects keyed by 0-based ID, with optional
/// per-object [`Pickable`] geometry and the current [`PickView`].
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct InteractiveContext {
    next_id: usize,
    states: HashMap<usize, ObjectState>,
    geometry: HashMap<usize, Pickable>,
    /// Transient runtime camera — not persisted (recomputed each frame
    /// by the renderer).
    #[serde(skip)]
    view: Option<PickView>,
}

impl InteractiveContext {
    /// Allocate a fresh ID and register it in the `Visible` state with
    /// no pickable geometry. Returns the new ID.
    pub fn display(&mut self) -> usize {
        let id = self.next_id;
        self.next_id += 1;
        self.states.insert(id, ObjectState::Visible);
        id
    }

    /// Allocate a fresh ID, register it `Visible`, and attach the
    /// given [`Pickable`] geometry so the `ais_select_*` pickers can
    /// ray-cast against it. Returns the new ID.
    pub fn display_geometry(&mut self, geometry: Pickable) -> usize {
        let id = self.display();
        self.geometry.insert(id, geometry);
        id
    }

    /// Look up the current state of `id`.
    pub fn state(&self, id: usize) -> Option<ObjectState> {
        self.states.get(&id).copied()
    }

    /// Borrow the pickable geometry of `id`, if any was registered.
    pub fn geometry(&self, id: usize) -> Option<&Pickable> {
        self.geometry.get(&id)
    }

    /// Force `id` into the given state. Returns the prior state, or
    /// `None` if `id` was never displayed.
    pub fn set_state(&mut self, id: usize, state: ObjectState) -> Option<ObjectState> {
        self.states.insert(id, state)
    }

    /// Install the current camera so picks can unproject. The renderer
    /// calls this every frame with `projection * view` and the
    /// viewport size.
    pub fn set_view(&mut self, view_proj: Matrix4<f64>, width: f32, height: f32) {
        self.view = Some(PickView {
            view_proj,
            width,
            height,
        });
    }

    /// The current [`PickView`], if [`set_view`](Self::set_view) has
    /// been called.
    pub fn view(&self) -> Option<PickView> {
        self.view
    }

    /// Unproject a screen pixel `(x, y)` (origin top-left, +y down)
    /// into a world-space [`Ray`].
    ///
    /// Returns `None` if no camera has been installed. The ray is
    /// built by inverse-transforming the near- and far-plane points
    /// of the pixel through the view-projection matrix — correct for
    /// both perspective and orthographic cameras.
    pub fn ray_at(&self, x: f32, y: f32) -> Option<Ray> {
        let view = self.view?;
        let inv = view.view_proj.try_inverse()?;
        // Pixel → NDC. NDC x in [-1, 1], y in [-1, 1] with +y up, so
        // flip the screen y.
        let ndc_x = (2.0 * x as f64 / view.width.max(1.0) as f64) - 1.0;
        let ndc_y = 1.0 - (2.0 * y as f64 / view.height.max(1.0) as f64);
        // Near point at NDC z = -1 (OpenGL-style) … but nalgebra's
        // perspective uses z in [-1, 1] too; use both ends.
        let unproject = |ndc_z: f64| -> Vector3<f64> {
            let clip = Vector4::new(ndc_x, ndc_y, ndc_z, 1.0);
            let world = inv * clip;
            if world.w.abs() < 1e-12 {
                Vector3::new(world.x, world.y, world.z)
            } else {
                Vector3::new(world.x / world.w, world.y / world.w, world.z / world.w)
            }
        };
        let near = unproject(-1.0);
        let far = unproject(1.0);
        let dir = far - near;
        let len = dir.norm();
        if len < 1e-12 {
            return None;
        }
        Some(Ray {
            origin: near,
            direction: dir / len,
        })
    }

    /// Ray-cast `ray` against every object that carries [`Pickable`]
    /// geometry and is not `Hidden`; return the nearest hit.
    ///
    /// `tolerance` is the pixel-scale proximity used for polyline /
    /// point picks (in world units — the caller derives it from the
    /// camera's units-per-pixel at the geometry depth).
    pub fn pick_nearest(&self, ray: &Ray, tolerance: f64) -> Option<PickHit> {
        let mut best: Option<PickHit> = None;
        for (&id, geom) in &self.geometry {
            if self.states.get(&id) == Some(&ObjectState::Hidden) {
                continue;
            }
            if let Some(t) = ray_pick_geometry(ray, geom, tolerance) {
                if best.is_none_or(|b| t < b.distance) {
                    best = Some(PickHit { id, distance: t });
                }
            }
        }
        best
    }

    /// Ray-cast `ray` against the geometry of `id` only, returning the
    /// nearest sub-element index hit (face / edge / vertex index for a
    /// `Tagged*` geometry, `0` for a plain one), or `None` on a miss.
    ///
    /// This is the subshape-selection primitive — `ais_select_face` /
    /// `ais_select_edge` / `ais_select_vertex` route through it.
    pub fn pick_subshape(&self, id: usize, ray: &Ray, tolerance: f64) -> Option<usize> {
        let geom = self.geometry.get(&id)?;
        ray_pick_subshape(ray, geom, tolerance).map(|(_, sub)| sub)
    }

    /// Project a world-space point to screen pixels `(x, y)` with the
    /// installed camera. Returns `None` if no camera is set or the
    /// point is behind the camera (clip `w <= 0`).
    pub fn project(&self, world: Vector3<f64>) -> Option<(f32, f32)> {
        let view = self.view?;
        let clip = view.view_proj * Vector4::new(world.x, world.y, world.z, 1.0);
        if clip.w <= 1e-9 {
            return None;
        }
        let ndc_x = clip.x / clip.w;
        let ndc_y = clip.y / clip.w;
        let sx = (ndc_x + 1.0) * 0.5 * view.width as f64;
        let sy = (1.0 - ndc_y) * 0.5 * view.height as f64;
        Some((sx as f32, sy as f32))
    }

    /// All IDs currently in the `Selected` state.
    pub fn selected_ids(&self) -> Vec<usize> {
        let mut out: Vec<usize> = self
            .states
            .iter()
            .filter_map(|(id, s)| (*s == ObjectState::Selected).then_some(*id))
            .collect();
        out.sort();
        out
    }

    /// IDs of every object carrying pickable geometry.
    pub fn geometry_ids(&self) -> Vec<usize> {
        let mut out: Vec<usize> = self.geometry.keys().copied().collect();
        out.sort();
        out
    }

    /// Total count of registered objects (regardless of state).
    pub fn len(&self) -> usize {
        self.states.len()
    }

    /// True iff no objects have been displayed yet.
    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }
}

/// Ray-cast against one [`Pickable`], returning the hit `t` (distance
/// along the ray) of the closest intersection, or `None`. Tagged
/// variants are picked as whole objects here — for the sub-element
/// index use [`ray_pick_subshape`].
pub(crate) fn ray_pick_geometry(ray: &Ray, geom: &Pickable, tolerance: f64) -> Option<f64> {
    ray_pick_subshape(ray, geom, tolerance).map(|(t, _)| t)
}

/// Ray-cast against one [`Pickable`], returning `(hit_t, sub_index)`.
///
/// For the tagged variants `sub_index` is the hit face / edge / vertex
/// index. For the plain variants it is `0` (the object has no
/// sub-elements). `None` when the ray misses.
pub(crate) fn ray_pick_subshape(
    ray: &Ray,
    geom: &Pickable,
    tolerance: f64,
) -> Option<(f64, usize)> {
    let v3 = |p: &[f64; 3]| Vector3::new(p[0], p[1], p[2]);
    match geom {
        Pickable::Mesh { triangles } => {
            let mut best: Option<f64> = None;
            for tri in triangles.chunks_exact(3) {
                if let Some(t) = ray_triangle(ray, v3(&tri[0]), v3(&tri[1]), v3(&tri[2])) {
                    if best.is_none_or(|x| t < x) {
                        best = Some(t);
                    }
                }
            }
            best.map(|t| (t, 0))
        }
        Pickable::Polyline { points } => {
            let mut best: Option<f64> = None;
            for seg in points.windows(2) {
                if let Some(t) = ray_segment(ray, v3(&seg[0]), v3(&seg[1]), tolerance) {
                    if best.is_none_or(|x| t < x) {
                        best = Some(t);
                    }
                }
            }
            best.map(|t| (t, 0))
        }
        Pickable::Point { position } => ray_point(ray, v3(position), tolerance).map(|t| (t, 0)),
        Pickable::TaggedMesh {
            triangles,
            face_tags,
        } => {
            let mut best: Option<(f64, usize)> = None;
            for (i, tri) in triangles.chunks_exact(3).enumerate() {
                if let Some(t) = ray_triangle(ray, v3(&tri[0]), v3(&tri[1]), v3(&tri[2])) {
                    let tag = face_tags.get(i).copied().unwrap_or(0);
                    if best.is_none_or(|(x, _)| t < x) {
                        best = Some((t, tag));
                    }
                }
            }
            best
        }
        Pickable::TaggedEdges { edges } => {
            let mut best: Option<(f64, usize)> = None;
            for (edge_idx, pts) in edges {
                for seg in pts.windows(2) {
                    if let Some(t) = ray_segment(ray, v3(&seg[0]), v3(&seg[1]), tolerance) {
                        if best.is_none_or(|(x, _)| t < x) {
                            best = Some((t, *edge_idx));
                        }
                    }
                }
            }
            best
        }
        Pickable::TaggedVertices { vertices } => {
            let mut best: Option<(f64, usize)> = None;
            for (vtx_idx, pos) in vertices {
                if let Some(t) = ray_point(ray, v3(pos), tolerance) {
                    if best.is_none_or(|(x, _)| t < x) {
                        best = Some((t, *vtx_idx));
                    }
                }
            }
            best
        }
    }
}

/// Möller-Trumbore ray-triangle intersection — returns the positive
/// hit distance, or `None`.
fn ray_triangle(ray: &Ray, a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>) -> Option<f64> {
    let e1 = b - a;
    let e2 = c - a;
    let pvec = ray.direction.cross(&e2);
    let det = e1.dot(&pvec);
    if det.abs() < 1e-12 {
        return None;
    }
    let inv = 1.0 / det;
    let tvec = ray.origin - a;
    let u = tvec.dot(&pvec) * inv;
    if !(-1e-9..=1.0 + 1e-9).contains(&u) {
        return None;
    }
    let qvec = tvec.cross(&e1);
    let v = ray.direction.dot(&qvec) * inv;
    if v < -1e-9 || u + v > 1.0 + 1e-9 {
        return None;
    }
    let t = e2.dot(&qvec) * inv;
    if t > 1e-9 {
        Some(t)
    } else {
        None
    }
}

/// Closest approach of a ray to a segment — a hit if the minimum
/// distance is within `tolerance`. Returns the ray `t` at closest
/// approach.
fn ray_segment(ray: &Ray, a: Vector3<f64>, b: Vector3<f64>, tolerance: f64) -> Option<f64> {
    // Standard segment-segment / line-line closest-point solve.
    let d1 = ray.direction;
    let d2 = b - a;
    let r = ray.origin - a;
    let a11 = d1.dot(&d1);
    let a12 = -d1.dot(&d2);
    let a22 = d2.dot(&d2);
    let b1 = -d1.dot(&r);
    let b2 = d2.dot(&r);
    let denom = a11 * a22 - a12 * a12;
    if denom.abs() < 1e-12 {
        return None; // parallel
    }
    let t_ray = (b1 * a22 - a12 * b2) / denom;
    let mut s_seg = (a11 * b2 - b1 * a12) / denom;
    if t_ray < 0.0 {
        return None; // behind the camera
    }
    s_seg = s_seg.clamp(0.0, 1.0);
    let p_ray = ray.origin + d1 * t_ray;
    let p_seg = a + d2 * s_seg;
    if (p_ray - p_seg).norm() <= tolerance {
        Some(t_ray)
    } else {
        None
    }
}

/// Closest approach of a ray to a point — a hit if within `tolerance`.
fn ray_point(ray: &Ray, p: Vector3<f64>, tolerance: f64) -> Option<f64> {
    let w = p - ray.origin;
    let t = w.dot(&ray.direction);
    if t < 0.0 {
        return None;
    }
    let closest = ray.origin + ray.direction * t;
    if (closest - p).norm() <= tolerance {
        Some(t)
    } else {
        None
    }
}

/// Construct a fresh, empty interactive context.
///
/// This op cannot fail — wraps [`InteractiveContext::default`] for API
/// consistency.
pub fn ais_interactive_context() -> Result<InteractiveContext, OcctVizError> {
    Ok(InteractiveContext::default())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A view-projection that looks down -Z at the origin from z=+10,
    /// with a simple perspective. Good enough to unproject a ray.
    fn test_view() -> Matrix4<f64> {
        let eye = Vector3::new(0.0, 0.0, 10.0);
        let view = Matrix4::look_at_rh(
            &nalgebra::Point3::from(eye),
            &nalgebra::Point3::origin(),
            &Vector3::y(),
        );
        let proj = Matrix4::new_perspective(1.0, 60f64.to_radians(), 0.1, 100.0);
        proj * view
    }

    /// A unit-ish triangle facing the camera, centred on the axis.
    fn axis_triangle() -> Pickable {
        Pickable::Mesh {
            triangles: vec![[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, 0.0]],
        }
    }

    #[test]
    fn fresh_context_is_empty() {
        let ctx = ais_interactive_context().unwrap();
        assert_eq!(ctx.len(), 0);
        assert!(ctx.is_empty());
        assert!(ctx.selected_ids().is_empty());
        assert!(ctx.geometry_ids().is_empty());
    }

    #[test]
    fn display_allocates_monotonic_ids() {
        let mut ctx = ais_interactive_context().unwrap();
        let a = ctx.display();
        let b = ctx.display();
        let c = ctx.display();
        assert_eq!(a, 0);
        assert_eq!(b, 1);
        assert_eq!(c, 2);
        assert_eq!(ctx.len(), 3);
    }

    #[test]
    fn display_geometry_attaches_pickable() {
        let mut ctx = ais_interactive_context().unwrap();
        let id = ctx.display_geometry(axis_triangle());
        assert!(matches!(ctx.geometry(id), Some(Pickable::Mesh { .. })));
        assert_eq!(ctx.geometry_ids(), vec![id]);
    }

    #[test]
    fn ray_at_needs_a_view() {
        let ctx = ais_interactive_context().unwrap();
        assert!(ctx.ray_at(50.0, 50.0).is_none());
    }

    #[test]
    fn center_ray_points_into_the_scene() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        // A ray through the screen centre starts near the eye (z≈+10)
        // and points toward -Z.
        let ray = ctx.ray_at(50.0, 50.0).unwrap();
        assert!(ray.origin.z > 0.0, "origin near the eye: {:?}", ray.origin);
        assert!(
            ray.direction.z < -0.5,
            "centre ray points into -Z: {:?}",
            ray.direction
        );
    }

    #[test]
    fn pick_nearest_hits_the_triangle_under_the_cursor() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(axis_triangle());
        // The centre ray hits the axis triangle at z=0.
        let ray = ctx.ray_at(50.0, 50.0).unwrap();
        let hit = ctx
            .pick_nearest(&ray, 0.1)
            .expect("centre ray hits triangle");
        assert_eq!(hit.id, id);
        // The hit is ~10 units away (eye at z=10, triangle at z=0).
        assert!((hit.distance - 10.0).abs() < 1.0, "t={}", hit.distance);
    }

    #[test]
    fn pick_nearest_returns_closest_of_two() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        // Far triangle at z=0, near triangle at z=5 — the near one wins.
        let _far = ctx.display_geometry(axis_triangle());
        let near = ctx.display_geometry(Pickable::Mesh {
            triangles: vec![[-1.0, -1.0, 5.0], [1.0, -1.0, 5.0], [0.0, 1.0, 5.0]],
        });
        let ray = ctx.ray_at(50.0, 50.0).unwrap();
        let hit = ctx.pick_nearest(&ray, 0.1).unwrap();
        assert_eq!(hit.id, near, "nearer triangle should win");
    }

    #[test]
    fn hidden_objects_are_not_picked() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        let id = ctx.display_geometry(axis_triangle());
        ctx.set_state(id, ObjectState::Hidden);
        let ray = ctx.ray_at(50.0, 50.0).unwrap();
        assert!(ctx.pick_nearest(&ray, 0.1).is_none());
    }

    #[test]
    fn ray_misses_an_off_axis_triangle() {
        let mut ctx = ais_interactive_context().unwrap();
        ctx.set_view(test_view(), 100.0, 100.0);
        // Triangle far off to the side — the centre ray misses it.
        ctx.display_geometry(Pickable::Mesh {
            triangles: vec![[50.0, 50.0, 0.0], [52.0, 50.0, 0.0], [51.0, 52.0, 0.0]],
        });
        let ray = ctx.ray_at(50.0, 50.0).unwrap();
        assert!(ctx.pick_nearest(&ray, 0.1).is_none());
    }
}
