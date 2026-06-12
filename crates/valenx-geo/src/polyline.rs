//! Procedural polyline + linear-extrusion primitives.
//!
//! ## Scope
//!
//! v0 of the "generate simple geometry without needing FreeCAD"
//! pipeline. Lets users programmatically build:
//!
//! - **`Polyline`** — an ordered list of 3D vertices, optionally
//!   marked as closed (last vertex implicitly connects back to the
//!   first).
//! - **`linear_extrude(polyline, direction)`** — sweeps a closed
//!   polyline along a direction vector, producing a Mesh of
//!   triangulated side walls + top / bottom caps. The caps are
//!   triangulated via fan triangulation from vertex 0 — robust for
//!   convex polygons, fragile for concave / self-intersecting ones
//!   (we'll do ear-clipping in a follow-up).
//!
//! ## Why
//!
//! Sweep workflows often need a tiny synthetic geometry for
//! validation (a cube, a wedge, a small box around an inlet), and
//! today the only path is via FreeCAD's Python script. This
//! module gives Rust callers a one-line way to build those
//! geometries directly into `valenx_mesh::Mesh`, ready to feed
//! the meshing / CFD pipeline.
//!
//! ## What's NOT here yet
//!
//! - Ear-clipping or constrained Delaunay for concave polygons
//! - Sweep along a curve (only linear extrude)
//! - Boolean ops between procedural shapes (use
//!   `valenx_mesh::boolean::union_concatenate` after meshing each
//!   piece)
//! - NURBS / B-Spline curves

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use valenx_mesh::{ElementBlock, ElementType, Mesh};

/// An ordered list of 3D vertices forming a polyline. `closed`
/// indicates the last vertex implicitly connects back to the first
/// — required for [`linear_extrude`] so the extrusion has well-
/// defined side walls + caps.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Polyline {
    pub vertices: Vec<Vector3<f64>>,
    pub closed: bool,
}

impl Polyline {
    /// Convenience: a closed quadrilateral in the XY plane.
    /// Vertices are (0,0), (w,0), (w,h), (0,h) at z=0.
    pub fn rectangle_xy(width: f64, height: f64) -> Self {
        Self {
            vertices: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(width, 0.0, 0.0),
                Vector3::new(width, height, 0.0),
                Vector3::new(0.0, height, 0.0),
            ],
            closed: true,
        }
    }

    /// Convenience: a closed regular polygon centred at the origin
    /// with `n` sides and circumradius `r`, lying in the XY plane.
    /// Useful for cylinders / hexagonal prisms when paired with
    /// [`linear_extrude`].
    ///
    /// Panics if `n < 3` or `n > MAX_REGULAR_POLYGON_SIDES`. The
    /// fallible [`try_regular_polygon_xy`](Self::try_regular_polygon_xy)
    /// returns a `Result` for hostile inputs (file-loaded sweep
    /// configs that haven't been bounds-checked yet).
    pub fn regular_polygon_xy(n: usize, r: f64) -> Self {
        Self::try_regular_polygon_xy(n, r).expect("valid n + r")
    }

    /// Fallible variant of [`regular_polygon_xy`](Self::regular_polygon_xy).
    /// Round-4 DoS hardening: caps `n` at
    /// [`MAX_REGULAR_POLYGON_SIDES`] so a hostile input can't drive
    /// `Vec::with_capacity(usize::MAX)` to OOM the host.
    pub fn try_regular_polygon_xy(n: usize, r: f64) -> Result<Self, PolygonError> {
        if n < 3 {
            return Err(PolygonError::TooFewSides { n });
        }
        if n > MAX_REGULAR_POLYGON_SIDES {
            return Err(PolygonError::TooManySides {
                n,
                max: MAX_REGULAR_POLYGON_SIDES,
            });
        }
        let mut v = Vec::with_capacity(n);
        for i in 0..n {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            v.push(Vector3::new(r * theta.cos(), r * theta.sin(), 0.0));
        }
        Ok(Self {
            vertices: v,
            closed: true,
        })
    }

    /// Total length of the polyline — the sum of consecutive segment lengths,
    /// plus the closing edge (last → first) when [`closed`](Self::closed) is set.
    /// Fewer than two vertices yields `0.0`. A length, distinct from a polygon area.
    pub fn length(&self) -> f64 {
        let n = self.vertices.len();
        if n < 2 {
            return 0.0;
        }
        let mut total: f64 = self.vertices.windows(2).map(|w| (w[1] - w[0]).norm()).sum();
        if self.closed {
            total += (self.vertices[0] - self.vertices[n - 1]).norm();
        }
        total
    }
}

/// Round-4 DoS hardening: upper bound on the side count of a regular
/// polygon. 100,000 is generous (any smooth-curve approximation that
/// needs more sides than this should use a NURBS curve, not a polygon)
/// and small enough to bound the worst-case allocation at ~2 MB.
pub const MAX_REGULAR_POLYGON_SIDES: usize = 100_000;

/// Errors raised by polygon constructors. Surfaces the DoS guard's
/// "too many sides" path as a structured error rather than a panic.
#[derive(Debug, Error)]
pub enum PolygonError {
    /// `n < 3` — a polygon needs at least 3 vertices to enclose area.
    #[error("regular polygon needs at least 3 sides; got {n}")]
    TooFewSides {
        /// The side count that was rejected.
        n: usize,
    },
    /// `n > MAX_REGULAR_POLYGON_SIDES` — DoS guard.
    #[error("regular polygon side count {n} exceeds the safety cap of {max}")]
    TooManySides {
        /// The side count that was rejected.
        n: usize,
        /// The cap that triggered the rejection.
        max: usize,
    },
}

/// Errors raised by [`linear_extrude`].
#[derive(Debug, Error)]
pub enum ExtrudeError {
    /// First and last vertices of the polyline are not coincident.
    #[error("polyline must be closed for extrusion")]
    NotClosed,
    /// Fewer than three distinct vertices — extruding wouldn't yield
    /// a closed surface.
    #[error("polyline must have at least 3 vertices to extrude")]
    TooFewVertices,
    /// The extrusion direction vector has zero length.
    #[error("extrusion direction must have non-zero length")]
    ZeroDirection,
}

/// Linearly extrude a closed polyline along `direction`. The output
/// is a triangle mesh ([`ElementType::Tri3`]) with:
///
/// - one side-wall quad per polyline edge (split into 2 triangles)
/// - one top cap triangulated via fan from vertex 0
/// - one bottom cap (same fan, reversed winding)
///
/// Total triangles: `2*N + 2*(N-2)` where N is the polyline vertex
/// count. For a triangle (N=3) that's 8 triangles; for a hex (N=6)
/// that's 20.
///
/// `id` becomes the output Mesh's id so audit logs / file dumps
/// retain a hint of the source.
pub fn linear_extrude(
    polyline: &Polyline,
    direction: Vector3<f64>,
    id: impl Into<String>,
) -> Result<Mesh, ExtrudeError> {
    if !polyline.closed {
        return Err(ExtrudeError::NotClosed);
    }
    let n = polyline.vertices.len();
    if n < 3 {
        return Err(ExtrudeError::TooFewVertices);
    }
    if direction.norm() < 1e-30 {
        return Err(ExtrudeError::ZeroDirection);
    }

    let mut mesh = Mesh::new(id);

    // Bottom and top vertex rings, in declaration order. Bottom
    // copies the polyline directly; top translates by `direction`.
    for v in &polyline.vertices {
        mesh.nodes.push(*v);
    }
    for v in &polyline.vertices {
        mesh.nodes.push(*v + direction);
    }

    let mut tris = ElementBlock::new(ElementType::Tri3);

    // Side walls: for each edge i -> (i+1)%n, emit two triangles
    // sharing the diagonal from bottom[i] to top[(i+1)%n].
    let bot = |i: usize| -> u32 { i as u32 };
    let top = |i: usize| -> u32 { (n + i) as u32 };
    for i in 0..n {
        let j = (i + 1) % n;
        // CCW winding from outside, assuming `direction` points
        // "up" relative to the polyline's standard orientation.
        // First triangle: bot[i], bot[j], top[j]
        tris.connectivity
            .extend_from_slice(&[bot(i), bot(j), top(j)]);
        // Second triangle: bot[i], top[j], top[i]
        tris.connectivity
            .extend_from_slice(&[bot(i), top(j), top(i)]);
    }

    // Bottom cap: fan triangulation from vertex 0 (CW from below
    // to keep the outward normal pointing down).
    for i in 1..(n - 1) {
        tris.connectivity
            .extend_from_slice(&[bot(0), bot(i + 1), bot(i)]);
    }

    // Top cap: same fan, reverse winding so the normal points up.
    for i in 1..(n - 1) {
        tris.connectivity
            .extend_from_slice(&[top(0), top(i), top(i + 1)]);
    }

    mesh.element_blocks.push(tris);
    mesh.recompute_stats();
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rectangle_xy_has_four_vertices_closed() {
        let r = Polyline::rectangle_xy(2.0, 3.0);
        assert_eq!(r.vertices.len(), 4);
        assert!(r.closed);
        assert_eq!(r.vertices[0], Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(r.vertices[2], Vector3::new(2.0, 3.0, 0.0));
    }

    #[test]
    fn regular_polygon_has_n_vertices_closed() {
        let p = Polyline::regular_polygon_xy(6, 1.0);
        assert_eq!(p.vertices.len(), 6);
        assert!(p.closed);
        // First vertex sits at (r, 0, 0).
        assert!((p.vertices[0].x - 1.0).abs() < 1e-12);
        assert!(p.vertices[0].y.abs() < 1e-12);
    }

    #[test]
    fn length_sums_segments_and_closing_edge() {
        // Closed 4×3 rectangle → perimeter 2·(4+3) = 14.
        let r = Polyline::rectangle_xy(4.0, 3.0);
        assert!((r.length() - 14.0).abs() < 1e-12);
        // Open 3-then-4 path (no closing edge) → 3 + 4 = 7.
        let verts = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(3.0, 0.0, 0.0),
            Vector3::new(3.0, 4.0, 0.0),
        ];
        let open = Polyline {
            vertices: verts.clone(),
            closed: false,
        };
        assert!((open.length() - 7.0).abs() < 1e-12);
        // The same vertices closed add the 3-4-5 closing edge → 7 + 5 = 12.
        let closed = Polyline {
            vertices: verts,
            closed: true,
        };
        assert!((closed.length() - 12.0).abs() < 1e-12);
        // Fewer than two vertices → 0.
        let single = Polyline {
            vertices: vec![Vector3::zeros()],
            closed: false,
        };
        assert_eq!(single.length(), 0.0);
    }

    #[test]
    fn extrude_rectangle_yields_box_mesh() {
        // 4-vertex polyline + extrude -> 8 nodes (4 bottom + 4 top)
        // and 12 triangles (4 quads x 2 + 2 caps x 2) = wait,
        // bottom cap fan = N-2 = 2; top cap = 2. Side walls = N x 2 = 8.
        // Total = 8 + 4 = 12.
        let r = Polyline::rectangle_xy(1.0, 1.0);
        let mesh = linear_extrude(&r, Vector3::new(0.0, 0.0, 0.5), "box").expect("extrude");
        assert_eq!(mesh.id, "box");
        assert_eq!(mesh.nodes.len(), 8);
        assert_eq!(mesh.element_blocks.len(), 1);
        assert_eq!(mesh.element_blocks[0].element_type, ElementType::Tri3);
        // 12 triangles * 3 indices each = 36 connectivity entries.
        assert_eq!(mesh.element_blocks[0].connectivity.len(), 36);
        // Top vertices sit at z = 0.5.
        for top_i in 4..8 {
            assert!((mesh.nodes[top_i].z - 0.5).abs() < 1e-12);
        }
    }

    #[test]
    fn extrude_hexagon_produces_expected_triangle_count() {
        // 6-vertex hex: side walls = 12 tris; caps = 2 * (6-2) = 8 tris.
        let p = Polyline::regular_polygon_xy(6, 1.0);
        let mesh = linear_extrude(&p, Vector3::new(0.0, 0.0, 1.0), "hex").expect("extrude");
        assert_eq!(mesh.nodes.len(), 12);
        assert_eq!(mesh.element_blocks[0].connectivity.len() / 3, 20);
    }

    #[test]
    fn extrude_rejects_open_polyline() {
        let p = Polyline {
            vertices: vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(1.0, 1.0, 0.0),
            ],
            closed: false,
        };
        let err = linear_extrude(&p, Vector3::new(0.0, 0.0, 1.0), "x").unwrap_err();
        assert!(matches!(err, ExtrudeError::NotClosed));
    }

    #[test]
    fn extrude_rejects_polyline_with_too_few_vertices() {
        let p = Polyline {
            vertices: vec![Vector3::zeros(), Vector3::new(1.0, 0.0, 0.0)],
            closed: true,
        };
        let err = linear_extrude(&p, Vector3::new(0.0, 0.0, 1.0), "x").unwrap_err();
        assert!(matches!(err, ExtrudeError::TooFewVertices));
    }

    #[test]
    fn extrude_rejects_zero_length_direction() {
        let r = Polyline::rectangle_xy(1.0, 1.0);
        let err = linear_extrude(&r, Vector3::zeros(), "x").unwrap_err();
        assert!(matches!(err, ExtrudeError::ZeroDirection));
    }

    #[test]
    fn extrude_recomputes_stats() {
        let r = Polyline::rectangle_xy(1.0, 1.0);
        let mesh = linear_extrude(&r, Vector3::new(0.0, 0.0, 0.5), "box").expect("extrude");
        assert_eq!(mesh.stats.node_count, 8);
        assert_eq!(mesh.stats.element_count, 12);
    }

    #[test]
    fn side_wall_connectivity_uses_bottom_and_top_vertices() {
        // Smoke test: every connectivity index lies in [0, 2N).
        let p = Polyline::regular_polygon_xy(5, 1.0);
        let mesh = linear_extrude(&p, Vector3::new(0.0, 0.0, 1.0), "pent").expect("extrude");
        let max_idx = mesh.element_blocks[0]
            .connectivity
            .iter()
            .copied()
            .max()
            .unwrap();
        assert!(max_idx < (2 * 5) as u32, "max idx {max_idx} >= 2N");
    }
}
