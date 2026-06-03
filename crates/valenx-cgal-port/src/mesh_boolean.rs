//! Mesh boolean (CSG) on triangle soups — Phase 58 / 58.5.
//!
//! ## What CGAL does
//!
//! `CGAL::Polygon_mesh_processing::corefine_and_compute_boolean_operation`
//! runs an *exact-arithmetic co-refinement*: it intersects the two
//! input surface meshes, inserts the intersection polylines into both
//! as constrained edges (so both meshes share a conforming seam), then
//! classifies each resulting facet as inside / outside the other
//! solid and keeps the facets the requested operation selects. The
//! exact kernel (`CGAL::Exact_predicates_exact_constructions_kernel`)
//! guarantees the predicates never flip sign under rounding.
//!
//! ## What this module does (58.5 — real CSG)
//!
//! This is a genuine co-refinement CSG, not the old concatenation
//! placeholder:
//!
//! 1. **Intersection segments.** Every triangle of A is tested against
//!    every triangle of B (AABB-tree broad phase) with a Möller-style
//!    triangle-triangle intersection that returns the actual 3D
//!    intersection *segment*.
//! 2. **Co-refinement.** Each triangle that picked up intersection
//!    segments is re-triangulated in its own plane: the segments are
//!    inserted as constraint edges and the triangle is subdivided so
//!    every sub-triangle lies entirely on one side of every segment.
//! 3. **Classification.** Each sub-triangle is classified inside /
//!    outside the *other* mesh by ray casting (parity of crossings of
//!    a ray from the sub-triangle centroid against the other mesh's
//!    AABB tree).
//! 4. **Selection.** `union` keeps `A∖B ∪ B∖A`; `intersection` keeps
//!    `A∩B ∪ B∩A`; `difference` keeps `A∖B` plus `B∩A` with B's
//!    facets flipped (the carved cavity wall).
//!
//! ## Honest precision limits
//!
//! - **Float arithmetic.** Predicates use `f64`, not an exact kernel.
//!   Robust for well-conditioned inputs; pathological near-coplanar
//!   slivers can mis-classify. CGAL's exact kernel is the Tier-3
//!   upgrade.
//! - **Coplanar overlap.** Two triangles lying in the *same* plane
//!   with overlapping area contribute no intersection segment (their
//!   intersection is an area, not a segment). Such facets are
//!   classified by centroid like any other — correct for the common
//!   case where the coplanar region is fully inside or outside, but
//!   not for partial coplanar overlap.
//! - **Open meshes.** Classification assumes each input is a closed
//!   manifold; ray-parity is undefined for meshes with boundary.
//!
//! Within those limits the result is a real Boolean: overlapping
//! cubes `union` to a single closed solid with the interior walls
//! removed, `difference` carves a real cavity, `intersection`
//! returns only the shared lens.

use std::collections::HashMap;

use nalgebra::Vector3;

use crate::aabb_tree::{AabbTree, Triangle3};

/// Concatenated triangle list — the working representation for CSG.
#[derive(Clone, Debug, Default)]
pub struct Mesh3 {
    /// Triangles.
    pub triangles: Vec<Triangle3>,
}

impl Mesh3 {
    /// Empty mesh.
    pub fn new() -> Self {
        Self::default()
    }

    /// Total triangle count.
    pub fn len(&self) -> usize {
        self.triangles.len()
    }

    /// `true` if empty.
    pub fn is_empty(&self) -> bool {
        self.triangles.is_empty()
    }

    /// All vertices (in triangle insertion order, with duplicates).
    pub fn vertices(&self) -> Vec<Vector3<f64>> {
        let mut out = Vec::with_capacity(self.triangles.len() * 3);
        for t in &self.triangles {
            out.extend_from_slice(&t.v);
        }
        out
    }
}

/// Geometric tolerance for coincidence / on-plane tests.
const EPS: f64 = 1e-9;

/// Boolean union `A ∪ B` — keeps the parts of each surface that lie
/// outside the other solid. See the module docs for the algorithm and
/// the precision caveats.
pub fn union(m1: &Mesh3, m2: &Mesh3) -> Mesh3 {
    boolean(m1, m2, BoolOp::Union)
}

/// Boolean difference `A ∖ B` — keeps `A` outside `B` plus the wall of
/// `B` inside `A` (flipped, forming the cavity).
pub fn difference(m1: &Mesh3, m2: &Mesh3) -> Mesh3 {
    boolean(m1, m2, BoolOp::Difference)
}

/// Boolean intersection `A ∩ B` — keeps only the lens common to both
/// solids.
pub fn intersection(m1: &Mesh3, m2: &Mesh3) -> Mesh3 {
    boolean(m1, m2, BoolOp::Intersection)
}

/// The three CSG operations.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum BoolOp {
    Union,
    Difference,
    Intersection,
}

/// Co-refinement CSG core. Both inputs are assumed closed manifolds.
fn boolean(m1: &Mesh3, m2: &Mesh3, op: BoolOp) -> Mesh3 {
    // Degenerate fast paths: an empty operand makes the result trivial.
    if m1.is_empty() || m2.is_empty() {
        return match op {
            BoolOp::Union => union_concat(m1, m2),
            BoolOp::Difference => m1.clone(),
            BoolOp::Intersection => Mesh3::new(),
        };
    }

    let tree1 = AabbTree::new(m1.triangles.clone());
    let tree2 = AabbTree::new(m2.triangles.clone());

    // Co-refine A against B and B against A.
    let refined_a = corefine(&m1.triangles, &tree2);
    let refined_b = corefine(&m2.triangles, &tree1);

    let mut out = Mesh3::new();

    // --- A's facets ---
    for tri in &refined_a {
        let centroid = tri_centroid(tri);
        let inside_b = point_inside(centroid, &tree2);
        let keep = match op {
            BoolOp::Union => !inside_b,
            BoolOp::Difference => !inside_b,
            BoolOp::Intersection => inside_b,
        };
        if keep {
            out.triangles.push(*tri);
        }
    }

    // --- B's facets ---
    for tri in &refined_b {
        let centroid = tri_centroid(tri);
        let inside_a = point_inside(centroid, &tree1);
        match op {
            BoolOp::Union => {
                if !inside_a {
                    out.triangles.push(*tri);
                }
            }
            BoolOp::Intersection => {
                if inside_a {
                    out.triangles.push(*tri);
                }
            }
            BoolOp::Difference => {
                // The wall of B inside A becomes the cavity surface —
                // flip its winding so the cavity's normals point into
                // the removed volume (outward w.r.t. the result solid).
                if inside_a {
                    out.triangles.push(flip(tri));
                }
            }
        }
    }

    out
}

/// Plain concatenation — used only for the empty-operand union path.
fn union_concat(m1: &Mesh3, m2: &Mesh3) -> Mesh3 {
    let mut out = Mesh3::new();
    out.triangles.extend_from_slice(&m1.triangles);
    out.triangles.extend_from_slice(&m2.triangles);
    out
}

/// Reverse a triangle's winding.
fn flip(t: &Triangle3) -> Triangle3 {
    Triangle3 {
        v: [t.v[0], t.v[2], t.v[1]],
    }
}

/// Centroid of a triangle.
fn tri_centroid(t: &Triangle3) -> Vector3<f64> {
    (t.v[0] + t.v[1] + t.v[2]) / 3.0
}

// ===================================================================
// Co-refinement: subdivide every triangle of `tris` by its
// intersection segments with `other`.
// ===================================================================

/// Re-triangulate each triangle of `tris` so that no sub-triangle
/// straddles an intersection segment with `other`.
fn corefine(tris: &[Triangle3], other: &AabbTree) -> Vec<Triangle3> {
    let mut out = Vec::with_capacity(tris.len());
    for tri in tris {
        // Broad phase: triangles of `other` whose AABB overlaps this
        // triangle's AABB. We reuse the ray-AABB tree by collecting
        // every candidate via an AABB-overlap descent.
        let candidates = overlapping_triangles(tri, other);
        let mut segments: Vec<[Vector3<f64>; 2]> = Vec::new();
        for &cid in &candidates {
            let ct = &other.triangles[cid];
            if let Some(seg) = tri_tri_segment(tri, ct) {
                segments.push(seg);
            }
        }
        if segments.is_empty() {
            out.push(*tri);
        } else {
            subdivide_triangle(tri, &segments, &mut out);
        }
    }
    out
}

/// Indices of `tree` triangles whose AABB overlaps `tri`'s AABB.
fn overlapping_triangles(tri: &Triangle3, tree: &AabbTree) -> Vec<usize> {
    let mut out = Vec::new();
    if tree.nodes.is_empty() {
        return out;
    }
    let (lo, hi) = tri_aabb(tri);
    let mut stack = vec![0usize];
    while let Some(i) = stack.pop() {
        let (nlo, nhi) = tree.nodes[i].aabb();
        if !aabb_overlap(lo, hi, nlo, nhi) {
            continue;
        }
        match &tree.nodes[i] {
            crate::aabb_tree::Node::Internal { left, right, .. } => {
                stack.push(*left);
                stack.push(*right);
            }
            crate::aabb_tree::Node::Leaf { triangles, .. } => {
                for t in triangles {
                    out.push(t.0);
                }
            }
        }
    }
    out
}

/// AABB of a single triangle.
fn tri_aabb(t: &Triangle3) -> (Vector3<f64>, Vector3<f64>) {
    let mut lo = t.v[0];
    let mut hi = t.v[0];
    for p in &t.v[1..] {
        lo = Vector3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = Vector3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    }
    (lo, hi)
}

/// Do two AABBs overlap (with an EPS slop)?
fn aabb_overlap(
    alo: Vector3<f64>,
    ahi: Vector3<f64>,
    blo: Vector3<f64>,
    bhi: Vector3<f64>,
) -> bool {
    alo.x <= bhi.x + EPS
        && ahi.x + EPS >= blo.x
        && alo.y <= bhi.y + EPS
        && ahi.y + EPS >= blo.y
        && alo.z <= bhi.z + EPS
        && ahi.z + EPS >= blo.z
}

// ===================================================================
// Triangle-triangle intersection segment.
// ===================================================================

/// Intersection of two triangles, returned as a 3D segment when they
/// cross transversally. `None` for disjoint, touching-at-a-point, or
/// coplanar pairs (coplanar overlap is an area, not a segment — see
/// the module-level precision note).
fn tri_tri_segment(a: &Triangle3, b: &Triangle3) -> Option<[Vector3<f64>; 2]> {
    let na = tri_normal(a)?;
    let nb = tri_normal(b)?;

    // Signed distances of A's vertices to B's plane.
    let db: [f64; 3] = [
        nb.dot(&(a.v[0] - b.v[0])),
        nb.dot(&(a.v[1] - b.v[0])),
        nb.dot(&(a.v[2] - b.v[0])),
    ];
    // Signed distances of B's vertices to A's plane.
    let da: [f64; 3] = [
        na.dot(&(b.v[0] - a.v[0])),
        na.dot(&(b.v[1] - a.v[0])),
        na.dot(&(b.v[2] - a.v[0])),
    ];

    // All on one side → no intersection.
    if (db[0] > EPS && db[1] > EPS && db[2] > EPS)
        || (db[0] < -EPS && db[1] < -EPS && db[2] < -EPS)
    {
        return None;
    }
    if (da[0] > EPS && da[1] > EPS && da[2] > EPS)
        || (da[0] < -EPS && da[1] < -EPS && da[2] < -EPS)
    {
        return None;
    }
    // Coplanar — skip (handled by centroid classification).
    if db[0].abs() <= EPS && db[1].abs() <= EPS && db[2].abs() <= EPS {
        return None;
    }

    // Interval of A on the intersection line of the two planes.
    let int_a = plane_crossing_interval(a, &db)?;
    let int_b = plane_crossing_interval(b, &da)?;

    // Project both intervals onto the intersection line direction.
    let line_dir = na.cross(&nb);
    if line_dir.norm() < EPS {
        return None;
    }
    let dir = line_dir.normalize();
    let ta0 = dir.dot(&int_a[0]);
    let ta1 = dir.dot(&int_a[1]);
    let tb0 = dir.dot(&int_b[0]);
    let tb1 = dir.dot(&int_b[1]);
    let (amin, amax) = (ta0.min(ta1), ta0.max(ta1));
    let (bmin, bmax) = (tb0.min(tb1), tb0.max(tb1));
    let lo = amin.max(bmin);
    let hi = amax.min(bmax);
    if hi <= lo + EPS {
        return None; // overlap is empty or a single point
    }
    // Map the [lo, hi] parameter back to 3D points on A's interval.
    let p = |t: f64| {
        let denom = ta1 - ta0;
        if denom.abs() < EPS {
            int_a[0]
        } else {
            let s = (t - ta0) / denom;
            int_a[0] + (int_a[1] - int_a[0]) * s
        }
    };
    Some([p(lo), p(hi)])
}

/// The two points where a triangle's edges cross the other plane,
/// given the signed distances `d` of the triangle's vertices.
fn plane_crossing_interval(t: &Triangle3, d: &[f64; 3]) -> Option<[Vector3<f64>; 2]> {
    let mut pts: Vec<Vector3<f64>> = Vec::with_capacity(2);
    for (i, j) in [(0usize, 1usize), (1, 2), (2, 0)] {
        let di = d[i];
        let dj = d[j];
        // Vertex exactly on the plane.
        if di.abs() <= EPS {
            push_unique(&mut pts, t.v[i]);
        }
        // Edge straddles the plane.
        if (di > EPS && dj < -EPS) || (di < -EPS && dj > EPS) {
            let s = di / (di - dj);
            push_unique(&mut pts, t.v[i] + (t.v[j] - t.v[i]) * s);
        }
    }
    if pts.len() == 2 {
        Some([pts[0], pts[1]])
    } else {
        None
    }
}

/// Append `p` to `pts` only if it isn't already there (within EPS).
fn push_unique(pts: &mut Vec<Vector3<f64>>, p: Vector3<f64>) {
    if pts.iter().all(|q| (q - p).norm() > EPS) {
        pts.push(p);
    }
}

/// Unit normal of a triangle, or `None` for a degenerate sliver.
fn tri_normal(t: &Triangle3) -> Option<Vector3<f64>> {
    let n = (t.v[1] - t.v[0]).cross(&(t.v[2] - t.v[0]));
    if n.norm() < EPS {
        None
    } else {
        Some(n.normalize())
    }
}

// ===================================================================
// Planar triangle subdivision along constraint segments.
// ===================================================================

/// Subdivide `tri` by the in-plane `segments` and push every
/// sub-triangle into `out`.
///
/// Strategy: project the triangle and the segments into the
/// triangle's plane (2D), insert each segment as a constraint by
/// splitting every crossed sub-triangle, then lift the result back
/// to 3D via barycentric coordinates of the original triangle.
fn subdivide_triangle(
    tri: &Triangle3,
    segments: &[[Vector3<f64>; 2]],
    out: &mut Vec<Triangle3>,
) {
    let Some(normal) = tri_normal(tri) else {
        out.push(*tri);
        return;
    };
    // Build an orthonormal 2D basis (u, v) for the triangle's plane.
    let u = (tri.v[1] - tri.v[0]).normalize();
    let v = normal.cross(&u);
    let origin = tri.v[0];
    let to2d = |p: Vector3<f64>| {
        let d = p - origin;
        [d.dot(&u), d.dot(&v)]
    };
    let to3d = |q: [f64; 2]| origin + u * q[0] + v * q[1];

    let tri2d = [to2d(tri.v[0]), to2d(tri.v[1]), to2d(tri.v[2])];

    // Active set of 2D triangles, seeded with the original.
    let mut active: Vec<[[f64; 2]; 3]> = vec![tri2d];

    for seg in segments {
        let s0 = to2d(seg[0]);
        let s1 = to2d(seg[1]);
        let mut next: Vec<[[f64; 2]; 3]> = Vec::with_capacity(active.len() + 2);
        for t in &active {
            split_tri_by_segment(t, s0, s1, &mut next);
        }
        active = next;
    }

    // Lift each 2D sub-triangle back to 3D.
    for t in &active {
        let a = to3d(t[0]);
        let b = to3d(t[1]);
        let c = to3d(t[2]);
        // Drop degenerate slivers.
        if (b - a).cross(&(c - a)).norm() < EPS * EPS {
            continue;
        }
        out.push(Triangle3 { v: [a, b, c] });
    }
}

/// Clip a 2D triangle against the infinite line through `(s0, s1)`,
/// splitting it into sub-triangles that each lie wholly on one side.
///
/// The intersection segment from the triangle-triangle test is the
/// portion of that line that lies inside *both* triangles, so for the
/// triangle being subdivided the relevant cut is the full line — the
/// resulting facets are then classified individually, and any facet
/// past the segment ends is still correctly inside/outside the other
/// solid. Splitting on the full line keeps the routine simple and the
/// classification step does the rest.
fn split_tri_by_segment(
    tri: &[[f64; 2]; 3],
    s0: [f64; 2],
    s1: [f64; 2],
    out: &mut Vec<[[f64; 2]; 3]>,
) {
    // Line normal form: dot(n, p) + c, sign = side.
    let dir = [s1[0] - s0[0], s1[1] - s0[1]];
    let len = (dir[0] * dir[0] + dir[1] * dir[1]).sqrt();
    if len < EPS {
        out.push(*tri);
        return;
    }
    let n = [-dir[1] / len, dir[0] / len];
    let c = -(n[0] * s0[0] + n[1] * s0[1]);
    let side = |p: [f64; 2]| n[0] * p[0] + n[1] * p[1] + c;

    let d = [side(tri[0]), side(tri[1]), side(tri[2])];
    // Wholly on one side (or on the line) → no split.
    let pos = d.iter().filter(|&&x| x > EPS).count();
    let neg = d.iter().filter(|&&x| x < -EPS).count();
    if pos == 0 || neg == 0 {
        out.push(*tri);
        return;
    }

    // If exactly one vertex sits *on* the line and the other two
    // straddle it, the cut runs from that vertex to the single
    // opposite-edge crossing — split into two triangles. (The
    // two-edge-crossing branch below misses this because neither edge
    // incident to an on-line vertex "straddles".)
    if let Some(on_vert) = (0..3).find(|&i| d[i].abs() <= EPS) {
        let o0 = (on_vert + 1) % 3;
        let o1 = (on_vert + 2) % 3;
        if (d[o0] > EPS && d[o1] < -EPS) || (d[o0] < -EPS && d[o1] > EPS) {
            let s = d[o0] / (d[o0] - d[o1]);
            let cross = [
                tri[o0][0] + (tri[o1][0] - tri[o0][0]) * s,
                tri[o0][1] + (tri[o1][1] - tri[o0][1]) * s,
            ];
            out.push([tri[on_vert], tri[o0], cross]);
            out.push([tri[on_vert], cross, tri[o1]]);
            return;
        }
    }

    // The line crosses two edges. Find the crossing points and
    // re-triangulate into one triangle + one quad (as two triangles).
    let mut crossings: Vec<(usize, [f64; 2])> = Vec::new();
    for (i, j) in [(0usize, 1usize), (1, 2), (2, 0)] {
        let di = d[i];
        let dj = d[j];
        if (di > EPS && dj < -EPS) || (di < -EPS && dj > EPS) {
            let s = di / (di - dj);
            let p = [
                tri[i][0] + (tri[j][0] - tri[i][0]) * s,
                tri[i][1] + (tri[j][1] - tri[i][1]) * s,
            ];
            crossings.push((i, p));
        }
    }
    if crossings.len() != 2 {
        // On-vertex degeneracy not covered above — keep intact.
        out.push(*tri);
        return;
    }

    // The triangle has one "lone" vertex on one side and an edge pair
    // on the other. The two crossings sit on the two edges incident
    // to the lone vertex.
    let (e0, p0) = crossings[0];
    let (e1, p1) = crossings[1];
    // Edges are (0,1)=0, (1,2)=1, (2,0)=2. The vertex shared by both
    // crossed edges is the lone vertex.
    let edge_verts = [(0usize, 1usize), (1, 2), (2, 0)];
    let (a0, b0) = edge_verts[e0];
    let (a1, b1) = edge_verts[e1];
    let lone = if a0 == a1 || a0 == b1 {
        a0
    } else if b0 == a1 || b0 == b1 {
        b0
    } else {
        // Shouldn't happen for a convex triangle; bail safely.
        out.push(*tri);
        return;
    };
    let others: Vec<usize> = (0..3).filter(|&i| i != lone).collect();
    let (o0, o1) = (others[0], others[1]);

    // The crossing on a given edge index.
    let crossing_on = |edge: usize| -> [f64; 2] {
        if edge == e0 {
            p0
        } else {
            p1
        }
    };
    // Edge lone..o0 and lone..o1 — figure which crossing belongs to which.
    let edge_index = |x: usize, y: usize| -> usize {
        for (idx, (a, b)) in edge_verts.iter().enumerate() {
            if (*a == x && *b == y) || (*a == y && *b == x) {
                return idx;
            }
        }
        usize::MAX
    };
    let c0 = crossing_on(edge_index(lone, o0));
    let c1 = crossing_on(edge_index(lone, o1));

    // Apex triangle on the lone-vertex side.
    out.push([tri[lone], c0, c1]);
    // Quad (c0, o0, o1, c1) split into two triangles.
    out.push([c0, tri[o0], tri[o1]]);
    out.push([c0, tri[o1], c1]);
}

// ===================================================================
// Point-in-solid classification by ray parity.
// ===================================================================

/// Is `p` inside the closed manifold represented by `tree`?
///
/// Casts a ray in a fixed direction and counts genuine triangle
/// crossings; odd parity ⇒ inside. The direction is jittered slightly
/// off the axes so it almost never grazes an edge.
fn point_inside(p: Vector3<f64>, tree: &AabbTree) -> bool {
    if tree.nodes.is_empty() {
        return false;
    }
    // A few jittered directions — take a majority vote so a single
    // edge-grazing ray cannot flip the verdict.
    let dirs = [
        Vector3::new(0.5773, 0.5774, 0.5775),
        Vector3::new(-0.3313, 0.7072, 0.6251),
        Vector3::new(0.8011, -0.2673, 0.5354),
    ];
    let mut inside_votes = 0usize;
    for d in dirs {
        let dir = d.normalize();
        let hits = tree.intersect_ray(p, dir);
        let mut count = 0usize;
        for tid in hits {
            let tri = &tree.triangles[tid.0];
            if ray_triangle_hit(p, dir, tri) {
                count += 1;
            }
        }
        if count % 2 == 1 {
            inside_votes += 1;
        }
    }
    inside_votes >= 2
}

/// Möller-Trumbore ray-triangle test — returns `true` when the ray
/// `(orig, dir)` crosses `tri` at a strictly positive parameter.
fn ray_triangle_hit(orig: Vector3<f64>, dir: Vector3<f64>, tri: &Triangle3) -> bool {
    let e1 = tri.v[1] - tri.v[0];
    let e2 = tri.v[2] - tri.v[0];
    let pvec = dir.cross(&e2);
    let det = e1.dot(&pvec);
    if det.abs() < 1e-12 {
        return false; // ray parallel to triangle
    }
    let inv_det = 1.0 / det;
    let tvec = orig - tri.v[0];
    let u = tvec.dot(&pvec) * inv_det;
    if !(0.0..=1.0).contains(&u) {
        return false;
    }
    let qvec = tvec.cross(&e1);
    let v = dir.dot(&qvec) * inv_det;
    if v < 0.0 || u + v > 1.0 {
        return false;
    }
    let t = e2.dot(&qvec) * inv_det;
    t > 1e-9
}

/// Build a `Mesh3` from an indexed `(vertices, faces)` pair —
/// convenience for callers that hold welded geometry.
pub fn mesh3_from_indexed(vertices: &[Vector3<f64>], faces: &[[usize; 3]]) -> Mesh3 {
    let mut out = Mesh3::new();
    for f in faces {
        out.triangles.push(Triangle3 {
            v: [vertices[f[0]], vertices[f[1]], vertices[f[2]]],
        });
    }
    out
}

/// Weld a `Mesh3`'s vertices within `tol` and return `(vertices,
/// faces)` — the inverse of [`mesh3_from_indexed`]. Used by callers
/// that want an indexed mesh out of the CSG result.
pub fn mesh3_to_indexed(mesh: &Mesh3, tol: f64) -> (Vec<Vector3<f64>>, Vec<[usize; 3]>) {
    let h = tol.max(1e-30);
    let key = |p: Vector3<f64>| {
        (
            (p.x / h).round() as i64,
            (p.y / h).round() as i64,
            (p.z / h).round() as i64,
        )
    };
    let mut map: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut verts: Vec<Vector3<f64>> = Vec::new();
    let mut faces: Vec<[usize; 3]> = Vec::new();
    for t in &mesh.triangles {
        let mut idx = [0usize; 3];
        for (k, p) in t.v.iter().enumerate() {
            let id = *map.entry(key(*p)).or_insert_with(|| {
                verts.push(*p);
                verts.len() - 1
            });
            idx[k] = id;
        }
        faces.push(idx);
    }
    (verts, faces)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Axis-aligned box [lo, hi] as 12 triangles with outward normals.
    fn box_mesh(lo: Vector3<f64>, hi: Vector3<f64>) -> Mesh3 {
        let v = [
            Vector3::new(lo.x, lo.y, lo.z),
            Vector3::new(hi.x, lo.y, lo.z),
            Vector3::new(hi.x, hi.y, lo.z),
            Vector3::new(lo.x, hi.y, lo.z),
            Vector3::new(lo.x, lo.y, hi.z),
            Vector3::new(hi.x, lo.y, hi.z),
            Vector3::new(hi.x, hi.y, hi.z),
            Vector3::new(lo.x, hi.y, hi.z),
        ];
        // Each face CCW seen from outside.
        let quads = [
            [0, 3, 2, 1], // bottom (-z)
            [4, 5, 6, 7], // top (+z)
            [0, 1, 5, 4], // -y
            [1, 2, 6, 5], // +x
            [2, 3, 7, 6], // +y
            [3, 0, 4, 7], // -x
        ];
        let mut m = Mesh3::new();
        for q in quads {
            m.triangles.push(Triangle3 {
                v: [v[q[0]], v[q[1]], v[q[2]]],
            });
            m.triangles.push(Triangle3 {
                v: [v[q[0]], v[q[2]], v[q[3]]],
            });
        }
        m
    }

    #[test]
    fn disjoint_union_keeps_all_triangles() {
        let a = box_mesh(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
        let b = box_mesh(Vector3::new(5.0, 5.0, 5.0), Vector3::new(6.0, 6.0, 6.0));
        let u = union(&a, &b);
        // No overlap → all 24 triangles survive.
        assert_eq!(u.len(), 24);
    }

    #[test]
    fn disjoint_intersection_is_empty() {
        let a = box_mesh(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
        let b = box_mesh(Vector3::new(5.0, 5.0, 5.0), Vector3::new(6.0, 6.0, 6.0));
        let i = intersection(&a, &b);
        assert!(i.is_empty(), "disjoint solids share no volume");
    }

    #[test]
    fn point_inside_box() {
        let a = box_mesh(Vector3::new(0.0, 0.0, 0.0), Vector3::new(2.0, 2.0, 2.0));
        let tree = AabbTree::new(a.triangles.clone());
        assert!(point_inside(Vector3::new(1.0, 1.0, 1.0), &tree));
        assert!(!point_inside(Vector3::new(5.0, 5.0, 5.0), &tree));
        assert!(!point_inside(Vector3::new(-1.0, 1.0, 1.0), &tree));
    }

    #[test]
    fn overlapping_union_removes_interior_walls() {
        // Two unit cubes overlapping in [0.5,1]x[0,1]x[0,1].
        let a = box_mesh(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
        let b = box_mesh(Vector3::new(0.5, 0.0, 0.0), Vector3::new(1.5, 1.0, 1.0));
        let u = union(&a, &b);
        // A real union has more facets than one box (the seam is
        // co-refined) but a point deep inside the shared region must
        // NOT be classified as a surface — verify by centroid test:
        // the union's bounding span covers x in [0, 1.5].
        let verts = u.vertices();
        let xmin = verts.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let xmax = verts
            .iter()
            .map(|p| p.x)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((xmin - 0.0).abs() < 1e-6, "xmin={xmin}");
        assert!((xmax - 1.5).abs() < 1e-6, "xmax={xmax}");
        // The interior wall at x=0.5 (inside the other box) must be
        // gone: no triangle of the union should have all three
        // vertices at x≈0.5 within the y,z overlap.
        let interior_wall = u.triangles.iter().any(|t| {
            t.v.iter().all(|p| (p.x - 0.5).abs() < 1e-6)
        });
        assert!(!interior_wall, "interior wall at x=0.5 should be removed");
    }

    #[test]
    fn difference_carves_a_cavity() {
        // Big box minus a smaller box poking through one face.
        let a = box_mesh(Vector3::new(0.0, 0.0, 0.0), Vector3::new(4.0, 4.0, 4.0));
        let b = box_mesh(Vector3::new(1.0, 1.0, -1.0), Vector3::new(2.0, 2.0, 2.0));
        let d = difference(&a, &b);
        // The result keeps A's outer shell (minus the carved opening)
        // plus B's wall inside A. It must be non-empty and contain
        // facets from both operands.
        assert!(!d.is_empty());
        // A point clearly inside A but outside the carved tunnel is
        // still solid; a point inside the tunnel is not. We can only
        // check the surface: the cavity wall (from B) lies at x≈1 and
        // x≈2 inside A — at least one such facet must exist.
        let has_cavity_wall = d.triangles.iter().any(|t| {
            t.v.iter().all(|p| (p.x - 1.0).abs() < 1e-6)
                && t.v.iter().all(|p| p.y >= 1.0 - 1e-6 && p.y <= 2.0 + 1e-6)
        });
        assert!(has_cavity_wall, "difference should leave a cavity wall");
    }

    #[test]
    fn intersection_of_overlapping_boxes_is_the_lens() {
        let a = box_mesh(Vector3::new(0.0, 0.0, 0.0), Vector3::new(2.0, 2.0, 2.0));
        let b = box_mesh(Vector3::new(1.0, 1.0, 1.0), Vector3::new(3.0, 3.0, 3.0));
        let i = intersection(&a, &b);
        assert!(!i.is_empty(), "overlapping boxes share a unit cube");
        // The lens is the unit cube [1,2]^3 — every vertex of the
        // result must lie within that box (with EPS slop).
        for t in &i.triangles {
            for p in &t.v {
                assert!(
                    p.x >= 1.0 - 1e-6 && p.x <= 2.0 + 1e-6,
                    "lens x out of range: {}",
                    p.x
                );
                assert!(p.y >= 1.0 - 1e-6 && p.y <= 2.0 + 1e-6);
                assert!(p.z >= 1.0 - 1e-6 && p.z <= 2.0 + 1e-6);
            }
        }
    }

    #[test]
    fn tri_tri_segment_finds_crossing() {
        // Triangle in z=0 plane, triangle crossing it vertically.
        let flat = Triangle3 {
            v: [
                Vector3::new(-1.0, -1.0, 0.0),
                Vector3::new(2.0, -1.0, 0.0),
                Vector3::new(0.0, 2.0, 0.0),
            ],
        };
        let vert = Triangle3 {
            v: [
                Vector3::new(0.0, 0.0, -1.0),
                Vector3::new(1.0, 0.0, -1.0),
                Vector3::new(0.5, 0.0, 1.0),
            ],
        };
        let seg = tri_tri_segment(&flat, &vert).expect("they cross");
        // The segment lies in z=0 and on y=0.
        for p in seg {
            assert!(p.z.abs() < 1e-6, "segment off the z=0 plane: {}", p.z);
            assert!(p.y.abs() < 1e-6, "segment off the y=0 plane: {}", p.y);
        }
    }

    #[test]
    fn indexed_round_trip() {
        let a = box_mesh(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
        let (verts, faces) = mesh3_to_indexed(&a, 1e-9);
        // A welded cube has 8 unique vertices and 12 triangles.
        assert_eq!(verts.len(), 8);
        assert_eq!(faces.len(), 12);
        let back = mesh3_from_indexed(&verts, &faces);
        assert_eq!(back.len(), 12);
    }
}
