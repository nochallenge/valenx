//! Plane-cut operations on a canonical [`Mesh`].
//!
//! Two complementary surface-level primitives:
//!
//! - [`intersect_plane`] — returns the **cross-section edges** where
//!   the cut plane hits each surface triangle. Used to draw the cut
//!   overlay in the viewport so the user can preview where the slice
//!   will land before they commit.
//! - [`slice()`] — actually trims the mesh, keeping only the elements
//!   whose centroids lie on the positive side of the plane
//!   (`(centroid - point) · normal >= 0`). Centroid-based culling is
//!   the pre-alpha approximation; a proper polygon-clipping pass that
//!   produces fresh boundary triangles along the cut is its own
//!   follow-up.
//!
//! ## Scope
//!
//! These operations target triangular surface meshes today — the
//! shape an imported STL turns into when it's promoted to canonical
//! form. Volume elements (Tet4, Hex8, …) are handled by
//! centroid-keep for [`slice()`], but [`intersect_plane`] only emits
//! edges from `ElementType::Tri3` blocks since that's all the
//! current overlay needs. Adding volume slicing is a planned
//! follow-up once the BRep kernel lands.

use nalgebra::Vector3;

use crate::element::{ElementBlock, ElementType};
use crate::mesh::Mesh;

/// One line segment of the surface's intersection with a cutting
/// plane — the two endpoints lie on the plane and on the original
/// triangle's edges.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LineSegment {
    pub a: Vector3<f64>,
    pub b: Vector3<f64>,
}

/// Sign-comparison tolerance shared by every plane-intersection
/// path. Vertices with `|signed distance| < EPS` count as lying *on*
/// the plane so floating-point noise doesn't spawn spurious slivers.
const EPS: f64 = 1e-12;

/// Intersect a single triangle (three world-space vertices) with the
/// plane defined by `point` / `normal`, appending at most one
/// [`LineSegment`] to `out`.
///
/// This is the per-triangle core shared by [`intersect_plane`] (which
/// walks a canonical [`Mesh`]'s Tri3 blocks) and
/// [`intersect_plane_triangles`] (which walks a raw triangle soup).
/// Keeping the geometry in one place means both entry points stay
/// bit-for-bit identical and only get tested once.
///
/// A triangle straddling the plane (vertices on both sides, or one
/// edge crossing while the opposite vertex sits on the plane) yields
/// exactly one segment. Triangles entirely on one side, or coincident
/// with the plane, contribute nothing.
fn intersect_triangle(
    v: &[Vector3<f64>; 3],
    point: Vector3<f64>,
    normal: Vector3<f64>,
    out: &mut Vec<LineSegment>,
) {
    let d: [f64; 3] = [
        (v[0] - point).dot(&normal),
        (v[1] - point).dot(&normal),
        (v[2] - point).dot(&normal),
    ];
    // All-same-side cases: skip.
    let positive = d.iter().filter(|x| **x > EPS).count();
    let negative = d.iter().filter(|x| **x < -EPS).count();
    if positive == 3 || negative == 3 {
        return;
    }
    // Coplanar triangle (all three on the plane): one full triangle
    // outline would technically be the intersection, but emitting
    // three degenerate edges from every coplanar triangle would spam
    // the overlay. Skip.
    if positive == 0 && negative == 0 {
        return;
    }
    // Walk the three triangle edges, recording the two crossings
    // between positive- and negative-side vertices.
    let mut crossings: [Option<Vector3<f64>>; 2] = [None, None];
    let mut idx = 0;
    for k in 0..3 {
        let a = v[k];
        let b = v[(k + 1) % 3];
        let da = d[k];
        let db = d[(k + 1) % 3];
        // Strict sign-change crossing — both endpoints off the plane,
        // on opposite sides.
        if (da > EPS && db < -EPS) || (da < -EPS && db > EPS) {
            let t = da / (da - db);
            let p = a + (b - a) * t;
            if idx < 2 {
                crossings[idx] = Some(p);
                idx += 1;
            }
        }
    }
    if let (Some(a), Some(b)) = (crossings[0], crossings[1]) {
        out.push(LineSegment { a, b });
    } else if idx == 1 {
        // One edge crosses, the opposite vertex sits on the plane —
        // emit segment between the crossing and that on-plane vertex.
        let on_plane = (0..3).find(|k| d[*k].abs() <= EPS).map(|k| v[k]);
        if let (Some(c), Some(o)) = (crossings[0], on_plane) {
            out.push(LineSegment { a: c, b: o });
        }
    }
}

/// Return every cross-section edge produced by intersecting the
/// mesh's surface triangles with the plane defined by `point` (any
/// point on the plane) and `normal` (its outward direction; need not
/// be unit length — the algorithm only compares signs of
/// `(vertex - point) · normal`).
///
/// Each triangle with vertices straddling the plane (one or two
/// vertices on the negative side and the rest on the non-negative
/// side) yields exactly one segment between the two edges that cross.
/// Triangles entirely on one side, or coincident with the plane, are
/// skipped. Floating-point ties (`|signed| < EPS`) are treated as
/// "on the plane" so we don't generate spurious tiny segments from
/// numerical noise.
pub fn intersect_plane(mesh: &Mesh, point: Vector3<f64>, normal: Vector3<f64>) -> Vec<LineSegment> {
    let mut out = Vec::new();
    // R34 S2 (defense-in-depth): `tri[k]` are connectivity values, so
    // an out-of-range index would panic `mesh.nodes[..]`. Use `.get()`
    // and skip any triangle citing a vertex past `nodes.len()` — the
    // overlay simply omits a degenerate triangle rather than crashing.
    // Backs the per-loader parse guards (OBJ/gmsh/netgen/PLY).
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let (Some(&v0), Some(&v1), Some(&v2)) = (
                mesh.nodes.get(tri[0] as usize),
                mesh.nodes.get(tri[1] as usize),
                mesh.nodes.get(tri[2] as usize),
            ) else {
                continue;
            };
            let v: [Vector3<f64>; 3] = [v0, v1, v2];
            intersect_triangle(&v, point, normal, &mut out);
        }
    }
    out
}

/// Plane-intersection over a raw triangle soup — same algorithm as
/// [`intersect_plane`] but takes `&[[[f64; 3]; 3]]` directly so
/// callers holding loose STL triangles don't need to round-trip
/// through a canonical [`Mesh`] (allocating node arrays, element
/// blocks, and a connectivity remap) just to preview a cut.
///
/// Each element of `triangles` is one triangle as three `[x, y, z]`
/// corner positions. Output and skip rules are identical to
/// [`intersect_plane`]: a straddling triangle yields one segment,
/// triangles on one side or coplanar yield nothing.
pub fn intersect_plane_triangles(
    triangles: &[[[f64; 3]; 3]],
    point: Vector3<f64>,
    normal: Vector3<f64>,
) -> Vec<LineSegment> {
    let mut out = Vec::new();
    for t in triangles {
        let v: [Vector3<f64>; 3] = [
            Vector3::new(t[0][0], t[0][1], t[0][2]),
            Vector3::new(t[1][0], t[1][1], t[1][2]),
            Vector3::new(t[2][0], t[2][1], t[2][2]),
        ];
        intersect_triangle(&v, point, normal, &mut out);
    }
    out
}

/// Return a fresh mesh containing only the elements whose centroid
/// satisfies `(centroid - point) · normal >= 0`. The original mesh
/// is unmodified.
///
/// Node arrays are deduplicated to just those referenced by surviving
/// elements; connectivity indices are remapped accordingly. `id` of
/// the returned mesh is `"<original>_cut"`. Regions and boundary
/// groups are dropped — they reference positional element indices
/// that don't survive the cull, and rebuilding them correctly would
/// need topological knowledge we don't have yet.
///
/// **Pre-alpha approximation**: centroid-keep means elements
/// straddling the plane are kept or discarded wholesale, not
/// clipped. A proper polygon-clip pass that emits fresh boundary
/// triangles along the cut is a planned follow-up — the simpler
/// path here covers the "drag the slider until the mesh looks
/// roughly right" interactive workflow today.
pub fn slice(mesh: &Mesh, point: Vector3<f64>, normal: Vector3<f64>) -> Mesh {
    let mut out = Mesh::new(format!("{}_cut", mesh.id));
    if mesh.nodes.is_empty() {
        return out;
    }

    // Per-element keep decision: centroid relative to the plane.
    let mut kept_nodes: Vec<u32> = vec![u32::MAX; mesh.nodes.len()];
    let mut next_idx: u32 = 0;

    for block in &mesh.element_blocks {
        let stride = block.element_type.nodes_per_element();
        if stride == 0 {
            continue;
        }
        let mut new_block = ElementBlock::new(block.element_type);
        'elem: for chunk in block.connectivity.chunks_exact(stride) {
            // R34 S2 (defense-in-depth): `chunk` holds connectivity
            // values that index `mesh.nodes` (centroid) and
            // `kept_nodes` (the remap). An out-of-range index would
            // panic either. Validate the WHOLE element up front and
            // drop it if any vertex is past `nodes.len()` — we must
            // not half-mutate `kept_nodes`/`out.nodes` for a corrupt
            // element. Backs the per-loader parse guards.
            for &idx in chunk {
                if (idx as usize) >= mesh.nodes.len() {
                    continue 'elem;
                }
            }
            // Compute the element centroid in the original node array.
            let mut centroid = Vector3::zeros();
            for &idx in chunk {
                centroid += mesh.nodes[idx as usize];
            }
            centroid /= stride as f64;
            if (centroid - point).dot(&normal) < 0.0 {
                continue; // discard
            }
            // Keep: emit the chunk into the new block, remapping
            // each node index through the kept_nodes table so we
            // only carry the nodes the surviving elements reach.
            for &idx in chunk {
                let slot = &mut kept_nodes[idx as usize];
                if *slot == u32::MAX {
                    *slot = next_idx;
                    next_idx += 1;
                    out.nodes.push(mesh.nodes[idx as usize]);
                }
                new_block.connectivity.push(*slot);
            }
        }
        if !new_block.connectivity.is_empty() {
            out.element_blocks.push(new_block);
        }
    }
    out.recompute_stats();
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::{ElementBlock, ElementType};

    fn pt(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    /// Two-triangle unit square on the z=0 plane.
    fn unit_square_tris() -> Mesh {
        let mut m = Mesh::new("square");
        m.nodes = vec![
            pt(0.0, 0.0, 0.0),
            pt(1.0, 0.0, 0.0),
            pt(1.0, 1.0, 0.0),
            pt(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2, 0, 2, 3];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        m
    }

    /// Cube surface mesh: 12 triangles, 8 unique nodes. Each face
    /// is two triangles; we wind them so face normals point outward.
    fn unit_cube_surface() -> Mesh {
        let mut m = Mesh::new("cube");
        m.nodes = vec![
            pt(0.0, 0.0, 0.0), // 0
            pt(1.0, 0.0, 0.0), // 1
            pt(1.0, 1.0, 0.0), // 2
            pt(0.0, 1.0, 0.0), // 3
            pt(0.0, 0.0, 1.0), // 4
            pt(1.0, 0.0, 1.0), // 5
            pt(1.0, 1.0, 1.0), // 6
            pt(0.0, 1.0, 1.0), // 7
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        // Bottom (z=0)
        blk.connectivity.extend_from_slice(&[0, 2, 1, 0, 3, 2]);
        // Top (z=1)
        blk.connectivity.extend_from_slice(&[4, 5, 6, 4, 6, 7]);
        // Front (y=0)
        blk.connectivity.extend_from_slice(&[0, 1, 5, 0, 5, 4]);
        // Back (y=1)
        blk.connectivity.extend_from_slice(&[2, 3, 7, 2, 7, 6]);
        // Left (x=0)
        blk.connectivity.extend_from_slice(&[0, 4, 7, 0, 7, 3]);
        // Right (x=1)
        blk.connectivity.extend_from_slice(&[1, 2, 6, 1, 6, 5]);
        m.element_blocks = vec![blk];
        m.recompute_stats();
        m
    }

    #[test]
    fn intersect_plane_far_away_returns_nothing() {
        // A plane sitting at z = 10 doesn't touch the z=0 square.
        let m = unit_square_tris();
        let segs = intersect_plane(&m, pt(0.0, 0.0, 10.0), pt(0.0, 0.0, 1.0));
        assert!(segs.is_empty());
    }

    #[test]
    fn intersect_plane_cuts_cube_through_midline() {
        // Plane x = 0.5 cuts the unit cube into two halves. Each of
        // the four side faces (front, back, left, right is unaffected
        // — wait, left and right ARE coplanar with the cut at x=0).
        // Cleaner: count we get at least 4 segments (one per face
        // the plane crosses without being coplanar).
        let m = unit_cube_surface();
        let segs = intersect_plane(&m, pt(0.5, 0.0, 0.0), pt(1.0, 0.0, 0.0));
        // 4 faces (top, bottom, front, back) each contribute 2 tris,
        // and each tri that straddles the plane produces 1 segment.
        // Per face: 2 tris straddle → 2 segments. 4 faces × 2 = 8.
        assert_eq!(segs.len(), 8);
    }

    #[test]
    fn slice_keeps_half_of_unit_square() {
        // Plane through (0.5, 0, 0) with +x normal keeps the right half.
        let m = unit_square_tris();
        let cut = slice(&m, pt(0.5, 0.0, 0.0), pt(1.0, 0.0, 0.0));
        // Tri (0, 1, 2) has centroid x ≈ 0.667 — kept.
        // Tri (0, 2, 3) has centroid x ≈ 0.333 — discarded.
        assert_eq!(cut.element_blocks.len(), 1);
        assert_eq!(cut.element_blocks[0].connectivity.len(), 3);
        // 3 unique nodes (0, 1, 2 of the original) survive.
        assert_eq!(cut.nodes.len(), 3);
    }

    #[test]
    fn slice_keeps_everything_when_plane_is_outside() {
        // Plane way below the square (z=-10) with +z normal: every
        // centroid is on the positive side → keep all.
        let m = unit_square_tris();
        let cut = slice(&m, pt(0.0, 0.0, -10.0), pt(0.0, 0.0, 1.0));
        assert_eq!(cut.element_blocks[0].connectivity.len(), 6);
        assert_eq!(cut.nodes.len(), 4);
    }

    #[test]
    fn slice_discards_everything_when_plane_is_above() {
        // Plane above the square with +z normal: every centroid is
        // on the negative side → discard all → empty mesh.
        let m = unit_square_tris();
        let cut = slice(&m, pt(0.0, 0.0, 10.0), pt(0.0, 0.0, 1.0));
        assert!(cut.element_blocks.is_empty() || cut.element_blocks[0].connectivity.is_empty());
        assert!(cut.nodes.is_empty());
    }

    #[test]
    fn slice_cube_in_half_keeps_six_triangles_or_more() {
        // x=0.5 plane bisects the 12-triangle cube. Centroid-keep:
        // bottom face (z=0) has 2 tris, each centroid is at x≈0.667
        // and x≈0.333 — one kept, one tossed. Same for top, front,
        // back. Left face's both centroids are at x≈0 (discard).
        // Right face's both centroids are at x=1 (keep both).
        // Total kept: 1 + 1 + 1 + 1 + 0 + 2 = 6.
        let m = unit_cube_surface();
        let cut = slice(&m, pt(0.5, 0.0, 0.0), pt(1.0, 0.0, 0.0));
        let kept_tris: usize = cut
            .element_blocks
            .iter()
            .map(|b| b.connectivity.len() / 3)
            .sum();
        assert_eq!(kept_tris, 6, "expected 6 kept triangles, got {kept_tris}");
        // Cached stats reflect the deduplicated node array.
        assert_eq!(cut.stats.node_count as usize, cut.nodes.len());
        assert_eq!(cut.stats.element_count as usize, kept_tris);
    }

    #[test]
    fn slice_id_carries_origin() {
        let m = unit_square_tris();
        let cut = slice(&m, pt(0.5, 0.0, 0.0), pt(1.0, 0.0, 0.0));
        assert_eq!(cut.id, "square_cut");
    }

    #[test]
    fn slice_remaps_node_indices_compactly() {
        // After culling, the surviving block's connectivity must
        // reference 0..cut.nodes.len() — no gaps, no out-of-range.
        let m = unit_cube_surface();
        let cut = slice(&m, pt(0.5, 0.0, 0.0), pt(1.0, 0.0, 0.0));
        let max_node = cut.nodes.len() as u32;
        for block in &cut.element_blocks {
            for &idx in &block.connectivity {
                assert!(idx < max_node, "stale index {idx} >= {max_node}");
            }
        }
    }

    #[test]
    fn slice_empty_mesh_returns_empty() {
        let m = Mesh::new("empty");
        let cut = slice(&m, pt(0.0, 0.0, 0.0), pt(0.0, 0.0, 1.0));
        assert!(cut.nodes.is_empty());
        assert!(cut.element_blocks.is_empty());
    }

    #[test]
    fn intersect_plane_triangles_matches_mesh_path() {
        // The triangle-soup entry point must produce exactly the same
        // segments as the Mesh path for the same geometry — both
        // funnel through `intersect_triangle`, this guards the wiring.
        let m = unit_cube_surface();
        let mesh_segs = intersect_plane(&m, pt(0.5, 0.0, 0.0), pt(1.0, 0.0, 0.0));
        // Rebuild the same 12 triangles as a raw soup.
        let mut soup: Vec<[[f64; 3]; 3]> = Vec::new();
        for block in &m.element_blocks {
            for tri in block.connectivity.chunks_exact(3) {
                soup.push([
                    [
                        m.nodes[tri[0] as usize].x,
                        m.nodes[tri[0] as usize].y,
                        m.nodes[tri[0] as usize].z,
                    ],
                    [
                        m.nodes[tri[1] as usize].x,
                        m.nodes[tri[1] as usize].y,
                        m.nodes[tri[1] as usize].z,
                    ],
                    [
                        m.nodes[tri[2] as usize].x,
                        m.nodes[tri[2] as usize].y,
                        m.nodes[tri[2] as usize].z,
                    ],
                ]);
            }
        }
        let soup_segs = intersect_plane_triangles(&soup, pt(0.5, 0.0, 0.0), pt(1.0, 0.0, 0.0));
        assert_eq!(soup_segs, mesh_segs);
        assert_eq!(soup_segs.len(), 8);
    }

    #[test]
    fn intersect_plane_triangles_single_tri_through_midline() {
        // One triangle in the z=0 plane spanning x∈[0,2]; a plane at
        // x=1 cuts it into exactly one segment between the two edges
        // that straddle.
        let soup = [[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 2.0, 0.0]]];
        let segs = intersect_plane_triangles(&soup, pt(1.0, 0.0, 0.0), pt(1.0, 0.0, 0.0));
        assert_eq!(segs.len(), 1);
        // Both endpoints land on the cutting plane x = 1.
        assert!((segs[0].a.x - 1.0).abs() < 1e-9);
        assert!((segs[0].b.x - 1.0).abs() < 1e-9);
    }

    #[test]
    fn intersect_plane_triangles_skips_when_plane_misses() {
        // A plane far above a flat triangle in z=0 touches nothing.
        let soup = [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]];
        let segs = intersect_plane_triangles(&soup, pt(0.0, 0.0, 5.0), pt(0.0, 0.0, 1.0));
        assert!(segs.is_empty());
        // Empty soup is also a no-op rather than a panic.
        let empty: [[[f64; 3]; 3]; 0] = [];
        assert!(intersect_plane_triangles(&empty, pt(0.0, 0.0, 0.0), pt(0.0, 0.0, 1.0)).is_empty());
    }

    /// R34 S2 (RED→GREEN): defense-in-depth sink seal. A mesh whose
    /// Tri3 connectivity cites a vertex past `nodes.len()` must NOT
    /// panic `intersect_plane` or `slice`. Pre-fix both did
    /// `mesh.nodes[tri[k] as usize]` / `mesh.nodes[idx as usize]` and
    /// panicked "index out of bounds". Post-fix the out-of-range
    /// triangle/element is dropped and a result is returned.
    #[test]
    fn out_of_range_connectivity_does_not_panic() {
        let mut m = Mesh::new("hostile");
        m.nodes = vec![pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0), pt(0.0, 1.0, 0.0)];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        // A triangle citing vertex 7 (out of range).
        blk.connectivity = vec![0, 1, 7];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        // Neither entry point may panic.
        let segs = intersect_plane(&m, pt(0.5, 0.0, 0.0), pt(1.0, 0.0, 0.0));
        assert!(segs.is_empty(), "the out-of-range triangle yields no segment");
        let cut = slice(&m, pt(0.0, 0.0, -10.0), pt(0.0, 0.0, 1.0));
        // The bad element is dropped, so no connectivity survives.
        let total: usize = cut.element_blocks.iter().map(|b| b.connectivity.len()).sum();
        assert_eq!(total, 0, "the out-of-range element must be dropped");
    }

    /// R34 S2: a valid triangle alongside an out-of-range one. The
    /// valid one drives the cut while the bad one is dropped, no panic.
    #[test]
    fn out_of_range_triangle_skipped_valid_kept() {
        let mut m = Mesh::new("mixed");
        m.nodes = vec![
            pt(0.0, 0.0, 0.0),
            pt(1.0, 0.0, 0.0),
            pt(1.0, 1.0, 0.0),
            pt(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        // First tri valid; second cites vertex 99.
        blk.connectivity = vec![0, 1, 2, 0, 2, 99];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        // Plane below the mesh with +z normal: the one valid triangle
        // is kept; the bad triangle is dropped.
        let cut = slice(&m, pt(0.0, 0.0, -10.0), pt(0.0, 0.0, 1.0));
        let kept_tris: usize = cut
            .element_blocks
            .iter()
            .map(|b| b.connectivity.len() / 3)
            .sum();
        assert_eq!(kept_tris, 1, "exactly the one valid triangle should survive");
    }

    #[test]
    fn intersect_plane_skips_non_tri3_blocks() {
        // A Hex8 block (volume) shouldn't yield surface intersections
        // — the intersect_plane helper only works on Tri3 surface
        // meshes today.
        let mut m = Mesh::new("hex");
        m.nodes = vec![
            pt(0.0, 0.0, 0.0),
            pt(1.0, 0.0, 0.0),
            pt(1.0, 1.0, 0.0),
            pt(0.0, 1.0, 0.0),
            pt(0.0, 0.0, 1.0),
            pt(1.0, 0.0, 1.0),
            pt(1.0, 1.0, 1.0),
            pt(0.0, 1.0, 1.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Hex8);
        blk.connectivity = vec![0, 1, 2, 3, 4, 5, 6, 7];
        m.element_blocks = vec![blk];
        let segs = intersect_plane(&m, pt(0.5, 0.0, 0.0), pt(1.0, 0.0, 0.0));
        assert!(segs.is_empty());
    }
}
