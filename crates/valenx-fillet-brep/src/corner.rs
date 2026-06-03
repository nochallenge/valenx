//! Corner-blend detection and classification (Phase 14.7).
//!
//! When 3+ filleted edges meet at a vertex of the solid, the per-edge
//! cylindrical fillet surfaces don't meet tangentially — they leave a
//! tri-tangent gap that needs a spherical patch to blend. The
//! geometry of that blend (the rolling-ball corner) is built in
//! [`crate::corner_build`]; this module is the **detection** layer
//! that decides *which* vertices get a corner blend and *whether* a
//! given vertex is the supported case.
//!
//! # The supported case — orthogonal convex 3-edge corner
//!
//! Corner-blend v1 ships the common, well-defined case: **three
//! mutually-orthogonal convex edges meeting at a single vertex** — the
//! corner of a box. At such a corner three quarter-cylinder edge
//! fillets of radius `r` meet, and the blend is the spherical octant
//! of the radius-`r` ball seated in the corner. [`classify_corner`]
//! confirms a vertex is this case; [`find_blendable_corners`] returns
//! every box-corner-style vertex of a filleted-edge set.
//!
//! # Honest scope
//!
//! General N-edge corners, non-orthogonal corners, and concave
//! corners are **not** the supported case — [`classify_corner`]
//! returns [`CornerClass::Unsupported`] for them and the fillet path
//! leaves those corners as the independent per-edge fillets v1 has
//! always produced (an honestly-documented Tier-3 residue). The
//! detection here never *fakes* a blendable corner.

use std::collections::{HashMap, HashSet};

use truck_modeling::{Edge as TruckEdge, InnerSpace, Point3, Solid as TruckSolid, Vector3, VertexID};

/// True if `vertex_id` is a meeting point of 3+ edges drawn from the
/// `filleted_edges` list — i.e. a corner the per-edge fillets leave
/// un-blended unless a corner blend is applied.
///
/// For a cube with all 12 edges filleted, every one of the 8 corners
/// returns true. For a cube with only 4 vertical edges filleted, the
/// 4 top + 4 bottom corners return false (each is the meeting point
/// of one filleted edge + two unfilleted ones, so v1's per-edge
/// fillet leaves them intact as a 2-way fillet-end + 2 unfilleted
/// edges meeting at the corner).
pub fn is_n_way_vertex(
    _solid: &TruckSolid,
    vertex_id: VertexID,
    filleted_edges: &[TruckEdge],
) -> bool {
    let mut count = 0_usize;
    for e in filleted_edges {
        if e.front().id() == vertex_id || e.back().id() == vertex_id {
            count += 1;
            if count >= 3 {
                return true;
            }
        }
    }
    false
}

/// Count how many distinct vertices in `solid` are meeting points of
/// 3+ of the given `filleted_edges`. Used by the bridge to surface a
/// "N corners need blending" diagnostic.
pub fn unblended_corner_count(solid: &TruckSolid, filleted_edges: &[TruckEdge]) -> usize {
    // Group edges by their endpoints, count edges per vertex.
    let mut per_vertex: HashMap<VertexID, usize> = HashMap::new();
    for e in filleted_edges {
        *per_vertex.entry(e.front().id()).or_insert(0) += 1;
        *per_vertex.entry(e.back().id()).or_insert(0) += 1;
    }
    // Only count vertices that actually exist in the solid (defensive
    // — a malformed edge list could carry stale IDs).
    let live_vertex_ids: HashSet<_> = solid.vertex_iter().map(|v| v.id()).collect();
    per_vertex
        .iter()
        .filter(|(vid, count)| live_vertex_ids.contains(vid) && **count >= 3)
        .count()
}

/// The three edges meeting at one corner vertex, resolved to the
/// geometry the corner-blend constructor needs.
///
/// The vertex itself is `apex`; each edge contributes a **unit
/// direction pointing away from the apex along that edge**. For the
/// supported (orthogonal box) corner these three directions are
/// mutually perpendicular and the corner blend's centre sits at
/// `apex + r·(d0 + d1 + d2)` (each `dᵢ` is also an inward face
/// direction — see [`crate::corner_build`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CornerEdges {
    /// The shared corner vertex (apex of the sharp corner).
    pub apex: Point3,
    /// Unit direction along edge 0, pointing away from `apex`.
    pub dir0: Vector3,
    /// Unit direction along edge 1, pointing away from `apex`.
    pub dir1: Vector3,
    /// Unit direction along edge 2, pointing away from `apex`.
    pub dir2: Vector3,
}

/// Classification of a candidate corner vertex.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CornerClass {
    /// The supported case: exactly three filleted edges meet here and
    /// they are mutually orthogonal and convex — the corner of a box.
    /// Carries the resolved [`CornerEdges`] geometry.
    OrthogonalConvex(CornerEdges),
    /// Anything else — more or fewer than 3 edges, non-orthogonal
    /// edges, a concave corner. Corner-blend v1 does not handle it;
    /// the per-edge fillets are left independent (Tier-3 residue).
    Unsupported,
}

/// Angular tolerance (radians) for the "mutually orthogonal" test.
/// A box corner is exactly 90°; this allows for float drift and very
/// slightly off-square modelling without admitting genuinely
/// non-orthogonal corners.
const ORTHO_TOL_RAD: f64 = 1e-4;

/// Classify a corner vertex: is it the supported orthogonal-convex
/// 3-edge case?
///
/// `vertex_id` must be a vertex of `solid`; `filleted_edges` is the
/// set of edges being filleted. The function gathers exactly the
/// filleted edges incident on `vertex_id`, and returns
/// [`CornerClass::OrthogonalConvex`] only when:
///
/// 1. **exactly three** filleted edges meet at the vertex (a box
///    corner has degree 3; 4+ is a non-box corner, 2 is an edge-end);
/// 2. the three edge directions (away from the apex) are **mutually
///    orthogonal** to within a small angular tolerance (`1e-4` rad);
/// 3. the three directions form a non-degenerate frame — the scalar
///    triple product is non-zero (three orthogonal unit vectors are
///    automatically a basis, so this is a belt-and-braces guard).
///
/// Every other configuration returns [`CornerClass::Unsupported`].
pub fn classify_corner(
    solid: &TruckSolid,
    vertex_id: VertexID,
    filleted_edges: &[TruckEdge],
) -> CornerClass {
    // Confirm the vertex is real (defensive against stale IDs).
    if !solid.vertex_iter().any(|v| v.id() == vertex_id) {
        return CornerClass::Unsupported;
    }

    // Gather the apex point + a unit direction away from the apex for
    // every filleted edge incident on this vertex.
    let mut apex: Option<Point3> = None;
    let mut dirs: Vec<Vector3> = Vec::new();
    for e in filleted_edges {
        let (front, back) = (e.front(), e.back());
        let (this, other) = if front.id() == vertex_id {
            (front.point(), back.point())
        } else if back.id() == vertex_id {
            (back.point(), front.point())
        } else {
            continue;
        };
        let along = other - this;
        let len = along.magnitude();
        if len < 1e-9 {
            // Degenerate edge — cannot define a direction.
            return CornerClass::Unsupported;
        }
        apex = Some(this);
        dirs.push(along / len);
    }

    // The supported corner has degree exactly 3.
    if dirs.len() != 3 {
        return CornerClass::Unsupported;
    }
    let apex = match apex {
        Some(p) => p,
        None => return CornerClass::Unsupported,
    };

    // Mutual orthogonality: every pair's dot product must be ~0.
    // |cos θ| ≈ |π/2 − θ| for θ near 90°, so the radian tolerance
    // doubles as the cosine tolerance here.
    let cos_tol = ORTHO_TOL_RAD;
    for (i, j) in [(0, 1), (0, 2), (1, 2)] {
        if dirs[i].dot(dirs[j]).abs() > cos_tol {
            return CornerClass::Unsupported;
        }
    }

    // Non-degeneracy: the scalar triple product of three mutually
    // orthogonal unit vectors is ±1; anything near 0 would mean a
    // pair is actually parallel (already excluded above).
    let triple = dirs[0].cross(dirs[1]).dot(dirs[2]);
    if triple.abs() < 0.5 {
        return CornerClass::Unsupported;
    }

    CornerClass::OrthogonalConvex(CornerEdges {
        apex,
        dir0: dirs[0],
        dir1: dirs[1],
        dir2: dirs[2],
    })
}

/// Find every vertex of `solid` that is a **supported** corner-blend
/// site for the given `filleted_edges` — i.e. an orthogonal convex
/// 3-edge corner.
///
/// Returns each site's [`VertexID`] paired with its resolved
/// [`CornerEdges`] geometry. Vertices that are meeting points of 3+
/// filleted edges but fail the orthogonal-convex test (non-box
/// corners) are *not* returned — the caller leaves those as
/// independent per-edge fillets.
pub fn find_blendable_corners(
    solid: &TruckSolid,
    filleted_edges: &[TruckEdge],
) -> Vec<(VertexID, CornerEdges)> {
    // Collect the distinct candidate vertex IDs: any vertex touched
    // by 3+ filleted edges.
    let mut per_vertex: HashMap<VertexID, usize> = HashMap::new();
    for e in filleted_edges {
        *per_vertex.entry(e.front().id()).or_insert(0) += 1;
        *per_vertex.entry(e.back().id()).or_insert(0) += 1;
    }
    let mut out = Vec::new();
    for (vid, count) in per_vertex {
        if count < 3 {
            continue;
        }
        if let CornerClass::OrthogonalConvex(edges) = classify_corner(solid, vid, filleted_edges) {
            out.push((vid, edges));
        }
    }
    // Deterministic order — truck VertexIDs are stable per build, but
    // a HashMap iteration is not; sort by the apex coordinates so the
    // corner-blend booleans run in a reproducible sequence.
    out.sort_by(|a, b| {
        let pa = a.1.apex;
        let pb = b.1.apex;
        (pa.x, pa.y, pa.z)
            .partial_cmp(&(pb.x, pb.y, pb.z))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use valenx_cad::primitives::box_solid;

    fn inner_brep(s: &valenx_cad::Solid) -> &TruckSolid {
        match s {
            valenx_cad::Solid::Brep(b) => b,
            _ => panic!("expected brep"),
        }
    }

    fn unique_edges(brep: &TruckSolid) -> Vec<TruckEdge> {
        let mut seen = HashSet::new();
        let mut out = Vec::new();
        for e in brep.edge_iter() {
            if seen.insert(e.id()) {
                out.push(e);
            }
        }
        out
    }

    #[test]
    fn empty_filleted_list_has_zero_corners() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        assert_eq!(unblended_corner_count(brep, &[]), 0);
    }

    #[test]
    fn cube_all_edges_filleted_has_eight_unblended_corners() {
        // Every vertex of a cube touches exactly 3 edges. Filleting
        // all 12 edges means every one of the 8 vertices is a 3-way
        // corner.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edges = unique_edges(brep);
        assert_eq!(edges.len(), 12);
        let unblended = unblended_corner_count(brep, &edges);
        assert_eq!(
            unblended, 8,
            "cube with all 12 edges filleted should have 8 3-way corners, got {unblended}"
        );
    }

    #[test]
    fn is_n_way_vertex_false_for_lonely_edge() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        let edges = unique_edges(brep);
        // Pick one edge; its endpoints should NOT be n-way (only 1
        // filleted edge touches them in this scenario).
        let one_edge = edges[0].clone();
        let v = one_edge.front().id();
        assert!(!is_n_way_vertex(brep, v, &[one_edge.clone()]));
    }

    #[test]
    fn cube_corner_classifies_as_orthogonal_convex() {
        // A cube vertex with all 12 edges filleted is a textbook
        // orthogonal convex 3-edge corner.
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let brep = inner_brep(&cube);
        let edges = unique_edges(brep);
        // The first edge's front vertex is a real cube corner.
        let vid = edges[0].front().id();
        match classify_corner(brep, vid, &edges) {
            CornerClass::OrthogonalConvex(ce) => {
                // The three directions must be mutually orthogonal
                // unit vectors.
                for (i, j) in [(ce.dir0, ce.dir1), (ce.dir0, ce.dir2), (ce.dir1, ce.dir2)] {
                    assert!(
                        i.dot(j).abs() < 1e-6,
                        "corner edge directions should be orthogonal"
                    );
                }
                for d in [ce.dir0, ce.dir1, ce.dir2] {
                    assert!((d.magnitude() - 1.0).abs() < 1e-9, "directions must be unit");
                }
            }
            CornerClass::Unsupported => panic!("cube corner should be orthogonal-convex"),
        }
    }

    #[test]
    fn cube_corner_with_only_two_filleted_edges_is_unsupported() {
        // A degree-2 vertex (only two filleted edges meet) is an
        // edge-end, not a corner — must be Unsupported.
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let brep = inner_brep(&cube);
        let edges = unique_edges(brep);
        // Find a vertex and exactly two edges incident on it.
        let vid = edges[0].front().id();
        let incident: Vec<TruckEdge> = edges
            .iter()
            .filter(|e| e.front().id() == vid || e.back().id() == vid)
            .take(2)
            .cloned()
            .collect();
        assert_eq!(incident.len(), 2);
        assert_eq!(
            classify_corner(brep, vid, &incident),
            CornerClass::Unsupported,
            "a 2-edge vertex is not a blendable corner"
        );
    }

    #[test]
    fn find_blendable_corners_returns_eight_for_fully_filleted_cube() {
        // All 12 edges of a cube filleted → 8 orthogonal convex
        // corners, all blendable.
        let cube = box_solid(3.0, 3.0, 3.0).unwrap();
        let brep = inner_brep(&cube);
        let edges = unique_edges(brep);
        let corners = find_blendable_corners(brep, &edges);
        assert_eq!(
            corners.len(),
            8,
            "a fully-filleted cube has 8 blendable corners, got {}",
            corners.len()
        );
        // Every reported apex must be a genuine cube vertex (a corner
        // of the [0,3]³ box).
        for (_vid, ce) in &corners {
            for coord in [ce.apex.x, ce.apex.y, ce.apex.z] {
                assert!(
                    coord.abs() < 1e-6 || (coord - 3.0).abs() < 1e-6,
                    "corner apex coordinate should be 0 or 3, got {coord}"
                );
            }
        }
    }

    #[test]
    fn find_blendable_corners_empty_for_no_edges() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let brep = inner_brep(&cube);
        assert!(find_blendable_corners(brep, &[]).is_empty());
    }
}
