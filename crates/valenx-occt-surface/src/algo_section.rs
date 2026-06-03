//! Phase 70 — `BRepAlgoAPI_Section`: intersection curves of two BReps.
//!
//! ## What OCCT does
//!
//! `BRepAlgoAPI_Section` computes the intersection of two `TopoDS_Shape`s
//! (typically two `TopoDS_Solid`s, but any pair of BReps with surface
//! topology) and returns the intersection curves as a `TopoDS_Shape`
//! made of `TopoDS_Edge`s. Each output edge lies on both input shapes;
//! callers usually drive 2D drafting projections or surface-on-surface
//! constraint solving with the result. Options include:
//!
//! - `ComputePCurveOn1(true)` — emit pcurves on the first shape.
//! - `ComputePCurveOn2(true)` — emit pcurves on the second shape.
//! - `Approximation(true)` — fit a B-spline curve to each intersection
//!   polyline rather than returning piecewise-linear edges.
//!
//! ## v1 status — real tessellated-section implementation
//!
//! This is the genuine section operator. Both inputs are tessellated
//! at a fixed chord-error budget; every triangle pair is tested with a
//! Möller-style triangle-triangle **segment** intersection that
//! returns the actual 3D segment endpoints (not just the boolean
//! crossing predicate). The resulting segment soup is the set of
//! section edges.
//!
//! [`algo_section`] returns the *unchained* segment list — one
//! `([f64; 3], [f64; 3])` per intersected triangle pair, the BRep-edge
//! analogue. [`cut_api_section`](crate::cut_api_section()) is the
//! companion that chains the same segments into ordered polylines and
//! is the place to look for the pcurve-aware variant. The B-spline
//! `Approximation` mode is a downstream refinement
//! (`valenx_surface::fit::nurbs_curve`) and is not applied here — the
//! caller can fit the polylines if needed.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// 3D polyline segment representation of an intersection curve.
///
/// Each tuple is `(start_point, end_point)` in world coordinates.
pub type SectionEdge = ([f64; 3], [f64; 3]);

/// Chord tolerance used for the section tessellation.
const SECTION_TESS_TOLERANCE: f64 = 0.25;

/// Compute the intersection edges of two solids.
///
/// Returns a list of `(start, end)` polyline segments — the section
/// where the boundary surfaces of `a` and `b` cross. The list is
/// unordered; chain it with [`cut_api_section`](crate::cut_api_section())
/// (or `valenx_surface` polyline tools) for connected curves.
///
/// # Errors
///
/// [`OcctSurfaceError::TruckLimit`] when either solid fails to
/// tessellate.
///
/// # Example
///
/// ```
/// use valenx_cad::box_solid;
/// use valenx_occt_surface::algo_section;
///
/// // Two boxes whose interiors overlap — their boundaries cross.
/// let a = box_solid(2.0, 2.0, 2.0).unwrap();
/// let b = box_solid(2.0, 2.0, 2.0).unwrap().translated(1.0, 1.0, 1.0).unwrap();
/// let edges = algo_section(&a, &b).unwrap();
/// assert!(!edges.is_empty());
/// ```
pub fn algo_section(a: &Solid, b: &Solid) -> Result<Vec<SectionEdge>, OcctSurfaceError> {
    let mesh_a = valenx_cad::solid_to_mesh(a, SECTION_TESS_TOLERANCE)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("section: tessellate a: {e:?}")))?;
    let mesh_b = valenx_cad::solid_to_mesh(b, SECTION_TESS_TOLERANCE)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("section: tessellate b: {e:?}")))?;

    let tris_a = triangles(&mesh_a);
    let tris_b = triangles(&mesh_b);

    let mut edges: Vec<SectionEdge> = Vec::new();
    for ta in &tris_a {
        for tb in &tris_b {
            if let Some(seg) = tri_tri_segment(ta, tb) {
                edges.push((seg[0], seg[1]));
            }
        }
    }
    Ok(edges)
}

/// Extract the world-space triangle list from a mesh.
fn triangles(mesh: &valenx_mesh::Mesh) -> Vec<[[f64; 3]; 3]> {
    let mut out = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type != valenx_mesh::ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let mut corners = [[0.0; 3]; 3];
            let mut ok = true;
            for (k, &idx) in tri.iter().enumerate() {
                match mesh.nodes.get(idx as usize) {
                    Some(p) => corners[k] = [p.x, p.y, p.z],
                    None => {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                out.push(corners);
            }
        }
    }
    out
}

/// Triangle-triangle intersection **segment** (Möller-style). Returns
/// the 3D endpoints of the shared segment, or `None` if the triangles
/// do not cross transversally.
fn tri_tri_segment(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 2]> {
    let n_a = cross(sub(a[1], a[0]), sub(a[2], a[0]));
    let n_b = cross(sub(b[1], b[0]), sub(b[2], b[0]));

    let d_a = -dot(n_a, a[0]);
    let db = [
        dot(n_a, b[0]) + d_a,
        dot(n_a, b[1]) + d_a,
        dot(n_a, b[2]) + d_a,
    ];
    if all_same_side(&db) {
        return None;
    }
    let d_b = -dot(n_b, b[0]);
    let da = [
        dot(n_b, a[0]) + d_b,
        dot(n_b, a[1]) + d_b,
        dot(n_b, a[2]) + d_b,
    ];
    if all_same_side(&da) {
        return None;
    }

    let line_dir = cross(n_a, n_b);
    if norm(line_dir) < 1e-14 {
        return None; // coplanar — no curve
    }

    let seg_a = plane_crossing_points(a, &da)?;
    let seg_b = plane_crossing_points(b, &db)?;

    let base = seg_a[0];
    let proj = |p: [f64; 3]| dot(sub(p, base), line_dir);
    let (a0, a1) = ordered(proj(seg_a[0]), proj(seg_a[1]));
    let (b0, b1) = ordered(proj(seg_b[0]), proj(seg_b[1]));
    let lo = a0.max(b0);
    let hi = a1.min(b1);
    if hi < lo - 1e-12 {
        return None;
    }
    let inv = 1.0 / dot(line_dir, line_dir);
    let pt_at = |t: f64| add(base, scale(line_dir, t * inv));
    let p_lo = pt_at(lo);
    let p_hi = pt_at(hi);
    if dist(p_lo, p_hi) < 1e-9 {
        return None; // touch at a point, not a segment
    }
    Some([p_lo, p_hi])
}

/// The two 3D points where a triangle's edges cross a plane.
fn plane_crossing_points(t: &[[f64; 3]; 3], d: &[f64; 3]) -> Option<[[f64; 3]; 2]> {
    let mut pts: Vec<[f64; 3]> = Vec::new();
    for k in 0..3 {
        let (i, j) = (k, (k + 1) % 3);
        let (di, dj) = (d[i], d[j]);
        if (di > 0.0) != (dj > 0.0) {
            let denom = di - dj;
            if denom.abs() < 1e-20 {
                continue;
            }
            let s = di / denom;
            pts.push(add(t[i], scale(sub(t[j], t[i]), s)));
        }
    }
    for k in 0..3 {
        if d[k].abs() < 1e-12 {
            pts.push(t[k]);
        }
    }
    if pts.len() < 2 {
        return None;
    }
    Some([pts[0], pts[pts.len() - 1]])
}

// --- vector helpers ---

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm(sub(a, b))
}
fn ordered(x: f64, y: f64) -> (f64, f64) {
    if x <= y {
        (x, y)
    } else {
        (y, x)
    }
}
fn all_same_side(d: &[f64; 3]) -> bool {
    (d[0] > 1e-12 && d[1] > 1e-12 && d[2] > 1e-12)
        || (d[0] < -1e-12 && d[1] < -1e-12 && d[2] < -1e-12)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn section_of_overlapping_cubes_is_non_empty() {
        // Two boxes whose interiors overlap — boundaries cross.
        let a = box_solid(2.0, 2.0, 2.0).unwrap();
        let b = box_solid(2.0, 2.0, 2.0)
            .unwrap()
            .translated(1.0, 1.0, 1.0)
            .unwrap();
        let edges = algo_section(&a, &b).unwrap();
        assert!(!edges.is_empty(), "overlapping cubes intersect");
        // Each edge is a real, non-degenerate segment.
        for (s, e) in &edges {
            let len = dist(*s, *e);
            assert!(len > 1e-9, "section edge is degenerate");
        }
    }

    #[test]
    fn section_of_disjoint_cubes_is_empty() {
        let a = box_solid(1.0, 1.0, 1.0).unwrap();
        let b = box_solid(1.0, 1.0, 1.0)
            .unwrap()
            .translated(10.0, 10.0, 10.0)
            .unwrap();
        let edges = algo_section(&a, &b).unwrap();
        assert!(edges.is_empty(), "far-apart cubes do not intersect");
    }

    #[test]
    fn section_edges_lie_on_both_solids_bounding_region() {
        // The section of two overlapping axis-aligned boxes must lie
        // within the overlap region's bounding box.
        let a = box_solid(4.0, 4.0, 4.0).unwrap();
        let b = box_solid(4.0, 4.0, 4.0)
            .unwrap()
            .translated(2.0, 2.0, 2.0)
            .unwrap();
        let edges = algo_section(&a, &b).unwrap();
        assert!(!edges.is_empty());
        // The overlap box is [2,4]^3.
        for (s, e) in &edges {
            for p in [s, e] {
                for c in 0..3 {
                    assert!(
                        p[c] >= 2.0 - 1e-6 && p[c] <= 4.0 + 1e-6,
                        "section point outside the overlap region: {p:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn tri_tri_segment_finds_crossing() {
        let a = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [0.0, 4.0, 0.0]];
        let b = [[1.0, 1.0, -1.0], [1.0, 1.0, 1.0], [2.0, 2.0, 0.0]];
        let seg = tri_tri_segment(&a, &b).expect("triangles cross");
        assert!(seg[0][2].abs() < 1e-9);
        assert!(seg[1][2].abs() < 1e-9);
    }
}
