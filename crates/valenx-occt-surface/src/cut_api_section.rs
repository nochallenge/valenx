//! Phase 95 — `BRepAlgoAPI_Section` variant returning curves on both
//! shapes (with pcurves).
//!
//! ## What OCCT does
//!
//! The pcurve-enabled mode of `BRepAlgoAPI_Section`: when configured
//! with `ComputePCurveOn1(true)` and `ComputePCurveOn2(true)`, the
//! operator returns not just the 3D intersection curves but also the
//! 2D parametric (`Geom2d_Curve`) representations on each input
//! surface. Pcurves are essential for re-importing the section into
//! either input as a trim curve.
//!
//! ## v1 status
//!
//! **Honest implementation of the 3D-section core** (Phase 95.5).
//! Both solids are tessellated; every triangle pair is tested with a
//! Möller-style triangle-triangle **segment** intersection (not just
//! the boolean predicate — the actual 3D segment endpoints are
//! computed); the resulting unordered segment soup is chained into
//! ordered polylines (`curves_3d`).
//!
//! **Pcurves are intentionally left empty.** A pcurve is a curve in
//! a surface's `(u, v)` parameter domain — it is only well-defined
//! when the input carries a parametric surface. valenx's `Solid`
//! enum, once tessellated for sectioning, is a triangle mesh with no
//! global parametrization, so there is no honest `(u, v)` to report.
//! Producing real pcurves needs a parametric-BRep section (the
//! Tier-2 `algo_section` BRep variant); until then `pcurves_a` /
//! `pcurves_b` stay empty rather than carrying fabricated values.

use valenx_cad::Solid;

use crate::error::OcctSurfaceError;

/// Output of [`cut_api_section`].
///
/// - `curves_3d` — i-th polyline is the i-th 3D intersection curve
///   in world coordinates.
/// - `pcurves_a` — i-th polyline is the corresponding 2D pcurve on
///   surface `a` in `(u, v)` parameter space.
/// - `pcurves_b` — same for surface `b`.
#[derive(Clone, Debug, Default)]
pub struct SectionWithPcurves {
    /// 3D intersection polyline curves.
    pub curves_3d: Vec<Vec<[f64; 3]>>,
    /// 2D pcurves on surface `a`. Empty for mesh-tessellated inputs —
    /// see the module docs.
    pub pcurves_a: Vec<Vec<[f64; 2]>>,
    /// 2D pcurves on surface `b`. Empty for mesh-tessellated inputs.
    pub pcurves_b: Vec<Vec<[f64; 2]>>,
}

/// Chord tolerance used for the section tessellation.
const SECTION_TESS_TOLERANCE: f64 = 0.25;

/// Points closer than this (world units) are treated as coincident
/// when chaining segments into polylines.
const WELD_TOLERANCE: f64 = 1.0e-6;

/// Both-shapes section: returns the 3D intersection curves where the
/// boundary surfaces of `a` and `b` cross.
///
/// The 3D curves are real geometry, computed by triangle-triangle
/// segment intersection over the tessellated boundaries.
/// `pcurves_a` / `pcurves_b` are empty (see the module docs — a
/// tessellated solid has no parameter domain).
///
/// # Errors
///
/// [`OcctSurfaceError::TruckLimit`] when either solid fails to
/// tessellate.
///
/// # Example
///
/// ```
/// use valenx_occt_surface::cut_api_section;
/// use valenx_cad::box_solid;
/// // Two unit cubes, the second shifted to overlap the first.
/// let a = box_solid(2.0, 2.0, 2.0).unwrap();
/// let b = box_solid(2.0, 2.0, 2.0).unwrap().translated(1.0, 1.0, 1.0).unwrap();
/// let sec = cut_api_section(&a, &b).unwrap();
/// assert!(!sec.curves_3d.is_empty(), "overlapping cubes intersect");
/// ```
pub fn cut_api_section(a: &Solid, b: &Solid) -> Result<SectionWithPcurves, OcctSurfaceError> {
    let mesh_a = valenx_cad::solid_to_mesh(a, SECTION_TESS_TOLERANCE)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("section: tessellate a: {e:?}")))?;
    let mesh_b = valenx_cad::solid_to_mesh(b, SECTION_TESS_TOLERANCE)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("section: tessellate b: {e:?}")))?;

    let tris_a = triangles(&mesh_a);
    let tris_b = triangles(&mesh_b);

    // Collect every triangle-triangle intersection segment.
    let mut segments: Vec<[[f64; 3]; 2]> = Vec::new();
    for ta in &tris_a {
        for tb in &tris_b {
            if let Some(seg) = tri_tri_segment(ta, tb) {
                segments.push(seg);
            }
        }
    }

    let curves_3d = chain_segments(segments);
    Ok(SectionWithPcurves {
        curves_3d,
        pcurves_a: Vec::new(),
        pcurves_b: Vec::new(),
    })
}

/// Extract the triangle list (as world-space corner triples) from a
/// mesh, expanding only `Tri3` blocks.
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
/// do not cross (or are coplanar — coplanar overlap is an area, not a
/// curve, so a section reports nothing for it).
fn tri_tri_segment(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 2]> {
    let n_a = cross(sub(a[1], a[0]), sub(a[2], a[0]));
    let n_b = cross(sub(b[1], b[0]), sub(b[2], b[0]));

    // Signed distances of B's vertices from A's plane and vice-versa.
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

    // Intersection line direction.
    let line_dir = cross(n_a, n_b);
    if norm(line_dir) < 1e-14 {
        // Coplanar — no curve.
        return None;
    }

    // Each triangle crosses the other's plane along a sub-segment of
    // the intersection line. Compute the 3D crossing points.
    let seg_a = plane_crossing_points(a, &da)?;
    let seg_b = plane_crossing_points(b, &db)?;

    // Parametrise all four crossing points along the line direction
    // and take the overlap of the two intervals.
    let base = seg_a[0];
    let proj = |p: [f64; 3]| dot(sub(p, base), line_dir);
    let (a0, a1) = ordered(proj(seg_a[0]), proj(seg_a[1]));
    let (b0, b1) = ordered(proj(seg_b[0]), proj(seg_b[1]));
    let lo = a0.max(b0);
    let hi = a1.min(b1);
    if hi < lo - 1e-12 {
        return None; // intervals disjoint
    }
    let inv = 1.0 / dot(line_dir, line_dir);
    let pt_at = |t: f64| add(base, scale(line_dir, t * inv));
    let p_lo = pt_at(lo);
    let p_hi = pt_at(hi);
    if dist(p_lo, p_hi) < WELD_TOLERANCE {
        return None; // touching at a point, not a segment
    }
    Some([p_lo, p_hi])
}

/// The two 3D points where a triangle's edges cross the plane the
/// signed distances `d` were measured against. Returns `None` for a
/// degenerate crossing (entirely on-plane).
fn plane_crossing_points(t: &[[f64; 3]; 3], d: &[f64; 3]) -> Option<[[f64; 3]; 2]> {
    let mut pts: Vec<[f64; 3]> = Vec::new();
    for k in 0..3 {
        let (i, j) = (k, (k + 1) % 3);
        let (di, dj) = (d[i], d[j]);
        // Edge straddles the plane.
        if (di > 0.0) != (dj > 0.0) {
            let denom = di - dj;
            if denom.abs() < 1e-20 {
                continue;
            }
            let s = di / denom;
            pts.push(add(t[i], scale(sub(t[j], t[i]), s)));
        }
    }
    // On-plane vertices also count as crossing points.
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

/// Chain an unordered segment soup into ordered polylines by greedily
/// linking segments whose endpoints coincide within [`WELD_TOLERANCE`].
fn chain_segments(segments: Vec<[[f64; 3]; 2]>) -> Vec<Vec<[f64; 3]>> {
    let mut remaining: Vec<[[f64; 3]; 2]> = segments;
    let mut polylines: Vec<Vec<[f64; 3]>> = Vec::new();

    while let Some(seed) = remaining.pop() {
        let mut chain: Vec<[f64; 3]> = vec![seed[0], seed[1]];
        // Extend forwards and backwards until nothing connects.
        let mut grew = true;
        while grew {
            grew = false;
            let front = *chain.first().unwrap();
            let back = *chain.last().unwrap();
            let mut found: Option<usize> = None;
            let mut prepend = false;
            let mut new_pt = [0.0; 3];
            for (idx, seg) in remaining.iter().enumerate() {
                if dist(seg[0], back) < WELD_TOLERANCE {
                    found = Some(idx);
                    new_pt = seg[1];
                    prepend = false;
                    break;
                } else if dist(seg[1], back) < WELD_TOLERANCE {
                    found = Some(idx);
                    new_pt = seg[0];
                    prepend = false;
                    break;
                } else if dist(seg[0], front) < WELD_TOLERANCE {
                    found = Some(idx);
                    new_pt = seg[1];
                    prepend = true;
                    break;
                } else if dist(seg[1], front) < WELD_TOLERANCE {
                    found = Some(idx);
                    new_pt = seg[0];
                    prepend = true;
                    break;
                }
            }
            if let Some(idx) = found {
                remaining.swap_remove(idx);
                if prepend {
                    chain.insert(0, new_pt);
                } else {
                    chain.push(new_pt);
                }
                grew = true;
            }
        }
        polylines.push(chain);
    }
    polylines
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
/// True when all three signed distances lie strictly on one side of
/// the plane (no crossing).
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
        // Two boxes whose interiors overlap — their boundary surfaces
        // genuinely cross, so the section curves must be non-empty.
        let a = box_solid(2.0, 2.0, 2.0).unwrap();
        let b = box_solid(2.0, 2.0, 2.0)
            .unwrap()
            .translated(1.0, 1.0, 1.0)
            .unwrap();
        let sec = cut_api_section(&a, &b).unwrap();
        assert!(
            !sec.curves_3d.is_empty(),
            "overlapping cubes must produce section curves"
        );
        // Each polyline has at least two points.
        for c in &sec.curves_3d {
            assert!(c.len() >= 2);
        }
        // Pcurves are honestly empty for mesh-tessellated input.
        assert!(sec.pcurves_a.is_empty());
        assert!(sec.pcurves_b.is_empty());
    }

    #[test]
    fn section_of_disjoint_cubes_is_empty() {
        let a = box_solid(1.0, 1.0, 1.0).unwrap();
        let b = box_solid(1.0, 1.0, 1.0)
            .unwrap()
            .translated(10.0, 10.0, 10.0)
            .unwrap();
        let sec = cut_api_section(&a, &b).unwrap();
        assert!(
            sec.curves_3d.is_empty(),
            "far-apart cubes must not intersect"
        );
    }

    #[test]
    fn tri_tri_segment_finds_crossing() {
        // Triangle A in the z=0 plane, triangle B straddling it
        // vertically — they cross along a segment.
        let a = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [0.0, 4.0, 0.0]];
        let b = [[1.0, 1.0, -1.0], [1.0, 1.0, 1.0], [2.0, 2.0, 0.0]];
        let seg = tri_tri_segment(&a, &b).expect("triangles cross");
        // Both endpoints lie on the z=0 plane.
        assert!(seg[0][2].abs() < 1e-9);
        assert!(seg[1][2].abs() < 1e-9);
    }

    #[test]
    fn tri_tri_segment_none_when_separated() {
        let a = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let b = [[0.0, 0.0, 5.0], [1.0, 0.0, 5.0], [0.0, 1.0, 5.0]];
        assert!(tri_tri_segment(&a, &b).is_none());
    }
}
