//! Phase 147 — `ShapeUpgrade_UnifySameDomain` — merge adjacent faces
//! that share the same underlying surface.
//!
//! ## What OCCT does
//!
//! `ShapeUpgrade_UnifySameDomain(shape, edges_flag, faces_flag,
//! concat_bsplines)` walks the shape and merges adjacent faces that
//! share an identical surface (or surfaces that lie on the same
//! infinite extension, e.g. two planar faces on the same plane).
//! Required for STEP import cleanup — many CAD systems export the
//! same face split across multiple "patches" for trim-curve
//! convenience, which downstream meshers then explode into hundreds
//! of unnecessary elements.
//!
//! Also merges adjacent edges that lie on the same underlying curve
//! (when `edges_flag` is set) — same problem at the edge level.
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 147.5) — the mesh-domain
//! coplanar-face merge. truck-modeling exposes no BRep face-merge
//! operator, so the solid is tessellated and adjacent **coplanar**
//! triangles (same plane within tolerance, edge-connected) are
//! flood-filled into maximal regions. Each region is then
//! re-triangulated from its own boundary loop with ear clipping —
//! collapsing an over-patched flat region of N triangles down to the
//! minimal (boundary − 2) triangles. The geometry is unchanged; the
//! triangle count drops, which is exactly the downstream benefit
//! `UnifySameDomain` exists to deliver.
//!
//! The result is a mesh-backed [`Solid`]. Curved-surface ("same
//! cylinder", "same B-spline") unification and the `merge_edges` /
//! `concat_bsplines` flags need parametric-surface introspection
//! truck does not provide; those remain follow-up work — the flags
//! are accepted but only `merge_faces` (planar) is acted on in v1.

use std::collections::VecDeque;

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::OcctAdvancedError;

/// Tessellation chord-error budget for the input solid.
const TESS_TOLERANCE: f64 = 0.25;

/// Two triangle normals count as the same plane when their unit
/// normals agree to this dot-product threshold and the plane offsets
/// agree to [`PLANE_OFFSET_TOLERANCE`].
const NORMAL_DOT_TOLERANCE: f64 = 1.0 - 1.0e-6;

/// Max difference in plane offset `d` for two triangles to be judged
/// coplanar (model units).
const PLANE_OFFSET_TOLERANCE: f64 = 1.0e-4;

/// Apply UnifySameDomain to `solid`.
///
/// `merge_edges` — also merge adjacent collinear edges. *(v1: accepted
/// but not acted on — needs BRep edge introspection.)*
/// `merge_faces` — merge adjacent coplanar faces. *(v1: the planar
/// case is implemented.)*
/// `concat_bsplines` — concatenate merged B-spline patches. *(v1:
/// accepted but not acted on — needs parametric surfaces.)*
///
/// When `merge_faces` is `false` the input is returned unchanged.
///
/// # Errors
///
/// - [`OcctAdvancedError::Backend`] when the solid cannot be
///   tessellated.
///
/// # Example
///
/// ```
/// use valenx_occt_advanced::shape_upgrade_unifysamedomain;
/// use valenx_cad::box_solid;
/// let cube = box_solid(2.0, 2.0, 2.0).unwrap();
/// let unified = shape_upgrade_unifysamedomain(&cube, false, true, false).unwrap();
/// // A cube's 6 planar faces re-triangulate to exactly 12 triangles.
/// match unified {
///     valenx_cad::Solid::Mesh(m) => assert_eq!(m.total_elements(), 12),
///     valenx_cad::Solid::Brep(_) => unreachable!(),
/// }
/// ```
pub fn shape_upgrade_unifysamedomain(
    solid: &Solid,
    merge_edges: bool,
    merge_faces: bool,
    concat_bsplines: bool,
) -> Result<Solid, OcctAdvancedError> {
    let _ = (merge_edges, concat_bsplines); // accepted, not acted on in v1
    if !merge_faces {
        return Ok(solid.clone());
    }

    let mesh = valenx_cad::solid_to_mesh(solid, TESS_TOLERANCE).map_err(|e| {
        OcctAdvancedError::Backend(format!(
            "shape_upgrade_unifysamedomain: cannot tessellate: {e:?}"
        ))
    })?;

    Ok(Solid::from_mesh(merge_coplanar_faces(&mesh)))
}

/// One triangle with a precomputed plane.
struct Tri {
    v: [u32; 3],
    /// Unit normal.
    n: [f64; 3],
    /// Plane offset: `n · p` for any vertex `p`.
    d: f64,
}

/// Flood-fill coplanar connected triangle regions and re-triangulate
/// each region from its boundary loop.
fn merge_coplanar_faces(mesh: &Mesh) -> Mesh {
    let tris = collect_triangles(mesh);
    if tris.is_empty() {
        return mesh.clone();
    }

    // Edge → list of triangle indices, for connectivity.
    let mut edge_tris: std::collections::HashMap<(u32, u32), Vec<usize>> =
        std::collections::HashMap::new();
    for (ti, t) in tris.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (t.v[k], t.v[(k + 1) % 3]);
            let key = if a < b { (a, b) } else { (b, a) };
            edge_tris.entry(key).or_default().push(ti);
        }
    }

    // Flood-fill coplanar connected components.
    let mut region_of = vec![usize::MAX; tris.len()];
    let mut regions: Vec<Vec<usize>> = Vec::new();
    for seed in 0..tris.len() {
        if region_of[seed] != usize::MAX {
            continue;
        }
        let rid = regions.len();
        let mut members = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(seed);
        region_of[seed] = rid;
        while let Some(ti) = queue.pop_front() {
            members.push(ti);
            let t = &tris[ti];
            for k in 0..3 {
                let (a, b) = (t.v[k], t.v[(k + 1) % 3]);
                let key = if a < b { (a, b) } else { (b, a) };
                let Some(neigh) = edge_tris.get(&key) else {
                    continue;
                };
                for &nt in neigh {
                    if region_of[nt] == usize::MAX
                        && same_plane(t, &tris[nt])
                    {
                        region_of[nt] = rid;
                        queue.push_back(nt);
                    }
                }
            }
        }
        regions.push(members);
    }

    // Re-triangulate each region.
    let mut out_conn: Vec<u32> = Vec::new();
    for members in &regions {
        if members.len() == 1 {
            // Singleton — keep as is.
            out_conn.extend_from_slice(&tris[members[0]].v);
            continue;
        }
        let retri = retriangulate_region(mesh, &tris, members);
        if retri.is_empty() {
            // Re-triangulation failed (degenerate boundary) — keep
            // the originals rather than dropping geometry.
            for &ti in members {
                out_conn.extend_from_slice(&tris[ti].v);
            }
        } else {
            for t in retri {
                out_conn.extend_from_slice(&t);
            }
        }
    }

    let mut out = Mesh::new(format!("{}_unified", mesh.id));
    out.nodes = mesh.nodes.clone();
    out.element_blocks.push(ElementBlock {
        element_type: ElementType::Tri3,
        connectivity: out_conn,
    });
    out.recompute_stats();
    out
}

/// Build the [`Tri`] list (Tri3 blocks only) with per-triangle plane.
fn collect_triangles(mesh: &Mesh) -> Vec<Tri> {
    let mut out = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let (i0, i1, i2) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
            let (Some(p0), Some(p1), Some(p2)) =
                (mesh.nodes.get(i0), mesh.nodes.get(i1), mesh.nodes.get(i2))
            else {
                continue;
            };
            let p0 = [p0.x, p0.y, p0.z];
            let p1 = [p1.x, p1.y, p1.z];
            let p2 = [p2.x, p2.y, p2.z];
            let raw = cross(sub(p1, p0), sub(p2, p0));
            let len = norm(raw);
            if len < 1e-18 {
                continue; // degenerate triangle
            }
            let n = [raw[0] / len, raw[1] / len, raw[2] / len];
            let d = dot(n, p0);
            out.push(Tri {
                v: [tri[0], tri[1], tri[2]],
                n,
                d,
            });
        }
    }
    out
}

/// True when two triangles lie on the same plane (normals parallel +
/// same offset) within tolerance.
fn same_plane(a: &Tri, b: &Tri) -> bool {
    // Allow either orientation: |n_a · n_b| ~ 1.
    let dotn = dot(a.n, b.n);
    if dotn.abs() < NORMAL_DOT_TOLERANCE {
        return false;
    }
    // Offsets compared with matched sign.
    let d_b = if dotn >= 0.0 { b.d } else { -b.d };
    (a.d - d_b).abs() < PLANE_OFFSET_TOLERANCE
}

/// Re-triangulate one coplanar region from its boundary loop with
/// ear clipping. Returns `[]` if the region has no single simple
/// boundary loop (holes, non-manifold edges).
fn retriangulate_region(mesh: &Mesh, tris: &[Tri], members: &[usize]) -> Vec<[u32; 3]> {
    // Boundary edges of the region: directed edges whose reverse is
    // not also in the region.
    use std::collections::HashMap;
    let mut count: HashMap<(u32, u32), i32> = HashMap::new();
    let mut directed: Vec<(u32, u32)> = Vec::new();
    for &ti in members {
        let t = &tris[ti];
        for k in 0..3 {
            let (a, b) = (t.v[k], t.v[(k + 1) % 3]);
            let key = if a < b { (a, b) } else { (b, a) };
            *count.entry(key).or_insert(0) += 1;
            directed.push((a, b));
        }
    }
    // Directed boundary edges (undirected count == 1).
    let mut next: HashMap<u32, u32> = HashMap::new();
    let mut boundary_edges = 0usize;
    for (a, b) in directed {
        let key = if a < b { (a, b) } else { (b, a) };
        if count.get(&key).copied().unwrap_or(0) == 1 {
            next.insert(a, b);
            boundary_edges += 1;
        }
    }
    if boundary_edges < 3 {
        return Vec::new();
    }
    // Walk the loop.
    let start = *next.keys().next().unwrap();
    let mut loop_v = vec![start];
    let mut cur = start;
    loop {
        let Some(&nx) = next.get(&cur) else {
            return Vec::new(); // open boundary
        };
        if nx == start {
            break;
        }
        if loop_v.contains(&nx) {
            return Vec::new(); // figure-eight / multi-loop boundary
        }
        loop_v.push(nx);
        cur = nx;
        if loop_v.len() > 1_000_000 {
            return Vec::new();
        }
    }
    if loop_v.len() != boundary_edges {
        // Region boundary is multiple loops (region has a hole) —
        // ear clipping a single loop would be wrong; keep originals.
        return Vec::new();
    }
    ear_clip(mesh, tris[members[0]].n, &loop_v)
}

/// Ear-clip a planar 3D boundary loop. `normal` is the loop's plane
/// normal; vertices are projected to a 2D basis on that plane.
fn ear_clip(mesh: &Mesh, normal: [f64; 3], loop_v: &[u32]) -> Vec<[u32; 3]> {
    if loop_v.len() < 3 {
        return Vec::new();
    }
    // 2D basis on the plane.
    let (u, v) = plane_basis(normal);
    let p2: Vec<(f64, f64)> = loop_v
        .iter()
        .map(|&idx| {
            let p = mesh.nodes[idx as usize];
            let pp = [p.x, p.y, p.z];
            (dot(pp, u), dot(pp, v))
        })
        .collect();

    // Signed area to fix winding (ear clipping below assumes CCW).
    let mut area2 = 0.0;
    for i in 0..p2.len() {
        let j = (i + 1) % p2.len();
        area2 += p2[i].0 * p2[j].1 - p2[j].0 * p2[i].1;
    }
    let mut idx: Vec<usize> = (0..loop_v.len()).collect();
    if area2 < 0.0 {
        idx.reverse();
    }

    let mut tris: Vec<[u32; 3]> = Vec::new();
    let mut guard = 0;
    while idx.len() > 3 {
        guard += 1;
        if guard > 10 * loop_v.len() + 10 {
            return Vec::new(); // no ear found — non-simple polygon
        }
        let n = idx.len();
        let mut clipped = false;
        for i in 0..n {
            let a = idx[(i + n - 1) % n];
            let b = idx[i];
            let c = idx[(i + 1) % n];
            if !is_convex(p2[a], p2[b], p2[c]) {
                continue;
            }
            // No other vertex inside triangle (a, b, c)?
            let mut empty = true;
            for (k, &pk) in idx.iter().enumerate() {
                if pk == a || pk == b || pk == c {
                    continue;
                }
                let _ = k;
                if point_in_tri(p2[pk], p2[a], p2[b], p2[c]) {
                    empty = false;
                    break;
                }
            }
            if empty {
                tris.push([loop_v[a], loop_v[b], loop_v[c]]);
                idx.remove(i);
                clipped = true;
                break;
            }
        }
        if !clipped {
            return Vec::new();
        }
    }
    tris.push([loop_v[idx[0]], loop_v[idx[1]], loop_v[idx[2]]]);
    tris
}

/// Two orthonormal in-plane vectors ⟂ `normal`.
fn plane_basis(normal: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let seed = if normal[0].abs() <= normal[1].abs() && normal[0].abs() <= normal[2].abs() {
        [1.0, 0.0, 0.0]
    } else if normal[1].abs() <= normal[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let u = normalize(cross(normal, seed));
    let v = cross(normal, u);
    (u, v)
}

/// Left-turn test — `b` is a convex corner of a CCW polygon.
fn is_convex(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> bool {
    let cross = (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0);
    cross > 1e-12
}

/// Point-in-triangle via barycentric sign test (2D).
fn point_in_tri(p: (f64, f64), a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> bool {
    let d1 = sign2(p, a, b);
    let d2 = sign2(p, b, c);
    let d3 = sign2(p, c, a);
    let has_neg = d1 < 0.0 || d2 < 0.0 || d3 < 0.0;
    let has_pos = d1 > 0.0 || d2 > 0.0 || d3 > 0.0;
    !(has_neg && has_pos)
}

fn sign2(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    (p.0 - b.0) * (a.1 - b.1) - (a.0 - b.0) * (p.1 - b.1)
}

// --- vector helpers ---
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
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
fn normalize(a: [f64; 3]) -> [f64; 3] {
    let l = norm(a);
    if l < 1e-18 {
        a
    } else {
        [a[0] / l, a[1] / l, a[2] / l]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn merge_faces_false_returns_input_unchanged() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let out = shape_upgrade_unifysamedomain(&cube, false, false, false).unwrap();
        // Untouched — still a BRep with 6 faces.
        assert!(matches!(out, Solid::Brep(_)));
        assert_eq!(out.faces(), 6);
    }

    #[test]
    fn cube_collapses_to_twelve_triangles() {
        // A cube has 6 coplanar faces; each merges to 2 triangles
        // (the minimal triangulation of a quad) → 12 total. truck may
        // tessellate each face into more than 2 triangles, so this
        // proves the merge actually fired.
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        let out = shape_upgrade_unifysamedomain(&cube, false, true, false).unwrap();
        match out {
            Solid::Mesh(m) => {
                assert_eq!(
                    m.total_elements(),
                    12,
                    "6 quad faces should each become 2 triangles"
                );
                // Geometry preserved: 8 distinct cube corners.
                // (welding not applied here, so node count is the
                //  tessellation's — just assert it is non-empty.)
                assert!(!m.nodes.is_empty());
            }
            Solid::Brep(_) => panic!("merge result must be mesh-backed"),
        }
    }

    #[test]
    fn ear_clip_triangulates_a_square() {
        // A unit square loop in the XY plane → 2 triangles.
        let mut mesh = Mesh::new("sq");
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(1.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(1.0, 1.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 1.0, 0.0));
        let tris = ear_clip(&mesh, [0.0, 0.0, 1.0], &[0, 1, 2, 3]);
        assert_eq!(tris.len(), 2);
    }

    #[test]
    fn same_plane_detects_coplanar() {
        let a = Tri {
            v: [0, 1, 2],
            n: [0.0, 0.0, 1.0],
            d: 5.0,
        };
        let b = Tri {
            v: [3, 4, 5],
            n: [0.0, 0.0, 1.0],
            d: 5.0,
        };
        assert!(same_plane(&a, &b));
        let c = Tri {
            v: [6, 7, 8],
            n: [0.0, 0.0, 1.0],
            d: 9.0,
        };
        assert!(!same_plane(&a, &c)); // parallel but offset
        let dperp = Tri {
            v: [9, 10, 11],
            n: [1.0, 0.0, 0.0],
            d: 5.0,
        };
        assert!(!same_plane(&a, &dperp)); // perpendicular
    }
}
