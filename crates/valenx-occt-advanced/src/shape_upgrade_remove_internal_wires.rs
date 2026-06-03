//! Phase 149 — `ShapeUpgrade_RemoveInternalWires` — delete wires
//! inside face boundaries.
//!
//! ## What OCCT does
//!
//! A face in OCCT's BRep has one outer wire + zero or more inner
//! ("hole") wires. The inner wires define holes in the face — useful
//! for a face with a real cut-out, defective when the wire is left
//! over from a partial boolean that should have removed it.
//! `ShapeUpgrade_RemoveInternalWires(face_or_shell)` walks every face
//! and discards inner wires whose area falls below a tolerance —
//! cleanly filling in spurious holes.
//!
//! Used in two scenarios: (1) IFC import where coplanar wall faces
//! arrive with leftover annotation wires; (2) post-boolean cleanup
//! where the operator overshoots and leaves a near-zero-area inner
//! wire that fails STEP export.
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 149.5) — the mesh-domain
//! equivalent. truck-modeling has no in-place BRep face mutation, so
//! the solid is tessellated; every open **boundary loop** of the
//! mesh (the mesh-domain analogue of an inner wire / hole) is
//! detected; each loop's enclosed area is measured by the projected
//! shoelace formula; and loops whose area falls **below
//! `area_tolerance`** are sealed by ear-clipping their boundary. A
//! genuine hole (a loop whose area is at or above the tolerance) is
//! deliberately kept — matching OCCT, which only removes the
//! *spurious* sub-tolerance wires.
//!
//! The result is a mesh-backed [`Solid`].

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::OcctAdvancedError;

/// Tessellation chord-error budget for the input solid.
const TESS_TOLERANCE: f64 = 0.25;

/// Walk `solid` and remove inner wires whose enclosed area is below
/// `area_tolerance`.
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for non-positive tolerance.
/// - [`OcctAdvancedError::Backend`] when the solid cannot be
///   tessellated.
///
/// # Example
///
/// ```
/// use valenx_occt_advanced::shape_upgrade_remove_internal_wires;
/// use valenx_cad::box_solid;
/// let cube = box_solid(1.0, 1.0, 1.0).unwrap();
/// // A closed cube has no open boundary loops — it round-trips
/// // through the mesh repair unchanged in triangle count.
/// let cleaned = shape_upgrade_remove_internal_wires(&cube, 1e-3).unwrap();
/// assert!(matches!(cleaned, valenx_cad::Solid::Mesh(_)));
/// ```
pub fn shape_upgrade_remove_internal_wires(
    solid: &Solid,
    area_tolerance: f64,
) -> Result<Solid, OcctAdvancedError> {
    if !area_tolerance.is_finite() || area_tolerance <= 0.0 {
        return Err(OcctAdvancedError::bad_input(
            "area_tolerance",
            "must be positive finite",
        ));
    }

    let mesh = valenx_cad::solid_to_mesh(solid, TESS_TOLERANCE).map_err(|e| {
        OcctAdvancedError::Backend(format!(
            "shape_upgrade_remove_internal_wires: cannot tessellate: {e:?}"
        ))
    })?;

    Ok(Solid::from_mesh(fill_small_loops(&mesh, area_tolerance)))
}

/// Fill every open boundary loop whose enclosed area is below
/// `area_tolerance`; leave larger (genuine-hole) loops open.
fn fill_small_loops(mesh: &Mesh, area_tolerance: f64) -> Mesh {
    let loops = valenx_mesh::boundary_loops(mesh);
    let mut out = mesh.clone();
    out.id = format!("{}_no_internal_wires", mesh.id);

    let mut new_tris: Vec<u32> = Vec::new();
    for lp in &loops {
        if lp.len() < 3 {
            continue;
        }
        let area = loop_area(mesh, lp);
        if area < area_tolerance {
            // Spurious wire — seal it.
            let tris = ear_clip(mesh, lp);
            for t in tris {
                new_tris.extend_from_slice(&t);
            }
        }
        // else: genuine hole — keep it open.
    }

    if !new_tris.is_empty() {
        // Append to the first Tri3 block (or create one).
        if let Some(block) = out
            .element_blocks
            .iter_mut()
            .find(|b| b.element_type == ElementType::Tri3)
        {
            block.connectivity.extend_from_slice(&new_tris);
        } else {
            out.element_blocks.push(ElementBlock {
                element_type: ElementType::Tri3,
                connectivity: new_tris,
            });
        }
    }
    out.recompute_stats();
    out
}

/// Enclosed area of a 3D boundary loop, via Newell's normal + the
/// projected shoelace formula.
fn loop_area(mesh: &Mesh, loop_v: &[u32]) -> f64 {
    // Newell's method gives a vector whose magnitude is twice the
    // planar polygon area (for a planar loop) and a robust average
    // normal for a slightly non-planar one.
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;
    for i in 0..loop_v.len() {
        let a = mesh.nodes[loop_v[i] as usize];
        let b = mesh.nodes[loop_v[(i + 1) % loop_v.len()] as usize];
        nx += (a.y - b.y) * (a.z + b.z);
        ny += (a.z - b.z) * (a.x + b.x);
        nz += (a.x - b.x) * (a.y + b.y);
    }
    0.5 * (nx * nx + ny * ny + nz * nz).sqrt()
}

/// Ear-clip a (near-)planar 3D boundary loop into triangles.
fn ear_clip(mesh: &Mesh, loop_v: &[u32]) -> Vec<[u32; 3]> {
    if loop_v.len() < 3 {
        return Vec::new();
    }
    // Average normal (Newell) → 2D basis.
    let mut nrm = [0.0; 3];
    for i in 0..loop_v.len() {
        let a = mesh.nodes[loop_v[i] as usize];
        let b = mesh.nodes[loop_v[(i + 1) % loop_v.len()] as usize];
        nrm[0] += (a.y - b.y) * (a.z + b.z);
        nrm[1] += (a.z - b.z) * (a.x + b.x);
        nrm[2] += (a.x - b.x) * (a.y + b.y);
    }
    let len = (nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]).sqrt();
    if len < 1e-18 {
        return Vec::new();
    }
    let n = [nrm[0] / len, nrm[1] / len, nrm[2] / len];
    let (u, v) = plane_basis(n);
    let p2: Vec<(f64, f64)> = loop_v
        .iter()
        .map(|&idx| {
            let p = mesh.nodes[idx as usize];
            (
                p.x * u[0] + p.y * u[1] + p.z * u[2],
                p.x * v[0] + p.y * v[1] + p.z * v[2],
            )
        })
        .collect();

    // Winding fix.
    let mut area2 = 0.0;
    for i in 0..p2.len() {
        let j = (i + 1) % p2.len();
        area2 += p2[i].0 * p2[j].1 - p2[j].0 * p2[i].1;
    }
    let mut idx: Vec<usize> = (0..loop_v.len()).collect();
    if area2 < 0.0 {
        idx.reverse();
    }

    let mut tris = Vec::new();
    let mut guard = 0;
    while idx.len() > 3 {
        guard += 1;
        if guard > 10 * loop_v.len() + 10 {
            return Vec::new();
        }
        let m = idx.len();
        let mut clipped = false;
        for i in 0..m {
            let a = idx[(i + m - 1) % m];
            let b = idx[i];
            let c = idx[(i + 1) % m];
            if !is_convex(p2[a], p2[b], p2[c]) {
                continue;
            }
            let mut empty = true;
            for &pk in &idx {
                if pk == a || pk == b || pk == c {
                    continue;
                }
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

fn plane_basis(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let seed = if n[0].abs() <= n[1].abs() && n[0].abs() <= n[2].abs() {
        [1.0, 0.0, 0.0]
    } else if n[1].abs() <= n[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let u = normalize(cross(n, seed));
    let v = cross(n, u);
    (u, v)
}

fn is_convex(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> bool {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0) > 1e-12
}

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

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize(a: [f64; 3]) -> [f64; 3] {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
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
    fn rejects_zero_tolerance() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = shape_upgrade_remove_internal_wires(&cube, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_negative_tolerance() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = shape_upgrade_remove_internal_wires(&cube, -1.0).unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn closed_cube_has_no_loops_to_fill() {
        // A closed cube has zero open boundary loops, so nothing is
        // added; the result is the unmodified tessellation.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let base = valenx_cad::solid_to_mesh(&cube, TESS_TOLERANCE).unwrap();
        let cleaned = shape_upgrade_remove_internal_wires(&cube, 1e-3).unwrap();
        match cleaned {
            Solid::Mesh(m) => assert_eq!(m.total_elements(), base.total_elements()),
            Solid::Brep(_) => panic!("result must be mesh-backed"),
        }
    }

    #[test]
    fn small_loop_is_filled_large_loop_is_kept() {
        // A flat sheet (2 triangles) with a small square hole punched
        // in the middle. With a tolerance above the hole's area the
        // hole gets sealed; below it, the hole stays open.
        let mut mesh = Mesh::new("holed");
        // Outer ring (10x10) corners 0..3, inner hole (1x1) corners 4..7.
        for (x, y) in [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)] {
            mesh.nodes.push(nalgebra::Vector3::new(x, y, 0.0));
        }
        for (x, y) in [(4.0, 4.0), (5.0, 4.0), (5.0, 5.0), (4.0, 5.0)] {
            mesh.nodes.push(nalgebra::Vector3::new(x, y, 0.0));
        }
        // Triangulate the ring-with-hole into 8 triangles (a quad
        // strip around the hole).
        let strip: [[u32; 3]; 8] = [
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ];
        let mut conn = Vec::new();
        for t in strip {
            conn.extend_from_slice(&t);
        }
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: conn,
        });
        mesh.recompute_stats();

        // Hole area = 1.0. Outer boundary area = 100.0.
        // Tolerance 2.0 > 1.0: the hole loop is filled. The outer
        // boundary (area 100) stays open.
        let filled = fill_small_loops(&mesh, 2.0);
        assert!(
            filled.total_elements() > mesh.total_elements(),
            "the sub-tolerance hole must be sealed"
        );

        // Tolerance 0.5 < 1.0: nothing is filled.
        let untouched = fill_small_loops(&mesh, 0.5);
        assert_eq!(
            untouched.total_elements(),
            mesh.total_elements(),
            "a hole above tolerance must be kept"
        );
    }

    #[test]
    fn loop_area_of_unit_square() {
        let mut mesh = Mesh::new("sq");
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(1.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(1.0, 1.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 1.0, 0.0));
        let area = loop_area(&mesh, &[0, 1, 2, 3]);
        assert!((area - 1.0).abs() < 1e-9, "unit square area should be 1.0");
    }
}
