//! Geometric measurement + validity checks for solids.
//!
//! Where [`crate::primitives`] / [`crate::boolean`] *construct* solids,
//! this module *measures* them — volume, surface area, centroid — and
//! checks structural validity (is the boundary a closed 2-manifold?).
//!
//! # Why measurement matters
//!
//! A CAD kernel is only trustworthy if you can prove its output. A
//! commercial kernel (Parasolid, ACIS, OCCT) ships a mass-properties
//! engine precisely so a downstream consumer can assert *"this solid
//! has the volume I expect"*. Valenx's validation suite leans on this
//! module to check primitives, booleans and fillets against their
//! exact analytic ground truth.
//!
//! # Method
//!
//! All three measurements are computed from a tessellation of the
//! solid's boundary, so the result carries the tessellation's
//! discretisation error — a coarser tolerance gives a coarser answer.
//! For *flat-faced* solids (box, prism, any boolean of them) the
//! tessellation is exact and so is the measurement; for *curved*
//! solids (cylinder, sphere, cone, torus, fillets) the answer
//! converges to the true value as the tolerance shrinks.
//!
//! - **Volume** — the divergence-theorem surface integral
//!   `∭ dV = ∯ x dy dz`, evaluated per triangle by
//!   [`truck_meshalgo`]'s `CalcVolume`. Exact for a closed flat-faced
//!   mesh; converges from *below* for a convex curved solid (the
//!   inscribed facets cut the corners).
//! - **Surface area** — the sum of every boundary triangle's area
//!   (`½‖(q−p)×(r−p)‖`). Converges from *below* for a convex curved
//!   surface, same reason.
//! - **Closed-solid check** — the boundary tessellation's vertices are
//!   *welded* (truck tessellates each face with an independent vertex
//!   array — see [`is_closed_solid_tol`]), zero-area pole slivers are
//!   dropped, and every remaining triangle edge is counted; a valid
//!   solid is a closed orientable 2-manifold, so each directed edge is
//!   used exactly once and its reverse exactly once.
//!
//! Mesh-backed solids ([`Solid::Mesh`]) are measured directly off
//! their cached triangle mesh — no BRep tessellation step.

use truck_meshalgo::analyzers::CalcVolume;
use truck_meshalgo::prelude::MeshableShape;
use truck_meshalgo::prelude::PolygonMesh;
use truck_meshalgo::tessellation::MeshedShape;

use crate::solid::{CadError, Solid};

/// Tessellation tolerance used by [`solid_volume`] / [`solid_area`] /
/// [`is_closed_solid`] when the caller doesn't supply one.
///
/// Tight enough (`1e-3` model units) that a unit-scale curved solid
/// measures to ~4 significant figures, which is what the validation
/// suite asserts against.
pub const DEFAULT_MEASURE_TOLERANCE: f64 = 1.0e-3;

/// Tessellate a BRep solid to a single merged [`PolygonMesh`].
///
/// Crate-internal helper shared by every measurement function. Returns
/// [`CadError::Tessellation`] if the tolerance is not strictly
/// positive and finite.
fn brep_polygon(brep: &truck_modeling::Solid, tol: f64) -> Result<PolygonMesh, CadError> {
    if !tol.is_finite() || tol <= 0.0 {
        return Err(CadError::Tessellation(format!(
            "measurement tolerance must be finite and > 0, got {tol}"
        )));
    }
    Ok(brep.triangulation(tol).to_polygon())
}

/// Build a [`PolygonMesh`] of a mesh-backed solid's cached triangles.
///
/// truck's `CalcVolume` / `Topology` traits are implemented on
/// `PolygonMesh`, so a [`Solid::Mesh`] is converted into one to share
/// the same measurement code path as a tessellated BRep.
fn mesh_to_polygon(mesh: &valenx_mesh::Mesh) -> PolygonMesh {
    use truck_meshalgo::prelude::{Faces, StandardAttributes, StandardVertex};
    use truck_modeling::Point3;

    let positions: Vec<Point3> = mesh
        .nodes
        .iter()
        .map(|n| Point3::new(n.x, n.y, n.z))
        .collect();
    let mut tri_faces: Vec<[StandardVertex; 3]> = Vec::new();
    for block in &mesh.element_blocks {
        // Only triangle blocks contribute a boundary surface; the CAD
        // pipeline only ever emits Tri3, but guard the stride anyway.
        let stride = block.element_type.nodes_per_element();
        if stride != 3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            tri_faces.push([
                StandardVertex::from(tri[0] as usize),
                StandardVertex::from(tri[1] as usize),
                StandardVertex::from(tri[2] as usize),
            ]);
        }
    }
    let attrs = StandardAttributes {
        positions,
        ..Default::default()
    };
    PolygonMesh::new(attrs, Faces::from_tri_and_quad_faces(tri_faces, Vec::new()))
}

/// Signed volume of a solid, in cubic model units.
///
/// The result is *signed* by the boundary orientation: a correctly
/// built solid (outward-facing normals) returns a **positive** volume;
/// a negative value means the boundary is inside-out.
///
/// # Accuracy
///
/// - Flat-faced solids (box, prism, booleans of them) — exact.
/// - Curved solids (cylinder, sphere, cone, torus) — converges to the
///   true volume from *below* as `tol → 0`; the inscribed triangle
///   facets undershoot a convex curved boundary.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly
/// positive.
pub fn solid_volume_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let poly = match solid {
        Solid::Brep(b) => brep_polygon(b, tol)?,
        Solid::Mesh(m) => mesh_to_polygon(m),
    };
    Ok(poly.volume())
}

/// Volume of a solid at [`DEFAULT_MEASURE_TOLERANCE`].
pub fn solid_volume(solid: &Solid) -> Result<f64, CadError> {
    solid_volume_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// Total boundary surface area of a solid, in square model units.
///
/// Computed as the sum of every boundary triangle's area. Like
/// [`solid_volume_tol`] this converges from below for curved solids
/// and is exact for flat-faced ones.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly
/// positive.
pub fn solid_area_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let poly = match solid {
        Solid::Brep(b) => brep_polygon(b, tol)?,
        Solid::Mesh(m) => mesh_to_polygon(m),
    };
    let positions = poly.positions();
    let mut area = 0.0;
    for tri in poly.tri_faces() {
        let p = positions[tri[0].pos];
        let q = positions[tri[1].pos];
        let r = positions[tri[2].pos];
        // ½‖(q−p)×(r−p)‖.
        let u = q - p;
        let v = r - p;
        let cross = truck_modeling::Vector3::new(
            u.y * v.z - u.z * v.y,
            u.z * v.x - u.x * v.z,
            u.x * v.y - u.y * v.x,
        );
        area += 0.5 * truck_modeling::InnerSpace::magnitude(cross);
    }
    // Quads: tessellation rarely emits them for a solid boundary, but
    // split any that appear into two triangles.
    for quad in poly.quad_faces() {
        for tri in [[quad[0], quad[1], quad[2]], [quad[0], quad[2], quad[3]]] {
            let p = positions[tri[0].pos];
            let q = positions[tri[1].pos];
            let r = positions[tri[2].pos];
            let u = q - p;
            let v = r - p;
            let cross = truck_modeling::Vector3::new(
                u.y * v.z - u.z * v.y,
                u.z * v.x - u.x * v.z,
                u.x * v.y - u.y * v.x,
            );
            area += 0.5 * truck_modeling::InnerSpace::magnitude(cross);
        }
    }
    Ok(area)
}

/// Surface area of a solid at [`DEFAULT_MEASURE_TOLERANCE`].
pub fn solid_area(solid: &Solid) -> Result<f64, CadError> {
    solid_area_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// Weld a polygon mesh's triangle indices onto a deduplicated vertex
/// set, so coincident positions collapse to one index.
///
/// truck's BRep tessellator emits one *independent* vertex array per
/// face — vertices on a shared BRep edge appear once per adjacent
/// face, never index-merged. [`Topology::shell_condition`] matches
/// boundary edges by vertex-index pairs, so on an un-welded mesh every
/// face looks like an isolated island with all four edges "open". This
/// helper snaps positions within `weld_tol` to a common index so the
/// shell-condition test sees the true topology.
///
/// Returns the welded triangle list as index triples into a
/// deduplicated `Vec<Point3>` (also returned).
fn weld_triangles(
    poly: &PolygonMesh,
    weld_tol: f64,
) -> (Vec<truck_modeling::Point3>, Vec<[usize; 3]>) {
    use std::collections::HashMap;
    use truck_modeling::Point3;

    let positions = poly.positions();
    let h = weld_tol.max(1e-30);
    let key = |p: &Point3| -> (i64, i64, i64) {
        (
            (p.x / h).round() as i64,
            (p.y / h).round() as i64,
            (p.z / h).round() as i64,
        )
    };
    let mut buckets: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    let mut deduped: Vec<Point3> = Vec::new();
    let mut remap: Vec<usize> = Vec::with_capacity(positions.len());
    let tol_sq = weld_tol * weld_tol;

    for p in positions.iter() {
        let center = key(p);
        let mut found: Option<usize> = None;
        'outer: for dx in -1..=1i64 {
            for dy in -1..=1i64 {
                for dz in -1..=1i64 {
                    let nb = (center.0 + dx, center.1 + dy, center.2 + dz);
                    if let Some(cands) = buckets.get(&nb) {
                        for &ci in cands {
                            let q = &deduped[ci];
                            let d = (p.x - q.x).powi(2)
                                + (p.y - q.y).powi(2)
                                + (p.z - q.z).powi(2);
                            if d <= tol_sq {
                                found = Some(ci);
                                break 'outer;
                            }
                        }
                    }
                }
            }
        }
        let idx = match found {
            Some(i) => i,
            None => {
                let i = deduped.len();
                deduped.push(*p);
                buckets.entry(center).or_default().push(i);
                i
            }
        };
        remap.push(idx);
    }

    let mut tris: Vec<[usize; 3]> = Vec::new();
    for tri in poly.tri_faces() {
        tris.push([remap[tri[0].pos], remap[tri[1].pos], remap[tri[2].pos]]);
    }
    for quad in poly.quad_faces() {
        let q = [
            remap[quad[0].pos],
            remap[quad[1].pos],
            remap[quad[2].pos],
            remap[quad[3].pos],
        ];
        tris.push([q[0], q[1], q[2]]);
        tris.push([q[0], q[2], q[3]]);
    }
    (deduped, tris)
}

/// Whether the solid's boundary is a **closed 2-manifold** — the
/// defining property of a valid solid.
///
/// A closed solid's boundary tessellates to a mesh where every edge is
/// shared by exactly two triangles with consistent winding. An open
/// shell, a non-manifold edge, or a flipped face all fail this check.
///
/// This is the cheapest honest "is this a valid solid?" test the
/// kernel can offer without a full BRep-topology audit: a boolean that
/// silently produced a self-intersecting or non-watertight result is
/// caught here.
///
/// # Method
///
/// The solid's boundary is tessellated and the triangle vertices are
/// **welded** (truck's per-face tessellation does not share vertices
/// across faces — see `weld_triangles`). **Degenerate triangles are
/// dropped first**: at a singular pole (a sphere's poles, a cone's
/// apex) truck's tessellator emits a handful of zero-area slivers
/// whose vertices weld to two distinct positions — they are
/// tessellation artifacts of the surface singularity, not part of the
/// boundary, and a pole that fans to one vertex is still a valid
/// manifold. After they are dropped, every remaining triangle edge is
/// counted: a closed orientable 2-manifold uses each directed edge
/// exactly once and its reverse exactly once. A boundary edge (open
/// shell), a non-manifold edge, or a flipped face all fail the check.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly
/// positive.
pub fn is_closed_solid_tol(solid: &Solid, tol: f64) -> Result<bool, CadError> {
    use std::collections::HashMap;

    let poly = match solid {
        Solid::Brep(b) => brep_polygon(b, tol)?,
        Solid::Mesh(m) => mesh_to_polygon(m),
    };
    // Weld within a fraction of the tessellation tolerance — tight
    // enough not to merge genuinely-distinct geometry, loose enough to
    // collapse the per-face duplicate vertices the tessellator emits.
    let weld_tol = (tol * 0.25).max(1e-9);
    let (_pts, tris) = weld_triangles(&poly, weld_tol);
    // Count each directed edge over the non-degenerate triangles. A
    // closed orientable manifold uses every undirected edge exactly
    // twice — once `(a,b)`, once `(b,a)`.
    let mut directed: HashMap<(usize, usize), i32> = HashMap::new();
    let mut real_tris = 0usize;
    for t in &tris {
        // Degenerate triangle (two welded vertices coincide) — a
        // zero-area pole sliver. Skip it; it contributes no boundary.
        if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] {
            continue;
        }
        real_tris += 1;
        for &(a, b) in &[(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            *directed.entry((a, b)).or_insert(0) += 1;
        }
    }
    if real_tris == 0 {
        return Ok(false);
    }
    for (&(a, b), &count) in &directed {
        // Each directed edge must appear exactly once, and its reverse
        // exactly once. count != 1 ⇒ a flipped face or a non-manifold
        // edge; missing reverse ⇒ an open boundary edge.
        if count != 1 {
            return Ok(false);
        }
        if directed.get(&(b, a)).copied().unwrap_or(0) != 1 {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether the solid is a closed 2-manifold, at
/// [`DEFAULT_MEASURE_TOLERANCE`].
pub fn is_closed_solid(solid: &Solid) -> Result<bool, CadError> {
    is_closed_solid_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The Euler characteristic `χ = V − E + F` of a BRep solid's
/// topology.
///
/// For a solid whose boundary is a single sphere-topology shell
/// (genus 0 — a box, cylinder, sphere, cone, any simply-connected
/// solid) the Euler–Poincaré formula gives `χ = 2`. A genus-`g`
/// solid (a torus is genus 1) gives `χ = 2 − 2g`.
///
/// Returns `None` for a mesh-backed solid — triangle counts are not
/// BRep topology and the relation would not hold.
pub fn euler_characteristic(solid: &Solid) -> Option<i64> {
    match solid {
        Solid::Brep(_) => {
            let v = solid.vertices() as i64;
            let e = solid.edges() as i64;
            let f = solid.faces() as i64;
            Some(v - e + f)
        }
        Solid::Mesh(_) => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{box_solid, cylinder, prism};
    use std::f64::consts::PI;

    #[test]
    fn unit_cube_volume_is_exact() {
        // A box is flat-faced — its tessellation is exact, so the
        // measured volume must hit 1.0 to machine precision.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let v = solid_volume(&cube).unwrap();
        assert!((v - 1.0).abs() < 1e-9, "unit cube volume {v} != 1.0");
    }

    #[test]
    fn box_volume_is_product_of_dims() {
        let b = box_solid(2.0, 3.0, 4.0).unwrap();
        let v = solid_volume(&b).unwrap();
        assert!((v - 24.0).abs() < 1e-9, "2×3×4 box volume {v} != 24");
    }

    #[test]
    fn unit_cube_area_is_exact() {
        // 6 unit faces → area 6.0, exact for a flat-faced solid.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let a = solid_area(&cube).unwrap();
        assert!((a - 6.0).abs() < 1e-9, "unit cube area {a} != 6.0");
    }

    #[test]
    fn box_is_a_closed_solid() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        assert!(is_closed_solid(&cube).unwrap(), "a box must be closed");
    }

    #[test]
    fn box_euler_characteristic_is_two() {
        // A box is genus-0: V−E+F = 8−12+6 = 2.
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        assert_eq!(euler_characteristic(&cube), Some(2));
    }

    #[test]
    fn cylinder_volume_converges_to_pi_r2_h() {
        // πr²h for r=1, h=2 → 2π. A curved solid tessellates from
        // below, so a fine tolerance must land within ~1% under.
        let cyl = cylinder(1.0, 2.0).unwrap();
        let v = solid_volume_tol(&cyl, 1e-3).unwrap();
        let exact = PI * 1.0 * 1.0 * 2.0;
        assert!(
            v > 0.0 && v <= exact + 1e-6,
            "cylinder volume {v} should not exceed exact {exact}"
        );
        assert!(
            (exact - v) / exact < 0.02,
            "cylinder volume {v} should be within 2% of {exact}"
        );
    }

    #[test]
    fn prism_volume_equals_base_area_times_height() {
        // A unit right-triangle prism: base area ½, height 3 → 1.5.
        let tri = prism(&[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], 3.0).unwrap();
        let v = solid_volume(&tri).unwrap();
        assert!((v - 1.5).abs() < 1e-9, "triangle prism volume {v} != 1.5");
    }

    #[test]
    fn measurement_rejects_bad_tolerance() {
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        assert!(matches!(
            solid_volume_tol(&cube, 0.0),
            Err(CadError::Tessellation(_))
        ));
        assert!(matches!(
            solid_area_tol(&cube, -1.0),
            Err(CadError::Tessellation(_))
        ));
        assert!(matches!(
            is_closed_solid_tol(&cube, f64::NAN),
            Err(CadError::Tessellation(_))
        ));
    }

    #[test]
    fn mesh_backed_solid_volume_uses_cached_mesh() {
        // A mesh-backed solid is measured off its triangles directly.
        // Build a closed unit tetrahedron mesh and confirm its volume.
        use valenx_mesh::{ElementBlock, ElementType, Mesh};
        let mut mesh = Mesh::new("tet");
        // A tetrahedron with one corner at origin and unit legs has
        // volume 1/6.
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(1.0, 0.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 1.0, 0.0));
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 0.0, 1.0));
        let mut block = ElementBlock::new(ElementType::Tri3);
        // Outward-facing windings for the 4 faces.
        block.connectivity = vec![
            0, 2, 1, // bottom (z=0), normal -z
            0, 1, 3, // front (y=0), normal -y
            0, 3, 2, // left (x=0), normal -x
            1, 2, 3, // slanted hypotenuse face
        ];
        mesh.element_blocks.push(block);
        mesh.recompute_stats();
        let solid = Solid::from_mesh(mesh);
        let v = solid_volume(&solid).unwrap();
        assert!(
            (v - 1.0 / 6.0).abs() < 1e-9,
            "unit tetrahedron volume {v} != 1/6"
        );
        assert!(
            is_closed_solid(&solid).unwrap(),
            "the tetrahedron mesh is closed"
        );
    }
}


