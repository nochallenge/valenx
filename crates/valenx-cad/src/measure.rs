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

/// The **volume centroid** (centre of mass at uniform density) of a solid,
/// `[x, y, z]` in model units, computed at tessellation tolerance `tol`.
///
/// This is the third mass-property measure promised in the module overview,
/// alongside [`solid_volume_tol`] and [`solid_area_tol`]. It is the
/// divergence-theorem volume centroid summed over the boundary triangulation:
/// each boundary triangle `(p, q, r)` spans a signed tetrahedron with the origin
/// of volume `vₜ/6` (with the scalar triple product `vₜ = p·(q×r)`) whose own
/// centroid is `(p + q + r)/4`, so the solid centroid is their volume-weighted
/// mean `C = Σ vₜ·(p + q + r) / (4·Σ vₜ)`. Like [`solid_volume_tol`] it is exact
/// for a flat-faced solid and converges as `tol → 0` for a curved one, and it is
/// translation-equivariant (moving the solid moves the centroid the same way).
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or if
/// the boundary encloses no volume (a degenerate / empty solid, where the centroid
/// is undefined).
pub fn solid_centroid_tol(solid: &Solid, tol: f64) -> Result<[f64; 3], CadError> {
    let poly = match solid {
        Solid::Brep(b) => brep_polygon(b, tol)?,
        Solid::Mesh(m) => mesh_to_polygon(m),
    };
    let positions = poly.positions();
    // Divergence-theorem volume centroid: accumulate 6× the signed tetrahedron
    // volume `vol6 = Σ p·(q×r)` and the first moment `Σ vₜ·(p+q+r)`; the centroid
    // is the moment over `4·vol6`.
    let mut vol6 = 0.0;
    let (mut mx, mut my, mut mz) = (0.0, 0.0, 0.0);
    for tri in poly.tri_faces() {
        let p = positions[tri[0].pos];
        let q = positions[tri[1].pos];
        let r = positions[tri[2].pos];
        let v = p.x * (q.y * r.z - q.z * r.y) + p.y * (q.z * r.x - q.x * r.z)
            + p.z * (q.x * r.y - q.y * r.x);
        vol6 += v;
        mx += v * (p.x + q.x + r.x);
        my += v * (p.y + q.y + r.y);
        mz += v * (p.z + q.z + r.z);
    }
    for quad in poly.quad_faces() {
        for tri in [[quad[0], quad[1], quad[2]], [quad[0], quad[2], quad[3]]] {
            let p = positions[tri[0].pos];
            let q = positions[tri[1].pos];
            let r = positions[tri[2].pos];
            let v = p.x * (q.y * r.z - q.z * r.y) + p.y * (q.z * r.x - q.x * r.z)
                + p.z * (q.x * r.y - q.y * r.x);
            vol6 += v;
            mx += v * (p.x + q.x + r.x);
            my += v * (p.y + q.y + r.y);
            mz += v * (p.z + q.z + r.z);
        }
    }
    if vol6.abs() < 1e-12 {
        return Err(CadError::Tessellation(
            "solid boundary encloses no volume; centroid undefined".to_string(),
        ));
    }
    Ok([mx / (4.0 * vol6), my / (4.0 * vol6), mz / (4.0 * vol6)])
}

/// Volume centroid of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_centroid_tol`].
pub fn solid_centroid(solid: &Solid) -> Result<[f64; 3], CadError> {
    solid_centroid_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **axis-aligned bounding box** of a solid — `(min, max)` corners `[x, y, z]`
/// in model units — computed at tessellation tolerance `tol`.
///
/// The spatial companion to [`solid_centroid_tol`] / [`solid_volume_tol`]: the
/// tightest box, aligned to the model axes, that encloses the solid's boundary. It
/// is the min and max of every coordinate over the boundary tessellation, so it is
/// exact for a flat-faced solid and converges inward (the inscribed facets cut the
/// corners by `O(tol)`) for a curved one. Useful for view-fitting, broad-phase
/// culling, and as a quick size readout. The box always contains the
/// [`solid_centroid_tol`].
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or if
/// the tessellation yields no boundary vertices (a degenerate / empty solid).
pub fn solid_bounding_box_tol(solid: &Solid, tol: f64) -> Result<([f64; 3], [f64; 3]), CadError> {
    let poly = match solid {
        Solid::Brep(b) => brep_polygon(b, tol)?,
        Solid::Mesh(m) => mesh_to_polygon(m),
    };
    let positions = poly.positions();
    if positions.is_empty() {
        return Err(CadError::Tessellation(
            "solid has no boundary vertices; bounding box undefined".to_string(),
        ));
    }
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for p in positions {
        min[0] = min[0].min(p.x);
        min[1] = min[1].min(p.y);
        min[2] = min[2].min(p.z);
        max[0] = max[0].max(p.x);
        max[1] = max[1].max(p.y);
        max[2] = max[2].max(p.z);
    }
    Ok((min, max))
}

/// Axis-aligned bounding box of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_bounding_box_tol`].
pub fn solid_bounding_box(solid: &Solid) -> Result<([f64; 3], [f64; 3]), CadError> {
    solid_bounding_box_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **centroidal radius of gyration** `k = √(⟨r²⟩ − |c|²)` (model units) of a
/// solid — the root-mean-square distance of the solid's volume from its centre of
/// mass, the length scale that sets its rotational inertia (`I ≈ V·k²` at unit
/// density). Computed at tessellation tolerance `tol`.
///
/// It completes the mass-property suite after [`solid_volume_tol`],
/// [`solid_area_tol`], [`solid_centroid_tol`] and [`solid_bounding_box_tol`]. The
/// divergence-theorem volume integral is extended with the diagonal second moments
/// `∫x², ∫y², ∫z² dV` (the tetrahedron covariance `(V/20)(Σpᵢ² + (Σpᵢ)²)` per axis),
/// then the centroid offset is removed by the parallel-axis theorem. Exact for a
/// flat-faced solid, converging for a curved one. For an `Lx×Ly×Lz` box it is
/// `√((Lx²+Ly²+Lz²)/12)`; for a solid sphere of radius `R`, `√(3/5)·R`.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or if
/// the boundary encloses no volume (a degenerate / empty solid).
pub fn solid_radius_of_gyration_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let poly = match solid {
        Solid::Brep(b) => brep_polygon(b, tol)?,
        Solid::Mesh(m) => mesh_to_polygon(m),
    };
    let positions = poly.positions();
    let mut vol = 0.0;
    let mut first = [0.0_f64; 3]; // ∫x, ∫y, ∫z dV (→ centroid)
    let mut second = [0.0_f64; 3]; // ∫x², ∫y², ∫z² dV about the origin
    for tri in poly.tri_faces() {
        let a = positions[tri[0].pos];
        let b = positions[tri[1].pos];
        let c = positions[tri[2].pos];
        let v = (a.x * (b.y * c.z - b.z * c.y) + a.y * (b.z * c.x - b.x * c.z)
            + a.z * (b.x * c.y - b.y * c.x))
            / 6.0;
        vol += v;
        first[0] += v * (a.x + b.x + c.x) / 4.0;
        first[1] += v * (a.y + b.y + c.y) / 4.0;
        first[2] += v * (a.z + b.z + c.z) / 4.0;
        second[0] += v / 20.0 * (a.x * a.x + b.x * b.x + c.x * c.x + (a.x + b.x + c.x).powi(2));
        second[1] += v / 20.0 * (a.y * a.y + b.y * b.y + c.y * c.y + (a.y + b.y + c.y).powi(2));
        second[2] += v / 20.0 * (a.z * a.z + b.z * b.z + c.z * c.z + (a.z + b.z + c.z).powi(2));
    }
    for quad in poly.quad_faces() {
        for tri in [[quad[0], quad[1], quad[2]], [quad[0], quad[2], quad[3]]] {
            let a = positions[tri[0].pos];
            let b = positions[tri[1].pos];
            let c = positions[tri[2].pos];
            let v = (a.x * (b.y * c.z - b.z * c.y) + a.y * (b.z * c.x - b.x * c.z)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0;
            vol += v;
            first[0] += v * (a.x + b.x + c.x) / 4.0;
            first[1] += v * (a.y + b.y + c.y) / 4.0;
            first[2] += v * (a.z + b.z + c.z) / 4.0;
            second[0] += v / 20.0 * (a.x * a.x + b.x * b.x + c.x * c.x + (a.x + b.x + c.x).powi(2));
            second[1] += v / 20.0 * (a.y * a.y + b.y * b.y + c.y * c.y + (a.y + b.y + c.y).powi(2));
            second[2] += v / 20.0 * (a.z * a.z + b.z * b.z + c.z * c.z + (a.z + b.z + c.z).powi(2));
        }
    }
    if vol.abs() < 1e-12 {
        return Err(CadError::Tessellation(
            "solid encloses no volume; radius of gyration undefined".to_string(),
        ));
    }
    // ⟨r²⟩ about the origin, minus the centroid offset (parallel-axis theorem).
    let mean_sq = (second[0] + second[1] + second[2]) / vol;
    let cen_sq = (first[0] * first[0] + first[1] * first[1] + first[2] * first[2]) / (vol * vol);
    Ok((mean_sq - cen_sq).max(0.0).sqrt())
}

/// Centroidal radius of gyration of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_radius_of_gyration_tol`].
pub fn solid_radius_of_gyration(solid: &Solid) -> Result<f64, CadError> {
    solid_radius_of_gyration_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The three **centroidal principal moments of inertia** `[I₁, I₂, I₃]` (unit
/// density, model units `u⁵`, sorted descending `I₁ ≥ I₂ ≥ I₃`) of a solid — the
/// eigenvalues of the centroidal inertia tensor, whose eigenvectors are the
/// principal axes. Computed at tessellation tolerance `tol`.
///
/// This is the full rotational mass property, completing the suite after
/// [`solid_volume_tol`], [`solid_area_tol`], [`solid_centroid_tol`],
/// [`solid_bounding_box_tol`] and [`solid_radius_of_gyration_tol`]. The
/// divergence-theorem volume integral is extended with the six second moments
/// `∫x², ∫y², ∫z², ∫xy, ∫xz, ∫yz dV` (the tetrahedron covariance
/// `(V/20)(Σpᵢpⱼ + (Σpᵢ)(Σpⱼ))`), every moment is shifted to the centroid by the
/// parallel-axis theorem, the symmetric inertia tensor `I = [[∫y²+z², −∫xy, −∫xz],
/// [−∫xy, ∫x²+z², −∫yz], [−∫xz, −∫yz, ∫x²+y²]]` is assembled, and its symmetric
/// eigendecomposition is taken. For an `Lx×Ly×Lz` box the values are
/// `(V/12)(Lⱼ²+Lₖ²)`; the trace relates to the radius of gyration `k` by
/// `I₁+I₂+I₃ = 2·V·k²`.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or if
/// the boundary encloses no volume (a degenerate / empty solid).
pub fn solid_principal_moments_tol(solid: &Solid, tol: f64) -> Result<[f64; 3], CadError> {
    let poly = match solid {
        Solid::Brep(b) => brep_polygon(b, tol)?,
        Solid::Mesh(m) => mesh_to_polygon(m),
    };
    let positions = poly.positions();
    let mut vol = 0.0;
    let mut first = [0.0_f64; 3]; // ∫x, ∫y, ∫z dV (→ centroid)
    let mut m = [0.0_f64; 6]; // origin second moments: xx, yy, zz, xy, xz, yz
    for tri in poly.tri_faces() {
        let a = positions[tri[0].pos];
        let b = positions[tri[1].pos];
        let c = positions[tri[2].pos];
        let (sx, sy, sz) = (a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z);
        let v = (a.x * (b.y * c.z - b.z * c.y) + a.y * (b.z * c.x - b.x * c.z)
            + a.z * (b.x * c.y - b.y * c.x))
            / 6.0;
        vol += v;
        first[0] += v * sx / 4.0;
        first[1] += v * sy / 4.0;
        first[2] += v * sz / 4.0;
        m[0] += v / 20.0 * (a.x * a.x + b.x * b.x + c.x * c.x + sx * sx);
        m[1] += v / 20.0 * (a.y * a.y + b.y * b.y + c.y * c.y + sy * sy);
        m[2] += v / 20.0 * (a.z * a.z + b.z * b.z + c.z * c.z + sz * sz);
        m[3] += v / 20.0 * (a.x * a.y + b.x * b.y + c.x * c.y + sx * sy);
        m[4] += v / 20.0 * (a.x * a.z + b.x * b.z + c.x * c.z + sx * sz);
        m[5] += v / 20.0 * (a.y * a.z + b.y * b.z + c.y * c.z + sy * sz);
    }
    for quad in poly.quad_faces() {
        for tri in [[quad[0], quad[1], quad[2]], [quad[0], quad[2], quad[3]]] {
            let a = positions[tri[0].pos];
            let b = positions[tri[1].pos];
            let c = positions[tri[2].pos];
            let (sx, sy, sz) = (a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z);
            let v = (a.x * (b.y * c.z - b.z * c.y) + a.y * (b.z * c.x - b.x * c.z)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0;
            vol += v;
            first[0] += v * sx / 4.0;
            first[1] += v * sy / 4.0;
            first[2] += v * sz / 4.0;
            m[0] += v / 20.0 * (a.x * a.x + b.x * b.x + c.x * c.x + sx * sx);
            m[1] += v / 20.0 * (a.y * a.y + b.y * b.y + c.y * c.y + sy * sy);
            m[2] += v / 20.0 * (a.z * a.z + b.z * b.z + c.z * c.z + sz * sz);
            m[3] += v / 20.0 * (a.x * a.y + b.x * b.y + c.x * c.y + sx * sy);
            m[4] += v / 20.0 * (a.x * a.z + b.x * b.z + c.x * c.z + sx * sz);
            m[5] += v / 20.0 * (a.y * a.z + b.y * b.z + c.y * c.z + sy * sz);
        }
    }
    if vol.abs() < 1e-12 {
        return Err(CadError::Tessellation(
            "solid encloses no volume; inertia undefined".to_string(),
        ));
    }
    let (cx, cy, cz) = (first[0] / vol, first[1] / vol, first[2] / vol);
    // Parallel-axis shift of every second moment to the centroid.
    let mxx = m[0] - vol * cx * cx;
    let myy = m[1] - vol * cy * cy;
    let mzz = m[2] - vol * cz * cz;
    let mxy = m[3] - vol * cx * cy;
    let mxz = m[4] - vol * cx * cz;
    let myz = m[5] - vol * cy * cz;
    let tensor = nalgebra::Matrix3::new(
        myy + mzz,
        -mxy,
        -mxz,
        -mxy,
        mxx + mzz,
        -myz,
        -mxz,
        -myz,
        mxx + myy,
    );
    let eig = tensor.symmetric_eigen().eigenvalues;
    let mut vals = [eig[0], eig[1], eig[2]];
    // Sort descending: I₁ ≥ I₂ ≥ I₃.
    vals.sort_by(|p, q| q.partial_cmp(p).unwrap_or(std::cmp::Ordering::Equal));
    Ok(vals)
}

/// Centroidal principal moments of inertia of a solid at
/// [`DEFAULT_MEASURE_TOLERANCE`]. See [`solid_principal_moments_tol`].
pub fn solid_principal_moments(solid: &Solid) -> Result<[f64; 3], CadError> {
    solid_principal_moments_tol(solid, DEFAULT_MEASURE_TOLERANCE)
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
    fn solid_principal_moments_match_the_box_inertia() {
        use crate::primitives::sphere;
        // For an Lx×Ly×Lz box (unit density) the centroidal principal moments are
        // I = (V/12)(Lⱼ²+Lₖ²): with V=48, I_x=4·(16+36)=208, I_y=4·(4+36)=160,
        // I_z=4·(4+16)=80 — exact for a flat-faced solid and independent of the
        // tensor/eigendecomposition path.
        let i = solid_principal_moments(&box_solid(2.0, 4.0, 6.0).unwrap()).unwrap();
        let expect = [208.0, 160.0, 80.0];
        for (got, want) in i.iter().zip(expect.iter()) {
            assert!((got - want).abs() / want < 1e-9, "principal moments {i:?} vs {expect:?}");
        }
        // Sorted descending.
        assert!(i[0] >= i[1] && i[1] >= i[2], "I₁ ≥ I₂ ≥ I₃");

        // THREADS solid_radius_of_gyration + solid_volume: the trace of the inertia
        // tensor is 2·V·⟨r²⟩, so √((I₁+I₂+I₃)/(2V)) = k.
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        let v = solid_volume(&b).unwrap();
        let k_from_inertia = ((i[0] + i[1] + i[2]) / (2.0 * v)).sqrt();
        let k = solid_radius_of_gyration(&b).unwrap();
        assert!((k_from_inertia - k).abs() / k < 1e-9, "√(ΣI/2V) = radius of gyration");

        // A sphere is isotropic: the three principal moments are equal.
        let s = solid_principal_moments(&sphere(2.0).unwrap()).unwrap();
        assert!((s[0] - s[2]).abs() / s[0] < 1e-3, "sphere principal moments equal: {s:?}");
    }

    #[test]
    fn solid_radius_of_gyration_matches_the_box_formula() {
        use crate::primitives::sphere;
        // For an Lx×Ly×Lz box the centroidal radius of gyration is
        // k = √((Lx²+Ly²+Lz²)/12) — exact for a flat-faced solid (textbook formula,
        // independent of the divergence-theorem integral the impl uses).
        let k = solid_radius_of_gyration(&box_solid(2.0, 4.0, 6.0).unwrap()).unwrap();
        let expect = ((2.0_f64 * 2.0 + 4.0 * 4.0 + 6.0 * 6.0) / 12.0).sqrt();
        assert!((k - expect).abs() / expect < 1e-9, "box k = √(56/12) = {expect}, got {k}");

        // A solid sphere of radius R has k = √(3/5)·R about its centre.
        let r = 2.0;
        let ks = solid_radius_of_gyration(&sphere(r).unwrap()).unwrap();
        let expect_s = (3.0_f64 / 5.0).sqrt() * r;
        assert!((ks - expect_s).abs() / expect_s < 1e-3, "sphere k = √(3/5)·R = {expect_s}, got {ks}");

        // k is a length: scaling the box ×2 scales k ×2.
        let k1 = solid_radius_of_gyration(&box_solid(1.0, 1.0, 1.0).unwrap()).unwrap();
        let k2 = solid_radius_of_gyration(&box_solid(2.0, 2.0, 2.0).unwrap()).unwrap();
        assert!((k2 - 2.0 * k1).abs() / (2.0 * k1) < 1e-9, "k ∝ length");
    }

    #[test]
    fn solid_bounding_box_brackets_the_geometry() {
        // A box from the origin to (2,4,6): exact AABB for a flat-faced solid.
        let (min, max) = solid_bounding_box(&box_solid(2.0, 4.0, 6.0).unwrap()).unwrap();
        let (emin, emax) = ([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        for k in 0..3 {
            assert!((min[k] - emin[k]).abs() < 1e-9, "min[{k}] = {} != {}", min[k], emin[k]);
            assert!((max[k] - emax[k]).abs() < 1e-9, "max[{k}] = {} != {}", max[k], emax[k]);
        }
        // Its centre is the box centroid (threads solid_centroid, #261).
        let c = solid_centroid(&box_solid(2.0, 4.0, 6.0).unwrap()).unwrap();
        for k in 0..3 {
            assert!(((min[k] + max[k]) / 2.0 - c[k]).abs() < 1e-9, "bbox centre = centroid at {k}");
        }

        // A cylinder (base disk on the origin in X-Y, axis +Z): z-extent is exactly [0,h].
        let h = 5.0;
        let (cmin, cmax) = solid_bounding_box(&cylinder(1.5, h).unwrap()).unwrap();
        assert!(cmin[2].abs() < 1e-9 && (cmax[2] - h).abs() < 1e-9, "cylinder z ∈ [0, h]");
        // The bounding box must contain the solid's centroid; min ≤ max on every axis.
        let cc = solid_centroid(&cylinder(1.5, h).unwrap()).unwrap();
        for k in 0..3 {
            assert!(cmin[k] <= cmax[k], "min ≤ max at {k}");
            assert!(
                cmin[k] - 1e-9 <= cc[k] && cc[k] <= cmax[k] + 1e-9,
                "centroid inside bbox at {k}"
            );
        }
    }

    #[test]
    fn solid_centroid_matches_analytic_geometry() {
        use crate::primitives::sphere;
        // A box runs from the origin to (dx, dy, dz), so its centroid is the
        // half-diagonal — exact for a flat-faced solid (tessellation is exact).
        let c = solid_centroid(&box_solid(2.0, 4.0, 6.0).unwrap()).unwrap();
        assert!((c[0] - 1.0).abs() < 1e-9, "box cx {} != 1", c[0]);
        assert!((c[1] - 2.0).abs() < 1e-9, "box cy {} != 2", c[1]);
        assert!((c[2] - 3.0).abs() < 1e-9, "box cz {} != 3", c[2]);

        // A cylinder's base disk is centred on the origin in the X-Y plane with its
        // axis along +Z, so the centroid sits on the axis at half height.
        let h = 5.0;
        let cc = solid_centroid(&cylinder(1.5, h).unwrap()).unwrap();
        assert!(cc[0].abs() < 1e-6 && cc[1].abs() < 1e-6, "cylinder off-axis: {cc:?}");
        assert!((cc[2] - h / 2.0).abs() < 1e-6, "cylinder cz {} != h/2", cc[2]);

        // An origin-centred sphere is centroid-symmetric about the origin (a looser
        // bound than the box/cylinder: the revolved-semicircle tessellation is not
        // perfectly pole-symmetric, leaving a ~1e-5 residual on a radius-2 sphere).
        let cs = solid_centroid(&sphere(2.0).unwrap()).unwrap();
        assert!(
            cs[0].abs() < 1e-4 && cs[1].abs() < 1e-4 && cs[2].abs() < 1e-4,
            "sphere centroid not at origin: {cs:?}"
        );

        // SANITY (threads solid_volume — same tessellation): the box encloses 48 u³.
        let v = solid_volume(&box_solid(2.0, 4.0, 6.0).unwrap()).unwrap();
        assert!((v - 48.0).abs() < 1e-9, "box volume {v} != 48");
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


