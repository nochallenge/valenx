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

/// The **volume-equivalent sphere radius** `r_eq = (3·V / 4π)^(1/3)` (model units) of a
/// solid — the radius of the sphere that has the same volume `V` (from
/// [`solid_volume_tol`]), evaluated at tessellation tolerance `tol`. It is the standard
/// **equivalent spherical diameter** size used to report particle / grain size from
/// sieving, sedimentation, or laser diffraction, and it equals the actual radius for a
/// sphere. Unlike the dimensionless [`solid_sphericity_tol`] (a shape ratio) this is a
/// length (a size).
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid
/// volume cannot be computed.
pub fn solid_equivalent_sphere_radius_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let v = solid_volume_tol(solid, tol)?;
    Ok((3.0 * v / (4.0 * std::f64::consts::PI)).cbrt())
}

/// Volume-equivalent sphere radius of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_equivalent_sphere_radius_tol`].
pub fn solid_equivalent_sphere_radius(solid: &Solid) -> Result<f64, CadError> {
    solid_equivalent_sphere_radius_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **volume-equivalent cube edge** `s = V^(1/3)` (model units) of a solid — the edge
/// length of the cube that has the same volume `V` (from [`solid_volume_tol`]), evaluated
/// at tessellation tolerance `tol`. It is the cube companion to the equivalent-sphere
/// diameter ([`solid_equivalent_sphere_radius_tol`]) — another standard nominal "size" of
/// an irregular body — and equals the actual edge for a cube.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid
/// volume cannot be computed.
pub fn solid_equivalent_cube_edge_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    Ok(solid_volume_tol(solid, tol)?.cbrt())
}

/// Volume-equivalent cube edge of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_equivalent_cube_edge_tol`].
pub fn solid_equivalent_cube_edge(solid: &Solid) -> Result<f64, CadError> {
    solid_equivalent_cube_edge_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **volume-equivalent sphere surface area** `A_eq = 4π·r_eq²` (model units²) of a
/// solid — the surface area of the sphere that has the same volume, with the
/// equivalent-sphere radius `r_eq` from [`solid_equivalent_sphere_radius_tol`], evaluated
/// at tessellation tolerance `tol`. By the isoperimetric inequality it is the **minimum
/// possible surface area** for the solid's volume, so `A_eq ≤ A` always; their ratio is
/// exactly the [`solid_sphericity_tol`]. Unlike that dimensionless sphericity this is an
/// area — the least surface a given volume can present (the diffusion / heat-transfer
/// floor).
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid
/// volume cannot be computed.
pub fn solid_equivalent_sphere_area_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let r = solid_equivalent_sphere_radius_tol(solid, tol)?;
    Ok(4.0 * std::f64::consts::PI * r * r)
}

/// Volume-equivalent sphere surface area of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_equivalent_sphere_area_tol`].
pub fn solid_equivalent_sphere_area(solid: &Solid) -> Result<f64, CadError> {
    solid_equivalent_sphere_area_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **equivalent spherical diameter** `d_eq = 2·r_eq` (model units) of a solid — the
/// diameter of the sphere with the same volume, evaluated at tessellation tolerance `tol`
/// from the [`solid_equivalent_sphere_radius_tol`]. It is the standard reported particle
/// size in particle sizing, powder metallurgy and sedimentology (the "ESD"). Its square
/// times π is the volume-equivalent sphere surface area [`solid_equivalent_sphere_area_tol`]
/// (`A = π·d²`).
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid
/// encloses no volume.
pub fn solid_equivalent_sphere_diameter_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    Ok(2.0 * solid_equivalent_sphere_radius_tol(solid, tol)?)
}

/// Equivalent spherical diameter of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_equivalent_sphere_diameter_tol`].
pub fn solid_equivalent_sphere_diameter(solid: &Solid) -> Result<f64, CadError> {
    solid_equivalent_sphere_diameter_tol(solid, DEFAULT_MEASURE_TOLERANCE)
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

/// The **Wadell sphericity** `Ψ = π^(1/3)·(6·V)^(2/3) / A` (dimensionless) of a solid
/// — the surface area of a sphere of the same volume divided by the solid's own
/// surface area, the standard measure of how *sphere-like* a shape is. It is exactly
/// `1` for a sphere (which minimises surface area for a given volume) and strictly
/// less than `1` for every other shape; it is widely used in sedimentation, powder
/// technology and granular-flow modelling. Computed at tessellation tolerance `tol`
/// from [`solid_volume_tol`] and [`solid_area_tol`].
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the
/// solid has no surface area (a degenerate / empty solid).
pub fn solid_sphericity_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let v = solid_volume_tol(solid, tol)?;
    let a = solid_area_tol(solid, tol)?;
    if a <= 0.0 {
        return Err(CadError::Tessellation(
            "solid has no surface area; sphericity undefined".to_string(),
        ));
    }
    Ok(std::f64::consts::PI.cbrt() * (6.0 * v).powf(2.0 / 3.0) / a)
}

/// Wadell sphericity of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_sphericity_tol`].
pub fn solid_sphericity(solid: &Solid) -> Result<f64, CadError> {
    solid_sphericity_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **specific surface area** `A/V` (units of inverse length) of a solid — its
/// total surface area per unit volume, computed at tessellation tolerance `tol` from
/// [`solid_area_tol`] and [`solid_volume_tol`]. It is `3/R` for a sphere of radius
/// `R`, `6/a` for a cube of side `a`, and in general decreases as `1/length` as a
/// shape is scaled up.
///
/// Unlike the dimensionless [`solid_sphericity_tol`], this is a *dimensional* measure
/// that captures size as well as shape: it is the governing quantity for surface-rate
/// processes — convective and radiative heat transfer, catalytic and dissolution
/// reaction rates, drying, and biological exchange — where flux scales with the
/// exposed area but capacity with the enclosed volume.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or if the
/// solid encloses no volume (a degenerate / empty solid, where `A/V` is undefined).
pub fn solid_specific_surface_area_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let v = solid_volume_tol(solid, tol)?;
    let a = solid_area_tol(solid, tol)?;
    if v <= 0.0 {
        return Err(CadError::Tessellation(
            "solid encloses no volume; specific surface area undefined".to_string(),
        ));
    }
    Ok(a / v)
}

/// Specific surface area `A/V` of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_specific_surface_area_tol`].
pub fn solid_specific_surface_area(solid: &Solid) -> Result<f64, CadError> {
    solid_specific_surface_area_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **mean chord length** `⟨ℓ⟩ = 4·V / A` (model units) of a solid — the average length
/// of a uniformly-random straight chord through it, evaluated at tessellation tolerance
/// `tol`. This is Cauchy's classic result for a convex body (the mean free path / mean
/// linear intercept of stereology and radiation transport): a sphere of radius `r` has mean
/// chord exactly `4r/3`. It is the reciprocal-scaled companion of the
/// [`solid_specific_surface_area_tol`] (`A/V`), since `⟨ℓ⟩ = 4 / (A/V)`.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, the solid is
/// empty, or its surface area is non-positive (degenerate).
pub fn solid_mean_chord_length_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let v = solid_volume_tol(solid, tol)?;
    let a = solid_area_tol(solid, tol)?;
    if a <= 0.0 {
        return Err(CadError::Tessellation(
            "degenerate solid (non-positive surface area); mean chord undefined".to_string(),
        ));
    }
    Ok(4.0 * v / a)
}

/// Mean chord length `4·V/A` of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_mean_chord_length_tol`].
pub fn solid_mean_chord_length(solid: &Solid) -> Result<f64, CadError> {
    solid_mean_chord_length_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **Sauter mean diameter** `d₃₂ = 6·V / A` (model units) of a solid — the diameter
/// of a sphere with the same volume-to-surface ratio, evaluated at tessellation tolerance
/// `tol`. It is the standard particle-sizing diameter in sprays, aerosols, powders, and
/// packed beds (the interfacial-area-weighted mean), and the surface-to-volume companion
/// to the Cauchy mean chord [`solid_mean_chord_length_tol`] (`4V/A`): `d₃₂ = 1.5·⟨ℓ⟩`. For
/// a sphere it returns the true diameter `2r`.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid is
/// degenerate (non-positive surface area).
pub fn solid_sauter_mean_diameter_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let v = solid_volume_tol(solid, tol)?;
    let a = solid_area_tol(solid, tol)?;
    if a <= 0.0 {
        return Err(CadError::Tessellation(
            "degenerate solid (non-positive surface area); Sauter mean diameter undefined"
                .to_string(),
        ));
    }
    Ok(6.0 * v / a)
}

/// Sauter mean diameter of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_sauter_mean_diameter_tol`].
pub fn solid_sauter_mean_diameter(solid: &Solid) -> Result<f64, CadError> {
    solid_sauter_mean_diameter_tol(solid, DEFAULT_MEASURE_TOLERANCE)
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

/// The **axis-aligned bounding-box center** `c = (min + max) / 2` (model units) of a solid
/// — the geometric mid-point of its AABB ([`solid_bounding_box_tol`]), evaluated at
/// tessellation tolerance `tol`. It is the natural pivot / camera-framing point. Unlike the
/// mass [`solid_centroid_tol`] it weights only the extremes of the geometry, so the two
/// coincide for a centro-symmetric body (a box, a sphere) but differ for a lopsided one
/// (the centroid leans toward the bulk, the box center toward the extent).
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid
/// is empty (no bounding box).
pub fn solid_bounding_box_center_tol(solid: &Solid, tol: f64) -> Result<[f64; 3], CadError> {
    let (min, max) = solid_bounding_box_tol(solid, tol)?;
    Ok([
        (min[0] + max[0]) * 0.5,
        (min[1] + max[1]) * 0.5,
        (min[2] + max[2]) * 0.5,
    ])
}

/// AABB center of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_bounding_box_center_tol`].
pub fn solid_bounding_box_center(solid: &Solid) -> Result<[f64; 3], CadError> {
    solid_bounding_box_center_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **axis-aligned bounding-box extents** `[Δx, Δy, Δz] = max − min` (model units) of a
/// solid — the per-axis width, height and depth of its AABB, evaluated at tessellation
/// tolerance `tol` from [`solid_bounding_box_tol`]. These are the raw dimensions the rest of
/// the AABB family is built from: the [`solid_bounding_box_diagonal_tol`] is their Euclidean
/// norm, the [`solid_bounding_box_volume_tol`] their product, and the
/// [`solid_bounding_box_aspect_ratio_tol`] the ratio of the largest to the smallest.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid is
/// empty (no bounding box).
pub fn solid_bounding_box_extents_tol(solid: &Solid, tol: f64) -> Result<[f64; 3], CadError> {
    let (min, max) = solid_bounding_box_tol(solid, tol)?;
    Ok([max[0] - min[0], max[1] - min[1], max[2] - min[2]])
}

/// AABB extents of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_bounding_box_extents_tol`].
pub fn solid_bounding_box_extents(solid: &Solid) -> Result<[f64; 3], CadError> {
    solid_bounding_box_extents_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **axis-aligned bounding-box elongation** = intermediate extent / longest extent
/// (dimensionless, in `(0, 1]`) of a solid — the Zingg `b/a` shape ratio, evaluated at
/// tessellation tolerance `tol` from the sorted [`solid_bounding_box_extents_tol`]. It
/// classifies a shape from rod-like (low, one axis dominates) to equant (`1`, all three
/// extents comparable); together with the flatness `c/b` it locates a particle on the Zingg
/// diagram. Distinct from the [`solid_bounding_box_aspect_ratio_tol`] (longest / shortest).
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid is
/// degenerate (a zero longest extent).
pub fn solid_bounding_box_elongation_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let mut s = solid_bounding_box_extents_tol(solid, tol)?;
    s.sort_by(|a, b| a.total_cmp(b));
    let (mid, longest) = (s[1], s[2]);
    if longest <= 0.0 {
        return Err(CadError::Tessellation(
            "degenerate (zero longest extent); elongation undefined".to_string(),
        ));
    }
    Ok(mid / longest)
}

/// AABB elongation of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_bounding_box_elongation_tol`].
pub fn solid_bounding_box_elongation(solid: &Solid) -> Result<f64, CadError> {
    solid_bounding_box_elongation_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **axis-aligned bounding-box flatness** = shortest extent / intermediate extent
/// (dimensionless, in `(0, 1]`) of a solid — the Zingg `c/b` shape ratio, evaluated at
/// tessellation tolerance `tol` from the sorted [`solid_bounding_box_extents_tol`]. It is
/// the second Zingg axis, the companion to the [`solid_bounding_box_elongation_tol`] `b/a`:
/// low for a flat / plate-like solid (the shortest axis collapses), `1` for an equant one
/// (all three extents comparable). Their product `flatness · elongation = c/a` is the
/// reciprocal of the [`solid_bounding_box_aspect_ratio_tol`] (longest / shortest).
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid is
/// degenerate (a zero intermediate extent).
pub fn solid_bounding_box_flatness_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let mut s = solid_bounding_box_extents_tol(solid, tol)?;
    s.sort_by(|a, b| a.total_cmp(b));
    let (shortest, mid) = (s[0], s[1]);
    if mid <= 0.0 {
        return Err(CadError::Tessellation(
            "degenerate (zero intermediate extent); flatness undefined".to_string(),
        ));
    }
    Ok(shortest / mid)
}

/// AABB flatness of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_bounding_box_flatness_tol`].
pub fn solid_bounding_box_flatness(solid: &Solid) -> Result<f64, CadError> {
    solid_bounding_box_flatness_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **Corey shape factor** `CSF = c / √(a·b)` (dimensionless, in `(0, 1]`) of a
/// solid, where `a ≥ b ≥ c` are the sorted AABB extents — the standard
/// sedimentology / hydraulics particle-shape descriptor (it governs the drag and
/// settling velocity of a non-spherical grain), evaluated at tessellation tolerance
/// `tol` from the sorted [`solid_bounding_box_extents_tol`]. It is `1` for an equant
/// cube or sphere and tends to `0` for a flat plate or a thin rod. It compounds the two
/// Zingg axes: `CSF =` [`solid_bounding_box_flatness_tol`] `· √(`[`solid_bounding_box_elongation_tol`]`)`.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid is
/// degenerate (a zero intermediate or longest extent).
pub fn solid_bounding_box_corey_shape_factor_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let mut s = solid_bounding_box_extents_tol(solid, tol)?;
    s.sort_by(|a, b| a.total_cmp(b));
    let (shortest, mid, longest) = (s[0], s[1], s[2]);
    if mid <= 0.0 || longest <= 0.0 {
        return Err(CadError::Tessellation(
            "degenerate (zero extent); Corey shape factor undefined".to_string(),
        ));
    }
    Ok(shortest / (mid * longest).sqrt())
}

/// AABB Corey shape factor of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_bounding_box_corey_shape_factor_tol`].
pub fn solid_bounding_box_corey_shape_factor(solid: &Solid) -> Result<f64, CadError> {
    solid_bounding_box_corey_shape_factor_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **Sneed–Folk maximum-projection sphericity** `ψ_p = ∛(c² / (a·b))`
/// (dimensionless, in `(0, 1]`) of a solid, where `a ≥ b ≥ c` are the sorted AABB
/// extents — a classic sedimentology particle-shape descriptor (the ratio of the maximum
/// projection area of the particle to that of a sphere of equal volume), evaluated at
/// tessellation tolerance `tol`. It is `1` for an equant cube or sphere and drops toward
/// `0` for a blade, rod, or disc. It is the cube root of the
/// [`solid_bounding_box_corey_shape_factor`] squared: `ψ_p = CSF^(2/3)`.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid is
/// degenerate (a zero intermediate or longest extent).
pub fn solid_bounding_box_max_projection_sphericity_tol(
    solid: &Solid,
    tol: f64,
) -> Result<f64, CadError> {
    let mut s = solid_bounding_box_extents_tol(solid, tol)?;
    s.sort_by(|a, b| a.total_cmp(b));
    let (shortest, mid, longest) = (s[0], s[1], s[2]);
    if mid <= 0.0 || longest <= 0.0 {
        return Err(CadError::Tessellation(
            "degenerate (zero extent); maximum projection sphericity undefined".to_string(),
        ));
    }
    Ok((shortest * shortest / (mid * longest)).cbrt())
}

/// AABB maximum-projection sphericity of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_bounding_box_max_projection_sphericity_tol`].
pub fn solid_bounding_box_max_projection_sphericity(solid: &Solid) -> Result<f64, CadError> {
    solid_bounding_box_max_projection_sphericity_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **axis-aligned bounding-box aspect ratio** = longest extent / shortest extent
/// (dimensionless, `≥ 1`) of a solid — the elongation (slenderness) of its AABB, computed
/// at tessellation tolerance `tol` from [`solid_bounding_box_tol`]. It is `1` for a cube
/// or a sphere and grows as a shape becomes more elongated or plate-like; a standard
/// descriptor for mesh-quality screening, particle / grain shape classification, and
/// camera view-fit.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid
/// is degenerate (a zero-extent bounding box, where the ratio is undefined).
pub fn solid_bounding_box_aspect_ratio_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let (min, max) = solid_bounding_box_tol(solid, tol)?;
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    let longest = dx.max(dy).max(dz);
    let shortest = dx.min(dy).min(dz);
    if shortest <= 0.0 {
        return Err(CadError::Tessellation(
            "degenerate (zero-extent) bounding box; aspect ratio undefined".to_string(),
        ));
    }
    Ok(longest / shortest)
}

/// AABB aspect ratio of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_bounding_box_aspect_ratio_tol`].
pub fn solid_bounding_box_aspect_ratio(solid: &Solid) -> Result<f64, CadError> {
    solid_bounding_box_aspect_ratio_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **axis-aligned bounding-box space diagonal** `d = √(Δx² + Δy² + Δz²)` (model
/// units) of a solid — the length of the main diagonal of its AABB, computed at
/// tessellation tolerance `tol` from [`solid_bounding_box_tol`]. It is the standard
/// single-number measure of a solid's overall size (level-of-detail selection, frustum
/// culling, a conservative bounding radius for broad-phase collision, and tolerance
/// scaling); it is `a√3` for a cube of side `a` and `2R√3` for a sphere of radius `R`,
/// and is always at least as large as the longest individual extent.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the
/// solid is empty (no bounding box).
pub fn solid_bounding_box_diagonal_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let (min, max) = solid_bounding_box_tol(solid, tol)?;
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    Ok((dx * dx + dy * dy + dz * dz).sqrt())
}

/// AABB space diagonal of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_bounding_box_diagonal_tol`].
pub fn solid_bounding_box_diagonal(solid: &Solid) -> Result<f64, CadError> {
    solid_bounding_box_diagonal_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **axis-aligned bounding-box volume** `V_bbox = Δx·Δy·Δz` (model units³) of a
/// solid — the product of its AABB extents, computed at tessellation tolerance `tol`
/// from [`solid_bounding_box_tol`]. It is the occupancy volume used by spatial indices
/// (octree / BVH node extents, broad-phase culling) and the denominator of the fill
/// ratio [`solid_bounding_box_fill_ratio_tol`] (`fill = V / V_bbox`). It is always at
/// least the solid's own volume, with equality for an axis-aligned box.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the
/// solid is empty (no bounding box).
pub fn solid_bounding_box_volume_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let (min, max) = solid_bounding_box_tol(solid, tol)?;
    Ok((max[0] - min[0]) * (max[1] - min[1]) * (max[2] - min[2]))
}

/// AABB volume of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_bounding_box_volume_tol`].
pub fn solid_bounding_box_volume(solid: &Solid) -> Result<f64, CadError> {
    solid_bounding_box_volume_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **axis-aligned bounding-box surface area** `2·(Δx·Δy + Δy·Δz + Δz·Δx)` (model
/// units²) of a solid — the total area of the six AABB faces, computed at tessellation
/// tolerance `tol` from [`solid_bounding_box_tol`]. It is the cost metric of the
/// **surface-area heuristic (SAH)** used to build BVH / k-d-tree / octree ray-tracing
/// and broad-phase collision acceleration structures, where a node's expected traversal
/// cost is proportional to the surface area of its bounding box. For an axis-aligned box
/// it equals the solid's own surface area.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the solid
/// is empty (no bounding box).
pub fn solid_bounding_box_surface_area_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let (min, max) = solid_bounding_box_tol(solid, tol)?;
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    Ok(2.0 * (dx * dy + dy * dz + dz * dx))
}

/// AABB surface area of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_bounding_box_surface_area_tol`].
pub fn solid_bounding_box_surface_area(solid: &Solid) -> Result<f64, CadError> {
    solid_bounding_box_surface_area_tol(solid, DEFAULT_MEASURE_TOLERANCE)
}

/// The **axis-aligned bounding-box fill ratio** `V / V_bbox` (dimensionless) of a
/// solid — the fraction of its AABB that the solid actually occupies, the 3-D analogue
/// of rectangularity / extent. It is `1` for an axis-aligned box (which fills its AABB
/// exactly), `π/6 ≈ 0.524` for a sphere, and `π/4 ≈ 0.785` for an axis-aligned
/// cylinder; a box maximises it. A standard packing / shape descriptor (collision
/// pre-tests, nesting, particle morphology). Computed at tessellation tolerance `tol`
/// from [`solid_volume_tol`] and [`solid_bounding_box_tol`].
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, or the
/// bounding box is degenerate (zero extent on some axis).
pub fn solid_bounding_box_fill_ratio_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let v = solid_volume_tol(solid, tol)?;
    let (min, max) = solid_bounding_box_tol(solid, tol)?;
    let bbox_vol = (max[0] - min[0]) * (max[1] - min[1]) * (max[2] - min[2]);
    if bbox_vol <= 0.0 {
        return Err(CadError::Tessellation(
            "degenerate bounding box; fill ratio undefined".to_string(),
        ));
    }
    Ok(v / bbox_vol)
}

/// AABB fill ratio of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_bounding_box_fill_ratio_tol`].
pub fn solid_bounding_box_fill_ratio(solid: &Solid) -> Result<f64, CadError> {
    solid_bounding_box_fill_ratio_tol(solid, DEFAULT_MEASURE_TOLERANCE)
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

/// The **inertia anisotropy** `I_max / I_min` (dimensionless, `≥ 1`) of a solid — the
/// ratio of its largest to smallest centroidal principal moment of inertia
/// ([`solid_principal_moments_tol`]), evaluated at tessellation tolerance `tol`. Unlike
/// the axis-aligned [`solid_bounding_box_aspect_ratio_tol`] this is **rotation-invariant**
/// (it depends only on the mass distribution, not the orientation): it is `1` for any
/// inertially isotropic body — a sphere *and* a cube — and grows large for a slender rod
/// or thin plate, so it is a robust elongation / shape descriptor for grain and particle
/// classification.
///
/// # Errors
///
/// [`CadError::Tessellation`] if `tol` is not finite and strictly positive, the solid is
/// degenerate, or a principal moment is non-positive.
pub fn solid_inertia_anisotropy_tol(solid: &Solid, tol: f64) -> Result<f64, CadError> {
    let m = solid_principal_moments_tol(solid, tol)?;
    let max = m[0].max(m[1]).max(m[2]);
    let min = m[0].min(m[1]).min(m[2]);
    if min <= 0.0 {
        return Err(CadError::Tessellation(
            "degenerate inertia (non-positive principal moment)".to_string(),
        ));
    }
    Ok(max / min)
}

/// Inertia anisotropy of a solid at [`DEFAULT_MEASURE_TOLERANCE`]. See
/// [`solid_inertia_anisotropy_tol`].
pub fn solid_inertia_anisotropy(solid: &Solid) -> Result<f64, CadError> {
    solid_inertia_anisotropy_tol(solid, DEFAULT_MEASURE_TOLERANCE)
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

/// The **topological genus** `g` of a closed B-rep solid — the number of
/// through-holes ("handles") in its boundary surface. From the Euler characteristic
/// `χ = 2 − 2g` of a closed orientable surface (see [`euler_characteristic`]),
/// `g = (2 − χ)/2`: a simply-connected solid (box, sphere, cylinder, cone) is genus
/// `0`, a torus genus `1`, a double torus genus `2`. `None` for a mesh-backed solid
/// (triangle counts are not B-rep topology), matching [`euler_characteristic`].
pub fn solid_genus(solid: &Solid) -> Option<i64> {
    euler_characteristic(solid).map(|chi| (2 - chi) / 2)
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
    fn solid_bounding_box_surface_area_is_the_aabb_face_area() {
        use crate::primitives::sphere;

        // box(2,4,6): SA = 2·(2·4 + 4·6 + 2·6) = 2·44 = 88, and a box's AABB surface IS
        // the box, so = solid_area.
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        assert!(
            (solid_bounding_box_surface_area(&b).unwrap() - 88.0).abs() < 1e-9,
            "box AABB SA = 88"
        );
        assert!(
            (solid_bounding_box_surface_area(&b).unwrap() - solid_area(&b).unwrap()).abs() < 1e-9,
            "box's AABB surface is the box"
        );

        // Threads solid_bounding_box: SA = 2·(ΔxΔy + ΔyΔz + ΔzΔx).
        for s in [box_solid(2.0, 4.0, 6.0).unwrap(), sphere(2.0).unwrap()] {
            let (mn, mx) = solid_bounding_box(&s).unwrap();
            let (dx, dy, dz) = (mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]);
            let expected = 2.0 * (dx * dy + dy * dz + dz * dx);
            assert!(
                (solid_bounding_box_surface_area(&s).unwrap() - expected).abs() / expected < 1e-9,
                "SA = 2(ΔxΔy + ΔyΔz + ΔzΔx)"
            );
        }

        // A cube of side a has AABB SA 6a²; a sphere's AABB is the cube of side 2R.
        assert!(
            (solid_bounding_box_surface_area(&box_solid(2.0, 2.0, 2.0).unwrap()).unwrap() - 24.0)
                .abs()
                < 1e-9,
            "cube side 2 → 24"
        );
        let sph = sphere(2.0).unwrap();
        assert!(
            (solid_bounding_box_surface_area(&sph).unwrap() - 96.0).abs() / 96.0 < 1e-2,
            "sphere R=2 → 96"
        );

        // The AABB surface exceeds the inscribed sphere's own area.
        assert!(
            solid_bounding_box_surface_area(&sph).unwrap() > solid_area(&sph).unwrap(),
            "AABB surface > sphere area"
        );
    }

    #[test]
    fn solid_bounding_box_volume_is_the_extent_product() {
        use crate::primitives::sphere;

        // box(2,4,6): V_bbox = 2·4·6 = 48, and a box fills its AABB so = solid_volume.
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        assert!((solid_bounding_box_volume(&b).unwrap() - 48.0).abs() < 1e-9, "box V_bbox = 48");
        assert!(
            (solid_bounding_box_volume(&b).unwrap() - solid_volume(&b).unwrap()).abs() < 1e-9,
            "box fills its AABB"
        );

        // Threads solid_bounding_box: V_bbox = product of extents.
        for s in [box_solid(2.0, 4.0, 6.0).unwrap(), sphere(2.0).unwrap()] {
            let (mn, mx) = solid_bounding_box(&s).unwrap();
            let expected = (mx[0] - mn[0]) * (mx[1] - mn[1]) * (mx[2] - mn[2]);
            assert!(
                (solid_bounding_box_volume(&s).unwrap() - expected).abs() / expected < 1e-9,
                "V_bbox = Π extent"
            );
        }

        // Double-threads fill_ratio + volume via the exact identity V_bbox = V / fill.
        for s in [box_solid(2.0, 4.0, 6.0).unwrap(), sphere(2.0).unwrap()] {
            let from_fill = solid_volume(&s).unwrap() / solid_bounding_box_fill_ratio(&s).unwrap();
            assert!(
                (solid_bounding_box_volume(&s).unwrap() - from_fill).abs() / from_fill < 1e-9,
                "V_bbox = V / fill_ratio"
            );
        }

        // A sphere's AABB is the cube of side 2R, and the AABB volume bounds the solid.
        let sph = sphere(2.0).unwrap();
        assert!(
            (solid_bounding_box_volume(&sph).unwrap() - 64.0).abs() / 64.0 < 1e-2,
            "sphere R=2 → 64"
        );
        assert!(
            solid_bounding_box_volume(&sph).unwrap() >= solid_volume(&sph).unwrap(),
            "AABB ⊇ solid"
        );
    }

    #[test]
    fn solid_bounding_box_elongation_is_the_zingg_mid_over_long() {
        use crate::primitives::sphere;

        // Worked: box(2,4,6) → sorted extents [2,4,6], elongation = 4/6 ≈ 0.6667.
        let bx = box_solid(2.0, 4.0, 6.0).unwrap();
        let e = solid_bounding_box_elongation(&bx).unwrap();
        assert!((e - 4.0 / 6.0).abs() <= 1e-9 * e, "elongation = mid/longest = 4/6");

        // Threads solid_bounding_box_extents (#363): mid/longest of the sorted extents.
        let mut s = solid_bounding_box_extents(&bx).unwrap();
        s.sort_by(|a, b| a.total_cmp(b));
        assert!((e - s[1] / s[2]).abs() <= 1e-9 * e, "= sorted[1]/sorted[2]");

        // Range 0 < elongation ≤ 1.
        assert!(e > 0.0 && e <= 1.0, "0 < elongation ≤ 1");

        // Equant: a cube and a sphere → ≈ 1 (all extents comparable).
        let cube = box_solid(3.0, 3.0, 3.0).unwrap();
        assert!((solid_bounding_box_elongation(&cube).unwrap() - 1.0).abs() <= 1e-9, "cube → 1");
        assert!(
            (solid_bounding_box_elongation(&sphere(2.0).unwrap()).unwrap() - 1.0).abs() <= 2e-2,
            "sphere → 1"
        );

        // The _tol wrapper agrees with the default (box exact at any tol).
        assert!(
            (solid_bounding_box_elongation(&bx).unwrap()
                - solid_bounding_box_elongation_tol(&bx, 0.1).unwrap())
            .abs()
                < 1e-9,
            "_tol wrapper agrees"
        );
    }

    #[test]
    fn solid_bounding_box_flatness_is_the_zingg_short_over_mid() {
        use crate::primitives::sphere;

        // Worked: box(2,4,6) → sorted extents [2,4,6], flatness = 2/4 = 0.5.
        let bx = box_solid(2.0, 4.0, 6.0).unwrap();
        let f = solid_bounding_box_flatness(&bx).unwrap();
        assert!((f - 2.0 / 4.0).abs() <= 1e-9 * f, "flatness = short/mid = 2/4");

        // Threads solid_bounding_box_extents (#363): shortest/intermediate of the sorted.
        let mut s = solid_bounding_box_extents(&bx).unwrap();
        s.sort_by(|a, b| a.total_cmp(b));
        assert!((f - s[0] / s[1]).abs() <= 1e-9 * f, "= sorted[0]/sorted[1]");

        // Range 0 < flatness ≤ 1.
        assert!(f > 0.0 && f <= 1.0, "0 < flatness ≤ 1");

        // Cross-check threading elongation (#375) + aspect_ratio: flatness·elongation
        // = (c/b)·(b/a) = c/a = 1/aspect_ratio.
        let elong = solid_bounding_box_elongation(&bx).unwrap();
        let aspect = solid_bounding_box_aspect_ratio(&bx).unwrap();
        assert!(
            (f * elong - 1.0 / aspect).abs() <= 1e-9,
            "flatness·elongation = 1/aspect_ratio"
        );

        // Equant: a cube and a sphere → ≈ 1.
        let cube = box_solid(3.0, 3.0, 3.0).unwrap();
        assert!((solid_bounding_box_flatness(&cube).unwrap() - 1.0).abs() <= 1e-9, "cube → 1");
        assert!(
            (solid_bounding_box_flatness(&sphere(2.0).unwrap()).unwrap() - 1.0).abs() <= 2e-2,
            "sphere → 1"
        );

        // The _tol wrapper agrees with the default (box exact at any tol).
        assert!(
            (solid_bounding_box_flatness(&bx).unwrap()
                - solid_bounding_box_flatness_tol(&bx, 0.1).unwrap())
            .abs()
                < 1e-9,
            "_tol wrapper agrees"
        );
    }

    #[test]
    fn solid_bounding_box_corey_shape_factor_is_c_over_sqrt_ab() {
        use crate::primitives::sphere;

        // Worked: box(2,4,6) → sorted [2,4,6], CSF = 2/√(4·6) = 2/√24 ≈ 0.40824829.
        let bx = box_solid(2.0, 4.0, 6.0).unwrap();
        let csf = solid_bounding_box_corey_shape_factor(&bx).unwrap();
        assert!((csf - 0.40824829046).abs() <= 1e-9, "CSF = 2/√24 ≈ 0.40825");

        // Threads solid_bounding_box_extents (#363): c/√(a·b) of the sorted extents.
        let mut s = solid_bounding_box_extents(&bx).unwrap();
        s.sort_by(|a, b| a.total_cmp(b));
        assert!((csf - s[0] / (s[1] * s[2]).sqrt()).abs() <= 1e-9 * csf, "= s0/√(s1·s2)");

        // Cross-check threading flatness (#381) + elongation (#375): CSF = flatness·√elongation.
        let f = solid_bounding_box_flatness(&bx).unwrap();
        let e = solid_bounding_box_elongation(&bx).unwrap();
        assert!((csf - f * e.sqrt()).abs() <= 1e-9 * csf, "CSF = flatness·√elongation");

        // Range 0 < CSF ≤ 1; equant cube and sphere → ≈ 1.
        assert!(csf > 0.0 && csf <= 1.0, "0 < CSF ≤ 1");
        let cube = box_solid(3.0, 3.0, 3.0).unwrap();
        assert!(
            (solid_bounding_box_corey_shape_factor(&cube).unwrap() - 1.0).abs() <= 1e-9,
            "cube → 1"
        );
        assert!(
            (solid_bounding_box_corey_shape_factor(&sphere(2.0).unwrap()).unwrap() - 1.0).abs()
                <= 3e-2,
            "sphere → 1"
        );

        // The _tol wrapper agrees with the default (box exact at any tol).
        assert!(
            (solid_bounding_box_corey_shape_factor(&bx).unwrap()
                - solid_bounding_box_corey_shape_factor_tol(&bx, 0.1).unwrap())
            .abs()
                < 1e-9,
            "_tol wrapper agrees"
        );
    }

    #[test]
    fn solid_bounding_box_max_projection_sphericity_is_cbrt_c2_over_ab() {
        use crate::primitives::sphere;

        // (a) WORKED: box(2,4,6) → sorted [2,4,6], ψ_p = ∛(2²/(4·6)) = ∛(1/6) ≈ 0.55032.
        let bx = box_solid(2.0, 4.0, 6.0).unwrap();
        let psi = solid_bounding_box_max_projection_sphericity(&bx).unwrap();
        assert!((psi - 0.5503212081).abs() <= 1e-9, "ψ_p = ∛(1/6) ≈ 0.55032");

        // (b) THREAD solid_bounding_box_extents (#363): ∛(s0²/(s1·s2)) of the sorted.
        let mut s = solid_bounding_box_extents(&bx).unwrap();
        s.sort_by(|a, b| a.total_cmp(b));
        assert!(
            (psi - (s[0] * s[0] / (s[1] * s[2])).cbrt()).abs() <= 1e-9 * psi,
            "= ∛(s0²/(s1·s2))"
        );

        // (c) CROSS-CHECK threading corey_shape_factor (#387) (non-tautological): ψ_p = CSF^(2/3).
        let csf = solid_bounding_box_corey_shape_factor(&bx).unwrap();
        assert!((psi - csf.powf(2.0 / 3.0)).abs() <= 1e-9 * psi, "ψ_p = CSF^(2/3)");

        // (d) Range 0 < ψ_p ≤ 1; equant cube and sphere → ≈ 1.
        assert!(psi > 0.0 && psi <= 1.0, "0 < ψ_p ≤ 1");
        let cube = box_solid(3.0, 3.0, 3.0).unwrap();
        assert!(
            (solid_bounding_box_max_projection_sphericity(&cube).unwrap() - 1.0).abs() <= 1e-9,
            "cube → 1"
        );
        assert!(
            (solid_bounding_box_max_projection_sphericity(&sphere(2.0).unwrap()).unwrap() - 1.0)
                .abs()
                <= 3e-2,
            "sphere → 1"
        );

        // (e) The _tol wrapper agrees with the default (box exact at any tol).
        assert!(
            (solid_bounding_box_max_projection_sphericity(&bx).unwrap()
                - solid_bounding_box_max_projection_sphericity_tol(&bx, 0.1).unwrap())
            .abs()
                < 1e-9,
            "_tol wrapper agrees"
        );
    }

    #[test]
    fn solid_bounding_box_aspect_ratio_is_the_extent_elongation() {
        use crate::primitives::sphere;

        // box(2,4,6): longest/shortest = 6/2 = 3.
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        assert!((solid_bounding_box_aspect_ratio(&b).unwrap() - 3.0).abs() < 1e-9, "box 2×4×6 → 3");

        // Threads solid_bounding_box: longest extent over shortest.
        let (mn, mx) = solid_bounding_box(&b).unwrap();
        let (dx, dy, dz) = (mx[0] - mn[0], mx[1] - mn[1], mx[2] - mn[2]);
        let expected = dx.max(dy).max(dz) / dx.min(dy).min(dz);
        assert!(
            (solid_bounding_box_aspect_ratio(&b).unwrap() - expected).abs() / expected < 1e-9,
            "aspect = longest/shortest extent"
        );

        // A cube and a sphere have aspect ratio 1 (isotropic AABB).
        assert!(
            (solid_bounding_box_aspect_ratio(&box_solid(2.0, 2.0, 2.0).unwrap()).unwrap() - 1.0)
                .abs()
                < 1e-9,
            "cube → 1"
        );
        assert!(
            (solid_bounding_box_aspect_ratio(&sphere(2.0).unwrap()).unwrap() - 1.0).abs() < 1e-2,
            "sphere → ≈ 1"
        );

        // The aspect ratio is always ≥ 1 (long over short).
        assert!(solid_bounding_box_aspect_ratio(&b).unwrap() >= 1.0, "ratio ≥ 1");
    }

    #[test]
    fn solid_bounding_box_center_is_the_aabb_midpoint() {
        use crate::primitives::sphere;

        // Threads solid_bounding_box: c = (min + max)/2.
        for s in [box_solid(2.0, 4.0, 6.0).unwrap(), sphere(2.0).unwrap()] {
            let (mn, mx) = solid_bounding_box(&s).unwrap();
            let c = solid_bounding_box_center(&s).unwrap();
            for k in 0..3 {
                assert!((c[k] - (mn[k] + mx[k]) * 0.5).abs() < 1e-9, "c = (min+max)/2");
            }
        }

        // Validates against solid_centroid (independent path): for a centro-symmetric
        // solid the AABB center equals the mass centroid.
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        let cb = solid_bounding_box_center(&b).unwrap();
        let gb = solid_centroid(&b).unwrap();
        for k in 0..3 {
            assert!((cb[k] - gb[k]).abs() < 1e-9, "box: bbox center = mass centroid");
        }
        let sph = sphere(2.0).unwrap();
        let cs = solid_bounding_box_center(&sph).unwrap();
        let gs = solid_centroid(&sph).unwrap();
        for k in 0..3 {
            assert!((cs[k] - gs[k]).abs() < 2e-2, "sphere: bbox center ≈ mass centroid");
        }
    }

    #[test]
    fn solid_bounding_box_extents_are_the_raw_dimensions() {
        use crate::primitives::sphere;

        // Worked: box(2,4,6) extents are {2,4,6} (sorted, axis-order-independent).
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        let mut sorted = solid_bounding_box_extents(&b).unwrap();
        sorted.sort_by(|a, c| a.total_cmp(c));
        assert!(
            (sorted[0] - 2.0).abs() < 1e-9
                && (sorted[1] - 4.0).abs() < 1e-9
                && (sorted[2] - 6.0).abs() < 1e-9,
            "box extents are {{2, 4, 6}}"
        );

        for s in [box_solid(2.0, 4.0, 6.0).unwrap(), sphere(2.0).unwrap()] {
            let e = solid_bounding_box_extents(&s).unwrap();

            // Threads solid_bounding_box: extents == [max − min].
            let (mn, mx) = solid_bounding_box(&s).unwrap();
            for k in 0..3 {
                assert!((e[k] - (mx[k] - mn[k])).abs() < 1e-9, "extent = max − min");
            }

            // Cross-check solid_bounding_box_diagonal: d = √Σextent².
            let diag = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
            assert!(
                (diag - solid_bounding_box_diagonal(&s).unwrap()).abs() <= 1e-9 * diag,
                "√Σe² = diagonal"
            );
            // Cross-check solid_bounding_box_volume: V = ∏extent.
            let vol = e[0] * e[1] * e[2];
            assert!(
                (vol - solid_bounding_box_volume(&s).unwrap()).abs() <= 1e-9 * vol,
                "∏e = bbox volume"
            );
            // Cross-check solid_bounding_box_aspect_ratio: max/min.
            let aspect = e[0].max(e[1]).max(e[2]) / e[0].min(e[1]).min(e[2]);
            assert!(
                (aspect - solid_bounding_box_aspect_ratio(&s).unwrap()).abs() <= 1e-9 * aspect,
                "max/min = aspect ratio"
            );
        }
    }

    #[test]
    fn solid_bounding_box_diagonal_matches_the_aabb_extents() {
        use crate::primitives::sphere;

        // box(2,4,6): AABB extents (2,4,6) → d = √(4+16+36) = √56 ≈ 7.4833148.
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        assert!(
            (solid_bounding_box_diagonal(&b).unwrap() - 56.0_f64.sqrt()).abs() < 1e-9,
            "box diagonal = √56"
        );

        // Threads solid_bounding_box: d = √(Σ extent²).
        for s in [box_solid(2.0, 4.0, 6.0).unwrap(), sphere(2.0).unwrap()] {
            let (mn, mx) = solid_bounding_box(&s).unwrap();
            let expected =
                ((mx[0] - mn[0]).powi(2) + (mx[1] - mn[1]).powi(2) + (mx[2] - mn[2]).powi(2)).sqrt();
            assert!(
                (solid_bounding_box_diagonal(&s).unwrap() - expected).abs() / expected < 1e-9,
                "d = √(Σ extent²)"
            );
        }

        // A cube of side a has diagonal a√3; a sphere's AABB is the cube of side 2R.
        assert!(
            (solid_bounding_box_diagonal(&box_solid(2.0, 2.0, 2.0).unwrap()).unwrap()
                - 2.0 * 3.0_f64.sqrt())
            .abs()
                < 1e-9,
            "cube side 2 → 2√3"
        );
        assert!(
            (solid_bounding_box_diagonal(&sphere(2.0).unwrap()).unwrap() - 4.0 * 3.0_f64.sqrt())
                .abs()
                < 1e-2,
            "sphere R=2 → 4√3"
        );

        // The space diagonal exceeds the largest single extent.
        assert!(solid_bounding_box_diagonal(&b).unwrap() > 6.0, "diagonal > z-extent");
    }

    #[test]
    fn solid_mean_chord_length_is_cauchys_four_v_over_a() {
        use crate::primitives::sphere;

        // Threads solid_specific_surface_area (A/V): ⟨ℓ⟩ = 4 / SSA (independent path).
        for s in [box_solid(2.0, 4.0, 6.0).unwrap(), sphere(2.0).unwrap()] {
            let mc = solid_mean_chord_length(&s).unwrap();
            let ssa = solid_specific_surface_area(&s).unwrap();
            assert!((mc - 4.0 / ssa).abs() <= 1e-9 * mc, "⟨ℓ⟩ = 4/SSA");
        }

        // Worked sphere (Cauchy): the mean chord of a sphere of radius r is exactly 4r/3.
        let sph = solid_mean_chord_length(&sphere(2.0).unwrap()).unwrap();
        assert!((sph - 4.0 * 2.0 / 3.0).abs() <= 2e-2 * (4.0 * 2.0 / 3.0), "sphere → 4r/3 = 8/3");

        // Worked box(2,4,6): V=48, A=88 → ⟨ℓ⟩ = 192/88 ≈ 2.1818 (flat-faced, exact at any
        // tolerance — also exercises the _tol wrapper directly).
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        assert!(
            (solid_mean_chord_length(&b).unwrap() - 192.0 / 88.0).abs() <= 1e-9 * (192.0 / 88.0),
            "box → 192/88"
        );
        assert!(
            (solid_mean_chord_length_tol(&b, 0.1).unwrap() - 192.0 / 88.0).abs()
                <= 1e-9 * (192.0 / 88.0),
            "_tol wrapper agrees"
        );
    }

    #[test]
    fn solid_sauter_mean_diameter_is_six_v_over_a() {
        use crate::primitives::sphere;

        // (a) WORKED: cube(3) → d₃₂ = 6·27/54 = 3 (= side); box(2,4,6): V=48, A=88 →
        // d₃₂ = 6·48/88 = 288/88 ≈ 3.2727 (flat-faced, exact).
        let cube = box_solid(3.0, 3.0, 3.0).unwrap();
        assert!(
            (solid_sauter_mean_diameter(&cube).unwrap() - 3.0).abs() <= 1e-9 * 3.0,
            "cube → side"
        );
        let bx = box_solid(2.0, 4.0, 6.0).unwrap();
        let d = solid_sauter_mean_diameter(&bx).unwrap();
        assert!((d - 6.0 * 48.0 / 88.0).abs() <= 1e-9 * d, "box → 6·48/88");

        // (b) THREAD solid_mean_chord_length (#357): d₃₂ = 1.5·⟨ℓ⟩ (6V/A = 1.5·4V/A).
        assert!(
            (d - 1.5 * solid_mean_chord_length(&bx).unwrap()).abs() <= 1e-9 * d,
            "d₃₂ = 1.5·mean_chord"
        );

        // (c) THREAD solid_specific_surface_area (A/V): d₃₂ = 6 / SSA (independent path).
        assert!(
            (d - 6.0 / solid_specific_surface_area(&bx).unwrap()).abs() <= 1e-9 * d,
            "d₃₂ = 6/SSA"
        );

        // (d) SPHERE → diameter: a sphere of radius 2 returns its true diameter 2r = 4.
        assert!(
            (solid_sauter_mean_diameter(&sphere(2.0).unwrap()).unwrap() - 4.0).abs() <= 3e-2 * 4.0,
            "sphere → 2r = 4"
        );

        // (e) The _tol wrapper agrees with the default (box exact at any tol).
        assert!(
            (d - solid_sauter_mean_diameter_tol(&bx, 0.1).unwrap()).abs() < 1e-9,
            "_tol wrapper agrees"
        );
    }

    #[test]
    fn solid_specific_surface_area_matches_area_over_volume() {
        use crate::primitives::sphere;

        // box(2,4,6): A = 2·(2·4+2·6+4·6) = 88, V = 48 → A/V = 88/48 = 11/6.
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        assert!(
            (solid_specific_surface_area(&b).unwrap() - 88.0 / 48.0).abs() < 1e-9,
            "box A/V = 88/48"
        );

        // Threads solid_area + solid_volume: A/V is exactly their ratio.
        for s in [box_solid(2.0, 4.0, 6.0).unwrap(), sphere(2.0).unwrap()] {
            let expected = solid_area(&s).unwrap() / solid_volume(&s).unwrap();
            assert!(
                (solid_specific_surface_area(&s).unwrap() - expected).abs() / expected < 1e-9,
                "A/V = area / volume"
            );
        }

        // A sphere of radius R has A/V = 3/R; a cube of side a has 6/a.
        assert!(
            (solid_specific_surface_area(&sphere(2.0).unwrap()).unwrap() - 3.0 / 2.0).abs() < 1e-2,
            "sphere → 3/R"
        );
        let cube = box_solid(2.0, 2.0, 2.0).unwrap();
        assert!(
            (solid_specific_surface_area(&cube).unwrap() - 3.0).abs() < 1e-9,
            "cube side 2 → 6/2 = 3"
        );

        // Dimensional: A/V scales as 1/length, so doubling the size halves it.
        let small = solid_specific_surface_area(&box_solid(1.0, 1.0, 1.0).unwrap()).unwrap();
        let big = solid_specific_surface_area(&box_solid(2.0, 2.0, 2.0).unwrap()).unwrap();
        assert!((big - 0.5 * small).abs() < 1e-9, "A/V ∝ 1/length");
    }

    #[test]
    fn solid_bounding_box_fill_ratio_is_one_for_a_box() {
        use crate::primitives::sphere;
        use std::f64::consts::PI;

        // A box fills its AABB exactly: V = 48, bbox = 2·4·6 = 48 → ratio 1.
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        assert!((solid_bounding_box_fill_ratio(&b).unwrap() - 1.0).abs() < 1e-9, "box fills its AABB");

        // Threads solid_volume + solid_bounding_box: ratio = V / (extent product).
        for s in [box_solid(2.0, 4.0, 6.0).unwrap(), sphere(2.0).unwrap()] {
            let (mn, mx) = solid_bounding_box(&s).unwrap();
            let bbox_vol = (mx[0] - mn[0]) * (mx[1] - mn[1]) * (mx[2] - mn[2]);
            let expected = solid_volume(&s).unwrap() / bbox_vol;
            assert!(
                (solid_bounding_box_fill_ratio(&s).unwrap() - expected).abs() / expected < 1e-9,
                "ratio = V / V_bbox"
            );
        }

        // A sphere fills π/6 of its cube AABB; a cylinder π/4 of its square-section AABB.
        assert!(
            (solid_bounding_box_fill_ratio(&sphere(2.0).unwrap()).unwrap() - PI / 6.0).abs() < 1e-2,
            "sphere → π/6"
        );
        assert!(
            (solid_bounding_box_fill_ratio(&cylinder(2.0, 5.0).unwrap()).unwrap() - PI / 4.0).abs()
                < 1e-2,
            "cylinder → π/4"
        );

        // The box maximises the fill ratio (and 0 < ratio ≤ 1).
        let r_sphere = solid_bounding_box_fill_ratio(&sphere(2.0).unwrap()).unwrap();
        assert!(r_sphere > 0.0 && r_sphere < 1.0, "sphere fill in (0,1)");
        assert!(
            solid_bounding_box_fill_ratio(&b).unwrap() > r_sphere,
            "box fills more of its AABB than a sphere"
        );
    }

    #[test]
    fn solid_equivalent_sphere_radius_is_the_volume_equivalent_size() {
        use crate::primitives::sphere;

        // box(2,4,6): V = 48 → r_eq = (3·48/4π)^(1/3) = (36/π)^(1/3).
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        let expected = ((3.0 * 48.0) / (4.0 * std::f64::consts::PI)).cbrt();
        assert!(
            (solid_equivalent_sphere_radius(&b).unwrap() - expected).abs() / expected < 1e-9,
            "r_eq = (3V/4π)^(1/3)"
        );

        // Threads solid_volume: r_eq recomputed from the measured volume (box + sphere).
        for s in [box_solid(2.0, 4.0, 6.0).unwrap(), sphere(2.0).unwrap()] {
            let from_v = (3.0 * solid_volume(&s).unwrap() / (4.0 * std::f64::consts::PI)).cbrt();
            assert!(
                (solid_equivalent_sphere_radius(&s).unwrap() - from_v).abs() / from_v < 1e-9,
                "r_eq from V"
            );
        }

        // Sphere round-trip: the equivalent-sphere radius of a sphere is its own radius.
        assert!(
            (solid_equivalent_sphere_radius(&sphere(2.0).unwrap()).unwrap() - 2.0).abs() < 2e-2,
            "sphere R=2 → r_eq ≈ 2 (tessellation)"
        );

        // Monotonic in size: a bigger box has a bigger equivalent radius.
        assert!(
            solid_equivalent_sphere_radius(&box_solid(4.0, 4.0, 4.0).unwrap()).unwrap()
                > solid_equivalent_sphere_radius(&box_solid(2.0, 2.0, 2.0).unwrap()).unwrap(),
            "bigger box → bigger r_eq"
        );
    }

    #[test]
    fn solid_equivalent_cube_edge_is_the_volume_cube_root() {
        use crate::primitives::sphere;

        // (a) WORKED: box(2,3,4) → V = 24 → s = 24^(1/3) ≈ 2.88450.
        let bx = box_solid(2.0, 3.0, 4.0).unwrap();
        let edge = solid_equivalent_cube_edge(&bx).unwrap();
        assert!((edge - 24.0_f64.cbrt()).abs() <= 1e-9 * edge, "s = V^(1/3)");

        // (b) THREAD solid_volume (non-tautological round-trip): s³ = V.
        assert!(
            (edge.powi(3) - solid_volume(&bx).unwrap()).abs() <= 1e-9 * solid_volume(&bx).unwrap(),
            "s³ = V"
        );

        // (c) CUBE recovers its own side: box(5,5,5) → V = 125 → s = 5.
        let cube = box_solid(5.0, 5.0, 5.0).unwrap();
        assert!(
            (solid_equivalent_cube_edge(&cube).unwrap() - 5.0).abs() <= 1e-9 * 5.0,
            "cube → s = 5"
        );

        // (d) SPHERE: V = (4/3)π·R³ → s = ((4/3)π)^(1/3)·R ≈ 3.224 for R = 2 (tessellation).
        let s_sphere = solid_equivalent_cube_edge(&sphere(2.0).unwrap()).unwrap();
        let expected = ((4.0 / 3.0) * std::f64::consts::PI).cbrt() * 2.0;
        assert!((s_sphere - expected).abs() / expected < 3e-2, "sphere → s ≈ 3.224");

        // (e) The _tol wrapper agrees with the default (box exact at any tol).
        assert!(
            (edge - solid_equivalent_cube_edge_tol(&bx, 0.1).unwrap()).abs() < 1e-9,
            "_tol wrapper agrees"
        );
    }

    #[test]
    fn solid_equivalent_sphere_diameter_is_twice_the_radius() {
        use crate::primitives::sphere;
        use std::f64::consts::PI;

        for s in [box_solid(2.0, 4.0, 6.0).unwrap(), sphere(2.0).unwrap()] {
            // Threads solid_equivalent_sphere_radius: d = 2·r_eq.
            let d = solid_equivalent_sphere_diameter(&s).unwrap();
            assert!(
                (d - 2.0 * solid_equivalent_sphere_radius(&s).unwrap()).abs() <= 1e-9 * d,
                "d = 2·r_eq"
            );
            // Threads solid_equivalent_sphere_area via the sphere-surface relation A = π·d²
            // (A_eq = 4πr² = π(2r)² = πd²) — area-via-radius vs π·diameter² (independent).
            assert!(
                (solid_equivalent_sphere_area(&s).unwrap() - PI * d.powi(2)).abs()
                    <= 1e-9 * solid_equivalent_sphere_area(&s).unwrap(),
                "A_eq = π·d²"
            );
        }

        // A sphere's equivalent sphere is itself: sphere(2) → d ≈ 4 (= 2·radius).
        assert!(
            (solid_equivalent_sphere_diameter(&sphere(2.0).unwrap()).unwrap() - 4.0).abs()
                <= 2e-2 * 4.0,
            "sphere → d = 4"
        );

        // The _tol wrapper agrees with the default (a flat-faced box is exact at any tol).
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        assert!(
            (solid_equivalent_sphere_diameter(&b).unwrap()
                - solid_equivalent_sphere_diameter_tol(&b, 0.1).unwrap())
            .abs()
                < 1e-9,
            "_tol wrapper agrees"
        );
    }

    #[test]
    fn solid_equivalent_sphere_area_is_the_isoperimetric_floor() {
        use crate::primitives::sphere;

        // Threads solid_equivalent_sphere_radius: A_eq = 4π·r_eq².
        for s in [box_solid(2.0, 4.0, 6.0).unwrap(), sphere(2.0).unwrap()] {
            let r = solid_equivalent_sphere_radius(&s).unwrap();
            let expected = 4.0 * std::f64::consts::PI * r * r;
            assert!(
                (solid_equivalent_sphere_area(&s).unwrap() - expected).abs() / expected < 1e-9,
                "A_eq = 4π·r_eq²"
            );
        }

        // Validates the sphericity definition: Ψ = A_eq / A_actual (a different code path
        // than sphericity's π^(1/3)(6V)^(2/3) numerator, so it cross-checks both).
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        assert!(
            (solid_sphericity(&b).unwrap()
                - solid_equivalent_sphere_area(&b).unwrap() / solid_area(&b).unwrap())
            .abs()
                < 1e-9,
            "Ψ = A_eq / A"
        );

        // Isoperimetric bound: the equivalent sphere has the least surface, so A_eq < A.
        assert!(
            solid_equivalent_sphere_area(&b).unwrap() < solid_area(&b).unwrap(),
            "A_eq ≤ A (sphere minimises surface)"
        );

        // A sphere is its own equivalent sphere → A_eq ≈ its actual area.
        let sph = sphere(2.0).unwrap();
        assert!(
            (solid_equivalent_sphere_area(&sph).unwrap() - solid_area(&sph).unwrap()).abs()
                / solid_area(&sph).unwrap()
                < 2e-2,
            "sphere → A_eq ≈ A"
        );
    }

    #[test]
    fn solid_sphericity_is_one_for_a_sphere() {
        use crate::primitives::sphere;
        use std::f64::consts::PI;

        // A unit cube: Ψ = π^(1/3)·6^(2/3)/6 ≈ 0.806 (computed independently).
        let cube = box_solid(1.0, 1.0, 1.0).unwrap();
        let psi_cube = solid_sphericity(&cube).unwrap();
        let expect_cube = PI.cbrt() * 6.0_f64.powf(2.0 / 3.0) / 6.0;
        assert!((psi_cube - expect_cube).abs() / expect_cube < 1e-9, "cube Ψ = {psi_cube}");

        // Threads solid_volume + solid_area: Ψ = π^(1/3)·(6V)^(2/3)/A.
        for s in [box_solid(1.0, 2.0, 4.0).unwrap(), sphere(2.0).unwrap()] {
            let psi = solid_sphericity(&s).unwrap();
            let from_va = PI.cbrt() * (6.0 * solid_volume(&s).unwrap()).powf(2.0 / 3.0)
                / solid_area(&s).unwrap();
            assert!((psi - from_va).abs() / from_va < 1e-9, "Ψ = π^⅓(6V)^⅔/A");
        }

        // The sphere maximises sphericity → Ψ ≈ 1.
        let psi_sphere = solid_sphericity(&sphere(2.0).unwrap()).unwrap();
        assert!((psi_sphere - 1.0).abs() < 1e-2, "sphere Ψ ≈ 1, got {psi_sphere}");

        // Dimensionless ⇒ scale-invariant: a 1-cube and a 3-cube share Ψ.
        let psi_big = solid_sphericity(&box_solid(3.0, 3.0, 3.0).unwrap()).unwrap();
        assert!((psi_cube - psi_big).abs() / psi_cube < 1e-9, "Ψ is scale-invariant");

        // A non-cube box is less sphere-like than the cube, and still below 1.
        let psi_slab = solid_sphericity(&box_solid(1.0, 2.0, 4.0).unwrap()).unwrap();
        assert!(psi_slab < psi_cube && psi_slab < 1.0, "slab less spherical: {psi_slab}");
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
    fn solid_genus_counts_handles() {
        use crate::primitives::{sphere, torus};
        // g = (2 − χ)/2 threads euler_characteristic for every solid; simply-connected
        // primitives (χ = 2) are genus 0.
        for s in [
            box_solid(1.0, 1.0, 1.0).unwrap(),
            cylinder(1.0, 2.0).unwrap(),
            sphere(1.0).unwrap(),
        ] {
            assert_eq!(solid_genus(&s), euler_characteristic(&s).map(|c| (2 - c) / 2));
            assert_eq!(solid_genus(&s), Some(0), "simply-connected → genus 0");
        }
        // A torus has one through-hole → genus 1 (χ = 0).
        let t = torus(2.0, 0.5).unwrap();
        assert_eq!(euler_characteristic(&t), Some(0), "torus χ = 0");
        assert_eq!(solid_genus(&t), Some(1), "torus → genus 1");
        // A mesh-backed solid has no B-rep topology → None (matches euler).
    }

    #[test]
    fn solid_inertia_anisotropy_is_rotation_invariant_shape() {
        use crate::primitives::sphere;

        // Inertially isotropic bodies → anisotropy 1: a cube has equal principal moments.
        assert!(
            (solid_inertia_anisotropy(&box_solid(2.0, 2.0, 2.0).unwrap()).unwrap() - 1.0).abs()
                < 1e-9,
            "cube → 1 (isotropic inertia)"
        );
        // ... and so does a sphere (tessellated, looser tol).
        assert!(
            (solid_inertia_anisotropy(&sphere(2.0).unwrap()).unwrap() - 1.0).abs() < 2e-2,
            "sphere → ≈ 1"
        );

        // box(2,4,6): principal moments [208, 160, 80] → I_max/I_min = 208/80 = 2.6.
        let b = box_solid(2.0, 4.0, 6.0).unwrap();
        assert!(
            (solid_inertia_anisotropy(&b).unwrap() - 2.6).abs() < 1e-9 * 2.6,
            "box 2×4×6 → 208/80 = 2.6"
        );

        // Threads solid_principal_moments: max/min of the moments.
        let m = solid_principal_moments(&b).unwrap();
        let expected = m[0].max(m[1]).max(m[2]) / m[0].min(m[1]).min(m[2]);
        assert!(
            (solid_inertia_anisotropy(&b).unwrap() - expected).abs() / expected < 1e-9,
            "= I_max/I_min"
        );

        // Always ≥ 1.
        assert!(solid_inertia_anisotropy(&b).unwrap() >= 1.0, "anisotropy ≥ 1");
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


