//! Closed-BRep sewing — stitch the per-region fitted surfaces of the
//! Phase-23-v2 pipeline into a watertight solid.
//!
//! # The gap this closes
//!
//! The Phase-23-v2 pipeline ([`crate::feature_detect`] +
//! [`crate::reconstruct`]) detects regions and fits a NURBS surface per
//! region with a tolerance report — but it stops there: the fitted
//! patches are an *open* set of disconnected faces, not a closed solid.
//! This module performs the **sew**: it identifies which fitted regions
//! are adjacent, joins them along their shared boundary, and produces a
//! single watertight [`Solid`].
//!
//! # Strategy — recognise the closed cases, fall back otherwise
//!
//! A general mesh→BRep sew is a Tier-3 kernel problem. This v1 takes
//! the honest route: it **recognises the two common closed shapes** and
//! reconstructs each as a real `Solid::Brep`, and for everything else
//! it produces a watertight **mesh-backed** shell.
//!
//! - **Fitted box** — six planar regions that pair into three
//!   mutually-perpendicular parallel-plane pairs. From the six fitted
//!   planes [`sew_regions`] recovers the box's three edge lengths and
//!   its orientation frame, builds a canonical `valenx_cad::box_solid`,
//!   and rigidly re-orients it onto the recovered frame — a true closed
//!   `Solid::Brep` with six faces, twelve shared edges, eight shared
//!   vertices.
//! - **Fitted cylinder** — one cylindrical region plus up to two planar
//!   cap regions. The cylinder fit gives the axis + radius; the caps
//!   give the height; [`sew_regions`] builds a `valenx_cad::cylinder`
//!   and orients it onto the recovered axis — a true closed
//!   `Solid::Brep` (side + 2 caps, shared circular edges).
//! - **Anything else** — the fitted regions' own triangles are welded
//!   (coincident-vertex merge, the mesh-domain edge share) into a
//!   single shell; the result is a mesh-backed [`Solid`]. The
//!   [`SewReport`] records whether that welded shell came out
//!   watertight (every edge shared by exactly two triangles) or stayed
//!   an open patch set.
//!
//! # Honest scope
//!
//! - The **box** and **cylinder** closed cases reconstruct a genuine
//!   `Solid::Brep`. Every other input sews in the mesh domain — a
//!   watertight mesh-backed `Solid` when the welded patches close, an
//!   open patch set otherwise. A general parametric trim-and-stitch of
//!   arbitrary fitted NURBS faces into a `Solid::Brep` shell stays the
//!   documented Tier-3 follow-up (it needs the parametric-BRep
//!   substrate, the same gate as the BRep fillet).
//! - The box recogniser requires the six planes to be genuinely
//!   axis-paired and perpendicular within tolerance; a sheared or
//!   incomplete box falls through to the mesh-domain weld.
//! - "Watertight" for the mesh fallback is the discrete 2-manifold
//!   test (every edge incident to exactly two triangles); it does not
//!   re-fit the welded geometry to the NURBS patches.

use nalgebra::Vector3;

use valenx_cad::Solid;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

use crate::feature_detect::{detect_cylinders, detect_planes, PlanarRegion};
use crate::reconstruct::ReconstructError;

/// What kind of closed solid the sew produced — surfaced so a caller
/// (and the tests) can tell a real BRep reconstruction from a
/// mesh-domain fallback.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SewOutcome {
    /// The fitted regions were recognised as a **box** and rebuilt as a
    /// true closed `Solid::Brep` (six planar faces).
    ClosedBrepBox,
    /// The fitted regions were recognised as a **cylinder** and rebuilt
    /// as a true closed `Solid::Brep` (lateral face + caps).
    ClosedBrepCylinder,
    /// The fitted regions were welded in the mesh domain and the result
    /// is **watertight** — a closed mesh-backed `Solid` (every edge is
    /// shared by exactly two triangles).
    WatertightMesh,
    /// The fitted regions were welded but the result is **not closed** —
    /// the sew returned the open patch set as a mesh-backed `Solid`.
    /// The caller has the patches but no closed solid.
    OpenPatchSet,
}

impl SewOutcome {
    /// True when the sew produced a genuinely closed solid (a BRep box,
    /// a BRep cylinder, or a watertight mesh) — i.e. anything but
    /// [`SewOutcome::OpenPatchSet`].
    pub fn is_closed(self) -> bool {
        !matches!(self, SewOutcome::OpenPatchSet)
    }

    /// True when the sew produced a real parametric `Solid::Brep`
    /// (the box or cylinder case) rather than a mesh-backed solid.
    pub fn is_brep(self) -> bool {
        matches!(
            self,
            SewOutcome::ClosedBrepBox | SewOutcome::ClosedBrepCylinder
        )
    }
}

/// The result of a sew — the produced solid plus a report of what
/// happened.
#[derive(Clone, Debug)]
pub struct SewResult {
    /// The sewn solid. A `Solid::Brep` for the recognised box /
    /// cylinder cases, a mesh-backed `Solid` otherwise.
    pub solid: Solid,
    /// What kind of result this is — see [`SewOutcome`].
    pub outcome: SewOutcome,
    /// Number of fitted regions the sew consumed.
    pub region_count: usize,
    /// Number of **free** (boundary) edges in the welded mesh — edges
    /// incident to exactly one triangle. Zero ⟺ the welded mesh is
    /// watertight. For the BRep box / cylinder cases this is `0` by
    /// construction (a primitive BRep is closed).
    pub free_edge_count: usize,
}

/// A diagnostic-only [`SewReport`] alias — the report is the
/// [`SewResult`] minus the solid, but callers usually want both, so
/// [`SewResult`] is the primary type.
pub type SewReport = SewResult;

/// Sew the regions detected in `mesh` into a watertight solid.
///
/// Runs region detection ([`detect_planes`] + [`detect_cylinders`]),
/// then:
///
/// 1. If the regions look like a **box** (six planar regions forming
///    three perpendicular parallel-plane pairs), reconstruct a true
///    closed `Solid::Brep` box.
/// 2. Else if they look like a **cylinder** (one cylindrical region,
///    optionally with planar caps), reconstruct a true closed
///    `Solid::Brep` cylinder.
/// 3. Else weld the regions' triangles into a single shell and return
///    it as a mesh-backed `Solid`, reporting whether the weld came out
///    watertight.
///
/// `normal_tolerance_deg` / `distance_tolerance` are the plane-region
/// detection tolerances (see [`detect_planes`]).
///
/// # Errors
///
/// [`ReconstructError::EmptyMesh`] if `mesh` has no `Tri3` triangles.
pub fn sew_regions(
    mesh: &Mesh,
    normal_tolerance_deg: f64,
    distance_tolerance: f64,
) -> Result<SewResult, ReconstructError> {
    if mesh
        .element_blocks
        .iter()
        .all(|b| b.element_type != ElementType::Tri3)
    {
        return Err(ReconstructError::EmptyMesh);
    }

    let planes = detect_planes(mesh, normal_tolerance_deg, distance_tolerance);

    // --- closed case 1: a fitted box ---
    if let Some(brep) = try_sew_box(&planes) {
        return Ok(SewResult {
            solid: brep,
            outcome: SewOutcome::ClosedBrepBox,
            region_count: planes.len(),
            free_edge_count: 0,
        });
    }

    // --- closed case 2: a fitted cylinder ---
    let cylinders = detect_cylinders(mesh, 0.05);
    if let Some(brep) = try_sew_cylinder(mesh, &cylinders, &planes) {
        return Ok(SewResult {
            solid: brep,
            outcome: SewOutcome::ClosedBrepCylinder,
            region_count: cylinders.len() + planes.len(),
            free_edge_count: 0,
        });
    }

    // --- fallback: weld the regions' triangles into a shell ---
    let welded = weld_regions_mesh(mesh, &planes, distance_tolerance.max(1e-6));
    let free = count_free_edges(&welded);
    let outcome = if free == 0 && welded.total_elements() > 0 {
        SewOutcome::WatertightMesh
    } else {
        SewOutcome::OpenPatchSet
    };
    Ok(SewResult {
        solid: Solid::from_mesh(welded),
        outcome,
        region_count: planes.len(),
        free_edge_count: free,
    })
}

/// Try to recognise six planar regions as a box and rebuild it as a
/// closed `Solid::Brep`.
///
/// A box has exactly six faces that pair into three parallel-plane
/// pairs, and the three pair normals are mutually perpendicular.
/// Returns `None` if the regions do not satisfy that — the caller
/// then falls through to the mesh-domain weld.
fn try_sew_box(planes: &[PlanarRegion]) -> Option<Solid> {
    if planes.len() != 6 {
        return None;
    }
    // Collect the six oriented plane normals.
    let normals: Vec<Vector3<f64>> = planes
        .iter()
        .map(|r| Vector3::from(r.normal).normalize())
        .collect();

    // Pair each plane with its anti-parallel partner. After pairing we
    // must have exactly three pairs.
    let mut paired = [false; 6];
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..6 {
        if paired[i] {
            continue;
        }
        let mut partner = None;
        for j in (i + 1)..6 {
            if paired[j] {
                continue;
            }
            // Anti-parallel normals: dot ≈ −1 (the two opposite faces
            // of a box face opposite ways).
            if normals[i].dot(&normals[j]) < -0.985 {
                partner = Some(j);
                break;
            }
        }
        let p = partner?;
        paired[i] = true;
        paired[p] = true;
        pairs.push((i, p));
    }
    if pairs.len() != 3 {
        return None;
    }

    // The three pair axes (one representative normal per pair) must be
    // mutually perpendicular.
    let axes: Vec<Vector3<f64>> = pairs.iter().map(|&(i, _)| normals[i]).collect();
    for a in 0..3 {
        for b in (a + 1)..3 {
            if axes[a].dot(&axes[b]).abs() > 0.05 {
                return None; // not perpendicular — a sheared shape
            }
        }
    }

    // Edge length along each axis = perpendicular distance between the
    // pair's two parallel planes, measured between their centroids.
    let mut lengths = [0.0f64; 3];
    let mut centre = Vector3::zeros();
    for (k, &(i, j)) in pairs.iter().enumerate() {
        let ci = Vector3::from(planes[i].centroid);
        let cj = Vector3::from(planes[j].centroid);
        // Distance along this axis between the two opposite faces.
        lengths[k] = (cj - ci).dot(&axes[k]).abs();
        centre += ci + cj;
        if lengths[k] < 1e-9 {
            return None; // a degenerate (zero-thickness) pair
        }
    }
    centre /= 6.0; // mean of the six face centroids = box centre

    // Build the canonical axis-aligned box, then re-orient it.
    //
    // `box_solid(dx, dy, dz)` puts a corner at the origin and the body
    // in the +octant. We want the box centred and oriented onto the
    // recovered (axes[0], axes[1], axes[2]) frame.
    let canonical = valenx_cad::box_solid(lengths[0], lengths[1], lengths[2]).ok()?;
    // Move the canonical box so its centre sits at the origin.
    let centred = canonical
        .translated(-lengths[0] / 2.0, -lengths[1] / 2.0, -lengths[2] / 2.0)
        .ok()?;
    // Orient: rotate the canonical X/Y/Z frame onto axes[0/1/2]. We
    // build the rotation as two successive minimal rotations (X→axis0,
    // then the residual Y→axis1) — robust and avoids assembling a raw
    // matrix through the public `valenx-cad` API.
    let oriented = orient_frame(&centred, &axes)?;
    // Translate the oriented, centred box to the recovered box centre.
    oriented.translated(centre.x, centre.y, centre.z).ok()
}

/// Re-orient a centred solid so its canonical `+X / +Y / +Z` axes map
/// onto `target[0] / target[1] / target[2]` (an orthonormal frame).
///
/// Done as two minimal rotations: first rotate `+X` onto `target[0]`;
/// then, in the rotated frame, rotate the (already-rotated) `+Y` onto
/// `target[1]`. The third axis follows for free since both frames are
/// orthonormal and right-handed-or-flipped consistently.
fn orient_frame(solid: &Solid, target: &[Vector3<f64>]) -> Option<Solid> {
    let x = Vector3::new(1.0, 0.0, 0.0);
    // --- rotation 1: +X → target[0] ---
    let (axis1, angle1) = minimal_rotation(x, target[0])?;
    let step1 = solid
        .rotated((0.0, 0.0, 0.0), (axis1.x, axis1.y, axis1.z), angle1)
        .ok()?;
    // Where +Y landed after rotation 1.
    let y = Vector3::new(0.0, 1.0, 0.0);
    let y_rotated = rotate_vector(y, axis1, angle1);
    // --- rotation 2: rotated +Y → target[1], about target[0] ---
    // Project both onto the plane ⊥ target[0] and find the angle.
    let n = target[0];
    let yr_perp = (y_rotated - n * y_rotated.dot(&n)).normalize();
    let t1_perp = (target[1] - n * target[1].dot(&n)).normalize();
    let cos = yr_perp.dot(&t1_perp).clamp(-1.0, 1.0);
    let sin = yr_perp.cross(&t1_perp).dot(&n);
    let angle2 = sin.atan2(cos);
    let step2 = step1
        .rotated((0.0, 0.0, 0.0), (n.x, n.y, n.z), angle2)
        .ok()?;
    Some(step2)
}

/// The minimal rotation `(axis, angle)` carrying unit vector `from`
/// onto unit vector `to`. Handles the parallel and antiparallel
/// degenerate cases.
fn minimal_rotation(from: Vector3<f64>, to: Vector3<f64>) -> Option<(Vector3<f64>, f64)> {
    let from = from.normalize();
    let to = to.normalize();
    let dot = from.dot(&to).clamp(-1.0, 1.0);
    if dot > 1.0 - 1e-12 {
        // Already aligned — zero rotation (any axis works).
        return Some((Vector3::new(1.0, 0.0, 0.0), 0.0));
    }
    if dot < -1.0 + 1e-12 {
        // Antiparallel — 180° about any axis ⊥ `from`.
        let perp = if from.x.abs() < 0.9 {
            from.cross(&Vector3::new(1.0, 0.0, 0.0))
        } else {
            from.cross(&Vector3::new(0.0, 1.0, 0.0))
        };
        return Some((perp.normalize(), std::f64::consts::PI));
    }
    let axis = from.cross(&to).normalize();
    Some((axis, dot.acos()))
}

/// Rotate `v` about unit `axis` by `angle` radians (Rodrigues).
fn rotate_vector(v: Vector3<f64>, axis: Vector3<f64>, angle: f64) -> Vector3<f64> {
    let (s, c) = angle.sin_cos();
    v * c + axis.cross(&v) * s + axis * axis.dot(&v) * (1.0 - c)
}

/// Try to recognise a cylindrical region (optionally with planar caps)
/// and rebuild it as a closed `Solid::Brep` cylinder.
///
/// Returns `None` if there is no dominant cylindrical region — the
/// caller then falls through to the mesh-domain weld.
fn try_sew_cylinder(
    mesh: &Mesh,
    cylinders: &[crate::feature_detect::CylindricalRegion],
    planes: &[PlanarRegion],
) -> Option<Solid> {
    let cyl = cylinders.first()?;
    let axis = Vector3::from(cyl.axis_direction).normalize();
    let radius = cyl.radius;
    if radius <= 1e-9 {
        return None;
    }

    // The cylinder height: the span of the lateral region's vertices
    // projected onto the axis. Caps (planar regions perpendicular to
    // the axis) confirm the closed form but the height is read from
    // the lateral geometry so it works cap-or-no-cap.
    let mut min_along = f64::INFINITY;
    let mut max_along = f64::NEG_INFINITY;
    let tri_conn = flatten_tri3(mesh);
    for &t in &cyl.triangle_indices {
        if let Some(tri) = tri_conn.get(t) {
            for &node in tri {
                if let Some(p) = mesh.nodes.get(node as usize) {
                    let along = p.dot(&axis);
                    min_along = min_along.min(along);
                    max_along = max_along.max(along);
                }
            }
        }
    }
    if !min_along.is_finite() {
        return None;
    }
    let height = max_along - min_along;
    if height <= 1e-9 {
        return None;
    }

    // The base-disk centre: the cylinder's axis origin slid along the
    // axis to the minimum-along cross-section.
    let axis_origin = Vector3::from(cyl.axis_origin);
    let along0 = axis_origin.dot(&axis);
    let base_centre = axis_origin + axis * (min_along - along0);

    // The caps are reported only as a confidence signal — a closed
    // cylinder usually shows two planar regions whose normals are
    // parallel to the axis. We do not *require* them (a cap-less tube
    // scan still reconstructs to a closed primitive cylinder), but if
    // planar regions exist they should be axis-aligned, not arbitrary.
    let cap_like = planes
        .iter()
        .filter(|r| Vector3::from(r.normal).normalize().dot(&axis).abs() > 0.9)
        .count();
    // If there are planar regions and none of them is a cap, the shape
    // is probably not a plain cylinder — bail to the mesh weld.
    if !planes.is_empty() && cap_like == 0 {
        return None;
    }

    // Build the canonical +Z cylinder and orient it onto the axis.
    let canonical = valenx_cad::cylinder(radius, height).ok()?;
    let (rot_axis, angle) = minimal_rotation(Vector3::new(0.0, 0.0, 1.0), axis)?;
    let oriented = canonical
        .rotated((0.0, 0.0, 0.0), (rot_axis.x, rot_axis.y, rot_axis.z), angle)
        .ok()?;
    oriented
        .translated(base_centre.x, base_centre.y, base_centre.z)
        .ok()
}

/// Weld every region's triangles into a single shell mesh.
///
/// Each planar region contributes its own source triangles; they are
/// concatenated and coincident vertices within `tolerance` are merged,
/// so triangles from adjacent regions that meet along a shared
/// boundary now reference the same vertices — the mesh-domain
/// equivalent of sharing an edge. This is the fallback sew for shapes
/// the box / cylinder recognisers do not match.
fn weld_regions_mesh(mesh: &Mesh, planes: &[PlanarRegion], tolerance: f64) -> Mesh {
    let tri_conn = flatten_tri3(mesh);
    let mut combined = Mesh::new("sewn_patches");
    let mut block = ElementBlock::new(ElementType::Tri3);
    for region in planes {
        for &t in &region.triangle_indices {
            if let Some(tri) = tri_conn.get(t) {
                let base = combined.nodes.len() as u32;
                let mut ok = true;
                let mut verts = [Vector3::zeros(); 3];
                for (k, &node) in tri.iter().enumerate() {
                    match mesh.nodes.get(node as usize) {
                        Some(p) => verts[k] = *p,
                        None => {
                            ok = false;
                            break;
                        }
                    }
                }
                if !ok {
                    continue;
                }
                for v in verts {
                    combined.nodes.push(v);
                }
                block
                    .connectivity
                    .extend_from_slice(&[base, base + 1, base + 2]);
            }
        }
    }
    if !block.connectivity.is_empty() {
        combined.element_blocks.push(block);
    }
    // Weld the seam: coincident boundary vertices collapse to one, so
    // adjacent patches now share their edge vertices.
    let mut welded = valenx_mesh::boolean::merge_coincident_nodes(&combined, tolerance);
    welded.recompute_stats();
    welded
}

/// Count the **free** edges of a triangle mesh — edges incident to
/// exactly one triangle. A watertight 2-manifold has zero free edges
/// (every edge is shared by two triangles).
fn count_free_edges(mesh: &Mesh) -> usize {
    use std::collections::HashMap;
    let mut edge_uses: HashMap<(u32, u32), u32> = HashMap::new();
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                // Canonical (low, high) edge key so the two winding
                // directions count as the same edge.
                let key = if a < b { (a, b) } else { (b, a) };
                *edge_uses.entry(key).or_insert(0) += 1;
            }
        }
    }
    edge_uses.values().filter(|&&c| c == 1).count()
}

/// Flatten a mesh's `Tri3` connectivity into a list of `[u32; 3]`
/// triangles, in declaration order — the indexing every `*Region`
/// type's `triangle_indices` uses.
fn flatten_tri3(mesh: &Mesh) -> Vec<[u32; 3]> {
    let mut out = Vec::new();
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            out.push([tri[0], tri[1], tri[2]]);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A unit cube mesh — 12 triangles, the canonical box test input.
    fn cube_mesh() -> Mesh {
        let mut m = Mesh::new("cube");
        let verts = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (1.0, 0.0, 1.0),
            (1.0, 1.0, 1.0),
            (0.0, 1.0, 1.0),
        ];
        for (x, y, z) in verts {
            m.nodes.push(Vector3::new(x, y, z));
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        // -Z, +Z, -Y, +Y, -X, +X — outward-wound.
        block.connectivity.extend_from_slice(&[0, 2, 1, 0, 3, 2]);
        block.connectivity.extend_from_slice(&[4, 5, 6, 4, 6, 7]);
        block.connectivity.extend_from_slice(&[0, 1, 5, 0, 5, 4]);
        block.connectivity.extend_from_slice(&[3, 7, 6, 3, 6, 2]);
        block.connectivity.extend_from_slice(&[0, 4, 7, 0, 7, 3]);
        block.connectivity.extend_from_slice(&[1, 2, 6, 1, 6, 5]);
        m.element_blocks.push(block);
        m
    }

    /// A non-cube box `dx × dy × dz` mesh.
    fn box_mesh(dx: f64, dy: f64, dz: f64) -> Mesh {
        let mut m = Mesh::new("box");
        let verts = [
            (0.0, 0.0, 0.0),
            (dx, 0.0, 0.0),
            (dx, dy, 0.0),
            (0.0, dy, 0.0),
            (0.0, 0.0, dz),
            (dx, 0.0, dz),
            (dx, dy, dz),
            (0.0, dy, dz),
        ];
        for (x, y, z) in verts {
            m.nodes.push(Vector3::new(x, y, z));
        }
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity.extend_from_slice(&[0, 2, 1, 0, 3, 2]);
        block.connectivity.extend_from_slice(&[4, 5, 6, 4, 6, 7]);
        block.connectivity.extend_from_slice(&[0, 1, 5, 0, 5, 4]);
        block.connectivity.extend_from_slice(&[3, 7, 6, 3, 6, 2]);
        block.connectivity.extend_from_slice(&[0, 4, 7, 0, 7, 3]);
        block.connectivity.extend_from_slice(&[1, 2, 6, 1, 6, 5]);
        m.element_blocks.push(block);
        m
    }

    /// A faceted closed cylinder mesh: `sides` facets, `radius`,
    /// `height`, with both end caps (fan-triangulated).
    fn closed_cylinder_mesh(radius: f64, height: f64, sides: usize) -> Mesh {
        let mut m = Mesh::new("cyl");
        // Two rings of `sides` vertices.
        for ring in 0..2 {
            let z = ring as f64 * height;
            for s in 0..sides {
                let a = std::f64::consts::TAU * s as f64 / sides as f64;
                m.nodes
                    .push(Vector3::new(radius * a.cos(), radius * a.sin(), z));
            }
        }
        // Two cap-centre vertices.
        let bottom_centre = m.nodes.len() as u32;
        m.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        let top_centre = m.nodes.len() as u32;
        m.nodes.push(Vector3::new(0.0, 0.0, height));
        let mut block = ElementBlock::new(ElementType::Tri3);
        for s in 0..sides {
            let next = (s + 1) % sides;
            let b0 = s as u32;
            let b1 = next as u32;
            let t0 = (sides + s) as u32;
            let t1 = (sides + next) as u32;
            // Lateral wall.
            block.connectivity.extend_from_slice(&[b0, b1, t1]);
            block.connectivity.extend_from_slice(&[b0, t1, t0]);
            // Bottom cap fan.
            block
                .connectivity
                .extend_from_slice(&[bottom_centre, b1, b0]);
            // Top cap fan.
            block.connectivity.extend_from_slice(&[top_centre, t0, t1]);
        }
        m.element_blocks.push(block);
        m
    }

    #[test]
    fn empty_mesh_errors() {
        let m = Mesh::new("empty");
        let err = sew_regions(&m, 1.0, 1e-3).unwrap_err();
        assert!(matches!(err, ReconstructError::EmptyMesh));
    }

    #[test]
    fn sews_a_unit_cube_into_a_closed_brep_box() {
        // The headline closed case: a cube's six planar regions sew
        // into a real closed Solid::Brep with six faces.
        let m = cube_mesh();
        let result = sew_regions(&m, 1.0, 1e-3).unwrap();
        assert_eq!(result.outcome, SewOutcome::ClosedBrepBox);
        assert!(result.outcome.is_closed());
        assert!(result.outcome.is_brep(), "a box must sew to a real BRep");
        assert_eq!(result.region_count, 6, "a cube is six planar regions");
        assert_eq!(result.free_edge_count, 0, "a BRep box is closed");
        // The result is a real BRep with the topology of a box.
        match &result.solid {
            Solid::Brep(_) => {
                assert_eq!(result.solid.faces(), 6, "box has 6 faces");
                assert_eq!(result.solid.edges(), 12, "box has 12 edges");
                assert_eq!(result.solid.vertices(), 8, "box has 8 vertices");
            }
            _ => panic!("box sew must yield a Solid::Brep"),
        }
    }

    #[test]
    fn sews_a_non_cube_box_recovering_its_dimensions() {
        // A 3 × 2 × 5 box: the sew must recover a closed BRep box whose
        // tessellated bounding box matches the input dimensions.
        let m = box_mesh(3.0, 2.0, 5.0);
        let result = sew_regions(&m, 1.0, 1e-3).unwrap();
        assert_eq!(result.outcome, SewOutcome::ClosedBrepBox);
        // Tessellate the sewn BRep and measure its extent.
        let mesh = valenx_cad::solid_to_mesh(&result.solid, 0.1).unwrap();
        let mut lo = Vector3::repeat(f64::INFINITY);
        let mut hi = Vector3::repeat(f64::NEG_INFINITY);
        for n in &mesh.nodes {
            lo = lo.inf(n);
            hi = hi.sup(n);
        }
        let ext = hi - lo;
        // The recovered box has the same three edge lengths (order may
        // differ since the axes are recovered from the planes).
        let mut got = [ext.x, ext.y, ext.z];
        got.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut want = [3.0, 2.0, 5.0];
        want.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (g, w) in got.iter().zip(want.iter()) {
            assert!((g - w).abs() < 0.2, "edge {g} should be ~{w}");
        }
    }

    #[test]
    fn sews_a_closed_cylinder_into_a_closed_brep_cylinder() {
        // A faceted closed cylinder sews into a real Solid::Brep
        // cylinder.
        let m = closed_cylinder_mesh(2.0, 5.0, 32);
        let result = sew_regions(&m, 1.0, 1e-3).unwrap();
        assert_eq!(
            result.outcome,
            SewOutcome::ClosedBrepCylinder,
            "a closed cylinder must sew to a BRep cylinder"
        );
        assert!(result.outcome.is_brep());
        assert!(result.outcome.is_closed());
        assert_eq!(result.free_edge_count, 0);
        match &result.solid {
            Solid::Brep(_) => {
                // A truck cylinder has a lateral face + 2 caps.
                assert!(result.solid.faces() >= 3, "cylinder has ≥3 faces");
            }
            _ => panic!("cylinder sew must yield a Solid::Brep"),
        }
    }

    #[test]
    fn sew_outcome_helpers_classify_correctly() {
        assert!(SewOutcome::ClosedBrepBox.is_brep());
        assert!(SewOutcome::ClosedBrepBox.is_closed());
        assert!(SewOutcome::ClosedBrepCylinder.is_brep());
        assert!(SewOutcome::WatertightMesh.is_closed());
        assert!(!SewOutcome::WatertightMesh.is_brep());
        assert!(!SewOutcome::OpenPatchSet.is_closed());
        assert!(!SewOutcome::OpenPatchSet.is_brep());
    }

    #[test]
    fn count_free_edges_zero_for_a_closed_cube_mesh() {
        // The cube mesh is a closed 2-manifold once its coincident
        // vertices are welded — every edge is shared by 2 triangles.
        let m = cube_mesh();
        let welded = valenx_mesh::boolean::merge_coincident_nodes(&m, 1e-6);
        assert_eq!(
            count_free_edges(&welded),
            0,
            "a welded closed cube has no free edges"
        );
    }

    #[test]
    fn count_free_edges_nonzero_for_an_open_patch() {
        // A single triangle is an open patch — all three of its edges
        // are free (incident to one triangle only).
        let mut m = Mesh::new("tri");
        m.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(1.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(0.0, 1.0, 0.0));
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity.extend_from_slice(&[0, 1, 2]);
        m.element_blocks.push(block);
        assert_eq!(count_free_edges(&m), 3, "a lone triangle has 3 free edges");
    }

    #[test]
    fn minimal_rotation_aligns_x_to_z() {
        // The minimal rotation carrying +X onto +Z, then applied,
        // should land +X on +Z.
        let from = Vector3::new(1.0, 0.0, 0.0);
        let to = Vector3::new(0.0, 0.0, 1.0);
        let (axis, angle) = minimal_rotation(from, to).unwrap();
        let rotated = rotate_vector(from, axis, angle);
        assert!((rotated - to).norm() < 1e-9, "got {rotated:?}");
    }

    #[test]
    fn minimal_rotation_handles_antiparallel() {
        // +X → −X is the antiparallel degenerate case: a 180° turn.
        let from = Vector3::new(1.0, 0.0, 0.0);
        let to = Vector3::new(-1.0, 0.0, 0.0);
        let (axis, angle) = minimal_rotation(from, to).unwrap();
        let rotated = rotate_vector(from, axis, angle);
        assert!((rotated - to).norm() < 1e-9, "antiparallel got {rotated:?}");
    }

    #[test]
    fn a_non_box_non_cylinder_mesh_falls_back_to_a_mesh_sew() {
        // A shape that is neither a box nor a cylinder — a single
        // open triangle pair forming an L — must fall through to the
        // mesh-domain weld and report honestly.
        let mut m = Mesh::new("L");
        // Two triangles sharing one edge, both in z=0 (one planar
        // region) — not a closed solid.
        m.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(1.0, 0.0, 0.0));
        m.nodes.push(Vector3::new(1.0, 1.0, 0.0));
        m.nodes.push(Vector3::new(0.0, 1.0, 0.0));
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity.extend_from_slice(&[0, 1, 2, 0, 2, 3]);
        m.element_blocks.push(block);
        let result = sew_regions(&m, 1.0, 1e-3).unwrap();
        // One planar region, open → not a closed solid.
        assert!(
            !result.outcome.is_brep(),
            "a flat quad is no BRep primitive"
        );
        assert_eq!(result.outcome, SewOutcome::OpenPatchSet);
        assert!(result.free_edge_count > 0, "an open quad has free edges");
        // It still returns a (mesh-backed) solid, never an error.
        assert!(matches!(result.solid, Solid::Mesh(_)));
    }
}
