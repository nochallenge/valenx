//! Phase 91 — `BRepOffsetAPI_MakeEvolved` (evolved surface).
//!
//! ## What OCCT does
//!
//! `BRepOffsetAPI_MakeEvolved(profile, spine, ...)` sweeps a 2D
//! profile along a 2D spine such that the profile is normal to the
//! spine at every point. The result is a 3D evolved surface (a
//! `TopoDS_Shell` or `TopoDS_Face`). The motivating use case is
//! gear-tooth involute generation and moulding-profile extrusion,
//! hence the name "evolved" — the geometry "evolves" the profile's
//! cross-section along the spine's curvature.
//!
//! Concretely: the spine is a planar curve in the XY plane. At each
//! point of the spine, a local frame is erected with `+x` along the
//! spine's in-plane normal and `+y` along world `+Z`; the profile
//! `(p, q)` is placed into that frame. As the spine bends, the
//! profile rolls with it, so the swept surface "evolves" the profile
//! around the path's curvature.
//!
//! ## v1 status — real mesh-domain evolved surface
//!
//! Honest implementation (Phase 91.5). Built on the shared
//! frame-transport primitives in [`crate::sweep_support`]:
//!
//! 1. The 2D spine `(x, y)` is lifted to the 3D XY plane.
//! 2. At every spine vertex a frame is built — `+x` is the unit
//!    in-plane normal of the spine (90° rotation of the planar
//!    tangent), `+y` is world `+Z`.
//! 3. The 2D profile `(p, q)` is mapped into that frame at every
//!    station: `point = spine_xy + p·normal + q·ẑ`.
//! 4. Consecutive profile rows are stitched into a quad strip
//!    (two triangles per quad).
//!
//! The result is a mesh-backed [`Solid`] carrying the evolved surface
//! as a triangle strip — an open surface, not a closed solid (an
//! evolved surface generally is not closed). It carries no BRep
//! topology, so apply it last in a feature chain.
//!
//! ### Honest scope
//!
//! The profile is treated as an **open** polyline (a moulding cross
//! section), so the surface is a strip, not a tube. The spine is also
//! treated as open. A closed profile or closed spine would simply
//! produce a degenerate seam; callers wanting a closed swept tube
//! should use [`sweep_api_pipe`](fn@crate::sweep_api_pipe) /
//! [`pipe_shell`](fn@crate::pipe_shell) instead. True G2 surface
//! fitting of the evolved patch (OCCT returns a B-spline surface) is
//! a Tier-3 follow-up — this v1 is the
//! faithful tessellated geometry.

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::OcctSurfaceError;
use crate::sweep_support::{normalize, Vec3};

/// Build an evolved surface by sweeping `profile_xy` along `spine_xy`,
/// keeping the profile in the plane normal to the spine.
///
/// Both inputs are 2D polylines. `profile_xy` is the cross-section
/// (treated as open — a moulding profile); `spine_xy` is the planar
/// path the profile is swept along. The result is a mesh-backed
/// [`Solid`] holding the evolved surface as a triangle strip.
///
/// The profile's first coordinate runs along the spine's in-plane
/// normal; its second coordinate runs along world `+Z`.
///
/// # Errors
///
/// - [`OcctSurfaceError::BadInput`] for a profile or spine with fewer
///   than 2 points, non-finite coordinates, or a spine with a
///   zero-length segment (no defined tangent).
pub fn sweep_api_evolved(
    profile_xy: &[(f64, f64)],
    spine_xy: &[(f64, f64)],
) -> Result<Solid, OcctSurfaceError> {
    if profile_xy.len() < 2 {
        return Err(OcctSurfaceError::bad_input(
            "profile_xy",
            "need at least 2 profile points",
        ));
    }
    if spine_xy.len() < 2 {
        return Err(OcctSurfaceError::bad_input(
            "spine_xy",
            "need at least 2 spine points",
        ));
    }
    for (x, y) in profile_xy.iter().chain(spine_xy.iter()) {
        if !x.is_finite() || !y.is_finite() {
            return Err(OcctSurfaceError::bad_input(
                "profile_xy/spine_xy",
                "input contains a non-finite coordinate",
            ));
        }
    }

    // Per-vertex in-plane normals of the planar spine.
    let normals = spine_inplane_normals(spine_xy)?;

    // Lay out one profile row per spine station.
    let rows = profile_xy.len();
    let cols = spine_xy.len();
    let mut nodes: Vec<nalgebra::Vector3<f64>> = Vec::with_capacity(rows * cols);
    let z_axis: Vec3 = [0.0, 0.0, 1.0];
    for (c, &(sx, sy)) in spine_xy.iter().enumerate() {
        let nrm = normals[c];
        for &(p, q) in profile_xy {
            // point = spine + p·normal + q·ẑ
            nodes.push(nalgebra::Vector3::new(
                sx + p * nrm[0] + q * z_axis[0],
                sy + p * nrm[1] + q * z_axis[1],
                p * nrm[2] + q * z_axis[2],
            ));
        }
    }

    // Stitch the grid into a triangle strip. Node index of
    // (profile row r, spine column c) is `c * rows + r`.
    let mut conn: Vec<u32> = Vec::new();
    for c in 0..cols - 1 {
        for r in 0..rows - 1 {
            let a = (c * rows + r) as u32;
            let b = ((c + 1) * rows + r) as u32;
            let a1 = a + 1;
            let b1 = b + 1;
            // Two triangles per quad.
            conn.extend_from_slice(&[a, b, b1]);
            conn.extend_from_slice(&[a, b1, a1]);
        }
    }

    let mut mesh = Mesh::new("evolved-surface");
    mesh.nodes = nodes;
    mesh.element_blocks.push(ElementBlock {
        element_type: ElementType::Tri3,
        connectivity: conn,
    });
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// Per-vertex unit in-plane normals of a planar (XY) spine.
///
/// The spine tangent is averaged at interior vertices; the normal is
/// the tangent rotated +90° in the XY plane (`(tx, ty) → (-ty, tx)`),
/// lifted to 3D with a zero Z component.
fn spine_inplane_normals(spine_xy: &[(f64, f64)]) -> Result<Vec<Vec3>, OcctSurfaceError> {
    let n = spine_xy.len();
    let mut normals = Vec::with_capacity(n);
    for i in 0..n {
        let incoming = if i > 0 {
            (
                spine_xy[i].0 - spine_xy[i - 1].0,
                spine_xy[i].1 - spine_xy[i - 1].1,
            )
        } else {
            (0.0, 0.0)
        };
        let outgoing = if i + 1 < n {
            (
                spine_xy[i + 1].0 - spine_xy[i].0,
                spine_xy[i + 1].1 - spine_xy[i].1,
            )
        } else {
            (0.0, 0.0)
        };
        let len2 = |a: (f64, f64)| (a.0 * a.0 + a.1 * a.1).sqrt();
        let li = len2(incoming);
        let lo = len2(outgoing);
        let tangent = match (li > 1e-12, lo > 1e-12) {
            (true, true) => (
                incoming.0 / li + outgoing.0 / lo,
                incoming.1 / li + outgoing.1 / lo,
            ),
            (true, false) => (incoming.0 / li, incoming.1 / li),
            (false, true) => (outgoing.0 / lo, outgoing.1 / lo),
            (false, false) => {
                return Err(OcctSurfaceError::bad_input(
                    "spine_xy",
                    "spine has a zero-length segment",
                ));
            }
        };
        // Lift the tangent to 3D and rotate +90° about Z to get the
        // in-plane normal.
        let t3: Vec3 = normalize([tangent.0, tangent.1, 0.0]);
        let nrm: Vec3 = [-t3[1], t3[0], 0.0];
        normals.push(nrm);
    }
    Ok(normals)
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::solid_to_mesh;

    #[test]
    fn evolved_rejects_short_profile() {
        let err = sweep_api_evolved(&[(0.0, 0.0)], &[(0.0, 0.0), (5.0, 0.0)])
            .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn evolved_rejects_short_spine() {
        let err = sweep_api_evolved(&[(0.0, 0.0), (1.0, 0.0)], &[(0.0, 0.0)])
            .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn evolved_rejects_non_finite() {
        let err = sweep_api_evolved(
            &[(0.0, 0.0), (f64::NAN, 1.0)],
            &[(0.0, 0.0), (5.0, 0.0)],
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn straight_spine_evolved_surface_is_a_flat_strip() {
        // A 2-point profile (a vertical line q ∈ [0, 1]) swept along a
        // straight spine running +X. The spine's in-plane normal is
        // +Y, so the profile's `p` coordinate maps to Y and `q` to Z.
        // The result is a flat rectangle in the XZ-offset plane.
        let profile = [(0.0, 0.0), (0.0, 1.0)]; // p=0, q∈[0,1]
        let spine = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)];
        let solid = sweep_api_evolved(&profile, &spine).unwrap();
        let Solid::Mesh(m) = &solid else {
            panic!("expected mesh solid");
        };
        // 2 profile rows × 3 spine columns = 6 nodes.
        assert_eq!(m.nodes.len(), 6);
        // 2 quads → 4 triangles.
        assert_eq!(m.total_elements(), 4);
        // q maps to Z: the surface spans z ∈ [0, 1].
        let zmin = m.nodes.iter().map(|n| n.z).fold(f64::INFINITY, f64::min);
        let zmax = m.nodes.iter().map(|n| n.z).fold(f64::NEG_INFINITY, f64::max);
        assert!((zmin - 0.0).abs() < 1e-9 && (zmax - 1.0).abs() < 1e-9);
    }

    #[test]
    fn profile_offset_follows_the_spine_normal() {
        // Profile with a non-zero `p` (offset along the spine normal).
        // For a +X spine the normal is +Y, so a profile point at p=2
        // must land at y = 2.
        let profile = [(2.0, 0.0), (2.0, 1.0)];
        let spine = [(0.0, 0.0), (1.0, 0.0)];
        let solid = sweep_api_evolved(&profile, &spine).unwrap();
        let mesh = solid_to_mesh(&solid, 0.1).unwrap();
        // Every node's Y must be ≈ 2 (the profile offset along +Y).
        for n in &mesh.nodes {
            assert!((n.y - 2.0).abs() < 1e-6, "node Y should be 2, got {}", n.y);
        }
    }

    #[test]
    fn curved_spine_evolves_the_profile_around_the_bend() {
        // An L-shaped spine: the in-plane normal differs on the two
        // legs, so the evolved profile rolls around the corner. The
        // surface must contain nodes whose normal-offset points in
        // distinctly different directions — i.e. it is genuinely 3D /
        // non-planar, not a flat strip.
        let profile = [(1.0, 0.0), (1.0, 2.0)];
        let spine = [(0.0, 0.0), (3.0, 0.0), (3.0, 3.0)];
        let solid = sweep_api_evolved(&profile, &spine).unwrap();
        let mesh = solid_to_mesh(&solid, 0.1).unwrap();
        // The first leg runs +X (normal +Y) and the last leg runs +Y
        // (normal -X). So some nodes are offset in +Y and others in
        // -X — the X spread and Y spread are both non-trivial.
        let xs: Vec<f64> = mesh.nodes.iter().map(|n| n.x).collect();
        let ys: Vec<f64> = mesh.nodes.iter().map(|n| n.y).collect();
        let xspread = xs.iter().cloned().fold(f64::MIN, f64::max)
            - xs.iter().cloned().fold(f64::MAX, f64::min);
        let yspread = ys.iter().cloned().fold(f64::MIN, f64::max)
            - ys.iter().cloned().fold(f64::MAX, f64::min);
        assert!(xspread > 1.0, "evolved surface should spread in X");
        assert!(yspread > 1.0, "evolved surface should spread in Y");
    }
}
