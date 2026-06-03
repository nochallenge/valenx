//! Phase 100 — `BRepFeat_MakeDFill` (draft on a face relative to a
//! neutral plane).
//!
//! ## What OCCT does
//!
//! Feature-based counterpart to [`crate::offset_api_draft_angle()`]:
//! applies a draft to selected faces of an existing solid, with the
//! difference that the result is *committed* as a feature in the
//! topology — so the feature can subsequently be modified, suppressed,
//! or rolled back through the parametric history. Maps to SolidWorks
//! "Draft", Inventor "Face Draft", FreeCAD Part Design "Draft".
//!
//! ## v1 status — real mesh-domain draft taper
//!
//! This is a genuine draft, implemented in the mesh domain. The solid
//! is tessellated and every vertex on the **release side** of the
//! neutral plane (the side `direction` points toward) is sheared
//! *radially outward* by `tan(angle) · h`, where `h` is the vertex's
//! signed distance past the neutral plane along `direction`. That is
//! exactly the wall-taper a mould draft produces: walls flush with
//! the neutral plane stay put; geometry further from it leans out so
//! the part releases from the mould.
//!
//! ### Honest scope
//!
//! - **Face selection.** OCCT's `BRepFeat` drafts a *selected set of
//!   `TopoDS_Face`s*. A tessellated mesh has no `TopoDS_Face`s, so
//!   `face_indices` cannot be honoured per-face — the mesh-domain
//!   draft tapers *all* wall geometry past the neutral plane. The
//!   argument is validated (non-empty) and retained for API parity;
//!   true per-face draft needs BRep face redirection (Tier 3).
//! - **Result** is a mesh-backed [`Solid`]; downstream booleans
//!   refuse it (apply the draft last in a feature chain).

use valenx_cad::Solid;
use valenx_mesh::Mesh;

use crate::error::OcctSurfaceError;

/// Tessellation chord tolerance for the mesh-domain draft.
const DRAFT_TESS_TOLERANCE: f64 = 0.05;

/// Apply a draft to `base`, measured against `neutral_plane_z`.
///
/// `direction` is the mould-release direction; `angle_rad` is the
/// draft angle; `neutral_plane_z` is the world-Z of the neutral plane.
/// `face_indices` is validated for API parity — see the module docs on
/// the mesh-domain face-selection scope.
///
/// # Errors
///
/// - [`OcctSurfaceError::BadInput`] for an empty face list, a zero
///   `direction`, or a non-finite / zero angle.
/// - [`OcctSurfaceError::TruckLimit`] if `base` fails to tessellate.
pub fn feat_make_draft(
    base: &Solid,
    face_indices: &[usize],
    direction: [f64; 3],
    angle_rad: f64,
    neutral_plane_z: f64,
) -> Result<Solid, OcctSurfaceError> {
    if face_indices.is_empty() {
        return Err(OcctSurfaceError::bad_input(
            "face_indices",
            "need at least one face to draft",
        ));
    }
    let dir_len = (direction[0].powi(2) + direction[1].powi(2) + direction[2].powi(2)).sqrt();
    if dir_len < f64::EPSILON {
        return Err(OcctSurfaceError::bad_input(
            "direction",
            "must be non-zero",
        ));
    }
    if !angle_rad.is_finite() || angle_rad.abs() < f64::EPSILON {
        return Err(OcctSurfaceError::bad_input(
            "angle_rad",
            "must be non-zero finite",
        ));
    }
    if !neutral_plane_z.is_finite() {
        return Err(OcctSurfaceError::bad_input(
            "neutral_plane_z",
            "must be finite",
        ));
    }

    // Tessellate, then taper.
    let mut mesh = valenx_cad::solid_to_mesh(base, DRAFT_TESS_TOLERANCE)
        .map_err(|e| OcctSurfaceError::TruckLimit(format!("draft: tessellate: {e:?}")))?;
    apply_draft(&mut mesh, direction, dir_len, angle_rad, neutral_plane_z);
    Ok(Solid::from_mesh(mesh))
}

/// Shear every vertex on the release side of the neutral plane
/// radially outward by `tan(angle) · h`.
fn apply_draft(
    mesh: &mut Mesh,
    direction: [f64; 3],
    dir_len: f64,
    angle_rad: f64,
    neutral_z: f64,
) {
    // Unit release direction.
    let d = [
        direction[0] / dir_len,
        direction[1] / dir_len,
        direction[2] / dir_len,
    ];
    // Centroid of the geometry — the radial-outward direction at each
    // vertex is measured from the centroid axis.
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for n in &mesh.nodes {
        cx += n.x;
        cy += n.y;
        cz += n.z;
    }
    let count = mesh.nodes.len().max(1) as f64;
    let centroid = [cx / count, cy / count, cz / count];
    let tan_a = angle_rad.tan();

    for n in &mut mesh.nodes {
        // Signed height past the neutral plane along the release
        // direction. The neutral plane is z = neutral_z; the release
        // side is where `d` points. We measure h as the component of
        // (vertex − neutral-point) along d.
        let to_vertex = [n.x - centroid[0], n.y - centroid[1], n.z - neutral_z];
        let h = to_vertex[0] * d[0] + to_vertex[1] * d[1] + to_vertex[2] * d[2];
        if h <= 0.0 {
            continue; // on / below the neutral plane — no taper
        }
        // Radial-outward direction in the plane perpendicular to `d`:
        // (vertex − centroid) projected off `d`.
        let rad = [n.x - centroid[0], n.y - centroid[1], n.z - centroid[2]];
        let rad_dot = rad[0] * d[0] + rad[1] * d[1] + rad[2] * d[2];
        let perp = [
            rad[0] - d[0] * rad_dot,
            rad[1] - d[1] * rad_dot,
            rad[2] - d[2] * rad_dot,
        ];
        let perp_len = (perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]).sqrt();
        if perp_len < 1e-12 {
            continue; // on the release axis — nothing to taper
        }
        let outward = [
            perp[0] / perp_len,
            perp[1] / perp_len,
            perp[2] / perp_len,
        ];
        let shift = tan_a * h;
        n.x += outward[0] * shift;
        n.y += outward[1] * shift;
        n.z += outward[2] * shift;
    }
    mesh.recompute_stats();
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::box_solid;

    #[test]
    fn feat_draft_rejects_empty_faces() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_draft(&base, &[], [0.0, 0.0, 1.0], 0.05, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn feat_draft_rejects_zero_direction() {
        let base = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = feat_make_draft(&base, &[0], [0.0; 3], 0.05, 0.0).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn feat_draft_tapers_geometry_above_the_neutral_plane() {
        // A 2x2x4 box drafted +5° with the neutral plane at its base
        // (z=0) and release direction +Z. The top of the box (z=4)
        // must end up *wider* than the base — the walls slant outward.
        let base = box_solid(2.0, 2.0, 4.0).unwrap();
        let drafted = feat_make_draft(
            &base,
            &[0, 1, 2, 3], // ignored per-face — see module docs
            [0.0, 0.0, 1.0],
            5f64.to_radians(),
            0.0,
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&drafted, 0.1).unwrap();
        // Max |x| among bottom vertices vs top vertices.
        let bottom = mesh
            .nodes
            .iter()
            .filter(|n| n.z < 0.5)
            .map(|n| (n.x - 1.0).abs().max((n.y - 1.0).abs()))
            .fold(0.0_f64, f64::max);
        let top = mesh
            .nodes
            .iter()
            .filter(|n| n.z > 3.5)
            .map(|n| (n.x - 1.0).abs().max((n.y - 1.0).abs()))
            .fold(0.0_f64, f64::max);
        assert!(
            top > bottom + 0.05,
            "drafted box should be wider at the top: top={top}, bottom={bottom}"
        );
    }

    #[test]
    fn feat_draft_leaves_the_neutral_plane_unchanged() {
        // Vertices exactly on the neutral plane don't move.
        let base = box_solid(2.0, 2.0, 2.0).unwrap();
        let drafted = feat_make_draft(
            &base,
            &[0],
            [0.0, 0.0, 1.0],
            10f64.to_radians(),
            0.0, // neutral plane at the box's base
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&drafted, 0.1).unwrap();
        // Bottom-face vertices keep their original half-width — for a
        // 2×2 box centred on the x=1 / y=1 axis that half-width is
        // exactly 1.0. The draft only tapers geometry above z=0, so a
        // vertex still on the neutral plane must not have moved
        // outward (or inward) from that 1.0 half-width.
        for n in &mesh.nodes {
            if n.z < 1e-6 {
                let half = (n.x - 1.0).abs().max((n.y - 1.0).abs());
                assert!(
                    (half - 1.0).abs() < 1e-6,
                    "neutral-plane vertex moved: {n:?} (half-width {half}, expected 1.0)"
                );
            }
        }
    }
}
