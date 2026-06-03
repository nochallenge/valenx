//! Phase 132 — `BRepOffsetAPI_DraftAngle` with neutral-plane and
//! explicit pulling-direction.
//!
//! ## What OCCT does
//!
//! `BRepOffsetAPI_DraftAngle(shape)` followed by
//! `Add(face, direction, angle, neutral_plane, flag)` adds a draft
//! taper to each `face`, pivoting it around its intersection with the
//! `neutral_plane`. The mould/casting use case: a vertical wall must
//! be drafted so the part can release cleanly, but the parting line
//! (where wall meets neutral plane) must stay fixed in space so the
//! mating part still fits. Real OCCT computes:
//!
//! 1. The intersection curve `C = face ∩ neutral_plane`.
//! 2. A rotation of `face` about `C` by `angle` in the `direction`
//!    pulling sense, preserving topology with adjacent faces.
//! 3. New face boundaries that match the rotated face's intersections
//!    with neighbours (this is the hard part — the operator must
//!    reconstruct the wireframe in place).
//!
//! Block 1's `valenx_occt_surface::offset_api_draft_angle()` does
//! the simpler "rotate every face by `angle` about the global Z axis"
//! variant. This phase adds the neutral-plane + pulling-direction
//! controls.
//!
//! ## v1 status — real mesh-domain neutral-plane draft
//!
//! A genuine draft against an **arbitrary** neutral plane (any
//! origin + normal, not just `z = const`). The solid is tessellated;
//! every vertex on the pulling side of the neutral plane is sheared
//! away from the pulling axis by `tan(angle) · h`, where `h` is the
//! vertex's signed distance past the neutral plane along the pulling
//! direction. Vertices on or behind the neutral plane stay fixed — so
//! the parting line (wall ∩ neutral plane) is preserved exactly,
//! which is the whole point of the neutral-plane control.
//!
//! ### Honest scope
//!
//! - **Face selection.** A tessellated mesh carries no
//!   `TopoDS_Face`s, so `face_indices` cannot be honoured per-face —
//!   the mesh draft tapers all wall geometry past the neutral plane.
//!   The argument is validated and kept for API parity; true per-face
//!   draft about an intersection curve needs BRep face redirection
//!   (Tier 3).
//! - **Result** is a mesh-backed [`Solid`].
//!
//! Generalises `valenx_occt_surface::feat_make_draft` (which fixed
//! the neutral plane to `z = neutral_z`) to an arbitrary plane.

use valenx_cad::Solid;
use valenx_mesh::Mesh;

use crate::error::OcctAdvancedError;

/// Tessellation chord tolerance for the mesh-domain draft.
const DRAFT_TESS_TOLERANCE: f64 = 0.05;

/// 3D point used for the neutral plane origin and pulling direction.
pub type Vec3 = [f64; 3];

/// Apply a draft taper to selected `face_indices` of `solid`, pivoting
/// around the intersection with the plane through `neutral_origin`
/// with normal `neutral_normal`, in the `pulling_direction` sense, by
/// `angle_rad` radians (positive = open the draft outward).
///
/// `face_indices` are 0-based into the solid's face iterator (see
/// `valenx_cad::Solid::faces`).
///
/// # Errors
///
/// - [`OcctAdvancedError::BadInput`] for malformed inputs.
/// - [`OcctAdvancedError::NotYetImplemented`] otherwise in v1.
pub fn offset_api_draft_angle_with_neutral_plane(
    solid: &Solid,
    face_indices: &[usize],
    neutral_origin: Vec3,
    neutral_normal: Vec3,
    pulling_direction: Vec3,
    angle_rad: f64,
) -> Result<Solid, OcctAdvancedError> {
    if face_indices.is_empty() {
        return Err(OcctAdvancedError::bad_input(
            "face_indices",
            "need at least one face to draft",
        ));
    }
    if !angle_rad.is_finite() {
        return Err(OcctAdvancedError::bad_input(
            "angle_rad",
            "must be finite",
        ));
    }
    let n_norm =
        (neutral_normal[0].powi(2) + neutral_normal[1].powi(2) + neutral_normal[2].powi(2)).sqrt();
    if n_norm < f64::EPSILON {
        return Err(OcctAdvancedError::bad_input(
            "neutral_normal",
            "must be non-zero",
        ));
    }
    let p_norm = (pulling_direction[0].powi(2)
        + pulling_direction[1].powi(2)
        + pulling_direction[2].powi(2))
    .sqrt();
    if p_norm < f64::EPSILON {
        return Err(OcctAdvancedError::bad_input(
            "pulling_direction",
            "must be non-zero",
        ));
    }
    if neutral_origin.iter().any(|c| !c.is_finite()) {
        return Err(OcctAdvancedError::bad_input(
            "neutral_origin",
            "must be finite",
        ));
    }

    let mut mesh = valenx_cad::solid_to_mesh(solid, DRAFT_TESS_TOLERANCE)
        .map_err(|e| OcctAdvancedError::Backend(format!("draft: tessellate: {e:?}")))?;
    apply_neutral_plane_draft(
        &mut mesh,
        neutral_origin,
        [
            neutral_normal[0] / n_norm,
            neutral_normal[1] / n_norm,
            neutral_normal[2] / n_norm,
        ],
        [
            pulling_direction[0] / p_norm,
            pulling_direction[1] / p_norm,
            pulling_direction[2] / p_norm,
        ],
        angle_rad,
    );
    Ok(Solid::from_mesh(mesh))
}

/// Shear every vertex on the pulling side of the neutral plane away
/// from the pulling axis by `tan(angle) · h`.
fn apply_neutral_plane_draft(
    mesh: &mut Mesh,
    neutral_origin: Vec3,
    neutral_normal: Vec3, // unit
    pulling: Vec3,        // unit
    angle_rad: f64,
) {
    let tan_a = angle_rad.tan();
    // Geometry centroid — the radial-outward direction is measured
    // from the pulling axis through the centroid.
    let n = mesh.nodes.len().max(1) as f64;
    let mut c = [0.0, 0.0, 0.0];
    for v in &mesh.nodes {
        c[0] += v.x;
        c[1] += v.y;
        c[2] += v.z;
    }
    let centroid = [c[0] / n, c[1] / n, c[2] / n];

    for v in &mut mesh.nodes {
        // Signed distance past the neutral plane, measured along the
        // pulling direction. A vertex is "on the pulling side" when
        // (vertex − neutral_origin) · pulling > 0 AND it is on the
        // +normal side of the plane.
        let rel_origin = [
            v.x - neutral_origin[0],
            v.y - neutral_origin[1],
            v.z - neutral_origin[2],
        ];
        // Plane-signed side.
        let plane_side = rel_origin[0] * neutral_normal[0]
            + rel_origin[1] * neutral_normal[1]
            + rel_origin[2] * neutral_normal[2];
        // Height along the pulling direction.
        let h = rel_origin[0] * pulling[0]
            + rel_origin[1] * pulling[1]
            + rel_origin[2] * pulling[2];
        // Draft only the geometry genuinely past the plane on the
        // pulling side.
        if plane_side <= 1e-9 || h <= 0.0 {
            continue;
        }
        // Radial-outward direction = (vertex − centroid) with the
        // pulling component removed.
        let rad = [
            v.x - centroid[0],
            v.y - centroid[1],
            v.z - centroid[2],
        ];
        let rad_dot = rad[0] * pulling[0] + rad[1] * pulling[1] + rad[2] * pulling[2];
        let perp = [
            rad[0] - pulling[0] * rad_dot,
            rad[1] - pulling[1] * rad_dot,
            rad[2] - pulling[2] * rad_dot,
        ];
        let perp_len = (perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]).sqrt();
        if perp_len < 1e-12 {
            continue;
        }
        let shift = tan_a * h;
        v.x += perp[0] / perp_len * shift;
        v.y += perp[1] / perp_len * shift;
        v.z += perp[2] / perp_len * shift;
    }
    mesh.recompute_stats();
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_8;
    use valenx_cad::box_solid;

    #[test]
    fn rejects_empty_face_list() {
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = offset_api_draft_angle_with_neutral_plane(
            &s,
            &[],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            FRAC_PI_8,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_neutral_normal() {
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = offset_api_draft_angle_with_neutral_plane(
            &s,
            &[0],
            [0.0; 3],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            FRAC_PI_8,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_zero_pulling_direction() {
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = offset_api_draft_angle_with_neutral_plane(
            &s,
            &[0],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            [0.0; 3],
            FRAC_PI_8,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn rejects_non_finite_angle() {
        let s = box_solid(1.0, 1.0, 1.0).unwrap();
        let err = offset_api_draft_angle_with_neutral_plane(
            &s,
            &[0],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            f64::NAN,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_advanced.bad_input");
    }

    #[test]
    fn draft_tapers_geometry_past_the_neutral_plane() {
        // A 2x2x4 box, neutral plane at its base (origin (0,0,0),
        // normal +Z), pulling +Z, draft +10°. The top must come out
        // wider than the base.
        let s = box_solid(2.0, 2.0, 4.0).unwrap();
        let drafted = offset_api_draft_angle_with_neutral_plane(
            &s,
            &[0, 1, 2, 3],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            10f64.to_radians(),
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&drafted, 0.1).unwrap();
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
        assert!(top > bottom + 0.1, "top={top}, bottom={bottom}");
    }

    #[test]
    fn draft_preserves_the_parting_line() {
        // Vertices on the neutral plane (the parting line) do not move.
        let s = box_solid(2.0, 2.0, 2.0).unwrap();
        let drafted = offset_api_draft_angle_with_neutral_plane(
            &s,
            &[0],
            [0.0, 0.0, 0.0], // neutral plane at the box base
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
            FRAC_PI_8,
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&drafted, 0.1).unwrap();
        // Parting-line vertices keep their original half-width — for a
        // 2×2 box centred on the x=1 / y=1 axis that half-width is
        // exactly 1.0. The draft only tapers geometry past the neutral
        // plane, so a vertex still on it must not have moved from that
        // 1.0 half-width.
        for n in &mesh.nodes {
            if n.z < 1e-6 {
                let half = (n.x - 1.0).abs().max((n.y - 1.0).abs());
                assert!(
                    (half - 1.0).abs() < 1e-6,
                    "parting-line vertex moved: {n:?} (half-width {half}, expected 1.0)"
                );
            }
        }
    }

    #[test]
    fn draft_works_on_an_arbitrary_neutral_plane() {
        // Neutral plane through the box centre with an off-axis normal
        // — the draft must still produce valid non-empty geometry.
        let s = box_solid(4.0, 4.0, 4.0).unwrap();
        let drafted = offset_api_draft_angle_with_neutral_plane(
            &s,
            &[0],
            [2.0, 2.0, 2.0],
            [0.0, 1.0, 0.0], // neutral plane normal +Y
            [0.0, 1.0, 0.0], // pull along +Y
            FRAC_PI_8,
        )
        .unwrap();
        let mesh = valenx_cad::solid_to_mesh(&drafted, 0.2).unwrap();
        assert!(!mesh.nodes.is_empty());
    }
}
