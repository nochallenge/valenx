//! Phase 89 — `BRepOffsetAPI_MakePipeShell` with auxiliary spine
//! support.
//!
//! ## What OCCT does
//!
//! Higher-level wrapper around [`BRepFill_PipeShell`](crate::pipe_shell())
//! exposed through the OffsetAPI namespace. The key extra feature is
//! `SetMode(auxiliary_spine, with_contact, with_correction)` which
//! constrains the swept frame to maintain contact (or near-contact)
//! with a secondary path. Used for:
//!
//! - **Ducting** — the duct must hug an enclosure wall along its
//!   length even when the centreline is curvy.
//! - **Tube heat exchangers** — tube wraps a finned core; the fin
//!   provides the auxiliary spine.
//!
//! ## v1 status — real mesh-domain sweep with auxiliary-spine roll
//!
//! Honest implementation (Phase 89.5), built on the same Bishop /
//! rotation-minimising frame as [`crate::sweep_api_pipe()`] (the
//! shared transport now lives in [`crate::sweep_support`]).
//!
//! - **No auxiliary spine.** Pure parallel-transport sweep —
//!   identical behaviour to `sweep_api_pipe`.
//! - **With an auxiliary spine.** At every primary-spine station the
//!   transported frame's roll is *corrected* so the profile's local
//!   `+u` axis points toward the matching point on the auxiliary
//!   spine. This is OCCT's `auxiliary_spine` mode: the cross-section
//!   no longer twists freely, it tracks the secondary path. The two
//!   spines are matched by **normalised arc length** so they need not
//!   have the same number of points.
//!
//! ### Honest scope
//!
//! `with_contact = false` (orient toward the aux spine) is fully
//! implemented. `with_contact = true` in OCCT additionally *scales /
//! shears* the profile so it physically touches the auxiliary spine;
//! that profile-deformation step is **not** done here — `with_contact`
//! is accepted and currently behaves the same as `false` (orientation
//! only). True contact deformation needs a per-station profile re-fit
//! and is a Tier-3 follow-up. The result is a mesh-backed [`Solid`];
//! like every mesh-domain sweep it carries no BRep topology, so apply
//! it last in a feature chain.

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::OcctSurfaceError;
use crate::sweep_support::{
    arc_length_param, cross, dot, norm, normalize, perp_basis, rotate_frame,
    sample_polyline_at, sub, vertex_tangents, Vec3,
};

/// Pipe-shell sweep with an optional auxiliary spine.
///
/// `profile` is a closed cross-section polygon — its points are taken
/// in their own local plane: the first three non-collinear points
/// define the plane, and the polygon is auto-closed. `spine` is the
/// primary sweep path (a 3D polyline). `auxiliary_spine`, when
/// present, is a secondary 3D polyline that the profile's `+u` axis is
/// rolled to face at every station. `with_contact` is accepted for
/// API parity (see the module docs — orientation is honoured, profile
/// contact-deformation is not).
///
/// # Errors
///
/// - [`OcctSurfaceError::BadInput`] for a profile with fewer than 3
///   points, a spine with fewer than 2 points, an auxiliary spine
///   with fewer than 2 points, non-finite coordinates, a degenerate
///   (zero-length) spine segment, or a profile whose points are all
///   collinear (no cross-section plane).
pub fn sweep_api_pipe_shell(
    profile: &[[f64; 3]],
    spine: &[[f64; 3]],
    auxiliary_spine: Option<&[[f64; 3]]>,
    with_contact: bool,
) -> Result<Solid, OcctSurfaceError> {
    let _ = with_contact; // orientation-only in v1 — see module docs.
    if profile.len() < 3 {
        return Err(OcctSurfaceError::bad_input(
            "profile",
            format!("need at least 3 profile points, got {}", profile.len()),
        ));
    }
    if spine.len() < 2 {
        return Err(OcctSurfaceError::bad_input(
            "spine",
            "spine needs at least two points",
        ));
    }
    for p in profile.iter().chain(spine.iter()) {
        if p.iter().any(|c| !c.is_finite()) {
            return Err(OcctSurfaceError::bad_input(
                "profile/spine",
                "input contains a non-finite coordinate",
            ));
        }
    }
    if let Some(aux) = auxiliary_spine {
        if aux.len() < 2 {
            return Err(OcctSurfaceError::bad_input(
                "auxiliary_spine",
                "auxiliary spine needs at least two points",
            ));
        }
        for p in aux {
            if p.iter().any(|c| !c.is_finite()) {
                return Err(OcctSurfaceError::bad_input(
                    "auxiliary_spine",
                    "auxiliary spine contains a non-finite coordinate",
                ));
            }
        }
    }

    // Express the profile in its own local 2D frame so it can be
    // re-placed at every spine station.
    let profile_local = project_profile_to_local(profile)?;

    // Per-vertex tangents + transported frames along the primary spine.
    let tangents = vertex_tangents(spine)?;
    let stations = build_stations(spine, &tangents, auxiliary_spine);

    let mesh = sweep_mesh(&profile_local, &stations);
    Ok(Solid::from_mesh(mesh))
}

/// One spine station: origin + an orthonormal `(u, v)` cross-section
/// frame (tangent = u × v).
struct Station {
    origin: Vec3,
    u: Vec3,
    v: Vec3,
}

/// Build the per-station frames. Without an auxiliary spine this is a
/// pure parallel transport; with one, each frame's roll is corrected
/// so `+u` points toward the matching auxiliary-spine point.
fn build_stations(
    spine: &[[f64; 3]],
    tangents: &[Vec3],
    auxiliary_spine: Option<&[[f64; 3]]>,
) -> Vec<Station> {
    let n = spine.len();
    let mut stations = Vec::with_capacity(n);

    // Seed frame at vertex 0.
    let (mut u, mut v) = perp_basis(tangents[0]);
    for i in 0..n {
        if i > 0 {
            // Parallel-transport the frame across this segment.
            let (ru, rv) = rotate_frame(u, v, tangents[i - 1], tangents[i]);
            u = ru;
            v = rv;
        }
        let origin: Vec3 = spine[i];
        let (fu, fv) = match auxiliary_spine {
            None => (u, v),
            Some(aux) => {
                // Match by normalised arc length so spines of
                // different resolutions still pair up.
                let s = arc_length_param(spine, i);
                let aim = sample_polyline_at(aux, s);
                roll_toward(origin, tangents[i], u, v, aim)
            }
        };
        stations.push(Station {
            origin,
            u: fu,
            v: fv,
        });
    }
    stations
}

/// Rotate the frame `(u, v)` about the station tangent so the `+u`
/// axis points (as closely as the tangent-plane allows) toward the
/// world point `aim`. Returns the corrected `(u, v)`.
///
/// This is the auxiliary-spine roll correction: the component of
/// `aim - origin` lying in the cross-section plane becomes the new
/// `+u`, and `+v` is rebuilt orthogonal to keep a right-handed frame.
fn roll_toward(origin: Vec3, tangent: Vec3, u: Vec3, v: Vec3, aim: Vec3) -> (Vec3, Vec3) {
    let to_aim = sub(aim, origin);
    // Project `to_aim` onto the cross-section plane (remove the
    // tangent component).
    let along = dot(to_aim, tangent);
    let in_plane = [
        to_aim[0] - tangent[0] * along,
        to_aim[1] - tangent[1] * along,
        to_aim[2] - tangent[2] * along,
    ];
    if norm(in_plane) < 1e-9 {
        // The auxiliary point is on the spine axis — no roll signal;
        // keep the transported frame.
        return (u, v);
    }
    let new_u = normalize(in_plane);
    // Right-handed: v = tangent × u.
    let new_v = cross(tangent, new_u);
    (new_u, new_v)
}

/// Build the swept triangle mesh: side walls between consecutive
/// cross-section rings plus triangulated end caps.
fn sweep_mesh(profile_local: &[(f64, f64)], stations: &[Station]) -> Mesh {
    let p = profile_local.len();
    let rings = stations.len();
    let mut nodes = Vec::with_capacity(p * rings);
    for st in stations {
        for (x, y) in profile_local {
            nodes.push(nalgebra::Vector3::new(
                st.origin[0] + st.u[0] * x + st.v[0] * y,
                st.origin[1] + st.u[1] * x + st.v[1] * y,
                st.origin[2] + st.u[2] * x + st.v[2] * y,
            ));
        }
    }
    let mut conn: Vec<u32> = Vec::new();
    for r in 0..rings - 1 {
        let a = (r * p) as u32;
        let b = ((r + 1) * p) as u32;
        for k in 0..p {
            let k1 = ((k + 1) % p) as u32;
            let k0 = k as u32;
            conn.extend_from_slice(&[a + k0, a + k1, b + k1]);
            conn.extend_from_slice(&[a + k0, b + k1, b + k0]);
        }
    }
    // End caps (fan-triangulate the first and last rings).
    let last = ((rings - 1) * p) as u32;
    for k in 1..p - 1 {
        conn.extend_from_slice(&[0, (k + 1) as u32, k as u32]);
        conn.extend_from_slice(&[last, last + k as u32, last + (k + 1) as u32]);
    }
    let mut mesh = Mesh::new("pipe-shell-sweep");
    mesh.nodes = nodes;
    mesh.element_blocks.push(ElementBlock {
        element_type: ElementType::Tri3,
        connectivity: conn,
    });
    mesh.recompute_stats();
    mesh
}

/// Express a 3D profile polygon in its own local 2D `(u, v)` frame.
///
/// The profile plane's normal is the first non-degenerate
/// `edge0 × edgeK` cross product; the `(u, v)` basis spans that plane.
/// Every profile point is projected onto that basis **relative to the
/// profile centroid**, so the local coordinates are centred on the
/// profile's own centre. `sweep_mesh` then places that centroid on
/// each spine station — a profile centred on the origin therefore
/// produces a cross-section centred on the spine, instead of riding
/// the spine by its first vertex (which offset the whole swept tube
/// by the vertex-0-to-centroid vector).
fn project_profile_to_local(
    profile: &[[f64; 3]],
) -> Result<Vec<(f64, f64)>, OcctSurfaceError> {
    let origin = profile[0];
    let e0 = sub(profile[1], origin);
    if norm(e0) < 1e-12 {
        return Err(OcctSurfaceError::bad_input(
            "profile",
            "profile's first edge is degenerate",
        ));
    }
    // Find a normal from the first non-collinear triple.
    let mut normal = [0.0; 3];
    for q in profile.iter().skip(2) {
        let ek = sub(*q, origin);
        let n = cross(e0, ek);
        if norm(n) > 1e-12 {
            normal = normalize(n);
            break;
        }
    }
    if norm(normal) < 1e-12 {
        return Err(OcctSurfaceError::bad_input(
            "profile",
            "profile points are collinear — no cross-section plane",
        ));
    }
    let u = normalize(e0);
    let v = cross(normal, u);
    // Centroid of the profile vertices — the spine pierces the
    // cross-section plane here, so this is the local-frame origin.
    let n_pts = profile.len() as f64;
    let centroid = profile.iter().fold([0.0; 3], |acc, p| {
        [acc[0] + p[0], acc[1] + p[1], acc[2] + p[2]]
    });
    let centroid = [centroid[0] / n_pts, centroid[1] / n_pts, centroid[2] / n_pts];
    let local: Vec<(f64, f64)> = profile
        .iter()
        .map(|p| {
            let d = sub(*p, centroid);
            (dot(d, u), dot(d, v))
        })
        .collect();
    Ok(local)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A unit square profile in the XY plane.
    fn square() -> Vec<[f64; 3]> {
        vec![
            [-0.5, -0.5, 0.0],
            [0.5, -0.5, 0.0],
            [0.5, 0.5, 0.0],
            [-0.5, 0.5, 0.0],
        ]
    }

    #[test]
    fn pipe_shell_validates_inputs() {
        let err = sweep_api_pipe_shell(&[], &[[0.0; 3], [1.0; 3]], None, false).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn pipe_shell_rejects_short_spine() {
        let err = sweep_api_pipe_shell(&square(), &[[0.0; 3]], None, false).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn pipe_shell_without_aux_sweeps_a_mesh() {
        let solid = sweep_api_pipe_shell(
            &square(),
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 2.0]],
            None,
            false,
        )
        .expect("pipe shell should sweep");
        match &solid {
            Solid::Mesh(m) => {
                // 4-gon profile, 3 stations → 2 ring gaps × 4 quads × 2
                // tris = 16 wall + 2 caps × 2 = 20 triangles.
                assert_eq!(m.total_elements(), 20);
                assert_eq!(m.nodes.len(), 12);
            }
            Solid::Brep(_) => panic!("expected a mesh-backed solid"),
        }
    }

    #[test]
    fn pipe_shell_rejects_collinear_profile() {
        // Three collinear points define no plane.
        let line = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let err = sweep_api_pipe_shell(
            &line,
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            None,
            false,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn pipe_shell_rejects_short_auxiliary_spine() {
        let err = sweep_api_pipe_shell(
            &square(),
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            Some(&[[1.0, 0.0, 0.0]]),
            false,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn auxiliary_spine_rolls_the_profile_to_face_it() {
        // A straight +Z spine. Without an aux spine the profile's roll
        // is arbitrary-but-fixed. With an aux spine running parallel
        // along +X offset, every station's +u axis must point toward
        // +X (the aux spine direction in the cross-section plane).
        let spine = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 2.0]];
        // Aux spine: parallel to the primary, offset +X by 5.
        let aux = vec![[5.0, 0.0, 0.0], [5.0, 0.0, 1.0], [5.0, 0.0, 2.0]];
        let solid =
            sweep_api_pipe_shell(&square(), &spine, Some(&aux), false).unwrap();
        let Solid::Mesh(m) = &solid else {
            panic!("expected mesh solid");
        };
        // Profile is a unit square centred on the station. With +u
        // toward +X, the ring at station 0 spans x ∈ [-0.5, 0.5].
        let first_ring = &m.nodes[0..4];
        let max_x = first_ring.iter().map(|n| n.x).fold(f64::MIN, f64::max);
        assert!(
            (max_x - 0.5).abs() < 1e-6,
            "profile +u should face +X, max_x = {max_x}"
        );
        // The whole first ring sits in the z = 0 cross-section plane.
        for n in first_ring {
            assert!(n.z.abs() < 1e-6, "first ring should lie at z=0");
        }
    }

    #[test]
    fn with_contact_flag_is_accepted() {
        // with_contact = true must not error (it is orientation-only
        // in v1 — see module docs).
        let solid = sweep_api_pipe_shell(
            &square(),
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            Some(&[[2.0, 0.0, 0.0], [2.0, 0.0, 1.0]]),
            true,
        )
        .expect("with_contact should be accepted");
        assert!(matches!(solid, Solid::Mesh(_)));
    }
}
