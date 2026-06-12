//! Phase 88 — `BRepOffsetAPI_MakePipe` (basic profile sweep along a
//! spine).
//!
//! ## What OCCT does
//!
//! `BRepOffsetAPI_MakePipe(spine, profile)` is the "simple pipe"
//! constructor — sweep a single profile along a spine wire with no
//! guide curves, no profile evolution, just rigid-body translation
//! of the cross-section along the path's tangent frame. The result
//! is a `TopoDS_Solid` if the profile is a closed face, otherwise a
//! `TopoDS_Shell`.
//!
//! Equivalent to a Frenet-frame sweep with constant cross-section,
//! the simplest sweep in the OCCT family.
//!
//! ## v1 status
//!
//! **Honest implementation** (Phase 88.5). The axis-aligned
//! straight-spine case still delegates to [`valenx_cad::prism`] (a
//! prism is the degenerate one-segment pipe). The **general curved
//! polyline spine** is now implemented directly: the profile is
//! placed at every spine vertex with a parallel-transport (Bishop)
//! frame so the cross-section does not twist, consecutive cross-
//! section rings are connected with quad side walls, and the two
//! ends are capped with triangle fans. The swept body is returned as
//! a mesh-backed [`Solid`] (`Solid::from_mesh`) — it is real
//! geometry, but carries no BRep topology, so downstream booleans
//! refuse it (apply the pipe last in a feature chain, same rule as
//! the fillet pipeline).

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::OcctSurfaceError;

/// Sweep `profile_xy` along a polyline `spine`.
///
/// `profile_xy` is a closed polygon in the local cross-section plane
/// (its own XY); it is auto-closed. The spine is a polyline of 3D
/// points. The cross-section is parallel-transported along the spine
/// so it stays perpendicular to the local tangent without
/// accumulating roll.
///
/// If `spine` is a single straight segment along +Z from the origin
/// the call delegates to [`valenx_cad::prism`] and returns a true
/// BRep solid. Every other (general) spine returns a mesh-backed
/// [`Solid`].
///
/// # Errors
///
/// - [`OcctSurfaceError::BadInput`] for fewer than 3 profile points,
///   fewer than 2 spine points, non-finite inputs, or a spine with a
///   zero-length segment.
/// - [`OcctSurfaceError::TruckLimit`] when the prism builder rejects
///   a straight-line case.
pub fn sweep_api_pipe(
    profile_xy: &[(f64, f64)],
    spine: &[[f64; 3]],
) -> Result<Solid, OcctSurfaceError> {
    if profile_xy.len() < 3 {
        return Err(OcctSurfaceError::bad_input(
            "profile_xy",
            format!("need at least 3 points, got {}", profile_xy.len()),
        ));
    }
    if spine.len() < 2 {
        return Err(OcctSurfaceError::bad_input(
            "spine",
            format!("spine needs at least two points, got {}", spine.len()),
        ));
    }
    for (x, y) in profile_xy {
        if !x.is_finite() || !y.is_finite() {
            return Err(OcctSurfaceError::bad_input(
                "profile_xy",
                "profile contains a non-finite coordinate",
            ));
        }
    }
    for p in spine {
        if p.iter().any(|c| !c.is_finite()) {
            return Err(OcctSurfaceError::bad_input(
                "spine",
                "spine contains a non-finite coordinate",
            ));
        }
    }

    // Easy case: a single +Z segment from the origin is a pure
    // extrusion — delegate to the real BRep prism.
    if spine.len() == 2
        && spine[0] == [0.0, 0.0, 0.0]
        && spine[1][0] == 0.0
        && spine[1][1] == 0.0
        && spine[1][2] > 0.0
    {
        return valenx_cad::prism(profile_xy, spine[1][2])
            .map_err(|e| OcctSurfaceError::TruckLimit(format!("pipe-as-prism: {e:?}")));
    }

    // General curved spine: parallel-transport the profile.
    let stations = transport_frames(spine)?;
    let mesh = sweep_mesh(profile_xy, &stations);
    Ok(Solid::from_mesh(mesh))
}

/// One spine station: the point plus an orthonormal (u, v) frame
/// spanning the cross-section plane (tangent = u × v).
struct Station {
    origin: [f64; 3],
    u: [f64; 3],
    v: [f64; 3],
}

/// Compute a parallel-transport (Bishop / rotation-minimising) frame
/// for every spine vertex. The first frame is seeded arbitrarily;
/// each subsequent frame is the previous one rotated by the minimal
/// rotation carrying the previous tangent onto the current tangent.
/// This avoids the unbounded twist a raw Frenet frame produces at
/// inflection points.
fn transport_frames(spine: &[[f64; 3]]) -> Result<Vec<Station>, OcctSurfaceError> {
    let n = spine.len();
    // Per-vertex tangents: segment direction averaged at interior
    // vertices, single-segment direction at the ends.
    let mut tangents = Vec::with_capacity(n);
    for i in 0..n {
        let incoming = if i > 0 {
            sub(spine[i], spine[i - 1])
        } else {
            [0.0; 3]
        };
        let outgoing = if i + 1 < n {
            sub(spine[i + 1], spine[i])
        } else {
            [0.0; 3]
        };
        let t = match (norm(incoming) > 1e-12, norm(outgoing) > 1e-12) {
            (true, true) => normalize(add(normalize(incoming), normalize(outgoing))),
            (true, false) => normalize(incoming),
            (false, true) => normalize(outgoing),
            (false, false) => {
                return Err(OcctSurfaceError::bad_input(
                    "spine",
                    "spine has a zero-length segment",
                ));
            }
        };
        tangents.push(t);
    }

    let mut stations = Vec::with_capacity(n);
    // Seed frame at vertex 0: pick any vector ⟂ tangent.
    let (mut u, mut v) = perp_basis(tangents[0]);
    stations.push(Station {
        origin: spine[0],
        u,
        v,
    });
    for i in 1..n {
        // Rotate (u, v) by the minimal rotation t[i-1] → t[i].
        let (ru, rv) = rotate_frame(u, v, tangents[i - 1], tangents[i]);
        u = ru;
        v = rv;
        stations.push(Station {
            origin: spine[i],
            u,
            v,
        });
    }
    Ok(stations)
}

/// Build the swept triangle mesh: side walls connecting consecutive
/// cross-section rings plus end caps.
fn sweep_mesh(profile_xy: &[(f64, f64)], stations: &[Station]) -> Mesh {
    let p = profile_xy.len();
    let rings = stations.len();

    // Lay out one ring of `p` nodes per station.
    let mut nodes = Vec::with_capacity(p * rings);
    for st in stations {
        for (x, y) in profile_xy {
            nodes.push(nalgebra::Vector3::new(
                st.origin[0] + st.u[0] * x + st.v[0] * y,
                st.origin[1] + st.u[1] * x + st.v[1] * y,
                st.origin[2] + st.u[2] * x + st.v[2] * y,
            ));
        }
    }

    let mut conn: Vec<u32> = Vec::new();
    // Side walls: quad between (ring r, ring r+1) and (corner k, k+1),
    // split into two triangles.
    for r in 0..rings - 1 {
        let a = (r * p) as u32;
        let b = ((r + 1) * p) as u32;
        for k in 0..p {
            let k1 = ((k + 1) % p) as u32;
            let k0 = k as u32;
            // tri 1: (a+k0, a+k1, b+k1)
            conn.extend_from_slice(&[a + k0, a + k1, b + k1]);
            // tri 2: (a+k0, b+k1, b+k0)
            conn.extend_from_slice(&[a + k0, b + k1, b + k0]);
        }
    }
    // End caps: fan-triangulate the first and last rings.
    let last = ((rings - 1) * p) as u32;
    for k in 1..p - 1 {
        // Start cap (reversed winding so it faces outward / -tangent).
        conn.extend_from_slice(&[0, (k + 1) as u32, k as u32]);
        // End cap.
        conn.extend_from_slice(&[last, last + k as u32, last + (k + 1) as u32]);
    }

    let block = ElementBlock {
        element_type: ElementType::Tri3,
        connectivity: conn,
    };
    let mut mesh = Mesh::new("pipe-sweep");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    mesh
}

// --- small vector helpers ---

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
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
    if l < 1e-20 {
        a
    } else {
        scale(a, 1.0 / l)
    }
}

/// Two orthonormal vectors spanning the plane ⟂ unit vector `t`.
fn perp_basis(t: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let seed = if t[0].abs() <= t[1].abs() && t[0].abs() <= t[2].abs() {
        [1.0, 0.0, 0.0]
    } else if t[1].abs() <= t[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let u = normalize(cross(t, seed));
    let v = cross(t, u);
    (u, v)
}

/// Rotate the orthonormal frame `(u, v)` by the minimal rotation that
/// carries unit tangent `t0` onto unit tangent `t1` (Rodrigues).
fn rotate_frame(u: [f64; 3], v: [f64; 3], t0: [f64; 3], t1: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let axis = cross(t0, t1);
    let sin_a = norm(axis);
    let cos_a = dot(t0, t1).clamp(-1.0, 1.0);
    if sin_a < 1e-12 {
        // Parallel (no rotation) or antiparallel (degenerate — keep
        // the frame; a 180° spine kink is unusual and the visual
        // error is local).
        return (u, v);
    }
    let k = scale(axis, 1.0 / sin_a);
    (rodrigues(u, k, cos_a, sin_a), rodrigues(v, k, cos_a, sin_a))
}

/// Rodrigues rotation of `x` about unit axis `k` by an angle whose
/// cosine/sine are `c`/`s`.
fn rodrigues(x: [f64; 3], k: [f64; 3], c: f64, s: f64) -> [f64; 3] {
    let kx = cross(k, x);
    let kdotx = dot(k, x);
    [
        x[0] * c + kx[0] * s + k[0] * kdotx * (1.0 - c),
        x[1] * c + kx[1] * s + k[1] * kdotx * (1.0 - c),
        x[2] * c + kx[2] * s + k[2] * kdotx * (1.0 - c),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pipe_along_z_axis_acts_as_prism() {
        let s = sweep_api_pipe(
            &[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)],
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]],
        )
        .expect("axis-aligned pipe should delegate to prism");
        assert_eq!(s.faces(), 5);
    }

    #[test]
    fn pipe_with_curved_spine_builds_mesh_solid() {
        let s = sweep_api_pipe(
            &[(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)],
            &[[0.0, 0.0, 0.0], [1.0, 0.0, 1.0], [2.0, 1.0, 2.0]],
        )
        .expect("curved spine pipe should build");
        // Mesh-backed solid: BRep face count is 0, but the cached mesh
        // carries real triangles.
        match &s {
            Solid::Mesh(m) => {
                // 4-sided profile, 3 stations → 2 ring gaps × 4 quads
                // × 2 tris = 16 wall tris + 2 caps × 2 tris = 20 tris.
                assert_eq!(m.total_elements(), 20);
                assert_eq!(m.nodes.len(), 12);
            }
            Solid::Brep(_) => panic!("curved spine should yield a mesh-backed solid"),
        }
    }

    #[test]
    fn pipe_rejects_too_few_spine_points() {
        let err = sweep_api_pipe(&[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)], &[[0.0; 3]]).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn pipe_rejects_degenerate_spine_segment() {
        let err = sweep_api_pipe(
            &[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)],
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn pipe_rejects_short_profile() {
        let err = sweep_api_pipe(
            &[(0.0, 0.0), (1.0, 0.0)],
            &[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }
}
