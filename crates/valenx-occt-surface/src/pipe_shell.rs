//! Phase 71 — `BRepFill_PipeShell`: profile swept along a path with
//! multiple guide curves.
//!
//! ## What OCCT does
//!
//! `BRepFill_PipeShell` builds a shell (BRep face collection) by
//! sweeping one or more `TopoDS_Wire` profiles along a `TopoDS_Wire`
//! spine. Unlike the basic [`crate::sweep_api_pipe()`], the shell
//! variant supports:
//!
//! - **Multiple profiles** — `Add(profile, ...)` blends between
//!   successive profile shapes along the spine (degree-3 G2 in OCCT
//!   v7+), so the cross-section evolves along the path.
//! - **Guide rails** — `SetMode(guide_wire, with_contact,
//!   with_correction)` constrains the swept profile's roll to track an
//!   auxiliary curve.
//! - **Frame transitions** — `SetTransitionMode(RightCorner |
//!   RoundCorner | Transformed)` for sharp spine angles.
//!
//! Used heavily by turbine-blade modelling and HVAC duct design —
//! anything whose cross-section evolves along an arbitrary 3D path.
//!
//! ## v1 status — real mesh-domain multi-profile sweep
//!
//! Honest implementation (Phase 71.5), built on the shared Bishop /
//! rotation-minimising frame in [`crate::sweep_support`].
//!
//! - **Single profile** — equivalent to [`crate::sweep_api_pipe()`]:
//!   the profile is parallel-transported along the spine.
//! - **Multiple profiles** — each profile is assigned a normalised
//!   arc-length station along the spine (profile 0 at `s = 0`, the
//!   last at `s = 1`, the rest evenly spaced). Every profile is
//!   resampled to a common ring size; at each spine station the
//!   active cross-section is the **linear blend** of the two
//!   bracketing profiles. The cross-section therefore evolves along
//!   the path.
//! - **Guide rails** — when `guide_points` carries at least one
//!   guide polyline, the transported frame's roll is corrected at
//!   every station so the cross-section's `+u` axis faces the matching
//!   point on the first guide (matched by normalised arc length).
//!   This is OCCT's `SetMode(guide, with_contact=false)`.
//!
//! ### Honest scope
//!
//! - Profile blending is **linear** in v1 (OCCT uses a degree-3 G2
//!   blend); good for visualisation and meshing, the smooth-spline
//!   blend is a follow-up.
//! - Only the **first** guide rail is used (for roll). Genuine
//!   multi-guide cross-section *warping* — where two or more rails
//!   independently deform the profile so it touches every rail — is a
//!   Tier-3 constrained-surface problem and is *not* done here.
//! - The result is a mesh-backed [`Solid`]; it carries no BRep
//!   topology, so apply it last in a feature chain.

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::OcctSurfaceError;
use crate::sweep_support::{
    arc_length_param, cross, dot, lerp, norm, normalize, perp_basis, resample_closed_polygon,
    rotate_frame, sample_polyline_at, sub, vertex_tangents, Vec3,
};

/// Ring sample count every profile is resampled to before stitching.
const PIPE_SHELL_RING: usize = 24;

/// Sweep one or more profiles along a spine with optional guide rails.
///
/// `profiles[0]` is the start cross-section, `profiles[N-1]` the end
/// cross-section; intermediate profiles interpolate the evolution.
/// Each profile is a closed 3D polygon. `spine_points` is the sweep
/// path as a polyline. `guide_points`, when present and non-empty,
/// supplies guide polylines whose first entry corrects the swept
/// frame's roll.
///
/// # Errors
///
/// - [`OcctSurfaceError::BadInput`] for an empty profile list, any
///   profile with fewer than 3 points, a spine with fewer than 2
///   points, non-finite coordinates, or a spine with a zero-length
///   segment.
pub fn pipe_shell(
    profiles: &[Vec<[f64; 3]>],
    spine_points: &[[f64; 3]],
    guide_points: Option<&[Vec<[f64; 3]>]>,
) -> Result<Solid, OcctSurfaceError> {
    if profiles.is_empty() {
        return Err(OcctSurfaceError::bad_input(
            "profiles",
            "need at least one cross-section profile",
        ));
    }
    if spine_points.len() < 2 {
        return Err(OcctSurfaceError::bad_input(
            "spine_points",
            "spine needs at least two points",
        ));
    }
    for (i, p) in profiles.iter().enumerate() {
        if p.len() < 3 {
            return Err(OcctSurfaceError::bad_input(
                "profiles",
                format!("profile {i} has fewer than 3 points"),
            ));
        }
        for pt in p {
            if pt.iter().any(|c| !c.is_finite()) {
                return Err(OcctSurfaceError::bad_input(
                    "profiles",
                    format!("profile {i} has a non-finite coordinate"),
                ));
            }
        }
    }
    for pt in spine_points {
        if pt.iter().any(|c| !c.is_finite()) {
            return Err(OcctSurfaceError::bad_input(
                "spine_points",
                "spine contains a non-finite coordinate",
            ));
        }
    }

    // --- normalise the profiles into a common local cross-section ---
    // Each profile is projected into its own 2D frame, resampled to a
    // common ring count, and tagged with its spine arc-length station.
    let n_profiles = profiles.len();
    let mut profile_rings: Vec<Vec<(f64, f64)>> = Vec::with_capacity(n_profiles);
    for p in profiles {
        let local3 = project_to_local_plane(p)?;
        let resampled = resample_closed_polygon(&local3, PIPE_SHELL_RING);
        // The resampled ring is still in the (u, v, 0) local frame:
        // drop the (numerically zero) third coordinate.
        let ring2: Vec<(f64, f64)> =
            resampled.iter().map(|q| (q[0], q[1])).collect();
        profile_rings.push(ring2);
    }
    // Profile stations: evenly spaced in normalised arc length.
    let profile_station = |idx: usize| -> f64 {
        if n_profiles < 2 {
            0.0
        } else {
            idx as f64 / (n_profiles - 1) as f64
        }
    };

    // --- spine frames ---
    let tangents = vertex_tangents(spine_points)?;
    let guide0: Option<&[[f64; 3]]> = guide_points
        .and_then(|g| g.first())
        .filter(|g| g.len() >= 2)
        .map(|g| g.as_slice());
    let stations = build_stations(spine_points, &tangents, guide0);

    // --- sweep: blend the cross-section at every spine station ---
    let rings = stations.len();
    let p = PIPE_SHELL_RING;
    let mut nodes: Vec<nalgebra::Vector3<f64>> = Vec::with_capacity(rings * p);
    for st in &stations {
        let s = st.arc;
        let section = blend_profile(&profile_rings, &profile_station, s);
        for &(x, y) in &section {
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
    // End caps.
    let last = ((rings - 1) * p) as u32;
    for k in 1..p - 1 {
        conn.extend_from_slice(&[0, (k + 1) as u32, k as u32]);
        conn.extend_from_slice(&[last, last + k as u32, last + (k + 1) as u32]);
    }

    let mut mesh = Mesh::new("pipe-shell");
    mesh.nodes = nodes;
    mesh.element_blocks.push(ElementBlock {
        element_type: ElementType::Tri3,
        connectivity: conn,
    });
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// One spine station: origin, an orthonormal `(u, v)` cross-section
/// frame, and the normalised arc-length parameter along the spine.
struct Station {
    origin: Vec3,
    u: Vec3,
    v: Vec3,
    arc: f64,
}

/// Build per-station frames: parallel transport, optionally roll-
/// corrected toward a guide rail.
fn build_stations(
    spine: &[[f64; 3]],
    tangents: &[Vec3],
    guide: Option<&[[f64; 3]]>,
) -> Vec<Station> {
    let n = spine.len();
    let mut stations = Vec::with_capacity(n);
    let (mut u, mut v) = perp_basis(tangents[0]);
    for i in 0..n {
        if i > 0 {
            let (ru, rv) = rotate_frame(u, v, tangents[i - 1], tangents[i]);
            u = ru;
            v = rv;
        }
        let origin: Vec3 = spine[i];
        let arc = arc_length_param(spine, i);
        let (fu, fv) = match guide {
            None => (u, v),
            Some(g) => {
                let aim = sample_polyline_at(g, arc);
                roll_toward(origin, tangents[i], u, v, aim)
            }
        };
        stations.push(Station {
            origin,
            u: fu,
            v: fv,
            arc,
        });
    }
    stations
}

/// Roll the frame about the tangent so `+u` faces `aim` (the
/// guide-rail point). Identical correction to the auxiliary spine in
/// [`crate::sweep_api_pipe_shell`].
fn roll_toward(origin: Vec3, tangent: Vec3, u: Vec3, v: Vec3, aim: Vec3) -> (Vec3, Vec3) {
    let to_aim = sub(aim, origin);
    let along = dot(to_aim, tangent);
    let in_plane = [
        to_aim[0] - tangent[0] * along,
        to_aim[1] - tangent[1] * along,
        to_aim[2] - tangent[2] * along,
    ];
    if norm(in_plane) < 1e-9 {
        return (u, v);
    }
    let new_u = normalize(in_plane);
    let new_v = cross(tangent, new_u);
    (new_u, new_v)
}

/// Blend the cross-section at spine arc-length `s` from the profile
/// stations. Picks the two bracketing profiles and linearly
/// interpolates them vertex-by-vertex (all profiles share the common
/// resampled ring size, so the blend is well-defined).
fn blend_profile(
    rings: &[Vec<(f64, f64)>],
    station_of: &impl Fn(usize) -> f64,
    s: f64,
) -> Vec<(f64, f64)> {
    let n = rings.len();
    if n == 1 {
        return rings[0].clone();
    }
    // Find the bracketing pair [lo, hi] with station_of(lo) <= s.
    let mut lo = 0usize;
    for i in 0..n - 1 {
        if station_of(i + 1) <= s + 1e-12 {
            lo = i + 1;
        }
    }
    let hi = (lo + 1).min(n - 1);
    if lo == hi {
        return rings[lo].clone();
    }
    let s_lo = station_of(lo);
    let s_hi = station_of(hi);
    let t = if (s_hi - s_lo).abs() > 1e-12 {
        ((s - s_lo) / (s_hi - s_lo)).clamp(0.0, 1.0)
    } else {
        0.0
    };
    rings[lo]
        .iter()
        .zip(rings[hi].iter())
        .map(|(&(ax, ay), &(bx, by))| {
            let a3 = lerp([ax, ay, 0.0], [bx, by, 0.0], t);
            (a3[0], a3[1])
        })
        .collect()
}

/// Project a 3D closed profile polygon into its own 2D `(u, v)` plane,
/// returned as `Vec3` with a zero third coordinate so it can be fed to
/// [`resample_closed_polygon`].
fn project_to_local_plane(profile: &[[f64; 3]]) -> Result<Vec<Vec3>, OcctSurfaceError> {
    let origin = profile[0];
    let e0 = sub(profile[1], origin);
    if norm(e0) < 1e-12 {
        return Err(OcctSurfaceError::bad_input(
            "profiles",
            "a profile's first edge is degenerate",
        ));
    }
    let mut normal = [0.0; 3];
    for q in profile.iter().skip(2) {
        let n = cross(e0, sub(*q, origin));
        if norm(n) > 1e-12 {
            normal = normalize(n);
            break;
        }
    }
    if norm(normal) < 1e-12 {
        return Err(OcctSurfaceError::bad_input(
            "profiles",
            "a profile's points are collinear — no cross-section plane",
        ));
    }
    let u = normalize(e0);
    let v = cross(normal, u);
    Ok(profile
        .iter()
        .map(|p| {
            let d = sub(*p, origin);
            [dot(d, u), dot(d, v), 0.0]
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::solid_to_mesh;

    /// A square profile of half-width `half` in the plane `z = z`.
    fn square(z: f64, half: f64) -> Vec<[f64; 3]> {
        vec![
            [-half, -half, z],
            [half, -half, z],
            [half, half, z],
            [-half, half, z],
        ]
    }

    #[test]
    fn pipe_shell_validates_empty_profiles() {
        let err = pipe_shell(&[], &[[0.0; 3], [1.0; 3]], None).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn pipe_shell_rejects_short_spine() {
        let err = pipe_shell(&[square(0.0, 1.0)], &[[0.0; 3]], None).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn pipe_shell_rejects_thin_profile() {
        let thin = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]; // only 2 points
        let err = pipe_shell(&[thin], &[[0.0; 3], [0.0, 0.0, 1.0]], None).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn single_profile_sweeps_a_constant_section_tube() {
        let solid = pipe_shell(
            &[square(0.0, 1.0)],
            &[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 2.0]],
            None,
        )
        .expect("single-profile pipe shell should sweep");
        let Solid::Mesh(m) = &solid else {
            panic!("expected mesh solid");
        };
        // 3 stations × PIPE_SHELL_RING nodes.
        assert_eq!(m.nodes.len(), 3 * PIPE_SHELL_RING);
        // It produced side-wall + cap triangles.
        assert!(m.total_elements() > PIPE_SHELL_RING);
    }

    #[test]
    fn two_profiles_evolve_the_cross_section_along_the_spine() {
        // A small square at the spine start, a big square at the end.
        // The swept tube must be narrow near s=0 and wide near s=1.
        let small = square(0.0, 0.5);
        let big = square(0.0, 3.0);
        let spine = vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 2.0],
            [0.0, 0.0, 4.0],
            [0.0, 0.0, 6.0],
        ];
        let solid = pipe_shell(&[small, big], &spine, None).unwrap();
        let mesh = solid_to_mesh(&solid, 0.1).unwrap();
        // Cross-section extent near the bottom (z≈0) vs the top (z≈6).
        let extent_at = |zlo: f64, zhi: f64| -> f64 {
            mesh.nodes
                .iter()
                .filter(|n| n.z >= zlo && n.z <= zhi)
                .map(|n| n.x.abs().max(n.y.abs()))
                .fold(0.0_f64, f64::max)
        };
        let bottom = extent_at(-0.5, 1.0);
        let top = extent_at(5.0, 6.5);
        assert!(
            top > bottom * 2.0,
            "section should widen along the spine: bottom {bottom}, top {top}"
        );
    }

    #[test]
    fn three_profiles_blend_through_an_intermediate_section() {
        // small → big → small. The middle of the spine must be the
        // widest part.
        let profiles = [square(0.0, 0.5), square(0.0, 3.0), square(0.0, 0.5)];
        let spine: Vec<[f64; 3]> =
            (0..=8).map(|k| [0.0, 0.0, k as f64]).collect();
        let solid = pipe_shell(&profiles, &spine, None).unwrap();
        let mesh = solid_to_mesh(&solid, 0.1).unwrap();
        let extent_at = |zlo: f64, zhi: f64| -> f64 {
            mesh.nodes
                .iter()
                .filter(|n| n.z >= zlo && n.z <= zhi)
                .map(|n| n.x.abs().max(n.y.abs()))
                .fold(0.0_f64, f64::max)
        };
        let middle = extent_at(3.5, 4.5);
        let ends = extent_at(-0.5, 0.5).max(extent_at(7.5, 8.5));
        assert!(
            middle > ends * 1.5,
            "the middle section should be widest: middle {middle}, ends {ends}"
        );
    }

    #[test]
    fn guide_rail_rolls_the_cross_section() {
        // A straight +Z spine with a guide rail offset along +X. The
        // first ring's max-X vertex should sit at the profile's
        // half-width along +X (the guide direction).
        let profile = square(0.0, 1.0);
        let spine = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 2.0]];
        let guide = vec![[4.0, 0.0, 0.0], [4.0, 0.0, 1.0], [4.0, 0.0, 2.0]];
        let solid = pipe_shell(&[profile], &spine, Some(&[guide])).unwrap();
        let Solid::Mesh(m) = &solid else {
            panic!("expected mesh solid");
        };
        let first_ring = &m.nodes[0..PIPE_SHELL_RING];
        let max_x = first_ring.iter().map(|n| n.x).fold(f64::MIN, f64::max);
        // The unit square resampled to a ring spans |u| ≤ 1 → max x ≈ 1.
        assert!(
            max_x > 0.9,
            "guide rail should roll +u toward +X, max_x = {max_x}"
        );
    }
}
