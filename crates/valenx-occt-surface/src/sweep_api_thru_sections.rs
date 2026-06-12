//! Phase 90 — `BRepOffsetAPI_ThruSections` (loft through profiles).
//!
//! ## What OCCT does
//!
//! `BRepOffsetAPI_ThruSections` builds a lofted surface (or solid)
//! through a sequence of profile wires. Unlike a sweep, there's no
//! spine — the algorithm interpolates between successive profiles
//! using a tensor-product B-spline surface (degree 3 by default). The
//! caller can configure:
//!
//! - `IsSolid(bool)` — return a `TopoDS_Solid` (closing the ends
//!   with planar caps) vs a `TopoDS_Shell`.
//! - `IsRuled(bool)` — use a ruled (linear-in-v) interpolation
//!   between adjacent sections, sacrificing smoothness for
//!   predictability.
//! - `SetCriteriumWeight(w1, w2, w3)` — energy-minimisation weights
//!   for the smooth interpolation case.
//!
//! Used as the bread-and-butter surfacing tool: hull plating,
//! aerofoil sections, swimming-pool form, anything described as a
//! pile of cross-sections.
//!
//! ## v1 status — real loft
//!
//! This is a genuine loft, not a stub. The profiles are closed 3D
//! polygons; the algorithm:
//!
//! 1. **Resamples** every profile to a common vertex count
//!    (`LOFT_RING_SAMPLES`) by arc-length sampling around the profile
//!    polygon — so profiles with different point counts still stitch.
//! 2. For `is_ruled = true` connects consecutive profile rings
//!    directly with quad side walls (ruled / linear-in-v).
//! 3. For `is_ruled = false` first **resamples in the loft (v)
//!    direction** through a Catmull-Rom spline so the side surface is
//!    smooth across the section joints, then stitches the denser ring
//!    set.
//! 4. For `is_solid = true` triangulates a planar cap at the first
//!    and last profile.
//!
//! The result is a watertight (when `is_solid`) mesh-backed
//! [`Solid`]. A mesh-backed solid carries no BRep topology, so
//! downstream booleans refuse it — apply the loft last in a feature
//! chain (the same rule as the sweep / fillet pipelines).

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::OcctSurfaceError;

/// Vertices each profile is resampled to before stitching.
pub const LOFT_RING_SAMPLES: usize = 24;

/// Intermediate rings inserted *between* each pair of profiles when
/// `is_ruled = false` (smooth loft) — higher = smoother side surface.
const SMOOTH_SUBDIVISIONS: usize = 4;

/// Loft a surface through a sequence of profile wires.
///
/// Each profile is a closed polygon in 3D. `is_solid` controls
/// whether the result has planar end caps; `is_ruled` controls the
/// interpolation style (linear vs. smooth Catmull-Rom).
///
/// # Errors
///
/// [`OcctSurfaceError::BadInput`] for fewer than 2 profiles, a
/// profile with fewer than 3 points, or non-finite coordinates.
pub fn sweep_api_thru_sections(
    profiles: &[Vec<[f64; 3]>],
    is_solid: bool,
    is_ruled: bool,
) -> Result<Solid, OcctSurfaceError> {
    if profiles.len() < 2 {
        return Err(OcctSurfaceError::bad_input(
            "profiles",
            format!(
                "need at least 2 profile wires to loft, got {}",
                profiles.len()
            ),
        ));
    }
    for (idx, p) in profiles.iter().enumerate() {
        if p.len() < 3 {
            return Err(OcctSurfaceError::bad_input(
                "profiles",
                format!("profile {idx} has fewer than 3 points"),
            ));
        }
        for pt in p {
            if pt.iter().any(|c| !c.is_finite()) {
                return Err(OcctSurfaceError::bad_input(
                    "profiles",
                    format!("profile {idx} has a non-finite coordinate"),
                ));
            }
        }
    }

    // Resample every profile to a common ring of LOFT_RING_SAMPLES
    // vertices so different point counts still stitch cleanly.
    let mut rings: Vec<Vec<[f64; 3]>> = profiles
        .iter()
        .map(|p| resample_closed_polygon(p, LOFT_RING_SAMPLES))
        .collect();

    // Smooth loft: insert Catmull-Rom-interpolated intermediate rings.
    if !is_ruled {
        rings = smooth_rings(&rings);
    }

    let mesh = loft_mesh(&rings, is_solid);
    Ok(Solid::from_mesh(mesh))
}

/// Resample a closed polygon to exactly `n` vertices spaced by equal
/// arc length around its perimeter.
fn resample_closed_polygon(poly: &[[f64; 3]], n: usize) -> Vec<[f64; 3]> {
    // Cumulative arc length around the closed loop.
    let m = poly.len();
    let mut seg_len = Vec::with_capacity(m);
    let mut total = 0.0;
    for i in 0..m {
        let a = poly[i];
        let b = poly[(i + 1) % m];
        let d = dist(a, b);
        seg_len.push(d);
        total += d;
    }
    if total < 1e-12 {
        // Degenerate profile — just repeat the first vertex.
        return vec![poly[0]; n];
    }
    let step = total / n as f64;
    let mut out = Vec::with_capacity(n);
    let mut seg = 0usize;
    let mut seg_start = 0.0;
    for k in 0..n {
        let target = k as f64 * step;
        // Advance to the segment containing `target`.
        while seg + 1 < m && seg_start + seg_len[seg] < target {
            seg_start += seg_len[seg];
            seg += 1;
        }
        let local = if seg_len[seg] > 1e-12 {
            (target - seg_start) / seg_len[seg]
        } else {
            0.0
        };
        let a = poly[seg];
        let b = poly[(seg + 1) % m];
        out.push(lerp(a, b, local.clamp(0.0, 1.0)));
    }
    out
}

/// Insert smooth intermediate rings between each pair of profile
/// rings via a per-vertex Catmull-Rom spline through the profiles.
fn smooth_rings(profiles: &[Vec<[f64; 3]>]) -> Vec<Vec<[f64; 3]>> {
    let n_profiles = profiles.len();
    let ring_size = profiles[0].len();
    let mut out: Vec<Vec<[f64; 3]>> = Vec::new();
    for seg in 0..n_profiles - 1 {
        // Control profiles for this Catmull-Rom segment.
        let p0 = &profiles[seg.saturating_sub(1)];
        let p1 = &profiles[seg];
        let p2 = &profiles[seg + 1];
        let p3 = &profiles[(seg + 2).min(n_profiles - 1)];
        // Push the start profile, then SMOOTH_SUBDIVISIONS interior
        // rings; the last segment additionally pushes the end profile.
        let steps = SMOOTH_SUBDIVISIONS + 1;
        for s in 0..steps {
            let t = s as f64 / steps as f64;
            let mut ring = Vec::with_capacity(ring_size);
            for v in 0..ring_size {
                ring.push(catmull_rom(p0[v], p1[v], p2[v], p3[v], t));
            }
            out.push(ring);
        }
    }
    // Close with the final profile.
    out.push(profiles[n_profiles - 1].clone());
    out
}

/// Catmull-Rom interpolation of one vertex across four control rings.
fn catmull_rom(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], t: f64) -> [f64; 3] {
    let t2 = t * t;
    let t3 = t2 * t;
    let mut out = [0.0; 3];
    for c in 0..3 {
        out[c] = 0.5
            * ((2.0 * p1[c])
                + (-p0[c] + p2[c]) * t
                + (2.0 * p0[c] - 5.0 * p1[c] + 4.0 * p2[c] - p3[c]) * t2
                + (-p0[c] + 3.0 * p1[c] - 3.0 * p2[c] + p3[c]) * t3);
    }
    out
}

/// Build the lofted triangle mesh from a stack of equal-size rings.
fn loft_mesh(rings: &[Vec<[f64; 3]>], is_solid: bool) -> Mesh {
    let mut mesh = Mesh::new("loft");
    let ring_size = rings[0].len();
    // Append all ring vertices.
    for ring in rings {
        for v in ring {
            mesh.nodes.push(nalgebra::Vector3::new(v[0], v[1], v[2]));
        }
    }
    let mut conn: Vec<u32> = Vec::new();
    // Side walls: stitch consecutive rings with quads (two triangles).
    for r in 0..rings.len() - 1 {
        let base_a = (r * ring_size) as u32;
        let base_b = ((r + 1) * ring_size) as u32;
        for k in 0..ring_size {
            let kn = ((k + 1) % ring_size) as u32;
            let a0 = base_a + k as u32;
            let a1 = base_a + kn;
            let b0 = base_b + k as u32;
            let b1 = base_b + kn;
            conn.extend_from_slice(&[a0, a1, b1]);
            conn.extend_from_slice(&[a0, b1, b0]);
        }
    }
    // End caps for a solid loft — triangle fans around each ring's
    // centroid.
    if is_solid {
        // First ring cap (wound so the normal faces outward, away from
        // the loft body).
        let c0 = ring_centroid(&rings[0]);
        let c0_idx = mesh.nodes.len() as u32;
        mesh.nodes.push(nalgebra::Vector3::new(c0[0], c0[1], c0[2]));
        for k in 0..ring_size {
            let kn = ((k + 1) % ring_size) as u32;
            conn.extend_from_slice(&[c0_idx, kn, k as u32]);
        }
        // Last ring cap.
        let last = rings.len() - 1;
        let cl = ring_centroid(&rings[last]);
        let cl_idx = mesh.nodes.len() as u32;
        mesh.nodes.push(nalgebra::Vector3::new(cl[0], cl[1], cl[2]));
        let base = (last * ring_size) as u32;
        for k in 0..ring_size {
            let kn = ((k + 1) % ring_size) as u32;
            conn.extend_from_slice(&[cl_idx, base + k as u32, base + kn]);
        }
    }
    mesh.element_blocks.push(ElementBlock {
        element_type: ElementType::Tri3,
        connectivity: conn,
    });
    mesh.recompute_stats();
    mesh
}

/// Centroid of a ring of points.
fn ring_centroid(ring: &[[f64; 3]]) -> [f64; 3] {
    let mut c = [0.0; 3];
    for p in ring {
        for k in 0..3 {
            c[k] += p[k];
        }
    }
    let n = ring.len().max(1) as f64;
    [c[0] / n, c[1] / n, c[2] / n]
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_cad::solid_to_mesh;

    /// A square profile at height `z`, side `2·half`.
    fn square(z: f64, half: f64) -> Vec<[f64; 3]> {
        vec![
            [-half, -half, z],
            [half, -half, z],
            [half, half, z],
            [-half, half, z],
        ]
    }

    #[test]
    fn loft_rejects_single_profile() {
        let err = sweep_api_thru_sections(
            &[vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]]],
            true,
            false,
        )
        .unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn loft_rejects_degenerate_profile() {
        let p0 = square(0.0, 1.0);
        let p1 = vec![[0.0, 0.0, 1.0], [1.0, 0.0, 1.0]]; // only 2 points
        let err = sweep_api_thru_sections(&[p0, p1], true, true).unwrap_err();
        assert_eq!(err.code(), "occt_surface.bad_input");
    }

    #[test]
    fn ruled_loft_of_two_squares_is_a_box_shell() {
        // A unit square at z=0 and z=2, ruled loft → a box-like shell.
        let p0 = square(0.0, 1.0);
        let p1 = square(2.0, 1.0);
        let solid = sweep_api_thru_sections(&[p0, p1], true, true).unwrap();
        let mesh = solid_to_mesh(&solid, 0.1).unwrap();
        // The mesh spans z ∈ [0, 2].
        let zmin = mesh.nodes.iter().map(|n| n.z).fold(f64::INFINITY, f64::min);
        let zmax = mesh
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((zmin - 0.0).abs() < 1e-6, "zmin={zmin}");
        assert!((zmax - 2.0).abs() < 1e-6, "zmax={zmax}");
        // It has side-wall + cap triangles.
        assert!(mesh.element_blocks[0].connectivity.len() / 3 > LOFT_RING_SAMPLES);
    }

    #[test]
    fn loft_widens_between_profiles() {
        // A small square at z=0 and a big square at z=4 — a frustum.
        // The mesh's XY extent at the top must exceed the bottom.
        let p0 = square(0.0, 1.0);
        let p1 = square(4.0, 3.0);
        let solid = sweep_api_thru_sections(&[p0, p1], false, true).unwrap();
        let mesh = solid_to_mesh(&solid, 0.1).unwrap();
        let bottom_extent = mesh
            .nodes
            .iter()
            .filter(|n| n.z < 0.5)
            .map(|n| n.x.abs())
            .fold(0.0_f64, f64::max);
        let top_extent = mesh
            .nodes
            .iter()
            .filter(|n| n.z > 3.5)
            .map(|n| n.x.abs())
            .fold(0.0_f64, f64::max);
        assert!(
            top_extent > bottom_extent * 1.5,
            "top {top_extent} should dwarf bottom {bottom_extent}"
        );
    }

    #[test]
    fn smooth_loft_inserts_intermediate_rings() {
        // Three profiles, smooth (non-ruled) loft → more rings than
        // the ruled case → strictly more triangles.
        let profiles = [square(0.0, 1.0), square(2.0, 2.0), square(4.0, 1.0)];
        let ruled = sweep_api_thru_sections(&profiles, false, true).unwrap();
        let smooth = sweep_api_thru_sections(&profiles, false, false).unwrap();
        let ruled_tris = solid_to_mesh(&ruled, 0.1).unwrap().element_blocks[0]
            .connectivity
            .len();
        let smooth_tris = solid_to_mesh(&smooth, 0.1).unwrap().element_blocks[0]
            .connectivity
            .len();
        assert!(
            smooth_tris > ruled_tris,
            "smooth loft ({smooth_tris}) should have more triangles than ruled ({ruled_tris})"
        );
    }

    #[test]
    fn resample_gives_requested_count() {
        let poly = square(0.0, 1.0);
        let rs = resample_closed_polygon(&poly, LOFT_RING_SAMPLES);
        assert_eq!(rs.len(), LOFT_RING_SAMPLES);
        // Every resampled point still lies on the square's perimeter
        // (max |x| or |y| == 1 for a unit square).
        for p in &rs {
            let on_edge = (p[0].abs() - 1.0).abs() < 1e-6 || (p[1].abs() - 1.0).abs() < 1e-6;
            assert!(on_edge, "resampled point off the perimeter: {p:?}");
        }
    }

    #[test]
    fn loft_handles_profiles_with_different_point_counts() {
        // A 4-gon and an 8-gon — resampling makes them stitchable.
        let square4 = square(0.0, 1.0);
        let octagon: Vec<[f64; 3]> = (0..8)
            .map(|k| {
                let a = std::f64::consts::TAU * k as f64 / 8.0;
                [a.cos() * 2.0, a.sin() * 2.0, 3.0]
            })
            .collect();
        let solid = sweep_api_thru_sections(&[square4, octagon], true, true).unwrap();
        let mesh = solid_to_mesh(&solid, 0.1).unwrap();
        assert!(!mesh.nodes.is_empty());
    }
}
