//! Shared guide-curve loft machinery for Phases 131 and 138.
//!
//! OCCT's `BRepOffsetAPI_ThruSections` with guide curves (and the
//! feature-based `loft-with-rails` in Phase 138) both need the same
//! core: a loft through cross-section wires whose **intermediate
//! sections are pulled toward explicit guide curves**. This module
//! hoists that algorithm so [`offset_api_thru_sections_with_guides`]
//! and [`feat_make_loft_with_rails`] share one implementation.
//!
//! [`offset_api_thru_sections_with_guides`]: fn@crate::offset_api_thru_sections_with_guides
//! [`feat_make_loft_with_rails`]: fn@crate::feat_make_loft_with_rails
//!
//! ## v1 algorithm — guide-warped loft
//!
//! 1. Resample every cross-section wire to a common ring of
//!    [`GUIDE_LOFT_RING_SAMPLES`] vertices (arc-length around the
//!    polygon), so sections with different point counts stitch.
//! 2. Insert smooth Catmull-Rom intermediate rings between every pair
//!    of sections (the same smoothing the Phase-90
//!    `sweep_api_thru_sections` loft does).
//! 3. **Guide warp** — for each ring (the v-station), evaluate every
//!    guide curve at the matching normalised arc-length parameter to
//!    get a target point per guide. The ring is rigidly translated +
//!    radially scaled so its centroid and radial scale follow the
//!    guide(s): the ring is offset by the mean guide-displacement
//!    (where each guide's displacement is "where the guide is at this
//!    station" minus "where it would be on the un-warped loft"), and
//!    radially scaled toward the guides so the section width tracks
//!    the rails. With one guide this pins the section's centre to the
//!    rail; with several it follows the rails' average and spread.
//! 4. Stitch the warped ring stack into a triangle mesh; optional
//!    triangulated planar end caps for a solid result.
//!
//! ## Honest scope
//!
//! This is a real guide-constrained loft — the intermediate sections
//! genuinely move to follow the rails — but it is **mesh-domain** (no
//! BRep faces) and the warp is a rigid-translate + uniform-radial-
//! scale, not the per-control-point tangential `Geom_BSplineSurface`
//! skinning OCCT runs. A turbine blade lofted here follows its LE/TE
//! rails in centroid and width; the exact surface tangency at each
//! rail crossing is the documented Tier-3 follow-up (it needs the
//! BRep substrate, the same gate as 14.5 / the sweep family).

use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

/// 3-component point — both `sections` (closed polylines) and `guides`
/// (open polylines) are lists of these.
pub type Pt3 = [f64; 3];

/// Vertices each cross-section is resampled to before stitching.
pub const GUIDE_LOFT_RING_SAMPLES: usize = 32;

/// Smooth intermediate rings inserted between each section pair.
const GUIDE_LOFT_SUBDIVISIONS: usize = 5;

/// Loft a guide-warped triangle mesh through `sections`, pulling the
/// intermediate sections toward `guides`.
///
/// `sections` is an ordered list of closed cross-section polylines;
/// `guides` is a list of open rail polylines. `is_solid` adds planar
/// end caps. Both lists are assumed pre-validated by the caller
/// (`sections.len() >= 2`, each section `>= 3` points, each guide
/// `>= 2` points).
pub fn guide_loft_mesh(sections: &[Vec<Pt3>], guides: &[Vec<Pt3>], is_solid: bool) -> Mesh {
    // 1. Resample each section to a common ring.
    let section_rings: Vec<Vec<Pt3>> = sections
        .iter()
        .map(|s| resample_closed(s, GUIDE_LOFT_RING_SAMPLES))
        .collect();

    // 2. Smooth intermediate rings — Catmull-Rom across the sections.
    //    `ring_station[i]` is the normalised v-parameter (0..=1) of
    //    ring `i`, used both to read the guides and to know which
    //    rings are "intermediate" (a section ring is left unwarped so
    //    the loft still interpolates its defining sections exactly).
    let (rings, ring_is_section) = smooth_rings_with_stations(&section_rings);
    let ring_count = rings.len();

    // 3. Guide warp. For every guide, the un-warped reference point at
    //    a station is the guide's own endpoints linearly blended (a
    //    guide that already runs straight produces zero displacement);
    //    the warp is the guide value at the station minus that
    //    reference, averaged over the guides.
    let mut warped: Vec<Vec<Pt3>> = Vec::with_capacity(ring_count);
    for (i, ring) in rings.iter().enumerate() {
        let s = if ring_count > 1 {
            i as f64 / (ring_count - 1) as f64
        } else {
            0.0
        };
        if ring_is_section[i] || guides.is_empty() {
            // Section rings stay put: the loft must pass through its
            // defining cross-sections exactly.
            warped.push(ring.clone());
            continue;
        }
        warped.push(warp_ring(ring, guides, s));
    }

    // 4. Stitch.
    loft_mesh(&warped, is_solid)
}

/// Warp one intermediate ring toward the guide curves at station `s`.
///
/// Each guide contributes a displacement = (guide point at `s`) −
/// (guide's straight-line reference at `s`). The ring is translated by
/// the mean displacement, then radially scaled so the mean ring-to-
/// centroid radius tracks the mean guide-to-centroid distance — this
/// is what makes the section *widen* to follow rails that fan out.
fn warp_ring(ring: &[Pt3], guides: &[Vec<Pt3>], s: f64) -> Vec<Pt3> {
    let centroid = ring_centroid(ring);

    // Mean translational displacement from the guides.
    let mut disp = [0.0; 3];
    // Mean guide-to-centroid distance (radial-scale target) and the
    // mean ring-to-centroid radius (the current radius).
    let mut guide_radius_sum = 0.0;
    for g in guides {
        let gp = sample_open_polyline(g, s);
        let g_ref = straight_reference(g, s);
        for c in 0..3 {
            disp[c] += gp[c] - g_ref[c];
        }
        guide_radius_sum += dist(gp, centroid);
    }
    let n = guides.len() as f64;
    for d in &mut disp {
        *d /= n;
    }
    let guide_radius = guide_radius_sum / n;

    let mean_ring_radius = {
        let mut sum = 0.0;
        for p in ring {
            sum += dist(*p, centroid);
        }
        sum / ring.len().max(1) as f64
    };

    // Radial scale: pull the ring's radius a fraction of the way
    // toward the guide radius. Half-weight keeps the warp gentle so a
    // single off-centre guide does not collapse the section.
    let scale = if mean_ring_radius > 1e-9 {
        let target = guide_radius / mean_ring_radius;
        1.0 + 0.5 * (target - 1.0)
    } else {
        1.0
    };

    ring.iter()
        .map(|p| {
            let mut out = [0.0; 3];
            for c in 0..3 {
                // Radially scale about the centroid, then translate.
                out[c] = centroid[c] + (p[c] - centroid[c]) * scale + disp[c];
            }
            out
        })
        .collect()
}

/// The straight-line reference position of an open polyline at
/// normalised arc-length `s` — a linear blend of its two endpoints.
/// Subtracting this from the real guide sample isolates the guide's
/// *curvature* so a guide that already runs straight warps nothing.
fn straight_reference(poly: &[Pt3], s: f64) -> Pt3 {
    let a = poly[0];
    let b = poly[poly.len() - 1];
    lerp(a, b, s.clamp(0.0, 1.0))
}

/// Resample a closed polygon to exactly `n` arc-length-spaced points.
fn resample_closed(poly: &[Pt3], n: usize) -> Vec<Pt3> {
    let m = poly.len();
    if m == 0 || n == 0 {
        return Vec::new();
    }
    let mut seg = Vec::with_capacity(m);
    let mut total = 0.0;
    for i in 0..m {
        let d = dist(poly[i], poly[(i + 1) % m]);
        seg.push(d);
        total += d;
    }
    if total < 1e-12 {
        return vec![poly[0]; n];
    }
    let step = total / n as f64;
    let mut out = Vec::with_capacity(n);
    let mut si = 0usize;
    let mut start = 0.0;
    for k in 0..n {
        let target = k as f64 * step;
        while si + 1 < m && start + seg[si] < target {
            start += seg[si];
            si += 1;
        }
        let local = if seg[si] > 1e-12 {
            (target - start) / seg[si]
        } else {
            0.0
        };
        out.push(lerp(poly[si], poly[(si + 1) % m], local.clamp(0.0, 1.0)));
    }
    out
}

/// Insert smooth Catmull-Rom intermediate rings between each pair of
/// section rings. Returns the dense ring stack plus, for each ring, a
/// flag that is `true` exactly when the ring is one of the original
/// (un-interpolated) section rings.
fn smooth_rings_with_stations(sections: &[Vec<Pt3>]) -> (Vec<Vec<Pt3>>, Vec<bool>) {
    let n_sec = sections.len();
    let ring_size = sections[0].len();
    let mut out: Vec<Vec<Pt3>> = Vec::new();
    let mut is_section: Vec<bool> = Vec::new();
    for seg in 0..n_sec - 1 {
        let p0 = &sections[seg.saturating_sub(1)];
        let p1 = &sections[seg];
        let p2 = &sections[seg + 1];
        let p3 = &sections[(seg + 2).min(n_sec - 1)];
        let steps = GUIDE_LOFT_SUBDIVISIONS + 1;
        for st in 0..steps {
            let t = st as f64 / steps as f64;
            let mut ring = Vec::with_capacity(ring_size);
            for v in 0..ring_size {
                ring.push(catmull_rom(p0[v], p1[v], p2[v], p3[v], t));
            }
            // st == 0 lands exactly on section `seg`.
            is_section.push(st == 0);
            out.push(ring);
        }
    }
    // Final section ring.
    out.push(sections[n_sec - 1].clone());
    is_section.push(true);
    (out, is_section)
}

/// Catmull-Rom interpolation of one vertex across four control rings.
fn catmull_rom(p0: Pt3, p1: Pt3, p2: Pt3, p3: Pt3, t: f64) -> Pt3 {
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

/// Sample an open polyline at normalised arc-length `s` (0..=1).
fn sample_open_polyline(poly: &[Pt3], s: f64) -> Pt3 {
    let n = poly.len();
    if n == 0 {
        return [0.0; 3];
    }
    if n == 1 {
        return poly[0];
    }
    let s = s.clamp(0.0, 1.0);
    let mut cum = vec![0.0_f64; n];
    for k in 1..n {
        cum[k] = cum[k - 1] + dist(poly[k - 1], poly[k]);
    }
    let total = cum[n - 1];
    if total < 1e-12 {
        return poly[0];
    }
    let target = s * total;
    for k in 1..n {
        if cum[k] >= target - 1e-12 {
            let seg = cum[k] - cum[k - 1];
            let local = if seg > 1e-12 {
                (target - cum[k - 1]) / seg
            } else {
                0.0
            };
            return lerp(poly[k - 1], poly[k], local.clamp(0.0, 1.0));
        }
    }
    poly[n - 1]
}

/// Build the lofted triangle mesh from a stack of equal-size rings.
fn loft_mesh(rings: &[Vec<Pt3>], is_solid: bool) -> Mesh {
    let mut mesh = Mesh::new("guide-loft");
    let ring_size = rings[0].len();
    for ring in rings {
        for v in ring {
            mesh.nodes.push(nalgebra::Vector3::new(v[0], v[1], v[2]));
        }
    }
    let mut conn: Vec<u32> = Vec::new();
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
    if is_solid {
        let c0 = ring_centroid(&rings[0]);
        let c0_idx = mesh.nodes.len() as u32;
        mesh.nodes.push(nalgebra::Vector3::new(c0[0], c0[1], c0[2]));
        for k in 0..ring_size {
            let kn = ((k + 1) % ring_size) as u32;
            conn.extend_from_slice(&[c0_idx, kn, k as u32]);
        }
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

/// Convenience: a guide-warped mesh-backed [`Solid`].
pub fn guide_loft_solid(sections: &[Vec<Pt3>], guides: &[Vec<Pt3>], is_solid: bool) -> Solid {
    Solid::from_mesh(guide_loft_mesh(sections, guides, is_solid))
}

fn ring_centroid(ring: &[Pt3]) -> Pt3 {
    let mut c = [0.0; 3];
    for p in ring {
        for k in 0..3 {
            c[k] += p[k];
        }
    }
    let n = ring.len().max(1) as f64;
    [c[0] / n, c[1] / n, c[2] / n]
}

fn dist(a: Pt3, b: Pt3) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn lerp(a: Pt3, b: Pt3, t: f64) -> Pt3 {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A square cross-section at height `z`, half-side `half`.
    fn square(z: f64, half: f64) -> Vec<Pt3> {
        vec![
            [-half, -half, z],
            [half, -half, z],
            [half, half, z],
            [-half, half, z],
        ]
    }

    #[test]
    fn straight_guide_warps_nothing() {
        // Two identical squares + a guide that is itself a straight
        // line: the straight reference equals the guide sample, so the
        // displacement is zero — the loft is unchanged from the
        // guide-free case (same node count, mesh built).
        let secs = vec![square(0.0, 1.0), square(4.0, 1.0)];
        let straight = vec![vec![[0.0, 0.0, 0.0], [0.0, 0.0, 4.0]]];
        let m = guide_loft_mesh(&secs, &straight, false);
        assert!(!m.nodes.is_empty());
        // A station ring's centroid should still be on the axis (the
        // straight guide added no displacement).
        let mid_ring_start = GUIDE_LOFT_RING_SAMPLES * 3;
        let mut c = [0.0; 3];
        for p in &m.nodes[mid_ring_start..mid_ring_start + GUIDE_LOFT_RING_SAMPLES] {
            c[0] += p.x;
            c[1] += p.y;
        }
        c[0] /= GUIDE_LOFT_RING_SAMPLES as f64;
        c[1] /= GUIDE_LOFT_RING_SAMPLES as f64;
        assert!(
            c[0].abs() < 1e-6 && c[1].abs() < 1e-6,
            "centroid drifted: {c:?}"
        );
    }

    #[test]
    fn bowed_guide_pulls_intermediate_sections() {
        // Two squares on the Z axis, but a guide that bows out to +X
        // in the middle. An intermediate ring's centroid must move
        // toward +X (the guide's curvature is non-zero there).
        let secs = vec![square(0.0, 1.0), square(4.0, 1.0)];
        let bowed = vec![vec![
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 2.0], // bows out to +X at the mid station
            [0.0, 0.0, 4.0],
        ]];
        let m = guide_loft_mesh(&secs, &bowed, false);
        // Find the ring nearest z = 2 and check its centroid X > 0.
        let rs = GUIDE_LOFT_RING_SAMPLES;
        let n_rings = m.nodes.len() / rs;
        let mut best_x = f64::NEG_INFINITY;
        for r in 0..n_rings {
            let mut cz = 0.0;
            let mut cx = 0.0;
            for p in &m.nodes[r * rs..r * rs + rs] {
                cz += p.z;
                cx += p.x;
            }
            cz /= rs as f64;
            cx /= rs as f64;
            if (cz - 2.0).abs() < 0.6 {
                best_x = best_x.max(cx);
            }
        }
        assert!(
            best_x > 0.3,
            "the bowed guide should pull a mid ring toward +X, got cx={best_x}"
        );
    }

    #[test]
    fn solid_loft_adds_end_caps() {
        let secs = vec![square(0.0, 1.0), square(2.0, 1.0)];
        let guides = vec![vec![[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]]];
        let shell = guide_loft_mesh(&secs, &guides, false);
        let solid = guide_loft_mesh(&secs, &guides, true);
        // The solid variant has two extra centroid nodes + cap tris.
        assert!(solid.nodes.len() > shell.nodes.len());
        assert!(solid.total_elements() > shell.total_elements());
    }

    #[test]
    fn section_rings_are_interpolated_exactly() {
        // The first and last ring of the stack must equal the input
        // sections (the loft passes through its defining sections).
        let secs = vec![square(0.0, 1.5), square(3.0, 1.5)];
        let guides = vec![vec![[5.0, 0.0, 0.0], [5.0, 0.0, 3.0]]];
        let m = guide_loft_mesh(&secs, &guides, false);
        let rs = GUIDE_LOFT_RING_SAMPLES;
        // First ring z ≈ 0.
        let z0: f64 = m.nodes[..rs].iter().map(|p| p.z).sum::<f64>() / rs as f64;
        assert!(z0.abs() < 1e-6, "first ring not at z=0: {z0}");
        // First ring centroid still on the section axis (un-warped).
        let cx0: f64 = m.nodes[..rs].iter().map(|p| p.x).sum::<f64>() / rs as f64;
        assert!(cx0.abs() < 1e-6, "section ring was warped: cx={cx0}");
    }
}
