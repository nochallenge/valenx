//! Loft evaluator — interpolates a surface between 2+ profile sketches.
//!
//! Phase 13B Task 12; **Phase 13.5** graduated the simple case to a
//! real BRep solid.
//!
//! ## Phase 13.5 — BRep fast path
//!
//! When the loft is the **simple two-profile case** — exactly two
//! profile sketches, `closed = false`, and both sketches' closed
//! wires have the **same number of edges** — the evaluator builds a
//! genuine [`valenx_cad::Solid::Brep`]: the side wall is a real
//! `truck-modeling` `try_wire_homotopy` shell (one NURBS face per
//! edge pair), and the two ends are real planar `try_attach_plane`
//! caps; the faces assemble into a closed `truck` `Solid`. This loft
//! round-trips through STEP/IGES and composes with downstream BRep
//! booleans.
//!
//! Every other case — 3+ profiles, a `closed` periodic loft, or two
//! profiles whose wires have mismatched edge counts — falls through
//! to the mesh-domain loft below (`try_wire_homotopy` requires equal
//! edge counts, and an N>2-profile loft would need a multi-section
//! BRep skin truck does not expose). The fall-through is silent and
//! lossless: the caller gets real swept geometry either way, only the
//! topology backend differs.
//!
//! ## Mesh-domain loft (the general path)
//!
//! Each profile is sampled into a closed polyline of N points
//! (matching across profiles); connecting "rungs" stitch
//! corresponding samples on adjacent profiles, with triangulated caps
//! at the first and last profile.
//!
//! ## Mesh-path limitations
//!
//! - **Profile sampling:** uses `extract_profile_lines` per profile and
//!   re-samples to a common count. Profiles with different vertex
//!   counts are normalized by linear interpolation. Real CAD lofts
//!   match by parametric arc-length.
//! - **Guide curves:** stored on the [`crate::feature::LoftParams`]
//!   but not consumed yet — Phase 14+ work.
//! - **Closed / ruled flags:** `closed` wraps the last profile back to
//!   the first (no caps); `ruled` is honored (straight rungs only —
//!   no smoothed/swept blend).
//! - **Mesh output:** downstream BRep ops on a mesh-domain loft fail
//!   with [`valenx_cad::CadError::MeshBackedSolid`].

use nalgebra::Vector3;
use truck_modeling::builder;
use truck_modeling::{Point3, Solid as TruckSolid, Wire};
use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_sketch::extrude::extract_profile_lines;

use crate::feature::LoftParams;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Common cross-section sample count — every profile is resampled to
/// this many vertices regardless of its native vertex count.
pub const LOFT_SAMPLES_PER_RING: usize = 32;

/// Evaluate a Loft: sample each profile, build connecting rings, cap
/// the ends if not `closed`.
///
/// **Phase 13.5:** the simple two-profile, equal-edge-count,
/// non-closed case returns a real [`Solid::Brep`] via
/// [`try_brep_loft`]; everything else uses the mesh-domain path.
pub(crate) fn evaluate(tree: &FeatureTree, p: &LoftParams) -> Result<Solid, FeatureError> {
    if p.profile_sketches.len() < 2 {
        return Err(FeatureError::BadParameter {
            name: "profile_sketches",
            reason: format!(
                "loft requires at least 2 profile sketches, got {}",
                p.profile_sketches.len()
            ),
        });
    }

    // Phase 13.5 BRep fast path — a genuine truck loft solid for the
    // simple two-profile case. `None` means the case did not qualify
    // (3+ profiles, closed, or mismatched edge counts) — fall through
    // to the mesh-domain loft.
    if !p.closed && p.profile_sketches.len() == 2 {
        if let Some(brep) = try_brep_loft(tree, p)? {
            return Ok(brep);
        }
    }

    // Sample each profile to LOFT_SAMPLES_PER_RING points around its
    // perimeter, distributed along the profile's z-elevation that
    // grows monotonically with profile index.
    let mut rings: Vec<Vec<Vector3<f64>>> = Vec::with_capacity(p.profile_sketches.len());
    for (i, sref) in p.profile_sketches.iter().enumerate() {
        let sketch = tree.get_sketch(*sref)?;
        let waypoints = extract_profile_lines(sketch, 1e-6)?;
        if waypoints.len() < 3 {
            return Err(FeatureError::EmptyProfile);
        }
        // v1: stack profiles along Z by their index (1 unit apart).
        let z = i as f64;
        rings.push(resample_ring_xy(&waypoints, LOFT_SAMPLES_PER_RING, z));
    }

    // Build the side mesh: for each pair of adjacent rings, connect
    // corresponding samples with two triangles per quad.
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut conn: Vec<u32> = Vec::new();
    let n = LOFT_SAMPLES_PER_RING;

    // Push every ring's nodes into the flat buffer, ring-major order.
    for ring in &rings {
        for p in ring {
            nodes.push(*p);
        }
    }

    let pairs: Vec<(usize, usize)> = if p.closed {
        // Connect each adjacent pair AND wrap last→first.
        let mut v: Vec<(usize, usize)> = (0..rings.len() - 1).map(|i| (i, i + 1)).collect();
        v.push((rings.len() - 1, 0));
        v
    } else {
        (0..rings.len() - 1).map(|i| (i, i + 1)).collect()
    };

    for (a, b) in pairs {
        let base_a = a * n;
        let base_b = b * n;
        for k in 0..n {
            let k1 = (k + 1) % n;
            // Triangle 1: a[k], b[k], a[k+1]
            conn.push((base_a + k) as u32);
            conn.push((base_b + k) as u32);
            conn.push((base_a + k1) as u32);
            // Triangle 2: a[k+1], b[k], b[k+1]
            conn.push((base_a + k1) as u32);
            conn.push((base_b + k) as u32);
            conn.push((base_b + k1) as u32);
        }
    }

    // Caps for non-closed lofts — fan-triangulate the first and last
    // ring around its centroid. Skips for `closed = true` because the
    // shape is periodic with no boundary.
    if !p.closed {
        // First ring (z = 0) cap.
        let first_centroid = centroid(&rings[0]);
        let first_idx = nodes.len() as u32;
        nodes.push(first_centroid);
        for k in 0..n {
            let k1 = (k + 1) % n;
            conn.push(first_idx);
            conn.push((k1) as u32);
            conn.push((k) as u32);
        }
        // Last ring cap.
        let last_centroid = centroid(rings.last().unwrap());
        let last_idx = nodes.len() as u32;
        nodes.push(last_centroid);
        let base = (rings.len() - 1) * n;
        for k in 0..n {
            let k1 = (k + 1) % n;
            conn.push(last_idx);
            conn.push((base + k) as u32);
            conn.push((base + k1) as u32);
        }
    }

    let mut mesh = Mesh::new(format!("loft_{}_profiles", rings.len()));
    mesh.nodes = nodes;
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = conn;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    // The `ruled` flag is honored implicitly by straight quadrilaterals
    // between rings; smoothed lofts would interpolate intermediate
    // rings via a B-spline of the rings. v1 always emits the ruled
    // form.
    let _ = p.ruled;
    Ok(Solid::from_mesh(mesh))
}

/// Try to build the simple two-profile loft as a real BRep solid.
///
/// Returns:
/// - `Ok(Some(solid))` — a genuine [`Solid::Brep`] loft.
/// - `Ok(None)` — the inputs do not qualify for the BRep path (the
///   two wires have different edge counts, or a profile is degenerate
///   in a way the mesh path can still handle by resampling). The
///   caller falls through to the mesh-domain loft.
/// - `Err(_)` — a hard error (unknown sketch, empty profile) that
///   both paths would reject anyway.
///
/// Requires `p.profile_sketches.len() == 2`.
///
/// The two profiles are placed at `z = 0` and `z = 1` (matching the
/// mesh path's index-stacking), each built as a closed `truck` wire.
/// `try_wire_homotopy` connects them into a NURBS-faced side shell;
/// `try_attach_plane` caps both ends; the faces assemble into a
/// closed [`TruckSolid`].
fn try_brep_loft(tree: &FeatureTree, p: &LoftParams) -> Result<Option<Solid>, FeatureError> {
    debug_assert_eq!(p.profile_sketches.len(), 2);

    // Extract both profiles as (x, y) waypoint loops.
    let mut profiles: Vec<Vec<(f64, f64)>> = Vec::with_capacity(2);
    for sref in &p.profile_sketches {
        let sketch = tree.get_sketch(*sref)?;
        let wp = extract_profile_lines(sketch, 1e-6)?;
        if wp.len() < 3 {
            return Err(FeatureError::EmptyProfile);
        }
        profiles.push(wp);
    }

    // The BRep homotopy needs the two wires to have the SAME edge
    // count. If they differ, the BRep path cannot apply — let the
    // mesh path resample them to a common ring instead.
    if profiles[0].len() != profiles[1].len() {
        return Ok(None);
    }

    // Build a closed truck wire for each profile at its z-elevation.
    let wire0 = closed_wire(&profiles[0], 0.0);
    let wire1 = closed_wire(&profiles[1], 1.0);

    // The side wall: a homotopy shell, one NURBS face per edge pair.
    let mut shell = match builder::try_wire_homotopy(&wire0, &wire1) {
        Ok(s) => s,
        // truck refused the homotopy (e.g. a degenerate edge) — fall
        // back to the mesh loft rather than erroring.
        Err(_) => return Ok(None),
    };

    // Cap both open boundary loops with planar faces. The homotopy
    // shell's two boundaries are the bottom and top profile wires.
    let boundaries = shell.extract_boundaries();
    if boundaries.len() != 2 {
        // Not a clean two-boundary tube — the mesh path is safer.
        return Ok(None);
    }
    for boundary in boundaries {
        // The cap face must be oriented opposite the boundary so the
        // closed solid has consistent outward normals; truck's
        // `try_attach_plane` builds the face from the wire as given,
        // and we invert the boundary so the cap seals the shell.
        let cap_wire: Wire = boundary.inverse();
        match builder::try_attach_plane(&[cap_wire]) {
            Ok(face) => shell.push(face),
            Err(_) => return Ok(None),
        }
    }

    // The capped shell should now be a closed 2-manifold — a closed
    // shell has no open boundary edges left.
    if !shell.extract_boundaries().is_empty() {
        return Ok(None);
    }
    let solid = TruckSolid::new(vec![shell]);
    // `ruled` is implicit in the homotopy (a degree-1-in-v
    // interpolation between the two profiles) — record the read.
    let _ = p.ruled;
    Ok(Some(Solid::from_truck(solid)))
}

/// Build a closed `truck` wire from an `(x, y)` waypoint loop lifted
/// to `z`. The loop is auto-closed (a final edge links the last
/// waypoint back to the first).
fn closed_wire(waypoints: &[(f64, f64)], z: f64) -> Wire {
    let verts: Vec<_> = waypoints
        .iter()
        .map(|&(x, y)| builder::vertex(Point3::new(x, y, z)))
        .collect();
    let mut edges = Vec::with_capacity(verts.len());
    for i in 0..verts.len() {
        let next = (i + 1) % verts.len();
        edges.push(builder::line(&verts[i], &verts[next]));
    }
    edges.into()
}

/// Resample a closed polyline of (x, y) waypoints onto exactly `n`
/// evenly-distributed points (by chord-length), lifted to `z`.
fn resample_ring_xy(waypoints: &[(f64, f64)], n: usize, z: f64) -> Vec<Vector3<f64>> {
    // Cumulative chord lengths.
    let mut cum: Vec<f64> = vec![0.0];
    let m = waypoints.len();
    for i in 1..m {
        let (px, py) = waypoints[i - 1];
        let (qx, qy) = waypoints[i];
        let d = ((qx - px).powi(2) + (qy - py).powi(2)).sqrt();
        cum.push(cum.last().unwrap() + d);
    }
    // Close the loop if last != first.
    let closing = {
        let (px, py) = waypoints[m - 1];
        let (qx, qy) = waypoints[0];
        ((qx - px).powi(2) + (qy - py).powi(2)).sqrt()
    };
    cum.push(cum.last().unwrap() + closing);
    let total = *cum.last().unwrap();
    if total < 1e-12 {
        return vec![Vector3::new(waypoints[0].0, waypoints[0].1, z); n];
    }

    let mut out = Vec::with_capacity(n);
    for k in 0..n {
        let target = total * (k as f64) / (n as f64);
        // Find the segment containing `target`.
        let mut seg = 0;
        while seg + 1 < cum.len() && cum[seg + 1] < target {
            seg += 1;
        }
        let seg_start = cum[seg];
        let seg_end = if seg + 1 < cum.len() {
            cum[seg + 1]
        } else {
            seg_start
        };
        let len = (seg_end - seg_start).max(1e-12);
        let t = ((target - seg_start) / len).clamp(0.0, 1.0);
        let (sx, sy) = waypoints[seg % m];
        let next_idx = (seg + 1) % m;
        let (ex, ey) = waypoints[next_idx];
        out.push(Vector3::new(sx + t * (ex - sx), sy + t * (ey - sy), z));
    }
    out
}

fn centroid(ring: &[Vector3<f64>]) -> Vector3<f64> {
    if ring.is_empty() {
        return Vector3::zeros();
    }
    let mut acc = Vector3::zeros();
    for p in ring {
        acc += *p;
    }
    acc / ring.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::SketchRef;

    fn square_sketch(half: f64) -> valenx_sketch::Sketch {
        let mut s = valenx_sketch::Sketch::new();
        let a = s.add_point(-half, -half);
        let b = s.add_point(half, -half);
        let c = s.add_point(half, half);
        let d = s.add_point(-half, half);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s.add_line(c, d).unwrap();
        s.add_line(d, a).unwrap();
        s
    }

    /// An `n`-gon sketch of radius `r` — used to exercise the
    /// mismatched-edge-count fall-through.
    fn ngon_sketch(n: usize, r: f64) -> valenx_sketch::Sketch {
        use std::f64::consts::TAU;
        let mut s = valenx_sketch::Sketch::new();
        let mut ids = Vec::with_capacity(n);
        for i in 0..n {
            let a = i as f64 / n as f64 * TAU;
            ids.push(s.add_point(r * a.cos(), r * a.sin()));
        }
        for i in 0..n {
            s.add_line(ids[i], ids[(i + 1) % n]).unwrap();
        }
        s
    }

    #[test]
    fn loft_between_two_squares_is_a_real_brep() {
        // Phase 13.5: two equal-edge-count profiles, not closed →
        // a genuine BRep solid (homotopy side + planar caps).
        let mut tree = FeatureTree::new();
        let s0 = tree.add_sketch(square_sketch(1.0));
        let s1 = tree.add_sketch(square_sketch(0.5));
        let params = LoftParams {
            profile_sketches: vec![s0, s1],
            guide_curves: vec![],
            closed: false,
            ruled: true,
        };
        let solid = evaluate(&tree, &params).expect("loft succeeds");
        assert!(
            matches!(solid, Solid::Brep(_)),
            "the simple two-profile loft should graduate to a BRep solid"
        );
        // A square-to-square loft has 4 side faces + 2 planar caps.
        assert_eq!(solid.faces(), 6, "expected 4 side + 2 cap faces");
        // It tessellates and spans z ∈ [0, 1].
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.1).unwrap();
        assert!(mesh.total_elements() > 0);
        let zmin = mesh.nodes.iter().map(|n| n.z).fold(f64::INFINITY, f64::min);
        let zmax = mesh
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(zmin.abs() < 1e-6 && (zmax - 1.0).abs() < 1e-6);
    }

    #[test]
    fn brep_loft_composes_with_a_boolean() {
        // A BRep loft must round-trip through a downstream boolean —
        // the whole point of graduating it off the mesh domain.
        let mut tree = FeatureTree::new();
        let s0 = tree.add_sketch(square_sketch(1.0));
        let s1 = tree.add_sketch(square_sketch(1.0));
        let params = LoftParams {
            profile_sketches: vec![s0, s1],
            guide_curves: vec![],
            closed: false,
            ruled: true,
        };
        let loft = evaluate(&tree, &params).expect("loft succeeds");
        assert!(matches!(loft, Solid::Brep(_)));
        // Union the loft with a box that overlaps its *interior*. The
        // box is positioned so it shares NO coplanar face with the loft
        // (the loft spans z ∈ [0, 1]; the box spans z ∈ [0.3, 1.5], its
        // caps clear of both loft caps). `truck_shapeops::or` returns
        // `None` whenever two operands share a coplanar face — see the
        // boolean-robustness note in `valenx_cad::boolean` — so a box
        // sharing the loft's z = 0 cap would trip that genuine truck
        // limitation rather than test loft↔boolean composition.
        let box_s = valenx_cad::box_solid(0.6, 0.6, 1.2)
            .unwrap()
            .translated(0.2, 0.2, 0.3)
            .unwrap();
        let fused = valenx_cad::union(&loft, &box_s).expect("BRep union of the loft");
        assert!(fused.faces() > 0, "the fused solid should have faces");
    }

    #[test]
    fn loft_with_mismatched_edge_counts_falls_back_to_mesh() {
        // A 4-gon and a 6-gon cannot homotopy (different edge counts)
        // — the loft must fall through to the mesh-domain path.
        let mut tree = FeatureTree::new();
        let s0 = tree.add_sketch(ngon_sketch(4, 1.0));
        let s1 = tree.add_sketch(ngon_sketch(6, 1.0));
        let params = LoftParams {
            profile_sketches: vec![s0, s1],
            guide_curves: vec![],
            closed: false,
            ruled: true,
        };
        let solid = evaluate(&tree, &params).expect("loft succeeds");
        assert!(
            matches!(solid, Solid::Mesh(_)),
            "mismatched edge counts should fall back to the mesh loft"
        );
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.1).unwrap();
        assert!(mesh.total_elements() > 0);
    }

    #[test]
    fn three_profile_loft_stays_mesh_domain() {
        // 3+ profiles need a multi-section BRep skin truck does not
        // expose — the mesh-domain loft handles it.
        let mut tree = FeatureTree::new();
        let s0 = tree.add_sketch(square_sketch(1.0));
        let s1 = tree.add_sketch(square_sketch(2.0));
        let s2 = tree.add_sketch(square_sketch(1.0));
        let params = LoftParams {
            profile_sketches: vec![s0, s1, s2],
            guide_curves: vec![],
            closed: false,
            ruled: true,
        };
        let solid = evaluate(&tree, &params).expect("loft succeeds");
        assert!(
            matches!(solid, Solid::Mesh(_)),
            "a 3-profile loft stays mesh-domain"
        );
    }

    #[test]
    fn loft_rejects_single_profile() {
        let mut tree = FeatureTree::new();
        let s0 = tree.add_sketch(square_sketch(1.0));
        let params = LoftParams {
            profile_sketches: vec![s0],
            guide_curves: vec![],
            closed: false,
            ruled: true,
        };
        let err = evaluate(&tree, &params).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter {
                name: "profile_sketches",
                ..
            }
        ));
    }

    #[test]
    fn loft_rejects_missing_sketch() {
        let tree = FeatureTree::new();
        let params = LoftParams {
            profile_sketches: vec![SketchRef(0), SketchRef(1)],
            guide_curves: vec![],
            closed: false,
            ruled: true,
        };
        let err = evaluate(&tree, &params).unwrap_err();
        assert_eq!(err.code(), "feature_tree.unknown_sketch");
    }
}
