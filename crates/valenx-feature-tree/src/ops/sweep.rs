//! Sweep evaluator — sweeps a 2D profile along a path.
//!
//! Phase 13B Task 15.
//!
//! ## BRep fast path — a straight path is a true extrusion
//!
//! When the sweep path is a single **straight segment** (exactly two
//! waypoints), the sweep is, geometrically, an *extrusion* of the
//! profile cross-section along the path direction — and that is
//! exactly what `truck_modeling::builder::tsweep` builds. The CAD-
//! depth pass graduates this case to a genuine [`Solid::Brep`]: the
//! profile is assembled as a planar face in the plane perpendicular to
//! the path direction and `tsweep`'d along it, producing a real closed
//! BRep tube that round-trips through STEP/IGES and composes with
//! downstream BRep booleans. See `try_brep_straight_sweep`.
//!
//! ## Mesh-domain path — curved / multi-segment paths
//!
//! A path with three or more waypoints (a polyline / curve) needs a
//! genuine general path-sweep, which `truck` 0.6 does **not** expose
//! (it ships only `tsweep` — linear — and `rsweep` — rotational). Such
//! a sweep stays mesh-domain: walk the path in N steps, place the
//! profile cross-section perpendicular to the path tangent at each
//! step, and stitch adjacent cross-sections with triangle strips. The
//! output is a [`Solid::Mesh`] and downstream BRep ops on it fail with
//! [`valenx_cad::CadError::MeshBackedSolid`] (same as Loft).

use nalgebra::Vector3;
use truck_modeling::{builder, Point3, Solid as TruckSolid, Vector3 as TruckVec3, Wire};
use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_sketch::extrude::extract_profile_lines;

use crate::feature::SweepParams;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// How many samples to take along the path. Matches Loft's per-ring
/// sample count for ring resolution.
pub const SWEEP_PATH_STEPS: usize = 32;

/// How many samples around the cross-section profile.
pub const SWEEP_PROFILE_SAMPLES: usize = 32;

/// Evaluate a Sweep: sample the profile + path, build cross-section
/// rings at each path step, stitch.
pub(crate) fn evaluate(tree: &FeatureTree, p: &SweepParams) -> Result<Solid, FeatureError> {
    let profile = tree.get_sketch(p.profile_sketch)?;
    let path = tree.get_sketch(p.path_sketch)?;

    let profile_wp = extract_profile_lines(profile, 1e-6)?;
    if profile_wp.len() < 3 {
        return Err(FeatureError::EmptyProfile);
    }
    let path_wp = extract_profile_lines(path, 1e-6)?;
    if path_wp.len() < 2 {
        return Err(FeatureError::BadParameter {
            name: "path_sketch",
            reason: format!(
                "sweep path needs at least 2 waypoints, got {}",
                path_wp.len()
            ),
        });
    }

    // BRep fast path — a straight (2-waypoint) path with no twist is a
    // genuine extrusion. Graduate it to a real `Solid::Brep` via
    // `tsweep`. A twist requires a rotational component the linear
    // `tsweep` cannot express, so a twisted sweep stays mesh-domain.
    if path_wp.len() == 2 && p.twist_angle.abs() < 1e-12 {
        if let Some(brep) = try_brep_straight_sweep(&profile_wp, &path_wp)? {
            return Ok(brep);
        }
    }

    // Resample both to fixed counts so we can stitch cleanly.
    let profile_ring = resample_xy(&profile_wp, SWEEP_PROFILE_SAMPLES, true);
    let path_samples = resample_xy(&path_wp, SWEEP_PATH_STEPS, false);

    // Build a ring at each path sample by placing the profile cross-
    // section in the plane *perpendicular* to the path tangent. The
    // profile sketch's local x maps to the in-XY-plane normal of the
    // path direction; its local y maps to world +Z. This makes the
    // swept result a genuine 3-D tube (an earlier revision kept every
    // ring flat at z = 0, collapsing the sweep into a degenerate
    // planar smear).
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let n_path = path_samples.len();
    let n_prof = profile_ring.len();
    for (i, pt) in path_samples.iter().enumerate() {
        let t_frac = if n_path > 1 {
            i as f64 / (n_path - 1) as f64
        } else {
            0.0
        };
        let twist = p.twist_angle * t_frac;
        let (cos_t, sin_t) = (twist.cos(), twist.sin());
        // Path tangent at this sample (central difference), in XY.
        let tangent = path_tangent_xy(&path_samples, i);
        // In-XY-plane normal of the tangent — the cross-section's
        // "width" axis.
        let (nx, ny) = (-tangent.1, tangent.0);
        for q in &profile_ring {
            // Optionally twist the profile in its own (x, y) plane.
            let (qx, qy) = if p.keep_profile_orientation {
                (q.x, q.y)
            } else {
                (q.x * cos_t - q.y * sin_t, q.x * sin_t + q.y * cos_t)
            };
            // Map profile-local (qx, qy) → world: qx along the in-plane
            // normal, qy along +Z. The cross-section therefore stands
            // perpendicular to the path.
            nodes.push(Vector3::new(pt.x + qx * nx, pt.y + qx * ny, qy));
        }
    }

    // Stitch rings: same pattern as Loft.
    let mut conn: Vec<u32> = Vec::new();
    for a in 0..n_path - 1 {
        let b = a + 1;
        let base_a = a * n_prof;
        let base_b = b * n_prof;
        for k in 0..n_prof {
            let k1 = (k + 1) % n_prof;
            conn.push((base_a + k) as u32);
            conn.push((base_b + k) as u32);
            conn.push((base_a + k1) as u32);
            conn.push((base_a + k1) as u32);
            conn.push((base_b + k) as u32);
            conn.push((base_b + k1) as u32);
        }
    }

    // Caps at first and last path sample.
    let first_centroid = centroid_slice(&nodes[..n_prof]);
    let first_idx = nodes.len() as u32;
    nodes.push(first_centroid);
    for k in 0..n_prof {
        let k1 = (k + 1) % n_prof;
        conn.push(first_idx);
        conn.push(k1 as u32);
        conn.push(k as u32);
    }
    let last_base = (n_path - 1) * n_prof;
    let last_slice = &nodes[last_base..last_base + n_prof].to_vec();
    let last_centroid = centroid_slice(last_slice);
    let last_idx = nodes.len() as u32;
    nodes.push(last_centroid);
    for k in 0..n_prof {
        let k1 = (k + 1) % n_prof;
        conn.push(last_idx);
        conn.push((last_base + k) as u32);
        conn.push((last_base + k1) as u32);
    }

    let mut mesh = Mesh::new("sweep");
    mesh.nodes = nodes;
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = conn;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// Try to build a straight-path sweep as a genuine BRep solid.
///
/// A straight (2-waypoint) path makes the sweep a pure extrusion of
/// the profile cross-section along the path direction — exactly what
/// `truck_modeling::builder::tsweep` produces. The profile cross-
/// section is assembled as a planar face standing **perpendicular** to
/// the path: the profile sketch's local x maps to the in-XY-plane
/// normal of the path direction, its local y maps to world +Z. That
/// face is then `tsweep`'d along the path vector.
///
/// Returns:
/// - `Ok(Some(brep))` — a genuine [`Solid::Brep`] tube.
/// - `Ok(None)` — the inputs do not qualify (degenerate path, or
///   `try_attach_plane` rejected the profile); the caller falls
///   through to the mesh-domain sweep.
/// - `Err(_)` — a hard error both paths would reject.
///
/// `path_wp` must have exactly two waypoints.
fn try_brep_straight_sweep(
    profile_wp: &[(f64, f64)],
    path_wp: &[(f64, f64)],
) -> Result<Option<Solid>, FeatureError> {
    debug_assert_eq!(path_wp.len(), 2);
    let (p0, p1) = (path_wp[0], path_wp[1]);
    let dx = p1.0 - p0.0;
    let dy = p1.1 - p0.1;
    let path_len = (dx * dx + dy * dy).sqrt();
    if path_len < 1e-9 {
        // Degenerate zero-length path — let the mesh path reject it.
        return Ok(None);
    }
    // In-XY-plane normal of the path direction — the cross-section's
    // width axis. Profile-local x maps to this, profile-local y to +Z.
    let (nx, ny) = (-dy / path_len, dx / path_len);

    // Build the profile face at the path's start point, standing
    // perpendicular to the path.
    let verts: Vec<_> = profile_wp
        .iter()
        .map(|&(qx, qy)| builder::vertex(Point3::new(p0.0 + qx * nx, p0.1 + qx * ny, qy)))
        .collect();
    let mut edges = Vec::with_capacity(verts.len());
    for i in 0..verts.len() {
        let next = (i + 1) % verts.len();
        edges.push(builder::line(&verts[i], &verts[next]));
    }
    let wire: Wire = edges.into();
    let face = match builder::try_attach_plane(&[wire]) {
        Ok(f) => f,
        // truck refused the profile face (self-intersecting / non-
        // planar wire) — fall through to the mesh sweep.
        Err(_) => return Ok(None),
    };
    // tsweep the face along the straight path vector → closed BRep.
    let solid: TruckSolid = builder::tsweep(&face, TruckVec3::new(dx, dy, 0.0));
    Ok(Some(Solid::from_truck(solid)))
}

/// Unit tangent of the path polyline at sample `i`, in the XY plane.
///
/// Uses a central difference on the interior samples and a one-sided
/// difference at the endpoints. Returns `(1, 0)` for a degenerate
/// (zero-length) neighbourhood so the caller never divides by zero.
///
/// Shared with [`super::pipe`], which places its cross-section
/// perpendicular to the path the same way.
pub(crate) fn path_tangent_xy(samples: &[Vector3<f64>], i: usize) -> (f64, f64) {
    let n = samples.len();
    if n < 2 {
        return (1.0, 0.0);
    }
    let (a, b) = if i == 0 {
        (samples[0], samples[1])
    } else if i == n - 1 {
        (samples[n - 2], samples[n - 1])
    } else {
        (samples[i - 1], samples[i + 1])
    };
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-12 {
        (1.0, 0.0)
    } else {
        (dx / len, dy / len)
    }
}

/// Resample an open or closed polyline to `n` points by chord-length.
/// `closed` adds a closing segment from last back to first.
pub(crate) fn resample_xy(waypoints: &[(f64, f64)], n: usize, closed: bool) -> Vec<Vector3<f64>> {
    let m = waypoints.len();
    if m == 0 {
        return Vec::new();
    }
    let mut cum: Vec<f64> = vec![0.0];
    for i in 1..m {
        let (px, py) = waypoints[i - 1];
        let (qx, qy) = waypoints[i];
        cum.push(cum.last().unwrap() + ((qx - px).powi(2) + (qy - py).powi(2)).sqrt());
    }
    if closed {
        let (px, py) = waypoints[m - 1];
        let (qx, qy) = waypoints[0];
        cum.push(cum.last().unwrap() + ((qx - px).powi(2) + (qy - py).powi(2)).sqrt());
    }
    let total = *cum.last().unwrap();
    if total < 1e-12 {
        return vec![Vector3::new(waypoints[0].0, waypoints[0].1, 0.0); n];
    }

    let mut out = Vec::with_capacity(n);
    let denom = if closed {
        n as f64
    } else {
        (n - 1).max(1) as f64
    };
    for k in 0..n {
        let target = total * (k as f64) / denom;
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
        let (sx, sy) = waypoints[seg.min(m - 1)];
        let next_idx = if closed {
            (seg + 1) % m
        } else {
            (seg + 1).min(m - 1)
        };
        let (ex, ey) = waypoints[next_idx];
        out.push(Vector3::new(sx + t * (ex - sx), sy + t * (ey - sy), 0.0));
    }
    out
}

fn centroid_slice(ring: &[Vector3<f64>]) -> Vector3<f64> {
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

    fn path_sketch_horizontal() -> valenx_sketch::Sketch {
        // Straight horizontal path along +X
        let mut s = valenx_sketch::Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(5.0, 0.0);
        s.add_line(a, b).unwrap();
        s
    }

    /// An L-shaped 3-waypoint path — a curved/multi-segment path that
    /// stays mesh-domain (no general path-sweep in truck).
    fn path_sketch_l_shape() -> valenx_sketch::Sketch {
        let mut s = valenx_sketch::Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(4.0, 0.0);
        let c = s.add_point(4.0, 4.0);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s
    }

    #[test]
    fn sweep_square_along_straight_path_graduates_to_brep() {
        // A straight-path sweep is a genuine extrusion → a real BRep
        // solid, not a mesh. The CAD-depth pass graduated this case.
        let mut tree = FeatureTree::new();
        let prof = tree.add_sketch(square_sketch(0.5));
        let path = tree.add_sketch(path_sketch_horizontal());
        let params = SweepParams {
            profile_sketch: prof,
            path_sketch: path,
            twist_angle: 0.0,
            keep_profile_orientation: true,
        };
        let solid = evaluate(&tree, &params).expect("sweep succeeds");
        assert!(
            matches!(solid, Solid::Brep(_)),
            "a straight-path sweep should graduate to a BRep solid"
        );
        // A 1×1 square swept 5 units along +X is a 1×1×5 prism →
        // volume 5 (flat-faced, so the measured volume is exact).
        let v = valenx_cad::measure::solid_volume_tol(&solid, 1e-3).unwrap();
        assert!(
            (v - 5.0).abs() < 1e-6,
            "1×1 square swept 5 along X → volume 5, got {v}"
        );
        assert!(
            valenx_cad::measure::is_closed_solid_tol(&solid, 1e-3).unwrap(),
            "the swept BRep must be a valid closed solid"
        );
    }

    #[test]
    fn brep_straight_sweep_composes_with_a_boolean() {
        // The graduated BRep sweep must round-trip through a downstream
        // BRep boolean — the point of graduating it off the mesh path.
        let mut tree = FeatureTree::new();
        let prof = tree.add_sketch(square_sketch(0.5));
        let path = tree.add_sketch(path_sketch_horizontal());
        let params = SweepParams {
            profile_sketch: prof,
            path_sketch: path,
            twist_angle: 0.0,
            keep_profile_orientation: true,
        };
        let swept = evaluate(&tree, &params).expect("sweep succeeds");
        assert!(matches!(swept, Solid::Brep(_)));
        // Intersect with a box that overlaps the tube's interior.
        let box_s = valenx_cad::box_solid(1.0, 2.0, 2.0)
            .unwrap()
            .translated(2.0, -1.0, -1.0)
            .unwrap();
        let inter =
            valenx_cad::intersection(&swept, &box_s).expect("BRep intersection of the swept tube");
        assert!(inter.faces() > 0, "the intersection should have faces");
    }

    #[test]
    fn sweep_along_curved_path_stays_mesh_domain() {
        // A 3-waypoint (L-shaped) path needs a general path-sweep that
        // truck does not expose — it stays mesh-domain.
        let mut tree = FeatureTree::new();
        let prof = tree.add_sketch(square_sketch(0.3));
        let path = tree.add_sketch(path_sketch_l_shape());
        let params = SweepParams {
            profile_sketch: prof,
            path_sketch: path,
            twist_angle: 0.0,
            keep_profile_orientation: true,
        };
        let solid = evaluate(&tree, &params).expect("sweep succeeds");
        assert!(
            matches!(solid, Solid::Mesh(_)),
            "a curved-path sweep stays mesh-domain"
        );
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.1).unwrap();
        assert!(mesh.total_elements() > 0);
    }

    #[test]
    fn twisted_straight_sweep_stays_mesh_domain() {
        // A twist requires a rotational component `tsweep` cannot
        // express, so even a straight path stays mesh-domain when
        // twisted.
        let mut tree = FeatureTree::new();
        let prof = tree.add_sketch(square_sketch(0.5));
        let path = tree.add_sketch(path_sketch_horizontal());
        let params = SweepParams {
            profile_sketch: prof,
            path_sketch: path,
            twist_angle: std::f64::consts::FRAC_PI_2,
            keep_profile_orientation: false,
        };
        let solid = evaluate(&tree, &params).expect("sweep succeeds");
        assert!(
            matches!(solid, Solid::Mesh(_)),
            "a twisted sweep stays mesh-domain (tsweep has no twist)"
        );
    }
}
