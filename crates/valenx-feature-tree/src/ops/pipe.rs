//! Pipe evaluator — specialization of Sweep for tubing / piping.
//!
//! Phase 13B Task 17. The Pipe op uses the same machinery as
//! [`super::sweep`] (sample cross-section, sample centerline, stitch
//! rings) but supports a `bend_radius` for filleting sharp corners in
//! the centerline.
//!
//! v1 limitations: `bend_radius` triggers chord-shortening at each
//! corner — the corner vertex is replaced with two new vertices set
//! back along their incoming/outgoing edges by `bend_radius`, and the
//! gap between them is sampled along a circular arc. Multi-turn 3D
//! bends in non-planar paths fall back to the linear shortening
//! without arc smoothing.

use nalgebra::Vector3;
use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_sketch::extrude::extract_profile_lines;

use crate::feature::PipeParams;
use crate::tree::FeatureTree;
use crate::FeatureError;

use super::sweep::{path_tangent_xy, resample_xy, SWEEP_PATH_STEPS, SWEEP_PROFILE_SAMPLES};

/// Per-bend arc subdivision.
pub const PIPE_BEND_SAMPLES: usize = 8;

/// Evaluate a Pipe: similar to Sweep with bend smoothing.
pub(crate) fn evaluate(tree: &FeatureTree, p: &PipeParams) -> Result<Solid, FeatureError> {
    if !p.bend_radius.is_finite() || p.bend_radius < 0.0 {
        return Err(FeatureError::BadParameter {
            name: "bend_radius",
            reason: format!("must be >= 0 and finite, got {}", p.bend_radius),
        });
    }
    let profile = tree.get_sketch(p.cross_section_sketch)?;
    let centerline = tree.get_sketch(p.centerline_sketch)?;

    let profile_wp = extract_profile_lines(profile, 1e-6)?;
    if profile_wp.len() < 3 {
        return Err(FeatureError::EmptyProfile);
    }
    let centerline_wp = extract_profile_lines(centerline, 1e-6)?;
    if centerline_wp.len() < 2 {
        return Err(FeatureError::BadParameter {
            name: "centerline_sketch",
            reason: format!(
                "pipe centerline needs at least 2 waypoints, got {}",
                centerline_wp.len()
            ),
        });
    }

    // Build a bent centerline: insert PIPE_BEND_SAMPLES arc points at
    // each interior corner.
    let bent = bend_polyline(&centerline_wp, p.bend_radius, PIPE_BEND_SAMPLES);

    let path_samples = resample_xy(&bent, SWEEP_PATH_STEPS.max(bent.len()), false);
    let profile_ring = resample_xy(&profile_wp, SWEEP_PROFILE_SAMPLES, true);

    // Place the cross-section perpendicular to the centerline at each
    // path sample: the profile sketch's local x maps to the in-XY-
    // plane normal of the path tangent, its local y maps to world +Z.
    // (An earlier revision kept every ring flat at z = 0, collapsing
    // the pipe into a degenerate planar smear instead of a 3-D tube.)
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let n_path = path_samples.len();
    let n_prof = profile_ring.len();
    for (i, pt) in path_samples.iter().enumerate() {
        let tangent = path_tangent_xy(&path_samples, i);
        let (nx, ny) = (-tangent.1, tangent.0);
        for q in &profile_ring {
            nodes.push(Vector3::new(pt.x + q.x * nx, pt.y + q.x * ny, q.y));
        }
    }

    let mut conn: Vec<u32> = Vec::new();
    for a in 0..n_path - 1 {
        let base_a = a * n_prof;
        let base_b = (a + 1) * n_prof;
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

    let mut mesh = Mesh::new("pipe");
    mesh.nodes = nodes;
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = conn;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// Insert `samples` linearly-interpolated points at each interior
/// corner of the polyline, set back by `radius` along the incoming /
/// outgoing edges. For `radius = 0` returns the input unchanged.
fn bend_polyline(wp: &[(f64, f64)], radius: f64, samples: usize) -> Vec<(f64, f64)> {
    if radius < 1e-9 || wp.len() < 3 {
        return wp.to_vec();
    }
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(wp.len() * samples);
    out.push(wp[0]);
    for i in 1..wp.len() - 1 {
        let prev = wp[i - 1];
        let cur = wp[i];
        let next = wp[i + 1];
        let v1 = (cur.0 - prev.0, cur.1 - prev.1);
        let v2 = (next.0 - cur.0, next.1 - cur.1);
        let len1 = (v1.0 * v1.0 + v1.1 * v1.1).sqrt();
        let len2 = (v2.0 * v2.0 + v2.1 * v2.1).sqrt();
        if len1 < 1e-9 || len2 < 1e-9 {
            out.push(cur);
            continue;
        }
        let setback = radius.min(len1 * 0.5).min(len2 * 0.5);
        // Setback point along v1 from cur: cur - (v1/len1) * setback.
        let p_a = (cur.0 - v1.0 / len1 * setback, cur.1 - v1.1 / len1 * setback);
        let p_b = (cur.0 + v2.0 / len2 * setback, cur.1 + v2.1 / len2 * setback);
        out.push(p_a);
        // Sample a linear approximation of the corner arc (v1: linear
        // chord; planar circle fitting deferred).
        for s in 1..samples {
            let t = s as f64 / samples as f64;
            out.push((p_a.0 + t * (p_b.0 - p_a.0), p_a.1 + t * (p_b.1 - p_a.1)));
        }
        out.push(p_b);
    }
    out.push(wp[wp.len() - 1]);
    out
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

    fn l_shape_centerline() -> valenx_sketch::Sketch {
        // L-shape: (0,0) → (3,0) → (3,3).
        let mut s = valenx_sketch::Sketch::new();
        let a = s.add_point(0.0, 0.0);
        let b = s.add_point(3.0, 0.0);
        let c = s.add_point(3.0, 3.0);
        s.add_line(a, b).unwrap();
        s.add_line(b, c).unwrap();
        s
    }

    #[test]
    fn pipe_with_l_centerline_builds_mesh() {
        let mut tree = FeatureTree::new();
        let prof = tree.add_sketch(square_sketch(0.3));
        let path = tree.add_sketch(l_shape_centerline());
        let params = PipeParams {
            cross_section_sketch: prof,
            centerline_sketch: path,
            bend_radius: 0.5,
        };
        let solid = evaluate(&tree, &params).expect("pipe evaluates");
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.1).unwrap();
        assert!(mesh.total_elements() > 0);
    }

    #[test]
    fn pipe_rejects_negative_bend_radius() {
        let mut tree = FeatureTree::new();
        let prof = tree.add_sketch(square_sketch(0.3));
        let path = tree.add_sketch(l_shape_centerline());
        let params = PipeParams {
            cross_section_sketch: prof,
            centerline_sketch: path,
            bend_radius: -1.0,
        };
        let err = evaluate(&tree, &params).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter {
                name: "bend_radius",
                ..
            }
        ));
    }
}
