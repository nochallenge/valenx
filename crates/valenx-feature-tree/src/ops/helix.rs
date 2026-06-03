//! Helix evaluator — sweeps a profile along a parametric helix path.
//!
//! Phase 13C Task 21. v1: discrete sampling — generate `samples_per_turn`
//!   `* turns` points along the helix axis, transform the profile to each
//!   sample, stitch the result.
//!
//! Mesh-backed output (same caveat as Loft / Sweep).

use std::f64::consts::TAU;

use nalgebra::Vector3;
use valenx_cad::Solid;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_sketch::extrude::extract_profile_lines;

use crate::feature::HelixParams;
use crate::tree::FeatureTree;
use crate::FeatureError;

/// Cross-section samples per helix turn.
pub const HELIX_SAMPLES_PER_TURN: usize = 24;

/// Cross-section profile samples around the perimeter.
pub const HELIX_PROFILE_SAMPLES: usize = 24;

/// Evaluate a Helix.
pub(crate) fn evaluate(tree: &FeatureTree, p: &HelixParams) -> Result<Solid, FeatureError> {
    if !p.pitch.is_finite() || p.pitch <= 0.0 {
        return Err(FeatureError::BadParameter {
            name: "pitch",
            reason: format!("must be > 0 and finite, got {}", p.pitch),
        });
    }
    if !p.turns.is_finite() || p.turns <= 0.0 {
        return Err(FeatureError::BadParameter {
            name: "turns",
            reason: format!("must be > 0 and finite, got {}", p.turns),
        });
    }
    let axis_len = p.axis_direction.norm();
    if axis_len < 1e-12 {
        return Err(FeatureError::BadParameter {
            name: "axis_direction",
            reason: format!("must have nonzero magnitude, got {:?}", p.axis_direction),
        });
    }
    let axis = p.axis_direction / axis_len;

    let profile = tree.get_sketch(p.profile_sketch)?;
    let profile_wp = extract_profile_lines(profile, 1e-6)?;
    if profile_wp.len() < 3 {
        return Err(FeatureError::EmptyProfile);
    }
    let profile_ring = crate::ops::sweep::resample_xy(&profile_wp, HELIX_PROFILE_SAMPLES, true);

    // Build an orthonormal basis (axis, u, v) for the helix frame.
    let (u, v) = orthonormal_basis(&axis);

    let total_samples = (HELIX_SAMPLES_PER_TURN as f64 * p.turns).ceil() as usize + 1;
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let n_prof = profile_ring.len();
    for i in 0..total_samples {
        let t = i as f64 / (total_samples - 1).max(1) as f64; // 0..=1 along helix
        let theta = if p.left_handed {
            -TAU * p.turns * t
        } else {
            TAU * p.turns * t
        };
        let along = p.pitch * p.turns * t;
        // Taper: radius grows linearly with `along`.
        let taper = p.taper_angle.to_radians().tan();
        // Center of the cross-section at this sample.
        let center = p.axis_origin + axis * along;
        let cos_t = theta.cos();
        let sin_t = theta.sin();
        for q in &profile_ring {
            // Profile (q.x, q.y) lives in the local (u, v) plane and
            // rotates by theta about the axis. Taper scales the radial
            // component by (1 + taper * along).
            let scale = 1.0 + taper * along;
            let radial = (q.x * cos_t - q.y * sin_t) * scale;
            let tangential = (q.x * sin_t + q.y * cos_t) * scale;
            let p_world = center + u * radial + v * tangential;
            nodes.push(p_world);
        }
    }

    // Stitch.
    let mut conn: Vec<u32> = Vec::new();
    for a in 0..total_samples - 1 {
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

    let mut mesh = Mesh::new("helix");
    mesh.nodes = nodes;
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = conn;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Ok(Solid::from_mesh(mesh))
}

/// Build any orthonormal pair (u, v) perpendicular to `axis`.
fn orthonormal_basis(axis: &Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
    let candidate = if axis.x.abs() < 0.9 {
        Vector3::new(1.0, 0.0, 0.0)
    } else {
        Vector3::new(0.0, 1.0, 0.0)
    };
    let u = (candidate - axis * candidate.dot(axis)).normalize();
    let v = axis.cross(&u);
    (u, v)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn circle_sketch_8(radius: f64) -> valenx_sketch::Sketch {
        use std::f64::consts::TAU;
        let mut s = valenx_sketch::Sketch::new();
        let n = 8;
        let mut ids = Vec::new();
        for i in 0..n {
            let a = (i as f64 / n as f64) * TAU;
            ids.push(s.add_point(radius * a.cos(), radius * a.sin()));
        }
        for i in 0..n {
            let j = (i + 1) % n;
            s.add_line(ids[i], ids[j]).unwrap();
        }
        s
    }

    #[test]
    fn helix_2_turns_builds_mesh() {
        let mut tree = FeatureTree::new();
        let prof = tree.add_sketch(circle_sketch_8(0.1));
        let params = HelixParams {
            profile_sketch: prof,
            pitch: 1.0,
            turns: 2.0,
            axis_origin: Vector3::zeros(),
            axis_direction: Vector3::new(0.0, 0.0, 1.0),
            taper_angle: 0.0,
            left_handed: false,
        };
        let solid = evaluate(&tree, &params).expect("helix succeeds");
        let mesh = valenx_cad::solid_to_mesh(&solid, 0.1).unwrap();
        assert!(mesh.total_elements() > 0);
    }

    #[test]
    fn helix_rejects_zero_pitch() {
        let mut tree = FeatureTree::new();
        let prof = tree.add_sketch(circle_sketch_8(0.1));
        let params = HelixParams {
            profile_sketch: prof,
            pitch: 0.0,
            turns: 1.0,
            axis_origin: Vector3::zeros(),
            axis_direction: Vector3::new(0.0, 0.0, 1.0),
            taper_angle: 0.0,
            left_handed: false,
        };
        let err = evaluate(&tree, &params).unwrap_err();
        assert!(matches!(
            err,
            FeatureError::BadParameter { name: "pitch", .. }
        ));
    }
}
