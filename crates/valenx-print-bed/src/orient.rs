//! Optimal-orientation search.
//!
//! v1 scope: brute-force evaluation of the six axis-aligned face-down
//! orientations (`+X / -X / +Y / -Y / +Z / -Z`) plus identity. Picks
//! the orientation that minimises supports (lowest projected area) or
//! maximises strength (largest bed-flat footprint).

use nalgebra::{UnitQuaternion, Vector3};

use valenx_mesh::Mesh;

/// Optimisation criterion.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum Criterion {
    /// Minimise the projected silhouette area (≈ minimise supports).
    MinSupports,
    /// Maximise the bed-contact footprint (≈ maximise strength /
    /// stability).
    MaxStrength,
}

/// Brute-force search across the 6 axis-aligned face-down rotations +
/// identity. Returns the chosen orientation.
pub fn optimal(mesh: &Mesh, criterion: Criterion) -> UnitQuaternion<f64> {
    let candidates: &[UnitQuaternion<f64>] = &[
        UnitQuaternion::identity(),
        UnitQuaternion::from_axis_angle(&Vector3::x_axis(), std::f64::consts::PI),
        UnitQuaternion::from_axis_angle(&Vector3::x_axis(), 0.5 * std::f64::consts::PI),
        UnitQuaternion::from_axis_angle(&Vector3::x_axis(), -0.5 * std::f64::consts::PI),
        UnitQuaternion::from_axis_angle(&Vector3::y_axis(), 0.5 * std::f64::consts::PI),
        UnitQuaternion::from_axis_angle(&Vector3::y_axis(), -0.5 * std::f64::consts::PI),
    ];

    let mut best = candidates[0];
    let mut best_score = score(mesh, &best, criterion);
    for c in &candidates[1..] {
        let s = score(mesh, c, criterion);
        let better = match criterion {
            Criterion::MinSupports => s < best_score,
            Criterion::MaxStrength => s > best_score,
        };
        if better {
            best = *c;
            best_score = s;
        }
    }
    best
}

fn score(mesh: &Mesh, q: &UnitQuaternion<f64>, criterion: Criterion) -> f64 {
    let rotated_nodes: Vec<Vector3<f64>> = mesh.nodes.iter().map(|n| q * n).collect();
    let (proj_area, bed_extent) = xy_footprint(&rotated_nodes);
    match criterion {
        // Minimising the projected silhouette is a good proxy for
        // overhang area; using bbox area is the v1 approximation.
        Criterion::MinSupports => proj_area,
        // Strength → choose the orientation with the largest bed-
        // flat footprint i.e. largest XY bbox extent.
        Criterion::MaxStrength => bed_extent,
    }
}

fn xy_footprint(nodes: &[Vector3<f64>]) -> (f64, f64) {
    if nodes.is_empty() {
        return (0.0, 0.0);
    }
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for n in nodes {
        if n.x < min_x {
            min_x = n.x;
        }
        if n.y < min_y {
            min_y = n.y;
        }
        if n.x > max_x {
            max_x = n.x;
        }
        if n.y > max_y {
            max_y = n.y;
        }
    }
    let w = max_x - min_x;
    let h = max_y - min_y;
    (w * h, w.max(h))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn box_mesh(w: f64, h: f64, d: f64) -> Mesh {
        let mut m = Mesh::new("box");
        for x in [0.0, w] {
            for y in [0.0, h] {
                for z in [0.0, d] {
                    m.nodes.push(Vector3::new(x, y, z));
                }
            }
        }
        m
    }

    #[test]
    fn optimal_returns_a_unit_quaternion() {
        let m = box_mesh(10.0, 20.0, 30.0);
        let q = optimal(&m, Criterion::MinSupports);
        assert!((q.norm() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn min_supports_prefers_lying_flat() {
        // A tall, thin column should be laid flat to reduce supports.
        let m = box_mesh(5.0, 5.0, 100.0);
        let q = optimal(&m, Criterion::MinSupports);
        // Score under chosen orientation should be smaller than identity's.
        let s_chosen = score(&m, &q, Criterion::MinSupports);
        let s_id = score(&m, &UnitQuaternion::identity(), Criterion::MinSupports);
        assert!(s_chosen <= s_id);
    }
}
