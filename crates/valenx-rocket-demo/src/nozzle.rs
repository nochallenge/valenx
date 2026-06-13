//! Parametric bell-nozzle geometry — a Rao parabolic-approximation contour
//! revolved into a watertight triangle [`Mesh`] for export (e.g. STL), so an
//! engine design can be turned into a printable model.
//!
//! Honest scope: a first-order **preliminary-design** geometry. The diverging
//! section is the thrust-optimized parabolic (Rao) *approximation* of the bell
//! — a quadratic Bézier between the throat- and exit-wall angles — revolved as
//! a closed solid of revolution (the engine's wetted bell contour). It is not
//! a method-of-characteristics exact contour, and not a wall-thickness shell.

use nalgebra::Vector3;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

/// A generated bell-nozzle contour plus its key radii and length.
#[derive(Debug, Clone)]
pub struct NozzleGeometry {
    /// Throat radius (m).
    pub throat_radius: f64,
    /// Exit radius (m).
    pub exit_radius: f64,
    /// Axial length of the diverging bell (m).
    pub length: f64,
    /// Contour stations `[axial_x, radius]`, throat (`x = 0`) → exit
    /// (`x = length`).
    pub contour: Vec<[f64; 2]>,
}

/// Build a Rao parabolic-approximation bell contour from the throat area and
/// expansion ratio.
///
/// `bell_fraction` is the bell length as a fraction of the equivalent 15°
/// cone (≈ 0.8 for an "80 % bell"); `theta_n_deg` / `theta_e_deg` are the
/// parabola's initial (just past the throat) and exit wall angles; `stations`
/// is the number of contour points (clamped to ≥ 2). Falls back to a straight
/// cone if the two angles are (near-)equal.
///
/// The exit radius follows from area conservation: `R_e = R_t · √ε`, so the
/// exit area is `ε` times the throat area by construction.
pub fn bell_nozzle_contour(
    throat_area: f64,
    expansion_ratio: f64,
    bell_fraction: f64,
    theta_n_deg: f64,
    theta_e_deg: f64,
    stations: usize,
) -> NozzleGeometry {
    let stations = stations.max(2);
    let r_t = (throat_area.max(0.0) / std::f64::consts::PI).sqrt();
    let r_e = r_t * expansion_ratio.max(1.0).sqrt();
    // Equivalent 15° conical length, scaled to the requested bell fraction.
    let cone_len = (r_e - r_t) / 15.0_f64.to_radians().tan();
    let length = bell_fraction.max(0.05) * cone_len;

    let theta_n = theta_n_deg.to_radians();
    let theta_e = theta_e_deg.to_radians();
    let n = [0.0, r_t];
    let e = [length, r_e];
    // Quadratic-Bézier control point Q is where the throat-angle line from N
    // meets the exit-angle line from E. det = sin(θe − θn) is well-defined
    // (nonzero) whenever the angles differ.
    let det = (theta_e - theta_n).sin();
    let contour: Vec<[f64; 2]> = if det.abs() > 1e-6 {
        let a = (length * theta_e.sin() - theta_e.cos() * (r_e - r_t)) / det;
        let q = [a * theta_n.cos(), r_t + a * theta_n.sin()];
        (0..stations)
            .map(|i| {
                let u = i as f64 / (stations - 1) as f64;
                let w0 = (1.0 - u) * (1.0 - u);
                let w1 = 2.0 * (1.0 - u) * u;
                let w2 = u * u;
                let x = w0 * n[0] + w1 * q[0] + w2 * e[0];
                let r = w0 * n[1] + w1 * q[1] + w2 * e[1];
                [x, r.max(0.0)]
            })
            .collect()
    } else {
        // Degenerate angles → straight cone from throat to exit.
        (0..stations)
            .map(|i| {
                let u = i as f64 / (stations - 1) as f64;
                [u * length, r_t + u * (r_e - r_t)]
            })
            .collect()
    };

    NozzleGeometry {
        throat_radius: r_t,
        exit_radius: r_e,
        length,
        contour,
    }
}

/// Revolve a default 80 %-bell contour for `(throat_area, expansion_ratio)`
/// into a watertight triangle [`Mesh`] (a closed solid of revolution along
/// +Z: throat at `z = 0`, exit at `z = length`, with flat throat and exit
/// caps). `segments` is the angular resolution (clamped to ≥ 3).
pub fn nozzle_mesh(throat_area: f64, expansion_ratio: f64, segments: usize) -> Mesh {
    let geom = bell_nozzle_contour(throat_area, expansion_ratio, 0.8, 22.0, 14.0, 40);
    revolve_contour(&geom, segments.max(3))
}

/// Revolve a [`NozzleGeometry`] contour around the +Z axis into a closed
/// (watertight) Tri3 [`Mesh`].
fn revolve_contour(geom: &NozzleGeometry, segments: usize) -> Mesh {
    let rings = geom.contour.len();
    let mut nodes: Vec<Vector3<f64>> = Vec::with_capacity(rings * segments + 2);
    // Ring vertices: one ring per contour station, `segments` around.
    for &[x, r] in &geom.contour {
        for j in 0..segments {
            let phi = std::f64::consts::TAU * j as f64 / segments as f64;
            nodes.push(Vector3::new(r * phi.cos(), r * phi.sin(), x));
        }
    }
    // Flat-cap centre vertices at the throat and exit planes.
    let throat_center = nodes.len() as u32;
    nodes.push(Vector3::new(0.0, 0.0, geom.contour[0][0]));
    let exit_center = nodes.len() as u32;
    nodes.push(Vector3::new(0.0, 0.0, geom.contour[rings - 1][0]));

    let idx = |i: usize, j: usize| -> u32 { (i * segments + (j % segments)) as u32 };
    let mut tris: Vec<u32> = Vec::with_capacity(6 * rings * segments);

    // Side wall: two triangles per (ring, segment) quad.
    for i in 0..rings - 1 {
        for j in 0..segments {
            let a = idx(i, j);
            let b = idx(i, j + 1);
            let c = idx(i + 1, j);
            let d = idx(i + 1, j + 1);
            tris.extend_from_slice(&[a, c, b, b, c, d]);
        }
    }
    // Throat cap — fan from the throat centre (wound to face −Z).
    for j in 0..segments {
        let a = idx(0, j);
        let b = idx(0, j + 1);
        tris.extend_from_slice(&[throat_center, b, a]);
    }
    // Exit cap — fan from the exit centre (wound to face +Z).
    for j in 0..segments {
        let a = idx(rings - 1, j);
        let b = idx(rings - 1, j + 1);
        tris.extend_from_slice(&[exit_center, a, b]);
    }

    let mut mesh = Mesh::new("nozzle");
    mesh.nodes = nodes;
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris;
    mesh.element_blocks.push(block);
    mesh
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn contour_runs_throat_to_exit_with_right_radii() {
        let g = bell_nozzle_contour(0.05, 16.0, 0.8, 22.0, 14.0, 30);
        let r_t = (0.05_f64 / std::f64::consts::PI).sqrt();
        assert!((g.throat_radius - r_t).abs() < 1e-12);
        // Exit radius = R_t·√ε ⇒ exit area = ε · throat area.
        assert!((g.exit_radius - r_t * 16.0_f64.sqrt()).abs() < 1e-9);
        assert!(g.length > 0.0);
        assert!((g.contour.first().unwrap()[1] - g.throat_radius).abs() < 1e-9);
        assert!((g.contour.last().unwrap()[1] - g.exit_radius).abs() < 1e-9);
        // The bell radius grows monotonically from throat to exit.
        for w in g.contour.windows(2) {
            assert!(w[1][1] >= w[0][1] - 1e-9, "bell radius should not shrink");
            assert!(w[1][0] >= w[0][0] - 1e-9, "axial station should advance");
        }
    }

    #[test]
    fn nozzle_mesh_is_watertight() {
        let m = nozzle_mesh(0.05, 16.0, 24);
        assert!(!m.nodes.is_empty());
        let block = &m.element_blocks[0];
        assert_eq!(block.element_type, ElementType::Tri3);
        let conn = &block.connectivity;
        assert_eq!(conn.len() % 3, 0);
        // Every connectivity index is in range.
        for &i in conn {
            assert!((i as usize) < m.nodes.len(), "index {i} out of range");
        }
        // Closed 2-manifold: every undirected edge is shared by exactly two
        // triangles.
        let mut edges: HashMap<(u32, u32), u32> = HashMap::new();
        for t in conn.chunks_exact(3) {
            for &(a, b) in &[(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *edges.entry(key).or_insert(0) += 1;
            }
        }
        let non_manifold = edges.values().filter(|&&c| c != 2).count();
        assert_eq!(non_manifold, 0, "{non_manifold} non-manifold edges");
    }

    #[test]
    fn mesh_exit_area_matches_expansion_ratio() {
        let throat_area = 0.05;
        let eps = 16.0;
        let m = nozzle_mesh(throat_area, eps, 48);
        // The largest ring radius is the exit radius.
        let max_r = m
            .nodes
            .iter()
            .map(|p| (p.x * p.x + p.y * p.y).sqrt())
            .fold(0.0_f64, f64::max);
        let exit_area = std::f64::consts::PI * max_r * max_r;
        assert!(
            (exit_area / throat_area - eps).abs() / eps < 0.02,
            "exit area ratio {} vs ε {}",
            exit_area / throat_area,
            eps
        );
    }

    #[test]
    fn degenerate_angles_fall_back_to_a_cone() {
        // Equal angles ⇒ no Bézier control point; a straight cone is produced
        // (still throat→exit, still the right radii).
        let g = bell_nozzle_contour(0.05, 16.0, 0.8, 15.0, 15.0, 20);
        assert!((g.contour.last().unwrap()[1] - g.exit_radius).abs() < 1e-9);
        // A cone's radius is exactly linear in the axial fraction.
        let mid = g.contour[g.contour.len() / 2];
        let frac = mid[0] / g.length;
        let expected = g.throat_radius + frac * (g.exit_radius - g.throat_radius);
        assert!((mid[1] - expected).abs() < 1e-9, "cone should be linear");
    }
}
