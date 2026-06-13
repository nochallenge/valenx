//! Procedural 3-D mesh of the **Valenx LV-1** launch vehicle.
//!
//! A triangulated surface of revolution — two cylindrical stages + a
//! conical nosecone — plus a bottom cap, a flared engine nozzle, and four
//! fins. Loaded into the central viewport so the designed rocket can be
//! inspected and orbited in 3-D (not just read off a flight chart).
//!
//! Dimensions are illustrative proportions in metres-ish units; the
//! vehicle stands along **+Z** with its base at `z = 0`, and the viewport
//! auto-frames it.

use std::f64::consts::TAU;

use nalgebra::Vector3;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

/// Build a triangulated 3-D surface mesh of the Valenx LV-1.
pub fn lv1_rocket_mesh() -> Mesh {
    let seg = 32usize;
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Axial profile (z, radius), bottom → top: stage 1, a step down to the
    // narrower stage 2, then the nosecone tapering to a point.
    let profile = [
        (0.0_f64, 0.9_f64),
        (14.0, 0.9),
        (14.0, 0.7),
        (20.0, 0.7),
        (26.0, 0.0),
    ];

    // One ring of `seg` vertices per profile point (a single apex vertex
    // when the radius is zero).
    let mut ring_start: Vec<usize> = Vec::with_capacity(profile.len());
    for &(z, r) in &profile {
        ring_start.push(nodes.len());
        if r <= 0.0 {
            nodes.push(Vector3::new(0.0, 0.0, z));
        } else {
            for k in 0..seg {
                let a = k as f64 / seg as f64 * TAU;
                nodes.push(Vector3::new(r * a.cos(), r * a.sin(), z));
            }
        }
    }

    // Stitch consecutive rings into the side surface.
    for w in 0..profile.len() - 1 {
        let s0 = ring_start[w];
        let s1 = ring_start[w + 1];
        if profile[w + 1].1 <= 0.0 {
            // Ring → apex (nosecone tip).
            let apex = s1;
            for k in 0..seg {
                tris.extend_from_slice(&[s0 + k, s0 + (k + 1) % seg, apex]);
            }
        } else {
            // Ring → ring band (cylinder / frustum): two triangles per quad.
            for k in 0..seg {
                let a0 = s0 + k;
                let a1 = s0 + (k + 1) % seg;
                let b0 = s1 + k;
                let b1 = s1 + (k + 1) % seg;
                tris.extend_from_slice(&[a0, a1, b1]);
                tris.extend_from_slice(&[a0, b1, b0]);
            }
        }
    }

    // Bottom cap (disk at z = 0), wound to face downward.
    let base = ring_start[0];
    let base_center = nodes.len();
    nodes.push(Vector3::new(0.0, 0.0, 0.0));
    for k in 0..seg {
        tris.extend_from_slice(&[base_center, base + (k + 1) % seg, base + k]);
    }

    // Engine nozzle — a short flared cone below the base.
    let noz_top = nodes.len();
    for k in 0..seg {
        let a = k as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(0.45 * a.cos(), 0.45 * a.sin(), 0.0));
    }
    let noz_lip = nodes.len();
    for k in 0..seg {
        let a = k as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(0.6 * a.cos(), 0.6 * a.sin(), -2.0));
    }
    for k in 0..seg {
        let a0 = noz_top + k;
        let a1 = noz_top + (k + 1) % seg;
        let b0 = noz_lip + k;
        let b1 = noz_lip + (k + 1) % seg;
        tris.extend_from_slice(&[a0, a1, b1]);
        tris.extend_from_slice(&[a0, b1, b0]);
    }

    // Four fins at the base — flat trapezoids, emitted with both windings
    // so each is visible from either side.
    for f in 0..4 {
        let ang = f as f64 / 4.0 * TAU;
        let (cx, cy) = (ang.cos(), ang.sin());
        let r_body = 0.9;
        let i = nodes.len();
        nodes.push(Vector3::new(cx * r_body, cy * r_body, 0.0)); // root low
        nodes.push(Vector3::new(cx * r_body, cy * r_body, 5.0)); // root high
        nodes.push(Vector3::new(cx * 2.4, cy * 2.4, -0.5)); // tip
        tris.extend_from_slice(&[i, i + 1, i + 2, i, i + 2, i + 1]);
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-lv1-rocket");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    mesh
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rocket_mesh_is_nonempty_and_triangulated() {
        let m = lv1_rocket_mesh();
        assert!(m.nodes.len() > 100, "node count {}", m.nodes.len());
        assert_eq!(m.element_blocks.len(), 1);
        let conn = &m.element_blocks[0].connectivity;
        assert!(!conn.is_empty());
        assert_eq!(conn.len() % 3, 0, "flat triangle connectivity");
        assert!(
            conn.iter().all(|&i| (i as usize) < m.nodes.len()),
            "every index is in range"
        );
        // Stands along +Z: base/nozzle below 0, nose tip near 26.
        let zmin = m.nodes.iter().map(|n| n.z).fold(f64::INFINITY, f64::min);
        let zmax = m
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(zmax > 25.0 && zmin < -1.0, "z span {zmin}..{zmax}");
    }
}
