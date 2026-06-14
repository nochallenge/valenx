//! Procedural 3-D mesh of the **Valenx LV-1** launch vehicle.
//!
//! A triangulated surface of revolution built from a realistic launcher
//! silhouette — a smooth **tangent-ogive payload fairing**, two cylindrical
//! stages joined by an interstage, a boat-tailed engine skirt, a curved
//! **engine bell**, and four swept fins. Loaded into the central viewport so
//! the designed rocket can be inspected and orbited in 3-D.
//!
//! The vehicle stands along **+Z** with its engine plane at `z = 0` and the
//! nose tip near `z = 30`; the viewport auto-frames it. Dimensions are
//! illustrative proportions in metres-ish units.

use std::f64::consts::TAU;

use nalgebra::Vector3;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

/// Revolve an axial `(z, radius)` profile into a surface of revolution,
/// appending vertices + triangles. A profile point with `radius ≤ 0` is an
/// apex (a single vertex), so the band below it is closed with a fan.
fn revolve(
    profile: &[(f64, f64)],
    seg: usize,
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
) {
    let mut ring_start = Vec::with_capacity(profile.len());
    for &(z, r) in profile {
        ring_start.push(nodes.len());
        if r <= 1e-9 {
            nodes.push(Vector3::new(0.0, 0.0, z));
        } else {
            for k in 0..seg {
                let a = k as f64 / seg as f64 * TAU;
                nodes.push(Vector3::new(r * a.cos(), r * a.sin(), z));
            }
        }
    }
    for w in 0..profile.len() - 1 {
        let s0 = ring_start[w];
        let s1 = ring_start[w + 1];
        if profile[w + 1].1 <= 1e-9 {
            let apex = s1;
            for k in 0..seg {
                tris.extend_from_slice(&[s0 + k, s0 + (k + 1) % seg, apex]);
            }
        } else if profile[w].1 <= 1e-9 {
            let apex = s0;
            for k in 0..seg {
                tris.extend_from_slice(&[apex, s1 + (k + 1) % seg, s1 + k]);
            }
        } else {
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
}

/// Build a triangulated 3-D surface mesh of the Valenx LV-1.
pub fn lv1_rocket_mesh() -> Mesh {
    let seg = 48usize;
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let r_body = 2.4_f64;

    // ── Body: engine skirt → stage 1 → interstage → stage 2 → fairing
    // shoulder, then a smooth tangent-ogive nose to the tip. ───────────────
    let mut profile: Vec<(f64, f64)> = vec![
        (0.0, 1.8),     // engine plane (necked for the bell)
        (1.4, r_body),  // boat-tail flare up to the body
        (16.0, r_body), // stage 1
        (16.7, 2.2),    // interstage step-down
        (23.0, 2.2),    // stage 2
        (23.6, r_body), // payload-fairing shoulder
    ];
    // Tangent-ogive fairing from the shoulder (r = r_body) to the tip (r = 0).
    let (nose_z0, nose_r, nose_len) = (23.6_f64, r_body, 8.5_f64);
    let rho = (nose_r * nose_r + nose_len * nose_len) / (2.0 * nose_r);
    let n_nose = 14;
    for i in 1..=n_nose {
        let frac = i as f64 / n_nose as f64; // 0 → shoulder, 1 → tip
        let x_from_tip = nose_len * (1.0 - frac);
        let r = (rho * rho - (nose_len - x_from_tip).powi(2))
            .max(0.0)
            .sqrt()
            - (rho - nose_r);
        profile.push((nose_z0 + nose_len * frac, r.max(0.0)));
    }
    revolve(&profile, seg, &mut nodes, &mut tris);

    // ── Curved engine bell below the engine plane (parabolic flare). ───────
    let bell: Vec<(f64, f64)> = (0..=6)
        .map(|i| {
            let t = i as f64 / 6.0; // 0 at the engine plane, 1 at the exit
            let z = -2.8 * t;
            let r = 1.0 + 1.1 * t * t; // throat-ish 1.0 → exit ~2.1
            (z, r)
        })
        .collect();
    revolve(&bell, seg, &mut nodes, &mut tris);
    // Close the bell exit with a disk so the base reads as a solid engine.
    let exit_z = bell.last().unwrap().0;
    let exit_start = nodes.len() - seg; // last revolved ring = bell exit
    let exit_center = nodes.len();
    nodes.push(Vector3::new(0.0, 0.0, exit_z));
    for k in 0..seg {
        tris.extend_from_slice(&[exit_center, exit_start + k, exit_start + (k + 1) % seg]);
    }

    // ── Four swept fins at the base — trapezoids, double-sided. ────────────
    for f in 0..4 {
        let ang = f as f64 / 4.0 * TAU;
        let (cx, cy) = (ang.cos(), ang.sin());
        let i = nodes.len();
        nodes.push(Vector3::new(cx * r_body, cy * r_body, 1.4)); // root fwd
        nodes.push(Vector3::new(cx * r_body, cy * r_body, 5.5)); // root aft-top
        nodes.push(Vector3::new(cx * 4.8, cy * 4.8, 3.6)); // tip-top (swept out)
        nodes.push(Vector3::new(cx * 4.8, cy * 4.8, -0.2)); // tip-aft (swept down)
                                                            // Two triangles (root-fwd, root-top, tip-top) + (root-fwd, tip-top,
                                                            // tip-aft), emitted with both windings so the flat fin shows on both
                                                            // faces.
        for &(a, b, c) in &[(i, i + 1, i + 2), (i, i + 2, i + 3)] {
            tris.extend_from_slice(&[a, b, c, a, c, b]);
        }
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-lv1-rocket");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    mesh
}

/// A procedurally **detailed liquid-rocket engine** — a combustion chamber, a
/// regeneratively cooled nozzle with a fluted "tube-bundle" surface (the
/// cooling channels that make a real engine look complex), and an injector
/// dome. Stands along +Z, exit at `z = 0`.
///
/// Honest scope: a far more engine-like *procedural* model than a bare bell,
/// but still research-grade geometry — the cooling channels are surface
/// flutes (not a real channel network), and there is no turbopump, injector
/// plate, gimbal or plumbing assembly. A true CAD/generative engine is the
/// documented next step.
pub fn detailed_engine_mesh() -> Mesh {
    let seg = 120usize; // high resolution to resolve the cooling channels
    let n_ch = 40.0_f64; // cooling channels around the circumference
    let amp = 0.07_f64; // channel ridge height
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let (exit_r, throat_r, chamber_r) = (2.1_f64, 0.5_f64, 1.25_f64);
    // Contour (z, base radius): diverging bell (exit → throat) → converging →
    // chamber.
    let mut contour: Vec<(f64, f64)> = Vec::new();
    for i in 0..=10 {
        let t = i as f64 / 10.0; // 0 at exit, 1 at throat
        contour.push((
            3.0 * t,
            throat_r + (exit_r - throat_r) * (1.0 - t).powf(1.6),
        ));
    }
    for i in 1..=5 {
        let t = i as f64 / 5.0;
        contour.push((3.0 + 1.3 * t, throat_r + (chamber_r - throat_r) * t));
    }
    contour.push((8.4, chamber_r));

    // Fluted revolve — the regen cooling channels (radius modulated by angle).
    let flute = |phi: f64| -> f64 {
        let ridge = 0.5 + 0.5 * (n_ch * phi).cos();
        ridge * ridge
    };
    let mut ring_start = Vec::with_capacity(contour.len());
    for &(z, rb) in &contour {
        ring_start.push(nodes.len());
        for k in 0..seg {
            let phi = k as f64 / seg as f64 * TAU;
            let r = rb + amp * flute(phi);
            nodes.push(Vector3::new(r * phi.cos(), r * phi.sin(), z));
        }
    }
    for w in 0..contour.len() - 1 {
        let (s0, s1) = (ring_start[w], ring_start[w + 1]);
        for k in 0..seg {
            let a0 = s0 + k;
            let a1 = s0 + (k + 1) % seg;
            let b0 = s1 + k;
            let b1 = s1 + (k + 1) % seg;
            tris.extend_from_slice(&[a0, a1, b1, a0, b1, b0]);
        }
    }

    // Injector dome — a smooth hemisphere closing the top of the chamber.
    let dome: Vec<(f64, f64)> = (0..=8)
        .map(|i| {
            let ang = (i as f64 / 8.0) * std::f64::consts::FRAC_PI_2;
            (8.4 + 1.5 * ang.sin(), chamber_r * ang.cos())
        })
        .collect();
    revolve(&dome, seg, &mut nodes, &mut tris);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-engine");
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
        assert!(m.nodes.len() > 200, "node count {}", m.nodes.len());
        assert_eq!(m.element_blocks.len(), 1);
        let conn = &m.element_blocks[0].connectivity;
        assert!(!conn.is_empty());
        assert_eq!(conn.len() % 3, 0, "flat triangle connectivity");
        assert!(
            conn.iter().all(|&i| (i as usize) < m.nodes.len()),
            "every index is in range"
        );
        // Stands along +Z: engine bell below 0, nose tip near 30.
        let zmin = m.nodes.iter().map(|n| n.z).fold(f64::INFINITY, f64::min);
        let zmax = m
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(zmax > 29.0 && zmin < -2.0, "z span {zmin}..{zmax}");
    }

    #[test]
    fn nose_tapers_smoothly_to_a_point() {
        // The tangent-ogive nose should reach a near-zero radius at the tip
        // (a smooth point, not a blunt cap).
        let m = lv1_rocket_mesh();
        let zmax = m
            .nodes
            .iter()
            .map(|n| n.z)
            .fold(f64::NEG_INFINITY, f64::max);
        let tip_r = m
            .nodes
            .iter()
            .filter(|n| n.z > zmax - 0.1)
            .map(|n| (n.x * n.x + n.y * n.y).sqrt())
            .fold(f64::INFINITY, f64::min);
        assert!(tip_r < 0.15, "nose tip radius {tip_r}");
    }
}
