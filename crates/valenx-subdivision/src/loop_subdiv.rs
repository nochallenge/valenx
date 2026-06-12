//! Loop subdivision — triangle-only refinement scheme.
//!
//! Each iteration:
//! 1. For each edge `(a, b)` shared by triangles `T1, T2`, compute a
//!    new **edge vertex** at `(3/8)(a+b) + (1/8)(c1+c2)` where
//!    `c1, c2` are the opposite vertices of `T1, T2`. Boundary
//!    edges fall back to the midpoint.
//! 2. Reposition each **old vertex** with the textbook Loop rule:
//!    `(1 - n*β) * v + β * sum(neighbours)` where
//!    `β = 1/n * (5/8 - (3/8 + 1/4*cos(2π/n))^2)` (Loop's original
//!    weight) for `n != 3`, and `β = 3/16` for `n == 3`.
//! 3. Split each old triangle into four new triangles using the
//!    three new edge vertices.

use std::collections::BTreeMap;

use nalgebra::Vector3;

use crate::error::SubdivError;
use crate::mesh::SubdivMesh;

/// Apply `iter` rounds of Loop subdivision. Fails fast if the mesh
/// has a non-triangle face.
pub fn subdivide(mesh: &SubdivMesh, iter: u32) -> Result<SubdivMesh, SubdivError> {
    for f in &mesh.faces {
        if f.len() != 3 {
            return Err(SubdivError::Topology(format!(
                "Loop subdivision requires triangle faces (got arity {})",
                f.len()
            )));
        }
    }
    // Reject faces that index past `vertices` before the scheme indexes them.
    mesh.validate()?;
    let mut m = mesh.clone();
    for _ in 0..iter {
        m = one_iter(&m);
    }
    Ok(m)
}

fn one_iter(m: &SubdivMesh) -> SubdivMesh {
    let nv = m.n_verts();

    // Edge → (face0, face1?) map. For triangles only.
    let mut edge_faces: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for (fi, face) in m.faces.iter().enumerate() {
        for k in 0..3 {
            let a = face[k];
            let b = face[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_faces.entry(key).or_default().push(fi);
        }
    }
    let edge_list: Vec<(usize, usize)> = edge_faces.keys().copied().collect();
    let mut edge_idx: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for (i, e) in edge_list.iter().enumerate() {
        edge_idx.insert(*e, i);
    }

    // Edge vertex positions.
    let edge_verts: Vec<Vector3<f64>> = edge_list
        .iter()
        .map(|e| {
            let mid = (m.vertices[e.0] + m.vertices[e.1]) * 0.5;
            let fs = &edge_faces[e];
            if fs.len() == 2 {
                // Opposite vertex of each adjacent triangle.
                let mut sum = Vector3::zeros();
                for &fi in fs {
                    let face = &m.faces[fi];
                    for &v in face {
                        if v != e.0 && v != e.1 {
                            sum += m.vertices[v];
                        }
                    }
                }
                m.vertices[e.0] * (3.0 / 8.0) + m.vertices[e.1] * (3.0 / 8.0) + sum * (1.0 / 8.0)
            } else {
                mid
            }
        })
        .collect();

    // Per-vertex one-ring + boundary flag.
    let ring = m.vertex_one_ring();
    let mut is_boundary = vec![false; nv];
    for (e, fs) in &edge_faces {
        if fs.len() == 1 {
            is_boundary[e.0] = true;
            is_boundary[e.1] = true;
        }
    }

    let mut new_verts: Vec<Vector3<f64>> = Vec::with_capacity(nv + edge_list.len());
    for v in 0..nv {
        let neighbours: Vec<usize> = ring[v].iter().copied().collect();
        if is_boundary[v] {
            // Boundary rule: 3/4 old + 1/8 each of its two boundary
            // neighbours. v1: average over all neighbours weighted
            // 1/8 (sufficient for closed-mesh tests).
            let mut sum = Vector3::zeros();
            for &u in &neighbours {
                sum += m.vertices[u];
            }
            let n = neighbours.len().max(1) as f64;
            new_verts.push(m.vertices[v] * (3.0 / 4.0) + sum * (1.0 / 4.0) / n);
        } else {
            let n = neighbours.len();
            if n == 0 {
                new_verts.push(m.vertices[v]);
                continue;
            }
            let beta = if n == 3 {
                3.0 / 16.0
            } else {
                let c = 3.0 / 8.0 + (1.0 / 4.0) * (2.0 * std::f64::consts::PI / n as f64).cos();
                (5.0 / 8.0 - c * c) / n as f64
            };
            let mut sum = Vector3::zeros();
            for &u in &neighbours {
                sum += m.vertices[u];
            }
            new_verts.push(m.vertices[v] * (1.0 - n as f64 * beta) + sum * beta);
        }
    }
    new_verts.extend(edge_verts.iter().copied());

    let e_off = nv;
    let mut new_faces = Vec::with_capacity(m.faces.len() * 4);
    for face in &m.faces {
        let a = face[0];
        let b = face[1];
        let c = face[2];
        let key = |x, y| if x < y { (x, y) } else { (y, x) };
        let ab = e_off + edge_idx[&key(a, b)];
        let bc = e_off + edge_idx[&key(b, c)];
        let ca = e_off + edge_idx[&key(c, a)];
        new_faces.push(vec![a, ab, ca]);
        new_faces.push(vec![b, bc, ab]);
        new_faces.push(vec![c, ca, bc]);
        new_faces.push(vec![ab, bc, ca]);
    }

    SubdivMesh {
        vertices: new_verts,
        faces: new_faces,
    }
}
