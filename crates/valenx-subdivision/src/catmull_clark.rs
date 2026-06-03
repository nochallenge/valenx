//! Catmull-Clark subdivision — quad-output scheme that works on
//! any polygon mesh.
//!
//! Each iteration:
//! 1. Compute a **face point** at the centroid of each face.
//! 2. Compute an **edge point** at `((v0 + v1 + f0 + f1) / 4)` where
//!    `f0, f1` are the face points of the two faces sharing the edge.
//!    Boundary edges (one adjacent face) use the midpoint of the
//!    edge plus the single face point, averaged.
//! 3. **Reposition each old vertex** with the textbook formula
//!    `(F + 2R + (n-3)P) / n` where `F` = average of adjacent face
//!    points, `R` = average of adjacent edge midpoints, `P` = old
//!    position, `n` = valence.
//! 4. Build the new face list: for every original face of arity `k`,
//!    emit `k` new quads `[edge_in, face_point, edge_out, vertex]`.

use std::collections::BTreeMap;

use nalgebra::Vector3;

use crate::mesh::SubdivMesh;

/// Apply `iter` rounds of Catmull-Clark subdivision.
pub fn subdivide(mesh: &SubdivMesh, iter: u32) -> SubdivMesh {
    let mut m = mesh.clone();
    for _ in 0..iter {
        m = one_iter(&m);
    }
    m
}

fn one_iter(m: &SubdivMesh) -> SubdivMesh {
    let nv = m.n_verts();
    let nf = m.n_faces();

    // Face points = centroid of each face.
    let face_points: Vec<Vector3<f64>> = m
        .faces
        .iter()
        .map(|f| {
            let s: Vector3<f64> = f.iter().map(|&i| m.vertices[i]).sum();
            s / f.len() as f64
        })
        .collect();

    // Build edge → (face0, face1?) map.
    let mut edge_faces: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for (fi, face) in m.faces.iter().enumerate() {
        for k in 0..face.len() {
            let a = face[k];
            let b = face[(k + 1) % face.len()];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_faces.entry(key).or_default().push(fi);
        }
    }

    // Edge points + edge index map.
    let edge_list: Vec<(usize, usize)> = edge_faces.keys().copied().collect();
    let mut edge_idx: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for (i, e) in edge_list.iter().enumerate() {
        edge_idx.insert(*e, i);
    }
    let edge_points: Vec<Vector3<f64>> = edge_list
        .iter()
        .map(|e| {
            let mid = (m.vertices[e.0] + m.vertices[e.1]) * 0.5;
            let fs = &edge_faces[e];
            if fs.len() == 2 {
                let fp_avg = (face_points[fs[0]] + face_points[fs[1]]) * 0.5;
                (mid + fp_avg) * 0.5
            } else {
                // Boundary edge — Catmull-Clark falls back to edge
                // midpoint (the canonical rule for an open patch).
                mid
            }
        })
        .collect();

    // Repositioned old vertices.
    let mut adj_faces: Vec<Vec<usize>> = vec![Vec::new(); nv];
    for (fi, face) in m.faces.iter().enumerate() {
        for &v in face {
            adj_faces[v].push(fi);
        }
    }
    let mut adj_edges: Vec<Vec<usize>> = vec![Vec::new(); nv];
    for (ei, e) in edge_list.iter().enumerate() {
        adj_edges[e.0].push(ei);
        adj_edges[e.1].push(ei);
    }
    let mut new_old_verts: Vec<Vector3<f64>> = Vec::with_capacity(nv);
    for v in 0..nv {
        let n = adj_faces[v].len().max(1) as f64;
        let f_avg: Vector3<f64> = adj_faces[v]
            .iter()
            .map(|&fi| face_points[fi])
            .sum::<Vector3<f64>>()
            / n;
        let r_avg: Vector3<f64> = adj_edges[v]
            .iter()
            .map(|&ei| {
                let e = edge_list[ei];
                (m.vertices[e.0] + m.vertices[e.1]) * 0.5
            })
            .sum::<Vector3<f64>>()
            / adj_edges[v].len().max(1) as f64;
        let p = m.vertices[v];
        new_old_verts.push((f_avg + r_avg * 2.0 + p * (n - 3.0)) / n);
    }

    // Assemble new vertex array: old verts, then face points, then
    // edge points. Index offsets:
    let f_off = nv;
    let e_off = nv + nf;

    let mut new_verts = Vec::with_capacity(nv + nf + edge_list.len());
    new_verts.extend(new_old_verts);
    new_verts.extend(face_points.iter().copied());
    new_verts.extend(edge_points.iter().copied());

    // Build new faces — each old face of arity k → k quads.
    let mut new_faces = Vec::with_capacity(m.faces.iter().map(|f| f.len()).sum());
    for (fi, face) in m.faces.iter().enumerate() {
        let k = face.len();
        for j in 0..k {
            let v_curr = face[j];
            let v_prev = face[(j + k - 1) % k];
            let v_next = face[(j + 1) % k];
            let e_in = {
                let key = if v_prev < v_curr {
                    (v_prev, v_curr)
                } else {
                    (v_curr, v_prev)
                };
                e_off + edge_idx[&key]
            };
            let e_out = {
                let key = if v_curr < v_next {
                    (v_curr, v_next)
                } else {
                    (v_next, v_curr)
                };
                e_off + edge_idx[&key]
            };
            let fp = f_off + fi;
            new_faces.push(vec![v_curr, e_out, fp, e_in]);
        }
    }

    SubdivMesh {
        vertices: new_verts,
        faces: new_faces,
    }
}
