//! Catmull-Clark with semi-sharp creases.
//!
//! For edges with `sharpness > 0`:
//! - the edge point becomes the **midpoint** (sharp rule),
//! - adjacent vertices use a tangent-line interpolation: average of
//!   the two crease-edge endpoints adjacent to the vertex, weighted
//!   by 6/8 with the vertex weighted 1/8 each side (the textbook
//!   sharp-vertex rule).
//!
//! Each iteration decrements all sharpness values by 1 (clamped at 0).

use std::collections::BTreeMap;

use nalgebra::Vector3;

use crate::crease::CreasedMesh;
use crate::mesh::SubdivMesh;

/// Apply `iter` rounds of crease-aware Catmull-Clark.
pub fn subdivide(creased: &CreasedMesh, iter: u32) -> SubdivMesh {
    // Guard against a malformed base mesh (a face index past `vertices`); this
    // entry is infallible, so return the base mesh unchanged rather than panic.
    if creased.mesh.validate().is_err() {
        return creased.mesh.clone();
    }
    let mut cur = creased.clone();
    for _ in 0..iter {
        cur = one_iter(&cur);
    }
    cur.mesh
}

fn one_iter(c: &CreasedMesh) -> CreasedMesh {
    let m = &c.mesh;
    let nv = m.n_verts();
    let nf = m.n_faces();

    // Face points.
    let face_points: Vec<Vector3<f64>> = m
        .faces
        .iter()
        .map(|f| {
            let s: Vector3<f64> = f.iter().map(|&i| m.vertices[i]).sum();
            s / f.len() as f64
        })
        .collect();

    // Edge → adjacent faces.
    let mut edge_faces: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    for (fi, face) in m.faces.iter().enumerate() {
        for k in 0..face.len() {
            let a = face[k];
            let b = face[(k + 1) % face.len()];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_faces.entry(key).or_default().push(fi);
        }
    }
    let edge_list: Vec<(usize, usize)> = edge_faces.keys().copied().collect();
    let mut edge_idx: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for (i, e) in edge_list.iter().enumerate() {
        edge_idx.insert(*e, i);
    }

    // Edge points — sharp rule for sharpness > 0, else smooth.
    let edge_points: Vec<Vector3<f64>> = edge_list
        .iter()
        .map(|e| {
            let mid = (m.vertices[e.0] + m.vertices[e.1]) * 0.5;
            let sharp = *c.creases.get(e).unwrap_or(&0.0);
            if sharp > 0.0 {
                mid
            } else {
                let fs = &edge_faces[e];
                if fs.len() == 2 {
                    let fp_avg = (face_points[fs[0]] + face_points[fs[1]]) * 0.5;
                    (mid + fp_avg) * 0.5
                } else {
                    mid
                }
            }
        })
        .collect();

    // Per-vertex crease tally — number of sharp edges incident.
    let mut sharp_edges_at: Vec<Vec<(usize, usize)>> = vec![Vec::new(); nv];
    for e in &edge_list {
        let sharp = *c.creases.get(e).unwrap_or(&0.0);
        if sharp > 0.0 {
            sharp_edges_at[e.0].push(*e);
            sharp_edges_at[e.1].push(*e);
        }
    }

    // Adjacency.
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

    // Vertex repositioning.
    let mut new_old_verts: Vec<Vector3<f64>> = Vec::with_capacity(nv);
    for v in 0..nv {
        let sharp_n = sharp_edges_at[v].len();
        match sharp_n {
            0 => {
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
            1 => {
                // Sharp endpoint — treat as smooth (the lone sharp
                // edge can't form a curve).
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
            2 => {
                // Crease vertex — (1/8) * neighbour_a + (1/8) * neighbour_b
                // + (3/4) * v, where neighbours are the far ends of the
                // two sharp edges incident here.
                let e1 = sharp_edges_at[v][0];
                let e2 = sharp_edges_at[v][1];
                let other = |e: (usize, usize)| if e.0 == v { e.1 } else { e.0 };
                let p = m.vertices[v];
                let a = m.vertices[other(e1)];
                let b = m.vertices[other(e2)];
                new_old_verts.push(p * (6.0 / 8.0) + a * (1.0 / 8.0) + b * (1.0 / 8.0));
            }
            _ => {
                // 3+ sharp edges meeting at a vertex — corner; keep fixed.
                new_old_verts.push(m.vertices[v]);
            }
        }
    }

    // Build new vertex array.
    let f_off = nv;
    let e_off = nv + nf;
    let mut new_verts = Vec::with_capacity(nv + nf + edge_list.len());
    new_verts.extend(new_old_verts);
    new_verts.extend(face_points.iter().copied());
    new_verts.extend(edge_points.iter().copied());

    // Build new faces (quads) the same way as smooth Catmull-Clark.
    let mut new_faces = Vec::with_capacity(m.faces.iter().map(|f| f.len()).sum());
    for (fi, face) in m.faces.iter().enumerate() {
        let k = face.len();
        for j in 0..k {
            let v_curr = face[j];
            let v_prev = face[(j + k - 1) % k];
            let v_next = face[(j + 1) % k];
            let key = |a, b| if a < b { (a, b) } else { (b, a) };
            let e_in = e_off + edge_idx[&key(v_prev, v_curr)];
            let e_out = e_off + edge_idx[&key(v_curr, v_next)];
            let fp = f_off + fi;
            new_faces.push(vec![v_curr, e_out, fp, e_in]);
        }
    }

    // Build new creases — map old (a,b) to two new sub-edges
    // (a, edge_point) and (edge_point, b), sharpness decremented.
    let mut new_creases: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    for e in &edge_list {
        let sharp = *c.creases.get(e).unwrap_or(&0.0);
        if sharp > 0.0 {
            let new_sharp = (sharp - 1.0).max(0.0);
            if new_sharp > 0.0 {
                let ep_idx = e_off + edge_idx[e];
                let k1 = if e.0 < ep_idx { (e.0, ep_idx) } else { (ep_idx, e.0) };
                let k2 = if ep_idx < e.1 { (ep_idx, e.1) } else { (e.1, ep_idx) };
                new_creases.insert(k1, new_sharp);
                new_creases.insert(k2, new_sharp);
            }
        }
    }

    CreasedMesh {
        mesh: SubdivMesh {
            vertices: new_verts,
            faces: new_faces,
        },
        creases: new_creases,
    }
}
