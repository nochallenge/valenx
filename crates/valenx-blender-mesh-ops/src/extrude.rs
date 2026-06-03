//! Region extrude — push a set of faces along a translation vector
//! and stitch the resulting walls.

use nalgebra::Vector3;

use crate::error::BlenderOpError;
use crate::mesh::Mesh;

/// Extrude the faces in `face_ids` by `vector`. Each input face is
/// replaced by a translated copy; new side-wall quads stitch the
/// boundary edges of the region to their translated counterparts.
pub fn region(
    mesh: &Mesh,
    face_ids: &[usize],
    vector: Vector3<f64>,
) -> Result<Mesh, BlenderOpError> {
    if face_ids.is_empty() {
        return Err(BlenderOpError::BadParameter {
            name: "face_ids",
            reason: "must not be empty".into(),
        });
    }
    for &fi in face_ids {
        if fi >= mesh.faces.len() {
            return Err(BlenderOpError::IndexOutOfRange {
                kind: "face",
                idx: fi,
                limit: mesh.faces.len(),
            });
        }
    }

    let mut out = mesh.clone();
    let mut vertex_map: std::collections::BTreeMap<usize, usize> = Default::default();

    // Per-region edge multiplicity — boundary edges appear once,
    // interior edges twice. Only the boundary needs walls.
    let mut edge_count: std::collections::BTreeMap<(usize, usize), usize> = Default::default();
    for &fi in face_ids {
        let f = &mesh.faces[fi];
        for k in 0..f.len() {
            let a = f[k];
            let b = f[(k + 1) % f.len()];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }

    // Translate every distinct vertex of the region exactly once.
    let mut all_verts = std::collections::BTreeSet::<usize>::new();
    for &fi in face_ids {
        for &v in &mesh.faces[fi] {
            all_verts.insert(v);
        }
    }
    for v in all_verts {
        let p = mesh.vertices[v];
        let new_idx = out.vertices.len();
        out.vertices.push(p + vector);
        vertex_map.insert(v, new_idx);
    }

    // Replace each input face with its translated cap (keep winding).
    for &fi in face_ids {
        let cap: Vec<usize> = mesh.faces[fi].iter().map(|v| vertex_map[v]).collect();
        out.faces[fi] = cap;
    }

    // For every BOUNDARY edge, add a side-wall quad linking the
    // original edge to its translated counterpart. Preserve winding
    // (the face's exterior side stays exterior).
    for (e, count) in &edge_count {
        if *count != 1 {
            continue;
        }
        let a_old = e.0;
        let b_old = e.1;
        let a_new = vertex_map[&a_old];
        let b_new = vertex_map[&b_old];
        out.faces.push(vec![a_old, b_old, b_new, a_new]);
    }

    Ok(out)
}
