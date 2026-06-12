//! Inset — for each face in `face_ids`, generate a smaller copy
//! pulled `distance` toward the face centroid, then stitch the
//! resulting ring of quads.

use crate::error::BlenderOpError;
use crate::mesh::{face_centroid, Mesh};

/// Inset the given faces by `distance` (in mesh units toward each
/// face's centroid).
pub fn faces(mesh: &Mesh, face_ids: &[usize], distance: f64) -> Result<Mesh, BlenderOpError> {
    if !distance.is_finite() || distance < 0.0 {
        return Err(BlenderOpError::BadParameter {
            name: "distance",
            reason: format!("must be finite and >= 0 (got {distance})"),
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
    for &fi in face_ids {
        let face = out.faces[fi].clone();
        let c = face_centroid(&out.vertices, &face);
        let base = out.vertices.len();
        // Create the inset ring vertices.
        for &v in &face {
            let p = out.vertices[v];
            let dir = c - p;
            let len = dir.norm();
            let off = if len > 1e-12 {
                dir / len * distance
            } else {
                dir
            };
            out.vertices.push(p + off);
        }
        // Inset face — replaces the original.
        let inset_ring: Vec<usize> = (0..face.len()).map(|i| base + i).collect();
        out.faces[fi] = inset_ring.clone();
        // Ring of quads between original and inset.
        let k = face.len();
        for i in 0..k {
            let a_old = face[i];
            let b_old = face[(i + 1) % k];
            let b_new = inset_ring[(i + 1) % k];
            let a_new = inset_ring[i];
            out.faces.push(vec![a_old, b_old, b_new, a_new]);
        }
    }
    Ok(out)
}
