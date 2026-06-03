//! Wings-style face extrude — push a face out along its normal
//! by `distance` and stitch the side walls.
//!
//! The original face is removed and replaced by:
//! - a translated copy (the "cap" — same arity, new vertex indices),
//! - one quad per edge of the original face linking the original
//!   vertices to the new cap vertices.

use nalgebra::Vector3;

use crate::error::SubdivError;
use crate::mesh::SubdivMesh;

/// Extrude a single face `face_id` along its area-weighted normal by
/// `distance`. Returns the new mesh.
pub fn extrude_face(
    mesh: &SubdivMesh,
    face_id: usize,
    distance: f64,
) -> Result<SubdivMesh, SubdivError> {
    if face_id >= mesh.faces.len() {
        return Err(SubdivError::IndexOutOfRange {
            kind: "face",
            idx: face_id,
            limit: mesh.faces.len(),
        });
    }
    if !distance.is_finite() {
        return Err(SubdivError::BadParameter {
            name: "distance",
            reason: "must be finite".into(),
        });
    }

    let mut out = mesh.clone();
    let face = out.faces[face_id].clone();
    if face.len() < 3 {
        return Err(SubdivError::Topology(format!(
            "face {face_id} has arity {} (< 3)",
            face.len()
        )));
    }
    let normal = face_normal(&out.vertices, &face);

    // Cap vertices — translated copies of the face vertices.
    let base = out.vertices.len();
    for &v in &face {
        let p = out.vertices[v];
        out.vertices.push(p + normal * distance);
    }

    // Replace the original face with the cap (same winding so the
    // outward normal still points outward after translation).
    let cap: Vec<usize> = (0..face.len()).map(|i| base + i).collect();
    out.faces[face_id] = cap.clone();

    // Side walls — one quad per edge of the original face. Wind the
    // quad so its outward normal points away from the cap's interior.
    let k = face.len();
    for i in 0..k {
        let a_old = face[i];
        let b_old = face[(i + 1) % k];
        let a_new = cap[i];
        let b_new = cap[(i + 1) % k];
        out.faces.push(vec![a_old, b_old, b_new, a_new]);
    }

    Ok(out)
}

/// Area-weighted face normal — Newell's method works for non-planar
/// polygons too.
pub(crate) fn face_normal(verts: &[Vector3<f64>], face: &[usize]) -> Vector3<f64> {
    let mut n: Vector3<f64> = Vector3::zeros();
    let k = face.len();
    for i in 0..k {
        let p = verts[face[i]];
        let q = verts[face[(i + 1) % k]];
        n.x += (p.y - q.y) * (p.z + q.z);
        n.y += (p.z - q.z) * (p.x + q.x);
        n.z += (p.x - q.x) * (p.y + q.y);
    }
    let len = n.norm();
    if len < 1e-12 {
        Vector3::<f64>::z()
    } else {
        n / len
    }
}
