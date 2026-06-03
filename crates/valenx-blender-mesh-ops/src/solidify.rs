//! Solidify — give a single-sided surface a thickness by extruding
//! every face outward and stitching the boundary.
//!
//! This is the Blender "Solidify" modifier; functionally identical
//! to `valenx-mesh::repair::shell` but exposed under the Blender
//! naming so the workbench mirrors the upstream UI.

use nalgebra::Vector3;

use crate::error::BlenderOpError;
use crate::mesh::{face_normal, Mesh};

/// Shell the mesh by `thickness`. Every face is duplicated on the
/// inside (offset by `-thickness * normal`), the original faces
/// stay on the outside, and boundary edges are stitched.
pub fn shell(mesh: &Mesh, thickness: f64) -> Result<Mesh, BlenderOpError> {
    if !thickness.is_finite() || thickness <= 0.0 {
        return Err(BlenderOpError::BadParameter {
            name: "thickness",
            reason: format!("must be finite and > 0 (got {thickness})"),
        });
    }
    let nv = mesh.vertices.len();
    // Inner vertices = outer minus vertex_normal * thickness.
    let mut vnorm = vec![Vector3::<f64>::zeros(); nv];
    for face in &mesh.faces {
        let n = face_normal(&mesh.vertices, face);
        for &v in face {
            vnorm[v] += n;
        }
    }
    let mut out = mesh.clone();
    let inner_off = out.vertices.len();
    for (v, n) in mesh.vertices.iter().zip(vnorm.iter()) {
        let mut nn = *n;
        if nn.norm() > 1e-12 {
            nn = nn.normalize();
        } else {
            nn = Vector3::z();
        }
        out.vertices.push(*v - nn * thickness);
    }
    // Inner faces — reverse winding so their normal flips inward.
    for face in &mesh.faces {
        let mut inner: Vec<usize> = face.iter().map(|i| inner_off + i).collect();
        inner.reverse();
        out.faces.push(inner);
    }
    // Stitch boundaries — every edge that appears in exactly one
    // face becomes a quad joining outer to inner.
    let mut edge_count: std::collections::BTreeMap<(usize, usize), usize> = Default::default();
    for face in &mesh.faces {
        for k in 0..face.len() {
            let a = face[k];
            let b = face[(k + 1) % face.len()];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    for (e, c) in &edge_count {
        if *c == 1 {
            let a_outer = e.0;
            let b_outer = e.1;
            let a_inner = inner_off + e.0;
            let b_inner = inner_off + e.1;
            out.faces.push(vec![a_outer, b_outer, b_inner, a_inner]);
        }
    }
    Ok(out)
}
