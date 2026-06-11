//! Edge bevel — for each edge in `edge_ids` (given as vertex pairs),
//! split the edge into `segments` pieces and pull each split vertex
//! `distance` along the average face normal at that vertex.
//!
//! v1 implementation: produces extra vertices along each input edge
//! that are offset along the *vertex normal*. The original mesh
//! topology is preserved (faces still reference original vertex
//! indices); the bevel state is captured in the new vertex list.
//! Subsequent ops can then read the new vertex positions.

use nalgebra::Vector3;

use crate::error::BlenderOpError;
use crate::mesh::{face_normal, Mesh};

/// Bevel the given edges. Returns the new mesh.
pub fn edges(
    mesh: &Mesh,
    edge_ids: &[(usize, usize)],
    distance: f64,
    segments: usize,
) -> Result<Mesh, BlenderOpError> {
    if !distance.is_finite() || distance < 0.0 {
        return Err(BlenderOpError::BadParameter {
            name: "distance",
            reason: format!("must be finite and >= 0 (got {distance})"),
        });
    }
    if segments == 0 {
        return Err(BlenderOpError::BadParameter {
            name: "segments",
            reason: "must be > 0".into(),
        });
    }
    let nv = mesh.vertices.len();
    for e in edge_ids {
        if e.0 >= nv {
            return Err(BlenderOpError::IndexOutOfRange {
                kind: "vertex",
                idx: e.0,
                limit: nv,
            });
        }
        if e.1 >= nv {
            return Err(BlenderOpError::IndexOutOfRange {
                kind: "vertex",
                idx: e.1,
                limit: nv,
            });
        }
    }

    let mut out = mesh.clone();

    // Vertex normal = sum of adjacent face normals, renormalised.
    let mut vnorm = vec![Vector3::<f64>::zeros(); mesh.vertices.len()];
    let mut vcount = vec![0u32; mesh.vertices.len()];
    for face in &mesh.faces {
        let n = face_normal(&mesh.vertices, face);
        for &v in face {
            vnorm[v] += n;
            vcount[v] += 1;
        }
    }
    for (n, c) in vnorm.iter_mut().zip(vcount.iter()) {
        if *c > 0 {
            *n /= f64::from(*c);
            let l = n.norm();
            if l > 1e-12 {
                *n /= l;
            }
        }
    }

    // For each edge: spawn `segments + 1` evenly-spaced points along
    // the edge, each offset along the interpolated vertex normal.
    for e in edge_ids {
        let a = mesh.vertices[e.0];
        let b = mesh.vertices[e.1];
        let na = vnorm[e.0];
        let nb = vnorm[e.1];
        for k in 0..=segments {
            let t = k as f64 / segments as f64;
            let p = a * (1.0 - t) + b * t;
            // Guard the normalize: an edge endpoint that belongs to no
            // face (vcount == 0) keeps a zero vertex normal. Guarding the
            // *interpolated* normal `raw` (rather than na/nb individually)
            // also covers the mid-edge case where na and nb cancel -- a
            // saddle edge with antiparallel endpoint normals makes `raw`
            // near-zero even when both endpoints are unit. `Vector3::
            // normalize` of a zero vector is NaN, which would poison the
            // emitted vertex; fall back to no offset (the on-edge point)
            // -- the same `> 1e-12` guard used for the per-vertex normals
            // above and in inset/solidify.
            let raw = na * (1.0 - t) + nb * t;
            let l = raw.norm();
            let n = if l > 1e-12 {
                raw / l
            } else {
                Vector3::zeros()
            };
            out.vertices.push(p + n * distance);
        }
    }

    Ok(out)
}
