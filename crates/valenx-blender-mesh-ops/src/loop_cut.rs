//! Loop cut — inserts `n_cuts` vertices along each edge of an
//! edge-loop and re-splits the faces. v1: just inserts the new
//! vertices along the edges and returns the mesh; the face split
//! is performed by [`crate::extrude::region`] downstream so the
//! op stays composable.
//!
//! Blender's actual loop cut threads through perpendicular face
//! diagonals; v1 here is the "vertex-only" mode (Ctrl+R then Esc),
//! which is what the modeler reaches for to add construction
//! geometry without committing to a face split.

use crate::error::BlenderOpError;
use crate::mesh::Mesh;

/// Insert `n_cuts` new vertices along each edge in `edge_loop`.
pub fn insert(
    mesh: &Mesh,
    edge_loop: &[(usize, usize)],
    n_cuts: usize,
) -> Result<Mesh, BlenderOpError> {
    if n_cuts == 0 {
        return Err(BlenderOpError::BadParameter {
            name: "n_cuts",
            reason: "must be > 0".into(),
        });
    }
    let nv = mesh.vertices.len();
    for e in edge_loop {
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
    for e in edge_loop {
        let a = mesh.vertices[e.0];
        let b = mesh.vertices[e.1];
        for k in 1..=n_cuts {
            let t = k as f64 / (n_cuts as f64 + 1.0);
            out.vertices.push(a * (1.0 - t) + b * t);
        }
    }
    Ok(out)
}
