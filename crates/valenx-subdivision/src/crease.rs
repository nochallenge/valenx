//! Semi-sharp creases — annotate a set of edges with sharpness
//! values so the Catmull-Clark scheme can preserve sharp features
//! while still smoothing everywhere else.
//!
//! Each subdivision iteration decrements all crease sharpness values
//! by 1; an edge with sharpness `> 0` uses the **infinitely sharp**
//! rule (edge point = midpoint, adjacent vertices clamp toward the
//! crease tangent). When sharpness reaches 0 the edge reverts to
//! smooth Catmull-Clark.

use serde::{Deserialize, Serialize};

use crate::error::SubdivError;
use crate::mesh::SubdivMesh;

/// Mesh + per-edge sharpness annotation.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct CreasedMesh {
    /// Underlying polygon mesh.
    pub mesh: SubdivMesh,
    /// `(low, high) -> sharpness >= 0.0`. Edges absent from the map
    /// are smooth.
    pub creases: std::collections::BTreeMap<(usize, usize), f64>,
}

impl CreasedMesh {
    /// New creased mesh with no sharp edges.
    pub fn from_mesh(m: SubdivMesh) -> Self {
        Self {
            mesh: m,
            creases: Default::default(),
        }
    }

    /// Mark the listed edges sharp at the given sharpness. Each edge
    /// is given as a vertex-index pair; we sort so the map key is
    /// (low, high) regardless of input order.
    pub fn set(&mut self, edges: &[(usize, usize)], sharpness: f64) -> Result<(), SubdivError> {
        if !sharpness.is_finite() || sharpness < 0.0 {
            return Err(SubdivError::BadParameter {
                name: "sharpness",
                reason: format!("must be finite and >= 0 (got {sharpness})"),
            });
        }
        let nv = self.mesh.n_verts();
        for (a, b) in edges {
            if *a >= nv {
                return Err(SubdivError::IndexOutOfRange {
                    kind: "vertex",
                    idx: *a,
                    limit: nv,
                });
            }
            if *b >= nv {
                return Err(SubdivError::IndexOutOfRange {
                    kind: "vertex",
                    idx: *b,
                    limit: nv,
                });
            }
            let key = if *a < *b { (*a, *b) } else { (*b, *a) };
            self.creases.insert(key, sharpness);
        }
        Ok(())
    }

    /// Lookup sharpness for an edge (returns 0.0 if smooth).
    pub fn sharpness(&self, a: usize, b: usize) -> f64 {
        let key = if a < b { (a, b) } else { (b, a) };
        *self.creases.get(&key).unwrap_or(&0.0)
    }
}

/// Top-level constructor — annotate `edges` on `mesh` with `sharpness`.
pub fn set_crease(
    mesh: &SubdivMesh,
    edges: &[(usize, usize)],
    sharpness: f64,
) -> Result<CreasedMesh, SubdivError> {
    let mut c = CreasedMesh::from_mesh(mesh.clone());
    c.set(edges, sharpness)?;
    Ok(c)
}
