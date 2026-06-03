//! Surface3d — ruled-surface descriptor.
//!
//! gCAD3D's ruled surface is the surface swept by a straight line
//! whose endpoints trace two input curves. v1 keeps the descriptor
//! abstract — the actual surface mesh / BRep is produced by
//! `valenx-surface::ruled` when invoked from the host.

use serde::{Deserialize, Serialize};

/// Ruled-surface descriptor referencing two curves by id.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RuledSurface {
    /// First curve id (free-form — opaque to this crate).
    pub curve1: String,
    /// Second curve id.
    pub curve2: String,
}

/// Build the descriptor — re-exposes the ruled-surface concept
/// under the gCAD3D namespace.
pub fn ruled(curve1: impl Into<String>, curve2: impl Into<String>) -> RuledSurface {
    RuledSurface {
        curve1: curve1.into(),
        curve2: curve2.into(),
    }
}
