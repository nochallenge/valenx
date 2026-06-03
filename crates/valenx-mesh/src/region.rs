//! Regions and boundary groups within a mesh.

use serde::{Deserialize, Serialize};

/// A named volumetric region inside the mesh (e.g. `"fluid"`,
/// `"solid"`, `"air"`, `"steel"`). Used to target material
/// assignment and per-region post-processing.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Region {
    pub name: String,
    /// Element indices (into the flat mesh element array) that
    /// belong to this region.
    pub element_indices: Vec<u32>,
    /// Optional stable ID shared with `RegionRef` on the fields
    /// side; defaults to `name` when not set.
    pub id: Option<String>,
}

/// A named surface group used for boundary-condition targeting
/// (e.g. `"inlet"`, `"outlet"`, `"walls"`).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BoundaryGroup {
    pub name: String,
    /// Element indices whose faces lie on this boundary. For CFD
    /// boundaries these are face elements; for FEA they may be node
    /// sets — use `kind` to disambiguate.
    pub element_indices: Vec<u32>,
    pub kind: BoundaryKind,
}

/// Whether a boundary group is element-face-based or node-based.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum BoundaryKind {
    /// Faces of surface elements (typical CFD).
    Faces,
    /// Node set (typical for FEA constraints / loads).
    Nodes,
}
