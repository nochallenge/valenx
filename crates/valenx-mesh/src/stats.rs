//! Cached mesh quality / size statistics. Computed once at mesh load
//! and surfaced in the UI + the browser tree.

use serde::{Deserialize, Serialize};

/// Lightweight statistics describing a mesh.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct MeshStats {
    pub node_count: u64,
    pub element_count: u64,
    pub region_count: u32,
    pub boundary_group_count: u32,
    /// Minimum element volume / area (sign preserved — negative
    /// values indicate inverted elements).
    pub min_element_size: Option<f64>,
    /// Maximum element aspect ratio.
    pub max_aspect_ratio: Option<f64>,
    /// Worst (largest) per-element equiangle skewness, range
    /// `[0, 1]`. `0` = all elements regular, approaches `1` as
    /// faces degenerate. See `crate::quality::equiangle_skewness`.
    pub max_skewness: Option<f64>,
    /// Minimum element orthogonality (cosine of angle between face
    /// normal and cell-centre vector, for CFD-style meshes).
    pub min_orthogonality: Option<f64>,
}
