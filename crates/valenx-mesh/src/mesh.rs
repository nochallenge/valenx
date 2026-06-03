//! The canonical `Mesh` type.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::element::ElementBlock;
use crate::region::{BoundaryGroup, Region};
use crate::stats::MeshStats;

/// A canonical finite-element / finite-volume mesh.
///
/// Layout:
/// - `nodes` is a flat vector of 3D coordinates (padded with zero-Z
///   for 2D meshes — coordinate dimension is always 3).
/// - `element_blocks` hold elements of one type each; mixed-type
///   meshes simply use several blocks.
/// - `regions` and `boundaries` index into the element arrays.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Mesh {
    pub id: String,
    pub nodes: Vec<Vector3<f64>>,
    pub element_blocks: Vec<ElementBlock>,
    pub regions: Vec<Region>,
    pub boundaries: Vec<BoundaryGroup>,
    /// Cached statistics — recomputed via [`Mesh::recompute_stats`].
    pub stats: MeshStats,
}

impl Mesh {
    /// Create an empty mesh with the given ID.
    pub fn new(id: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            ..Default::default()
        }
    }

    /// Total element count across all blocks.
    pub fn total_elements(&self) -> usize {
        self.element_blocks.iter().map(|b| b.count()).sum()
    }

    /// Refresh cached counts in `stats`. Cheap walk over the
    /// existing element blocks — does NOT compute geometric quality
    /// metrics. Use [`Mesh::recompute_quality_stats`] to also fill
    /// in `min_element_size` / `max_aspect_ratio` / `max_skewness` /
    /// `min_orthogonality`.
    pub fn recompute_stats(&mut self) {
        self.stats.node_count = self.nodes.len() as u64;
        self.stats.element_count = self.total_elements() as u64;
        self.stats.region_count = self.regions.len() as u32;
        self.stats.boundary_group_count = self.boundaries.len() as u32;
    }

    /// Recompute counts AND geometric quality scalars
    /// (`min_element_size`, `max_aspect_ratio`, `max_skewness`,
    /// `min_orthogonality`) into `stats`, and return the full
    /// [`crate::quality::QualityReport`] so callers can also
    /// inspect / cache / surface it without re-running the analysis.
    ///
    /// Calls [`crate::quality::report`] internally — O(elements)
    /// for the per-element scalars plus O(faces) for orthogonality.
    /// Adapters that import a mesh from disk should call this
    /// instead of the bare [`Mesh::recompute_stats`] when they want
    /// quality fields populated.
    pub fn recompute_quality_stats(&mut self) -> crate::quality::QualityReport {
        self.recompute_stats();
        let report = crate::quality::report(self);
        self.stats.min_element_size = report.min_size;
        self.stats.max_aspect_ratio = report.max_aspect_ratio;
        self.stats.max_skewness = report.max_skewness;
        self.stats.min_orthogonality = report.min_orthogonality;
        report
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::{ElementBlock, ElementType};

    #[test]
    fn empty_mesh_stats() {
        let mut m = Mesh::new("test");
        m.recompute_stats();
        assert_eq!(m.stats.node_count, 0);
        assert_eq!(m.stats.element_count, 0);
    }

    #[test]
    fn tri_mesh_stats() {
        let mut m = Mesh::new("triangle");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        m.element_blocks.push(block);
        m.recompute_stats();
        assert_eq!(m.stats.node_count, 3);
        assert_eq!(m.stats.element_count, 1);
    }

    #[test]
    fn recompute_quality_stats_populates_aspect_skew_size() {
        // Right-isoceles triangle: skew=0.25, aspect=sqrt(2),
        // size (area) = 0.5. recompute_quality_stats should fill
        // these into the cached stats AND return the report.
        let mut m = Mesh::new("right-iso");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        m.element_blocks.push(block);
        let report = m.recompute_quality_stats();
        // Counts populated as before.
        assert_eq!(m.stats.node_count, 3);
        assert_eq!(m.stats.element_count, 1);
        // Quality scalars now populated too.
        assert!((m.stats.min_element_size.unwrap() - 0.5).abs() < 1e-12);
        assert!((m.stats.max_aspect_ratio.unwrap() - (2.0_f64).sqrt()).abs() < 1e-12);
        assert!((m.stats.max_skewness.unwrap() - 0.25).abs() < 1e-12);
        // No interior faces -> orthogonality stays None.
        assert_eq!(m.stats.min_orthogonality, None);
        // Returned report mirrors the cached scalars.
        assert_eq!(report.element_count, 1);
        assert_eq!(report.max_aspect_ratio, m.stats.max_aspect_ratio);
        assert_eq!(report.max_skewness, m.stats.max_skewness);
    }

    #[test]
    fn recompute_quality_stats_populates_orthogonality_when_interior_faces_exist() {
        // Two stacked unit cubes share one face; orthogonality 1.0.
        let mut m = Mesh::new("two-stacked-hexes");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
            Vector3::new(0.0, 0.0, 2.0),
            Vector3::new(1.0, 0.0, 2.0),
            Vector3::new(1.0, 1.0, 2.0),
            Vector3::new(0.0, 1.0, 2.0),
        ];
        let mut block = ElementBlock::new(ElementType::Hex8);
        block.connectivity = vec![0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11];
        m.element_blocks.push(block);
        let report = m.recompute_quality_stats();
        assert!((m.stats.min_orthogonality.unwrap() - 1.0).abs() < 1e-12);
        // Returned report has the same value (no recomputation needed).
        assert_eq!(report.min_orthogonality, m.stats.min_orthogonality);
    }
}
