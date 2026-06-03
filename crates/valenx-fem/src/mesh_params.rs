//! Mesh-generation parameters.

use serde::{Deserialize, Serialize};

/// Element-order enum.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum ElementOrder {
    /// 4-node linear tetrahedron.
    Tet4,
    /// 10-node quadratic tetrahedron.
    Tet10,
    /// 8-node linear hexahedron (brick).
    Hex8,
    /// 20-node quadratic hexahedron.
    Hex20,
}

impl ElementOrder {
    /// CalculiX / Abaqus element-type code.
    pub fn ccx_code(&self) -> &'static str {
        match self {
            ElementOrder::Tet4 => "C3D4",
            ElementOrder::Tet10 => "C3D10",
            ElementOrder::Hex8 => "C3D8",
            ElementOrder::Hex20 => "C3D20",
        }
    }
}

/// Region of the solid where the user wants finer / coarser mesh.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RefinementRegion {
    /// 0-based face indices whose neighbourhood gets refined.
    pub face_indices: Vec<usize>,
    /// Target element size in this region (model units).
    pub element_size: f64,
}

/// Mesh parameters fed to the upstream meshing adapter (gmsh /
/// Netgen).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FemMeshParams {
    /// Global target element size in model units (metres typically).
    pub element_size: f64,
    /// Element order + topology.
    pub element_type: ElementOrder,
    /// Optional local refinement regions.
    pub refinement_regions: Vec<RefinementRegion>,
}

impl Default for FemMeshParams {
    fn default() -> Self {
        Self {
            element_size: 0.005,
            element_type: ElementOrder::Tet10,
            refinement_regions: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ccx_codes_are_standard() {
        assert_eq!(ElementOrder::Tet4.ccx_code(), "C3D4");
        assert_eq!(ElementOrder::Tet10.ccx_code(), "C3D10");
        assert_eq!(ElementOrder::Hex8.ccx_code(), "C3D8");
        assert_eq!(ElementOrder::Hex20.ccx_code(), "C3D20");
    }

    #[test]
    fn default_is_quadratic_tet() {
        let p = FemMeshParams::default();
        assert_eq!(p.element_type, ElementOrder::Tet10);
        assert!(p.element_size > 0.0);
        assert!(p.refinement_regions.is_empty());
    }
}
