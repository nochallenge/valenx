//! Element topology types.

use serde::{Deserialize, Serialize};

/// Supported element topologies. Covers the overwhelming majority of
/// practical CFD / FEA / EM meshes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ElementType {
    /// 1D line (2 nodes).
    Line2,
    /// 2D triangle (3 nodes).
    Tri3,
    /// 2D quadrilateral (4 nodes).
    Quad4,
    /// 3D tetrahedron (4 nodes).
    Tet4,
    /// 3D pyramid (5 nodes).
    Pyr5,
    /// 3D prism / wedge (6 nodes).
    Prism6,
    /// 3D hexahedron (8 nodes).
    Hex8,
    /// 2D triangle — second order (6 nodes).
    Tri6,
    /// 3D tetrahedron — second order (10 nodes).
    Tet10,
    /// 3D hexahedron — second order (20 nodes).
    Hex20,
}

impl ElementType {
    /// Number of nodes per element of this type (e.g. `3` for
    /// [`ElementType::Tri3`], `20` for [`ElementType::Hex20`]).
    pub fn nodes_per_element(self) -> usize {
        match self {
            ElementType::Line2 => 2,
            ElementType::Tri3 => 3,
            ElementType::Quad4 | ElementType::Tet4 => 4,
            ElementType::Pyr5 => 5,
            ElementType::Prism6 | ElementType::Tri6 => 6,
            ElementType::Hex8 => 8,
            ElementType::Tet10 => 10,
            ElementType::Hex20 => 20,
        }
    }

    /// Topological dimension (1/2/3).
    pub fn dim(self) -> u8 {
        match self {
            ElementType::Line2 => 1,
            ElementType::Tri3 | ElementType::Quad4 | ElementType::Tri6 => 2,
            ElementType::Tet4
            | ElementType::Pyr5
            | ElementType::Prism6
            | ElementType::Hex8
            | ElementType::Tet10
            | ElementType::Hex20 => 3,
        }
    }
}

/// A contiguous block of elements of a single type. A `Mesh` holds
/// one or more of these, letting us efficiently store mixed-element
/// meshes without a type-per-element tag.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ElementBlock {
    pub element_type: ElementType,
    /// Flat connectivity: `nodes_per_element * count` node indices.
    pub connectivity: Vec<u32>,
}

impl ElementBlock {
    /// New empty block (no elements). The element type is fixed at
    /// construction.
    pub fn new(element_type: ElementType) -> Self {
        Self {
            element_type,
            connectivity: Vec::new(),
        }
    }

    /// Number of elements in this block.
    pub fn count(&self) -> usize {
        let n = self.element_type.nodes_per_element();
        if n == 0 {
            0
        } else {
            self.connectivity.len() / n
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn node_counts() {
        assert_eq!(ElementType::Tri3.nodes_per_element(), 3);
        assert_eq!(ElementType::Tet4.nodes_per_element(), 4);
        assert_eq!(ElementType::Hex8.nodes_per_element(), 8);
        assert_eq!(ElementType::Hex20.nodes_per_element(), 20);
    }

    #[test]
    fn dim_lookup() {
        assert_eq!(ElementType::Line2.dim(), 1);
        assert_eq!(ElementType::Tri3.dim(), 2);
        assert_eq!(ElementType::Tet4.dim(), 3);
    }

    #[test]
    fn block_count() {
        let mut b = ElementBlock::new(ElementType::Tri3);
        b.connectivity = vec![0, 1, 2, 2, 1, 3];
        assert_eq!(b.count(), 2);
    }
}
