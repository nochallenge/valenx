//! Cell-face / cell-edge adjacency for finite-volume / FEA quality
//! metrics.
//!
//! Builds an index of every (d-1)-dimensional sub-entity in the mesh
//! and groups elements that share each one. Interior entries have
//! two owners; boundary entries have one. This is the data CFD-style
//! mesh quality metrics (orthogonality, non-orthogonality angle,
//! skewness-by-face) need to walk neighbouring cells, and is what
//! boundary-detection / mesh-dual / interface-mass-conservation
//! checks build on top of.
//!
//! - **Face adjacency** ([`build_face_adjacency`]) covers 3D elements
//!   (Tet4 / Hex8 / Pyr5 / Prism6 + their quadratic counterparts).
//!   Each face is a 2D polygon with 3 or 4 nodes.
//! - **Edge adjacency** ([`build_edge_adjacency`]) covers 2D elements
//!   (Tri3 / Quad4 + Tri6). Each edge is a 1D line with 2 nodes.
//!   Used for shell mesh quality, 2D-mesh boundary extraction, etc.

use std::collections::HashMap;

use crate::{ElementType, Mesh};

/// A reference to one element-and-one-of-its-faces. The pair
/// `(global_element_index, local_face_index)` is enough to recover
/// the face's nodes (via [`face_nodes_for`]) and the element's
/// type / centroid.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ElementFaceRef {
    /// Globally-unique element index. Elements are numbered by
    /// walking `mesh.element_blocks` in order; the first element of
    /// the second block has the global index of (count of block 0).
    pub global_element: usize,
    /// Which face of that element — 0..faces_per_element(type).
    pub local_face: u8,
}

/// Cell-face adjacency table. Each unique face in the mesh appears
/// exactly once, owned by either one element (boundary face) or two
/// (interior face).
#[derive(Clone, Debug, Default)]
pub struct FaceAdjacency {
    interior: Vec<InteriorFace>,
    boundary: Vec<BoundaryFace>,
}

/// A face shared by exactly two elements.
#[derive(Clone, Debug)]
pub struct InteriorFace {
    /// Sorted node indices identifying the face. Always sorted
    /// ascending so two adjacent cells produce the same key.
    pub sorted_nodes: Vec<u32>,
    pub left: ElementFaceRef,
    pub right: ElementFaceRef,
}

/// A face that belongs to a single element (i.e. on the mesh boundary).
#[derive(Clone, Debug)]
pub struct BoundaryFace {
    pub sorted_nodes: Vec<u32>,
    pub owner: ElementFaceRef,
}

impl FaceAdjacency {
    /// Number of interior (shared-by-two-elements) faces.
    pub fn interior_face_count(&self) -> usize {
        self.interior.len()
    }

    /// Number of boundary (single-owner) faces.
    pub fn boundary_face_count(&self) -> usize {
        self.boundary.len()
    }

    /// Borrow the interior face list.
    pub fn interior_faces(&self) -> &[InteriorFace] {
        &self.interior
    }

    /// Borrow the boundary face list.
    pub fn boundary_faces(&self) -> &[BoundaryFace] {
        &self.boundary
    }
}

/// Build the face adjacency for `mesh`. 3D-only: 1D / 2D blocks
/// (Line2 / Tri3 / Quad4 / Tri6) contribute no faces today and are
/// silently skipped. Quadratic 3D types (Tet10 / Hex20) reduce to
/// their corner subset, mirroring the quality-metric convention.
pub fn build_face_adjacency(mesh: &Mesh) -> FaceAdjacency {
    // Bucket: sorted-node-key -> Vec<ElementFaceRef> (1 or 2 entries
    // for well-formed meshes; >2 means non-manifold input which we
    // surface as separate boundary faces rather than crashing).
    let mut faces: HashMap<Vec<u32>, Vec<ElementFaceRef>> = HashMap::new();
    let mut global_element: usize = 0;
    for block in &mesh.element_blocks {
        let face_table = face_nodes_for(block.element_type);
        let npe = block.element_type.nodes_per_element();
        if npe == 0 || face_table.is_empty() {
            // Skip element types without a face decomposition (Line2,
            // 2D, etc.) — they advance the global counter but
            // contribute no entries.
            global_element += block.connectivity.len() / npe.max(1);
            continue;
        }
        let element_count = block.connectivity.len() / npe;
        for i in 0..element_count {
            let start = i * npe;
            let conn = &block.connectivity[start..start + npe];
            for (local_face, face_indices) in face_table.iter().enumerate() {
                let mut sorted: Vec<u32> = face_indices.iter().map(|&idx| conn[idx]).collect();
                sorted.sort_unstable();
                let entry = faces.entry(sorted).or_default();
                entry.push(ElementFaceRef {
                    global_element,
                    local_face: local_face as u8,
                });
            }
            global_element += 1;
        }
    }

    let mut adj = FaceAdjacency::default();
    for (sorted_nodes, owners) in faces {
        match owners.as_slice() {
            [solo] => adj.boundary.push(BoundaryFace {
                sorted_nodes,
                owner: *solo,
            }),
            [a, b] => adj.interior.push(InteriorFace {
                sorted_nodes,
                left: *a,
                right: *b,
            }),
            // Non-manifold (>2 owners): surface every owner as its
            // own boundary entry rather than dropping the data. This
            // lets quality reporters flag the mesh as bad without
            // panicking on malformed input.
            multiple => {
                for o in multiple {
                    adj.boundary.push(BoundaryFace {
                        sorted_nodes: sorted_nodes.clone(),
                        owner: *o,
                    });
                }
            }
        }
    }
    // Stable ordering for reproducibility. Sort interior by left's
    // global element, boundary by owner's global element.
    adj.interior
        .sort_by_key(|f| (f.left.global_element, f.left.local_face));
    adj.boundary
        .sort_by_key(|f| (f.owner.global_element, f.owner.local_face));
    adj
}

/// Local node indices for each face of an element. Indices are 0..npe
/// — actual node IDs come from indexing the element's connectivity
/// with these. Empty slice for types we don't decompose into faces
/// (1D / 2D today, plus unsupported types).
pub fn face_nodes_for(t: ElementType) -> &'static [&'static [usize]] {
    match t {
        ElementType::Tet4 => TET4_FACES,
        ElementType::Hex8 => HEX8_FACES,
        ElementType::Pyr5 => PYR5_FACES,
        ElementType::Prism6 => PRISM6_FACES,
        // Quadratic 3D: corner subset matches linear face decomposition.
        ElementType::Tet10 => TET4_FACES,
        ElementType::Hex20 => HEX8_FACES,
        // Line / 2D / unsupported: no faces.
        ElementType::Line2 | ElementType::Tri3 | ElementType::Quad4 | ElementType::Tri6 => &[],
    }
}

/// A reference to one element-and-one-of-its-edges. Edge analog of
/// [`ElementFaceRef`] for 2D meshes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ElementEdgeRef {
    pub global_element: usize,
    pub local_edge: u8,
}

/// Cell-edge adjacency table for 2D meshes. Mirrors
/// [`FaceAdjacency`] but for edges of 2D elements (Tri3 / Quad4 /
/// Tri6) rather than faces of 3D elements.
#[derive(Clone, Debug, Default)]
pub struct EdgeAdjacency {
    interior: Vec<InteriorEdge>,
    boundary: Vec<BoundaryEdge>,
}

/// An edge shared by exactly two 2D elements.
#[derive(Clone, Debug)]
pub struct InteriorEdge {
    /// Sorted node indices identifying the edge (always 2 nodes).
    pub sorted_nodes: Vec<u32>,
    pub left: ElementEdgeRef,
    pub right: ElementEdgeRef,
}

/// An edge that belongs to a single 2D element (mesh boundary).
#[derive(Clone, Debug)]
pub struct BoundaryEdge {
    pub sorted_nodes: Vec<u32>,
    pub owner: ElementEdgeRef,
}

impl EdgeAdjacency {
    /// Number of interior (shared-by-two-elements) edges.
    pub fn interior_edge_count(&self) -> usize {
        self.interior.len()
    }

    /// Number of boundary (single-owner) edges.
    pub fn boundary_edge_count(&self) -> usize {
        self.boundary.len()
    }

    /// Borrow the interior edge list.
    pub fn interior_edges(&self) -> &[InteriorEdge] {
        &self.interior
    }

    /// Borrow the boundary edge list.
    pub fn boundary_edges(&self) -> &[BoundaryEdge] {
        &self.boundary
    }
}

/// Build the edge adjacency for `mesh`. 2D-only: 3D blocks have
/// faces, not edges, so they're silently skipped (advance the
/// global counter, contribute nothing). Tri6 reduces to its corner
/// subset (Tri3) so quadratic 2D elements share the linear edge
/// topology.
pub fn build_edge_adjacency(mesh: &Mesh) -> EdgeAdjacency {
    let mut edges: HashMap<Vec<u32>, Vec<ElementEdgeRef>> = HashMap::new();
    let mut global_element: usize = 0;
    for block in &mesh.element_blocks {
        let edge_table = edge_nodes_for(block.element_type);
        let npe = block.element_type.nodes_per_element();
        if npe == 0 || edge_table.is_empty() {
            global_element += block.connectivity.len() / npe.max(1);
            continue;
        }
        let element_count = block.connectivity.len() / npe;
        for i in 0..element_count {
            let start = i * npe;
            let conn = &block.connectivity[start..start + npe];
            for (local_edge, edge_indices) in edge_table.iter().enumerate() {
                let mut sorted: Vec<u32> = edge_indices.iter().map(|&idx| conn[idx]).collect();
                sorted.sort_unstable();
                let entry = edges.entry(sorted).or_default();
                entry.push(ElementEdgeRef {
                    global_element,
                    local_edge: local_edge as u8,
                });
            }
            global_element += 1;
        }
    }

    let mut adj = EdgeAdjacency::default();
    for (sorted_nodes, owners) in edges {
        match owners.as_slice() {
            [solo] => adj.boundary.push(BoundaryEdge {
                sorted_nodes,
                owner: *solo,
            }),
            [a, b] => adj.interior.push(InteriorEdge {
                sorted_nodes,
                left: *a,
                right: *b,
            }),
            multiple => {
                // Non-manifold (T-junctions, dangling overlaps): split
                // each owner into its own boundary entry rather than
                // panic. Same policy as build_face_adjacency.
                for o in multiple {
                    adj.boundary.push(BoundaryEdge {
                        sorted_nodes: sorted_nodes.clone(),
                        owner: *o,
                    });
                }
            }
        }
    }
    adj.interior
        .sort_by_key(|e| (e.left.global_element, e.left.local_edge));
    adj.boundary
        .sort_by_key(|e| (e.owner.global_element, e.owner.local_edge));
    adj
}

/// Local node-pair indices for each edge of a 2D element. Empty
/// for 1D / 3D elements (their relevant adjacency is face-level
/// or doesn't exist).
pub fn edge_nodes_for(t: ElementType) -> &'static [&'static [usize]] {
    match t {
        ElementType::Tri3 => TRI3_EDGES,
        ElementType::Quad4 => QUAD4_EDGES,
        // Quadratic 2D: corner subset matches linear edges (mid-edge
        // nodes don't change topological adjacency).
        ElementType::Tri6 => TRI3_EDGES,
        // 1D / 3D: no edges.
        _ => &[],
    }
}

const TRI3_EDGES: &[&[usize]] = &[&[0, 1], &[1, 2], &[2, 0]];
const QUAD4_EDGES: &[&[usize]] = &[&[0, 1], &[1, 2], &[2, 3], &[3, 0]];

const TET4_FACES: &[&[usize]] = &[&[0, 1, 2], &[0, 1, 3], &[0, 2, 3], &[1, 2, 3]];

const HEX8_FACES: &[&[usize]] = &[
    &[0, 1, 2, 3], // bottom (-z)
    &[4, 5, 6, 7], // top    (+z)
    &[0, 1, 5, 4], // front  (-y)
    &[1, 2, 6, 5], // right  (+x)
    &[2, 3, 7, 6], // back   (+y)
    &[3, 0, 4, 7], // left   (-x)
];

const PYR5_FACES: &[&[usize]] = &[
    &[0, 1, 2, 3], // base (quad)
    &[0, 1, 4],
    &[1, 2, 4],
    &[2, 3, 4],
    &[3, 0, 4],
];

const PRISM6_FACES: &[&[usize]] = &[
    &[0, 1, 2], // bottom cap
    &[3, 4, 5], // top cap
    &[0, 1, 4, 3],
    &[1, 2, 5, 4],
    &[2, 0, 3, 5],
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ElementBlock;
    use nalgebra::Vector3;

    #[test]
    fn single_tet_has_four_boundary_faces_and_no_interior() {
        let mut m = Mesh::new("single-tet");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tet4);
        block.connectivity = vec![0, 1, 2, 3];
        m.element_blocks.push(block);

        let adj = build_face_adjacency(&m);
        assert_eq!(adj.interior_face_count(), 0);
        assert_eq!(adj.boundary_face_count(), 4);
    }

    #[test]
    fn single_hex_has_six_boundary_faces() {
        let mut m = Mesh::new("hex");
        m.nodes = unit_cube_nodes();
        let mut block = ElementBlock::new(ElementType::Hex8);
        block.connectivity = vec![0, 1, 2, 3, 4, 5, 6, 7];
        m.element_blocks.push(block);
        let adj = build_face_adjacency(&m);
        assert_eq!(adj.boundary_face_count(), 6);
        assert_eq!(adj.interior_face_count(), 0);
    }

    #[test]
    fn single_pyramid_has_one_quad_base_and_four_tri_sides() {
        let mut m = Mesh::new("pyr");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.5, 0.5, 1.0),
        ];
        let mut block = ElementBlock::new(ElementType::Pyr5);
        block.connectivity = vec![0, 1, 2, 3, 4];
        m.element_blocks.push(block);
        let adj = build_face_adjacency(&m);
        assert_eq!(adj.boundary_face_count(), 5);
        assert_eq!(adj.interior_face_count(), 0);
        // Exactly one face has 4 nodes (quad base); the other 4 are tris.
        let quad_count = adj
            .boundary_faces()
            .iter()
            .filter(|f| f.sorted_nodes.len() == 4)
            .count();
        assert_eq!(quad_count, 1);
    }

    #[test]
    fn single_prism_has_two_tri_caps_and_three_quad_sides() {
        let mut m = Mesh::new("prism");
        let h = (3.0_f64).sqrt() / 2.0;
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.5, h, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(0.5, h, 1.0),
        ];
        let mut block = ElementBlock::new(ElementType::Prism6);
        block.connectivity = vec![0, 1, 2, 3, 4, 5];
        m.element_blocks.push(block);
        let adj = build_face_adjacency(&m);
        assert_eq!(adj.boundary_face_count(), 5);
        let tri_count = adj
            .boundary_faces()
            .iter()
            .filter(|f| f.sorted_nodes.len() == 3)
            .count();
        let quad_count = adj
            .boundary_faces()
            .iter()
            .filter(|f| f.sorted_nodes.len() == 4)
            .count();
        assert_eq!(tri_count, 2);
        assert_eq!(quad_count, 3);
    }

    #[test]
    fn two_stacked_hexes_share_one_interior_face() {
        // Hex A: unit cube z in [0,1]. Hex B: unit cube z in [1,2]
        // sharing the top face of A as its bottom face.
        let mut m = Mesh::new("two-hexes");
        m.nodes = vec![
            // Bottom of A (z=0): 0..3
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            // Middle layer (z=1, shared): 4..7
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
            // Top of B (z=2): 8..11
            Vector3::new(0.0, 0.0, 2.0),
            Vector3::new(1.0, 0.0, 2.0),
            Vector3::new(1.0, 1.0, 2.0),
            Vector3::new(0.0, 1.0, 2.0),
        ];
        let mut block = ElementBlock::new(ElementType::Hex8);
        block.connectivity = vec![
            0, 1, 2, 3, 4, 5, 6, 7, // hex A
            4, 5, 6, 7, 8, 9, 10, 11, // hex B
        ];
        m.element_blocks.push(block);
        let adj = build_face_adjacency(&m);
        assert_eq!(adj.interior_face_count(), 1);
        // 12 total faces - 2 (the shared one is now interior, not in
        // boundary count for either) = 10 boundary.
        assert_eq!(adj.boundary_face_count(), 10);
        let interior = &adj.interior_faces()[0];
        let mut expected = vec![4u32, 5, 6, 7];
        expected.sort_unstable();
        assert_eq!(interior.sorted_nodes, expected);
    }

    #[test]
    fn line_and_2d_blocks_contribute_no_faces() {
        let mut m = Mesh::new("2d-1d");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut tri = ElementBlock::new(ElementType::Tri3);
        tri.connectivity = vec![0, 1, 2];
        m.element_blocks.push(tri);
        let mut line = ElementBlock::new(ElementType::Line2);
        line.connectivity = vec![0, 1];
        m.element_blocks.push(line);
        let adj = build_face_adjacency(&m);
        assert_eq!(adj.interior_face_count(), 0);
        assert_eq!(adj.boundary_face_count(), 0);
    }

    #[test]
    fn tet10_reduces_to_tet4_corner_subset_for_face_count() {
        // A Tet10 with arbitrary mid-edge nodes should still yield 4
        // boundary triangular faces — same as the linear corner tet.
        let mut m = Mesh::new("tet10");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            // Mid-edge nodes 4..9 — irrelevant for face topology.
            Vector3::new(0.5, 0.0, 0.0),
            Vector3::new(0.5, 0.5, 0.0),
            Vector3::new(0.0, 0.5, 0.0),
            Vector3::new(0.0, 0.0, 0.5),
            Vector3::new(0.5, 0.0, 0.5),
            Vector3::new(0.0, 0.5, 0.5),
        ];
        let mut block = ElementBlock::new(ElementType::Tet10);
        block.connectivity = (0..10).collect();
        m.element_blocks.push(block);
        let adj = build_face_adjacency(&m);
        assert_eq!(adj.boundary_face_count(), 4);
        assert_eq!(adj.interior_face_count(), 0);
        // Each boundary face uses 3 corner nodes (no mid-edges).
        for f in adj.boundary_faces() {
            assert_eq!(f.sorted_nodes.len(), 3);
            assert!(f.sorted_nodes.iter().all(|&n| n < 4));
        }
    }

    #[test]
    fn single_tri_has_three_boundary_edges_and_no_interior() {
        let mut m = Mesh::new("tri");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2];
        m.element_blocks.push(block);
        let adj = build_edge_adjacency(&m);
        assert_eq!(adj.interior_edge_count(), 0);
        assert_eq!(adj.boundary_edge_count(), 3);
    }

    #[test]
    fn single_quad_has_four_boundary_edges() {
        let mut m = Mesh::new("quad");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Quad4);
        block.connectivity = vec![0, 1, 2, 3];
        m.element_blocks.push(block);
        let adj = build_edge_adjacency(&m);
        assert_eq!(adj.interior_edge_count(), 0);
        assert_eq!(adj.boundary_edge_count(), 4);
    }

    #[test]
    fn two_tris_sharing_an_edge_have_one_interior_and_four_boundary() {
        // Tri A = nodes 0-1-2. Tri B = nodes 1-3-2 (shares edge 1-2).
        let mut m = Mesh::new("two-tris");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2, 1, 3, 2];
        m.element_blocks.push(block);
        let adj = build_edge_adjacency(&m);
        assert_eq!(adj.interior_edge_count(), 1);
        assert_eq!(adj.boundary_edge_count(), 4);
        let f = &adj.interior_edges()[0];
        assert_eq!(f.sorted_nodes, vec![1, 2]);
        assert_ne!(f.left.global_element, f.right.global_element);
    }

    #[test]
    fn edge_adjacency_skips_3d_blocks() {
        // 3D elements have faces, not edges — edge_adjacency should
        // ignore them and return empty when only Tet4 is present.
        let mut m = Mesh::new("tet-only");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tet4);
        block.connectivity = vec![0, 1, 2, 3];
        m.element_blocks.push(block);
        let adj = build_edge_adjacency(&m);
        assert_eq!(adj.interior_edge_count(), 0);
        assert_eq!(adj.boundary_edge_count(), 0);
    }

    #[test]
    fn tri6_reduces_to_tri3_corner_subset_for_edge_count() {
        // A Tri6 with mid-edge nodes should still produce 3 edges
        // using only the 3 corner nodes — same as the linear Tri3.
        let mut m = Mesh::new("tri6");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.5, 0.0, 0.0),
            Vector3::new(0.5, 0.5, 0.0),
            Vector3::new(0.0, 0.5, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri6);
        block.connectivity = (0..6).collect();
        m.element_blocks.push(block);
        let adj = build_edge_adjacency(&m);
        assert_eq!(adj.boundary_edge_count(), 3);
        for e in adj.boundary_edges() {
            assert_eq!(e.sorted_nodes.len(), 2);
            // Edges only reference corner nodes 0..3.
            assert!(e.sorted_nodes.iter().all(|&n| n < 3));
        }
    }

    #[test]
    fn mixed_blocks_get_distinct_global_element_indices() {
        // One Tet4 + one Hex8 in separate blocks. Global element
        // indices should be 0 (tet) and 1 (hex). Each contributes
        // boundary-only faces (no shared faces by construction).
        let mut m = Mesh::new("mixed");
        let mut nodes = vec![
            // Tet (0..3)
            Vector3::new(10.0, 0.0, 0.0),
            Vector3::new(11.0, 0.0, 0.0),
            Vector3::new(10.0, 1.0, 0.0),
            Vector3::new(10.0, 0.0, 1.0),
        ];
        nodes.extend_from_slice(&unit_cube_nodes());
        m.nodes = nodes;
        let mut tet = ElementBlock::new(ElementType::Tet4);
        tet.connectivity = vec![0, 1, 2, 3];
        m.element_blocks.push(tet);
        let mut hex = ElementBlock::new(ElementType::Hex8);
        hex.connectivity = vec![4, 5, 6, 7, 8, 9, 10, 11];
        m.element_blocks.push(hex);
        let adj = build_face_adjacency(&m);
        assert_eq!(adj.interior_face_count(), 0);
        assert_eq!(adj.boundary_face_count(), 4 + 6);
        // Boundary faces include both element 0 (tet) and element 1 (hex).
        let element_set: std::collections::HashSet<usize> = adj
            .boundary_faces()
            .iter()
            .map(|f| f.owner.global_element)
            .collect();
        assert_eq!(element_set, [0usize, 1].iter().copied().collect());
    }

    fn unit_cube_nodes() -> Vec<Vector3<f64>> {
        vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ]
    }

    #[test]
    fn two_tets_sharing_a_face_have_one_interior_and_six_boundary() {
        // Tet A = (0,1,2,3) with face [1,2,3] facing +x.
        // Tet B = (4,1,2,3) where node 4 = (-1, 0.33, 0.33) sits on
        // the opposite side of face [1,2,3]. Shared face [1,2,3] is
        // the only interior face. Each tet contributes 3 boundary
        // faces -> 6 total.
        let mut m = Mesh::new("two-tets");
        m.nodes = vec![
            Vector3::new(1.0, 0.0, 0.0),    // 0
            Vector3::new(0.0, 0.0, 0.0),    // 1
            Vector3::new(0.0, 1.0, 0.0),    // 2
            Vector3::new(0.0, 0.0, 1.0),    // 3
            Vector3::new(-1.0, 0.33, 0.33), // 4 - opposite side
        ];
        let mut block = ElementBlock::new(ElementType::Tet4);
        block.connectivity = vec![0, 1, 2, 3, 4, 1, 2, 3];
        m.element_blocks.push(block);

        let adj = build_face_adjacency(&m);
        assert_eq!(adj.interior_face_count(), 1);
        assert_eq!(adj.boundary_face_count(), 6);

        // The interior face's left/right must reference distinct
        // global elements (0 and 1).
        let f = &adj.interior_faces()[0];
        assert_ne!(f.left.global_element, f.right.global_element);
        assert!(
            (f.left.global_element == 0 && f.right.global_element == 1)
                || (f.left.global_element == 1 && f.right.global_element == 0)
        );
        assert_eq!(f.sorted_nodes, vec![1, 2, 3]);
    }
}
