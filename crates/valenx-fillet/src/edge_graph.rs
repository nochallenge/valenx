//! Edge-graph construction for triangle meshes.
//!
//! Wraps a [`valenx_mesh::Mesh`] with adjacency from each undirected
//! edge `(min_vertex, max_vertex)` to the indices of every [`Tri3`]
//! triangle that touches it. Boundary edges have exactly one
//! triangle; interior manifold edges have two. Non-manifold edges
//! (3+ touching triangles) can occur in degenerate input and are
//! tolerated here but rejected later by the fillet/chamfer ops.
//!
//! [`Tri3`]: valenx_mesh::ElementType::Tri3

use std::collections::HashMap;

use nalgebra::Vector3;
use valenx_mesh::{ElementType, Mesh};

/// Canonical undirected edge identifier: stores the two endpoint
/// vertex indices sorted so that `(a, b)` and `(b, a)` hash equal.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub struct EdgeKey(
    /// Smaller vertex index.
    pub usize,
    /// Larger vertex index.
    pub usize,
);

impl EdgeKey {
    /// Canonicalize: store the smaller index first.
    pub fn new(a: usize, b: usize) -> Self {
        if a <= b {
            EdgeKey(a, b)
        } else {
            EdgeKey(b, a)
        }
    }
}

/// One Tri3 triangle by its three vertex indices, in the original
/// winding order from the mesh's element block. Index `tri_idx` into
/// [`EdgeGraph::triangles`] is the same index used in
/// [`EdgeGraph::adjacency`]'s value vectors.
#[derive(Copy, Clone, Debug)]
pub struct Tri {
    /// Three node indices, in mesh winding order.
    pub v: [usize; 3],
}

/// Adjacency table: for each undirected edge, which triangles touch it.
///
/// Use [`EdgeGraph::from_mesh`] to build, then [`EdgeGraph::edges`] to
/// walk every edge.
#[derive(Clone, Debug, Default)]
pub struct EdgeGraph {
    /// `EdgeKey -> Vec<tri_idx>`. 1-element vec = boundary edge,
    /// 2-element vec = interior manifold edge, 3+ = non-manifold.
    pub adjacency: HashMap<EdgeKey, Vec<usize>>,
    /// Flattened list of every Tri3 triangle in the source mesh, in
    /// the order encountered across all Tri3 element blocks.
    pub triangles: Vec<Tri>,
}

impl EdgeGraph {
    /// True if the edge has exactly one adjacent triangle — i.e. it
    /// lies on the mesh boundary (a hole, or the outer perimeter of
    /// an open shell). Unknown / never-seen edges return `false`.
    pub fn is_boundary(&self, key: EdgeKey) -> bool {
        self.adjacency.get(&key).is_some_and(|v| v.len() == 1)
    }

    /// True if the edge has exactly two adjacent triangles — a
    /// well-formed manifold interior edge. Excludes both boundaries
    /// (1 triangle) and non-manifold edges (3+).
    pub fn is_manifold_interior(&self, key: EdgeKey) -> bool {
        self.adjacency.get(&key).is_some_and(|v| v.len() == 2)
    }

    /// Iterate every edge in the graph in arbitrary order. The
    /// resulting reference is stable for the lifetime of `self`.
    pub fn edges(&self) -> impl Iterator<Item = &EdgeKey> {
        self.adjacency.keys()
    }

    /// Compute the outward face normal of triangle `tri_idx`.
    /// The normal direction follows the right-hand rule on the
    /// triangle's stored winding `(v0, v1, v2)` —
    /// `n = (v1 - v0) cross (v2 - v0)`, normalized.
    ///
    /// Panics if `tri_idx` is out of bounds. Returns the zero vector
    /// for degenerate (collinear) triangles.
    pub fn triangle_normal(&self, mesh: &Mesh, tri_idx: usize) -> Vector3<f64> {
        let v = self.triangles[tri_idx].v;
        let p0 = mesh.nodes[v[0]];
        let p1 = mesh.nodes[v[1]];
        let p2 = mesh.nodes[v[2]];
        let n = (p1 - p0).cross(&(p2 - p0));
        let len = n.norm();
        if len < 1e-30 {
            Vector3::zeros()
        } else {
            n / len
        }
    }

    /// Compute the centroid (arithmetic mean of vertices) of
    /// triangle `tri_idx`.
    pub fn triangle_centroid(&self, mesh: &Mesh, tri_idx: usize) -> Vector3<f64> {
        let v = self.triangles[tri_idx].v;
        (mesh.nodes[v[0]] + mesh.nodes[v[1]] + mesh.nodes[v[2]]) / 3.0
    }

    /// Dihedral angle (radians) at the given interior edge, defined
    /// as the angle between the two adjacent face normals.
    ///
    /// Returns:
    /// - `Some(0.0)` for a perfectly flat edge (coplanar triangles).
    /// - `Some(π/2)` for a 90° fold (e.g. two cube faces sharing an
    ///   edge).
    /// - `None` for boundary or non-manifold edges (≠2 triangles).
    pub fn dihedral_angle(&self, mesh: &Mesh, key: EdgeKey) -> Option<f64> {
        let tris = self.adjacency.get(&key)?;
        if tris.len() != 2 {
            return None;
        }
        let n1 = self.triangle_normal(mesh, tris[0]);
        let n2 = self.triangle_normal(mesh, tris[1]);
        let cos_theta = n1.dot(&n2).clamp(-1.0, 1.0);
        Some(cos_theta.acos())
    }

    /// True if the edge is *convex* — the surface bulges outward
    /// across the edge, like the outside corner of a cube. Concave
    /// edges (inside of an L-bracket) return `false`.
    ///
    /// **Requires consistent outward-pointing normals on the input
    /// mesh.** Meshes produced by `valenx_cad::solid_to_mesh`
    /// satisfy this; hand-built test meshes may need their winding
    /// fixed.
    ///
    /// Heuristic: find tri 1's "opposite vertex" `q2` (the vertex
    /// of tri 1 not on the shared edge) and check which side of
    /// tri 0's plane it sits on. For a **convex** corner, `q2` is
    /// behind tri 0's outward normal (toward the interior bulk),
    /// so `n1 · (q2 - va) < 0`. For a **concave** corner, `q2`
    /// is on the positive normal side.
    ///
    /// Returns `None` for boundary or non-manifold edges.
    pub fn is_convex(&self, mesh: &Mesh, key: EdgeKey) -> Option<bool> {
        let tris = self.adjacency.get(&key)?;
        if tris.len() != 2 {
            return None;
        }
        let n1 = self.triangle_normal(mesh, tris[0]);
        let va = mesh.nodes[key.0];
        // Opposite vertex of tri 1 (the one not on the shared edge).
        let t1 = &self.triangles[tris[1]];
        let q2_idx =
            t1.v.iter()
                .copied()
                .find(|&v| v != key.0 && v != key.1)
                .expect("manifold triangle has 3 distinct vertices");
        let q2 = mesh.nodes[q2_idx];
        Some(n1.dot(&(q2 - va)) < 0.0)
    }

    /// Return every interior manifold edge that is convex **and**
    /// whose dihedral angle exceeds `threshold_rad`.
    ///
    /// `threshold_rad` is the minimum *deviation from flat* — so
    /// `0.0` returns every convex edge, while a typical filleting
    /// threshold like 45° (~0.785 rad) skips nearly-flat seams that
    /// shouldn't be rounded.
    ///
    /// Output is in arbitrary order.
    pub fn filletable_edges(&self, mesh: &Mesh, threshold_rad: f64) -> Vec<EdgeKey> {
        let mut out = Vec::new();
        for &key in self.adjacency.keys() {
            let Some(angle) = self.dihedral_angle(mesh, key) else {
                continue;
            };
            if angle < threshold_rad {
                continue;
            }
            let Some(convex) = self.is_convex(mesh, key) else {
                continue;
            };
            if !convex {
                continue;
            }
            out.push(key);
        }
        out
    }

    /// Walk every Tri3 element of `mesh` (skipping non-triangle
    /// element types) and register all 3 edges of each triangle into
    /// the adjacency map.
    pub fn from_mesh(mesh: &Mesh) -> Self {
        let mut g = EdgeGraph::default();
        for block in &mesh.element_blocks {
            if block.element_type != ElementType::Tri3 {
                continue;
            }
            for tri_nodes in block.connectivity.chunks_exact(3) {
                let v = [
                    tri_nodes[0] as usize,
                    tri_nodes[1] as usize,
                    tri_nodes[2] as usize,
                ];
                let tri_idx = g.triangles.len();
                g.triangles.push(Tri { v });
                for (a, b) in [(v[0], v[1]), (v[1], v[2]), (v[2], v[0])] {
                    g.adjacency
                        .entry(EdgeKey::new(a, b))
                        .or_default()
                        .push(tri_idx);
                }
            }
        }
        g
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use nalgebra::Vector3;
    use valenx_mesh::{ElementBlock, ElementType, Mesh};

    use super::*;

    /// Two right-isoceles triangles sharing the diagonal of a unit
    /// square in the XY plane.
    fn two_tri_quad() -> Mesh {
        let mut m = Mesh::new("quad");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![0, 1, 2, 0, 2, 3];
        m.element_blocks.push(block);
        m
    }

    #[test]
    fn edge_key_canonicalizes() {
        assert_eq!(EdgeKey::new(2, 5), EdgeKey(2, 5));
        assert_eq!(EdgeKey::new(5, 2), EdgeKey(2, 5));
        assert_eq!(EdgeKey::new(7, 7), EdgeKey(7, 7));
    }

    #[test]
    fn quad_has_five_edges_one_interior() {
        let m = two_tri_quad();
        let g = EdgeGraph::from_mesh(&m);
        assert_eq!(g.triangles.len(), 2);
        // 4 boundary edges (the perimeter) + 1 shared diagonal = 5.
        assert_eq!(g.adjacency.len(), 5);
        let interior_count = g.adjacency.values().filter(|v| v.len() == 2).count();
        let boundary_count = g.adjacency.values().filter(|v| v.len() == 1).count();
        assert_eq!(
            interior_count, 1,
            "quad should have 1 interior edge (the diagonal)"
        );
        assert_eq!(boundary_count, 4, "quad should have 4 perimeter edges");
    }

    #[test]
    fn shared_edge_lists_both_triangles() {
        let m = two_tri_quad();
        let g = EdgeGraph::from_mesh(&m);
        let diag = g
            .adjacency
            .get(&EdgeKey::new(0, 2))
            .expect("diagonal exists");
        assert_eq!(diag.len(), 2);
        // The two triangles in winding order are 0 and 1.
        let mut sorted = diag.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1]);
    }

    #[test]
    fn is_boundary_and_is_manifold_interior() {
        let m = two_tri_quad();
        let g = EdgeGraph::from_mesh(&m);
        // Diagonal (0, 2) is interior.
        assert!(g.is_manifold_interior(EdgeKey::new(0, 2)));
        assert!(!g.is_boundary(EdgeKey::new(0, 2)));
        // Perimeter edge (0, 1) is boundary.
        assert!(g.is_boundary(EdgeKey::new(0, 1)));
        assert!(!g.is_manifold_interior(EdgeKey::new(0, 1)));
        // Unknown edge.
        assert!(!g.is_boundary(EdgeKey::new(99, 100)));
        assert!(!g.is_manifold_interior(EdgeKey::new(99, 100)));
    }

    #[test]
    fn edges_iterates_all() {
        let m = two_tri_quad();
        let g = EdgeGraph::from_mesh(&m);
        let count = g.edges().count();
        assert_eq!(count, 5);
    }

    /// Two triangles meeting at a 90-degree fold along the X axis.
    /// Both face normals point toward the +Y+Z side; the dihedral
    /// (angle between the normals) is π/2. With outward-normal
    /// convention this represents a *concave* (interior corner)
    /// fold — see [`ninety_degree_convex_corner`] for the convex
    /// counterpart.
    fn ninety_degree_fold() -> Mesh {
        let mut m = Mesh::new("fold");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0), // 0
            Vector3::new(1.0, 0.0, 0.0), // 1
            Vector3::new(0.0, 1.0, 0.0), // 2 — XY plane
            Vector3::new(0.0, 0.0, 1.0), // 3 — XZ plane
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        // Tri 0: (0, 1, 2) — normal +Z.
        // Tri 1: (0, 3, 1) — normal +Y.
        block.connectivity = vec![0, 1, 2, 0, 3, 1];
        m.element_blocks.push(block);
        m
    }

    /// Two triangles meeting at a 90-degree convex outer corner along
    /// the X axis. Like the corner of the +X+Y+Z octant of a cube
    /// viewed from outside — both face normals point AWAY from the
    /// bulk (toward -Y and -Z respectively).
    fn ninety_degree_convex_corner() -> Mesh {
        let mut m = Mesh::new("convex-corner");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0), // 0
            Vector3::new(1.0, 0.0, 0.0), // 1 — shared edge with 0 along X axis
            Vector3::new(0.0, 1.0, 0.0), // 2 — opposite vertex of tri 0
            Vector3::new(0.0, 0.0, 1.0), // 3 — opposite vertex of tri 1
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        // Tri 0: (0, 2, 1) — winding gives normal -Z (outward, away
        // from the +Z hemisphere of the bulk).
        // Tri 1: (0, 1, 3) — winding gives normal -Y (outward).
        block.connectivity = vec![0, 2, 1, 0, 1, 3];
        m.element_blocks.push(block);
        m
    }

    #[test]
    fn flat_quad_has_zero_dihedral() {
        let m = two_tri_quad();
        let g = EdgeGraph::from_mesh(&m);
        let angle = g.dihedral_angle(&m, EdgeKey::new(0, 2)).unwrap();
        assert!(
            angle.abs() < 1e-9,
            "expected ≈0 for a flat quad, got {angle}"
        );
    }

    #[test]
    fn ninety_degree_fold_has_half_pi_dihedral() {
        let m = ninety_degree_fold();
        let g = EdgeGraph::from_mesh(&m);
        let angle = g.dihedral_angle(&m, EdgeKey::new(0, 1)).unwrap();
        assert!(
            (angle - std::f64::consts::FRAC_PI_2).abs() < 1e-9,
            "expected π/2, got {angle}"
        );
    }

    #[test]
    fn dihedral_returns_none_for_boundary() {
        let m = two_tri_quad();
        let g = EdgeGraph::from_mesh(&m);
        // (0, 1) is a perimeter edge.
        assert!(g.dihedral_angle(&m, EdgeKey::new(0, 1)).is_none());
    }

    #[test]
    fn triangle_normal_unit_length() {
        let m = two_tri_quad();
        let g = EdgeGraph::from_mesh(&m);
        let n = g.triangle_normal(&m, 0);
        assert!((n.norm() - 1.0).abs() < 1e-9);
        // The XY-plane triangle's normal should be ±Z.
        assert!(n.z.abs() > 0.999);
    }

    #[test]
    fn convex_corner_classifies_as_convex() {
        let m = ninety_degree_convex_corner();
        let g = EdgeGraph::from_mesh(&m);
        assert_eq!(g.is_convex(&m, EdgeKey::new(0, 1)), Some(true));
    }

    #[test]
    fn inward_normal_fold_classifies_as_concave() {
        // ninety_degree_fold's normals point INTO the +Y+Z region.
        // Under outward-normal convention this is a concave fold.
        let m = ninety_degree_fold();
        let g = EdgeGraph::from_mesh(&m);
        assert_eq!(g.is_convex(&m, EdgeKey::new(0, 1)), Some(false));
    }

    #[test]
    fn is_convex_returns_none_for_boundary() {
        let m = two_tri_quad();
        let g = EdgeGraph::from_mesh(&m);
        assert!(g.is_convex(&m, EdgeKey::new(0, 1)).is_none());
    }

    /// Construct a unit cube (corners at (0,0,0) and (1,1,1))
    /// tessellated as 2 triangles per face = 12 triangles. All face
    /// normals point outward (away from cube center).
    pub(crate) fn unit_cube() -> Mesh {
        let mut m = Mesh::new("unit-cube");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0), // 0
            Vector3::new(1.0, 0.0, 0.0), // 1
            Vector3::new(1.0, 1.0, 0.0), // 2
            Vector3::new(0.0, 1.0, 0.0), // 3
            Vector3::new(0.0, 0.0, 1.0), // 4
            Vector3::new(1.0, 0.0, 1.0), // 5
            Vector3::new(1.0, 1.0, 1.0), // 6
            Vector3::new(0.0, 1.0, 1.0), // 7
        ];
        let mut block = ElementBlock::new(ElementType::Tri3);
        block.connectivity = vec![
            // bottom (z=0) — normal -Z
            0, 2, 1, 0, 3, 2, // top (z=1) — normal +Z
            4, 5, 6, 4, 6, 7, // front (y=0) — normal -Y
            0, 1, 5, 0, 5, 4, // back (y=1) — normal +Y
            3, 7, 6, 3, 6, 2, // left (x=0) — normal -X
            0, 4, 7, 0, 7, 3, // right (x=1) — normal +X
            1, 2, 6, 1, 6, 5,
        ];
        m.element_blocks.push(block);
        m
    }

    #[test]
    fn unit_cube_has_twelve_filletable_edges() {
        let m = unit_cube();
        let g = EdgeGraph::from_mesh(&m);
        // 12 cube edges (between faces, 90° fold) + 6 face diagonals
        // (between coplanar triangles of the same face, 0° fold).
        assert_eq!(g.adjacency.len(), 12 + 6);
        // Threshold = 45° (~0.785 rad). All 12 cube edges have
        // dihedral π/2 ≈ 1.571 rad, all 6 face diagonals have
        // dihedral 0. So exactly 12 pass.
        let filletable = g.filletable_edges(&m, std::f64::consts::FRAC_PI_4);
        assert_eq!(
            filletable.len(),
            12,
            "expected 12 cube edges to be filletable, got {}",
            filletable.len()
        );
    }

    #[test]
    fn flat_quad_has_no_filletable_edges() {
        // A coplanar two-triangle quad has one interior edge with
        // dihedral 0; no edges should exceed even a 1° threshold.
        let m = two_tri_quad();
        let g = EdgeGraph::from_mesh(&m);
        let filletable = g.filletable_edges(&m, 1.0_f64.to_radians());
        assert!(filletable.is_empty());
    }

    #[test]
    fn concave_fold_is_not_filletable() {
        // Using the inward-normal fold which classifies as concave.
        let m = ninety_degree_fold();
        let g = EdgeGraph::from_mesh(&m);
        // Dihedral is π/2 so it passes the angle threshold, but
        // concave so it's not filletable.
        let filletable = g.filletable_edges(&m, 0.1);
        assert!(filletable.is_empty());
    }

    #[test]
    fn skips_non_tri3_blocks() {
        let mut m = Mesh::new("mixed");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut tet_block = ElementBlock::new(ElementType::Tet4);
        tet_block.connectivity = vec![0, 1, 2, 3];
        m.element_blocks.push(tet_block);
        let g = EdgeGraph::from_mesh(&m);
        assert!(g.triangles.is_empty(), "Tet4 should be skipped");
        assert!(g.adjacency.is_empty());
    }
}
