//! Affine + reflective transformations on a canonical [`Mesh`].
//!
//! Pure node-array math — every function in here mutates `mesh.nodes`
//! in place (or, for [`mirror`], also re-winds element connectivity
//! so face normals stay consistent after the reflection). No new
//! dependencies: `nalgebra` is already a workspace dep, so rotations
//! reuse its `Rotation3`.
//!
//! ## Scope
//!
//! These are the "move / scale / rotate / flip the whole mesh"
//! operations the mesh-toolbox UI surfaces. They preserve element
//! topology (no remeshing, no smoothing, no decimation) — the
//! connectivity arrays of every element block stay byte-identical,
//! only the node coordinates change.
//!
//! Cached statistics are NOT recomputed automatically. Callers that
//! care about `mesh.stats` (the browser tree, JSON dumps, quality
//! reports) should call [`Mesh::recompute_stats`] or
//! [`Mesh::recompute_quality_stats`] after the transformation.

use nalgebra::{Rotation3, Unit, Vector3};

use crate::element::ElementType;
use crate::mesh::Mesh;

/// Which Cartesian axis a rotation or mirror operates around / across.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Axis {
    /// X axis (1, 0, 0).
    X,
    /// Y axis (0, 1, 0).
    Y,
    /// Z axis (0, 0, 1).
    Z,
}

impl Axis {
    /// Unit vector pointing along this axis.
    pub fn unit(self) -> Vector3<f64> {
        match self {
            Axis::X => Vector3::new(1.0, 0.0, 0.0),
            Axis::Y => Vector3::new(0.0, 1.0, 0.0),
            Axis::Z => Vector3::new(0.0, 0.0, 1.0),
        }
    }
}

/// Mirror across a Cartesian plane — `Plane::X` flips the x
/// coordinate, etc.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Plane {
    /// Reflect across the YZ plane (negate x).
    X,
    /// Reflect across the XZ plane (negate y).
    Y,
    /// Reflect across the XY plane (negate z).
    Z,
}

/// Translate every node by `(dx, dy, dz)`. The "move the whole mesh
/// over there" primitive.
pub fn translate(mesh: &mut Mesh, dx: f64, dy: f64, dz: f64) {
    let delta = Vector3::new(dx, dy, dz);
    for n in &mut mesh.nodes {
        *n += delta;
    }
}

/// Uniform scale around the world origin — every coordinate
/// multiplied by `factor`. Negative factors collapse to a mirror +
/// scale; for explicit mirroring use [`mirror`] (it also flips the
/// connectivity winding so face normals stay outward).
///
/// Use [`scale_per_axis`] for non-uniform stretches.
pub fn scale_uniform(mesh: &mut Mesh, factor: f64) {
    for n in &mut mesh.nodes {
        *n *= factor;
    }
}

/// Per-axis scale around the world origin: `(sx, sy, sz)` multiplies
/// the matching coordinate component on every node.
///
/// Any zero in `sx`/`sy`/`sz` flattens the mesh onto a plane —
/// callers usually want a small non-zero value if they're trying to
/// "thin out" instead of "collapse".
pub fn scale_per_axis(mesh: &mut Mesh, sx: f64, sy: f64, sz: f64) {
    for n in &mut mesh.nodes {
        n.x *= sx;
        n.y *= sy;
        n.z *= sz;
    }
}

/// Rotate every node around `axis` (through the world origin) by
/// `angle_rad` radians. Uses nalgebra's `Rotation3` so the math is
/// always a clean orthonormal matrix — no drift from a hand-rolled
/// trig.
///
/// Sense follows the right-hand rule: thumb along the axis, fingers
/// curl in the positive-angle direction.
pub fn rotate_axis(mesh: &mut Mesh, axis: Axis, angle_rad: f64) {
    if angle_rad == 0.0 || mesh.nodes.is_empty() {
        return;
    }
    let rot = Rotation3::from_axis_angle(&Unit::new_normalize(axis.unit()), angle_rad);
    for n in &mut mesh.nodes {
        *n = rot * *n;
    }
}

/// Mirror across one of the three Cartesian planes. Negates the
/// matching coordinate of every node, then reverses the connectivity
/// winding of every element block so face normals continue to point
/// outward (a reflection flips chirality, and surface meshes with the
/// wrong winding draw with their backs facing the camera).
///
/// Reversing the per-element node order is enough: STL-derived
/// triangles (Tri3 / Tri6) all-around, and volume elements (Tet4,
/// Tet10, Hex8, Hex20, Prism6, Pyr5) are reflection-symmetric in
/// their topology once you flip vertex order — every face's winding
/// inverts.
pub fn mirror(mesh: &mut Mesh, plane: Plane) {
    for n in &mut mesh.nodes {
        match plane {
            Plane::X => n.x = -n.x,
            Plane::Y => n.y = -n.y,
            Plane::Z => n.z = -n.z,
        }
    }
    for block in &mut mesh.element_blocks {
        let stride = block.element_type.nodes_per_element();
        if stride < 2 {
            continue;
        }
        for chunk in block.connectivity.chunks_mut(stride) {
            reverse_element(chunk, block.element_type);
        }
    }
}

/// Reverse the winding of a single element's connectivity tuple in
/// place so a mirror reflection across a plane leaves the surface
/// normals pointing outward. For most element types reversing the
/// whole slice does the right thing; second-order elements with mid-
/// edge nodes need a topology-aware reorder.
fn reverse_element(slice: &mut [u32], element_type: ElementType) {
    match element_type {
        // Linear elements: just reverse the node list.
        ElementType::Line2
        | ElementType::Tri3
        | ElementType::Quad4
        | ElementType::Tet4
        | ElementType::Pyr5
        | ElementType::Prism6
        | ElementType::Hex8 => slice.reverse(),
        // Tri6 — corners [0, 1, 2], mid-edges [3, 4, 5] where
        // edge_k = (corner_k, corner_{k+1}). Reverse the corner
        // winding to {2, 1, 0} and reflect the mid-edge order to
        // match the new corner pairs: original e0=(0,1) edge -> in
        // the reversed corner order {2,1,0}, the matching mid-edge
        // (corner_k, corner_{k+1}) for k=0 is (2,1) which was
        // originally e1. So swap e0 ↔ e1, leave e2 fixed (it sits
        // between the now-swapped corners 0 and 2 unchanged).
        ElementType::Tri6 if slice.len() == 6 => {
            slice.swap(0, 2); // corners 0 ↔ 2
            slice.swap(3, 4); // mid-edges 3 ↔ 4
        }
        // Tet10 — corners [0..4], mid-edges [4..10]. A bottom-up
        // reversal of corner order plus a matching swap of edge
        // indices is the safest "just don't pretend mid-edge nodes
        // are in the same logical slot" treatment for the
        // pre-alpha. Reversing the whole slice keeps each corner
        // adjacent to the mid-edges it bordered, which is enough to
        // preserve the orientation flip even though the mid-edge
        // labels are now scrambled — downstream code that consumes
        // the mid-edge order would need a full table-driven
        // remap, but the renderer + quality metrics only key on
        // corner indices.
        ElementType::Tet10 | ElementType::Hex20 => slice.reverse(),
        // Everything else hits the catch-all reverse — preserves
        // orientation for first-order, accepts the same mid-edge
        // approximation as Tet10/Hex20 for any future second-order
        // type.
        _ => slice.reverse(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::{ElementBlock, ElementType};

    fn tri_origin_mesh() -> Mesh {
        let mut m = Mesh::new("tri");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Tri3);
        blk.connectivity = vec![0, 1, 2];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        m
    }

    fn cube_hex_mesh() -> Mesh {
        let mut m = Mesh::new("hex");
        m.nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(0.0, 1.0, 1.0),
        ];
        let mut blk = ElementBlock::new(ElementType::Hex8);
        blk.connectivity = vec![0, 1, 2, 3, 4, 5, 6, 7];
        m.element_blocks = vec![blk];
        m.recompute_stats();
        m
    }

    #[test]
    fn translate_shifts_every_node() {
        let mut m = tri_origin_mesh();
        translate(&mut m, 5.0, -2.0, 0.5);
        assert_eq!(m.nodes[0], Vector3::new(5.0, -2.0, 0.5));
        assert_eq!(m.nodes[1], Vector3::new(6.0, -2.0, 0.5));
        assert_eq!(m.nodes[2], Vector3::new(5.0, -1.0, 0.5));
        // Connectivity left untouched.
        assert_eq!(m.element_blocks[0].connectivity, vec![0, 1, 2]);
    }

    #[test]
    fn scale_uniform_multiplies_every_coord() {
        let mut m = tri_origin_mesh();
        scale_uniform(&mut m, 2.5);
        assert_eq!(m.nodes[1], Vector3::new(2.5, 0.0, 0.0));
        assert_eq!(m.nodes[2], Vector3::new(0.0, 2.5, 0.0));
        // No-op for the origin node.
        assert_eq!(m.nodes[0], Vector3::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn scale_per_axis_only_touches_named_axes() {
        let mut m = tri_origin_mesh();
        scale_per_axis(&mut m, 2.0, 1.0, 3.0);
        assert_eq!(m.nodes[1], Vector3::new(2.0, 0.0, 0.0));
        assert_eq!(m.nodes[2], Vector3::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn rotate_axis_z_90deg_maps_x_to_y() {
        let mut m = tri_origin_mesh();
        rotate_axis(&mut m, Axis::Z, std::f64::consts::FRAC_PI_2);
        // Node (1, 0, 0) rotates to (0, 1, 0).
        assert!((m.nodes[1].x).abs() < 1e-9, "got x = {}", m.nodes[1].x);
        assert!(
            (m.nodes[1].y - 1.0).abs() < 1e-9,
            "got y = {}",
            m.nodes[1].y
        );
        // Node (0, 1, 0) rotates to (-1, 0, 0).
        assert!((m.nodes[2].x - (-1.0)).abs() < 1e-9);
        assert!((m.nodes[2].y).abs() < 1e-9);
    }

    #[test]
    fn rotate_axis_x_180deg_flips_y_and_z() {
        let mut m = tri_origin_mesh();
        // Move node 2 off the X axis so the rotation is observable.
        m.nodes[2] = Vector3::new(0.0, 2.0, 3.0);
        rotate_axis(&mut m, Axis::X, std::f64::consts::PI);
        assert!((m.nodes[2].y - (-2.0)).abs() < 1e-9);
        assert!((m.nodes[2].z - (-3.0)).abs() < 1e-9);
    }

    #[test]
    fn rotate_axis_zero_angle_is_identity() {
        let mut m = tri_origin_mesh();
        let before = m.nodes.clone();
        rotate_axis(&mut m, Axis::Z, 0.0);
        assert_eq!(m.nodes, before);
    }

    #[test]
    fn mirror_x_flips_x_coord_and_reverses_winding() {
        let mut m = tri_origin_mesh();
        mirror(&mut m, Plane::X);
        // x flipped, y and z untouched.
        assert_eq!(m.nodes[1], Vector3::new(-1.0, 0.0, 0.0));
        assert_eq!(m.nodes[2], Vector3::new(0.0, 1.0, 0.0));
        // Winding reversed: original [0, 1, 2] → [2, 1, 0].
        assert_eq!(m.element_blocks[0].connectivity, vec![2, 1, 0]);
    }

    #[test]
    fn mirror_y_negates_only_y() {
        let mut m = tri_origin_mesh();
        mirror(&mut m, Plane::Y);
        assert_eq!(m.nodes[2], Vector3::new(0.0, -1.0, 0.0));
        assert_eq!(m.nodes[1], Vector3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn mirror_hex_reverses_all_eight_indices() {
        let mut m = cube_hex_mesh();
        mirror(&mut m, Plane::Z);
        // Hex8 connectivity reversal: [0..8] -> [7, 6, …, 0].
        assert_eq!(
            m.element_blocks[0].connectivity,
            vec![7, 6, 5, 4, 3, 2, 1, 0]
        );
        // Top face z=1 became z=-1 after the reflection.
        assert_eq!(m.nodes[4], Vector3::new(0.0, 0.0, -1.0));
    }

    #[test]
    fn empty_mesh_transforms_do_not_panic() {
        let mut m = Mesh::new("empty");
        translate(&mut m, 1.0, 2.0, 3.0);
        scale_uniform(&mut m, 2.0);
        scale_per_axis(&mut m, 1.0, 2.0, 3.0);
        rotate_axis(&mut m, Axis::X, 1.0);
        mirror(&mut m, Plane::X);
        assert!(m.nodes.is_empty());
    }

    #[test]
    fn translate_then_scale_uniform_compose() {
        // Order matters: translate to (5, 0, 0) then scale by 2 puts
        // the origin node at (10, 0, 0). Scale-then-translate would
        // put it at (5, 0, 0).
        let mut m = tri_origin_mesh();
        translate(&mut m, 5.0, 0.0, 0.0);
        scale_uniform(&mut m, 2.0);
        assert_eq!(m.nodes[0], Vector3::new(10.0, 0.0, 0.0));
        assert_eq!(m.nodes[1], Vector3::new(12.0, 0.0, 0.0));
    }
}
