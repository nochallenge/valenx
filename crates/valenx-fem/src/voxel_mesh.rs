//! Voxel → hexahedral FEM meshing (unstructured-3D-meshing roadmap
//! gap).
//!
//! Converts a 3-D binary occupancy grid — a voxelised solid or an
//! image-segmentation label volume — into a conforming **8-node
//! hexahedral** finite-element mesh: one [`ElementType::Hex8`] brick
//! per "on" voxel, with the shared corner nodes deduplicated on the
//! integer lattice.
//!
//! This is the standard *voxel-to-hex* (a.k.a. "dexel" / blocky-hex)
//! approach used by image-based FEA (Sandia's `automesh`,
//! `ScanIP`-style segmentation meshers, micro-CT bone/material
//! pipelines). It is exact and robust — no surface reconstruction, no
//! marching-cubes ambiguity — at the cost of a stair-stepped boundary.
//! For a labelled CT/segmentation grid every solid voxel becomes a
//! well-shaped (in fact, for cubic spacing, *perfect*) hex, so the mesh
//! drops straight into the native Hex8 solvers
//! ([`crate::native_solver::solve_linear_static_mixed`],
//! [`crate::modal_solver::solve_modal_mixed`], …) with no remeshing.
//!
//! # Why in-house (not the `automesh` crate)
//!
//! Sandia's `automesh` is the reference voxel→hex mesher, but it is
//! **GPL-3.0** licensed. valenx ships under permissive MIT / Apache-2.0,
//! so vendoring `automesh` would virally relicense this crate. The core
//! algorithm is small (one hex per voxel + lattice-node dedup), so it is
//! implemented here directly — no new dependency, license-clean, and
//! using the same [`valenx_mesh::Mesh`] type the rest of the workbench
//! consumes.
//!
//! # Node ordering
//!
//! Each voxel `(i,j,k)` emits its eight corners in the **standard hex
//! order** [`crate::elements::Hex8`] expects and that
//! [`crate::meshgen::structured_hex_mesh`] already uses: the `z=k`
//! (bottom) face CCW —
//! `(i,j) (i+1,j) (i+1,j+1) (i,j+1)` — then the `z=k+1` (top) face in the
//! same order. So a fully-occupied grid produces a mesh **identical** to
//! `structured_hex_mesh`, and the generated bricks have positive Jacobian
//! in the Hex8 element.
//!
//! # Coordinates
//!
//! Corner node `(I,J,K)` (a *lattice* corner, `0..=nx` etc.) maps to the
//! physical point `(I·dx, J·dy, K·dz)` from the grid [`spacing`]. An
//! optional origin offset is supported via [`VoxelGrid::origin`].
//!
//! [`spacing`]: VoxelGrid::spacing

use std::collections::HashMap;

use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

/// Maximum number of corner nodes the voxel mesher will allocate.
///
/// A fully-occupied grid has `(nx+1)·(ny+1)·(nz+1)` corner nodes; this
/// cap (matching [`crate::meshgen`]'s `MAX_NODES`) corresponds to a
/// multi-GiB mesh, well beyond any realistic image-based grid. Past it,
/// the implied allocation is in the same class as a hostile recipe
/// trying to OOM the host process, so the build is rejected rather than
/// attempted.
pub const MAX_VOXEL_NODES: usize = 200_000_000;

/// Errors from voxel-grid construction or meshing.
#[derive(Debug, thiserror::Error)]
pub enum VoxelMeshError {
    /// A grid dimension was zero, an `occupancy` length mismatch, a
    /// non-finite / non-positive spacing, or an arithmetic overflow on
    /// the implied node / element counts.
    #[error("invalid voxel grid: {reason}")]
    InvalidGrid {
        /// Human-readable explanation.
        reason: String,
    },
}

/// A 3-D binary occupancy grid — the input to voxel→hex meshing.
///
/// `occupancy` is row-major with **`x` fastest, then `y`, then `z`**:
/// voxel `(i,j,k)` (with `0 ≤ i < nx`, `0 ≤ j < ny`, `0 ≤ k < nz`) is at
/// linear index `i + nx·j + nx·ny·k`. A `true` entry is solid material
/// that becomes one hex element; a `false` entry is void.
///
/// `spacing` is the physical edge length of one voxel along each axis
/// (`dx, dy, dz`); `origin` shifts the whole lattice (defaults to the
/// origin). Together they place corner lattice node `(I,J,K)` at
/// `origin + (I·dx, J·dy, K·dz)`.
#[derive(Clone, Debug)]
pub struct VoxelGrid {
    /// Voxel count along `x`.
    pub nx: usize,
    /// Voxel count along `y`.
    pub ny: usize,
    /// Voxel count along `z`.
    pub nz: usize,
    /// Solid mask, length `nx·ny·nz`, indexed `i + nx·j + nx·ny·k`.
    pub occupancy: Vec<bool>,
    /// Physical voxel edge length per axis `(dx, dy, dz)`.
    pub spacing: [f64; 3],
    /// Physical position of lattice corner `(0,0,0)`.
    pub origin: [f64; 3],
}

impl VoxelGrid {
    /// Build a grid from explicit dimensions, an occupancy mask, and an
    /// (isotropic or anisotropic) spacing, with the lattice rooted at
    /// the origin.
    ///
    /// # Errors
    ///
    /// Returns [`VoxelMeshError::InvalidGrid`] if any dimension is zero,
    /// if `occupancy.len() != nx·ny·nz` (or that product overflows
    /// `usize`), or if any spacing is non-finite or non-positive.
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        occupancy: Vec<bool>,
        spacing: [f64; 3],
    ) -> Result<Self, VoxelMeshError> {
        Self::with_origin(nx, ny, nz, occupancy, spacing, [0.0, 0.0, 0.0])
    }

    /// Like [`VoxelGrid::new`] but with an explicit lattice `origin`.
    ///
    /// # Errors
    ///
    /// See [`VoxelGrid::new`]; additionally the `origin` components must
    /// be finite.
    pub fn with_origin(
        nx: usize,
        ny: usize,
        nz: usize,
        occupancy: Vec<bool>,
        spacing: [f64; 3],
        origin: [f64; 3],
    ) -> Result<Self, VoxelMeshError> {
        if nx == 0 || ny == 0 || nz == 0 {
            return Err(VoxelMeshError::InvalidGrid {
                reason: format!("dimensions must be positive (got nx={nx}, ny={ny}, nz={nz})"),
            });
        }
        let cell_count = nx
            .checked_mul(ny)
            .and_then(|v| v.checked_mul(nz))
            .ok_or_else(|| VoxelMeshError::InvalidGrid {
                reason: format!("cell count nx*ny*nz overflowed usize (nx={nx}, ny={ny}, nz={nz})"),
            })?;
        if occupancy.len() != cell_count {
            return Err(VoxelMeshError::InvalidGrid {
                reason: format!(
                    "occupancy length {} != nx*ny*nz = {cell_count}",
                    occupancy.len()
                ),
            });
        }
        if !spacing.iter().all(|s| s.is_finite()) {
            return Err(VoxelMeshError::InvalidGrid {
                reason: format!("spacing must be finite (got {spacing:?})"),
            });
        }
        if !spacing.iter().all(|&s| s > 0.0) {
            return Err(VoxelMeshError::InvalidGrid {
                reason: format!("spacing must be positive (got {spacing:?})"),
            });
        }
        if !origin.iter().all(|o| o.is_finite()) {
            return Err(VoxelMeshError::InvalidGrid {
                reason: format!("origin must be finite (got {origin:?})"),
            });
        }
        Ok(Self {
            nx,
            ny,
            nz,
            occupancy,
            spacing,
            origin,
        })
    }

    /// Linear index of voxel `(i,j,k)` into [`occupancy`], `x` fastest.
    ///
    /// [`occupancy`]: VoxelGrid::occupancy
    #[inline]
    fn voxel_index(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.nx * j + self.nx * self.ny * k
    }

    /// Whether voxel `(i,j,k)` is solid. Out-of-range coordinates read
    /// as void (`false`).
    #[inline]
    pub fn is_solid(&self, i: usize, j: usize, k: usize) -> bool {
        if i >= self.nx || j >= self.ny || k >= self.nz {
            return false;
        }
        self.occupancy[self.voxel_index(i, j, k)]
    }

    /// Number of solid (occupied) voxels — the hex element count of the
    /// generated mesh.
    pub fn solid_count(&self) -> usize {
        self.occupancy.iter().filter(|&&b| b).count()
    }

    /// Convert this occupancy grid into a hexahedral mesh: one 8-node
    /// hex per solid voxel, corner nodes deduplicated on the lattice.
    ///
    /// Only the corner nodes actually *used* by a solid voxel are
    /// emitted, so a sparse / hollow occupancy does not allocate the
    /// full dense corner lattice. Node ordering per hex matches
    /// [`crate::elements::Hex8`].
    ///
    /// # Errors
    ///
    /// Returns [`VoxelMeshError::InvalidGrid`] if the dense corner
    /// lattice `(nx+1)·(ny+1)·(nz+1)` would overflow `usize` or exceed
    /// [`MAX_VOXEL_NODES`]. (The cap is checked against the dense
    /// lattice, the worst case; the actual emitted node count is ≤ that.)
    pub fn to_hex_mesh(&self) -> Result<HexMesh, VoxelMeshError> {
        // Bound the dense corner lattice even though we emit lazily: a
        // hostile (nx,ny,nz) must not be able to imply an unbounded
        // lattice-key space.
        let nx1 = self.nx + 1; // nx>0 checked at construction, +1 cannot overflow a usize that held nx
        let ny1 = self.ny + 1;
        let nz1 = self.nz + 1;
        let lattice_nodes = nx1
            .checked_mul(ny1)
            .and_then(|v| v.checked_mul(nz1))
            .ok_or_else(|| VoxelMeshError::InvalidGrid {
                reason: "corner lattice (nx+1)*(ny+1)*(nz+1) overflowed usize".to_string(),
            })?;
        if lattice_nodes > MAX_VOXEL_NODES {
            return Err(VoxelMeshError::InvalidGrid {
                reason: format!(
                    "voxel mesh corner lattice would have {lattice_nodes} nodes, \
                     cap is {MAX_VOXEL_NODES}"
                ),
            });
        }

        // Lattice-corner index `(I,J,K)` → emitted node id. Corners are
        // 0..=nx etc.; the key uses the dense stride so equal corners
        // from adjacent voxels collapse to one node.
        let corner_key =
            |ci: usize, cj: usize, ck: usize| -> usize { ci + nx1 * cj + nx1 * ny1 * ck };

        let [dx, dy, dz] = self.spacing;
        let [ox, oy, oz] = self.origin;

        let mut node_of_corner: HashMap<usize, usize> = HashMap::new();
        let mut nodes: Vec<[f64; 3]> = Vec::new();
        let mut elements: Vec<[usize; 8]> = Vec::with_capacity(self.solid_count());

        // The eight corner offsets of a voxel in standard Hex8 order:
        // bottom (z=k) face CCW, then top (z=k+1) face.
        const CORNERS: [(usize, usize, usize); 8] = [
            (0, 0, 0),
            (1, 0, 0),
            (1, 1, 0),
            (0, 1, 0),
            (0, 0, 1),
            (1, 0, 1),
            (1, 1, 1),
            (0, 1, 1),
        ];

        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    if !self.occupancy[self.voxel_index(i, j, k)] {
                        continue;
                    }
                    let mut hex = [0usize; 8];
                    for (slot, &(di, dj, dk)) in CORNERS.iter().enumerate() {
                        let (ci, cj, ck) = (i + di, j + dj, k + dk);
                        let key = corner_key(ci, cj, ck);
                        let node = *node_of_corner.entry(key).or_insert_with(|| {
                            let id = nodes.len();
                            nodes.push([
                                ox + ci as f64 * dx,
                                oy + cj as f64 * dy,
                                oz + ck as f64 * dz,
                            ]);
                            id
                        });
                        hex[slot] = node;
                    }
                    elements.push(hex);
                }
            }
        }

        Ok(HexMesh { nodes, elements })
    }
}

/// A solid hexahedral mesh: deduplicated nodes plus 8-node hex
/// connectivity.
///
/// `nodes[n]` is the `[x,y,z]` coordinate of node `n`; `elements[e]` is
/// the eight node indices of hex `e`, in
/// [`crate::elements::Hex8`] order (bottom face CCW, then top). Convert
/// to the workbench's canonical [`valenx_mesh::Mesh`] (a single
/// [`ElementType::Hex8`] block) with [`HexMesh::to_mesh`] to feed the
/// native solvers.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HexMesh {
    /// Unique node coordinates.
    pub nodes: Vec<[f64; 3]>,
    /// Hex elements as 8 node indices each, standard Hex8 order.
    pub elements: Vec<[usize; 8]>,
}

impl HexMesh {
    /// Number of unique nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of hexahedral elements.
    pub fn element_count(&self) -> usize {
        self.elements.len()
    }

    /// Whether the mesh is empty (no elements — e.g. an all-void grid).
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// Convert to the canonical [`valenx_mesh::Mesh`] — a single
    /// [`ElementType::Hex8`] block — so the mesh drops directly into the
    /// native Hex8 assembly / solvers.
    ///
    /// Node indices are cast to `u32` (the [`ElementBlock`] connectivity
    /// width). Image-based grids are far below `u32::MAX` nodes; the
    /// [`MAX_VOXEL_NODES`] cap (200 M) sits under `u32::MAX` so no node
    /// id can truncate.
    pub fn to_mesh(&self, id: impl Into<String>) -> Mesh {
        let mut mesh = Mesh::new(id);
        mesh.nodes = self
            .nodes
            .iter()
            .map(|&[x, y, z]| Vector3::new(x, y, z))
            .collect();
        let mut conn: Vec<u32> = Vec::with_capacity(self.elements.len() * 8);
        for hex in &self.elements {
            for &n in hex {
                conn.push(n as u32);
            }
        }
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Hex8,
            connectivity: conn,
        });
        mesh.recompute_stats();
        mesh
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elements::{Hex8, SolidElement};

    /// Helper: a fully-solid `nx×ny×nz` grid with unit spacing.
    fn full_grid(nx: usize, ny: usize, nz: usize) -> VoxelGrid {
        VoxelGrid::new(nx, ny, nz, vec![true; nx * ny * nz], [1.0, 1.0, 1.0]).unwrap()
    }

    #[test]
    fn single_voxel_gives_eight_nodes_one_hex() {
        let mesh = full_grid(1, 1, 1).to_hex_mesh().unwrap();
        assert_eq!(
            mesh.node_count(),
            8,
            "single voxel must have 8 corner nodes"
        );
        assert_eq!(mesh.element_count(), 1, "single voxel must be one hex");
    }

    #[test]
    fn two_by_one_by_one_block_gives_twelve_nodes_two_hex_four_shared() {
        // 2×1×1: two cubes sharing one internal face (4 nodes). Total
        // corners = 2·(2·2·2) − 4 shared = 16 − 4 = 12.
        let mesh = full_grid(2, 1, 1).to_hex_mesh().unwrap();
        assert_eq!(mesh.element_count(), 2, "must be two hexes");
        assert_eq!(
            mesh.node_count(),
            12,
            "2x1x1 block must dedup the shared face to 12 nodes"
        );
    }

    #[test]
    fn node_count_matches_unique_corner_count() {
        // For a dense (fully solid) grid, every lattice corner is used,
        // so the node count must equal (nx+1)(ny+1)(nz+1).
        for (nx, ny, nz) in [(1, 1, 1), (2, 1, 1), (3, 2, 2), (4, 4, 1), (2, 3, 4)] {
            let mesh = full_grid(nx, ny, nz).to_hex_mesh().unwrap();
            let expected_nodes = (nx + 1) * (ny + 1) * (nz + 1);
            assert_eq!(
                mesh.node_count(),
                expected_nodes,
                "dense {nx}x{ny}x{nz}: nodes should equal corner lattice"
            );
            assert_eq!(
                mesh.element_count(),
                nx * ny * nz,
                "dense {nx}x{ny}x{nz}: one hex per voxel"
            );
        }
    }

    #[test]
    fn empty_grid_gives_empty_mesh() {
        let grid = VoxelGrid::new(2, 2, 2, vec![false; 8], [1.0, 1.0, 1.0]).unwrap();
        let mesh = grid.to_hex_mesh().unwrap();
        assert!(mesh.is_empty(), "all-void grid must yield no elements");
        assert_eq!(mesh.node_count(), 0, "all-void grid must yield no nodes");
        assert_eq!(mesh.element_count(), 0);
    }

    #[test]
    fn nodes_are_actually_deduplicated_no_coordinate_repeats() {
        // No two emitted nodes may share a coordinate (the lattice dedup
        // is exact because coords come from integer corner indices).
        let mesh = full_grid(3, 2, 2).to_hex_mesh().unwrap();
        let mut seen = std::collections::HashSet::new();
        for p in &mesh.nodes {
            let key = (
                (p[0] * 1e9).round() as i64,
                (p[1] * 1e9).round() as i64,
                (p[2] * 1e9).round() as i64,
            );
            assert!(seen.insert(key), "duplicate node coordinate {p:?}");
        }
        assert_eq!(seen.len(), mesh.node_count());
    }

    #[test]
    fn sparse_grid_only_emits_used_corners() {
        // Two diagonally opposite voxels in a 2×1×1 grid share nothing,
        // so 2 separate cubes → 16 nodes, 2 hexes; the in-between dense
        // lattice corners that no voxel touches are NOT emitted.
        // (Here use a 3×1×1 with the middle voxel void.)
        let grid = VoxelGrid::new(3, 1, 1, vec![true, false, true], [1.0, 1.0, 1.0]).unwrap();
        let mesh = grid.to_hex_mesh().unwrap();
        assert_eq!(mesh.element_count(), 2, "two solid voxels → two hexes");
        // Two disjoint cubes (they do not touch — separated by the void
        // voxel) → 8 + 8 = 16 distinct corner nodes, none shared.
        assert_eq!(
            mesh.node_count(),
            16,
            "disjoint voxels must not share any corner node"
        );
    }

    #[test]
    fn generated_hexes_have_positive_volume_and_correct_total() {
        // Every emitted hex must have positive Jacobian in the Hex8
        // element (correct node ordering), and a dense grid's hex
        // volumes must sum to the box volume.
        let (nx, ny, nz) = (2, 3, 2);
        let (dx, dy, dz) = (0.5, 0.25, 0.75);
        let grid = VoxelGrid::new(nx, ny, nz, vec![true; nx * ny * nz], [dx, dy, dz]).unwrap();
        let mesh = grid.to_hex_mesh().unwrap();
        let mut total = 0.0;
        for hex in &mesh.elements {
            let mut c = [Vector3::zeros(); 8];
            for (slot, &n) in hex.iter().enumerate() {
                let p = mesh.nodes[n];
                c[slot] = Vector3::new(p[0], p[1], p[2]);
            }
            let v = Hex8::new(c).volume().unwrap();
            assert!(v > 0.0, "voxel hex has non-positive volume {v}");
            total += v;
        }
        let expected = (nx as f64 * dx) * (ny as f64 * dy) * (nz as f64 * dz);
        assert!(
            (total - expected).abs() < 1e-9,
            "voxel hex volume sum {total} != box {expected}"
        );
    }

    #[test]
    fn anisotropic_spacing_and_origin_place_nodes_correctly() {
        let grid = VoxelGrid::with_origin(1, 1, 1, vec![true], [2.0, 3.0, 5.0], [10.0, 20.0, 30.0])
            .unwrap();
        let mesh = grid.to_hex_mesh().unwrap();
        // The far corner (1,1,1) sits at origin + spacing.
        let far = [10.0 + 2.0, 20.0 + 3.0, 30.0 + 5.0];
        assert!(
            mesh.nodes.iter().any(|p| (p[0] - far[0]).abs() < 1e-12
                && (p[1] - far[1]).abs() < 1e-12
                && (p[2] - far[2]).abs() < 1e-12),
            "far corner not placed at origin+spacing; nodes={:?}",
            mesh.nodes
        );
        // The near corner sits exactly at the origin.
        assert!(
            mesh.nodes.iter().any(|p| (p[0] - 10.0).abs() < 1e-12
                && (p[1] - 20.0).abs() < 1e-12
                && (p[2] - 30.0).abs() < 1e-12),
            "near corner not at origin"
        );
    }

    #[test]
    fn dense_voxel_mesh_matches_structured_hex_mesh() {
        // A fully-occupied voxel grid must reproduce the canonical
        // structured hex mesh (same node count, same element count).
        let (nx, ny, nz) = (3, 2, 2);
        let voxel = full_grid(nx, ny, nz).to_hex_mesh().unwrap().to_mesh("v");
        let structured =
            crate::meshgen::structured_hex_mesh(nx as f64, ny as f64, nz as f64, nx, ny, nz)
                .unwrap();
        assert_eq!(voxel.nodes.len(), structured.nodes.len());
        assert_eq!(
            voxel.element_blocks[0].connectivity.len(),
            structured.element_blocks[0].connectivity.len()
        );
        assert_eq!(voxel.element_blocks[0].element_type, ElementType::Hex8);
    }

    #[test]
    fn to_mesh_produces_single_hex8_block() {
        let mesh = full_grid(2, 2, 2)
            .to_hex_mesh()
            .unwrap()
            .to_mesh("voxel-test");
        assert_eq!(mesh.element_blocks.len(), 1);
        assert_eq!(mesh.element_blocks[0].element_type, ElementType::Hex8);
        assert_eq!(mesh.element_blocks[0].connectivity.len(), 8 * 8); // 8 voxels
        assert_eq!(mesh.nodes.len(), 27); // 3×3×3 corner lattice
    }

    // --- Validation / hostile-input guards ---

    #[test]
    fn rejects_zero_dimension() {
        let err = VoxelGrid::new(0, 1, 1, vec![], [1.0, 1.0, 1.0]).expect_err("nx=0");
        assert!(format!("{err}").contains("positive"));
    }

    #[test]
    fn rejects_occupancy_length_mismatch() {
        let err = VoxelGrid::new(2, 2, 2, vec![true; 7], [1.0, 1.0, 1.0])
            .expect_err("wrong occupancy length");
        assert!(format!("{err}").contains("occupancy length"));
    }

    #[test]
    fn rejects_non_finite_spacing() {
        let err =
            VoxelGrid::new(1, 1, 1, vec![true], [f64::NAN, 1.0, 1.0]).expect_err("NaN spacing");
        assert!(format!("{err}").contains("finite"));
    }

    #[test]
    fn rejects_non_positive_spacing() {
        let err = VoxelGrid::new(1, 1, 1, vec![true], [0.0, 1.0, 1.0]).expect_err("zero spacing");
        assert!(format!("{err}").contains("positive"));
    }

    #[test]
    fn is_solid_reads_void_outside_bounds() {
        let grid = full_grid(2, 2, 2);
        assert!(grid.is_solid(0, 0, 0));
        assert!(!grid.is_solid(2, 0, 0), "out-of-range must read as void");
        assert!(!grid.is_solid(0, 9, 0));
    }
}
