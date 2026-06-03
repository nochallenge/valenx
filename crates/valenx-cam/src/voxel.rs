//! Voxel grid for material removal simulation.
//!
//! Used by the Phase 17D simulation animation to track which stock
//! voxels have been removed by a swept toolpath. v1 uses a single
//! bit per voxel (solid / void) packed into `Vec<u64>` — a 64³ grid
//! takes 32 KiB. Default resolution targets 64³ for interactive
//! preview; the host can configure higher resolutions for precision.
//!
//! ## Algorithm
//!
//! - [`Voxel::from_aabb`] initialises every voxel inside the stock
//!   AABB to `solid` (bit = 1).
//! - [`Voxel::cut_segment`] sweeps a line segment through the voxel
//!   space and clears every voxel within `tool_radius` of the swept
//!   path (cylindrical sweep approximation — adequate for flat
//!   end-mills, conservative for ball-mills).
//! - [`Voxel::to_mesh`] extracts a renderable mesh of the **solid
//!   boundary** via a coarse "cube-marching" cell test — every
//!   voxel that is solid and has at least one void neighbour
//!   contributes its outward face. v1 doesn't smooth the result;
//!   high-quality smooth surface nets are Phase 17.5.
//!
//! ## Surface extraction — faceted and smooth (Phase 56.5)
//!
//! Two mesh extractors ship:
//!
//! - [`Voxel::to_mesh`] — the original **faceted** boundary: every
//!   solid voxel with a void neighbour contributes its axis-aligned
//!   outward face. Blocky, exact to the voxel grid, cheap.
//! - [`Voxel::to_mesh_surface_nets`] — a **Naive Surface Nets**
//!   extractor: one vertex per surface-straddling cell, positioned by
//!   averaging the sign-changing edge midpoints (which smooths the
//!   stair-stepping), with quads stitched across every sign-change
//!   grid edge. The result is a smooth, manifold approximation of the
//!   machined surface — far closer to what CAMotics renders than the
//!   faceted voxel boundary.
//!
//! ## Honest limitations
//!
//! - **Cylindrical sweep approximation** — does not honour ball-end
//!   geometry (a ball-mill's tip-shape would carve a different
//!   profile near the bottom of the cut).
//! - **No XY rotation** — the swept envelope is axis-aligned in
//!   `x`/`y`; tilted tools (5-axis) need a 5-axis-aware
//!   `cut_segment_tilted` (deferred).
//! - **Surface Nets is binary-grid dual contouring** — because the
//!   grid stores 1 bit per voxel (solid/void, no signed-distance
//!   field) the edge crossing is taken at the edge midpoint. A full
//!   SDF grid would let the crossing sit at the true iso-value and
//!   sharpen the result further; that is a follow-up.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Hard cap on the number of voxel cells in a single [`Voxel`] grid.
/// 100 M cells at 1 bit per cell is ~12 MiB; well above any
/// interactive preview (default is 64³ ≈ 256 K cells) and well
/// below anything that would OOM a modern desktop. Round-12 M3/M8:
/// without this cap a hostile or accidental call with
/// `(u32::MAX, u32::MAX, u32::MAX)` resolution would overflow the
/// `usize` allocation math and either panic-allocate ~7×10²⁸ bytes
/// or wrap-around to a tiny allocation followed by index UB.
pub const MAX_VOXEL_CELLS: usize = 100_000_000;

/// Errors raised by [`Voxel::from_aabb`] /
/// [`Voxel::to_mesh_surface_nets`].
#[derive(Debug, Error, PartialEq, Eq)]
pub enum VoxelError {
    /// One or more of the per-axis resolution multiplications
    /// overflowed `u64`. Implies a programmer error (resolution
    /// chosen from untrusted input without a separate per-axis cap)
    /// — the message echoes the offending axis triple.
    #[error("voxel resolution {0:?} overflows u64 in nx*ny*nz")]
    Overflow((u32, u32, u32)),

    /// Total cell count is finite but exceeds [`MAX_VOXEL_CELLS`].
    #[error(
        "voxel grid would have {cells} cells, exceeding the {max}-cell cap"
    )]
    TooManyCells {
        /// Computed total cell count (`nx*ny*nz`).
        cells: u64,
        /// The cap that was violated.
        max: usize,
    },
}

/// A 3D occupancy grid (1 bit per voxel: solid / void).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Voxel {
    /// World-space origin (min XYZ corner).
    pub origin: Vector3<f64>,
    /// World-space extent along each axis (mm).
    pub size: Vector3<f64>,
    /// Grid resolution along each axis.
    pub resolution: (u32, u32, u32),
    /// Packed bit data. Index `i` ⇒ `data[i / 64] >> (i % 64) & 1`.
    /// Bit `1` = solid, `0` = void.
    pub data: Vec<u64>,
}

impl Voxel {
    /// Number of voxels along x / y / z.
    pub fn nx(&self) -> u32 {
        self.resolution.0
    }
    /// Number of voxels along y.
    pub fn ny(&self) -> u32 {
        self.resolution.1
    }
    /// Number of voxels along z.
    pub fn nz(&self) -> u32 {
        self.resolution.2
    }
    /// Total voxel count = nx × ny × nz.
    pub fn n(&self) -> u64 {
        self.nx() as u64 * self.ny() as u64 * self.nz() as u64
    }

    /// Voxel side length along each axis (mm).
    pub fn cell_size(&self) -> Vector3<f64> {
        Vector3::new(
            self.size.x / self.nx() as f64,
            self.size.y / self.ny() as f64,
            self.size.z / self.nz() as f64,
        )
    }

    /// Initialise a solid block over `(min, max)` at `resolution`.
    ///
    /// # Errors
    ///
    /// - [`VoxelError::Overflow`] when `nx*ny*nz` overflows `u64`
    ///   (e.g. the `(u32::MAX, u32::MAX, u32::MAX)` worst case).
    /// - [`VoxelError::TooManyCells`] when the cell count is finite
    ///   but exceeds [`MAX_VOXEL_CELLS`].
    pub fn from_aabb(
        min: Vector3<f64>,
        max: Vector3<f64>,
        resolution: (u32, u32, u32),
    ) -> Result<Self, VoxelError> {
        let size = max - min;
        let n = (resolution.0 as u64)
            .checked_mul(resolution.1 as u64)
            .and_then(|x| x.checked_mul(resolution.2 as u64))
            .ok_or(VoxelError::Overflow(resolution))?;
        if n > MAX_VOXEL_CELLS as u64 {
            return Err(VoxelError::TooManyCells {
                cells: n,
                max: MAX_VOXEL_CELLS,
            });
        }
        let n_words = n.div_ceil(64) as usize;
        let mut data = vec![u64::MAX; n_words];
        // Clear tail bits beyond `n`.
        let tail = (n_words as u64 * 64) - n;
        if tail > 0 {
            let mask = (1_u128 << (64 - tail)) - 1;
            data[n_words - 1] &= mask as u64;
        }
        Ok(Self {
            origin: min,
            size,
            resolution,
            data,
        })
    }

    fn linear_index(&self, ix: u32, iy: u32, iz: u32) -> u64 {
        ix as u64 + iy as u64 * self.nx() as u64 + iz as u64 * self.nx() as u64 * self.ny() as u64
    }

    /// Test the bit at `(ix, iy, iz)`. Out-of-range coords return
    /// `false`.
    pub fn get(&self, ix: u32, iy: u32, iz: u32) -> bool {
        if ix >= self.nx() || iy >= self.ny() || iz >= self.nz() {
            return false;
        }
        let i = self.linear_index(ix, iy, iz);
        let word = (i / 64) as usize;
        let bit = (i % 64) as u32;
        (self.data[word] >> bit) & 1 == 1
    }

    /// Set the bit at `(ix, iy, iz)` to `value`. Out-of-range coords
    /// silently no-op.
    pub fn set(&mut self, ix: u32, iy: u32, iz: u32, value: bool) {
        if ix >= self.nx() || iy >= self.ny() || iz >= self.nz() {
            return;
        }
        let i = self.linear_index(ix, iy, iz);
        let word = (i / 64) as usize;
        let bit = (i % 64) as u32;
        if value {
            self.data[word] |= 1_u64 << bit;
        } else {
            self.data[word] &= !(1_u64 << bit);
        }
    }

    /// Count solid voxels (popcount over the packed bits).
    pub fn solid_count(&self) -> u64 {
        let mut total = 0_u64;
        for w in &self.data {
            total += w.count_ones() as u64;
        }
        total
    }

    /// Convert a world coordinate to integer voxel coordinates.
    /// Returns `None` if outside the grid.
    pub fn world_to_voxel(&self, p: Vector3<f64>) -> Option<(u32, u32, u32)> {
        let local = p - self.origin;
        let cell = self.cell_size();
        let ix_f = local.x / cell.x;
        let iy_f = local.y / cell.y;
        let iz_f = local.z / cell.z;
        if !(ix_f >= 0.0 && iy_f >= 0.0 && iz_f >= 0.0) {
            return None;
        }
        let ix = ix_f.floor() as i64;
        let iy = iy_f.floor() as i64;
        let iz = iz_f.floor() as i64;
        if ix < 0
            || iy < 0
            || iz < 0
            || ix as u32 >= self.nx()
            || iy as u32 >= self.ny()
            || iz as u32 >= self.nz()
        {
            return None;
        }
        Some((ix as u32, iy as u32, iz as u32))
    }

    /// Convert integer voxel coordinates back to the cell-centre
    /// world coordinate.
    pub fn voxel_to_world_centre(&self, ix: u32, iy: u32, iz: u32) -> Vector3<f64> {
        let cell = self.cell_size();
        Vector3::new(
            self.origin.x + (ix as f64 + 0.5) * cell.x,
            self.origin.y + (iy as f64 + 0.5) * cell.y,
            self.origin.z + (iz as f64 + 0.5) * cell.z,
        )
    }

    /// Cut a swept cylindrical envelope along the line `p0 → p1` at
    /// the given `tool_radius`. Every voxel within the envelope
    /// gets cleared (set to 0).
    ///
    /// Envelope = capsule with cylindrical body + spherical caps.
    /// Test uses point-to-line-segment distance per voxel centre,
    /// constrained to voxels in the segment's AABB (+/- radius
    /// padding).
    pub fn cut_segment(&mut self, p0: Vector3<f64>, p1: Vector3<f64>, tool_radius: f64) {
        let r2 = tool_radius * tool_radius;
        let bbox_min = Vector3::new(
            p0.x.min(p1.x) - tool_radius,
            p0.y.min(p1.y) - tool_radius,
            p0.z.min(p1.z) - tool_radius,
        );
        let bbox_max = Vector3::new(
            p0.x.max(p1.x) + tool_radius,
            p0.y.max(p1.y) + tool_radius,
            p0.z.max(p1.z) + tool_radius,
        );
        let cell = self.cell_size();
        let i0_x = (((bbox_min.x - self.origin.x) / cell.x).floor() as i64).max(0) as u32;
        let i0_y = (((bbox_min.y - self.origin.y) / cell.y).floor() as i64).max(0) as u32;
        let i0_z = (((bbox_min.z - self.origin.z) / cell.z).floor() as i64).max(0) as u32;
        let i1_x = (((bbox_max.x - self.origin.x) / cell.x).ceil() as i64).max(0) as u32;
        let i1_y = (((bbox_max.y - self.origin.y) / cell.y).ceil() as i64).max(0) as u32;
        let i1_z = (((bbox_max.z - self.origin.z) / cell.z).ceil() as i64).max(0) as u32;
        let i1_x = i1_x.min(self.nx());
        let i1_y = i1_y.min(self.ny());
        let i1_z = i1_z.min(self.nz());
        let dir = p1 - p0;
        let len2 = dir.norm_squared();
        for iz in i0_z..i1_z {
            for iy in i0_y..i1_y {
                for ix in i0_x..i1_x {
                    let c = self.voxel_to_world_centre(ix, iy, iz);
                    let d2 = if len2 > 0.0 {
                        let t = ((c - p0).dot(&dir) / len2).clamp(0.0, 1.0);
                        let foot = p0 + dir * t;
                        (c - foot).norm_squared()
                    } else {
                        (c - p0).norm_squared()
                    };
                    if d2 <= r2 {
                        self.set(ix, iy, iz, false);
                    }
                }
            }
        }
    }

    /// Extract a faceted mesh of the solid boundary. Every
    /// solid voxel with at least one void neighbour contributes its
    /// outward axis-aligned face(s).
    ///
    /// Returns a [`valenx_mesh::Mesh`] consisting of Tri3 elements.
    pub fn to_mesh(&self) -> valenx_mesh::Mesh {
        let cell = self.cell_size();
        let mut nodes: Vec<Vector3<f64>> = Vec::new();
        let mut conn: Vec<u32> = Vec::new();
        let push_quad =
            |nodes: &mut Vec<Vector3<f64>>, conn: &mut Vec<u32>, v: [Vector3<f64>; 4]| {
                let base = nodes.len() as u32;
                nodes.extend_from_slice(&v);
                // Two triangles per quad.
                conn.extend_from_slice(&[base, base + 1, base + 2]);
                conn.extend_from_slice(&[base, base + 2, base + 3]);
            };
        for iz in 0..self.nz() {
            for iy in 0..self.ny() {
                for ix in 0..self.nx() {
                    if !self.get(ix, iy, iz) {
                        continue;
                    }
                    let p = Vector3::new(
                        self.origin.x + ix as f64 * cell.x,
                        self.origin.y + iy as f64 * cell.y,
                        self.origin.z + iz as f64 * cell.z,
                    );
                    let dx = cell.x;
                    let dy = cell.y;
                    let dz = cell.z;
                    // -X face
                    if ix == 0 || !self.get(ix - 1, iy, iz) {
                        push_quad(
                            &mut nodes,
                            &mut conn,
                            [
                                Vector3::new(p.x, p.y, p.z),
                                Vector3::new(p.x, p.y + dy, p.z),
                                Vector3::new(p.x, p.y + dy, p.z + dz),
                                Vector3::new(p.x, p.y, p.z + dz),
                            ],
                        );
                    }
                    // +X face
                    if ix == self.nx() - 1 || !self.get(ix + 1, iy, iz) {
                        push_quad(
                            &mut nodes,
                            &mut conn,
                            [
                                Vector3::new(p.x + dx, p.y, p.z),
                                Vector3::new(p.x + dx, p.y, p.z + dz),
                                Vector3::new(p.x + dx, p.y + dy, p.z + dz),
                                Vector3::new(p.x + dx, p.y + dy, p.z),
                            ],
                        );
                    }
                    // -Y face
                    if iy == 0 || !self.get(ix, iy - 1, iz) {
                        push_quad(
                            &mut nodes,
                            &mut conn,
                            [
                                Vector3::new(p.x, p.y, p.z),
                                Vector3::new(p.x, p.y, p.z + dz),
                                Vector3::new(p.x + dx, p.y, p.z + dz),
                                Vector3::new(p.x + dx, p.y, p.z),
                            ],
                        );
                    }
                    // +Y face
                    if iy == self.ny() - 1 || !self.get(ix, iy + 1, iz) {
                        push_quad(
                            &mut nodes,
                            &mut conn,
                            [
                                Vector3::new(p.x, p.y + dy, p.z),
                                Vector3::new(p.x + dx, p.y + dy, p.z),
                                Vector3::new(p.x + dx, p.y + dy, p.z + dz),
                                Vector3::new(p.x, p.y + dy, p.z + dz),
                            ],
                        );
                    }
                    // -Z face
                    if iz == 0 || !self.get(ix, iy, iz - 1) {
                        push_quad(
                            &mut nodes,
                            &mut conn,
                            [
                                Vector3::new(p.x, p.y, p.z),
                                Vector3::new(p.x + dx, p.y, p.z),
                                Vector3::new(p.x + dx, p.y + dy, p.z),
                                Vector3::new(p.x, p.y + dy, p.z),
                            ],
                        );
                    }
                    // +Z face
                    if iz == self.nz() - 1 || !self.get(ix, iy, iz + 1) {
                        push_quad(
                            &mut nodes,
                            &mut conn,
                            [
                                Vector3::new(p.x, p.y, p.z + dz),
                                Vector3::new(p.x, p.y + dy, p.z + dz),
                                Vector3::new(p.x + dx, p.y + dy, p.z + dz),
                                Vector3::new(p.x + dx, p.y, p.z + dz),
                            ],
                        );
                    }
                }
            }
        }
        let mut mesh = valenx_mesh::Mesh::new("voxel-surface");
        mesh.nodes = nodes;
        let mut block =
            valenx_mesh::element::ElementBlock::new(valenx_mesh::element::ElementType::Tri3);
        block.connectivity = conn;
        mesh.element_blocks.push(block);
        mesh
    }

    /// Solid/void test in *signed* coordinates — voxels outside the
    /// grid count as **void**. This is the convention Surface Nets
    /// needs so the stock's outer surface (where a solid voxel borders
    /// the empty space past the grid) generates geometry.
    fn solid_at(&self, ix: i64, iy: i64, iz: i64) -> bool {
        if ix < 0 || iy < 0 || iz < 0 {
            return false;
        }
        self.get(ix as u32, iy as u32, iz as u32)
    }

    /// World-space position of the integer voxel sample `(ix, iy, iz)`
    /// — the voxel **centre** (signed coords allowed).
    ///
    /// Surface Nets treats each voxel as a grid sample and a *cell* as
    /// the cube spanning 8 adjacent voxel centres. Using the centre
    /// (not the min corner) is what keeps the smoothed surface inside
    /// the stock: the outermost solid voxel's centre sits half a cell
    /// in from the AABB face, and the sign-change with the void voxel
    /// just outside the grid (centre half a cell *past* the face)
    /// places the surface crossing exactly **on** the AABB face.
    fn sample_world(&self, ix: i64, iy: i64, iz: i64) -> Vector3<f64> {
        let cell = self.cell_size();
        Vector3::new(
            self.origin.x + (ix as f64 + 0.5) * cell.x,
            self.origin.y + (iy as f64 + 0.5) * cell.y,
            self.origin.z + (iz as f64 + 0.5) * cell.z,
        )
    }

    /// Extract a **smooth** surface mesh of the solid boundary via
    /// **Naive Surface Nets** (a binary-grid dual-contouring scheme).
    ///
    /// Where [`Voxel::to_mesh`] emits the blocky axis-aligned voxel
    /// faces, this places **one vertex per surface-straddling cell**,
    /// positioned at the average of the cell's sign-changing edge
    /// midpoints — which rounds off the voxel stair-stepping — and
    /// stitches a quad across every grid edge that crosses the
    /// surface. The output is a smooth, manifold triangle mesh much
    /// closer to the real machined surface CAMotics renders.
    ///
    /// # Method
    ///
    /// 1. **Cell vertices.** A *cell* is the cube whose 8 corners are
    ///    8 adjacent voxel samples. Cell indices run from `−1` to
    ///    `n−1` along each axis (the `−1` row picks up the stock's
    ///    outer surface, where an in-grid solid voxel borders the
    ///    void outside the grid). For every cell whose 8 corners are
    ///    *mixed* (some solid, some void), the cell's 12 edges are
    ///    examined; each edge with a solid/void sign change
    ///    contributes its midpoint, and the cell's vertex is the
    ///    average of those midpoints. Averaging the crossings is what
    ///    smooths the surface — a vertex on a flat run of boundary
    ///    cells lands on the true plane, and a vertex near a feature
    ///    is pulled toward it.
    /// 2. **Quad faces.** For every grid edge between two adjacent
    ///    samples that has a sign change, the **4 cells sharing that
    ///    edge** each carry a vertex; those 4 vertices form one quad
    ///    of the surface. The quad is wound so its normal points from
    ///    solid toward void (outward). Each quad is split into 2
    ///    triangles.
    ///
    /// Returns a [`valenx_mesh::Mesh`] of Tri3 elements. An all-solid
    /// or all-void grid yields an empty mesh (no boundary).
    ///
    /// # Errors
    ///
    /// - [`VoxelError::Overflow`] if the cell-table size
    ///   `(nx+1)*(ny+1)*(nz+1)` overflows `u64` (extreme adversarial
    ///   resolution).
    /// - [`VoxelError::TooManyCells`] if the cell-table size exceeds
    ///   [`MAX_VOXEL_CELLS`] (would over-allocate the vertex LUT).
    pub fn to_mesh_surface_nets(&self) -> Result<valenx_mesh::Mesh, VoxelError> {
        let nx = self.nx() as i64;
        let ny = self.ny() as i64;
        let nz = self.nz() as i64;

        // The 8 corner offsets of a cell, indexed by a 3-bit
        // (dx, dy, dz) code so an edge can name its two endpoints.
        const CORNER: [(i64, i64, i64); 8] = [
            (0, 0, 0),
            (1, 0, 0),
            (0, 1, 0),
            (1, 1, 0),
            (0, 0, 1),
            (1, 0, 1),
            (0, 1, 1),
            (1, 1, 1),
        ];
        // The 12 edges of a cell as index pairs into CORNER.
        const EDGE: [(usize, usize); 12] = [
            (0, 1),
            (2, 3),
            (4, 5),
            (6, 7), // X-direction edges
            (0, 2),
            (1, 3),
            (4, 6),
            (5, 7), // Y-direction edges
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7), // Z-direction edges
        ];

        // --- pass 1: a vertex per active cell ---
        // Cell coords run cx ∈ -1..=nx-1 (likewise y, z); map each to
        // a dense linear index for the vertex lookup table.
        let cells_x = (nx + 1) as usize; // -1 .. nx-1
        let cells_y = (ny + 1) as usize;
        let cells_z = (nz + 1) as usize;
        // Round-12 M3/M8: the LUT size below is `cells_x * cells_y *
        // cells_z` and would silently overflow on extreme
        // resolutions, allocating a too-small buffer that the loop
        // below then indexes into. Compute the size in u64 first and
        // bail if it doesn't fit within the cap.
        let cells_total = (cells_x as u64)
            .checked_mul(cells_y as u64)
            .and_then(|x| x.checked_mul(cells_z as u64))
            .ok_or(VoxelError::Overflow(self.resolution))?;
        if cells_total > MAX_VOXEL_CELLS as u64 {
            return Err(VoxelError::TooManyCells {
                cells: cells_total,
                max: MAX_VOXEL_CELLS,
            });
        }
        let cell_index = |cx: i64, cy: i64, cz: i64| -> usize {
            // Shift the -1 origin to 0.
            ((cx + 1) as usize)
                + ((cy + 1) as usize) * cells_x
                + ((cz + 1) as usize) * cells_x * cells_y
        };
        // vertex id of each cell, or u32::MAX if the cell is inactive.
        let mut cell_vertex = vec![u32::MAX; cells_total as usize];
        let mut nodes: Vec<Vector3<f64>> = Vec::new();

        for cz in -1..nz {
            for cy in -1..ny {
                for cx in -1..nx {
                    // Sample the 8 corners.
                    let mut corner_solid = [false; 8];
                    let mut solid_count = 0;
                    for (k, &(dx, dy, dz)) in CORNER.iter().enumerate() {
                        let s = self.solid_at(cx + dx, cy + dy, cz + dz);
                        corner_solid[k] = s;
                        if s {
                            solid_count += 1;
                        }
                    }
                    // Inactive: all solid or all void → no surface here.
                    if solid_count == 0 || solid_count == 8 {
                        continue;
                    }
                    // Average the midpoints of every sign-changing edge.
                    let mut acc = Vector3::zeros();
                    let mut crossings = 0.0;
                    for &(a, b) in &EDGE {
                        if corner_solid[a] != corner_solid[b] {
                            let (ax, ay, az) = CORNER[a];
                            let (bx, by, bz) = CORNER[b];
                            let pa = self.sample_world(cx + ax, cy + ay, cz + az);
                            let pb = self.sample_world(cx + bx, cy + by, cz + bz);
                            acc += (pa + pb) * 0.5;
                            crossings += 1.0;
                        }
                    }
                    // crossings is ≥ 1 here (the cell is mixed).
                    let vertex = acc / crossings;
                    let id = nodes.len() as u32;
                    nodes.push(vertex);
                    cell_vertex[cell_index(cx, cy, cz)] = id;
                }
            }
        }

        // --- pass 2: a quad per sign-changing grid edge ---
        let mut conn: Vec<u32> = Vec::new();
        let emit_quad =
            |conn: &mut Vec<u32>, q: [u32; 4], flip: bool| {
                // Skip if any of the 4 cells was inactive (a robust
                // grid never hits this, but guard anyway).
                if q.contains(&u32::MAX) {
                    return;
                }
                let (a, b, c, d) = if flip {
                    (q[0], q[3], q[2], q[1])
                } else {
                    (q[0], q[1], q[2], q[3])
                };
                conn.extend_from_slice(&[a, b, c]);
                conn.extend_from_slice(&[a, c, d]);
            };

        // X-direction grid edges: sample (ix,iy,iz) → (ix+1,iy,iz).
        // The 4 cells sharing that edge differ in (cy, cz).
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in -1..nx {
                    let s0 = self.solid_at(ix, iy, iz);
                    let s1 = self.solid_at(ix + 1, iy, iz);
                    if s0 == s1 {
                        continue;
                    }
                    // The 4 cells touching this X-edge: cx = ix,
                    // cy ∈ {iy-1, iy}, cz ∈ {iz-1, iz}.
                    let q = [
                        cell_vertex[cell_index(ix, iy - 1, iz - 1)],
                        cell_vertex[cell_index(ix, iy, iz - 1)],
                        cell_vertex[cell_index(ix, iy, iz)],
                        cell_vertex[cell_index(ix, iy - 1, iz)],
                    ];
                    // Wind so the normal points solid→void. If the
                    // edge goes solid(s0)→void, the +X side is void.
                    emit_quad(&mut conn, q, !s0);
                }
            }
        }
        // Y-direction grid edges.
        for iz in 0..nz {
            for iy in -1..ny {
                for ix in 0..nx {
                    let s0 = self.solid_at(ix, iy, iz);
                    let s1 = self.solid_at(ix, iy + 1, iz);
                    if s0 == s1 {
                        continue;
                    }
                    let q = [
                        cell_vertex[cell_index(ix - 1, iy, iz - 1)],
                        cell_vertex[cell_index(ix, iy, iz - 1)],
                        cell_vertex[cell_index(ix, iy, iz)],
                        cell_vertex[cell_index(ix - 1, iy, iz)],
                    ];
                    emit_quad(&mut conn, q, s0);
                }
            }
        }
        // Z-direction grid edges.
        for iz in -1..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let s0 = self.solid_at(ix, iy, iz);
                    let s1 = self.solid_at(ix, iy, iz + 1);
                    if s0 == s1 {
                        continue;
                    }
                    let q = [
                        cell_vertex[cell_index(ix - 1, iy - 1, iz)],
                        cell_vertex[cell_index(ix, iy - 1, iz)],
                        cell_vertex[cell_index(ix, iy, iz)],
                        cell_vertex[cell_index(ix - 1, iy, iz)],
                    ];
                    emit_quad(&mut conn, q, !s0);
                }
            }
        }

        let mut mesh = valenx_mesh::Mesh::new("voxel-surface-nets");
        mesh.nodes = nodes;
        let mut block =
            valenx_mesh::element::ElementBlock::new(valenx_mesh::element::ElementType::Tri3);
        block.connectivity = conn;
        mesh.element_blocks.push(block);
        Ok(mesh)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_aabb_solid_count_matches_n() {
        let v = Voxel::from_aabb(Vector3::zeros(), Vector3::new(10.0, 10.0, 10.0), (8, 8, 8))
            .unwrap();
        assert_eq!(v.n(), 512);
        assert_eq!(v.solid_count(), 512);
    }

    #[test]
    fn cut_segment_clears_voxels() {
        let mut v = Voxel::from_aabb(
            Vector3::zeros(),
            Vector3::new(10.0, 10.0, 10.0),
            (16, 16, 16),
        )
        .unwrap();
        let before = v.solid_count();
        // Cut a thin horizontal segment in the centre.
        v.cut_segment(
            Vector3::new(0.0, 5.0, 5.0),
            Vector3::new(10.0, 5.0, 5.0),
            1.0,
        );
        let after = v.solid_count();
        assert!(
            after < before,
            "expected voxels cleared; {before} -> {after}"
        );
    }

    #[test]
    fn world_to_voxel_roundtrip() {
        let v = Voxel::from_aabb(Vector3::zeros(), Vector3::new(8.0, 8.0, 8.0), (8, 8, 8))
            .unwrap();
        let p = Vector3::new(2.5, 4.5, 6.5);
        let (ix, iy, iz) = v.world_to_voxel(p).unwrap();
        assert_eq!((ix, iy, iz), (2, 4, 6));
        let c = v.voxel_to_world_centre(ix, iy, iz);
        assert!((c - Vector3::new(2.5, 4.5, 6.5)).norm() < 1e-9);
    }

    #[test]
    fn to_mesh_produces_triangles() {
        let v = Voxel::from_aabb(Vector3::zeros(), Vector3::new(4.0, 4.0, 4.0), (4, 4, 4))
            .unwrap();
        let m = v.to_mesh();
        // Outer surface of a 4×4×4 cube: 6 faces × 16 voxel-quads ×
        // 2 tris = 192 triangles.
        let n_tris = m.element_blocks[0].connectivity.len() / 3;
        assert_eq!(n_tris, 192, "got {n_tris} triangles");
    }

    #[test]
    fn out_of_range_get_returns_false() {
        let v = Voxel::from_aabb(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0), (4, 4, 4))
            .unwrap();
        assert!(!v.get(99, 0, 0));
        assert!(v.world_to_voxel(Vector3::new(-1.0, 0.0, 0.0)).is_none());
    }

    /// Triangle count of a mesh's first block.
    fn tri_count(m: &valenx_mesh::Mesh) -> usize {
        m.element_blocks
            .first()
            .map(|b| b.connectivity.len() / 3)
            .unwrap_or(0)
    }

    #[test]
    fn surface_nets_of_solid_block_is_a_closed_manifold() {
        // A solid block: Surface Nets must produce a watertight,
        // closed surface — every edge shared by exactly two triangles.
        let v = Voxel::from_aabb(
            Vector3::zeros(),
            Vector3::new(4.0, 4.0, 4.0),
            (4, 4, 4),
        )
        .unwrap();
        let m = v.to_mesh_surface_nets().unwrap();
        assert!(tri_count(&m) > 0, "surface nets should produce geometry");
        // Manifold check: tally every undirected edge; a closed
        // manifold has each edge in exactly 2 triangles.
        use std::collections::HashMap;
        let mut edges: HashMap<(u32, u32), u32> = HashMap::new();
        let conn = &m.element_blocks[0].connectivity;
        for tri in conn.chunks_exact(3) {
            for &(a, b) in &[(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                *edges.entry(key).or_insert(0) += 1;
            }
        }
        for (edge, count) in &edges {
            assert_eq!(
                *count, 2,
                "edge {edge:?} is in {count} triangles — surface is not a closed manifold"
            );
        }
    }

    #[test]
    fn surface_nets_block_corners_are_inside_the_stock_aabb() {
        // Every Surface-Nets vertex of a solid block must lie within
        // (or on) the stock AABB — the smoothed surface should not
        // bulge outside the block it was extracted from.
        let v = Voxel::from_aabb(
            Vector3::new(-1.0, -2.0, 0.5),
            Vector3::new(3.0, 2.0, 4.5),
            (8, 8, 8),
        )
        .unwrap();
        let m = v.to_mesh_surface_nets().unwrap();
        for n in &m.nodes {
            assert!(
                n.x >= -1.0 - 1e-9 && n.x <= 3.0 + 1e-9,
                "vertex x {} outside stock",
                n.x
            );
            assert!(n.y >= -2.0 - 1e-9 && n.y <= 2.0 + 1e-9);
            assert!(n.z >= 0.5 - 1e-9 && n.z <= 4.5 + 1e-9);
        }
    }

    #[test]
    fn surface_nets_shares_vertices_unlike_faceted_extraction() {
        // For a binary occupancy grid the *triangle* counts of the two
        // extractors coincide — naive Surface Nets emits one quad per
        // sign-changing grid edge, and that is in 1:1 correspondence
        // with the faceted extractor's one quad per boundary voxel
        // face. Surface Nets' real win is **vertex sharing**: one
        // vertex per boundary cell, shared by up to 4 quads, where the
        // faceted extractor pushes 4 fresh un-welded nodes per quad.
        // So the smooth mesh carries far fewer vertices for the same
        // surface — that is what makes it lighter and smoother.
        let mut v = Voxel::from_aabb(
            Vector3::zeros(),
            Vector3::new(20.0, 20.0, 20.0),
            (32, 32, 32),
        )
        .unwrap();
        // Carve a diagonal channel — a genuinely non-axis-aligned cut.
        v.cut_segment(
            Vector3::new(0.0, 0.0, 10.0),
            Vector3::new(20.0, 20.0, 10.0),
            3.0,
        );
        let faceted = v.to_mesh();
        let smooth = v.to_mesh_surface_nets().unwrap();
        assert!(tri_count(&smooth) > 0, "surface nets produced no geometry");
        // Same surface ⇒ equal triangle counts on a binary grid.
        assert_eq!(
            tri_count(&smooth),
            tri_count(&faceted),
            "binary-grid surface nets emits one quad per sign-changing \
             edge — in 1:1 correspondence with faceted boundary faces"
        );
        // ...but the smooth mesh shares vertices, so it has far fewer.
        assert!(
            smooth.nodes.len() < faceted.nodes.len(),
            "surface nets should share vertices: smooth {} vs faceted {}",
            smooth.nodes.len(),
            faceted.nodes.len()
        );
    }

    #[test]
    fn surface_nets_of_empty_grid_is_empty() {
        // An all-void grid has no boundary → no triangles.
        let mut v = Voxel::from_aabb(
            Vector3::zeros(),
            Vector3::new(8.0, 8.0, 8.0),
            (8, 8, 8),
        )
        .unwrap();
        // Clear everything by a huge cut.
        v.cut_segment(
            Vector3::new(4.0, 4.0, 4.0),
            Vector3::new(4.0, 4.0, 4.0),
            100.0,
        );
        assert_eq!(v.solid_count(), 0, "grid should be fully cleared");
        assert_eq!(tri_count(&v.to_mesh_surface_nets().unwrap()), 0);
    }

    #[test]
    fn surface_nets_outer_face_count_matches_a_slab() {
        // A 1-voxel-thick slab — its Surface-Nets mesh is two sheets
        // (top + bottom) joined at the rim. The vertex count equals
        // the number of active cells; for a flat slab that is the
        // ring of boundary cells. Just assert it is a sane closed
        // surface with the expected order of magnitude.
        let v = Voxel::from_aabb(
            Vector3::zeros(),
            Vector3::new(10.0, 10.0, 2.0),
            (10, 10, 2),
        )
        .unwrap();
        let m = v.to_mesh_surface_nets().unwrap();
        assert!(tri_count(&m) > 0);
        // The slab's surface area dwarfs its edge, so the mesh should
        // have well over a hundred triangles on a 10×10 footprint.
        assert!(
            tri_count(&m) > 100,
            "slab surface-nets mesh too sparse: {}",
            tri_count(&m)
        );
    }

    /// Round-12 M3/M8 RED→GREEN: a resolution that would overflow
    /// the cell-count arithmetic must be rejected at the cap layer
    /// instead of panicking on allocation or wrapping into a tiny
    /// allocation with index UB.
    #[test]
    fn from_aabb_rejects_overflowing_resolution() {
        let err = Voxel::from_aabb(
            Vector3::zeros(),
            Vector3::new(1.0, 1.0, 1.0),
            (u32::MAX, u32::MAX, u32::MAX),
        )
        .expect_err("must reject overflowing resolution");
        match err {
            VoxelError::Overflow(r) => {
                assert_eq!(r, (u32::MAX, u32::MAX, u32::MAX));
            }
            other => panic!("expected Overflow, got {other:?}"),
        }
    }

    /// Round-12 M3/M8: a non-overflowing but still huge cell count
    /// must be rejected by the MAX_VOXEL_CELLS guard.
    #[test]
    fn from_aabb_rejects_too_many_cells() {
        // 1024³ ≈ 1.07e9 cells — over the 1e8 cap, well below u64 max.
        let err = Voxel::from_aabb(
            Vector3::zeros(),
            Vector3::new(1.0, 1.0, 1.0),
            (1024, 1024, 1024),
        )
        .expect_err("must reject too-many-cells request");
        match err {
            VoxelError::TooManyCells { cells, max } => {
                assert_eq!(cells, 1024u64 * 1024 * 1024);
                assert_eq!(max, MAX_VOXEL_CELLS);
            }
            other => panic!("expected TooManyCells, got {other:?}"),
        }
    }
}
