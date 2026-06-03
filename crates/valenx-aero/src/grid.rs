//! The 3-D Cartesian staggered (marker-and-cell) grid and its field
//! storage.
//!
//! # Why a staggered grid
//!
//! A collocated 3-D grid — every variable at the same point — admits
//! the notorious **checkerboard pressure** mode: a zig-zag pressure
//! field the central-difference gradient cannot "see", so the solver
//! converges to a physically wrong oscillating answer. Harlow &
//! Welch's staggered grid cures it by storing each variable where its
//! governing flux is naturally defined:
//!
//! ```text
//!   pressure  p   — at cell centres            (nx   · ny   · nz  )
//!   u-velocity    — at  x-normal cell faces    ((nx+1)· ny   · nz  )
//!   v-velocity    — at  y-normal cell faces    (nx   ·(ny+1)· nz  )
//!   w-velocity    — at  z-normal cell faces    (nx   · ny   ·(nz+1))
//! ```
//!
//! With this layout the pressure gradient driving a face velocity is
//! an *exact* difference of the two adjacent cell pressures — no
//! interpolation, no checkerboard. It is the standard discretisation
//! behind the SIMPLE family of solvers, here extended to three
//! dimensions for external vehicle / aircraft aerodynamics.
//!
//! # Indexing
//!
//! Every field is a flat `Vec<f64>` in **x-fastest, then y, then z**
//! order. A pressure cell `(i, j, k)` is at linear index
//! `i + nx·(j + ny·k)`; the [`Field3`] type wraps the index
//! arithmetic and bounds checks so the solver reads `p.at(i, j, k)`.

/// Geometry of a uniform Cartesian staggered grid over a rectangular
/// box `[0, lx] × [0, ly] × [0, lz]`.
///
/// `nx · ny · nz` is the count of **pressure cells**; each cell is
/// `dx · dy · dz` with `dx = lx/nx`, etc. The grid is uniform — one
/// cell size per axis — which keeps the immersed-boundary voxelization
/// and the geometric-multigrid coarsening simple and robust.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid3 {
    /// Number of pressure cells along x.
    pub nx: usize,
    /// Number of pressure cells along y.
    pub ny: usize,
    /// Number of pressure cells along z.
    pub nz: usize,
    /// Physical domain length along x (m).
    pub lx: f64,
    /// Physical domain length along y (m).
    pub ly: f64,
    /// Physical domain length along z (m).
    pub lz: f64,
    /// World-space `x` of the domain's minimum corner (m).
    pub x0: f64,
    /// World-space `y` of the domain's minimum corner (m).
    pub y0: f64,
    /// World-space `z` of the domain's minimum corner (m).
    pub z0: f64,
}

impl Grid3 {
    /// Build a grid of `nx · ny · nz` cells over the box
    /// `[x0, x0+lx] × [y0, y0+ly] × [z0, z0+lz]`.
    ///
    /// # Panics
    ///
    /// Panics on a zero cell count or a non-positive / non-finite
    /// domain length — those are programmer errors, not runtime input.
    /// Callers that build a grid from external data should validate
    /// first ([`crate::WindTunnel`] does).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        lx: f64,
        ly: f64,
        lz: f64,
        x0: f64,
        y0: f64,
        z0: f64,
    ) -> Grid3 {
        assert!(
            nx > 0 && ny > 0 && nz > 0,
            "grid must have at least one cell per axis"
        );
        assert!(
            lx > 0.0 && ly > 0.0 && lz > 0.0,
            "domain dimensions must be positive"
        );
        assert!(
            lx.is_finite() && ly.is_finite() && lz.is_finite(),
            "domain dimensions must be finite"
        );
        Grid3 {
            nx,
            ny,
            nz,
            lx,
            ly,
            lz,
            x0,
            y0,
            z0,
        }
    }

    /// A grid anchored at the origin — convenience for unit / test use.
    pub fn at_origin(nx: usize, ny: usize, nz: usize, lx: f64, ly: f64, lz: f64) -> Grid3 {
        Grid3::new(nx, ny, nz, lx, ly, lz, 0.0, 0.0, 0.0)
    }

    /// Cell width along x — `lx / nx`.
    #[inline]
    pub fn dx(&self) -> f64 {
        self.lx / self.nx as f64
    }

    /// Cell width along y — `ly / ny`.
    #[inline]
    pub fn dy(&self) -> f64 {
        self.ly / self.ny as f64
    }

    /// Cell width along z — `lz / nz`.
    #[inline]
    pub fn dz(&self) -> f64 {
        self.lz / self.nz as f64
    }

    /// Total pressure-cell count.
    #[inline]
    pub fn cell_count(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// World-space `(x, y, z)` of the centre of pressure cell
    /// `(i, j, k)`.
    #[inline]
    pub fn cell_centre(&self, i: usize, j: usize, k: usize) -> (f64, f64, f64) {
        (
            self.x0 + (i as f64 + 0.5) * self.dx(),
            self.y0 + (j as f64 + 0.5) * self.dy(),
            self.z0 + (k as f64 + 0.5) * self.dz(),
        )
    }

    /// A zero-initialised cell-centred scalar field — `nx · ny · nz`.
    pub fn scalar_field(&self) -> Field3 {
        Field3::zeros(self.nx, self.ny, self.nz)
    }

    /// A zero-initialised staggered `u`-velocity field —
    /// `(nx+1) · ny · nz`, one value per x-normal cell face.
    pub fn u_field(&self) -> Field3 {
        Field3::zeros(self.nx + 1, self.ny, self.nz)
    }

    /// A zero-initialised staggered `v`-velocity field —
    /// `nx · (ny+1) · nz`, one value per y-normal cell face.
    pub fn v_field(&self) -> Field3 {
        Field3::zeros(self.nx, self.ny + 1, self.nz)
    }

    /// A zero-initialised staggered `w`-velocity field —
    /// `nx · ny · (nz+1)`, one value per z-normal cell face.
    pub fn w_field(&self) -> Field3 {
        Field3::zeros(self.nx, self.ny, self.nz + 1)
    }
}

/// A flat 3-D scalar array with `(i, j, k)` indexing.
///
/// Used for every field on the staggered grid — the pressure, the
/// three velocity components, the turbulence scalars, and the solver
/// work arrays. It carries its own `(nx, ny, nz)` extent so the index
/// arithmetic and the bounds checks live in one place.
#[derive(Clone, Debug, PartialEq)]
pub struct Field3 {
    /// Extent along x.
    pub nx: usize,
    /// Extent along y.
    pub ny: usize,
    /// Extent along z.
    pub nz: usize,
    /// `nx · ny · nz` values, x-fastest then y then z.
    pub data: Vec<f64>,
}

impl Field3 {
    /// A zero-filled field of the given extent.
    pub fn zeros(nx: usize, ny: usize, nz: usize) -> Field3 {
        Field3 {
            nx,
            ny,
            nz,
            data: vec![0.0; nx * ny * nz],
        }
    }

    /// A field filled with a constant value.
    pub fn filled(nx: usize, ny: usize, nz: usize, value: f64) -> Field3 {
        Field3 {
            nx,
            ny,
            nz,
            data: vec![value; nx * ny * nz],
        }
    }

    /// Linear index of `(i, j, k)`. Debug-asserts the bounds.
    #[inline]
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        debug_assert!(
            i < self.nx && j < self.ny && k < self.nz,
            "field index ({i},{j},{k}) out of bounds ({},{},{})",
            self.nx,
            self.ny,
            self.nz
        );
        i + self.nx * (j + self.ny * k)
    }

    /// Read the value at `(i, j, k)`.
    #[inline]
    pub fn at(&self, i: usize, j: usize, k: usize) -> f64 {
        self.data[self.index(i, j, k)]
    }

    /// Set the value at `(i, j, k)`.
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, k: usize, value: f64) {
        let idx = self.index(i, j, k);
        self.data[idx] = value;
    }

    /// Add to the value at `(i, j, k)`.
    #[inline]
    pub fn add(&mut self, i: usize, j: usize, k: usize, delta: f64) {
        let idx = self.index(i, j, k);
        self.data[idx] += delta;
    }

    /// Fill the whole field with one value.
    pub fn fill(&mut self, value: f64) {
        self.data.iter_mut().for_each(|v| *v = value);
    }

    /// The largest absolute value across the field.
    pub fn abs_max(&self) -> f64 {
        self.data.iter().fold(0.0, |m, &v| m.max(v.abs()))
    }

    /// The L2 (root-mean-square) norm of the field.
    pub fn l2_norm(&self) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.data.iter().map(|&v| v * v).sum();
        (sum / self.data.len() as f64).sqrt()
    }

    /// The arithmetic mean of the field.
    pub fn mean(&self) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        self.data.iter().sum::<f64>() / self.data.len() as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_cell_sizes_divide_the_domain() {
        let g = Grid3::at_origin(4, 5, 8, 2.0, 10.0, 16.0);
        assert!((g.dx() - 0.5).abs() < 1e-12);
        assert!((g.dy() - 2.0).abs() < 1e-12);
        assert!((g.dz() - 2.0).abs() < 1e-12);
        assert_eq!(g.cell_count(), 4 * 5 * 8);
    }

    #[test]
    fn staggered_field_dimensions_are_correct() {
        // For an nx·ny·nz grid: p is nx·ny·nz, u is (nx+1)·ny·nz,
        // v is nx·(ny+1)·nz, w is nx·ny·(nz+1) — the defining shape
        // of a 3-D MAC grid.
        let g = Grid3::at_origin(8, 6, 4, 1.0, 1.0, 1.0);
        let p = g.scalar_field();
        assert_eq!((p.nx, p.ny, p.nz), (8, 6, 4));
        let u = g.u_field();
        assert_eq!((u.nx, u.ny, u.nz), (9, 6, 4));
        let v = g.v_field();
        assert_eq!((v.nx, v.ny, v.nz), (8, 7, 4));
        let w = g.w_field();
        assert_eq!((w.nx, w.ny, w.nz), (8, 6, 5));
    }

    #[test]
    fn cell_centres_are_offset_by_half_a_cell_and_respect_origin() {
        let g = Grid3::new(2, 2, 2, 2.0, 2.0, 2.0, 10.0, 20.0, 30.0);
        let (x, y, z) = g.cell_centre(0, 0, 0);
        assert!((x - 10.5).abs() < 1e-12);
        assert!((y - 20.5).abs() < 1e-12);
        assert!((z - 30.5).abs() < 1e-12);
    }

    #[test]
    fn field_index_is_x_fastest() {
        let f = Field3::zeros(3, 2, 2);
        assert_eq!(f.index(0, 0, 0), 0);
        assert_eq!(f.index(2, 0, 0), 2); // x advances fastest
        assert_eq!(f.index(0, 1, 0), 3); // next y row
        assert_eq!(f.index(0, 0, 1), 6); // next z slab
    }

    #[test]
    fn field_get_set_add_round_trip() {
        let mut f = Field3::zeros(4, 4, 4);
        f.set(2, 3, 1, 7.5);
        assert_eq!(f.at(2, 3, 1), 7.5);
        f.add(2, 3, 1, 0.5);
        assert_eq!(f.at(2, 3, 1), 8.0);
    }

    #[test]
    fn field_norms_report_magnitude() {
        let mut f = Field3::zeros(2, 1, 2);
        f.set(0, 0, 0, 3.0);
        f.set(1, 0, 1, -4.0);
        assert_eq!(f.abs_max(), 4.0);
        // L2 of [3,0,0,-4] / 4 = sqrt(25/4) = 2.5.
        assert!((f.l2_norm() - 2.5).abs() < 1e-12);
        // mean = -0.25.
        assert!((f.mean() - (-0.25)).abs() < 1e-12);
    }
}
