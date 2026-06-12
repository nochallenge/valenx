//! The staggered (marker-and-cell) grid and its field storage.
//!
//! # Why a staggered grid
//!
//! A collocated grid — every variable at the same point — admits a
//! notorious **checkerboard pressure** mode: a zig-zag pressure field
//! that the central-difference pressure gradient cannot "see", so the
//! solver happily converges to a physically wrong oscillating answer.
//! Harlow & Welch's 1965 staggered grid cures it by storing each
//! variable where its governing flux is naturally defined:
//!
//! ```text
//!   pressure  p   — at cell centres          (nx × ny  values)
//!   u-velocity    — at vertical   cell faces ((nx+1) × ny values)
//!   v-velocity    — at horizontal cell faces (nx × (ny+1) values)
//! ```
//!
//! With this layout the pressure gradient driving a face velocity is an
//! *exact* difference of the two adjacent cell pressures — no
//! interpolation, no checkerboard. It is the standard discretisation
//! behind the SIMPLE family of solvers.
//!
//! # Indexing
//!
//! Every field is a flat `Vec<f64>` in **row-major, x-fastest** order.
//! A pressure cell `(i, j)` (`i` along x, `j` along y) is at linear
//! index `i + j·nx`; the [`Field`] type wraps the index arithmetic so
//! the solver reads `p[(i, j)]`.

/// Geometry of a uniform rectangular staggered grid.
///
/// `nx` × `ny` is the count of **pressure cells**; the domain is the
/// rectangle `[0, lx] × [0, ly]`, so each cell is `dx × dy` with
/// `dx = lx/nx`, `dy = ly/ny`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid {
    /// Number of pressure cells along x.
    pub nx: usize,
    /// Number of pressure cells along y.
    pub ny: usize,
    /// Physical domain length along x.
    pub lx: f64,
    /// Physical domain length along y.
    pub ly: f64,
}

impl Grid {
    /// Build a grid of `nx × ny` cells over the `lx × ly` domain.
    ///
    /// # Panics
    ///
    /// Panics on a zero cell count or a non-positive domain length —
    /// those are programmer errors, not runtime inputs.
    pub fn new(nx: usize, ny: usize, lx: f64, ly: f64) -> Grid {
        assert!(
            nx > 0 && ny > 0,
            "grid must have at least one cell per axis"
        );
        assert!(
            lx > 0.0 && ly > 0.0 && lx.is_finite() && ly.is_finite(),
            "domain dimensions must be positive and finite"
        );
        Grid { nx, ny, lx, ly }
    }

    /// Cell width along x — `lx / nx`.
    #[inline]
    pub fn dx(&self) -> f64 {
        self.lx / self.nx as f64
    }

    /// Cell height along y — `ly / ny`.
    #[inline]
    pub fn dy(&self) -> f64 {
        self.ly / self.ny as f64
    }

    /// World-space `(x, y)` of the centre of pressure cell `(i, j)`.
    #[inline]
    pub fn cell_centre(&self, i: usize, j: usize) -> (f64, f64) {
        ((i as f64 + 0.5) * self.dx(), (j as f64 + 0.5) * self.dy())
    }

    /// A zero-initialised [`Field`] for the cell-centred pressure
    /// (`nx × ny`).
    pub fn pressure_field(&self) -> Field {
        Field::zeros(self.nx, self.ny)
    }

    /// A zero-initialised [`Field`] for the staggered `u`-velocity
    /// component — `(nx+1) × ny`, one value per vertical cell face.
    pub fn u_field(&self) -> Field {
        Field::zeros(self.nx + 1, self.ny)
    }

    /// A zero-initialised [`Field`] for the staggered `v`-velocity
    /// component — `nx × (ny+1)`, one value per horizontal cell face.
    pub fn v_field(&self) -> Field {
        Field::zeros(self.nx, self.ny + 1)
    }
}

/// A flat 2-D scalar array with `(i, j)` indexing.
///
/// Used for every field on the staggered grid — the pressure, the two
/// velocity components, and the various solver work arrays. It carries
/// its own `(width, height)` so the index arithmetic and bounds checks
/// live in one place.
#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    /// Number of columns (x extent).
    pub width: usize,
    /// Number of rows (y extent).
    pub height: usize,
    /// `width · height` values, row-major (x-fastest).
    pub data: Vec<f64>,
}

impl Field {
    /// A zero-filled field of the given dimensions.
    pub fn zeros(width: usize, height: usize) -> Field {
        Field {
            width,
            height,
            data: vec![0.0; width * height],
        }
    }

    /// A field filled with a constant value.
    pub fn filled(width: usize, height: usize, value: f64) -> Field {
        Field {
            width,
            height,
            data: vec![value; width * height],
        }
    }

    /// Linear index of `(i, j)`. Debug-asserts the bounds; in release
    /// an out-of-range index simply wraps via the `Vec` panic.
    #[inline]
    pub fn index(&self, i: usize, j: usize) -> usize {
        debug_assert!(
            i < self.width && j < self.height,
            "field index ({i}, {j}) out of bounds ({}, {})",
            self.width,
            self.height
        );
        i + j * self.width
    }

    /// Read the value at `(i, j)`.
    #[inline]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        self.data[self.index(i, j)]
    }

    /// Set the value at `(i, j)`.
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, value: f64) {
        let idx = self.index(i, j);
        self.data[idx] = value;
    }

    /// Add to the value at `(i, j)`.
    #[inline]
    pub fn add(&mut self, i: usize, j: usize, delta: f64) {
        let idx = self.index(i, j);
        self.data[idx] += delta;
    }

    /// Fill the whole field with one value.
    pub fn fill(&mut self, value: f64) {
        self.data.iter_mut().for_each(|v| *v = value);
    }

    /// The largest absolute value across the field — a convenient
    /// residual / convergence norm.
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_cell_sizes_divide_the_domain() {
        let g = Grid::new(4, 5, 2.0, 10.0);
        assert!((g.dx() - 0.5).abs() < 1e-12);
        assert!((g.dy() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn staggered_field_dimensions_are_correct() {
        // For an nx × ny grid: p is nx×ny, u is (nx+1)×ny, v is
        // nx×(ny+1) — the defining shape of a MAC grid.
        let g = Grid::new(8, 6, 1.0, 1.0);
        assert_eq!(
            (g.pressure_field().width, g.pressure_field().height),
            (8, 6)
        );
        assert_eq!((g.u_field().width, g.u_field().height), (9, 6));
        assert_eq!((g.v_field().width, g.v_field().height), (8, 7));
    }

    #[test]
    fn cell_centres_are_offset_by_half_a_cell() {
        let g = Grid::new(2, 2, 2.0, 2.0);
        // Cell (0,0) centre is at (0.5, 0.5) for dx = dy = 1.
        let (x, y) = g.cell_centre(0, 0);
        assert!((x - 0.5).abs() < 1e-12 && (y - 0.5).abs() < 1e-12);
        let (x1, y1) = g.cell_centre(1, 1);
        assert!((x1 - 1.5).abs() < 1e-12 && (y1 - 1.5).abs() < 1e-12);
    }

    #[test]
    fn field_index_is_row_major_x_fastest() {
        let f = Field::zeros(3, 2);
        assert_eq!(f.index(0, 0), 0);
        assert_eq!(f.index(2, 0), 2); // x advances fastest
        assert_eq!(f.index(0, 1), 3); // next row
        assert_eq!(f.index(2, 1), 5);
    }

    #[test]
    fn field_get_set_round_trip() {
        let mut f = Field::zeros(4, 4);
        f.set(2, 3, 7.5);
        assert_eq!(f.at(2, 3), 7.5);
        f.add(2, 3, 0.5);
        assert_eq!(f.at(2, 3), 8.0);
    }

    #[test]
    fn field_norms_report_magnitude() {
        let mut f = Field::zeros(2, 2);
        f.set(0, 0, 3.0);
        f.set(1, 1, -4.0);
        assert_eq!(f.abs_max(), 4.0);
        // L2 of [3, 0, 0, -4] / 4 = sqrt(25/4) = 2.5.
        assert!((f.l2_norm() - 2.5).abs() < 1e-12);
    }

    #[test]
    fn filled_field_has_the_constant_value() {
        let f = Field::filled(3, 3, 2.0);
        assert!(f.data.iter().all(|&v| v == 2.0));
        assert_eq!(f.abs_max(), 2.0);
    }
}
