//! A uniform-grid spatial hash for neighbour search.
//!
//! SPH needs, for each particle, the list of other particles within the kernel
//! support radius `h`. Done naively that is an `O(N²)` all-pairs scan. The
//! standard acceleration structure is a **uniform grid**: space is partitioned
//! into cubic cells of side `cell_size` (chosen `≥ h`), each particle is binned
//! into the cell containing it, and a neighbour query for a particle only has to
//! inspect the particles in its own cell and the 26 surrounding cells (a 3×3×3
//! block). With a roughly uniform particle density this brings neighbour search
//! down to `O(N)`.
//!
//! This is a *hashed* grid: rather than allocating a dense 3-D array (which would
//! be wasteful for a sparse or unbounded domain) the integer cell coordinate
//! `(i, j, k)` is the key of a hash map whose value is the list of particle
//! indices in that cell. The domain therefore need not be known in advance and
//! particles may occupy negative coordinates.
//!
//! The query returns *candidate* indices — every particle in the 27-cell block —
//! because a cell may be as large as the support; the caller still applies the
//! exact `r ≤ h` test (the kernels in [`crate::kernels`] additionally return `0`
//! beyond `h`, so an over-inclusive candidate is harmless either way).

use std::collections::HashMap;

use nalgebra::Vector3;

use crate::error::FluidError;
use crate::particle::ParticleSystem;

/// Integer cell coordinate (a hash-map key).
type Cell = (i64, i64, i64);

/// A uniform-grid spatial hash over a set of particle positions.
///
/// Rebuilt each step from the current particle positions (particles move, so the
/// binning is only valid for the configuration it was built from).
#[derive(Debug, Clone)]
pub struct SpatialHash {
    cell_size: f64,
    inv_cell_size: f64,
    cells: HashMap<Cell, Vec<usize>>,
}

impl SpatialHash {
    /// Build a hash with the given cell size, leaving it empty.
    ///
    /// `cell_size` should be `≥` the kernel support `h` so that all neighbours
    /// within `h` fall in the 3×3×3 block around a particle's cell.
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] if `cell_size` is not finite and strictly
    /// positive.
    pub fn new(cell_size: f64) -> Result<Self, FluidError> {
        if !(cell_size.is_finite() && cell_size > 0.0) {
            return Err(FluidError::InvalidConfig(format!(
                "grid cell size must be finite and > 0, got {cell_size}"
            )));
        }
        Ok(Self {
            cell_size,
            inv_cell_size: 1.0 / cell_size,
            cells: HashMap::new(),
        })
    }

    /// Build a hash and immediately bin a particle system into it.
    ///
    /// # Errors
    /// [`FluidError::InvalidConfig`] if `cell_size` is invalid (see
    /// [`Self::new`]).
    pub fn build(cell_size: f64, system: &ParticleSystem) -> Result<Self, FluidError> {
        let mut grid = Self::new(cell_size)?;
        grid.rebuild(system);
        Ok(grid)
    }

    /// The cell side length.
    #[must_use]
    pub fn cell_size(&self) -> f64 {
        self.cell_size
    }

    /// The integer cell coordinate containing a position.
    #[must_use]
    pub fn cell_of(&self, p: Vector3<f64>) -> Cell {
        (
            (p.x * self.inv_cell_size).floor() as i64,
            (p.y * self.inv_cell_size).floor() as i64,
            (p.z * self.inv_cell_size).floor() as i64,
        )
    }

    /// Clear and re-bin every particle in `system` by its current position.
    pub fn rebuild(&mut self, system: &ParticleSystem) {
        self.cells.clear();
        for (idx, p) in system.particles().iter().enumerate() {
            let cell = self.cell_of(p.position);
            self.cells.entry(cell).or_default().push(idx);
        }
    }

    /// Clear and re-bin from a slice of positions directly.
    ///
    /// Used by solvers (e.g. PCISPH) that need to bin *predicted* positions that
    /// are not stored in a [`ParticleSystem`]. The index of each position in the
    /// slice is what neighbour queries return.
    pub fn rebuild_from(&mut self, positions: &[Vector3<f64>]) {
        self.cells.clear();
        for (idx, &pos) in positions.iter().enumerate() {
            let cell = self.cell_of(pos);
            self.cells.entry(cell).or_default().push(idx);
        }
    }

    /// The number of occupied cells.
    #[must_use]
    pub fn occupied_cells(&self) -> usize {
        self.cells.len()
    }

    /// Append the candidate-neighbour indices for `position` (its own cell plus
    /// the 26 surrounding cells) into `out`.
    ///
    /// `out` is cleared first. The returned indices are *candidates*: the caller
    /// must still apply the exact `‖r_i − r_j‖ ≤ h` distance test. Reusing one
    /// `out` buffer across queries avoids per-query allocation in the hot loop.
    pub fn neighbors_into(&self, position: Vector3<f64>, out: &mut Vec<usize>) {
        out.clear();
        let (ci, cj, ck) = self.cell_of(position);
        for di in -1..=1 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    if let Some(bucket) = self.cells.get(&(ci + di, cj + dj, ck + dk)) {
                        out.extend_from_slice(bucket);
                    }
                }
            }
        }
    }

    /// Convenience wrapper around [`Self::neighbors_into`] that allocates a
    /// fresh `Vec` (prefer the buffer-reusing variant in hot loops).
    #[must_use]
    pub fn neighbors(&self, position: Vector3<f64>) -> Vec<usize> {
        let mut out = Vec::new();
        self.neighbors_into(position, &mut out);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::particle::Particle;

    #[test]
    fn rejects_bad_cell_size() {
        assert!(SpatialHash::new(0.0).is_err());
        assert!(SpatialHash::new(-1.0).is_err());
        assert!(SpatialHash::new(f64::NAN).is_err());
        assert!(SpatialHash::new(0.5).is_ok());
    }

    #[test]
    fn cell_of_handles_negative_coordinates() {
        let g = SpatialHash::new(1.0).unwrap();
        assert_eq!(g.cell_of(Vector3::new(0.5, 0.5, 0.5)), (0, 0, 0));
        assert_eq!(g.cell_of(Vector3::new(-0.5, -0.5, -0.5)), (-1, -1, -1));
        assert_eq!(g.cell_of(Vector3::new(2.3, -2.3, 0.0)), (2, -3, 0));
    }

    #[test]
    fn neighbors_finds_only_nearby_particles() {
        let mut sys = ParticleSystem::new();
        // Three particles close together near the origin...
        sys.push(Particle::at(Vector3::new(0.0, 0.0, 0.0))).unwrap();
        sys.push(Particle::at(Vector3::new(0.1, 0.0, 0.0))).unwrap();
        sys.push(Particle::at(Vector3::new(-0.1, 0.1, 0.0)))
            .unwrap();
        // ...and one far away (many cells over).
        sys.push(Particle::at(Vector3::new(50.0, 50.0, 50.0)))
            .unwrap();

        let grid = SpatialHash::build(0.5, &sys).unwrap();

        // A query at the origin sees the three local particles but not the
        // distant one (it is far outside the 3×3×3 block).
        let mut near = grid.neighbors(Vector3::new(0.0, 0.0, 0.0));
        near.sort_unstable();
        assert_eq!(near, vec![0, 1, 2]);

        // A query at the far particle sees only itself.
        let far = grid.neighbors(Vector3::new(50.0, 50.0, 50.0));
        assert_eq!(far, vec![3]);
    }

    #[test]
    fn neighbors_spans_adjacent_cells() {
        // Two particles either side of a cell boundary are still neighbours
        // because the query inspects the 3×3×3 block, not just one cell.
        let mut sys = ParticleSystem::new();
        sys.push(Particle::at(Vector3::new(0.49, 0.0, 0.0)))
            .unwrap(); // cell 0
        sys.push(Particle::at(Vector3::new(0.51, 0.0, 0.0)))
            .unwrap(); // cell 1
        let grid = SpatialHash::build(0.5, &sys).unwrap();
        let mut n = grid.neighbors(Vector3::new(0.49, 0.0, 0.0));
        n.sort_unstable();
        assert_eq!(n, vec![0, 1]);
    }

    #[test]
    fn rebuild_reflects_moved_particles() {
        let mut sys = ParticleSystem::new();
        sys.push(Particle::at(Vector3::new(0.0, 0.0, 0.0))).unwrap();
        let mut grid = SpatialHash::build(0.5, &sys).unwrap();
        assert_eq!(grid.neighbors(Vector3::new(0.0, 0.0, 0.0)), vec![0]);

        // Move the particle far away and rebuild.
        sys.particles_mut()[0].position = Vector3::new(100.0, 0.0, 0.0);
        grid.rebuild(&sys);
        assert!(grid.neighbors(Vector3::new(0.0, 0.0, 0.0)).is_empty());
        assert_eq!(grid.neighbors(Vector3::new(100.0, 0.0, 0.0)), vec![0]);
    }
}
