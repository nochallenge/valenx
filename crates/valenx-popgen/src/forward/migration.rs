//! Migration models for structured populations.
//!
//! When a population is split into `d` *demes* (sub-populations),
//! migration moves individuals between them each generation. A
//! [`MigrationModel`] is a row-stochastic `d x d` matrix `M`: `M[i][j]`
//! is the probability that an individual sampled to be in deme `i` this
//! generation actually descends from a parent in deme `j` (a
//! *backward* migration rate, the convention used in coalescent
//! theory).
//!
//! Two standard topologies have constructors:
//!
//! - **Island model** ([`MigrationModel::island`]) — every deme
//!   exchanges migrants with every other deme at the same rate `m`.
//!   `M[i][j] = m / (d - 1)` off the diagonal, `1 - m` on it.
//! - **Stepping-stone model** ([`MigrationModel::stepping_stone`]) —
//!   migrants move only to the two nearest-neighbour demes (a 1-D
//!   lattice). Optionally circular (a ring).

use crate::error::{PopgenError, Result};
use crate::rng::Rng;
use serde::{Deserialize, Serialize};

/// A backward migration matrix over `d` demes.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct MigrationModel {
    /// Row-stochastic `d x d` backward-migration matrix.
    matrix: Vec<Vec<f64>>,
}

impl MigrationModel {
    /// Builds a model from an explicit `d x d` matrix.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] on an empty matrix or a negative entry;
    /// [`PopgenError::Dimension`] if the matrix is not square;
    /// [`PopgenError::Invalid`] if a row does not sum to ~1.
    pub fn new(matrix: Vec<Vec<f64>>) -> Result<Self> {
        let d = matrix.len();
        if d == 0 {
            return Err(PopgenError::invalid("matrix", "no demes"));
        }
        for (i, row) in matrix.iter().enumerate() {
            if row.len() != d {
                return Err(PopgenError::dimension(d, row.len(), "migration row"));
            }
            if row.iter().any(|&v| v < 0.0) {
                return Err(PopgenError::invalid(
                    "matrix",
                    format!("row {i} has a negative rate"),
                ));
            }
            let sum: f64 = row.iter().sum();
            if (sum - 1.0).abs() > 1e-6 {
                return Err(PopgenError::invalid(
                    "matrix",
                    format!("row {i} sums to {sum}, not 1"),
                ));
            }
        }
        Ok(MigrationModel { matrix })
    }

    /// Builds an island model: `d` demes, total migration rate `m`
    /// spread evenly over all other demes.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `d < 2` or `m` is outside `[0, 1]`.
    pub fn island(d: usize, m: f64) -> Result<Self> {
        if d < 2 {
            return Err(PopgenError::invalid("d", "need at least two demes"));
        }
        if !(0.0..=1.0).contains(&m) {
            return Err(PopgenError::invalid(
                "m",
                "migration rate must lie in [0, 1]",
            ));
        }
        let off = m / (d - 1) as f64;
        let matrix = (0..d)
            .map(|i| {
                (0..d)
                    .map(|j| if i == j { 1.0 - m } else { off })
                    .collect()
            })
            .collect();
        Ok(MigrationModel { matrix })
    }

    /// Builds a 1-D stepping-stone model: each deme exchanges migrants
    /// with its immediate neighbours at total rate `m`. If `circular`
    /// the lattice wraps into a ring.
    ///
    /// # Errors
    /// [`PopgenError::Invalid`] if `d < 2` or `m` is outside `[0, 1]`.
    pub fn stepping_stone(d: usize, m: f64, circular: bool) -> Result<Self> {
        if d < 2 {
            return Err(PopgenError::invalid("d", "need at least two demes"));
        }
        if !(0.0..=1.0).contains(&m) {
            return Err(PopgenError::invalid(
                "m",
                "migration rate must lie in [0, 1]",
            ));
        }
        let mut matrix = vec![vec![0.0; d]; d];
        // The index `i` addresses arbitrary neighbour columns of the
        // matrix, so an enumerate-iterator does not apply here.
        #[allow(clippy::needless_range_loop)]
        for i in 0..d {
            // Neighbour indices.
            let mut neighbours = Vec::new();
            if i > 0 {
                neighbours.push(i - 1);
            } else if circular {
                neighbours.push(d - 1);
            }
            if i + 1 < d {
                neighbours.push(i + 1);
            } else if circular {
                neighbours.push(0);
            }
            let share = if neighbours.is_empty() {
                0.0
            } else {
                m / neighbours.len() as f64
            };
            for &n in &neighbours {
                matrix[i][n] += share;
            }
            // The diagonal absorbs whatever did not migrate.
            let migrated: f64 = matrix[i].iter().sum();
            matrix[i][i] = 1.0 - migrated;
        }
        Ok(MigrationModel { matrix })
    }

    /// Number of demes.
    pub fn deme_count(&self) -> usize {
        self.matrix.len()
    }

    /// The backward-migration matrix.
    pub fn matrix(&self) -> &[Vec<f64>] {
        &self.matrix
    }

    /// Samples the source deme a deme-`i` offspring's parent came from.
    ///
    /// Returns `i` itself if `i` is out of range.
    pub fn sample_source(&self, i: usize, rng: &mut Rng) -> usize {
        if i >= self.matrix.len() {
            return i;
        }
        rng.weighted_index(&self.matrix[i])
    }

    /// Fraction of deme `i` that is non-resident (the realised
    /// immigration rate `1 - M[i][i]`).
    pub fn immigration_rate(&self, i: usize) -> f64 {
        self.matrix.get(i).map(|r| 1.0 - r[i]).unwrap_or(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn island_rows_are_stochastic() {
        let m = MigrationModel::island(4, 0.1).unwrap();
        assert_eq!(m.deme_count(), 4);
        for row in m.matrix() {
            let sum: f64 = row.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12);
        }
        // Off-diagonal entry = m / (d-1) = 0.1/3.
        assert!((m.matrix()[0][1] - 0.1 / 3.0).abs() < 1e-12);
        assert!((m.immigration_rate(0) - 0.1).abs() < 1e-12);
    }

    #[test]
    fn stepping_stone_linear_endpoints_have_one_neighbour() {
        let m = MigrationModel::stepping_stone(5, 0.2, false).unwrap();
        // Deme 0: only neighbour is deme 1 -> all migration goes there.
        assert!((m.matrix()[0][1] - 0.2).abs() < 1e-12);
        assert!(m.matrix()[0][2].abs() < 1e-12);
        // Interior deme 2 splits 0.2 between demes 1 and 3.
        assert!((m.matrix()[2][1] - 0.1).abs() < 1e-12);
        assert!((m.matrix()[2][3] - 0.1).abs() < 1e-12);
        for row in m.matrix() {
            assert!((row.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn stepping_stone_circular_has_no_endpoints() {
        let m = MigrationModel::stepping_stone(4, 0.2, true).unwrap();
        // Every deme has two neighbours in a ring.
        for i in 0..4 {
            assert!((m.immigration_rate(i) - 0.2).abs() < 1e-12);
        }
    }

    #[test]
    fn sample_source_respects_the_matrix() {
        // Deme 0 with no migration always stays in deme 0.
        let m = MigrationModel::island(3, 0.0).unwrap();
        let mut rng = Rng::new(1);
        for _ in 0..100 {
            assert_eq!(m.sample_source(0, &mut rng), 0);
        }
    }

    #[test]
    fn rejects_bad_input() {
        assert!(MigrationModel::island(1, 0.1).is_err());
        assert!(MigrationModel::island(3, 1.5).is_err());
        assert!(MigrationModel::stepping_stone(1, 0.1, false).is_err());
        // A non-square matrix.
        assert!(MigrationModel::new(vec![vec![1.0, 0.0], vec![1.0]]).is_err());
        // A row that does not sum to 1.
        assert!(
            MigrationModel::new(vec![vec![0.5, 0.0], vec![0.0, 1.0]]).is_err()
        );
    }
}
