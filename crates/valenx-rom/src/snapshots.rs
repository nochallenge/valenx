//! The snapshot container shared by every method in this crate.
//!
//! A [`Snapshots`] wraps a column-major state-in-time matrix: each **column**
//! is the full state vector at one time sample, and the columns run forward in
//! time. This is the universal input to POD, DMD, and Operator Inference.
//!
//! The type is a thin, validated wrapper around an [`nalgebra::DMatrix<f64>`].
//! Construction is fail-loud — an empty matrix or one holding a non-finite
//! entry is rejected at the door, so downstream solvers never have to second
//! guess their input.

use nalgebra::DMatrix;
use ndarray::Array2;

use crate::error::RomError;

/// A column-major snapshot matrix: `rows` = state dimension, `cols` = time.
///
/// Column `k` is the system state at time sample `k`; columns are ordered in
/// time. Validated on construction to be non-empty and finite.
#[derive(Debug, Clone, PartialEq)]
pub struct Snapshots {
    data: DMatrix<f64>,
}

impl Snapshots {
    /// Wrap an existing [`DMatrix`] (columns = time).
    ///
    /// # Errors
    /// - [`RomError::Empty`] if the matrix has zero rows or zero columns.
    /// - [`RomError::NonFinite`] if any entry is `NaN` or infinite.
    pub fn from_matrix(data: DMatrix<f64>) -> Result<Self, RomError> {
        let (rows, cols) = (data.nrows(), data.ncols());
        if rows == 0 || cols == 0 {
            return Err(RomError::Empty { rows, cols });
        }
        if data.iter().any(|v| !v.is_finite()) {
            return Err(RomError::NonFinite { what: "snapshots" });
        }
        Ok(Self { data })
    }

    /// Build from an [`ndarray::Array2<f64>`] whose columns are time samples.
    ///
    /// Provided for callers that assemble snapshots in `ndarray` (the workspace
    /// snapshot-container crate). The data is copied into the [`nalgebra`]
    /// representation the solvers use.
    ///
    /// # Errors
    /// Same as [`Snapshots::from_matrix`].
    pub fn from_ndarray(arr: &Array2<f64>) -> Result<Self, RomError> {
        let (rows, cols) = arr.dim();
        if rows == 0 || cols == 0 {
            return Err(RomError::Empty { rows, cols });
        }
        // ndarray is row-major; build the nalgebra matrix explicitly so the
        // (row, col) mapping is unambiguous regardless of memory order.
        let data = DMatrix::from_fn(rows, cols, |r, c| arr[(r, c)]);
        Self::from_matrix(data)
    }

    /// Build column-by-column from per-time-sample state vectors.
    ///
    /// `columns[k]` is the state at time `k`; every column must share the same
    /// length (the state dimension).
    ///
    /// # Errors
    /// - [`RomError::Empty`] if there are no columns or the columns are empty.
    /// - [`RomError::DimensionMismatch`] if columns differ in length.
    /// - [`RomError::NonFinite`] if any entry is non-finite.
    pub fn from_columns(columns: &[Vec<f64>]) -> Result<Self, RomError> {
        let cols = columns.len();
        if cols == 0 {
            return Err(RomError::Empty { rows: 0, cols: 0 });
        }
        let rows = columns[0].len();
        if rows == 0 {
            return Err(RomError::Empty { rows: 0, cols });
        }
        for (k, c) in columns.iter().enumerate() {
            if c.len() != rows {
                return Err(RomError::DimensionMismatch {
                    what: "snapshot column",
                    expected: rows,
                    got: c.len(),
                });
            }
            // record which column to keep the error specific; index unused
            let _ = k;
        }
        let data = DMatrix::from_fn(rows, cols, |r, c| columns[c][r]);
        Self::from_matrix(data)
    }

    /// The state dimension (number of rows).
    pub fn state_dim(&self) -> usize {
        self.data.nrows()
    }

    /// The number of time samples (number of columns).
    pub fn n_time(&self) -> usize {
        self.data.ncols()
    }

    /// Borrow the underlying column-major matrix (columns = time).
    pub fn matrix(&self) -> &DMatrix<f64> {
        &self.data
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_empty() {
        let e = Snapshots::from_matrix(DMatrix::<f64>::zeros(0, 3)).unwrap_err();
        assert_eq!(e.code(), "empty");
        let e = Snapshots::from_columns(&[]).unwrap_err();
        assert_eq!(e.code(), "empty");
    }

    #[test]
    fn rejects_non_finite() {
        let m = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, f64::NAN, 4.0]);
        assert_eq!(Snapshots::from_matrix(m).unwrap_err().code(), "non_finite");
    }

    #[test]
    fn rejects_ragged_columns() {
        let cols = vec![vec![1.0, 2.0], vec![3.0]];
        assert_eq!(
            Snapshots::from_columns(&cols).unwrap_err().code(),
            "dimension_mismatch"
        );
    }

    #[test]
    fn from_columns_is_column_major() {
        let cols = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let s = Snapshots::from_columns(&cols).unwrap();
        assert_eq!(s.state_dim(), 3);
        assert_eq!(s.n_time(), 2);
        assert_eq!(s.matrix()[(0, 0)], 1.0);
        assert_eq!(s.matrix()[(2, 1)], 6.0);
    }

    #[test]
    fn from_ndarray_matches_matrix() {
        let arr = Array2::from_shape_vec((2, 3), vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).unwrap();
        let s = Snapshots::from_ndarray(&arr).unwrap();
        assert_eq!(s.state_dim(), 2);
        assert_eq!(s.n_time(), 3);
        // row 0 = [1,2,3], row 1 = [4,5,6]
        assert_eq!(s.matrix()[(0, 2)], 3.0);
        assert_eq!(s.matrix()[(1, 0)], 4.0);
    }
}
