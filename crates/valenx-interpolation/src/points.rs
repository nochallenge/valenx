//! Validated container of interpolation knots.
//!
//! [`DataPoints`] holds `(x, y)` samples whose abscissae are finite
//! and strictly ascending. Every interpolator in this crate is built
//! from a `DataPoints`, so the "sorted-x required", "no duplicate x",
//! and "all coordinates finite" invariants are checked exactly once,
//! at construction, and never re-validated downstream.

use serde::{Deserialize, Serialize};

use crate::error::InterpolationError;

/// A set of interpolation knots with strictly ascending, finite
/// abscissae.
///
/// Construct with [`DataPoints::new`]; the invariants are guaranteed
/// for any value of this type, so the interpolators can index into
/// the data without re-checking.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DataPoints {
    xs: Vec<f64>,
    ys: Vec<f64>,
}

impl DataPoints {
    /// Build from parallel `x` / `y` slices.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolationError::EmptyData`] if no points are
    /// given, [`InterpolationError::LengthMismatch`] if the slices
    /// differ in length, [`InterpolationError::NonFinite`] for any
    /// NaN / infinite coordinate, [`InterpolationError::DuplicateAbscissa`]
    /// for a repeated `x`, and [`InterpolationError::Unsorted`] for any
    /// `x` that is not strictly greater than its predecessor.
    pub fn new(xs: &[f64], ys: &[f64]) -> Result<Self, InterpolationError> {
        if xs.is_empty() {
            return Err(InterpolationError::EmptyData { needed: 1 });
        }
        if xs.len() != ys.len() {
            return Err(InterpolationError::LengthMismatch {
                x_len: xs.len(),
                y_len: ys.len(),
            });
        }
        for (i, &x) in xs.iter().enumerate() {
            if !x.is_finite() {
                return Err(InterpolationError::NonFinite {
                    axis: "x",
                    index: i,
                });
            }
        }
        for (i, &y) in ys.iter().enumerate() {
            if !y.is_finite() {
                return Err(InterpolationError::NonFinite {
                    axis: "y",
                    index: i,
                });
            }
        }
        for i in 1..xs.len() {
            if xs[i] == xs[i - 1] {
                return Err(InterpolationError::DuplicateAbscissa {
                    index: i,
                    prev_index: i - 1,
                    x: xs[i],
                });
            }
            if xs[i] < xs[i - 1] {
                return Err(InterpolationError::Unsorted {
                    index: i,
                    prev_index: i - 1,
                    x: xs[i],
                    prev: xs[i - 1],
                });
            }
        }
        Ok(Self {
            xs: xs.to_vec(),
            ys: ys.to_vec(),
        })
    }

    /// Number of knots.
    pub fn len(&self) -> usize {
        self.xs.len()
    }

    /// Whether the set is empty.
    ///
    /// Always `false` for a value produced by [`DataPoints::new`],
    /// which rejects empty input; provided for API completeness so
    /// Clippy does not flag a `len` without an `is_empty`.
    pub fn is_empty(&self) -> bool {
        self.xs.is_empty()
    }

    /// The abscissae, in ascending order.
    pub fn xs(&self) -> &[f64] {
        &self.xs
    }

    /// The ordinates, aligned with [`DataPoints::xs`].
    pub fn ys(&self) -> &[f64] {
        &self.ys
    }

    /// First (smallest) abscissa, i.e. the lower domain bound.
    pub fn x_first(&self) -> f64 {
        self.xs[0]
    }

    /// Last (largest) abscissa, i.e. the upper domain bound.
    pub fn x_last(&self) -> f64 {
        self.xs[self.xs.len() - 1]
    }

    /// Reject a query abscissa that lies outside `[x_first, x_last]`.
    ///
    /// A small relative slack absorbs round-off at the endpoints so a
    /// query exactly at `x_last` never spuriously fails.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolationError::OutOfDomain`] when `at` is below
    /// the first knot or above the last (beyond the slack tolerance).
    pub(crate) fn check_domain(&self, at: f64) -> Result<(), InterpolationError> {
        if !at.is_finite() {
            return Err(InterpolationError::OutOfDomain {
                at,
                lo: self.x_first(),
                hi: self.x_last(),
            });
        }
        let lo = self.x_first();
        let hi = self.x_last();
        let span = (hi - lo).abs().max(1.0);
        let slack = span * 1e-12;
        if at < lo - slack || at > hi + slack {
            return Err(InterpolationError::OutOfDomain { at, lo, hi });
        }
        Ok(())
    }

    /// Index `i` of the segment `[xs[i], xs[i + 1]]` that contains
    /// `at` (clamped into `0..len-1`). Assumes `at` is already known
    /// to lie within the domain.
    pub(crate) fn segment_index(&self, at: f64) -> usize {
        // `partition_point` returns the count of leading knots that
        // are <= at; subtract one to land on the left knot of the
        // bracketing segment, clamped to a valid segment.
        let n = self.xs.len();
        let p = self.xs.partition_point(|&x| x <= at);
        if p == 0 {
            0
        } else if p >= n {
            n - 2
        } else {
            p - 1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_empty() {
        let err = DataPoints::new(&[], &[]).unwrap_err();
        assert_eq!(err.code(), "interpolation.empty_data");
    }

    #[test]
    fn rejects_length_mismatch() {
        let err = DataPoints::new(&[0.0, 1.0], &[0.0]).unwrap_err();
        assert_eq!(err.code(), "interpolation.length_mismatch");
    }

    #[test]
    fn rejects_unsorted() {
        let err = DataPoints::new(&[0.0, 2.0, 1.0], &[0.0, 0.0, 0.0]).unwrap_err();
        assert_eq!(err.code(), "interpolation.unsorted");
    }

    #[test]
    fn rejects_duplicate_abscissa() {
        let err = DataPoints::new(&[0.0, 1.0, 1.0], &[0.0, 0.0, 0.0]).unwrap_err();
        assert_eq!(err.code(), "interpolation.duplicate_abscissa");
    }

    #[test]
    fn rejects_non_finite() {
        let err = DataPoints::new(&[0.0, f64::NAN], &[0.0, 1.0]).unwrap_err();
        assert_eq!(err.code(), "interpolation.non_finite");
    }

    #[test]
    fn accepts_sorted_distinct() {
        let dp = DataPoints::new(&[0.0, 1.0, 3.0], &[2.0, 4.0, 8.0]).unwrap();
        assert_eq!(dp.len(), 3);
        assert!(!dp.is_empty());
        let df = dp.x_first() - 0.0;
        let dl = dp.x_last() - 3.0;
        assert!(df.abs() < 1e-15);
        assert!(dl.abs() < 1e-15);
    }

    #[test]
    fn segment_index_brackets_correctly() {
        let dp = DataPoints::new(&[0.0, 1.0, 2.0, 3.0], &[0.0; 4]).unwrap();
        assert_eq!(dp.segment_index(0.0), 0);
        assert_eq!(dp.segment_index(0.5), 0);
        assert_eq!(dp.segment_index(1.0), 1);
        assert_eq!(dp.segment_index(2.5), 2);
        assert_eq!(dp.segment_index(3.0), 2);
    }

    #[test]
    fn check_domain_rejects_outside() {
        let dp = DataPoints::new(&[0.0, 1.0], &[0.0, 1.0]).unwrap();
        assert!(dp.check_domain(0.5).is_ok());
        assert!(dp.check_domain(0.0).is_ok());
        assert!(dp.check_domain(1.0).is_ok());
        let err = dp.check_domain(1.5).unwrap_err();
        assert_eq!(err.code(), "interpolation.out_of_domain");
        assert!(dp.check_domain(-0.5).is_err());
    }
}
