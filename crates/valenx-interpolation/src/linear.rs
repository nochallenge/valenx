//! Piecewise-linear interpolation.
//!
//! On each segment `[x_i, x_{i+1}]` the interpolant is the straight
//! line joining `(x_i, y_i)` and `(x_{i+1}, y_{i+1})`. It passes
//! through every knot exactly and is `C0` (continuous) but not
//! generally `C1`.

use crate::error::InterpolationError;
use crate::points::DataPoints;

/// A piecewise-linear interpolant through a set of knots.
///
/// Build with [`LinearInterpolator::new`]; query with
/// [`LinearInterpolator::value_at`].
#[derive(Clone, Debug, PartialEq)]
pub struct LinearInterpolator {
    data: DataPoints,
}

impl LinearInterpolator {
    /// Build a linear interpolant from validated [`DataPoints`].
    ///
    /// A single point is allowed: the interpolant is then the
    /// constant `y0` on the degenerate domain `[x0, x0]`.
    pub fn new(data: DataPoints) -> Self {
        Self { data }
    }

    /// Build directly from `x` / `y` slices.
    ///
    /// # Errors
    ///
    /// Propagates any [`DataPoints::new`] validation error.
    pub fn from_points(xs: &[f64], ys: &[f64]) -> Result<Self, InterpolationError> {
        Ok(Self::new(DataPoints::new(xs, ys)?))
    }

    /// The underlying knots.
    pub fn data(&self) -> &DataPoints {
        &self.data
    }

    /// Evaluate the interpolant at `at`.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolationError::OutOfDomain`] if `at` lies
    /// outside `[x_first, x_last]` (extrapolation is rejected).
    pub fn value_at(&self, at: f64) -> Result<f64, InterpolationError> {
        self.data.check_domain(at)?;
        let xs = self.data.xs();
        let ys = self.data.ys();
        if xs.len() == 1 {
            return Ok(ys[0]);
        }
        let i = self.data.segment_index(at);
        let x0 = xs[i];
        let x1 = xs[i + 1];
        let y0 = ys[i];
        let y1 = ys[i + 1];
        let t = (at - x0) / (x1 - x0);
        Ok(y0 + t * (y1 - y0))
    }

    /// Slope of the segment containing `at`.
    ///
    /// The linear interpolant is piecewise-affine, so its derivative
    /// is the constant slope of the bracketing segment. At an interior
    /// knot this returns the slope of the segment to the right.
    ///
    /// # Errors
    ///
    /// Returns [`InterpolationError::OutOfDomain`] if `at` lies
    /// outside the domain, or [`InterpolationError::TooFewPoints`] if
    /// there is only a single knot (no segment, hence no slope).
    pub fn slope_at(&self, at: f64) -> Result<f64, InterpolationError> {
        self.data.check_domain(at)?;
        let xs = self.data.xs();
        let ys = self.data.ys();
        if xs.len() == 1 {
            return Err(InterpolationError::TooFewPoints { needed: 2, got: 1 });
        }
        let i = self.data.segment_index(at);
        Ok((ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < EPS, "expected {a} ~= {b}");
    }

    #[test]
    fn passes_through_every_knot() {
        let xs = [0.0, 1.0, 2.5, 4.0];
        let ys = [3.0, -1.0, 2.0, 7.0];
        let li = LinearInterpolator::from_points(&xs, &ys).unwrap();
        for (&x, &y) in xs.iter().zip(ys.iter()) {
            close(li.value_at(x).unwrap(), y);
        }
    }

    #[test]
    fn segment_is_the_connecting_line() {
        // Knots (0,0) and (2,10): midpoint must be (1, 5), quarter (0.5, 2.5).
        let li = LinearInterpolator::from_points(&[0.0, 2.0], &[0.0, 10.0]).unwrap();
        close(li.value_at(1.0).unwrap(), 5.0);
        close(li.value_at(0.5).unwrap(), 2.5);
        close(li.value_at(1.5).unwrap(), 7.5);
    }

    #[test]
    fn matches_independent_line_formula_across_a_segment() {
        // Segment from (1,2) to (3,-4): y = 2 + ((-4-2)/(3-1))*(x-1).
        let li = LinearInterpolator::from_points(&[1.0, 3.0], &[2.0, -4.0]).unwrap();
        for k in 0..=10 {
            let x = 1.0 + 2.0 * (k as f64 / 10.0);
            let expected = 2.0 + (-3.0) * (x - 1.0);
            close(li.value_at(x).unwrap(), expected);
        }
    }

    #[test]
    fn constant_data_gives_constant_output() {
        let li = LinearInterpolator::from_points(&[0.0, 1.0, 2.0, 3.0], &[5.0; 4]).unwrap();
        for k in 0..=12 {
            let x = 3.0 * (k as f64 / 12.0);
            close(li.value_at(x).unwrap(), 5.0);
        }
    }

    #[test]
    fn slope_is_segment_slope() {
        let li = LinearInterpolator::from_points(&[0.0, 2.0, 3.0], &[0.0, 10.0, 10.0]).unwrap();
        close(li.slope_at(1.0).unwrap(), 5.0); // first segment slope
        close(li.slope_at(2.5).unwrap(), 0.0); // flat second segment
    }

    #[test]
    fn single_point_is_constant() {
        let li = LinearInterpolator::from_points(&[4.0], &[9.0]).unwrap();
        close(li.value_at(4.0).unwrap(), 9.0);
        assert!(li.slope_at(4.0).is_err());
    }

    #[test]
    fn rejects_extrapolation() {
        let li = LinearInterpolator::from_points(&[0.0, 1.0], &[0.0, 1.0]).unwrap();
        assert_eq!(
            li.value_at(2.0).unwrap_err().code(),
            "interpolation.out_of_domain"
        );
    }
}
