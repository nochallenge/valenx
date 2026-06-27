//! A polynomial response-surface surrogate.
//!
//! When the underlying model `f` is expensive, a *surrogate* (metamodel) is a
//! cheap function `f̂` fitted to a modest set of `(input, output)` samples and
//! evaluated in `f`'s place. This module provides the simplest robust choice:
//! a **multivariate polynomial response surface** of configurable low degree
//! (`0`, `1`, or `2`), fitted by ordinary least squares.
//!
//! The basis is the full multivariate polynomial up to the chosen degree:
//!
//! * degree 0 — `{1}` (a constant);
//! * degree 1 — `{1, x_1, …, x_d}` (a hyperplane);
//! * degree 2 — adds the squares `x_i²` and all distinct cross terms
//!   `x_i·x_j` (`i < j`) — a full quadratic with interactions.
//!
//! The coefficients minimise `‖X β − y‖²` for the design matrix `X` of basis
//! evaluations; the system is solved by SVD (via `nalgebra`), which returns
//! the minimum-norm least-squares solution and is robust to mild
//! ill-conditioning. Goodness of fit is reported by the coefficient of
//! determination `R²`.
//!
//! ## Honesty note — low-order polynomial, not GP/kriging
//!
//! A polynomial response surface captures **smooth global trends**. It does
//! **not** interpolate the training points (unless the data are genuinely a
//! polynomial of the fitted degree), it extrapolates poorly outside the sampled
//! box, and a fixed low degree cannot represent sharp local features or strong
//! non-linearity. A **Gaussian-process / kriging** surrogate — which
//! interpolates the data and supplies a predictive variance — is the natural
//! next step and is a documented future extension; it is **not** implemented
//! here.

use crate::error::UqError;
use nalgebra::{DMatrix, DVector};

/// Maximum polynomial degree this surrogate supports.
pub const MAX_DEGREE: usize = 2;

/// A fitted polynomial response-surface surrogate over `d` inputs.
///
/// Construct one with [`PolynomialSurrogate::fit`], then call
/// [`PolynomialSurrogate::predict`] to evaluate it. [`PolynomialSurrogate::r_squared`]
/// reports the in-sample goodness of fit.
#[derive(Debug, Clone, PartialEq)]
pub struct PolynomialSurrogate {
    n_inputs: usize,
    degree: usize,
    /// Least-squares coefficients, aligned with the basis term order produced
    /// by [`Self::basis`].
    coefficients: Vec<f64>,
    /// In-sample coefficient of determination, cached at fit time.
    r_squared: f64,
}

impl PolynomialSurrogate {
    /// Fit a degree-`degree` polynomial surrogate to `samples` → `values` by
    /// least squares.
    ///
    /// `samples` is a slice of input vectors (all the same length `d`), and
    /// `values[k]` is the scalar response observed at `samples[k]`.
    ///
    /// # Errors
    /// * [`UqError::OutOfRange`] if `degree > MAX_DEGREE`.
    /// * [`UqError::EmptyInput`] if `samples` is empty.
    /// * [`UqError::DimensionMismatch`] if `samples.len() != values.len()`, the
    ///   sample rows are ragged, or there are fewer samples than basis terms
    ///   (an under-determined fit).
    /// * [`UqError::LinearAlgebra`] if the SVD of the design matrix fails.
    pub fn fit(samples: &[Vec<f64>], values: &[f64], degree: usize) -> Result<Self, UqError> {
        if degree > MAX_DEGREE {
            return Err(UqError::OutOfRange(format!(
                "polynomial degree must be <= {MAX_DEGREE} (got {degree})"
            )));
        }
        if samples.is_empty() {
            return Err(UqError::EmptyInput("surrogate fit with no samples".into()));
        }
        if samples.len() != values.len() {
            return Err(UqError::DimensionMismatch(format!(
                "samples ({}) and values ({}) must have equal length",
                samples.len(),
                values.len()
            )));
        }

        let n_inputs = samples[0].len();
        if let Some(bad) = samples.iter().position(|s| s.len() != n_inputs) {
            return Err(UqError::DimensionMismatch(format!(
                "sample row {bad} has length {} but row 0 has length {n_inputs}",
                samples[bad].len()
            )));
        }

        let n_terms = n_basis_terms(n_inputs, degree);
        let n_rows = samples.len();
        if n_rows < n_terms {
            return Err(UqError::DimensionMismatch(format!(
                "need at least {n_terms} samples to fit a degree-{degree} polynomial \
                 in {n_inputs} inputs (got {n_rows})"
            )));
        }

        // Design matrix X (n_rows × n_terms) and response vector y.
        let mut x = DMatrix::<f64>::zeros(n_rows, n_terms);
        for (r, sample) in samples.iter().enumerate() {
            let row = basis(sample, degree);
            for (c, term) in row.into_iter().enumerate() {
                x[(r, c)] = term;
            }
        }
        let y = DVector::<f64>::from_row_slice(values);

        // Minimum-norm least-squares solution via SVD: β = X⁺ y.
        let svd = x.clone().svd(true, true);
        let beta = svd
            .solve(&y, 1e-12)
            .map_err(|e| UqError::LinearAlgebra(format!("SVD least-squares solve failed: {e}")))?;
        let coefficients: Vec<f64> = beta.iter().copied().collect();

        // Coefficient of determination R² = 1 - SS_res / SS_tot.
        let predictions = &x * &beta;
        let ss_res: f64 = values
            .iter()
            .zip(predictions.iter())
            .map(|(&obs, &pred)| (obs - pred) * (obs - pred))
            .sum();
        let mean_y: f64 = values.iter().sum::<f64>() / n_rows as f64;
        let ss_tot: f64 = values.iter().map(|&v| (v - mean_y) * (v - mean_y)).sum();
        // If the response is constant (SS_tot == 0) the fit is exact iff the
        // residual is also zero; report 1.0 in that perfectly-explained case
        // and 0.0 otherwise rather than dividing by zero.
        let r_squared = if ss_tot > 0.0 {
            1.0 - ss_res / ss_tot
        } else if ss_res <= f64::EPSILON {
            1.0
        } else {
            0.0
        };

        Ok(Self {
            n_inputs,
            degree,
            coefficients,
            r_squared,
        })
    }

    /// Predict the response at input point `x`.
    ///
    /// `x` must have length [`Self::n_inputs`]; a wrong length yields
    /// [`UqError::DimensionMismatch`].
    ///
    /// # Errors
    /// [`UqError::DimensionMismatch`] if `x.len() != self.n_inputs()`.
    pub fn predict(&self, x: &[f64]) -> Result<f64, UqError> {
        if x.len() != self.n_inputs {
            return Err(UqError::DimensionMismatch(format!(
                "predict expected {} inputs, got {}",
                self.n_inputs,
                x.len()
            )));
        }
        let row = basis(x, self.degree);
        Ok(row
            .iter()
            .zip(self.coefficients.iter())
            .map(|(t, c)| t * c)
            .sum())
    }

    /// The in-sample coefficient of determination `R² ∈ (-∞, 1]`.
    ///
    /// `R² = 1` is a perfect fit; `R²` near 1 is a good fit; lower (or
    /// negative) values mean the polynomial explains little of the response
    /// variance. Because it is computed on the *training* sample it measures
    /// fit, not predictive accuracy.
    #[must_use]
    pub fn r_squared(&self) -> f64 {
        self.r_squared
    }

    /// The number of input dimensions the surrogate was fitted over.
    #[must_use]
    pub fn n_inputs(&self) -> usize {
        self.n_inputs
    }

    /// The polynomial degree of the surrogate.
    #[must_use]
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// The fitted least-squares coefficients, in basis-term order
    /// (constant, then linear, then — for degree 2 — squares and cross terms).
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }
}

/// Number of basis terms for a degree-`degree` polynomial in `d` inputs.
fn n_basis_terms(d: usize, degree: usize) -> usize {
    match degree {
        0 => 1,
        1 => 1 + d,
        // 1 (const) + d (linear) + d (squares) + d(d-1)/2 (cross terms).
        _ => 1 + 2 * d + d * d.saturating_sub(1) / 2,
    }
}

/// Evaluate the polynomial basis at point `x` for the given `degree`, in a
/// fixed canonical order: `1`, then `x_i`, then (degree 2) `x_i²`, then the
/// cross terms `x_i·x_j` for `i < j`.
fn basis(x: &[f64], degree: usize) -> Vec<f64> {
    let d = x.len();
    let mut terms = Vec::with_capacity(n_basis_terms(d, degree));
    terms.push(1.0); // constant
    if degree >= 1 {
        terms.extend_from_slice(x); // linear
    }
    if degree >= 2 {
        for &xi in x {
            terms.push(xi * xi); // squares
        }
        for i in 0..d {
            for j in (i + 1)..d {
                terms.push(x[i] * x[j]); // cross terms
            }
        }
    }
    terms
}
