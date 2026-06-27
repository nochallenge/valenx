//! Polynomial Chaos Expansion (PCE) — a spectral surrogate for low dimension.
//!
//! Polynomial chaos represents a model's output as a series in **orthogonal
//! polynomials of the input random variables**:
//!
//! ```text
//!   y = f(ξ) ≈ Σ_α  c_α · Ψ_α(ξ)
//! ```
//!
//! where each `Ψ_α` is a (tensor-product) orthogonal polynomial chosen to match
//! the input distribution (the *Wiener–Askey* scheme): **probabilists' Hermite**
//! polynomials for `N(0, 1)` inputs and **Legendre** polynomials for
//! `Uniform(−1, 1)` inputs. Once the coefficients `c_α` are known, the output
//! moments fall out *analytically* from orthogonality — no extra sampling:
//!
//! * **mean** `= c_0` (the coefficient of the constant basis term `Ψ_0 ≡ 1`);
//! * **variance** `= Σ_{α ≠ 0} c_α² · ‖Ψ_α‖²`,
//!
//! with the squared norms `‖Ψ_α‖² = E[Ψ_α(ξ)²]` known in closed form
//! (`n!` for Hermite, `1/(2n+1)` for Legendre).
//!
//! ## What this module fits
//!
//! A **single-input** (univariate) PCE of configurable degree, fitted to
//! `(input, output)` pairs by least-squares regression
//! ([`Pce::fit_regression`]). For an output that is genuinely a polynomial of
//! the input up to the fitted degree, regression recovers the exact
//! coefficients (to round-off) and the analytic mean/variance are exact.
//!
//! ## Honesty / scope caveats — the curse of dimensionality
//!
//! * PCE is provided here for **one input dimension**. A full multivariate PCE
//!   uses the tensor / total-degree product basis, whose term count grows as
//!   `C(d + p, p)` for `d` inputs and degree `p` — the **curse of
//!   dimensionality** that makes naive PCE impractical beyond a handful of
//!   inputs without sparse/adaptive truncation. That multivariate extension is
//!   deliberately out of scope; for many inputs use the global
//!   [`crate::surrogate::PolynomialSurrogate`] or sampling-based
//!   [`crate::statistics`] instead.
//! * The expansion is only as good as its **truncation**: if the true response
//!   has spectral content above the fitted degree, the omitted terms bias both
//!   the surrogate and the variance (which is then an *under*-estimate, since
//!   every dropped `c_α²‖Ψ_α‖²` is non-negative).
//! * Hermite assumes the input is **standard** normal `N(0, 1)` and Legendre
//!   assumes **`Uniform(−1, 1)`**. Other parameters must be mapped to the
//!   standard variable by the caller (`ξ = (x − μ)/σ`, or
//!   `ξ = 2(x − lo)/(hi − lo) − 1`) before fitting; the recovered moments are
//!   then in the standardised variable's terms unless the response itself was
//!   expressed in `x`.

use crate::error::UqError;
use nalgebra::{DMatrix, DVector};

/// Which orthogonal-polynomial family the expansion uses, chosen to match the
/// input distribution under the Wiener–Askey scheme.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PolyBasis {
    /// Probabilists' Hermite polynomials `He_n`. Orthogonal w.r.t. the standard
    /// normal weight: `E[He_m He_n] = n! · δ_mn` for `ξ ~ N(0, 1)`. Use for
    /// **normal** inputs (standardise to `N(0, 1)` first).
    Hermite,
    /// Legendre polynomials `P_n`. Orthogonal w.r.t. the uniform weight on
    /// `[−1, 1]`: `E[P_m P_n] = δ_mn / (2n + 1)` for `ξ ~ Uniform(−1, 1)`. Use
    /// for **uniform** inputs (map to `[−1, 1]` first).
    Legendre,
}

impl PolyBasis {
    /// Evaluate the basis polynomial of degree `n` at `x` via its three-term
    /// recurrence.
    ///
    /// * Hermite (probabilists'): `He_0 = 1`, `He_1 = x`,
    ///   `He_{n+1} = x·He_n − n·He_{n−1}`.
    /// * Legendre: `P_0 = 1`, `P_1 = x`,
    ///   `(n+1)·P_{n+1} = (2n+1)·x·P_n − n·P_{n−1}`.
    #[must_use]
    pub fn eval_poly(self, n: usize, x: f64) -> f64 {
        match self {
            PolyBasis::Hermite => {
                if n == 0 {
                    return 1.0;
                }
                if n == 1 {
                    return x;
                }
                let mut p0 = 1.0; // He_0
                let mut p1 = x; // He_1
                for k in 1..n {
                    let p2 = x * p1 - (k as f64) * p0;
                    p0 = p1;
                    p1 = p2;
                }
                p1
            }
            PolyBasis::Legendre => {
                if n == 0 {
                    return 1.0;
                }
                if n == 1 {
                    return x;
                }
                let mut p0 = 1.0; // P_0
                let mut p1 = x; // P_1
                for k in 1..n {
                    let kf = k as f64;
                    let p2 = ((2.0 * kf + 1.0) * x * p1 - kf * p0) / (kf + 1.0);
                    p0 = p1;
                    p1 = p2;
                }
                p1
            }
        }
    }

    /// The squared norm `‖Ψ_n‖² = E[Ψ_n(ξ)²]` under the matching weight.
    ///
    /// * Hermite: `n!`.
    /// * Legendre: `1 / (2n + 1)`.
    #[must_use]
    pub fn squared_norm(self, n: usize) -> f64 {
        match self {
            PolyBasis::Hermite => factorial(n),
            PolyBasis::Legendre => 1.0 / (2.0 * n as f64 + 1.0),
        }
    }
}

/// A fitted univariate Polynomial Chaos Expansion.
///
/// Build one with [`Pce::fit_regression`], then read [`Pce::mean`] /
/// [`Pce::variance`] for the analytic output moments or [`Pce::predict`] to
/// evaluate the surrogate at a (standardised) input.
#[derive(Debug, Clone, PartialEq)]
pub struct Pce {
    basis: PolyBasis,
    /// Coefficients `c_0 … c_degree`, index = polynomial degree.
    coefficients: Vec<f64>,
}

impl Pce {
    /// Fit a degree-`degree` univariate PCE to `(inputs, outputs)` by
    /// least-squares regression.
    ///
    /// `inputs[k]` is the **standardised** input value `ξ_k` (already mapped to
    /// `N(0, 1)` for [`PolyBasis::Hermite`] or to `[−1, 1]` for
    /// [`PolyBasis::Legendre`]) and `outputs[k]` is the observed response. The
    /// design matrix has column `j` equal to `Ψ_j(ξ_k)`; the normal equations
    /// are solved by SVD (minimum-norm least squares), matching the surrogate
    /// module.
    ///
    /// # Errors
    /// * [`UqError::EmptyInput`] if `inputs` is empty.
    /// * [`UqError::DimensionMismatch`] if `inputs.len() != outputs.len()` or
    ///   there are fewer samples than the `degree + 1` basis terms.
    /// * [`UqError::LinearAlgebra`] if the SVD solve fails.
    pub fn fit_regression(
        inputs: &[f64],
        outputs: &[f64],
        degree: usize,
        basis: PolyBasis,
    ) -> Result<Self, UqError> {
        if inputs.is_empty() {
            return Err(UqError::EmptyInput("PCE fit with no samples".into()));
        }
        if inputs.len() != outputs.len() {
            return Err(UqError::DimensionMismatch(format!(
                "inputs ({}) and outputs ({}) must have equal length",
                inputs.len(),
                outputs.len()
            )));
        }
        let n_terms = degree + 1;
        let n_rows = inputs.len();
        if n_rows < n_terms {
            return Err(UqError::DimensionMismatch(format!(
                "need at least {n_terms} samples to fit a degree-{degree} PCE (got {n_rows})"
            )));
        }

        // Design matrix X[k, j] = Ψ_j(ξ_k).
        let mut x = DMatrix::<f64>::zeros(n_rows, n_terms);
        for (r, &xi) in inputs.iter().enumerate() {
            for c in 0..n_terms {
                x[(r, c)] = basis.eval_poly(c, xi);
            }
        }
        let y = DVector::<f64>::from_row_slice(outputs);

        let svd = x.svd(true, true);
        let beta = svd
            .solve(&y, 1e-12)
            .map_err(|e| UqError::LinearAlgebra(format!("PCE SVD solve failed: {e}")))?;

        Ok(Self {
            basis,
            coefficients: beta.iter().copied().collect(),
        })
    }

    /// The expansion's mean = `c_0` (the constant-mode coefficient).
    #[must_use]
    pub fn mean(&self) -> f64 {
        self.coefficients[0]
    }

    /// The expansion's variance `= Σ_{n ≥ 1} c_n² · ‖Ψ_n‖²`.
    ///
    /// Each summand is non-negative, so a truncated expansion can only
    /// *under*-estimate the true variance (any omitted high-degree mode is
    /// dropped).
    #[must_use]
    pub fn variance(&self) -> f64 {
        self.coefficients
            .iter()
            .enumerate()
            .skip(1)
            .map(|(n, &c)| c * c * self.basis.squared_norm(n))
            .sum()
    }

    /// Evaluate the surrogate at the **standardised** input `xi`.
    #[must_use]
    pub fn predict(&self, xi: f64) -> f64 {
        self.coefficients
            .iter()
            .enumerate()
            .map(|(n, &c)| c * self.basis.eval_poly(n, xi))
            .sum()
    }

    /// The fitted coefficients `c_0 … c_degree`, indexed by polynomial degree.
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// The orthogonal-polynomial family this expansion uses.
    #[must_use]
    pub fn basis(&self) -> PolyBasis {
        self.basis
    }
}

/// `n!` as an `f64` (exact for `n ≤ 18`; PCE degrees are far smaller).
fn factorial(n: usize) -> f64 {
    (1..=n).map(|k| k as f64).product()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::SplitMix64;

    #[test]
    fn hermite_recurrence_values() {
        // He_2 = x²−1, He_3 = x³−3x, He_4 = x⁴−6x²+3.
        let b = PolyBasis::Hermite;
        let x = 2.0;
        assert!((b.eval_poly(0, x) - 1.0).abs() < 1e-12);
        assert!((b.eval_poly(1, x) - 2.0).abs() < 1e-12);
        assert!((b.eval_poly(2, x) - (x * x - 1.0)).abs() < 1e-12);
        assert!((b.eval_poly(3, x) - (x.powi(3) - 3.0 * x)).abs() < 1e-12);
        assert!((b.eval_poly(4, x) - (x.powi(4) - 6.0 * x * x + 3.0)).abs() < 1e-12);
    }

    #[test]
    fn legendre_recurrence_values() {
        // P_2 = (3x²−1)/2, P_3 = (5x³−3x)/2.
        let b = PolyBasis::Legendre;
        let x = 0.3;
        assert!((b.eval_poly(0, x) - 1.0).abs() < 1e-12);
        assert!((b.eval_poly(1, x) - x).abs() < 1e-12);
        assert!((b.eval_poly(2, x) - (3.0 * x * x - 1.0) / 2.0).abs() < 1e-12);
        assert!((b.eval_poly(3, x) - (5.0 * x.powi(3) - 3.0 * x) / 2.0).abs() < 1e-12);
    }

    #[test]
    fn squared_norms() {
        let h = PolyBasis::Hermite;
        assert_eq!(h.squared_norm(0), 1.0);
        assert_eq!(h.squared_norm(1), 1.0);
        assert_eq!(h.squared_norm(2), 2.0);
        assert_eq!(h.squared_norm(3), 6.0);
        let l = PolyBasis::Legendre;
        assert_eq!(l.squared_norm(0), 1.0);
        assert!((l.squared_norm(1) - 1.0 / 3.0).abs() < 1e-15);
        assert!((l.squared_norm(2) - 1.0 / 5.0).abs() < 1e-15);
    }

    /// PCE of a **known polynomial response** with `N(0,1)` inputs (Hermite)
    /// recovers the analytic mean and variance to < 1e-6.
    ///
    /// Response: y = x². In the Hermite basis x² = He_2(x) + 1, so the exact
    /// coefficients are c_0 = 1, c_2 = 1, and the analytic moments are
    ///   mean  = c_0                = 1   (= E[x²]   for x ~ N(0,1)),
    ///   var   = c_2² · 2!          = 2   (= Var[x²] for x ~ N(0,1)).
    #[test]
    fn pce_matches_analytic_moments_quadratic() {
        let mut rng = SplitMix64::new(0x1357_9BDF);
        let n = 2000;
        let inputs: Vec<f64> = (0..n).map(|_| rng.next_standard_normal()).collect();
        let outputs: Vec<f64> = inputs.iter().map(|&x| x * x).collect();

        let pce = Pce::fit_regression(&inputs, &outputs, 2, PolyBasis::Hermite).unwrap();

        // Coefficients: c_0 ≈ 1, c_1 ≈ 0, c_2 ≈ 1.
        assert!((pce.coefficients()[0] - 1.0).abs() < 1e-6);
        assert!(pce.coefficients()[1].abs() < 1e-6);
        assert!((pce.coefficients()[2] - 1.0).abs() < 1e-6);

        assert!((pce.mean() - 1.0).abs() < 1e-6, "mean = {}", pce.mean());
        assert!(
            (pce.variance() - 2.0).abs() < 1e-6,
            "variance = {}",
            pce.variance()
        );
    }

    /// A cubic response y = x³ − 2x exercises a higher-degree Hermite term.
    /// x³ − 2x = He_3(x) + 3x − 2x = He_3(x) + He_1(x), so c_1 = 1, c_3 = 1:
    ///   mean = 0,  var = c_1²·1! + c_3²·3! = 1 + 6 = 7.
    #[test]
    fn pce_matches_analytic_moments_cubic() {
        let mut rng = SplitMix64::new(0x2468_ACE0);
        let n = 4000;
        let inputs: Vec<f64> = (0..n).map(|_| rng.next_standard_normal()).collect();
        let outputs: Vec<f64> = inputs.iter().map(|&x| x.powi(3) - 2.0 * x).collect();

        let pce = Pce::fit_regression(&inputs, &outputs, 3, PolyBasis::Hermite).unwrap();
        assert!((pce.mean() - 0.0).abs() < 1e-6);
        assert!(
            (pce.variance() - 7.0).abs() < 1e-6,
            "var = {}",
            pce.variance()
        );
    }

    /// Legendre PCE on a uniform input. Response y = x² on Uniform(−1,1).
    /// x² = (2/3)·P_2(x) + 1/3, so c_0 = 1/3, c_2 = 2/3:
    ///   mean = 1/3                    (= E[x²]   = 1/3),
    ///   var  = (2/3)²·(1/5) = 4/45   (= Var[x²] = 4/45 for U(−1,1)).
    #[test]
    fn pce_legendre_uniform_moments() {
        let mut rng = SplitMix64::new(0x0F0F_0F0F);
        let n = 4000;
        // Uniform(−1, 1) draws.
        let inputs: Vec<f64> = (0..n).map(|_| 2.0 * rng.next_f64() - 1.0).collect();
        let outputs: Vec<f64> = inputs.iter().map(|&x| x * x).collect();

        let pce = Pce::fit_regression(&inputs, &outputs, 2, PolyBasis::Legendre).unwrap();
        assert!(
            (pce.mean() - 1.0 / 3.0).abs() < 1e-6,
            "mean = {}",
            pce.mean()
        );
        assert!(
            (pce.variance() - 4.0 / 45.0).abs() < 1e-6,
            "var = {}",
            pce.variance()
        );
    }

    #[test]
    fn predict_reproduces_polynomial() {
        let mut rng = SplitMix64::new(42);
        let n = 500;
        let inputs: Vec<f64> = (0..n).map(|_| rng.next_standard_normal()).collect();
        let outputs: Vec<f64> = inputs.iter().map(|&x| 3.0 + x * x).collect();
        let pce = Pce::fit_regression(&inputs, &outputs, 2, PolyBasis::Hermite).unwrap();
        for &xi in &[-1.5, 0.0, 0.7, 2.1] {
            assert!((pce.predict(xi) - (3.0 + xi * xi)).abs() < 1e-6);
        }
    }

    #[test]
    fn fit_rejects_bad_input() {
        assert!(Pce::fit_regression(&[], &[], 1, PolyBasis::Hermite).is_err());
        assert!(Pce::fit_regression(&[1.0, 2.0], &[1.0], 1, PolyBasis::Hermite).is_err());
        // Fewer samples than basis terms (degree 3 → 4 terms, only 2 samples).
        assert!(Pce::fit_regression(&[1.0, 2.0], &[1.0, 2.0], 3, PolyBasis::Hermite).is_err());
    }
}
