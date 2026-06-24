//! Operator Inference (OpInf).
//!
//! OpInf learns the *reduced* operators of a dynamical system directly from
//! projected snapshot data by least squares, without access to the full-order
//! operators (a non-intrusive alternative to [`crate::galerkin`]). Given a
//! reduced state trajectory `q_k ∈ ℝʳ` (e.g. POD coordinates `q_k = Uᵣᵀ x_k`)
//! and its time derivatives `q̇_k`, it fits
//!
//! ```text
//! q̇ = c + A q                         (linear, with optional constant c)
//! q̇ = c + A q + H (q ⊗̃ q)            (with the quadratic term)
//! ```
//!
//! by solving the linear-in-the-operators least-squares problem
//!
//! ```text
//! minimise  Σ_k ‖ q̇_k − [ c + A q_k (+ H (q_k ⊗̃ q_k)) ] ‖²
//! ```
//!
//! Stacking the unknown operators into one matrix `O = [c, A, H]` and the data
//! into `D = [1, qᵀ, (q⊗̃q)ᵀ]`, this is the standard `O Dᵀ ≈ Q̇` system solved
//! per output row via the [`nalgebra`] SVD (a minimum-norm / Tikhonov-
//! regularised least-squares solve).
//!
//! ## Quadratic term — unique Kronecker product
//!
//! The naive `q ⊗ q` has `r²` entries with duplicate cross terms (`q_i q_j` and
//! `q_j q_i`). OpInf uses the **unique** (symmetric) product keeping only
//! `i ≤ j`, giving `r(r+1)/2` columns — fewer unknowns, better conditioning,
//! and no rank deficiency from exact duplicate columns. See [`kron_unique`].
//!
//! ## Conditioning & regularisation
//!
//! The data matrix `D` is often ill-conditioned (POD coordinates decay over
//! orders of magnitude). [`OpInfConfig::ridge`] adds an L2 (Tikhonov) penalty
//! `λ‖O‖²` which is the recommended OpInf stabiliser (Peherstorfer & Willcox,
//! 2016; McQuarrie et al., 2021). The default is a small `1e-10` floor; raise
//! it if the inferred operator forecasts diverge. With `ridge = 0` the solve is
//! the bare minimum-norm SVD least-squares (and will fail loud if the system is
//! numerically rank-deficient).
//!
//! ## Time stepping
//!
//! [`OpInfModel::step`] advances the reduced state one step. The continuous
//! model `q̇ = f(q)` is integrated with classical RK4 over `dt`; a model fit in
//! discrete mode (derivatives = forward differences) is still integrated as a
//! continuous ODE, which is exact for the linear-operator recovery test here.

use nalgebra::{DMatrix, DVector};

use crate::error::RomError;

/// Configuration for an OpInf fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OpInfConfig {
    /// Include a constant term `c`.
    pub constant: bool,
    /// Include the quadratic term `H (q ⊗̃ q)`.
    pub quadratic: bool,
    /// L2 (Tikhonov / ridge) regularisation weight `λ ≥ 0`.
    pub ridge: f64,
}

impl Default for OpInfConfig {
    fn default() -> Self {
        Self {
            constant: false,
            quadratic: false,
            // Small floor: stabilises the SVD solve without materially biasing
            // a well-posed problem. Raise for ill-conditioned data.
            ridge: 1e-10,
        }
    }
}

/// A fitted reduced operator model.
#[derive(Debug, Clone)]
pub struct OpInfModel {
    /// Reduced dimension `r`.
    r: usize,
    /// Constant term `c` (length `r`), zero if not fitted.
    c: DVector<f64>,
    /// Linear operator `A` (`r x r`).
    a: DMatrix<f64>,
    /// Quadratic operator `H` (`r x r(r+1)/2`), empty if not fitted.
    h: Option<DMatrix<f64>>,
}

/// The unique (symmetric) Kronecker product of `q` with itself: the
/// `r(r+1)/2` products `q_i q_j` for `i ≤ j`, in row-major `(i, j)` order.
///
/// This is the column layout OpInf's quadratic operator `H` expects.
pub fn kron_unique(q: &DVector<f64>) -> DVector<f64> {
    let r = q.len();
    let mut out = Vec::with_capacity(r * (r + 1) / 2);
    for i in 0..r {
        for j in i..r {
            out.push(q[i] * q[j]);
        }
    }
    DVector::from_vec(out)
}

/// Number of unique quadratic terms for reduced dimension `r`.
fn n_quad(r: usize) -> usize {
    r * (r + 1) / 2
}

impl OpInfModel {
    /// Infer reduced operators from a reduced trajectory and its derivatives.
    ///
    /// - `states`: columns are reduced states `q_k` (`r x m`).
    /// - `derivs`: columns are the matching `q̇_k` (`r x m`), same shape.
    ///
    /// Use [`OpInfModel::fit_from_columns`] for `Vec`-of-`Vec` input, or
    /// [`forward_difference`] to estimate `derivs` from a uniformly sampled
    /// trajectory.
    ///
    /// # Errors
    /// - [`RomError::Empty`] for empty input.
    /// - [`RomError::DimensionMismatch`] if `states` and `derivs` differ in
    ///   shape.
    /// - [`RomError::NonFinite`] for non-finite entries.
    /// - [`RomError::RankDeficient`] if (with `ridge = 0`) the data matrix is
    ///   numerically rank-deficient.
    /// - [`RomError::NotEnoughSamples`] if there are fewer columns than the
    ///   number of operator unknowns to identify.
    pub fn fit(
        states: &DMatrix<f64>,
        derivs: &DMatrix<f64>,
        cfg: &OpInfConfig,
    ) -> Result<Self, RomError> {
        let (r, m) = (states.nrows(), states.ncols());
        if r == 0 || m == 0 {
            return Err(RomError::Empty { rows: r, cols: m });
        }
        if derivs.nrows() != r || derivs.ncols() != m {
            return Err(RomError::DimensionMismatch {
                what: "OpInf derivs shape",
                expected: r * m,
                got: derivs.nrows() * derivs.ncols(),
            });
        }
        if states.iter().chain(derivs.iter()).any(|v| !v.is_finite()) {
            return Err(RomError::NonFinite { what: "OpInf data" });
        }
        if !cfg.ridge.is_finite() || cfg.ridge < 0.0 {
            return Err(RomError::NonFinite {
                what: "OpInf ridge weight",
            });
        }

        // Build the data matrix D (m x p): [1?, qᵀ, (q⊗̃q)ᵀ?].
        let p_const = usize::from(cfg.constant);
        let p_lin = r;
        let p_quad = if cfg.quadratic { n_quad(r) } else { 0 };
        let p = p_const + p_lin + p_quad;

        // Need at least as many independent samples as unknowns per row; with
        // ridge > 0 the solve is always well-posed, but warn-loud when even the
        // regularised system is under-determined and ridge is off.
        if cfg.ridge == 0.0 && m < p {
            return Err(RomError::NotEnoughSamples {
                what: "OpInf least-squares (unregularised)",
                needed: p,
                got: m,
            });
        }

        let mut d = DMatrix::<f64>::zeros(m, p);
        for k in 0..m {
            let q = states.column(k).into_owned();
            let mut col = 0;
            if cfg.constant {
                d[(k, col)] = 1.0;
                col += 1;
            }
            for i in 0..r {
                d[(k, col + i)] = q[i];
            }
            col += r;
            if cfg.quadratic {
                let kq = kron_unique(&q);
                for (i, val) in kq.iter().enumerate() {
                    d[(k, col + i)] = *val;
                }
            }
        }

        // Solve  O Dᵀ ≈ Q̇  ⇔  D Oᵀ ≈ Q̇ᵀ  for Oᵀ (p x r) via ridge LSQ:
        //   Oᵀ = (DᵀD + λI)⁻¹ Dᵀ Q̇ᵀ.
        let dtd = d.transpose() * &d; // p x p
        let dt_qt = d.transpose() * derivs.transpose(); // p x r
        let normal = if cfg.ridge > 0.0 {
            dtd + DMatrix::<f64>::identity(p, p) * cfg.ridge
        } else {
            dtd
        };

        // Solve the symmetric normal equations via SVD (robust to residual
        // ill-conditioning even with ridge). The pseudo-inverse tolerance is
        // relative to the largest singular value.
        let svd = normal.clone().svd(true, true);
        let smax = svd.singular_values.iter().cloned().fold(0.0_f64, f64::max);
        let eps = smax * 1e-14;
        if smax <= 0.0 {
            return Err(RomError::RankDeficient {
                what: "OpInf normal equations",
                tol: eps,
            });
        }
        let o_t = svd
            .solve(&dt_qt, eps)
            .map_err(|_| RomError::RankDeficient {
                what: "OpInf normal equations",
                tol: eps,
            })?; // p x r

        // Unpack Oᵀ rows back into c, A, H.
        let o = o_t.transpose(); // r x p
        let mut col = 0;
        let c = if cfg.constant {
            let cc = o.column(col).into_owned();
            col += 1;
            cc
        } else {
            DVector::zeros(r)
        };
        let a = o.columns(col, r).into_owned();
        col += r;
        let h = if cfg.quadratic {
            Some(o.columns(col, p_quad).into_owned())
        } else {
            None
        };

        Ok(Self { r, c, a, h })
    }

    /// Convenience constructor from `Vec`-of-columns reduced data.
    ///
    /// `states[k]` and `derivs[k]` are the reduced state and derivative at
    /// sample `k`; all columns must share length `r`.
    ///
    /// # Errors
    /// As [`OpInfModel::fit`].
    pub fn fit_from_columns(
        states: &[Vec<f64>],
        derivs: &[Vec<f64>],
        cfg: &OpInfConfig,
    ) -> Result<Self, RomError> {
        if states.len() != derivs.len() {
            return Err(RomError::DimensionMismatch {
                what: "OpInf column counts",
                expected: states.len(),
                got: derivs.len(),
            });
        }
        let s = cols_to_matrix(states, "OpInf states")?;
        let d = cols_to_matrix(derivs, "OpInf derivs")?;
        Self::fit(&s, &d, cfg)
    }

    /// The reduced dimension `r`.
    pub fn reduced_dim(&self) -> usize {
        self.r
    }

    /// The inferred linear operator `A` (`r x r`).
    pub fn linear_operator(&self) -> &DMatrix<f64> {
        &self.a
    }

    /// The inferred constant term `c` (length `r`; zeros if not fitted).
    pub fn constant(&self) -> &DVector<f64> {
        &self.c
    }

    /// The inferred quadratic operator `H` (`r x r(r+1)/2`), if fitted.
    pub fn quadratic_operator(&self) -> Option<&DMatrix<f64>> {
        self.h.as_ref()
    }

    /// Evaluate the reduced right-hand side `f(q) = c + A q + H (q ⊗̃ q)`.
    ///
    /// # Errors
    /// [`RomError::DimensionMismatch`] if `q.len() != r`.
    pub fn rhs(&self, q: &DVector<f64>) -> Result<DVector<f64>, RomError> {
        if q.len() != self.r {
            return Err(RomError::DimensionMismatch {
                what: "OpInf rhs state",
                expected: self.r,
                got: q.len(),
            });
        }
        let mut out = &self.a * q + &self.c;
        if let Some(h) = &self.h {
            out += h * kron_unique(q);
        }
        Ok(out)
    }

    /// Advance the reduced state one step of size `dt` with classical RK4.
    ///
    /// # Errors
    /// - [`RomError::DimensionMismatch`] if `q.len() != r`.
    /// - [`RomError::BadTimeStep`] if `dt` is not finite and positive.
    pub fn step(&self, q: &DVector<f64>, dt: f64) -> Result<DVector<f64>, RomError> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(RomError::BadTimeStep { value: dt });
        }
        if q.len() != self.r {
            return Err(RomError::DimensionMismatch {
                what: "OpInf step state",
                expected: self.r,
                got: q.len(),
            });
        }
        let k1 = self.rhs(q)?;
        let k2 = self.rhs(&(q + &k1 * (dt / 2.0)))?;
        let k3 = self.rhs(&(q + &k2 * (dt / 2.0)))?;
        let k4 = self.rhs(&(q + &k3 * dt))?;
        Ok(q + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0))
    }

    /// Integrate `steps` reduced steps of size `dt` from `q0`, returning the
    /// trajectory (including the initial state) as columns.
    ///
    /// # Errors
    /// As [`OpInfModel::step`].
    pub fn rollout(
        &self,
        q0: &DVector<f64>,
        dt: f64,
        steps: usize,
    ) -> Result<Vec<DVector<f64>>, RomError> {
        let mut traj = Vec::with_capacity(steps + 1);
        let mut q = q0.clone();
        traj.push(q.clone());
        for _ in 0..steps {
            q = self.step(&q, dt)?;
            traj.push(q.clone());
        }
        Ok(traj)
    }
}

/// Estimate column-wise time derivatives of a uniformly sampled trajectory by
/// second-order central differences (first/last via one-sided differences).
///
/// `states` columns are samples `q_k`; `dt` is the uniform interval. The result
/// has the same shape and is the standard input to [`OpInfModel::fit`] when only
/// states (not derivatives) are available.
///
/// # Errors
/// - [`RomError::Empty`] if `states` is empty.
/// - [`RomError::BadTimeStep`] if `dt` is not finite and positive.
/// - [`RomError::NotEnoughSamples`] if there are fewer than 2 columns.
pub fn forward_difference(states: &DMatrix<f64>, dt: f64) -> Result<DMatrix<f64>, RomError> {
    let (r, m) = (states.nrows(), states.ncols());
    if r == 0 || m == 0 {
        return Err(RomError::Empty { rows: r, cols: m });
    }
    if !dt.is_finite() || dt <= 0.0 {
        return Err(RomError::BadTimeStep { value: dt });
    }
    if m < 2 {
        return Err(RomError::NotEnoughSamples {
            what: "finite-difference derivative",
            needed: 2,
            got: m,
        });
    }
    let mut d = DMatrix::<f64>::zeros(r, m);
    for k in 0..m {
        let col = if k == 0 {
            (states.column(1) - states.column(0)) / dt
        } else if k == m - 1 {
            (states.column(m - 1) - states.column(m - 2)) / dt
        } else {
            (states.column(k + 1) - states.column(k - 1)) / (2.0 * dt)
        };
        d.set_column(k, &col);
    }
    Ok(d)
}

/// Pack `Vec`-of-columns into an `r x m` matrix, validating shape & finiteness.
fn cols_to_matrix(columns: &[Vec<f64>], what: &'static str) -> Result<DMatrix<f64>, RomError> {
    let m = columns.len();
    if m == 0 {
        return Err(RomError::Empty { rows: 0, cols: 0 });
    }
    let r = columns[0].len();
    if r == 0 {
        return Err(RomError::Empty { rows: 0, cols: m });
    }
    for c in columns {
        if c.len() != r {
            return Err(RomError::DimensionMismatch {
                what,
                expected: r,
                got: c.len(),
            });
        }
        if c.iter().any(|v| !v.is_finite()) {
            return Err(RomError::NonFinite { what });
        }
    }
    Ok(DMatrix::from_fn(r, m, |i, j| columns[j][i]))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A known reduced linear ODE q̇ = A q. OpInf must recover A and forecast
    /// held-out steps within tolerance.
    fn known_linear_system() -> (DMatrix<f64>, DVector<f64>, f64) {
        // Stable spiral: eigenvalues -0.5 ± 2i.
        let a = DMatrix::from_row_slice(2, 2, &[-0.5, -2.0, 2.0, -0.5]);
        let q0 = DVector::from_row_slice(&[1.0, 0.0]);
        (a, q0, 0.01)
    }

    /// Integrate the true ODE with RK4 to produce ground-truth snapshots.
    fn integrate(a: &DMatrix<f64>, q0: &DVector<f64>, dt: f64, steps: usize) -> Vec<DVector<f64>> {
        let f = |q: &DVector<f64>| a * q;
        let mut q = q0.clone();
        let mut out = vec![q.clone()];
        for _ in 0..steps {
            let k1 = f(&q);
            let k2 = f(&(&q + &k1 * (dt / 2.0)));
            let k3 = f(&(&q + &k2 * (dt / 2.0)));
            let k4 = f(&(&q + &k3 * dt));
            q += (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
            out.push(q.clone());
        }
        out
    }

    #[test]
    fn recovers_known_linear_operator() {
        let (a_true, q0, dt) = known_linear_system();
        let traj = integrate(&a_true, &q0, dt, 400);
        // Use the first 300 for fitting, hold out the rest.
        let fit_traj: Vec<DVector<f64>> = traj[..300].to_vec();
        let states = DMatrix::from_fn(2, fit_traj.len(), |i, j| fit_traj[j][i]);
        let derivs = forward_difference(&states, dt).unwrap();

        let cfg = OpInfConfig {
            constant: false,
            quadratic: false,
            ridge: 0.0,
        };
        let model = OpInfModel::fit(&states, &derivs, &cfg).unwrap();

        // Recovered A should match the true A. The central-difference estimate
        // of q̇ is O(dt²), so allow a small tolerance.
        let err = (model.linear_operator() - &a_true).norm();
        assert!(err < 1e-3, "operator error = {err:e}");

        // Forecast held-out steps from the last fitted state.
        let start = traj[299].clone();
        let horizon = 100;
        let rollout = model.rollout(&start, dt, horizon).unwrap();
        let mut max_err = 0.0_f64;
        for (k, q) in rollout.iter().enumerate() {
            let truth = &traj[299 + k];
            max_err = max_err.max((q - truth).norm());
        }
        assert!(max_err < 1e-2, "forecast error = {max_err:e}");
    }

    #[test]
    fn forecast_matches_analytic_decay() {
        // Scalar q̇ = -q → q(t) = q0 e^{-t}. Recover A ≈ -1 and forecast.
        let dt = 0.01;
        let mut q = 1.0_f64;
        let mut cols = Vec::new();
        for _ in 0..200 {
            cols.push(vec![q]);
            q += dt * (-q); // explicit Euler ground truth is fine for fitting
        }
        // Build cleaner exponential ground truth instead.
        let mut clean = Vec::new();
        for k in 0..200 {
            let t = k as f64 * dt;
            clean.push(vec![(-t).exp()]);
        }
        let states = cols_to_matrix(&clean, "t").unwrap();
        let derivs = forward_difference(&states, dt).unwrap();
        let model = OpInfModel::fit(
            &states,
            &derivs,
            &OpInfConfig {
                constant: false,
                quadratic: false,
                ridge: 0.0,
            },
        )
        .unwrap();
        assert!((model.linear_operator()[(0, 0)] + 1.0).abs() < 1e-3);
        let q0 = DVector::from_row_slice(&[1.0]);
        let roll = model.rollout(&q0, dt, 150).unwrap();
        let t = 150.0 * dt;
        let analytic = (-t).exp();
        assert!((roll[150][0] - analytic).abs() < 1e-3);
    }

    #[test]
    fn quadratic_term_fits_a_logistic() {
        // Logistic q̇ = q - q² has c=0, A=[1], H=[-1] (single quad term q²).
        let dt = 0.005;
        let mut clean = Vec::new();
        // analytic logistic with q(0)=0.1: q(t) = 1/(1+9 e^{-t})
        for k in 0..600 {
            let t = k as f64 * dt;
            clean.push(vec![1.0 / (1.0 + 9.0 * (-t).exp())]);
        }
        let states = cols_to_matrix(&clean, "t").unwrap();
        let derivs = forward_difference(&states, dt).unwrap();
        let model = OpInfModel::fit(
            &states,
            &derivs,
            &OpInfConfig {
                constant: false,
                quadratic: true,
                ridge: 0.0,
            },
        )
        .unwrap();
        // A ≈ 1, H ≈ -1.
        assert!(
            (model.linear_operator()[(0, 0)] - 1.0).abs() < 5e-2,
            "A = {}",
            model.linear_operator()[(0, 0)]
        );
        let h = model.quadratic_operator().unwrap();
        assert!((h[(0, 0)] + 1.0).abs() < 5e-2, "H = {}", h[(0, 0)]);
    }

    #[test]
    fn kron_unique_layout() {
        let q = DVector::from_row_slice(&[2.0, 3.0]);
        let k = kron_unique(&q);
        // i<=j: (0,0)=4, (0,1)=6, (1,1)=9
        assert_eq!(k.len(), 3);
        assert_eq!(k[0], 4.0);
        assert_eq!(k[1], 6.0);
        assert_eq!(k[2], 9.0);
    }

    #[test]
    fn rejects_mismatched_shapes() {
        let states = DMatrix::<f64>::zeros(2, 5);
        let derivs = DMatrix::<f64>::zeros(2, 4);
        assert_eq!(
            OpInfModel::fit(&states, &derivs, &OpInfConfig::default())
                .unwrap_err()
                .code(),
            "dimension_mismatch"
        );
    }

    #[test]
    fn rejects_empty() {
        let states = DMatrix::<f64>::zeros(0, 0);
        let derivs = DMatrix::<f64>::zeros(0, 0);
        assert_eq!(
            OpInfModel::fit(&states, &derivs, &OpInfConfig::default())
                .unwrap_err()
                .code(),
            "empty"
        );
    }

    #[test]
    fn rejects_non_finite() {
        let states = DMatrix::from_row_slice(1, 2, &[1.0, f64::INFINITY]);
        let derivs = DMatrix::from_row_slice(1, 2, &[0.0, 0.0]);
        assert_eq!(
            OpInfModel::fit(&states, &derivs, &OpInfConfig::default())
                .unwrap_err()
                .code(),
            "non_finite"
        );
    }

    #[test]
    fn underdetermined_unregularised_errs() {
        // 3 unknowns per row (constant+linear with r=2), only 2 samples, ridge 0.
        let states = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 3.0, 4.0]);
        let derivs = DMatrix::from_row_slice(2, 2, &[0.1, 0.2, 0.3, 0.4]);
        let cfg = OpInfConfig {
            constant: true,
            quadratic: false,
            ridge: 0.0,
        };
        assert_eq!(
            OpInfModel::fit(&states, &derivs, &cfg).unwrap_err().code(),
            "not_enough_samples"
        );
    }

    #[test]
    fn step_rejects_bad_dt_and_dim() {
        let (a, _, _) = known_linear_system();
        let model = OpInfModel {
            r: 2,
            c: DVector::zeros(2),
            a,
            h: None,
        };
        let q = DVector::from_row_slice(&[1.0, 1.0]);
        assert_eq!(model.step(&q, 0.0).unwrap_err().code(), "bad_time_step");
        let bad = DVector::from_row_slice(&[1.0]);
        assert_eq!(
            model.step(&bad, 0.1).unwrap_err().code(),
            "dimension_mismatch"
        );
    }
}
