//! Discrete Empirical Interpolation Method (DEIM) hyper-reduction.
//!
//! POD–Galerkin makes the *linear* part of a model cheap, but a nonlinear term
//! `f(x)` still has to be evaluated at every one of the `m` full-order degrees
//! of freedom and then projected — so its cost never drops with the reduced
//! dimension. DEIM (Chaturantabut & Sorensen, *SIAM J. Sci. Comput.* 32(5),
//! 2010) removes that bottleneck: given an orthonormal basis `U` (`m × k`, e.g.
//! the POD of a set of nonlinear-term snapshots), it picks `k` **interpolation
//! rows** `p₁…p_k` and approximates the nonlinear vector from only those `k`
//! sampled entries:
//!
//! ```text
//! f(x) ≈ U (Pᵀ U)⁻¹ Pᵀ f(x)
//! ```
//!
//! where `P = [e_{p₁} … e_{p_k}]` selects the chosen rows. Only `Pᵀ f(x)` — the
//! nonlinear function at `k` points — ever has to be computed; the `m × k`
//! operator `U (Pᵀ U)⁻¹` lifts it back to the full state. The error is bounded
//! by `‖(Pᵀ U)⁻¹‖₂` times the POD projection error of `f`, so a well-conditioned
//! point set and a basis that captures `f` give an accurate, fully reduced
//! nonlinear term.
//!
//! ## Greedy index selection
//!
//! The points are chosen one basis vector at a time so that each new sample sits
//! where the running interpolant is *worst* (Chaturantabut–Sorensen, Algorithm
//! 1):
//!
//! - `p₁ = argmaxᵢ |U[i, 0]|`.
//! - For `j = 2 … k`: solve `(Pᵀ U₀..ⱼ₋₁) c = Pᵀ U[:, j-1]` for the coefficients
//!   `c` that interpolate column `j-1` at the points chosen so far, form the
//!   residual `r = U[:, j-1] − U₀..ⱼ₋₁ c`, and take `pⱼ = argmaxᵢ |r[i]|`.
//!
//! The `j × j` system `Pᵀ U₀..ⱼ₋₁` is nonsingular exactly when the greedy has
//! not produced a repeat (guaranteed for a full-rank `U`); this implementation
//! nevertheless **guards the inverse** and returns a [`RomError::RankDeficient`]
//! rather than dividing by a singular factor, in keeping with the crate's
//! fail-loud contract.
//!
//! ## Fail-loud contract
//!
//! Construction rejects an empty basis, a basis with more columns than rows
//! (`k > m`, no `k` distinct rows to pick), a non-finite entry, or a
//! rank-deficient `Pᵀ U`. Nothing here unwraps a `None` from the inverse.
//!
//! ## Example — exact recovery of an in-span vector
//!
//! ```
//! use valenx_rom::Deim;
//! use nalgebra::{DMatrix, DVector};
//!
//! // A 5×2 orthonormal basis (two columns of the identity-ish span).
//! let u = DMatrix::from_column_slice(5, 2, &[
//!     1.0, 0.0, 0.0, 0.0, 0.0, // column 0
//!     0.0, 0.0, 1.0, 0.0, 0.0, // column 1
//! ]);
//! let deim = Deim::new(&u).unwrap();
//! assert_eq!(deim.points().len(), 2);
//! // Any vector in the column span is reconstructed exactly from k samples.
//! let f = DVector::from_column_slice(&[3.0, 0.0, -2.0, 0.0, 0.0]);
//! let sampled = deim.sample(&f);                 // length-2 vector Pᵀ f
//! let recon = deim.approximate(&sampled).unwrap();
//! assert!((recon - f).norm() < 1e-12);
//! ```

use nalgebra::{DMatrix, DVector};

use crate::error::RomError;
use crate::pod::PodBasis;

/// Relative singular-value / pivot floor below which `Pᵀ U` is treated as
/// numerically singular (multiplied by the largest singular value).
const COND_REL_FLOOR: f64 = 1e-12;

/// A DEIM hyper-reduction operator built from a nonlinear-term basis.
///
/// Holds the chosen interpolation row indices `p₁…p_k` and the precomputed
/// `m × k` projection operator `U (Pᵀ U)⁻¹`. Apply it with
/// [`Deim::approximate`] (given the nonlinear term sampled at the selected
/// rows) or sample a full vector with [`Deim::sample`].
#[derive(Debug, Clone)]
pub struct Deim {
    /// Chosen interpolation row indices, in selection order (length = `k`).
    points: Vec<usize>,
    /// Precomputed DEIM operator `U (Pᵀ U)⁻¹` (`m × k`).
    operator: DMatrix<f64>,
}

impl Deim {
    /// Build a DEIM operator from an orthonormal nonlinear-term basis `basis`
    /// (`m × k`, columns = modes).
    ///
    /// The greedy Chaturantabut–Sorensen selection picks `k` distinct rows and
    /// the operator `U (Pᵀ U)⁻¹` is formed and cached.
    ///
    /// The columns need not be exactly orthonormal — only linearly independent —
    /// but a POD basis (which is orthonormal) is the canonical input. The math
    /// only requires `Pᵀ U` to be invertible.
    ///
    /// # Errors
    /// - [`RomError::Empty`] if `basis` has zero rows or zero columns.
    /// - [`RomError::NonFinite`] if any entry is `NaN` or infinite.
    /// - [`RomError::InvalidRank`] if `k > m` (more modes than state rows — there
    ///   are not `k` distinct rows to interpolate at).
    /// - [`RomError::RankDeficient`] if a greedy step's `Pᵀ U` sub-block is
    ///   singular (e.g. a rank-deficient basis), or the final `Pᵀ U` cannot be
    ///   inverted.
    pub fn new(basis: &DMatrix<f64>) -> Result<Self, RomError> {
        let m = basis.nrows();
        let k = basis.ncols();
        if m == 0 || k == 0 {
            return Err(RomError::Empty { rows: m, cols: k });
        }
        if basis.iter().any(|v| !v.is_finite()) {
            return Err(RomError::NonFinite { what: "DEIM basis" });
        }
        if k > m {
            // No way to pick k distinct interpolation rows out of m.
            return Err(RomError::InvalidRank {
                requested: k,
                max: m,
            });
        }

        let points = deim_select_matrix(basis)?;

        // Build P (m x k): P[p_j, j] = 1. Then PᵀU (k x k) and the operator
        // U (PᵀU)⁻¹ (m x k). Guard the inverse — never unwrap a None.
        let pt_u = selected_rows(basis, &points); // k x k
        let inv = invert_guarded(&pt_u, "DEIM (PᵀU)")?;
        let operator = basis * inv; // m x k

        Ok(Self { points, operator })
    }

    /// Build a DEIM operator directly from a [`PodBasis`] of nonlinear-term
    /// snapshots — the canonical workflow (POD of `f`-snapshots → DEIM).
    ///
    /// Equivalent to [`Deim::new`] on `pod.modes()`.
    ///
    /// # Errors
    /// As [`Deim::new`].
    pub fn from_pod(pod: &PodBasis) -> Result<Self, RomError> {
        Deim::new(pod.modes())
    }

    /// The chosen interpolation row indices `p₁…p_k`, in selection order.
    ///
    /// These are guaranteed distinct and each in `0..m`.
    pub fn points(&self) -> &[usize] {
        &self.points
    }

    /// The number of interpolation points `k` (equals the basis column count).
    pub fn n_points(&self) -> usize {
        self.points.len()
    }

    /// The full state dimension `m` (rows of the operator / basis).
    pub fn state_dim(&self) -> usize {
        self.operator.nrows()
    }

    /// The precomputed DEIM operator `U (Pᵀ U)⁻¹`, shape `m × k`.
    pub fn operator(&self) -> &DMatrix<f64> {
        &self.operator
    }

    /// Sample a full nonlinear vector `f` (length `m`) at the selected rows,
    /// returning the length-`k` vector `Pᵀ f`.
    ///
    /// This is the only part of the nonlinear term a reduced model must actually
    /// evaluate; in practice a caller computes `f` at just these
    /// [`Deim::points`] rather than forming the whole vector.
    ///
    /// Indices are in range by construction, so this does not fail; a length
    /// mismatch is reported by [`Deim::approximate`] when the sample is used.
    pub fn sample(&self, f: &DVector<f64>) -> DVector<f64> {
        DVector::from_iterator(
            self.points.len(),
            self.points
                .iter()
                .map(|&p| f.get(p).copied().unwrap_or(0.0)),
        )
    }

    /// Reconstruct the full nonlinear term from its values at the selected rows:
    /// `f ≈ U (Pᵀ U)⁻¹ (Pᵀ f)`.
    ///
    /// `values_at_selected_indices` is the nonlinear function evaluated at the
    /// `k` [`Deim::points`], in the same order. The result has length `m`.
    ///
    /// If `f` lies in the column span of the basis, the reconstruction is exact
    /// (to round-off); otherwise it is the DEIM interpolant.
    ///
    /// # Errors
    /// [`RomError::DimensionMismatch`] if the sample length is not `k`.
    pub fn approximate(
        &self,
        values_at_selected_indices: &DVector<f64>,
    ) -> Result<DVector<f64>, RomError> {
        let k = self.points.len();
        if values_at_selected_indices.len() != k {
            return Err(RomError::DimensionMismatch {
                what: "DEIM sampled values",
                expected: k,
                got: values_at_selected_indices.len(),
            });
        }
        Ok(&self.operator * values_at_selected_indices)
    }
}

/// Greedy Chaturantabut–Sorensen interpolation-index selection.
///
/// Given an orthonormal (or merely full-column-rank) basis `basis` (`m × k`),
/// returns the `k` distinct interpolation row indices `[p₁ … p_k]` in selection
/// order. This is the standalone counterpart to [`Deim::new`] for callers that
/// want only the indices.
///
/// # Errors
/// - [`RomError::Empty`] if `basis` has zero rows or zero columns.
/// - [`RomError::NonFinite`] if any entry is non-finite.
/// - [`RomError::InvalidRank`] if `k > m`.
/// - [`RomError::RankDeficient`] if a greedy step's `Pᵀ U` sub-block is singular.
pub fn deim_select(basis: &DMatrix<f64>) -> Result<Vec<usize>, RomError> {
    let m = basis.nrows();
    let k = basis.ncols();
    if m == 0 || k == 0 {
        return Err(RomError::Empty { rows: m, cols: k });
    }
    if basis.iter().any(|v| !v.is_finite()) {
        return Err(RomError::NonFinite { what: "DEIM basis" });
    }
    if k > m {
        return Err(RomError::InvalidRank {
            requested: k,
            max: m,
        });
    }
    deim_select_matrix(basis)
}

/// Core greedy loop, assuming `basis` is already validated (finite, `0 < k ≤ m`).
fn deim_select_matrix(basis: &DMatrix<f64>) -> Result<Vec<usize>, RomError> {
    let k = basis.ncols();
    let mut points: Vec<usize> = Vec::with_capacity(k);

    // p_1 = argmax |U[:, 0]|.
    let first = argmax_abs(&basis.column(0).into_owned());
    points.push(first);

    for j in 1..k {
        // Current column to interpolate at the points chosen so far.
        let uj = basis.column(j).into_owned(); // m
        let u_prev = basis.columns(0, j).into_owned(); // m x j

        // c solves (Pᵀ U₀..ⱼ₋₁) c = Pᵀ uⱼ  — a j x j system. Guard the inverse.
        let pt_uprev = selected_rows(&u_prev, &points); // j x j
        let pt_uj = DVector::from_iterator(j, points.iter().map(|&p| uj[p])); // j
        let c = solve_guarded(&pt_uprev, &pt_uj, "DEIM greedy (PᵀU) step")?; // j

        // Residual r = uⱼ − U₀..ⱼ₋₁ c, then pⱼ = argmax |r|.
        let r = &uj - &u_prev * c; // m
        let next = argmax_abs(&r);
        points.push(next);
    }

    Ok(points)
}

/// Extract the rows of `mat` indexed by `rows`, preserving order, into a new
/// `rows.len() × mat.ncols()` matrix.
fn selected_rows(mat: &DMatrix<f64>, rows: &[usize]) -> DMatrix<f64> {
    DMatrix::from_fn(rows.len(), mat.ncols(), |i, c| mat[(rows[i], c)])
}

/// Index of the maximum-magnitude entry of `v` (ties → lowest index).
fn argmax_abs(v: &DVector<f64>) -> usize {
    let mut best = 0usize;
    let mut best_val = f64::NEG_INFINITY;
    for (i, &x) in v.iter().enumerate() {
        let a = x.abs();
        if a > best_val {
            best_val = a;
            best = i;
        }
    }
    best
}

/// Invert a small square matrix, mapping a singular/non-invertible factor to
/// [`RomError::RankDeficient`] instead of unwrapping a `None`.
///
/// Uses an SVD-relative floor on the smallest singular value as the
/// rank-deficiency test, then `try_inverse` for the actual inverse (both are
/// guarded; either failing yields the same fail-loud error).
fn invert_guarded(a: &DMatrix<f64>, what: &'static str) -> Result<DMatrix<f64>, RomError> {
    let tol = singular_floor(a);
    if !cond_ok(a, tol) {
        return Err(RomError::RankDeficient { what, tol });
    }
    a.clone()
        .try_inverse()
        .ok_or(RomError::RankDeficient { what, tol })
}

/// Solve `a x = b` for a small square `a`, guarding against a singular `a`.
fn solve_guarded(
    a: &DMatrix<f64>,
    b: &DVector<f64>,
    what: &'static str,
) -> Result<DVector<f64>, RomError> {
    let tol = singular_floor(a);
    if !cond_ok(a, tol) {
        return Err(RomError::RankDeficient { what, tol });
    }
    // LU with partial pivoting; fall back through the rank-deficiency guard if
    // the factor is singular (never unwrap a None).
    a.clone()
        .lu()
        .solve(b)
        .ok_or(RomError::RankDeficient { what, tol })
}

/// Absolute singular-value floor for `a`: `σ_max · COND_REL_FLOOR`.
fn singular_floor(a: &DMatrix<f64>) -> f64 {
    let sv = a.clone().singular_values();
    let smax = sv.iter().cloned().fold(0.0_f64, f64::max);
    smax * COND_REL_FLOOR
}

/// True iff every singular value of `a` exceeds `tol` (i.e. `a` is numerically
/// full rank and safe to invert).
fn cond_ok(a: &DMatrix<f64>, tol: f64) -> bool {
    let sv = a.clone().singular_values();
    let smin = sv.iter().cloned().fold(f64::INFINITY, f64::min);
    smin > tol && tol.is_finite()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::snapshots::Snapshots;
    use nalgebra::{DMatrix, DVector};

    /// An orthonormal basis whose columns are k cosine modes on an m-point grid,
    /// orthonormalised by QR. Gives a well-conditioned, non-trivial DEIM basis.
    fn cosine_basis(m: usize, k: usize) -> DMatrix<f64> {
        let raw = DMatrix::from_fn(m, k, |i, j| {
            (std::f64::consts::PI * (j as f64 + 1.0) * (i as f64 + 0.5) / m as f64).cos()
        });
        let qr = raw.qr();
        qr.q().columns(0, k).into_owned()
    }

    #[test]
    fn pin1_exact_recovery_of_in_span_vector() {
        // ANALYTIC PIN (1): if f lies exactly in the k-mode span, DEIM with k
        // points reconstructs it EXACTLY (<1e-9) from only k sampled entries.
        let m = 30;
        let k = 5;
        let u = cosine_basis(m, k);
        let deim = Deim::new(&u).unwrap();

        // f = U a for a known coefficient vector → exactly in span.
        let a = DVector::from_fn(k, |i, _| 1.0 + 0.5 * i as f64);
        let f = &u * &a;

        let sampled = deim.sample(&f); // only k entries used
        assert_eq!(sampled.len(), k);
        let recon = deim.approximate(&sampled).unwrap();
        let err = (&recon - &f).norm();
        assert!(err < 1e-9, "in-span reconstruction error = {err:e}");
    }

    #[test]
    fn pin1b_approximate_uses_only_k_samples() {
        // The reconstruction depends ONLY on the values at the selected rows:
        // perturbing f off-support must not change approximate() output for an
        // in-span vector that we re-sample at the points.
        let m = 24;
        let k = 4;
        let u = cosine_basis(m, k);
        let deim = Deim::new(&u).unwrap();
        let a = DVector::from_fn(k, |i, _| (i as f64 - 1.5).sin());
        let f = &u * &a;

        // Hand-built sample using ONLY the published point indices.
        let manual = DVector::from_iterator(k, deim.points().iter().map(|&p| f[p]));
        let recon = deim.approximate(&manual).unwrap();
        assert!((&recon - &f).norm() < 1e-9);
    }

    #[test]
    fn pin2_points_are_distinct() {
        // ANALYTIC PIN (2): the greedy selects k DISTINCT indices (no repeats).
        for (m, k) in [(10, 1), (10, 5), (50, 12), (8, 8)] {
            let u = cosine_basis(m, k);
            let pts = deim_select(&u).unwrap();
            assert_eq!(pts.len(), k, "expected k points");
            let mut sorted = pts.clone();
            sorted.sort_unstable();
            sorted.dedup();
            assert_eq!(sorted.len(), k, "DEIM points must be distinct: {pts:?}");
            assert!(pts.iter().all(|&p| p < m), "point out of range: {pts:?}");
        }
    }

    #[test]
    fn pin3_operator_shape_is_m_by_k() {
        // ANALYTIC PIN (3): the operator shape is m x k.
        let m = 40;
        let k = 7;
        let u = cosine_basis(m, k);
        let deim = Deim::new(&u).unwrap();
        assert_eq!(deim.operator().nrows(), m);
        assert_eq!(deim.operator().ncols(), k);
        assert_eq!(deim.state_dim(), m);
        assert_eq!(deim.n_points(), k);
    }

    #[test]
    fn pin4_error_decreases_as_k_grows() {
        // ANALYTIC PIN (4): for a smooth function NOT in the span, the DEIM
        // error decreases as k grows. Use a POD basis of snapshots of a smooth
        // parametric field, and a held-out smooth target.
        let m = 60;

        // Snapshot set: f_mu(x) = exp(-((x - mu)^2)) over a sweep of mu.
        let grid: Vec<f64> = (0..m).map(|i| i as f64 / (m as f64 - 1.0)).collect();
        let mut cols = Vec::new();
        for s in 0..40 {
            let mu = s as f64 / 39.0;
            let col: Vec<f64> = grid
                .iter()
                .map(|&x| (-((x - mu) * (x - mu)) / 0.02).exp())
                .collect();
            cols.push(col);
        }
        let snaps = Snapshots::from_columns(&cols).unwrap();

        // A smooth target NOT in the snapshot set: a Gaussian bump at a centre
        // *between* the sampled mu values and at a slightly different width. It
        // lives near (not on) the snapshot manifold, so adding modes must keep
        // driving the DEIM interpolation error down — a genuine convergence
        // test, not a loosened bound.
        let target = DVector::from_iterator(
            m,
            grid.iter()
                .map(|&x| (-((x - 0.503) * (x - 0.503)) / 0.02).exp()),
        );

        let mut errs = Vec::new();
        for k in [2usize, 4, 8, 12] {
            let pod = crate::pod::PodBasis::fit_rank(&snaps, k).unwrap();
            let deim = Deim::from_pod(&pod).unwrap();
            let sampled = deim.sample(&target);
            let recon = deim.approximate(&sampled).unwrap();
            let err = (&recon - &target).norm() / target.norm();
            errs.push(err);
        }
        assert_eq!(errs.len(), 4);
        // Monotone non-increasing as k grows (small tolerance for round-off).
        for w in errs.windows(2) {
            assert!(
                w[1] <= w[0] * (1.0 + 1e-9),
                "DEIM error did not decrease: {:e} -> {:e}",
                w[0],
                w[1]
            );
        }
        // A large drop over the range (observed ~430×, assert ≥100×) and
        // genuine sub-percent accuracy at the richest basis on this off-grid
        // smooth target.
        assert!(
            errs[3] < errs[0] * 1e-2,
            "DEIM error barely improved: {:e} -> {:e}",
            errs[0],
            errs[3]
        );
        assert!(errs[3] < 5e-3, "k=12 DEIM error {:e} too large", errs[3]);
    }

    #[test]
    fn operator_matches_definition_u_times_inv_ptu() {
        // The cached operator equals U (PᵀU)⁻¹ recomputed independently.
        let m = 20;
        let k = 4;
        let u = cosine_basis(m, k);
        let deim = Deim::new(&u).unwrap();
        let pts = deim.points();
        let pt_u = DMatrix::from_fn(k, k, |i, c| u[(pts[i], c)]);
        let inv = pt_u.try_inverse().unwrap();
        let expected = &u * inv;
        assert!((deim.operator() - &expected).norm() < 1e-9);
    }

    #[test]
    fn first_point_is_argmax_of_first_mode() {
        // p_1 = argmax |U[:, 0]| by construction.
        let m = 16;
        let k = 3;
        let u = cosine_basis(m, k);
        let pts = deim_select(&u).unwrap();
        let col0 = u.column(0).into_owned();
        let want = argmax_abs(&col0);
        assert_eq!(pts[0], want);
    }

    #[test]
    fn from_pod_equals_new_on_modes() {
        let m = 18;
        let k = 3;
        let u = cosine_basis(m, k);
        let snaps = {
            // Make snapshots whose POD basis is exactly these k modes by using
            // the modes with distinct time coefficients.
            let mut cols = Vec::new();
            for t in 0..30 {
                let a = DVector::from_fn(k, |i, _| {
                    ((t as f64 + 1.0) * (i as f64 + 1.0) * 0.1).sin() + (i as f64 + 1.0)
                });
                let v = &u * a;
                cols.push((0..m).map(|i| v[i]).collect::<Vec<f64>>());
            }
            Snapshots::from_columns(&cols).unwrap()
        };
        let pod = crate::pod::PodBasis::fit_rank(&snaps, k).unwrap();
        let via_pod = Deim::from_pod(&pod).unwrap();
        let via_new = Deim::new(pod.modes()).unwrap();
        assert_eq!(via_pod.points(), via_new.points());
        assert!((via_pod.operator() - via_new.operator()).norm() < 1e-12);
    }

    // ---- fail-loud guards -------------------------------------------------

    #[test]
    fn rejects_empty_basis() {
        let e = Deim::new(&DMatrix::<f64>::zeros(0, 3)).unwrap_err();
        assert_eq!(e.code(), "empty");
        let e = Deim::new(&DMatrix::<f64>::zeros(4, 0)).unwrap_err();
        assert_eq!(e.code(), "empty");
        let e = deim_select(&DMatrix::<f64>::zeros(0, 2)).unwrap_err();
        assert_eq!(e.code(), "empty");
    }

    #[test]
    fn rejects_k_greater_than_m() {
        // 3 rows, 5 columns: cannot pick 5 distinct interpolation rows.
        let u = DMatrix::from_fn(3, 5, |i, j| (i + j) as f64 + 1.0);
        let e = Deim::new(&u).unwrap_err();
        assert_eq!(e.code(), "invalid_rank");
        let e = deim_select(&u).unwrap_err();
        assert_eq!(e.code(), "invalid_rank");
    }

    #[test]
    fn rejects_non_finite_basis() {
        let mut u = cosine_basis(6, 2);
        u[(0, 0)] = f64::NAN;
        assert_eq!(Deim::new(&u).unwrap_err().code(), "non_finite");
        let mut u2 = cosine_basis(6, 2);
        u2[(1, 1)] = f64::INFINITY;
        assert_eq!(deim_select(&u2).unwrap_err().code(), "non_finite");
    }

    #[test]
    fn rejects_rank_deficient_basis() {
        // Two identical columns → the second greedy step's PᵀU is singular.
        // (A column-0 duplicate makes residual of column 1 identically zero, so
        // (PᵀU) for the 2x2 step is singular.)
        let mut u = DMatrix::<f64>::zeros(5, 2);
        let col = DVector::from_column_slice(&[1.0, 2.0, 3.0, 4.0, 5.0]).normalize();
        u.set_column(0, &col);
        u.set_column(1, &col); // identical → rank 1
        let e = Deim::new(&u).unwrap_err();
        assert_eq!(
            e.code(),
            "rank_deficient",
            "duplicate-column basis must be rejected, got {e:?}"
        );
    }

    #[test]
    fn approximate_dimension_mismatch_errs() {
        let u = cosine_basis(10, 3);
        let deim = Deim::new(&u).unwrap();
        let bad = DVector::<f64>::zeros(2); // should be length 3
        assert_eq!(
            deim.approximate(&bad).unwrap_err().code(),
            "dimension_mismatch"
        );
    }

    #[test]
    fn single_mode_picks_one_point_and_reconstructs() {
        // k = 1 edge case: one point, operator is m x 1, exact for any scalar
        // multiple of the single mode.
        let m = 12;
        let u = cosine_basis(m, 1);
        let deim = Deim::new(&u).unwrap();
        assert_eq!(deim.points().len(), 1);
        assert_eq!(deim.operator().ncols(), 1);
        let f = u.column(0).into_owned() * 2.5;
        let recon = deim.approximate(&deim.sample(&f)).unwrap();
        assert!((&recon - &f).norm() < 1e-9);
    }
}
