//! Verification & Validation (V&V) harness — order-of-accuracy, grid
//! convergence, manufactured solutions, and error norms.
//!
//! *Verification* asks "are we solving the equations right?" — i.e. is the
//! discrete solver converging to the exact solution of its governing equations
//! at the rate the scheme promises? This module gathers the standard tools for
//! that question, kept deliberately small and analytic so each can be pinned to
//! a textbook value:
//!
//! * [`ConvergenceStudy`] — fit the **observed order of accuracy** as the
//!   least-squares slope of `log|error|` vs `log h` over a refinement sweep.
//! * [`ConvergenceStudy::gci`] — Roache's **Grid Convergence Index**, a
//!   standardised, safety-factored error band from three systematically-refined
//!   grids.
//! * [`Mms`] — a tiny **Method of Manufactured Solutions** framework: pick an
//!   analytic solution, apply the differential operator, and the leftover is
//!   the forcing term that makes that analytic field an *exact* solution of the
//!   forced equation — turning any solver into a verifiable one.
//! * [`l2_norm`] / [`linf_norm`] / [`rms_norm`] — discrete error norms between a
//!   numeric and an exact field, length-checked (fail-loud on a mismatch).
//!
//! ## Honesty / scope caveats
//!
//! * **Observed order needs asymptotic data.** The fitted slope only equals the
//!   scheme's formal order if every grid in the sweep is in the *asymptotic
//!   range* — fine enough that the leading truncation term dominates. Too-coarse
//!   grids, round-off-dominated too-fine grids, or a non-smooth solution all
//!   corrupt the estimate; this code reports the slope it sees, it cannot tell
//!   you the data were non-asymptotic.
//! * **GCI assumes** three grids in the asymptotic range with a *constant*
//!   refinement ratio `r > 1`, monotone convergence (the three values not
//!   bracketing an oscillation), and a known/observed order `p`. Roache's
//!   recommended safety factor is `1.25` for a careful three-grid study (`3.0`
//!   for a two-grid estimate). Outside those assumptions the band is not
//!   trustworthy.
//! * MMS verifies the **discretisation**, not the physics: it confirms the code
//!   solves its equations correctly, not that those equations model reality
//!   (that is *validation*, which needs experimental data).

use crate::error::UqError;

/// A grid-refinement convergence study: pairs of (grid spacing `h`, error).
///
/// The error is whatever scalar error measure the caller computed at each `h`
/// (e.g. an [`l2_norm`] against an exact or manufactured solution). Build it
/// with [`ConvergenceStudy::new`], then query [`ConvergenceStudy::observed_order`]
/// for the fitted order of accuracy or [`ConvergenceStudy::gci`] for a Grid
/// Convergence Index.
#[derive(Debug, Clone, PartialEq)]
pub struct ConvergenceStudy {
    /// `(h, error)` pairs, in caller order. Errors are stored as given; the
    /// order fit uses their absolute values.
    points: Vec<(f64, f64)>,
}

impl ConvergenceStudy {
    /// Build a study from `(h, error)` pairs.
    ///
    /// # Errors
    /// * [`UqError::EmptyInput`] if fewer than two points are supplied (a slope
    ///   needs at least two).
    /// * [`UqError::OutOfRange`] if any `h <= 0` or any value is non-finite, or
    ///   if any error is exactly zero (its `log` is undefined — a zero-error
    ///   grid carries no order information).
    pub fn new(points: &[(f64, f64)]) -> Result<Self, UqError> {
        if points.len() < 2 {
            return Err(UqError::EmptyInput(format!(
                "convergence study needs >= 2 (h, error) points (got {})",
                points.len()
            )));
        }
        for &(h, e) in points {
            if !h.is_finite() || !e.is_finite() {
                return Err(UqError::OutOfRange(
                    "convergence (h, error) values must be finite".into(),
                ));
            }
            if h <= 0.0 {
                return Err(UqError::OutOfRange(format!(
                    "grid spacing h must be > 0 (got {h})"
                )));
            }
            if e == 0.0 {
                return Err(UqError::OutOfRange(
                    "error must be non-zero (log of zero error is undefined)".into(),
                ));
            }
        }
        Ok(Self {
            points: points.to_vec(),
        })
    }

    /// The **observed order of accuracy** `p`: the least-squares slope of
    /// `log|error|` against `log h`.
    ///
    /// For a scheme of formal order `p` the error scales as `error ≈ C·hᵖ`, so
    /// `log|error| = log C + p·log h` — a straight line whose slope is `p`. The
    /// fit uses every point; the returned slope is positive when the error
    /// shrinks with `h` (the normal case).
    ///
    /// # Errors
    /// [`UqError::LinearAlgebra`] if every `h` is identical (zero variance in
    /// the abscissa makes the slope undefined).
    pub fn observed_order(&self) -> Result<f64, UqError> {
        let xs: Vec<f64> = self.points.iter().map(|&(h, _)| h.ln()).collect();
        let ys: Vec<f64> = self.points.iter().map(|&(_, e)| e.abs().ln()).collect();
        least_squares_slope(&xs, &ys)
    }

    /// Roache's **Grid Convergence Index** between the two finest grids of a
    /// three-grid study, using the supplied `safety_factor` `F_s`.
    ///
    /// Requires exactly three points; they are sorted by `h` so `f1` is the
    /// finest. With a constant refinement ratio `r = h2/h1 = h3/h2` and observed
    /// order
    /// ```text
    ///   p = ln(|f3 − f2| / |f2 − f1|) / ln r
    /// ```
    /// the index is
    /// ```text
    ///   GCI_21 = F_s · |ε21| / (rᵖ − 1),   ε21 = (f2 − f1) / f1.
    /// ```
    /// `F_s = 1.25` is Roache's recommendation for a careful three-grid study.
    /// The result is a relative error band (fraction; multiply by 100 for a
    /// percentage).
    ///
    /// # Errors
    /// * [`UqError::DimensionMismatch`] if the study does not have exactly three
    ///   points.
    /// * [`UqError::OutOfRange`] if `safety_factor <= 0`, the grids are not
    ///   (within a small tolerance) a constant refinement ratio `> 1`, the
    ///   finest value `f1` is zero, or the differences are degenerate (equal
    ///   grids, or `(f3−f2)`/`(f2−f1)` non-positive ⇒ non-monotone, so the
    ///   order is undefined).
    pub fn gci(&self, safety_factor: f64) -> Result<f64, UqError> {
        if self.points.len() != 3 {
            return Err(UqError::DimensionMismatch(format!(
                "GCI needs exactly 3 grids (got {})",
                self.points.len()
            )));
        }
        if !safety_factor.is_finite() || safety_factor <= 0.0 {
            return Err(UqError::OutOfRange(format!(
                "GCI safety factor must be a finite value > 0 (got {safety_factor})"
            )));
        }

        // Sort by h ascending: index 0 = finest grid.
        let mut pts = self.points.clone();
        pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let (h1, f1) = pts[0];
        let (h2, f2) = pts[1];
        let (h3, f3) = pts[2];

        let r21 = h2 / h1;
        let r32 = h3 / h2;
        if r21 <= 1.0 || (r21 - r32).abs() > 1e-6 * r21 {
            return Err(UqError::OutOfRange(format!(
                "GCI requires a constant refinement ratio r > 1 \
                 (got r21={r21}, r32={r32})"
            )));
        }
        if f1 == 0.0 {
            return Err(UqError::OutOfRange(
                "GCI: finest-grid value f1 is zero (relative error undefined)".into(),
            ));
        }

        let e21 = f2 - f1;
        let e32 = f3 - f2;
        // Monotone-convergence requirement: the ratio must be positive so its
        // log is real.
        let ratio = e32 / e21;
        if !ratio.is_finite() || ratio <= 0.0 {
            return Err(UqError::OutOfRange(
                "GCI: non-monotone or degenerate grid differences (observed order undefined)"
                    .into(),
            ));
        }

        let p = ratio.ln() / r21.ln();
        let denom = r21.powf(p) - 1.0;
        if denom.abs() <= f64::EPSILON {
            return Err(UqError::OutOfRange(
                "GCI: rᵖ − 1 ≈ 0 (degenerate order/ratio)".into(),
            ));
        }
        let rel_err = (e21 / f1).abs();
        Ok(safety_factor * rel_err / denom)
    }

    /// The observed order implied by the three-grid ratio (the `p` used inside
    /// [`ConvergenceStudy::gci`]). Convenience accessor for reporting alongside
    /// the GCI.
    ///
    /// # Errors
    /// Same conditions as [`ConvergenceStudy::gci`] regarding the three-grid
    /// requirement and a well-defined, monotone ratio.
    pub fn observed_order_three_grid(&self) -> Result<f64, UqError> {
        if self.points.len() != 3 {
            return Err(UqError::DimensionMismatch(format!(
                "three-grid order needs exactly 3 grids (got {})",
                self.points.len()
            )));
        }
        let mut pts = self.points.clone();
        pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let (h1, f1) = pts[0];
        let (h2, f2) = pts[1];
        let (h3, f3) = pts[2];
        let r = h2 / h1;
        let r32 = h3 / h2;
        if r <= 1.0 || (r - r32).abs() > 1e-6 * r {
            return Err(UqError::OutOfRange(
                "three-grid order requires a constant refinement ratio r > 1".into(),
            ));
        }
        let ratio = (f3 - f2) / (f2 - f1);
        if !ratio.is_finite() || ratio <= 0.0 {
            return Err(UqError::OutOfRange(
                "three-grid order: non-monotone / degenerate differences".into(),
            ));
        }
        Ok(ratio.ln() / r.ln())
    }
}

/// Least-squares slope of `y` on `x` (ordinary linear regression slope).
fn least_squares_slope(x: &[f64], y: &[f64]) -> Result<f64, UqError> {
    let n = x.len() as f64;
    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for (&xi, &yi) in x.iter().zip(y.iter()) {
        let dx = xi - mean_x;
        sxx += dx * dx;
        sxy += dx * (yi - mean_y);
    }
    if sxx <= f64::EPSILON {
        return Err(UqError::LinearAlgebra(
            "slope undefined: all abscissa values are identical".into(),
        ));
    }
    Ok(sxy / sxx)
}

/// The discrete **L² (root-mean-square-times-√n) error norm**
/// `√(Σ (a_i − b_i)²)` between a numeric field `numeric` and an exact field
/// `exact`.
///
/// This is the Euclidean norm of the error vector. (For a *grid-independent*
/// measure prefer [`rms_norm`], which divides out the sample count.)
///
/// # Errors
/// * [`UqError::EmptyInput`] if the fields are empty.
/// * [`UqError::DimensionMismatch`] if the two fields differ in length.
pub fn l2_norm(numeric: &[f64], exact: &[f64]) -> Result<f64, UqError> {
    check_pair(numeric, exact, "l2_norm")?;
    let ss: f64 = numeric
        .iter()
        .zip(exact.iter())
        .map(|(&a, &b)| (a - b) * (a - b))
        .sum();
    Ok(ss.sqrt())
}

/// The discrete **L∞ (maximum) error norm** `max_i |a_i − b_i|`.
///
/// # Errors
/// * [`UqError::EmptyInput`] if the fields are empty.
/// * [`UqError::DimensionMismatch`] if the two fields differ in length.
pub fn linf_norm(numeric: &[f64], exact: &[f64]) -> Result<f64, UqError> {
    check_pair(numeric, exact, "linf_norm")?;
    Ok(numeric
        .iter()
        .zip(exact.iter())
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0_f64, f64::max))
}

/// The **root-mean-square error norm** `√(Σ (a_i − b_i)² / n)` — the L² norm
/// normalised by the field length, so it does not grow merely because the mesh
/// is refined.
///
/// # Errors
/// * [`UqError::EmptyInput`] if the fields are empty.
/// * [`UqError::DimensionMismatch`] if the two fields differ in length.
pub fn rms_norm(numeric: &[f64], exact: &[f64]) -> Result<f64, UqError> {
    check_pair(numeric, exact, "rms_norm")?;
    let n = numeric.len() as f64;
    let ss: f64 = numeric
        .iter()
        .zip(exact.iter())
        .map(|(&a, &b)| (a - b) * (a - b))
        .sum();
    Ok((ss / n).sqrt())
}

/// Shared length/emptiness check for the error norms (fail-loud).
fn check_pair(a: &[f64], b: &[f64], who: &str) -> Result<(), UqError> {
    if a.is_empty() || b.is_empty() {
        return Err(UqError::EmptyInput(format!("{who}: empty field")));
    }
    if a.len() != b.len() {
        return Err(UqError::DimensionMismatch(format!(
            "{who}: field lengths differ ({} vs {})",
            a.len(),
            b.len()
        )));
    }
    Ok(())
}

/// A **Method of Manufactured Solutions** problem on a 1-D point set.
///
/// MMS turns *any* differential operator `L[·]` into a verifiable test: choose
/// an arbitrary smooth analytic field `u(x)`, evaluate `L[u]` analytically (or
/// to high accuracy), and **define** the forcing term `s(x) = L[u](x)`. Then
/// `u` is, by construction, the *exact* solution of the forced equation
/// `L[U] = s`. A correct solver fed `s` must reproduce `u` to within
/// discretisation error, so the error → 0 at the scheme's order under
/// refinement.
///
/// This is a tiny framework: the caller supplies the analytic solution `u(x)`
/// and the operator-applied-to-the-solution `L[u](x)` as closures, both
/// sampled on a set of node coordinates. [`Mms::forcing`] returns `s` at the
/// nodes; [`Mms::residual`] checks that a candidate discrete operator output
/// reproduces that forcing (≈ 0 residual) — i.e. that the manufactured field is
/// reproduced exactly.
pub struct Mms<U, L>
where
    U: Fn(f64) -> f64,
    L: Fn(f64) -> f64,
{
    solution: U,
    operator_on_solution: L,
}

impl<U, L> Mms<U, L>
where
    U: Fn(f64) -> f64,
    L: Fn(f64) -> f64,
{
    /// Create an MMS problem from the analytic solution `u(x)` and the analytic
    /// value of the differential operator applied to it, `L[u](x)`.
    pub fn new(solution: U, operator_on_solution: L) -> Self {
        Self {
            solution,
            operator_on_solution,
        }
    }

    /// The exact manufactured field `u(x)` sampled at `nodes`.
    #[must_use]
    pub fn exact_solution(&self, nodes: &[f64]) -> Vec<f64> {
        nodes.iter().map(|&x| (self.solution)(x)).collect()
    }

    /// The forcing term `s(x) = L[u](x)` sampled at `nodes` — the source that
    /// makes `u` an exact solution of `L[U] = s`.
    #[must_use]
    pub fn forcing(&self, nodes: &[f64]) -> Vec<f64> {
        nodes
            .iter()
            .map(|&x| (self.operator_on_solution)(x))
            .collect()
    }

    /// Residual `‖L_applied − s‖∞` between a candidate discrete-operator output
    /// `operator_output` (the result of applying the *discretised* operator to
    /// the manufactured field) and the manufactured forcing `s` at `nodes`.
    ///
    /// For an exact (or exactly-reproducing) operator this is `≈ 0`, the defining
    /// MMS verification check.
    ///
    /// # Errors
    /// [`UqError::DimensionMismatch`] if `operator_output.len() != nodes.len()`,
    /// or [`UqError::EmptyInput`] if `nodes` is empty.
    pub fn residual(&self, nodes: &[f64], operator_output: &[f64]) -> Result<f64, UqError> {
        let s = self.forcing(nodes);
        linf_norm(operator_output, &s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Synthetic second-order data: error = C·h². The fitted observed order
    /// must be ≈ 2.0.
    #[test]
    fn observed_order_second_order() {
        let c = 0.37;
        let points: Vec<(f64, f64)> = [1.0, 0.5, 0.25, 0.125, 0.0625]
            .iter()
            .map(|&h| (h, c * h * h))
            .collect();
        let study = ConvergenceStudy::new(&points).unwrap();
        let p = study.observed_order().unwrap();
        assert!((p - 2.0).abs() < 1e-3, "observed order = {p}");
    }

    /// Synthetic first-order data: error = C·h¹. Observed order ≈ 1.0.
    #[test]
    fn observed_order_first_order() {
        let c = 1.9;
        let points: Vec<(f64, f64)> = [1.0, 0.5, 0.25, 0.125, 0.0625]
            .iter()
            .map(|&h| (h, c * h))
            .collect();
        let study = ConvergenceStudy::new(&points).unwrap();
        let p = study.observed_order().unwrap();
        assert!((p - 1.0).abs() < 1e-3, "observed order = {p}");
    }

    /// Roache's textbook three-grid GCI worked example
    /// (Roache 1994/1998; reproduced in the NASA/NPARC and AIAA verification
    /// guides). Grids refine by r = 2 with values
    ///   f1 = 0.97050 (fine), f2 = 0.96854 (medium), f3 = 0.96178 (coarse).
    /// The published results are observed order p ≈ 1.786 and
    /// GCI_fine (Fs = 1.25) ≈ 0.00103 (0.103 %).
    #[test]
    fn gci_reproduces_roache_textbook_example() {
        // h's just need ratio 2; finest first is fine (smallest h).
        let points = [(1.0, 0.97050), (2.0, 0.96854), (4.0, 0.96178)];
        let study = ConvergenceStudy::new(&points).unwrap();

        let p = study.observed_order_three_grid().unwrap();
        assert!((p - 1.786_30).abs() < 1e-2, "observed order = {p}");

        let gci = study.gci(1.25).unwrap();
        assert!(
            (gci - 0.001_03).abs() < 5e-5,
            "GCI = {gci} (expected ≈ 0.00103)"
        );
    }

    #[test]
    fn error_norms_match_hand_computed() {
        // numeric − exact = [3, 4] → L2 = 5, L∞ = 4, RMS = 5/√2.
        let numeric = [4.0, 6.0];
        let exact = [1.0, 2.0];
        let l2 = l2_norm(&numeric, &exact).unwrap();
        let linf = linf_norm(&numeric, &exact).unwrap();
        let rms = rms_norm(&numeric, &exact).unwrap();
        assert!((l2 - 5.0).abs() < 1e-12, "L2 = {l2}");
        assert!((linf - 4.0).abs() < 1e-12, "Linf = {linf}");
        assert!(
            (rms - 5.0 / std::f64::consts::SQRT_2).abs() < 1e-12,
            "RMS = {rms}"
        );
    }

    #[test]
    fn error_norms_fail_loud_on_mismatch() {
        assert!(l2_norm(&[1.0, 2.0], &[1.0]).is_err());
        assert!(linf_norm(&[], &[]).is_err());
        assert!(rms_norm(&[1.0], &[1.0, 2.0]).is_err());
    }

    /// MMS worked example for the 1-D Poisson operator `L[u] = −u''`.
    ///
    /// Manufacture the solution u(x) = sin(πx). Then −u'' = π²·sin(πx), so the
    /// forcing is s(x) = π²·sin(πx). Feeding that exact operator output back in
    /// must give a residual of (essentially) zero — the manufactured field is
    /// reproduced exactly.
    #[test]
    fn mms_forcing_makes_solution_exact() {
        use std::f64::consts::PI;
        let mms = Mms::new(
            |x: f64| (PI * x).sin(),           // u(x)
            |x: f64| PI * PI * (PI * x).sin(), // L[u] = −u''
        );
        let nodes: Vec<f64> = (0..=10).map(|i| i as f64 / 10.0).collect();

        // The "discrete operator output" here is the exact operator applied to
        // the exact field — for a perfect operator it equals the forcing, so the
        // residual must vanish.
        let exact = mms.exact_solution(&nodes);
        let operator_output: Vec<f64> = exact
            .iter()
            .zip(nodes.iter())
            .map(|(_, &x)| PI * PI * (PI * x).sin())
            .collect();

        let residual = mms.residual(&nodes, &operator_output).unwrap();
        assert!(residual < 1e-12, "MMS residual = {residual}");

        // And the manufactured forcing equals s(x) = π² sin(πx) node-wise.
        let s = mms.forcing(&nodes);
        for (&x, &si) in nodes.iter().zip(s.iter()) {
            assert!((si - PI * PI * (PI * x).sin()).abs() < 1e-12);
        }
    }

    /// A genuine second-order central-difference discretisation of −u'' fed the
    /// MMS forcing recovers the manufactured field to second order: the L² error
    /// drops by ≈ 4× when h is halved.
    #[test]
    fn mms_central_difference_converges_second_order() {
        use std::f64::consts::PI;
        // Solve −u'' = s on (0,1) with Dirichlet BCs from u(x)=sin(πx), via the
        // standard tridiagonal central-difference scheme, on two grids.
        let solve_error = |n: usize| -> f64 {
            let h = 1.0 / n as f64;
            // Interior unknowns i = 1..n (nodes x_i = i h), u_0 = u_n = 0.
            let m = n - 1;
            // RHS s_i = π² sin(π x_i); the tridiagonal system is
            // (-u_{i-1} + 2u_i - u_{i+1})/h² = s_i.
            let mut a = vec![0.0; m]; // sub-diagonal
            let mut b = vec![0.0; m]; // diagonal
            let mut c = vec![0.0; m]; // super-diagonal
            let mut d = vec![0.0; m]; // rhs
            for i in 0..m {
                let x = (i + 1) as f64 * h;
                a[i] = -1.0 / (h * h);
                b[i] = 2.0 / (h * h);
                c[i] = -1.0 / (h * h);
                d[i] = PI * PI * (PI * x).sin();
            }
            // Thomas algorithm.
            for i in 1..m {
                let w = a[i] / b[i - 1];
                b[i] -= w * c[i - 1];
                d[i] -= w * d[i - 1];
            }
            let mut u = vec![0.0; m];
            u[m - 1] = d[m - 1] / b[m - 1];
            for i in (0..m - 1).rev() {
                u[i] = (d[i] - c[i] * u[i + 1]) / b[i];
            }
            // Compare to exact at interior nodes.
            let exact: Vec<f64> = (0..m).map(|i| (PI * (i + 1) as f64 * h).sin()).collect();
            rms_norm(&u, &exact).unwrap()
        };

        let e_coarse = solve_error(20);
        let e_fine = solve_error(40);
        let rate = e_coarse / e_fine;
        // Halving h should cut the error ~4× for a 2nd-order scheme.
        assert!(rate > 3.5 && rate < 4.5, "refinement ratio = {rate}");
    }

    #[test]
    fn convergence_study_rejects_bad_input() {
        assert!(ConvergenceStudy::new(&[(1.0, 0.1)]).is_err()); // < 2 points
        assert!(ConvergenceStudy::new(&[(1.0, 0.1), (-0.5, 0.2)]).is_err()); // h <= 0
        assert!(ConvergenceStudy::new(&[(1.0, 0.0), (0.5, 0.1)]).is_err()); // zero error
    }

    #[test]
    fn gci_rejects_non_three_grid_and_bad_ratio() {
        let two = ConvergenceStudy::new(&[(1.0, 0.1), (0.5, 0.05)]).unwrap();
        assert!(two.gci(1.25).is_err()); // not 3 grids
                                         // Non-constant refinement ratio.
        let bad = ConvergenceStudy::new(&[(1.0, 0.97), (2.0, 0.96), (5.0, 0.95)]).unwrap();
        assert!(bad.gci(1.25).is_err());
    }
}
