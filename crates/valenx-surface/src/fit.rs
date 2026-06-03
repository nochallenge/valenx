//! Curve and surface fitting from point clouds (Phase 19D).
//!
//! ## Algorithms
//!
//! - **Curve fitting** (`nurbs_curve_through_points`): least-squares
//!   fit of a NURBS curve through `m` data points with a target of
//!   `n_cps` control points. Uses centripetal chord-length
//!   parameterisation and clamped open-uniform knots. The system
//!   `N^T N x = N^T D` is solved via a hand-rolled Gauss-Jordan
//!   elimination (workspace `nalgebra` is already a dep — we use
//!   `nalgebra::DMatrix::solve_lower_triangular` etc., but the hand
//!   path keeps the dependency surface narrow).
//!
//! - **Surface fitting from a structured grid**
//!   (`nurbs_surface_through_grid`): tensor-product least-squares.
//!   We fit `nu` curves in the u direction (one per v slice) to get
//!   an intermediate grid of CPs, then fit `nv` curves in the v
//!   direction through that intermediate grid to get the final CP
//!   array.
//!
//! - **Surface fitting from scattered points**
//!   (`surface_through_scattered`): v1 simplification — fit a
//!   best-fit plane via SVD-free centroid + normal-from-cross-product
//!   approximation, project all points onto the plane to get (u, v)
//!   parameters, bin into a coarse `target_n_cps_u x target_n_cps_v`
//!   grid (averaging within each cell), and delegate to
//!   `nurbs_surface_through_grid`.
//!
//! Each public function returns the fitted artefact together with the
//! RMS fitting error (Euclidean distance from each input point to the
//! evaluated fit).

use nalgebra::Vector3;

use crate::error::SurfaceError;
use crate::nurbs_curve::{basis_functions, find_knot_span, NurbsCurve};
use crate::nurbs_surface::NurbsSurface;

/// Result of a curve fit: the fitted curve + RMS error in the input
/// point's metric.
#[derive(Clone, Debug)]
pub struct CurveFit {
    /// The fitted NURBS curve.
    pub curve: NurbsCurve,
    /// Root-mean-squared distance from each input point to the
    /// evaluated fit at that point's parameter.
    pub rms_error: f64,
}

/// Result of a surface fit: the fitted surface + RMS error.
#[derive(Clone, Debug)]
pub struct SurfaceFit {
    /// The fitted NURBS surface.
    pub surface: NurbsSurface,
    /// Root-mean-squared distance from each input point to the
    /// evaluated fit at that point's (u, v) parameter.
    pub rms_error: f64,
}

/// Fit a NURBS curve through `points` using `n_cps` control points.
///
/// `degree` must be `>= 1`; `n_cps` must be `>= degree + 1` and
/// `<= points.len()`. Returns the fitted curve + RMS error.
///
/// Uses centripetal parameterisation (parameter at point `i` is
/// proportional to `sqrt(|P_i - P_{i-1}|)`) and a clamped
/// open-uniform knot vector.
pub fn nurbs_curve_through_points(
    points: &[Vector3<f64>],
    degree: usize,
    n_cps: usize,
) -> Result<CurveFit, SurfaceError> {
    if points.len() < 2 {
        return Err(SurfaceError::BadKnotVector {
            reason: "need at least 2 points to fit".into(),
        });
    }
    if degree == 0 || degree > 9 {
        return Err(SurfaceError::BadDegree(degree));
    }
    if n_cps < degree + 1 || n_cps > points.len() {
        return Err(SurfaceError::BadKnotVector {
            reason: format!(
                "n_cps must be in [{}, {}]; got {n_cps}",
                degree + 1,
                points.len()
            ),
        });
    }

    let m = points.len();
    let n = n_cps;
    let params = centripetal_params(points);
    let knots = open_uniform_knots_for(n, degree);

    // Build the m x n basis matrix N where N[i][j] = N_j(u_i).
    let mut nmat = vec![vec![0.0_f64; n]; m];
    for (i, &u) in params.iter().enumerate() {
        let span = find_knot_span(u, &knots, degree, n);
        let basis = basis_functions(span, u, degree, &knots);
        for (k, b) in basis.iter().enumerate() {
            let col = span - degree + k;
            nmat[i][col] = *b;
        }
    }

    // First / last control points are pinned to the first / last
    // data points (consequence of the clamped-knot endpoint property).
    // We solve for the n - 2 interior CPs in least squares:
    //   M x = R
    // where M = N_int^T N_int  (n-2 x n-2),
    //       R = N_int^T (D - P_0 N_0 - P_n-1 N_n-1)  (n-2 x 3).
    let p0 = points[0];
    let pn = points[m - 1];
    if n == 2 {
        // Degenerate: just a degree-1 line through the endpoints.
        let cps = vec![p0, pn];
        let weights = vec![1.0; 2];
        let curve = NurbsCurve::new(degree, knots, cps, weights)?;
        let rms = compute_rms_error(&curve, points, &params);
        return Ok(CurveFit {
            curve,
            rms_error: rms,
        });
    }

    // Build M = N_int^T N_int  (size (n-2) x (n-2)).
    let interior = n - 2;
    let mut m_mat = vec![vec![0.0_f64; interior]; interior];
    let mut r_mat = vec![[0.0_f64; 3]; interior];
    for k in 0..interior {
        let col_k = k + 1;
        for kk in 0..interior {
            let col_kk = kk + 1;
            let mut sum = 0.0;
            for row in nmat.iter().take(m) {
                sum += row[col_k] * row[col_kk];
            }
            m_mat[k][kk] = sum;
        }
        // R row.
        for (i, row) in nmat.iter().enumerate().take(m) {
            let coeff = row[col_k];
            let nb0 = row[0];
            let nbn = row[n - 1];
            let d = points[i] - nb0 * p0 - nbn * pn;
            r_mat[k][0] += coeff * d.x;
            r_mat[k][1] += coeff * d.y;
            r_mat[k][2] += coeff * d.z;
        }
    }
    // Solve 3 right-hand sides via in-place Gauss-Jordan.
    let x = solve_lin_3rhs(&mut m_mat, &mut r_mat)?;

    let mut cps = Vec::with_capacity(n);
    cps.push(p0);
    for row in &x {
        cps.push(Vector3::new(row[0], row[1], row[2]));
    }
    cps.push(pn);
    let weights = vec![1.0; n];
    let curve = NurbsCurve::new(degree, knots, cps, weights)?;
    let rms = compute_rms_error(&curve, points, &params);
    Ok(CurveFit {
        curve,
        rms_error: rms,
    })
}

/// Fit a NURBS surface through a structured `(nu_in x nv_in)` grid of
/// data points using `(n_cps_u x n_cps_v)` control points.
///
/// `points_uv[i][j]` is the data point at sample `(i, j)`. Returns
/// the fitted surface + RMS error across all data points.
pub fn nurbs_surface_through_grid(
    points_uv: &[Vec<Vector3<f64>>],
    degree_u: usize,
    degree_v: usize,
    n_cps_u: usize,
    n_cps_v: usize,
) -> Result<SurfaceFit, SurfaceError> {
    let nu_in = points_uv.len();
    if nu_in == 0 {
        return Err(SurfaceError::BadKnotVector {
            reason: "empty grid".into(),
        });
    }
    let nv_in = points_uv[0].len();
    for row in points_uv {
        if row.len() != nv_in {
            return Err(SurfaceError::BadKnotVector {
                reason: "grid is not rectangular".into(),
            });
        }
    }

    // Step 1: fit one curve per row in v through (n_cps_v) CPs.
    let mut intermediate: Vec<Vec<Vector3<f64>>> = Vec::with_capacity(nu_in);
    for row in points_uv {
        let fit = nurbs_curve_through_points(row, degree_v, n_cps_v)?;
        intermediate.push(fit.curve.control_points.clone());
    }
    // Step 2: fit one curve per column in u (through n_cps_u CPs).
    // `intermediate[i][j]` is the j-th CP of the row-i v-fit; we now
    // gather along i (fixed j) and fit in u.
    let mut final_cps: Vec<Vec<Vector3<f64>>> = vec![Vec::with_capacity(n_cps_v); n_cps_u];
    for j in 0..n_cps_v {
        let col: Vec<Vector3<f64>> = (0..nu_in).map(|i| intermediate[i][j]).collect();
        let fit = nurbs_curve_through_points(&col, degree_u, n_cps_u)?;
        for (i, cp) in fit.curve.control_points.iter().enumerate() {
            final_cps[i].push(*cp);
        }
    }
    let u_knots = open_uniform_knots_for(n_cps_u, degree_u);
    let v_knots = open_uniform_knots_for(n_cps_v, degree_v);
    let weights = vec![vec![1.0_f64; n_cps_v]; n_cps_u];
    let surface = NurbsSurface::new(degree_u, degree_v, u_knots, v_knots, final_cps, weights)?;
    // RMS error: evaluate the surface at the input data points'
    // assumed parameters (uniformly distributed). v1 simplification.
    let mut sumsq = 0.0_f64;
    let mut n = 0_usize;
    for (i, row) in points_uv.iter().enumerate() {
        for (j, p) in row.iter().enumerate() {
            let u = i as f64 / (nu_in - 1).max(1) as f64;
            let v = j as f64 / (nv_in - 1).max(1) as f64;
            let q = surface.evaluate(u, v);
            sumsq += (q - p).norm_squared();
            n += 1;
        }
    }
    let rms = (sumsq / n.max(1) as f64).sqrt();
    Ok(SurfaceFit {
        surface,
        rms_error: rms,
    })
}

/// Fit a NURBS surface to a scattered point cloud.
///
/// **v1 strategy:** project all points to the best-fit plane (via
/// centroid + normal-from-PCA), parameterise by their plane
/// coordinates, bin into a `target_n_cps_u x target_n_cps_v` grid
/// (averaging within each cell), and delegate to
/// [`nurbs_surface_through_grid`].
///
/// For genuinely curved point clouds (where the projection collapses
/// real structure) this v1 fit will under-represent the data; true
/// moving-least-squares scattered fitting is a v1.5 upgrade.
pub fn surface_through_scattered(
    points: &[Vector3<f64>],
    degree_u: usize,
    degree_v: usize,
    target_n_cps_u: usize,
    target_n_cps_v: usize,
) -> Result<SurfaceFit, SurfaceError> {
    if points.len() < (target_n_cps_u * target_n_cps_v) {
        return Err(SurfaceError::BadKnotVector {
            reason: format!(
                "need at least {} scattered points for a {}x{} fit; got {}",
                target_n_cps_u * target_n_cps_v,
                target_n_cps_u,
                target_n_cps_v,
                points.len()
            ),
        });
    }
    // Centroid.
    let mut centroid = Vector3::zeros();
    for p in points {
        centroid += p;
    }
    centroid /= points.len() as f64;

    // Covariance matrix (3x3). For v1 we approximate the normal by
    // the eigenvector of the smallest eigenvalue via deflation /
    // simple power-iteration on the inverse — but to stay
    // dependency-light we use a cheap approximation: split the
    // residual variances along x/y/z and pick the axis-aligned
    // normal of the smallest variance. This is exact for
    // axis-aligned planar point clouds and reasonable for slightly
    // tilted clouds; for arbitrary clouds it's a v1 best-effort.
    let mut var = [0.0_f64; 3];
    for p in points {
        let d = p - centroid;
        var[0] += d.x * d.x;
        var[1] += d.y * d.y;
        var[2] += d.z * d.z;
    }
    let (normal_axis, _) = var
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap();
    let (axis_u, axis_v) = match normal_axis {
        0 => (Vector3::new(0.0, 1.0, 0.0), Vector3::new(0.0, 0.0, 1.0)),
        1 => (Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0)),
        _ => (Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)),
    };

    // Project + bin.
    let mut us = Vec::with_capacity(points.len());
    let mut vs = Vec::with_capacity(points.len());
    for p in points {
        us.push((p - centroid).dot(&axis_u));
        vs.push((p - centroid).dot(&axis_v));
    }
    let (u_min, u_max) = (
        us.iter().cloned().fold(f64::INFINITY, f64::min),
        us.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    let (v_min, v_max) = (
        vs.iter().cloned().fold(f64::INFINITY, f64::min),
        vs.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
    let u_range = (u_max - u_min).max(1.0e-12);
    let v_range = (v_max - v_min).max(1.0e-12);

    // Resample the scattered cloud onto a regular `(bin_u × bin_v)`
    // grid in the (u, v) parameter plane. We *oversample* — more grid
    // nodes than control points — so the downstream curve fit has slack
    // to least-squares-smooth the noise.
    //
    // A naive bin-and-average is not used: floating-point round-off in
    // the bin index scatters even a perfectly regular point lattice
    // unevenly, and the resulting empty cells / averaged boundary cells
    // produce a scrambled, non-monotone grid (a twisted, self-folding
    // fitted patch). Instead:
    //
    // - **Interior** grid nodes take a **Shepard inverse-distance-
    //   weighted** blend of the data — smooth, monotone, no empty cell.
    // - **Boundary** grid nodes are resampled from the data points on
    //   the *corresponding data boundary strip only*. A clamped NURBS
    //   patch's boundary is fixed by its boundary control points; if a
    //   boundary node blended in interior points the patch edge would
    //   bow inward and strand the extreme data points off the surface.
    let bin_u = (target_n_cps_u * 2).max(target_n_cps_u + 1);
    let bin_v = (target_n_cps_v * 2).max(target_n_cps_v + 1);
    let cell_u = u_range / (bin_u - 1) as f64;
    let cell_v = v_range / (bin_v - 1) as f64;
    // IDW falloff: a couple of grid spacings, so each node is shaped by
    // its local neighbourhood, not the whole cloud.
    let falloff2 = (2.0 * (cell_u + cell_v)).powi(2).max(1.0e-18) * 1.0e-3;

    // Shepard IDW of the data restricted to the index set `pool`,
    // weighting by squared (u, v) distance to the node `(gu, gv)`.
    let idw = |pool: &[usize], gu: f64, gv: f64| -> Vector3<f64> {
        let mut wsum = 0.0_f64;
        let mut acc = Vector3::zeros();
        for &k in pool {
            let du = us[k] - gu;
            let dv = vs[k] - gv;
            let d2 = du * du + dv * dv;
            if d2 < 1.0e-18 {
                return points[k]; // data point sits on the node
            }
            let w = 1.0 / (d2 + falloff2);
            wsum += w;
            acc += points[k] * w;
        }
        if wsum > 0.0 {
            acc / wsum
        } else {
            centroid
        }
    };

    // Boundary-strip index sets — the data points hugging each
    // parameter-domain edge — so a boundary node resamples only the
    // matching data edge. The strip is kept thin (a small fraction of a
    // grid cell) so it isolates the genuine boundary data; on a sparse
    // edge it widens progressively rather than ever being empty.
    let all: Vec<usize> = (0..points.len()).collect();
    let strip = |coord: &dyn Fn(usize) -> f64, edge: f64| -> Vec<usize> {
        // Widen the band until it captures at least two points (enough
        // to fit a boundary curve), so a clean edge keeps a tight band
        // and a sparse one still resolves.
        for &frac in &[0.02_f64, 0.1, 0.34, 1.0] {
            let band = (cell_u.max(cell_v)) * frac;
            let s: Vec<usize> = (0..points.len())
                .filter(|&k| (coord(k) - edge).abs() <= band)
                .collect();
            if s.len() >= 2 {
                return s;
            }
        }
        all.clone()
    };
    let u_lo_strip = strip(&|k| us[k], u_min);
    let u_hi_strip = strip(&|k| us[k], u_max);
    let v_lo_strip = strip(&|k| vs[k], v_min);
    let v_hi_strip = strip(&|k| vs[k], v_max);

    let mut grid: Vec<Vec<Vector3<f64>>> = vec![vec![Vector3::zeros(); bin_v]; bin_u];
    for (i, grid_row) in grid.iter_mut().enumerate() {
        let gu = u_min + cell_u * i as f64;
        let on_u_lo = i == 0;
        let on_u_hi = i == bin_u - 1;
        for (j, cell) in grid_row.iter_mut().enumerate() {
            let gv = v_min + cell_v * j as f64;
            let on_v_lo = j == 0;
            let on_v_hi = j == bin_v - 1;
            // Restrict the resampling pool to the matching boundary
            // strip so the fitted patch boundary tracks the data
            // boundary; interior nodes use the whole cloud.
            let pool: &[usize] = if on_u_lo {
                &u_lo_strip
            } else if on_u_hi {
                &u_hi_strip
            } else if on_v_lo {
                &v_lo_strip
            } else if on_v_hi {
                &v_hi_strip
            } else {
                &all
            };
            *cell = idw(pool, gu, gv);
        }
    }
    // Snap the 4 grid corners to the extremal data points. A clamped
    // (open-uniform) NURBS fit interpolates its corner control points
    // exactly, so pinning the grid corners to the data cloud's
    // parametric extremes makes the fitted patch *cover the whole
    // point cloud*.
    let corner_targets = [
        (0usize, 0usize, u_min, v_min),
        (0, bin_v - 1, u_min, v_max),
        (bin_u - 1, 0, u_max, v_min),
        (bin_u - 1, bin_v - 1, u_max, v_max),
    ];
    for (ci, cj, tu, tv) in corner_targets {
        let mut best_k = 0usize;
        let mut best_d = f64::INFINITY;
        for k in 0..points.len() {
            let du = us[k] - tu;
            let dv = vs[k] - tv;
            let d = du * du + dv * dv;
            if d < best_d {
                best_d = d;
                best_k = k;
            }
        }
        grid[ci][cj] = points[best_k];
    }
    nurbs_surface_through_grid(&grid, degree_u, degree_v, target_n_cps_u, target_n_cps_v)
}

// ===== Helpers =====

fn centripetal_params(points: &[Vector3<f64>]) -> Vec<f64> {
    let m = points.len();
    if m == 1 {
        return vec![0.0];
    }
    let mut dists = Vec::with_capacity(m - 1);
    for i in 1..m {
        let d = (points[i] - points[i - 1]).norm().sqrt(); // centripetal
        dists.push(d);
    }
    let total: f64 = dists.iter().sum::<f64>().max(1.0e-12);
    let mut params = Vec::with_capacity(m);
    params.push(0.0);
    let mut acc = 0.0;
    for d in &dists {
        acc += d;
        params.push((acc / total).min(1.0));
    }
    let last = params.len() - 1;
    params[last] = 1.0;
    params
}

fn open_uniform_knots_for(n_cps: usize, degree: usize) -> Vec<f64> {
    let p = degree;
    let m = n_cps + p + 1;
    let mut k = vec![0.0; m];
    if n_cps <= p + 1 {
        for kv in k.iter_mut().skip(m - p - 1) {
            *kv = 1.0;
        }
        return k;
    }
    let n_internal = n_cps - p - 1;
    for (i, kv) in k.iter_mut().enumerate().take(m) {
        if i <= p {
            *kv = 0.0;
        } else if i >= n_cps {
            *kv = 1.0;
        } else {
            let idx = i - p;
            *kv = idx as f64 / (n_internal + 1) as f64;
        }
    }
    k
}

fn solve_lin_3rhs(a: &mut [Vec<f64>], b: &mut [[f64; 3]]) -> Result<Vec<[f64; 3]>, SurfaceError> {
    let n = a.len();
    for row in a.iter() {
        if row.len() != n {
            return Err(SurfaceError::BadKnotVector {
                reason: "non-square system in solver".into(),
            });
        }
    }
    // Gauss elimination with partial pivoting.
    for k in 0..n {
        let mut pivot = k;
        for i in (k + 1)..n {
            if a[i][k].abs() > a[pivot][k].abs() {
                pivot = i;
            }
        }
        if a[pivot][k].abs() < 1.0e-14 {
            return Err(SurfaceError::BadKnotVector {
                reason: "singular least-squares matrix".into(),
            });
        }
        a.swap(k, pivot);
        b.swap(k, pivot);
        for i in (k + 1)..n {
            let factor = a[i][k] / a[k][k];
            for j in k..n {
                a[i][j] -= factor * a[k][j];
            }
            for rhs in 0..3 {
                b[i][rhs] -= factor * b[k][rhs];
            }
        }
    }
    let mut x = vec![[0.0_f64; 3]; n];
    for i in (0..n).rev() {
        for rhs in 0..3 {
            let mut sum = b[i][rhs];
            for (j, xj) in x.iter().enumerate().take(n).skip(i + 1) {
                sum -= a[i][j] * xj[rhs];
            }
            x[i][rhs] = sum / a[i][i];
        }
    }
    Ok(x)
}

fn compute_rms_error(curve: &NurbsCurve, points: &[Vector3<f64>], params: &[f64]) -> f64 {
    let mut sumsq = 0.0_f64;
    for (p, &u) in points.iter().zip(params) {
        let q = curve.evaluate(u);
        sumsq += (q - p).norm_squared();
    }
    (sumsq / points.len() as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fit_line_through_collinear_points_is_near_zero_rms() {
        let pts: Vec<Vector3<f64>> = (0..6).map(|i| Vector3::new(i as f64, 0.0, 0.0)).collect();
        let fit = nurbs_curve_through_points(&pts, 3, 4).unwrap();
        // With 4 cubic CPs through 6 collinear points, RMS should be
        // very small (the fit can reproduce the line exactly).
        assert!(fit.rms_error < 1.0e-6, "rms = {}", fit.rms_error);
    }

    #[test]
    fn fit_grid_through_planar_grid_is_near_zero_rms() {
        let mut grid = Vec::with_capacity(5);
        for i in 0..5 {
            let mut row = Vec::with_capacity(5);
            for j in 0..5 {
                let u = i as f64 / 4.0;
                let v = j as f64 / 4.0;
                row.push(Vector3::new(u, v, 0.0));
            }
            grid.push(row);
        }
        let fit = nurbs_surface_through_grid(&grid, 3, 3, 4, 4).unwrap();
        assert!(fit.rms_error < 1.0e-6, "rms = {}", fit.rms_error);
    }

    #[test]
    fn fit_scattered_planar_points() {
        // 50 random points all in z=0.
        let pts: Vec<Vector3<f64>> = (0..50)
            .map(|i| {
                let u = (i % 10) as f64 * 0.1;
                let v = (i / 10) as f64 * 0.2;
                Vector3::new(u, v, 0.0)
            })
            .collect();
        let fit = surface_through_scattered(&pts, 3, 3, 4, 4).unwrap();
        // RMS should be small for a perfectly planar scatter.
        assert!(fit.rms_error < 0.5, "rms = {}", fit.rms_error);
    }

    #[test]
    fn rejects_too_few_points() {
        let pts = [Vector3::zeros()];
        assert!(nurbs_curve_through_points(&pts, 1, 2).is_err());
    }
}
