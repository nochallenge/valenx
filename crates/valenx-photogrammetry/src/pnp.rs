//! Camera resectioning by PnP: registering a new view from 2D–3D
//! correspondences (Stage 4 — the heart of the incremental mapper).
//!
//! Stage 3 ([`crate::twoview`]) bootstraps a reconstruction from an image
//! *pair*: two views, their relative pose `(R, t)`, and a first set of
//! triangulated 3-D points. The **incremental mapper** then grows that
//! reconstruction one image at a time. The core operation each time a new
//! image is added is **resectioning** — recovering the new camera's pose from
//! correspondences between its 2-D feature points and the already-triangulated
//! 3-D scene points:
//!
//! ```text
//!  (already in the model)      (new image)
//!   3-D points  Xᵢ      ↔       2-D features  xᵢ
//!                       │
//!                       └─ PnP (resectioning) ──►  new camera pose (R, t)
//! ```
//!
//! This is the **Perspective-n-Point (PnP)** problem: given `n` known 3-D
//! points, their observed image projections, and the camera intrinsics `K`,
//! find the camera pose `(R, t)` with `X_cam = R·X_world + t`. Unlike the
//! two-view stage, PnP recovers a pose at **metric scale** — the 3-D points
//! already live in the model's coordinate frame, so the new camera slots into
//! the *same* frame (no separate scale ambiguity; see the honesty note below).
//!
//! ## Method — linear DLT-PnP (Hartley & Zisserman §7.1)
//!
//! We solve the **calibrated** PnP linearly by the Direct Linear Transform.
//! First normalize each image point to a calibrated ray `x̂ = K⁻¹ x`. The
//! projection of a world point is `x̂ ∝ [R | t] X` (in homogeneous form), so
//! the cross-product constraint `x̂ᵢ × ([R | t] Xᵢ) = 0` holds exactly. Each
//! correspondence contributes two independent rows to a homogeneous
//! `2n × 12` system `A p = 0`, where `p` is the row-major vectorization of the
//! `3×4` matrix `M = [R | t]`. The right singular vector of `A`'s smallest
//! singular value is reshaped back to a `3×4` `M`.
//!
//! The recovered `M` is only *approximately* of the form `[R | t]` (its
//! left `3×3` block is some scaled, slightly non-orthogonal matrix), so we:
//!
//! 1. **Orthonormalize** the `3×3` block by its own SVD `B = U Σ Vᵀ` →
//!    `R = U Vᵀ` (the nearest proper rotation in Frobenius norm), forcing
//!    `det R = +1` by flipping the sign of the smallest-singular-value
//!    direction if the SVD handed back a reflection.
//! 2. **Rescale** the translation by the same factor the orthonormalization
//!    removed from the block — the mean singular value `s̄ = (σ₁+σ₂+σ₃)/3` —
//!    so `t` stays metrically consistent with the cleaned `R`
//!    (`t = t_raw / s̄`).
//! 3. **Fix the overall sign** so the reconstructed points sit *in front of*
//!    the camera (positive `z = (R·X + t)_z`): the homogeneous solution is
//!    defined only up to sign, and the two signs correspond to the scene
//!    being in front of vs. behind the camera.
//!
//! See R. Hartley & A. Zisserman, *Multiple View Geometry in Computer
//! Vision*, 2nd ed., §7.1 (the DLT camera-resectioning algorithm).
//!
//! ## Robustness — RANSAC over reprojection error
//!
//! Real 2D–3D match sets are contaminated by mismatches. [`solve_pnp_ransac`]
//! wraps the linear solver in RANSAC (Fischler & Bolles 1981): it repeatedly
//! samples a minimal set of **6** correspondences, fits a pose with
//! [`solve_pnp`], scores **all** correspondences by **reprojection error**
//! (the pixel distance between the observed point and the projection of its
//! 3-D point through `K[R | t]`), keeps the largest consensus set, and refits
//! [`solve_pnp`] on those inliers. The sampling is deterministic for a given
//! [`PnpRansacParams::seed`] (the in-crate `SplitMix64`, no `rand`
//! dependency).
//!
//! ## Honesty notes
//!
//! - **Linear, not optimal.** DLT-PnP minimizes an *algebraic* residual, not
//!   the geometric reprojection error, so its accuracy is below that of EPnP
//!   or an iterative (Levenberg–Marquardt) refinement on noisy data. It is a
//!   solid, fast initializer — exactly what the incremental mapper needs
//!   before the global **bundle adjustment** stage (built later) polishes all
//!   poses and points together. On *clean* correspondences it is essentially
//!   exact (recovered to roundoff), which the tests assert.
//! - **Six-point minimum.** The `2n × 12` system has 11 degrees of freedom
//!   (the `3×4` `M` up to scale), so `n ≥ 6` correspondences (12 rows) are
//!   required. The classic *minimal* calibrated PnP is P3P (3 points, up to
//!   four discrete solutions); the linear DLT instead needs six but yields a
//!   single closed-form answer with no disambiguation step. We use six.
//! - **Coplanar degeneracy.** If all the 3-D points are **coplanar** (or
//!   collinear), the DLT system is rank-deficient and the linear method
//!   breaks down (a planar scene needs a homography-based pose instead). This
//!   stage targets the general non-coplanar case; a coplanar configuration
//!   returns a meaningless pose or, when the rank collapse defeats the SVD /
//!   sign / depth checks, [`None`]. Callers feeding planar scenes should use a
//!   planar-PnP path (not implemented here).
//! - **Scale.** Because the 3-D points are supplied in a fixed world frame,
//!   the recovered `(R, t)` is at that frame's metric scale — there is *no*
//!   free scale factor here (contrast the two-view stage, where translation
//!   was a unit direction only).

use crate::twoview::CameraIntrinsics;
use nalgebra::{DMatrix, Matrix3, Matrix3x4, Vector3};

/// Minimal number of 2D–3D correspondences for the linear DLT-PnP solver.
///
/// The `3×4` projection block `M = [R | t]` has 11 degrees of freedom (12
/// entries, homogeneous), and each correspondence yields two independent
/// equations, so six points (12 equations) are the linear minimum.
const MIN_POINTS: usize = 6;

/// A recovered camera pose: the rotation and translation of a view.
///
/// The convention is **`X_cam = R · X_world + t`** — `rotation` maps a point
/// from the world (model) frame into this camera's frame, and `translation`
/// is the world origin's position in the camera frame. The camera's
/// projection matrix is therefore `P = K [R | t]`, and the camera centre in
/// world coordinates is `C = −Rᵀ t`.
///
/// `rotation` is a proper rotation (`Rᵀ R = I`, `det R = +1`); see
/// [`solve_pnp`] for how it is orthonormalized out of the raw linear solution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CameraPose {
    /// Rotation `R` mapping world points into the camera frame
    /// (`X_cam = R · X_world + t`). A proper rotation: `Rᵀ R = I`,
    /// `det R = +1`.
    pub rotation: Matrix3<f64>,
    /// Translation `t` (the world origin expressed in the camera frame);
    /// `X_cam = R · X_world + t`. At the metric scale of the world points.
    pub translation: Vector3<f64>,
}

/// The result of a robust (RANSAC) PnP fit: the recovered pose and the
/// indices of the input correspondences that agreed with it.
#[derive(Debug, Clone)]
pub struct PnpResult {
    /// The recovered camera pose, refit on the consensus (inlier) set.
    pub pose: CameraPose,
    /// Indices (into the input `points3d` / `points2d` slices) of the
    /// correspondences whose reprojection error fell below the threshold —
    /// the consensus set that supports [`Self::pose`].
    pub inliers: Vec<usize>,
}

/// Parameters controlling the [`solve_pnp_ransac`] robust search.
#[derive(Debug, Clone, Copy)]
pub struct PnpRansacParams {
    /// Number of RANSAC iterations (random minimal samples to try). Each
    /// iteration draws six correspondences and fits a pose; more iterations
    /// raise the chance of hitting an all-inlier sample, at linear cost.
    pub iterations: usize,
    /// Inlier threshold on the **reprojection error**, in pixels. A
    /// correspondence counts as an inlier when the pixel distance between its
    /// observed 2-D point and the projection of its 3-D point through the
    /// candidate pose is below this value.
    pub threshold: f64,
    /// Seed for the deterministic `SplitMix64` sample-index PRNG, so a given
    /// input + params always yields the same result (reproducible builds /
    /// tests). No external `rand` dependency is used.
    pub seed: u64,
}

impl Default for PnpRansacParams {
    /// Sensible defaults: 1000 iterations, a 4 px reprojection threshold, and
    /// a fixed seed.
    fn default() -> Self {
        Self {
            iterations: 1000,
            threshold: 4.0,
            seed: 0x9E37_79B9_7F4A_7C15,
        }
    }
}

/// Project a world point through `K, R, t` to a pixel, guarding `z ≈ 0`.
///
/// Computes `x_cam = R · X + t`, then the pixel `(fx·x/z + cx, fy·y/z + cy)`
/// via `K`. When the camera-frame depth `z` is at (or near) zero — the point
/// lies on the camera's principal plane, where the perspective divide is
/// undefined — this returns [`None`] rather than dividing by zero. A point
/// that reprojects to `None` is simply never counted as a RANSAC inlier.
#[must_use]
pub fn project_point(
    k: &CameraIntrinsics,
    r: &Matrix3<f64>,
    t: &Vector3<f64>,
    x: &Vector3<f64>,
) -> Option<(f64, f64)> {
    let cam = r * x + t;
    // Guard the perspective divide: a point on/near the principal plane
    // (z ≈ 0) has no finite image; report it as unprojectable.
    if cam.z.abs() < 1e-12 {
        return None;
    }
    let u = k.fx * (cam.x / cam.z) + k.cx;
    let v = k.fy * (cam.y / cam.z) + k.cy;
    Some((u, v))
}

/// Solve **calibrated PnP** linearly by the DLT (Hartley & Zisserman §7.1):
/// recover the camera pose `(R, t)` from `n ≥ 6` 2D–3D correspondences and
/// the intrinsics `K`.
///
/// `points3d[i]` is a world-frame 3-D point and `points2d[i]` its observed
/// pixel `(u, v)` in the image; the two slices are matched element-wise. The
/// returned [`CameraPose`] follows the `X_cam = R · X_world + t` convention.
///
/// # Method
///
/// 1. Normalize each pixel to a calibrated ray `x̂ = K⁻¹ x`.
/// 2. From `x̂ᵢ × ([R | t] Xᵢ) = 0`, build the homogeneous `2n × 12` system
///    `A p = 0` (two independent rows per correspondence; `p` is the
///    row-major vectorization of `M = [R | t]`).
/// 3. SVD of `A`; the right singular vector of the smallest singular value is
///    reshaped into the `3×4` `M`.
/// 4. **Orthonormalize** `M`'s left `3×3` block via its SVD (`R = U Vᵀ`,
///    forced to `det = +1`), and **rescale** `t` by the block's mean singular
///    value so it stays metrically consistent with the cleaned `R`.
/// 5. **Fix the sign** so the points have positive camera-frame depth.
///
/// # Returns
///
/// [`None`] when there are fewer than six correspondences, when the two
/// slices differ in length, or when the configuration is degenerate (a failed
/// SVD, a rank-deficient / coplanar system that yields a non-finite or
/// non-invertible block, or a solution with no consistent positive-depth
/// sign). Never panics.
///
/// See the [module docs](self) for the linear-vs-iterative accuracy note, the
/// six-point minimum, and the coplanar-degeneracy caveat.
#[must_use]
pub fn solve_pnp(
    points3d: &[Vector3<f64>],
    points2d: &[(f64, f64)],
    k: &CameraIntrinsics,
) -> Option<CameraPose> {
    if points3d.len() != points2d.len() || points3d.len() < MIN_POINTS {
        return None;
    }
    let n = points3d.len();

    // K⁻¹ to normalize pixels to calibrated rays. With zero skew this is
    // closed-form, but we guard the (degenerate) singular-K case anyway.
    let k_inv = k.matrix().try_inverse()?;

    // Assemble the 2n × 12 DLT system. For a calibrated ray
    // x̂ = (a, b, c) and M = [m₀ᵀ; m₁ᵀ; m₂ᵀ] (rows of the 3×4 matrix), the
    // cross product x̂ × (M X) = 0 gives three rows, of which two are
    // independent; we use the first two (the standard DLT choice):
    //
    //   row r0:  b·(m₂·X) − c·(m₁·X) = 0
    //   row r1:  c·(m₀·X) − a·(m₂·X) = 0
    //
    // Each X is the homogeneous 4-vector (X, Y, Z, 1); the unknown p stacks
    // m₀, m₁, m₂ (12 entries, row-major).
    let mut a = DMatrix::<f64>::zeros(2 * n, 12);
    for i in 0..n {
        let p3 = points3d[i];
        let (u, v) = points2d[i];
        // Calibrated ray x̂ = K⁻¹ [u, v, 1]ᵀ.
        let ray = k_inv * Vector3::new(u, v, 1.0);
        let (xa, xb, xc) = (ray.x, ray.y, ray.z);
        // Homogeneous world point as four scalars.
        let xw = [p3.x, p3.y, p3.z, 1.0];

        let r0 = 2 * i;
        let r1 = 2 * i + 1;
        for j in 0..4 {
            // Columns 0..4 = m₀, 4..8 = m₁, 8..12 = m₂.
            // Row r0:  b·m₂·X − c·m₁·X.
            a[(r0, 4 + j)] = -xc * xw[j]; // −c · (m₁ · X)
            a[(r0, 8 + j)] = xb * xw[j]; //  b · (m₂ · X)
                                         // Row r1:  c·m₀·X − a·m₂·X.
            a[(r1, j)] = xc * xw[j]; //  c · (m₀ · X)
            a[(r1, 8 + j)] = -xa * xw[j]; // −a · (m₂ · X)
        }
    }

    // Smallest-singular-value right vector = last row of Vᵀ.
    let svd = a.svd(false, true);
    let v_t = svd.v_t?;

    // Degeneracy guard. A well-posed calibrated PnP constrains 11 of the 12
    // DOF of the homogeneous `M`, so the system has exactly a **one**-
    // dimensional solution null space: only the smallest singular value is
    // ≈ 0, and the second-smallest is clearly non-zero. If the
    // second-smallest singular value is also negligible relative to the
    // largest, the null space is higher-dimensional — the configuration is
    // rank-deficient (all-coincident, collinear, or coplanar points) and the
    // "solution" the SVD returns is an arbitrary null-space direction, i.e.
    // meaningless. Reject such configurations rather than return a bogus pose.
    // (`singular_values` are produced in descending order; for a `2n×12`
    // matrix there are 12 of them.)
    let sv_a = &svd.singular_values;
    let smax = sv_a[0];
    let second_smallest = sv_a[sv_a.len() - 2];
    // `smax <= 0` (or NaN) means a collapsed system; a small relative
    // second-smallest means a >1-D null space (degenerate configuration).
    if smax <= 0.0 || !smax.is_finite() || second_smallest / smax < 1e-9 {
        return None;
    }

    // The 12-vector solution (row-major M = [R | t]).
    let p = v_t.row(v_t.nrows() - 1);

    // Reshape into the 3×4 M (row-major: p = [m₀ | m₁ | m₂]).
    let mut m = Matrix3x4::<f64>::zeros();
    for r in 0..3 {
        for c in 0..4 {
            m[(r, c)] = p[r * 4 + c];
        }
    }

    // Split into the raw rotation block B and raw translation column.
    let b = m.fixed_view::<3, 3>(0, 0).into_owned();
    let t_raw = m.fixed_view::<3, 1>(0, 3).into_owned();

    // Bail out on a non-finite block (a fully degenerate / rank-collapsed
    // system can produce NaNs through the SVD).
    if b.iter().any(|x| !x.is_finite()) || t_raw.iter().any(|x| !x.is_finite()) {
        return None;
    }

    // Orthonormalize the rotation block: nearest proper rotation R = U Vᵀ.
    let bsvd = b.svd(true, true);
    let (Some(u), Some(vt)) = (bsvd.u, bsvd.v_t) else {
        return None;
    };
    let sv = bsvd.singular_values;
    // Mean singular value = the isotropic scale the block carried. This is the
    // factor orthonormalization strips from R, so we divide t by it to keep
    // t metrically consistent with the cleaned R. Guard a collapsed block.
    let s_mean = (sv[0] + sv[1] + sv[2]) / 3.0;
    // Reject a collapsed block (and NaN): the scale must be safely positive.
    if s_mean < 1e-12 || !s_mean.is_finite() {
        return None;
    }

    let mut r = u * vt;
    let mut t = t_raw / s_mean;
    // Force a proper rotation (det = +1). If the SVD product is a reflection,
    // flipping the sign of the last column of U flips det(R) while keeping it
    // orthogonal — the standard "nearest rotation" correction.
    if r.determinant() < 0.0 {
        let mut u_fixed = u;
        let last = u_fixed.ncols() - 1;
        for row in 0..3 {
            u_fixed[(row, last)] = -u_fixed[(row, last)];
        }
        r = u_fixed * vt;
    }

    // Fix the overall sign (the homogeneous solution is defined up to ±): pick
    // the sign that puts the scene in front of the camera (positive mean
    // depth z = (R·X + t)_z over the correspondences). Flipping the sign of
    // BOTH R and t is a valid alternative solution of A p = 0; only one of the
    // two has the points in front. Note flipping R's sign would break
    // det = +1, so we instead test depth under the proper R and, if the scene
    // is behind, negate t and rotate by R's "twin" — but the clean way that
    // preserves det = +1 is: the sign ambiguity of the DLT vector flips
    // (R, t) → (−R, −t); −R is improper, so the physically meaningful
    // disambiguation is carried entirely by t's sign together with the global
    // vector sign. We therefore evaluate depth for the current (R, t) and for
    // (R, t) with the *whole* solution negated re-orthonormalized — which, for
    // the rotation, returns the same R — leaving t's sign as the only free
    // choice. Concretely: choose the sign of the solution vector that yields
    // positive mean depth.
    let mean_depth = |rot: &Matrix3<f64>, tr: &Vector3<f64>| -> f64 {
        let mut acc = 0.0;
        for p3 in points3d {
            acc += (rot * p3 + tr).z;
        }
        acc / n as f64
    };

    // Negating the entire 12-vector p flips both B and t_raw. Re-running the
    // orthonormalization on −B yields the same R when det was already fixed,
    // but t flips sign. So the two physical candidates are (R, +t) and
    // (R, −t) — exactly the global-sign freedom. Pick the one in front.
    if mean_depth(&r, &t) < 0.0 {
        t = -t;
        // With t negated we must also negate the rotation's handling of the
        // raw block sign so depth is consistent; re-derive R from −B.
        let nb = -b;
        let nbsvd = nb.svd(true, true);
        if let (Some(nu), Some(nvt)) = (nbsvd.u, nbsvd.v_t) {
            let mut nr = nu * nvt;
            if nr.determinant() < 0.0 {
                let mut nu_fixed = nu;
                let last = nu_fixed.ncols() - 1;
                for row in 0..3 {
                    nu_fixed[(row, last)] = -nu_fixed[(row, last)];
                }
                nr = nu_fixed * nvt;
            }
            r = nr;
        }
    }

    // Final sanity: everything finite and R still a proper rotation.
    if r.iter().any(|x| !x.is_finite())
        || t.iter().any(|x| !x.is_finite())
        || (r.determinant() - 1.0).abs() > 1e-6
    {
        return None;
    }

    Some(CameraPose {
        rotation: r,
        translation: t,
    })
}

/// Robustly solve PnP with **RANSAC** over reprojection error: register a new
/// view from 2D–3D correspondences that may contain gross outliers.
///
/// `points3d` / `points2d` are matched element-wise (a world point and its
/// observed pixel). The algorithm repeatedly samples six correspondences,
/// fits a pose with [`solve_pnp`], scores **all** correspondences by
/// reprojection error (pixel distance below [`PnpRansacParams::threshold`]),
/// keeps the largest consensus set, and finally refits [`solve_pnp`] on that
/// inlier set. Sampling is deterministic for a given
/// [`PnpRansacParams::seed`].
///
/// # Returns
///
/// [`None`] when there are fewer than six correspondences, the slices differ
/// in length, or no sampled minimal set ever yields a pose with at least six
/// supporting inliers (a hopelessly degenerate or outlier-dominated set).
/// Never panics.
///
/// The pose accuracy is that of the linear [`solve_pnp`] refit on the
/// inliers — a strong initializer for the later bundle-adjustment stage, not
/// itself a maximum-likelihood estimate (see the [module docs](self)).
#[must_use]
pub fn solve_pnp_ransac(
    points3d: &[Vector3<f64>],
    points2d: &[(f64, f64)],
    k: &CameraIntrinsics,
    params: &PnpRansacParams,
) -> Option<PnpResult> {
    if points3d.len() != points2d.len() || points3d.len() < MIN_POINTS {
        return None;
    }
    let n = points3d.len();
    let thresh2 = params.threshold * params.threshold;

    let mut rng = SplitMix64::new(params.seed);
    let mut best_inliers: Vec<usize> = Vec::new();

    let iterations = params.iterations.max(1);
    for _ in 0..iterations {
        // Draw six distinct correspondence indices for the minimal fit.
        let sample = sample_indices(&mut rng, n, MIN_POINTS);
        let s3d: Vec<Vector3<f64>> = sample.iter().map(|&i| points3d[i]).collect();
        let s2d: Vec<(f64, f64)> = sample.iter().map(|&i| points2d[i]).collect();

        let Some(pose) = solve_pnp(&s3d, &s2d, k) else {
            continue;
        };

        // Score ALL correspondences by reprojection error.
        let inliers = consensus(points3d, points2d, k, &pose, thresh2);
        if inliers.len() > best_inliers.len() {
            best_inliers = inliers;
        }
    }

    // Need at least the minimal count to refit a pose.
    if best_inliers.len() < MIN_POINTS {
        return None;
    }

    // Refit on the full consensus set for a better pose.
    let in3d: Vec<Vector3<f64>> = best_inliers.iter().map(|&i| points3d[i]).collect();
    let in2d: Vec<(f64, f64)> = best_inliers.iter().map(|&i| points2d[i]).collect();
    let refined = solve_pnp(&in3d, &in2d, k)?;

    // Recompute the consensus under the refined pose so the reported inlier
    // set is consistent with the returned pose (the refit can only help).
    let final_inliers = consensus(points3d, points2d, k, &refined, thresh2);
    let inliers = if final_inliers.len() >= best_inliers.len() {
        final_inliers
    } else {
        best_inliers
    };

    Some(PnpResult {
        pose: refined,
        inliers,
    })
}

/// Indices of correspondences whose squared reprojection error under `pose`
/// is below `thresh2`. A point that projects behind/through the camera
/// (`project_point` returns `None`) is not an inlier.
fn consensus(
    points3d: &[Vector3<f64>],
    points2d: &[(f64, f64)],
    k: &CameraIntrinsics,
    pose: &CameraPose,
    thresh2: f64,
) -> Vec<usize> {
    let mut inliers = Vec::new();
    for i in 0..points3d.len() {
        let Some((pu, pv)) = project_point(k, &pose.rotation, &pose.translation, &points3d[i])
        else {
            continue;
        };
        let (ou, ov) = points2d[i];
        let du = pu - ou;
        let dv = pv - ov;
        if du * du + dv * dv < thresh2 {
            inliers.push(i);
        }
    }
    inliers
}

/// Draw `k` distinct indices from `0..n` using the in-crate `SplitMix64`
/// PRNG (partial Fisher–Yates over a scratch index list). Requires `k <= n`,
/// which the single caller guarantees (`n >= MIN_POINTS == k`).
fn sample_indices(rng: &mut SplitMix64, n: usize, k: usize) -> Vec<usize> {
    let mut pool: Vec<usize> = (0..n).collect();
    for i in 0..k {
        // Pick j uniformly in [i, n), swap into position i.
        let j = i + (rng.next_u64() as usize) % (n - i);
        pool.swap(i, j);
    }
    pool.truncate(k);
    pool
}

/// Minimal SplitMix64 PRNG, mirroring the one in [`crate::geometry`] /
/// [`crate::descriptor`].
///
/// Used here only to draw deterministic RANSAC minimal-sample indices, so the
/// crate needs no `rand` dependency and results are reproducible across runs
/// and machines. Not used for any security purpose.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// Seed the generator.
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Next raw 64-bit output.
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_intrinsics() -> CameraIntrinsics {
        CameraIntrinsics::new(800.0, 800.0, 320.0, 240.0)
    }

    /// A known, generic camera pose (yaw + pitch rotation, generic
    /// translation) used as the ground truth across the tests.
    fn known_pose() -> (Matrix3<f64>, Vector3<f64>) {
        let yaw = 0.30_f64;
        let pitch = -0.17_f64;
        let ry = Matrix3::new(
            yaw.cos(),
            0.0,
            yaw.sin(),
            0.0,
            1.0,
            0.0,
            -yaw.sin(),
            0.0,
            yaw.cos(),
        );
        let rx = Matrix3::new(
            1.0,
            0.0,
            0.0,
            0.0,
            pitch.cos(),
            -pitch.sin(),
            0.0,
            pitch.sin(),
            pitch.cos(),
        );
        let r = ry * rx;
        // Translation places the world origin a few units in front of the
        // camera so the scene below has positive depth.
        let t = Vector3::new(0.4, -0.25, 6.0);
        (r, t)
    }

    /// A fixed cloud of **non-coplanar** 3-D world points spanning a range of
    /// depths (so the DLT system is well-conditioned).
    fn scene_points() -> Vec<Vector3<f64>> {
        let raw = [
            (-0.6, -0.4, 0.3),
            (0.5, -0.3, -0.4),
            (0.2, 0.6, 0.7),
            (-0.4, 0.45, -0.6),
            (0.7, 0.5, 0.2),
            (-0.55, -0.25, -0.5),
            (0.3, -0.6, 0.55),
            (-0.15, 0.2, 0.9),
            (0.1, -0.1, -0.8),
            (0.6, 0.1, 0.45),
            (-0.7, 0.35, 0.15),
            (0.25, 0.3, -0.7),
        ];
        raw.iter().map(|&(x, y, z)| Vector3::new(x, y, z)).collect()
    }

    /// Project a world point through `K, R, t` to a pixel (test ground-truth
    /// projector; mirrors [`project_point`] but without the `None` guard).
    fn project(
        k: &Matrix3<f64>,
        r: &Matrix3<f64>,
        t: &Vector3<f64>,
        x: &Vector3<f64>,
    ) -> (f64, f64) {
        let cam = r * x + t;
        let px = k * cam;
        (px.x / px.z, px.y / px.z)
    }

    // ----- Test 1: solve_pnp recovers a known pose exactly. -----
    #[test]
    fn solve_pnp_recovers_known_pose() {
        let (r_true, t_true) = known_pose();
        let k = test_intrinsics();
        let km = k.matrix();
        let pts = scene_points();

        let img: Vec<(f64, f64)> = pts
            .iter()
            .map(|x| project(&km, &r_true, &t_true, x))
            .collect();

        let pose = solve_pnp(&pts, &img, &k).expect("PnP should solve clean ≥6-point input");

        let r_err = (pose.rotation - r_true).norm();
        assert!(r_err < 1e-6, "recovered R off by {r_err}");
        let t_err = (pose.translation - t_true).norm();
        assert!(
            t_err < 1e-6,
            "recovered t {:?} off truth {t_true:?} by {t_err}",
            pose.translation
        );
    }

    // ----- Test 2: the recovered R is a proper orthonormal rotation. -----
    #[test]
    fn recovered_rotation_is_orthonormal() {
        let (r_true, t_true) = known_pose();
        let k = test_intrinsics();
        let km = k.matrix();
        let pts = scene_points();
        let img: Vec<(f64, f64)> = pts
            .iter()
            .map(|x| project(&km, &r_true, &t_true, x))
            .collect();

        let pose = solve_pnp(&pts, &img, &k).expect("PnP should solve");

        // RᵀR ≈ I.
        let rtr = pose.rotation.transpose() * pose.rotation;
        let ortho_err = (rtr - Matrix3::identity()).norm();
        assert!(ortho_err < 1e-9, "RᵀR deviates from I by {ortho_err}");
        // det R ≈ +1.
        let det = pose.rotation.determinant();
        assert!((det - 1.0).abs() < 1e-9, "det R = {det}, expected +1");
    }

    // ----- Test 3: RANSAC recovers the pose + inliers despite ~30% outliers. -----
    #[test]
    fn ransac_is_robust_to_outliers() {
        let (r_true, t_true) = known_pose();
        let k = test_intrinsics();
        let km = k.matrix();

        // A larger clean point cloud for a meaningful inlier/outlier split.
        let mut clean3d: Vec<Vector3<f64>> = scene_points();
        // Add a few more so we have 20 inliers.
        let extra = [
            (-0.3, -0.55, 0.6),
            (0.45, 0.2, -0.35),
            (-0.2, -0.2, 0.75),
            (0.15, 0.5, 0.5),
            (-0.65, 0.1, -0.45),
            (0.55, -0.4, 0.25),
            (0.05, 0.35, -0.65),
            (-0.45, -0.35, 0.4),
        ];
        for &(x, y, z) in &extra {
            clean3d.push(Vector3::new(x, y, z));
        }
        let n_clean = clean3d.len();

        // True projections for the clean points.
        let mut pts3d = clean3d.clone();
        let mut pts2d: Vec<(f64, f64)> = clean3d
            .iter()
            .map(|x| project(&km, &r_true, &t_true, x))
            .collect();

        // Mark which indices are the true inliers.
        let true_inliers: std::collections::BTreeSet<usize> = (0..n_clean).collect();

        // Add ~30% gross outliers: real 3-D points paired with WRONG pixels
        // (a deterministic large pixel offset puts them far from any correct
        // reprojection).
        let n_out = (n_clean as f64 * 0.3).ceil() as usize; // ~6 of 20
        let mut rng = SplitMix64::new(0x1234_5678_9ABC_DEF0);
        for _ in 0..n_out {
            // A genuine-looking 3-D point.
            let gx = ((rng.next_u64() % 1000) as f64 / 1000.0) - 0.5;
            let gy = ((rng.next_u64() % 1000) as f64 / 1000.0) - 0.5;
            let gz = ((rng.next_u64() % 1000) as f64 / 1000.0) - 0.5;
            pts3d.push(Vector3::new(gx, gy, gz));
            // A wildly wrong pixel (offset by hundreds of px from any sane
            // projection), guaranteed to exceed the reprojection threshold.
            let bu = (rng.next_u64() % 640) as f64;
            let bv = (rng.next_u64() % 480) as f64;
            pts2d.push((bu + 1000.0, bv + 1000.0));
        }

        let params = PnpRansacParams {
            iterations: 500,
            threshold: 2.0,
            seed: 0xABCD_0123_4567_89AB,
        };
        let result =
            solve_pnp_ransac(&pts3d, &pts2d, &k, &params).expect("RANSAC PnP should succeed");

        // Pose matches the truth.
        let r_err = (result.pose.rotation - r_true).norm();
        assert!(r_err < 1e-5, "RANSAC recovered R off by {r_err}");
        let t_err = (result.pose.translation - t_true).norm();
        assert!(t_err < 1e-5, "RANSAC recovered t off by {t_err}");

        // The inlier set is essentially the true clean set: every reported
        // inlier is a true one, and we recovered (nearly) all of them.
        let found: std::collections::BTreeSet<usize> = result.inliers.iter().copied().collect();
        for &idx in &found {
            assert!(
                true_inliers.contains(&idx),
                "index {idx} flagged inlier but is a planted outlier"
            );
        }
        assert!(
            found.len() >= n_clean - 1,
            "expected ~{n_clean} inliers, got {}",
            found.len()
        );
    }

    // ----- Test 4: the recovered pose reprojects inliers within threshold. -----
    #[test]
    fn reprojection_within_threshold() {
        let (r_true, t_true) = known_pose();
        let k = test_intrinsics();
        let km = k.matrix();
        let pts = scene_points();
        let img: Vec<(f64, f64)> = pts
            .iter()
            .map(|x| project(&km, &r_true, &t_true, x))
            .collect();

        let pose = solve_pnp(&pts, &img, &k).expect("PnP should solve");

        // Every point reprojects to within a tight pixel tolerance under the
        // recovered pose (clean data → essentially exact).
        for (x, &(ou, ov)) in pts.iter().zip(img.iter()) {
            let (pu, pv) =
                project_point(&k, &pose.rotation, &pose.translation, x).expect("point in front");
            let err = ((pu - ou).powi(2) + (pv - ov).powi(2)).sqrt();
            assert!(err < 1e-6, "reprojection error {err} px exceeds tolerance");
        }
    }

    // ----- Test 5a: degenerate counts / mismatched lengths -> None, no panic. -----
    #[test]
    fn too_few_points_returns_none() {
        let k = test_intrinsics();
        let (r_true, t_true) = known_pose();
        let km = k.matrix();
        let pts = scene_points();
        let img: Vec<(f64, f64)> = pts
            .iter()
            .map(|x| project(&km, &r_true, &t_true, x))
            .collect();

        // Fewer than six points.
        assert!(
            solve_pnp(&pts[..5], &img[..5], &k).is_none(),
            "5 points must yield None"
        );
        assert!(solve_pnp(&[], &[], &k).is_none(), "empty must yield None");

        // Mismatched slice lengths.
        assert!(
            solve_pnp(&pts, &img[..pts.len() - 1], &k).is_none(),
            "length mismatch must yield None"
        );

        // RANSAC with too few points likewise None.
        assert!(
            solve_pnp_ransac(&pts[..5], &img[..5], &k, &PnpRansacParams::default()).is_none(),
            "RANSAC with <6 points must yield None"
        );
    }

    // ----- Test 5b: all-coincident 3-D points -> None, no panic. -----
    #[test]
    fn coincident_points_return_none() {
        let k = test_intrinsics();
        // Eight identical world points (rank-deficient: cannot fix a pose).
        let p = Vector3::new(0.1, -0.2, 0.5);
        let pts3d = vec![p; 8];
        // Project them through some pose; all pixels coincide too.
        let (r_true, t_true) = known_pose();
        let km = k.matrix();
        let px = project(&km, &r_true, &t_true, &p);
        let pts2d = vec![px; 8];

        // Must not panic; the rank-deficient system yields no valid pose.
        let out = solve_pnp(&pts3d, &pts2d, &k);
        assert!(out.is_none(), "all-coincident points must yield None");
    }

    // ----- Test 5b': coplanar 3-D points are degenerate -> None, no panic. ---
    #[test]
    fn coplanar_points_return_none() {
        // Honest caveat from the module docs made a tested guarantee: a
        // coplanar point set (here all on the world plane z = 0) is
        // rank-deficient for the linear DLT and must be rejected, not return a
        // bogus pose. (A planar scene needs a homography-based pose, which is
        // out of scope for this stage.)
        let k = test_intrinsics();
        let (r_true, t_true) = known_pose();
        let km = k.matrix();
        let planar = [
            (-0.6, -0.4, 0.0),
            (0.5, -0.3, 0.0),
            (0.2, 0.6, 0.0),
            (-0.4, 0.45, 0.0),
            (0.7, 0.5, 0.0),
            (-0.55, -0.25, 0.0),
            (0.3, -0.6, 0.0),
            (-0.15, 0.2, 0.0),
        ];
        let pts3d: Vec<Vector3<f64>> = planar
            .iter()
            .map(|&(x, y, z)| Vector3::new(x, y, z))
            .collect();
        let pts2d: Vec<(f64, f64)> = pts3d
            .iter()
            .map(|x| project(&km, &r_true, &t_true, x))
            .collect();

        // Must not panic; the rank-deficient planar system yields no valid pose.
        assert!(
            solve_pnp(&pts3d, &pts2d, &k).is_none(),
            "coplanar points must yield None (planar degeneracy)"
        );
    }

    // ----- Test 5c: project_point with z ≈ 0 does not divide by zero. -----
    #[test]
    fn project_point_guards_zero_depth() {
        let k = test_intrinsics();
        let r = Matrix3::identity();
        let t = Vector3::zeros();
        // A point exactly on the camera centre / principal plane (z = 0).
        let on_plane = Vector3::new(1.0, 1.0, 0.0);
        assert!(
            project_point(&k, &r, &t, &on_plane).is_none(),
            "z≈0 point must return None, not divide by zero"
        );
        // A point with tiny but non-zero depth still projects to a finite px.
        let near = Vector3::new(0.01, 0.01, 1e-6);
        if let Some((u, v)) = project_point(&k, &r, &t, &near) {
            assert!(
                u.is_finite() && v.is_finite(),
                "near-plane px must be finite"
            );
        }
        // A normal point in front projects fine.
        let front = Vector3::new(0.0, 0.0, 5.0);
        let (u, v) = project_point(&k, &r, &t, &front).expect("front point projects");
        assert!((u - k.cx).abs() < 1e-9 && (v - k.cy).abs() < 1e-9);
    }

    // ----- Bonus: RANSAC returns a non-inlier for a z≈0 / behind point. -----
    #[test]
    fn ransac_excludes_unprojectable_points() {
        // Build a clean set, then add one correspondence whose 3-D point sits
        // behind the camera for the true pose; it must not be an inlier.
        let (r_true, t_true) = known_pose();
        let k = test_intrinsics();
        let km = k.matrix();
        let mut pts3d = scene_points();
        let mut pts2d: Vec<(f64, f64)> = pts3d
            .iter()
            .map(|x| project(&km, &r_true, &t_true, x))
            .collect();

        // A point far behind the camera, paired with a plausible pixel.
        let behind = Vector3::new(0.0, 0.0, -100.0);
        pts3d.push(behind);
        pts2d.push((320.0, 240.0));

        let result = solve_pnp_ransac(&pts3d, &pts2d, &k, &PnpRansacParams::default())
            .expect("RANSAC should still solve from the clean majority");
        // The behind-camera correspondence (last index) is not an inlier.
        let last = pts3d.len() - 1;
        assert!(
            !result.inliers.contains(&last),
            "behind-camera point must be excluded from inliers"
        );
    }
}
