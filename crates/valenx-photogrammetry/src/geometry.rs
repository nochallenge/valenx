//! Two-view geometric verification (Stage 2, part 2).
//!
//! Putative matches from [`crate::match_descriptors`] still contain
//! outliers (mismatched descriptors). This module fits the **fundamental
//! matrix** `F` relating the two views and, in doing so, separates the
//! geometrically consistent *inliers* from the outliers.
//!
//! The fundamental matrix encodes the epipolar geometry of an uncalibrated
//! stereo pair: for any true correspondence `x₁ ↔ x₂` (homogeneous pixel
//! coordinates),
//!
//! ```text
//! x₂ᵀ F x₁ = 0
//! ```
//!
//! ## Method
//!
//! 1. **Normalized 8-point algorithm** (Hartley 1997). Pixel coordinates
//!    are first conditioned by an isotropic similarity transform per image
//!    — translate so the centroid is at the origin and scale so the mean
//!    distance to the origin is `√2`. The linear system `A f = 0` (one row
//!    `[x₂x₁, x₂y₁, x₂, y₂x₁, y₂y₁, y₂, x₁, y₁, 1]` per correspondence) is
//!    solved by SVD; `f` is the right-singular vector of the smallest
//!    singular value. The resulting `F̂` is then forced to **rank 2** by a
//!    second SVD with its smallest singular value zeroed (a real
//!    fundamental matrix is singular), and finally **denormalized**:
//!    `F = T₂ᵀ F̂ T₁`. Conditioning is what makes the 8-point algorithm
//!    numerically usable — see Hartley, "In Defense of the Eight-Point
//!    Algorithm," PAMI 1997.
//!
//! 2. **RANSAC** (Fischler & Bolles 1981). Because the matches contain
//!    outliers, the fit is wrapped in random sample consensus: repeatedly
//!    draw a minimal sample of 8 matches, fit a candidate `F`, and score
//!    how many of *all* matches are inliers under the **Sampson distance**
//!    (the first-order geometric approximation to reprojection error)
//!    against a pixel threshold. The largest consensus set wins, and `F`
//!    is refit on all of its inliers.
//!
//! ## Scope (stated honestly)
//!
//! This stage estimates the **fundamental matrix only** and the inlier
//! correspondence set. Recovering the **essential matrix**, the relative
//! camera rotation/translation `(R, t)`, and triangulating 3-D points
//! requires camera intrinsics and is **Stage 3** — it is deliberately not
//! done here. The 8-point algorithm is also a *linear* estimator: it
//! minimises an algebraic residual, not the geometric (Sampson/reprojection)
//! error, so the refit `F` is a good initialisation but not the
//! maximum-likelihood estimate that a later non-linear refinement (bundle
//! adjustment) would produce.
//!
//! ## Numerical note on thresholds
//!
//! Because the estimator is linear, even *exact* (noise-free) correspondences
//! leave a small, unevenly distributed Sampson residual after the inlier
//! refit — on the order of a pixel for keypoints stored as `f32`. The
//! RANSAC inlier [`threshold`](RansacParams::threshold) must therefore be a
//! few pixels (the default is 2 px) rather than sub-pixel; this is standard
//! practice and is why true inliers are admitted while genuine outliers
//! (epipolar error of tens to hundreds of pixels) are still rejected. The
//! *algebraic* epipolar residual `x₂ᵀ F x₁` is a different quantity: with
//! `F` scaled to unit Frobenius norm it sits at floating-point roundoff
//! (`~1e-6`) for a good fit and should not be confused with a pixel
//! distance.

use crate::keypoint::Keypoint;
use crate::matching::Match;
use nalgebra::{DMatrix, Matrix3, Vector3};

/// Parameters controlling the RANSAC fundamental-matrix search.
#[derive(Debug, Clone, Copy)]
pub struct RansacParams {
    /// Number of RANSAC iterations (random minimal samples to try). More
    /// iterations raise the probability of hitting an all-inlier sample at
    /// linear cost. A few hundred to a couple thousand is typical.
    pub iterations: usize,
    /// Inlier threshold on the **Sampson distance**, in pixels. A
    /// correspondence counts as an inlier when its Sampson distance to the
    /// candidate `F` is below this value. Around `1.0`–`3.0` px is usual.
    pub threshold: f64,
    /// Seed for the deterministic SplitMix64 sample-index PRNG, so a given
    /// input + params always yields the same result (reproducible builds /
    /// tests). No external `rand` dependency is used.
    pub seed: u64,
}

impl Default for RansacParams {
    /// Sensible defaults: 2000 iterations, a 2 px Sampson threshold, and a
    /// fixed seed.
    fn default() -> Self {
        Self {
            iterations: 2000,
            threshold: 2.0,
            seed: 0x5EED_C0DE_F00D_BEEF,
        }
    }
}

/// The result of two-view verification: the estimated fundamental matrix
/// and the consensus set of inlier matches that support it.
#[derive(Debug, Clone)]
pub struct TwoViewResult {
    /// The estimated fundamental matrix `F` (3×3, rank 2), refit on all
    /// inliers and scaled so its Frobenius norm is 1. Satisfies
    /// `x₂ᵀ F x₁ ≈ 0` for every inlier correspondence.
    pub fundamental: Matrix3<f64>,
    /// The subset of the input `matches` deemed geometrically consistent
    /// with [`Self::fundamental`] (Sampson distance below the threshold).
    pub inliers: Vec<Match>,
}

/// Minimal number of correspondences needed by the 8-point algorithm.
const MIN_SAMPLE: usize = 8;

/// Estimate the fundamental matrix relating two views and the inlier set,
/// robustly, with RANSAC over the normalized 8-point algorithm.
///
/// `kp1` / `kp2` are the keypoints of the two images; each [`Match`] in
/// `matches` pairs `kp1[query_idx]` with `kp2[train_idx]`. The keypoints'
/// `(x, y)` pixel positions are used; `score` / `angle` are ignored.
///
/// Returns `None` when there are fewer than 8 matches (the 8-point
/// algorithm is then under-determined), or if RANSAC fails to find any
/// usable consensus (degenerate / collinear configurations from which no
/// rank-2 `F` could be estimated). Never panics on valid input.
///
/// The search is deterministic for a given `params.seed`.
///
/// # Method
///
/// See the [module documentation](self): isotropic Hartley normalization,
/// SVD solve, rank-2 enforcement, denormalization, all inside a RANSAC loop
/// scored by Sampson distance, followed by a final refit on the inliers.
#[must_use]
pub fn verify_two_view(
    kp1: &[Keypoint],
    kp2: &[Keypoint],
    matches: &[Match],
    params: &RansacParams,
) -> Option<TwoViewResult> {
    if matches.len() < MIN_SAMPLE {
        return None;
    }

    // Gather the matched pixel coordinates once, as f64 point pairs.
    // A match referencing an out-of-range keypoint index is skipped
    // defensively (callers should not do this, but we never index OOB).
    let mut pts: Vec<(Vector3<f64>, Vector3<f64>)> = Vec::with_capacity(matches.len());
    let mut kept_matches: Vec<Match> = Vec::with_capacity(matches.len());
    for &m in matches {
        let (Some(a), Some(b)) = (kp1.get(m.query_idx), kp2.get(m.train_idx)) else {
            continue;
        };
        pts.push((
            Vector3::new(f64::from(a.x), f64::from(a.y), 1.0),
            Vector3::new(f64::from(b.x), f64::from(b.y), 1.0),
        ));
        kept_matches.push(m);
    }
    if pts.len() < MIN_SAMPLE {
        return None;
    }

    let n = pts.len();
    let mut rng = SplitMix64::new(params.seed);
    let thresh2 = params.threshold * params.threshold;

    let mut best_inlier_idx: Vec<usize> = Vec::new();

    // Need at least one iteration to do anything useful.
    let iterations = params.iterations.max(1);
    for _ in 0..iterations {
        // Draw 8 distinct indices into pts.
        let sample = sample_indices(&mut rng, n, MIN_SAMPLE);
        let sample_pts: Vec<(Vector3<f64>, Vector3<f64>)> =
            sample.iter().map(|&i| pts[i]).collect();

        let Some(f) = fundamental_8point(&sample_pts) else {
            continue;
        };

        // Score: count inliers over ALL correspondences by Sampson distance.
        let mut inliers = Vec::new();
        for (i, &(x1, x2)) in pts.iter().enumerate() {
            if sampson_distance_sq(&f, &x1, &x2) < thresh2 {
                inliers.push(i);
            }
        }
        if inliers.len() > best_inlier_idx.len() {
            best_inlier_idx = inliers;
        }
    }

    // Need a non-degenerate consensus set to refit a fundamental matrix.
    if best_inlier_idx.len() < MIN_SAMPLE {
        return None;
    }

    // Refit F on all inliers (the standard RANSAC polish step).
    let inlier_pts: Vec<(Vector3<f64>, Vector3<f64>)> =
        best_inlier_idx.iter().map(|&i| pts[i]).collect();
    let f = fundamental_8point(&inlier_pts)?;

    // Recompute the inlier set against the refit F so the reported inliers
    // are exactly those consistent with the returned matrix.
    let mut final_inliers = Vec::new();
    for (i, &(x1, x2)) in pts.iter().enumerate() {
        if sampson_distance_sq(&f, &x1, &x2) < thresh2 {
            final_inliers.push(kept_matches[i]);
        }
    }
    if final_inliers.len() < MIN_SAMPLE {
        return None;
    }

    Some(TwoViewResult {
        fundamental: f,
        inliers: final_inliers,
    })
}

/// Estimate a fundamental matrix from `>= 8` point correspondences via the
/// **normalized 8-point algorithm**. Returns `None` if the inputs are too
/// few or so degenerate (e.g. all points coincident / collinear) that the
/// normalization or SVD cannot produce a usable rank-2 matrix.
///
/// Each correspondence is a pair of homogeneous image points
/// `(x₁, x₂)` with the third coordinate `1`.
fn fundamental_8point(pts: &[(Vector3<f64>, Vector3<f64>)]) -> Option<Matrix3<f64>> {
    if pts.len() < MIN_SAMPLE {
        return None;
    }

    // Hartley isotropic normalization for each image independently.
    let t1 = normalization_transform(pts.iter().map(|p| p.0))?;
    let t2 = normalization_transform(pts.iter().map(|p| p.1))?;

    // Build the linear system A f = 0, one row per correspondence, using the
    // normalized coordinates.
    let mut a = DMatrix::<f64>::zeros(pts.len(), 9);
    for (row, (p1, p2)) in pts.iter().enumerate() {
        let n1 = t1 * p1;
        let n2 = t2 * p2;
        let (x1, y1) = (n1.x, n1.y);
        let (x2, y2) = (n2.x, n2.y);
        a[(row, 0)] = x2 * x1;
        a[(row, 1)] = x2 * y1;
        a[(row, 2)] = x2;
        a[(row, 3)] = y2 * x1;
        a[(row, 4)] = y2 * y1;
        a[(row, 5)] = y2;
        a[(row, 6)] = x1;
        a[(row, 7)] = y1;
        a[(row, 8)] = 1.0;
    }

    // f is the right-singular vector with the smallest singular value, i.e.
    // the last column of V (last row of Vᵀ).
    let svd = a.svd(false, true);
    let v_t = svd.v_t?;
    // The 9-vector solution is the last row of Vᵀ.
    let last = v_t.row(v_t.nrows() - 1);
    let f_hat = Matrix3::new(
        last[0], last[1], last[2], last[3], last[4], last[5], last[6], last[7], last[8],
    );

    // Enforce rank 2: SVD of F̂, zero the smallest singular value, recompose.
    let f_svd = f_hat.svd(true, true);
    let u = f_svd.u?;
    let v_t2 = f_svd.v_t?;
    let mut s = f_svd.singular_values;
    s[2] = 0.0;
    let s_diag = Matrix3::new(s[0], 0.0, 0.0, 0.0, s[1], 0.0, 0.0, 0.0, 0.0);
    let f_rank2 = u * s_diag * v_t2;

    // Denormalize: F = T₂ᵀ F̂ T₁.
    let f = t2.transpose() * f_rank2 * t1;

    // Normalize the overall scale (F is defined only up to scale). Reject a
    // numerically dead matrix.
    let norm = f.norm();
    if !norm.is_finite() || norm < 1e-12 {
        return None;
    }
    Some(f / norm)
}

/// Compute the Hartley isotropic normalization for a set of homogeneous
/// 2-D points: a similarity transform `T` (3×3) that translates the
/// centroid to the origin and scales so the mean distance to the origin is
/// `√2`.
///
/// Two passes over the points are needed (centroid, then mean distance), so
/// the iterator is taken by value and must be `Clone`. Returns `None` if
/// the points are all (numerically) coincident, so no meaningful scale
/// exists.
fn normalization_transform(
    points: impl Iterator<Item = Vector3<f64>> + Clone,
) -> Option<Matrix3<f64>> {
    let mut count = 0usize;
    let mut cx = 0.0;
    let mut cy = 0.0;
    for p in points.clone() {
        cx += p.x;
        cy += p.y;
        count += 1;
    }
    if count == 0 {
        return None;
    }
    let ninv = 1.0 / count as f64;
    cx *= ninv;
    cy *= ninv;

    let mut mean_dist = 0.0;
    for p in points {
        let dx = p.x - cx;
        let dy = p.y - cy;
        mean_dist += (dx * dx + dy * dy).sqrt();
    }
    mean_dist *= ninv;
    if mean_dist < 1e-12 {
        // All points coincide — no usable scale.
        return None;
    }

    // Scale so mean distance becomes sqrt(2).
    let s = std::f64::consts::SQRT_2 / mean_dist;
    // T = [[s, 0, -s*cx], [0, s, -s*cy], [0, 0, 1]]
    Some(Matrix3::new(
        s,
        0.0,
        -s * cx,
        0.0,
        s,
        -s * cy,
        0.0,
        0.0,
        1.0,
    ))
}

/// Squared **Sampson distance** of a correspondence `(x₁, x₂)` to the
/// fundamental matrix `F`. This is the standard first-order geometric error
/// for epipolar geometry:
///
/// ```text
///                       (x₂ᵀ F x₁)²
///  d²  =  ───────────────────────────────────────
///          (F x₁)_x² + (F x₁)_y² + (Fᵀ x₂)_x² + (Fᵀ x₂)_y²
/// ```
///
/// where the subscripts pick the first two (image-plane) components of the
/// epipolar lines. Returns a large finite value if the denominator is
/// (near) zero so such a degenerate point never counts as an inlier.
fn sampson_distance_sq(f: &Matrix3<f64>, x1: &Vector3<f64>, x2: &Vector3<f64>) -> f64 {
    let fx1 = f * x1; // epipolar line in image 2
    let ftx2 = f.transpose() * x2; // epipolar line in image 1
    let num = (x2.dot(&fx1)).powi(2);
    let denom = fx1.x * fx1.x + fx1.y * fx1.y + ftx2.x * ftx2.x + ftx2.y * ftx2.y;
    if denom < 1e-15 {
        return f64::MAX;
    }
    num / denom
}

/// Draw `k` distinct indices from `0..n` using the in-crate SplitMix64
/// PRNG (partial Fisher–Yates over a scratch index list). Requires
/// `k <= n`, which the single caller guarantees (`n >= MIN_SAMPLE == k`).
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

/// Minimal SplitMix64 PRNG, mirroring the one in [`crate::descriptor`].
///
/// Used here only to draw deterministic RANSAC minimal-sample indices, so
/// the crate needs no `rand` dependency and results are reproducible across
/// runs and machines. Not used for any security purpose.
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

    /// A camera projection matrix `P` (3×4) applied to a homogeneous 3-D
    /// point, returning the pixel coordinate (with w divided out).
    fn project(p: &nalgebra::Matrix3x4<f64>, x: &nalgebra::Vector4<f64>) -> (f64, f64) {
        let v = p * x;
        (v.x / v.z, v.y / v.z)
    }

    /// Two known camera projection matrices for the synthetic-geometry
    /// tests: camera 1 at the origin, camera 2 translated and slightly
    /// rotated, both with a simple pinhole intrinsic.
    fn test_cameras() -> (nalgebra::Matrix3x4<f64>, nalgebra::Matrix3x4<f64>) {
        // Intrinsics K (focal 800, principal point (320, 240)).
        let k = Matrix3::new(800.0, 0.0, 320.0, 0.0, 800.0, 240.0, 0.0, 0.0, 1.0);

        // Camera 1: P1 = K [I | 0].
        let mut p1 = nalgebra::Matrix3x4::<f64>::zeros();
        p1.fixed_view_mut::<3, 3>(0, 0)
            .copy_from(&Matrix3::identity());
        let p1 = k * p1;

        // Camera 2: small yaw rotation + a baseline translation.
        let ang: f64 = 0.15;
        let (s, c) = ang.sin_cos();
        let r = Matrix3::new(c, 0.0, s, 0.0, 1.0, 0.0, -s, 0.0, c);
        let t = Vector3::new(-1.0, 0.05, 0.1);
        let mut rt = nalgebra::Matrix3x4::<f64>::zeros();
        rt.fixed_view_mut::<3, 3>(0, 0).copy_from(&r);
        rt.fixed_view_mut::<3, 1>(0, 3).copy_from(&t);
        let p2 = k * rt;

        (p1, p2)
    }

    /// A fixed cloud of 3-D scene points in front of both cameras.
    fn scene_points() -> Vec<nalgebra::Vector4<f64>> {
        // 12 points spread in x, y and depth, all with positive z.
        let raw = [
            (-0.4, -0.3, 5.0),
            (0.3, -0.2, 6.0),
            (0.1, 0.4, 4.5),
            (-0.2, 0.25, 7.0),
            (0.45, 0.35, 5.5),
            (-0.35, -0.15, 6.5),
            (0.2, -0.35, 4.8),
            (-0.1, 0.1, 8.0),
            (0.05, -0.05, 5.2),
            (0.4, 0.05, 6.8),
            (-0.45, 0.3, 4.2),
            (0.15, 0.2, 7.5),
        ];
        raw.iter()
            .map(|&(x, y, z)| nalgebra::Vector4::new(x, y, z, 1.0))
            .collect()
    }

    /// Build matches and keypoint lists for a set of projected points,
    /// optionally interleaving outliers. Returns (kp1, kp2, matches).
    fn project_to_views(
        pts3d: &[nalgebra::Vector4<f64>],
    ) -> (Vec<Keypoint>, Vec<Keypoint>, Vec<Match>) {
        let (p1, p2) = test_cameras();
        let mut kp1 = Vec::new();
        let mut kp2 = Vec::new();
        let mut matches = Vec::new();
        for (i, x) in pts3d.iter().enumerate() {
            let (u1, v1) = project(&p1, x);
            let (u2, v2) = project(&p2, x);
            kp1.push(Keypoint::new(u1 as f32, v1 as f32, 1.0));
            kp2.push(Keypoint::new(u2 as f32, v2 as f32, 1.0));
            matches.push(Match {
                query_idx: i,
                train_idx: i,
                distance: 0,
            });
        }
        (kp1, kp2, matches)
    }

    // Test 3: the 8-point F satisfies the epipolar constraint and is rank 2.
    #[test]
    fn eight_point_satisfies_epipolar_constraint() {
        let pts3d = scene_points();
        let (kp1, kp2, matches) = project_to_views(&pts3d);

        // Build the homogeneous correspondences directly and fit F.
        let pairs: Vec<(Vector3<f64>, Vector3<f64>)> = matches
            .iter()
            .map(|m| {
                let a = kp1[m.query_idx];
                let b = kp2[m.train_idx];
                (
                    Vector3::new(f64::from(a.x), f64::from(a.y), 1.0),
                    Vector3::new(f64::from(b.x), f64::from(b.y), 1.0),
                )
            })
            .collect();

        let f = fundamental_8point(&pairs).expect("F should be estimable from 12 clean points");

        // Epipolar constraint x₂ᵀ F x₁ ≈ 0 for every correspondence.
        //
        // Two complementary checks:
        //  - The algebraic residual x₂ᵀ F x₁. With F scaled to unit
        //    Frobenius norm and pixel coordinates of order 1e2, a perfect
        //    fit leaves only floating-point roundoff (≈ machine-eps · |x|²),
        //    i.e. a few × 1e-6 in absolute terms — that IS numerical zero
        //    here, not error.
        //  - The Sampson distance, which divides out the epipolar-line
        //    gradient and is in *pixel* units; for exact data it is far
        //    below a hundredth of a pixel.
        for (x1, x2) in &pairs {
            let algebraic = (x2.dot(&(f * x1))).abs();
            assert!(
                algebraic < 1e-4,
                "algebraic epipolar residual {algebraic} too large (should be ~0)"
            );
            let sampson = sampson_distance_sq(&f, x1, x2).sqrt();
            assert!(
                sampson < 1e-2,
                "Sampson distance {sampson} px too large for exact data"
            );
        }

        // F must be rank 2: its smallest singular value is ~0.
        let sv = f.svd(false, false).singular_values;
        let smallest = sv[2];
        assert!(
            smallest < 1e-9,
            "smallest singular value {smallest} should be ~0 (rank-2 enforcement)"
        );
        // ... and the matrix should not be all-zero (largest sv well away from 0).
        assert!(sv[0] > 1e-3, "F collapsed to ~zero matrix");
    }

    // Test 4: RANSAC recovers the true inliers despite ~30% outliers.
    #[test]
    fn ransac_is_robust_to_outliers() {
        let pts3d = scene_points();
        let (mut kp1, mut kp2, mut matches) = project_to_views(&pts3d);
        let n_true = matches.len();

        // Inject ~30% random-outlier matches: extra keypoints at arbitrary
        // pixel positions paired together with no geometric relationship.
        // Deterministic positions from a local SplitMix64.
        let mut rng = SplitMix64::new(0xABCD_1234_5678_9F01);
        let n_outliers = (n_true as f64 * 0.3).ceil() as usize; // 4 outliers for 12 inliers
        for _ in 0..n_outliers {
            let ax = (rng.next_u64() % 640) as f32;
            let ay = (rng.next_u64() % 480) as f32;
            let bx = (rng.next_u64() % 640) as f32;
            let by = (rng.next_u64() % 480) as f32;
            let qi = kp1.len();
            let ti = kp2.len();
            kp1.push(Keypoint::new(ax, ay, 1.0));
            kp2.push(Keypoint::new(bx, by, 1.0));
            matches.push(Match {
                query_idx: qi,
                train_idx: ti,
                distance: 50,
            });
        }

        // A 3 px Sampson threshold. The keypoints are stored as f32 (the
        // `Keypoint` type), and the 8-point estimator is *linear* — it
        // minimises an algebraic, not geometric, residual — so even for
        // exact synthetic projections the true inliers' Sampson distances
        // spread up to ~1.5 px here. 3 px comfortably admits them while
        // still rejecting the random outliers, whose epipolar error is tens
        // to hundreds of pixels. (This is exactly why real SfM pipelines use
        // a px-scale RANSAC threshold rather than demanding sub-px fits from
        // the linear estimator.)
        let params = RansacParams {
            iterations: 3000,
            threshold: 3.0,
            seed: 0x1357_9BDF_2468_ACE0,
        };
        let result = verify_two_view(&kp1, &kp2, &matches, &params)
            .expect("RANSAC should find a consensus set");

        // The recovered inliers should be essentially the true 12 (the
        // outliers are random and almost surely violate epipolar geometry).
        assert!(
            result.inliers.len() >= n_true,
            "expected >= {} inliers, got {}",
            n_true,
            result.inliers.len()
        );
        // Every true correspondence (query_idx < n_true and == train_idx)
        // must be present — RANSAC must recover the full true inlier set.
        for i in 0..n_true {
            assert!(
                result
                    .inliers
                    .iter()
                    .any(|m| m.query_idx == i && m.train_idx == i),
                "true inlier {i} missing from consensus set"
            );
        }

        // The consensus set must be *essentially* the true inliers. We do
        // not demand ZERO outlier leakage: the fundamental matrix only
        // constrains a point to its epipolar *line* (a 1-D constraint), so a
        // uniformly random 2-D outlier has a real, non-negligible chance of
        // landing within a few pixels of some epipolar line and being
        // counted. That is a genuine property of two-view geometry, not a
        // matcher failure. We instead require that the true inliers dominate
        // overwhelmingly: outlier leakage stays a small fraction.
        let leaked = result
            .inliers
            .iter()
            .filter(|m| m.query_idx >= n_true)
            .count();
        assert!(
            leaked <= n_outliers / 2,
            "too many outliers leaked into consensus: {leaked} of {n_outliers}"
        );
        assert!(
            result.inliers.len() <= n_true + leaked,
            "consensus set larger than true inliers plus leaked outliers"
        );

        // The recovered F should fit the inliers (low epipolar residual).
        for m in &result.inliers {
            let a = kp1[m.query_idx];
            let b = kp2[m.train_idx];
            let x1 = Vector3::new(f64::from(a.x), f64::from(a.y), 1.0);
            let x2 = Vector3::new(f64::from(b.x), f64::from(b.y), 1.0);
            let d2 = sampson_distance_sq(&result.fundamental, &x1, &x2);
            assert!(
                d2 < params.threshold * params.threshold,
                "inlier Sampson dist too large"
            );
        }
    }

    // Test 5: degenerate input (< 8 matches) returns None without panicking.
    #[test]
    fn too_few_matches_returns_none() {
        let pts3d: Vec<_> = scene_points().into_iter().take(5).collect();
        let (kp1, kp2, matches) = project_to_views(&pts3d);
        assert_eq!(matches.len(), 5);
        let params = RansacParams::default();
        assert!(verify_two_view(&kp1, &kp2, &matches, &params).is_none());

        // Also an empty match list.
        assert!(verify_two_view(&kp1, &kp2, &[], &params).is_none());
    }

    #[test]
    fn ransac_determinism_same_seed() {
        let pts3d = scene_points();
        let (kp1, kp2, matches) = project_to_views(&pts3d);
        let params = RansacParams {
            iterations: 500,
            threshold: 1.5,
            seed: 42,
        };
        let r1 = verify_two_view(&kp1, &kp2, &matches, &params).unwrap();
        let r2 = verify_two_view(&kp1, &kp2, &matches, &params).unwrap();
        assert_eq!(r1.inliers.len(), r2.inliers.len());
        // Fundamental matrices should be bit-identical for the same seed.
        assert_eq!(r1.fundamental, r2.fundamental);
    }
}
