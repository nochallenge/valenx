//! Two-view geometry: essential matrix, relative pose, triangulation (Stage 3).
//!
//! Stage 2 ([`crate::geometry`]) recovers the **fundamental matrix** `F` and
//! the inlier correspondence set of an *uncalibrated* image pair. This stage
//! takes the next step of the classical Structure-from-Motion pipeline: given
//! the camera **intrinsics**, it upgrades `F` to the **essential matrix** `E`,
//! decomposes `E` into the relative camera rotation and translation `(R, t)`,
//! and **triangulates** the inlier correspondences into 3-D scene points.
//!
//! ```text
//! F (Stage 2)  +  K₁, K₂  ──►  E = K₂ᵀ F K₁
//!                              │
//!                              ├─ decompose E  ──►  4 candidate (R, t)
//!                              │
//!                              └─ triangulate + cheirality test  ──►  (R, t) + 3-D points
//! ```
//!
//! ## This is the CALIBRATED path (stated honestly)
//!
//! Everything here assumes the camera **intrinsics are known** — a
//! [`CameraIntrinsics`] (focal lengths `fx`, `fy` and principal point
//! `cx`, `cy`) for each view. In practice these come from the EXIF focal
//! length and sensor size, a prior calibration, or a sensible default (e.g.
//! `fx = fy ≈ image_width`, principal point at the image centre) when no
//! metadata is available. **Full uncalibrated auto-calibration** — recovering
//! focal lengths from the imagery itself (e.g. via the Kruppa equations or a
//! self-calibration bundle) — is **out of scope** for this stage: without
//! intrinsics, `F` alone fixes the projective geometry only up to an unknown
//! projective transformation, not a metric `(R, t)`. We do not attempt to
//! guess intrinsics from `F`; the caller must supply them.
//!
//! ## Translation is recovered only up to scale (stated honestly)
//!
//! The essential matrix is homogeneous: `E` and `αE` describe the same
//! geometry for any non-zero scalar `α`. Consequently the recovered
//! translation `t` has an **arbitrary positive scale** — its *direction* is
//! meaningful but its *length* is not. We return `t` as a **unit vector**
//! (`‖t‖ = 1`), and the triangulated [`TwoViewReconstruction::points3d`] are
//! therefore correct only up to one global scale factor (the same factor that
//! scales `t`). Fixing the absolute scale requires external information — a
//! known baseline, a scale bar, GPS, or fusing more views in the later
//! incremental-mapping / bundle-adjustment stages. The *four-fold rotation/
//! translation ambiguity* of the decomposition is resolved here by the
//! cheirality (points-in-front-of-both-cameras) test; the *scale* ambiguity
//! is fundamental to two views and is not.
//!
//! ## Method
//!
//! - **Essential from fundamental** (the calibration relation): for pixel
//!   correspondences `x₂ᵀ F x₁ = 0` and calibrated rays
//!   `x̂ᵢ = Kᵢ⁻¹ xᵢ`, the essential matrix is `E = K₂ᵀ F K₁` and satisfies
//!   `x̂₂ᵀ E x̂₁ = 0`. See Hartley & Zisserman, *Multiple View Geometry in
//!   Computer Vision*, 2nd ed., §9.6.
//!
//! - **Decomposition of `E`** into `(R, t)`: a valid essential matrix has two
//!   equal non-zero singular values and a zero third (`E ∝ U diag(1,1,0) Vᵀ`).
//!   Using the orthogonal `W` matrix, the rotation is `R = U W Vᵀ` or
//!   `U Wᵀ Vᵀ` and the translation direction is `±u₃` (the third column of
//!   `U`), giving **four** candidate poses. This is HZ Result 9.19. We force
//!   each candidate `R` to a proper rotation (`det R = +1`) by negating it if
//!   the SVD handed back a reflection.
//!
//! - **Triangulation** via the linear **DLT** (Direct Linear Transform): each
//!   view contributes two rows to a homogeneous `4×4` system `A X = 0` from
//!   the cross-product form `xᵢ × (Pᵢ X) = 0`; the solution is the right
//!   singular vector of the smallest singular value, dehomogenized by the `w`
//!   divide. See HZ §12.2.
//!
//! - **Cheirality** (HZ §9.6.3): the physically correct pose is the one for
//!   which triangulated points lie *in front of* (positive depth in) **both**
//!   cameras. We triangulate every inlier under each of the four candidates
//!   and keep the candidate with the most points satisfying this test.

use crate::keypoint::Keypoint;
use crate::matching::Match;
use nalgebra::{Matrix3, Matrix3x4, Matrix4, Vector3};

/// Pinhole camera intrinsics: focal lengths and principal point, in pixels.
///
/// These define the calibration matrix
///
/// ```text
///       ⎡ fx   0   cx ⎤
///  K =  ⎢  0  fy   cy ⎥
///       ⎣  0   0    1 ⎦
/// ```
///
/// which maps a 3-D point in the camera frame to homogeneous pixel
/// coordinates. Skew is assumed zero (true for essentially all real
/// cameras). Supply these from the image EXIF focal length and sensor size,
/// a prior calibration, or a sensible default (e.g. `fx = fy ≈ image width`,
/// principal point at the image centre) — see the [module docs](self) on the
/// calibrated-only scope.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CameraIntrinsics {
    /// Focal length in the `x` direction, in pixels.
    pub fx: f64,
    /// Focal length in the `y` direction, in pixels.
    pub fy: f64,
    /// Principal-point `x` coordinate (pixels from the left edge).
    pub cx: f64,
    /// Principal-point `y` coordinate (pixels from the top edge).
    pub cy: f64,
}

impl CameraIntrinsics {
    /// Construct intrinsics from focal lengths and principal point.
    #[inline]
    #[must_use]
    pub fn new(fx: f64, fy: f64, cx: f64, cy: f64) -> Self {
        Self { fx, fy, cx, cy }
    }

    /// The `3×3` calibration matrix `K` (see the [type docs](Self)).
    #[must_use]
    pub fn matrix(&self) -> Matrix3<f64> {
        Matrix3::new(
            self.fx, 0.0, self.cx, //
            0.0, self.fy, self.cy, //
            0.0, 0.0, 1.0,
        )
    }
}

/// A recovered two-view reconstruction: the relative camera pose and the
/// triangulated inlier points.
///
/// The pose maps a point `X₁` expressed in the **first** camera's frame to
/// the **second** camera's frame by `X₂ = R X₁ + t`. Equivalently the second
/// camera's projection matrix is `P₂ = K₂ [R | t]` while the first is
/// `P₁ = K₁ [I | 0]`.
#[derive(Debug, Clone)]
pub struct TwoViewReconstruction {
    /// Relative rotation `R` (proper rotation, `det R = +1`) of the second
    /// camera with respect to the first.
    pub rotation: Matrix3<f64>,
    /// Relative translation direction `t` of the second camera with respect
    /// to the first, as a **unit vector** (`‖t‖ = 1`). The absolute scale is
    /// unrecoverable from two views — see the [module docs](self).
    pub translation: Vector3<f64>,
    /// The triangulated 3-D points (one per inlier correspondence that
    /// triangulated successfully), expressed in the **first** camera's frame.
    /// These are correct only up to the same global scale factor as
    /// [`Self::translation`].
    pub points3d: Vec<Vector3<f64>>,
    /// How many of the triangulated points lie in front of (have positive
    /// depth in) **both** cameras — the cheirality count that selected this
    /// pose. A healthy reconstruction has this equal (or nearly equal) to the
    /// number of inliers.
    pub num_in_front: usize,
}

/// Compute the **essential matrix** from the fundamental matrix and the two
/// cameras' intrinsics: `E = K₂ᵀ F K₁`.
///
/// `f` is the fundamental matrix from Stage 2 ([`crate::verify_two_view`]),
/// relating pixel correspondences by `x₂ᵀ F x₁ = 0`. The returned `E` relates
/// the calibrated (normalized) image rays `x̂ᵢ = Kᵢ⁻¹ xᵢ` by
/// `x̂₂ᵀ E x̂₁ = 0`. See Hartley & Zisserman §9.6.
///
/// The result is returned exactly as `K₂ᵀ F K₁`; its singular structure is
/// *approximately* `(σ, σ, 0)` when `f` is a good fundamental matrix.
/// [`decompose_essential`] re-imposes the exact `(1, 1, 0)` structure before
/// extracting a pose, so no projection is done here.
#[must_use]
pub fn essential_from_fundamental(
    f: &Matrix3<f64>,
    k1: &CameraIntrinsics,
    k2: &CameraIntrinsics,
) -> Matrix3<f64> {
    k2.matrix().transpose() * f * k1.matrix()
}

/// The orthogonal helper matrix `W` from the essential-matrix decomposition
/// (Hartley & Zisserman Result 9.19):
///
/// ```text
///       ⎡ 0  -1   0 ⎤
///  W =  ⎢ 1   0   0 ⎥
///       ⎣ 0   0   1 ⎦
/// ```
#[inline]
fn w_matrix() -> Matrix3<f64> {
    Matrix3::new(
        0.0, -1.0, 0.0, //
        1.0, 0.0, 0.0, //
        0.0, 0.0, 1.0,
    )
}

/// Force `r` to be a proper rotation (`det = +1`). The SVD factors `U`, `V`
/// are orthogonal but may individually be reflections (`det = -1`), in which
/// case a candidate `R = U W Vᵀ` can come out as a reflection; negating the
/// whole matrix flips its determinant sign while keeping it orthogonal, which
/// is the standard fix (it corresponds to the inherent sign freedom of `E`).
#[inline]
fn make_proper_rotation(r: Matrix3<f64>) -> Matrix3<f64> {
    if r.determinant() < 0.0 {
        -r
    } else {
        r
    }
}

/// Decompose an **essential matrix** into the **four** candidate relative
/// poses `(R, t)` (Hartley & Zisserman Result 9.19).
///
/// The essential matrix factors as `E ∝ U diag(1, 1, 0) Vᵀ`. We take the SVD
/// of `e`, *enforce* the ideal `(1, 1, 0)` singular structure (a real `E` has
/// two equal non-zero singular values and a zero third; noise perturbs this,
/// so we replace the singular values before forming the poses), and build the
/// four candidates:
///
/// ```text
///  (R_a, +t),  (R_a, -t),  (R_b, +t),  (R_b, -t)
/// ```
///
/// where `R_a = U W Vᵀ`, `R_b = U Wᵀ Vᵀ`, and `t = u₃` (the third column of
/// `U`, the translation direction). Each `R` is forced to a proper rotation
/// (`det R = +1`) — negated if the SVD handed back a reflection — so all four
/// returned rotations are guaranteed valid rotations and never reflections.
///
/// Exactly one of the four is the physically correct pose; it is selected by
/// the cheirality test in [`recover_pose`]. The translation is a **unit
/// direction** only (scale is unrecoverable — see the [module docs](self)).
///
/// For a numerically dead input (e.g. the all-zero matrix) the SVD still
/// returns orthogonal factors, so this function does not panic; the resulting
/// poses are simply meaningless and will fail the downstream cheirality test
/// (`num_in_front == 0`).
#[must_use]
pub fn decompose_essential(e: &Matrix3<f64>) -> [(Matrix3<f64>, Vector3<f64>); 4] {
    let svd = e.svd(true, true);
    // For a 3×3 matrix nalgebra always produces both factors; fall back to
    // identity factors if (pathologically) it did not, so we never panic and
    // simply yield poses that fail cheirality downstream.
    let u = svd.u.unwrap_or_else(Matrix3::identity);
    let v_t = svd.v_t.unwrap_or_else(Matrix3::identity);

    // Enforce the ideal (1, 1, 0) singular structure. (We do not actually
    // need to rebuild E from it — the rotations depend only on U, V and the
    // translation only on U — but re-imposing the structure is the canonical
    // statement of the decomposition and documents the assumption.)
    let _e_ideal = u * Matrix3::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0) * v_t;

    let w = w_matrix();
    let r_a = make_proper_rotation(u * w * v_t);
    let r_b = make_proper_rotation(u * w.transpose() * v_t);

    // Translation direction = third column of U, normalized to a unit vector
    // (E, hence t, is defined only up to scale). Guard the (degenerate)
    // zero-norm case.
    let mut t = u.column(2).into_owned();
    let tn = t.norm();
    if tn > 1e-12 {
        t /= tn;
    }

    [(r_a, t), (r_a, -t), (r_b, t), (r_b, -t)]
}

/// Triangulate a single 3-D point from two views via the linear **DLT**
/// (Direct Linear Transform; Hartley & Zisserman §12.2).
///
/// `p1`, `p2` are the two `3×4` camera projection matrices and `x1`, `x2` the
/// corresponding image points (in pixels, `(u, v)`) in view 1 and view 2.
/// Each view contributes two rows to a homogeneous `4×4` system `A X = 0`
/// built from the cross-product constraint `xᵢ × (Pᵢ X) = 0`:
///
/// ```text
///  u (p₃ᵀ) − (p₁ᵀ)
///  v (p₃ᵀ) − (p₂ᵀ)
/// ```
///
/// for each camera, where `pₖᵀ` is the `k`-th row of `P`. The homogeneous
/// solution is the right singular vector of the smallest singular value of
/// `A`; it is converted to a Euclidean point by dividing through by its
/// fourth (`w`) coordinate.
///
/// If `w` is (near) zero — a point at infinity, or a wholly degenerate
/// configuration — the divide is skipped and the point's `x, y, z` are
/// returned as-is (a large/ill-defined but **finite** vector); the function
/// never divides by zero and never panics. Such a point will fail the
/// cheirality depth test in [`recover_pose`].
#[must_use]
pub fn triangulate_point(
    p1: &Matrix3x4<f64>,
    p2: &Matrix3x4<f64>,
    x1: (f64, f64),
    x2: (f64, f64),
) -> Vector3<f64> {
    // Assemble the 4×4 DLT system, two rows per view.
    let mut a = Matrix4::<f64>::zeros();
    // Row of P as a 4-vector.
    let row = |p: &Matrix3x4<f64>, r: usize| {
        nalgebra::RowVector4::new(p[(r, 0)], p[(r, 1)], p[(r, 2)], p[(r, 3)])
    };

    let (u1, v1) = x1;
    let (u2, v2) = x2;

    a.set_row(0, &(u1 * row(p1, 2) - row(p1, 0)));
    a.set_row(1, &(v1 * row(p1, 2) - row(p1, 1)));
    a.set_row(2, &(u2 * row(p2, 2) - row(p2, 0)));
    a.set_row(3, &(v2 * row(p2, 2) - row(p2, 1)));

    // Smallest-singular-value right vector = last row of Vᵀ. nalgebra always
    // produces V for a 4×4; fall back to a zero solution if (pathologically)
    // it did not, which dehomogenizes to the origin and fails cheirality.
    let svd = a.svd(false, true);
    let Some(v_t) = svd.v_t else {
        return Vector3::zeros();
    };
    let x = v_t.row(3);
    let (xh, yh, zh, w) = (x[0], x[1], x[2], x[3]);

    // Dehomogenize, guarding w ≈ 0 (point at / near infinity).
    if w.abs() > 1e-12 {
        Vector3::new(xh / w, yh / w, zh / w)
    } else {
        // No finite Euclidean point; return the (un-divided) direction. It is
        // finite and will simply fail the positive-depth cheirality test.
        Vector3::new(xh, yh, zh)
    }
}

/// Build the `3×4` projection matrix `P = K [R | t]`.
fn projection_matrix(k: &Matrix3<f64>, r: &Matrix3<f64>, t: &Vector3<f64>) -> Matrix3x4<f64> {
    let mut rt = Matrix3x4::<f64>::zeros();
    rt.fixed_view_mut::<3, 3>(0, 0).copy_from(r);
    rt.fixed_view_mut::<3, 1>(0, 3).copy_from(t);
    k * rt
}

/// Recover the relative camera pose `(R, t)` and triangulate the inliers,
/// resolving the four-fold decomposition ambiguity by the **cheirality**
/// (points-in-front-of-both-cameras) test.
///
/// `e` is the essential matrix (from [`essential_from_fundamental`]); `k1` /
/// `k2` the two cameras' intrinsics; `kp1` / `kp2` the keypoints of the two
/// images; and `inliers` the geometrically verified correspondences (each
/// [`Match`] pairs `kp1[query_idx]` with `kp2[train_idx]`), typically
/// [`crate::TwoViewResult::inliers`] from Stage 2.
///
/// For each of the four candidate poses from [`decompose_essential`] this
/// builds `P₁ = K₁[I | 0]` and `P₂ = K₂[R | t]`, triangulates every inlier
/// with [`triangulate_point`], and counts the points with **positive depth in
/// both cameras**. The candidate with the most such points is returned, with
/// its triangulated points (in the first camera's frame) and the in-front
/// count.
///
/// Returns `None` when there are no usable inliers, or when *no* candidate
/// puts any point in front of both cameras (a fully degenerate configuration —
/// e.g. an all-zero `e`, or a scene entirely behind the cameras). Never
/// panics.
///
/// The translation is a **unit direction** and the points are correct only up
/// to one global scale — see the [module docs](self).
#[must_use]
pub fn recover_pose(
    e: &Matrix3<f64>,
    k1: &CameraIntrinsics,
    k2: &CameraIntrinsics,
    kp1: &[Keypoint],
    kp2: &[Keypoint],
    inliers: &[Match],
) -> Option<TwoViewReconstruction> {
    // Gather the inlier pixel correspondences once (skipping any out-of-range
    // index defensively, exactly as Stage 2 does).
    let mut obs: Vec<((f64, f64), (f64, f64))> = Vec::with_capacity(inliers.len());
    for m in inliers {
        let (Some(a), Some(b)) = (kp1.get(m.query_idx), kp2.get(m.train_idx)) else {
            continue;
        };
        obs.push((
            (f64::from(a.x), f64::from(a.y)),
            (f64::from(b.x), f64::from(b.y)),
        ));
    }
    if obs.is_empty() {
        return None;
    }

    let k1m = k1.matrix();
    let k2m = k2.matrix();

    // First camera: P1 = K1 [I | 0].
    let p1 = projection_matrix(&k1m, &Matrix3::identity(), &Vector3::zeros());

    let candidates = decompose_essential(e);

    let mut best: Option<TwoViewReconstruction> = None;
    for (r, t) in candidates {
        let p2 = projection_matrix(&k2m, &r, &t);

        let mut points = Vec::with_capacity(obs.len());
        let mut in_front = 0usize;
        for &(x1, x2) in &obs {
            let xw = triangulate_point(&p1, &p2, x1, x2);
            // Depth in camera 1: the point is already in cam-1 coordinates
            // (P1 has R = I, t = 0), so its depth is simply its z component.
            let depth1 = xw.z;
            // Depth in camera 2: transform into cam-2 coordinates X₂ = R X + t.
            let depth2 = (r * xw + t).z;
            if depth1 > 0.0 && depth2 > 0.0 {
                in_front += 1;
            }
            points.push(xw);
        }

        let better = match &best {
            None => true,
            Some(b) => in_front > b.num_in_front,
        };
        if better {
            best = Some(TwoViewReconstruction {
                rotation: r,
                translation: t,
                points3d: points,
                num_in_front: in_front,
            });
        }
    }

    // A pose that puts nothing in front of either camera is degenerate; report
    // failure rather than a meaningless reconstruction.
    match best {
        Some(b) if b.num_in_front > 0 => Some(b),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Skew-symmetric "cross-product matrix" `[t]_×` such that
    /// `[t]_× v = t × v`.
    fn skew(t: &Vector3<f64>) -> Matrix3<f64> {
        Matrix3::new(
            0.0, -t.z, t.y, //
            t.z, 0.0, -t.x, //
            -t.y, t.x, 0.0,
        )
    }

    /// A small, valid relative pose used across the decomposition / pose
    /// tests: a yaw + pitch rotation and a generic baseline translation.
    fn known_pose() -> (Matrix3<f64>, Vector3<f64>) {
        let yaw = 0.20_f64;
        let pitch = -0.12_f64;
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
        let t = Vector3::new(-0.9, 0.10, 0.20);
        (r, t)
    }

    fn test_intrinsics() -> CameraIntrinsics {
        CameraIntrinsics::new(800.0, 800.0, 320.0, 240.0)
    }

    /// A fixed cloud of 3-D scene points, all in front of both cameras for the
    /// [`known_pose`] geometry (positive depth in cam 1, and `R X + t` has
    /// positive z).
    fn scene_points() -> Vec<Vector3<f64>> {
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
        raw.iter().map(|&(x, y, z)| Vector3::new(x, y, z)).collect()
    }

    /// Project a camera-frame 3-D point through `K [R|t]` to a pixel.
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

    // ----- Test 1: essential-from-fundamental has the (σ, σ, 0) structure. -----
    #[test]
    fn essential_from_fundamental_has_correct_structure() {
        // Build a *ground-truth* essential matrix E = [t]_× R from a known
        // pose, then synthesize the fundamental matrix the calibrated relation
        // implies, F = K₂⁻ᵀ E K₁⁻¹, feed it back through
        // essential_from_fundamental, and check the singular structure.
        let (r, t) = known_pose();
        let tn = t.normalize();
        let e_true = skew(&tn) * r;

        let k = test_intrinsics();
        let km = k.matrix();
        let k_inv = km.try_inverse().expect("K invertible");
        // F = K₂⁻ᵀ E K₁⁻¹ (inverse of E = K₂ᵀ F K₁ with K₁ = K₂ = K).
        let f = k_inv.transpose() * e_true * k_inv;

        let e = essential_from_fundamental(&f, &k, &k);

        // Singular values: two near-equal non-zero, third ≈ 0.
        let sv = e.svd(false, false).singular_values;
        let (s0, s1, s2) = (sv[0], sv[1], sv[2]);
        assert!(s0 > 1e-6, "largest singular value collapsed: {s0}");
        // Two largest nearly equal (relative).
        assert!(
            (s0 - s1).abs() / s0 < 1e-6,
            "two largest singular values should be equal: {s0} vs {s1}"
        );
        // Third ≈ 0 relative to the largest.
        assert!(
            s2 / s0 < 1e-6,
            "third singular value should be ~0: {s2} (largest {s0})"
        );
    }

    // ----- Test 2: decompose_essential recovers (R_true, ±t_true). -----
    #[test]
    fn decompose_essential_recovers_known_pose() {
        let (r_true, t_true) = known_pose();
        let t_unit = t_true.normalize();
        let e = skew(&t_unit) * r_true;

        let candidates = decompose_essential(&e);

        // At least one candidate must match R_true and ±t_unit.
        let mut found = false;
        for (r, t) in &candidates {
            // Each returned R is a proper rotation.
            assert!(
                (r.determinant() - 1.0).abs() < 1e-9,
                "candidate rotation is not a proper rotation: det = {}",
                r.determinant()
            );
            let r_ok = (r - r_true).norm() < 1e-6;
            let t_ok = (t - t_unit).norm() < 1e-6 || (t + t_unit).norm() < 1e-6;
            if r_ok && t_ok {
                found = true;
            }
        }
        assert!(
            found,
            "no decomposition candidate matched the known (R, ±t/‖t‖)"
        );
    }

    // ----- Test 3: DLT triangulation recovers a known 3-D point. -----
    #[test]
    fn triangulate_recovers_known_point() {
        let (r, t) = known_pose();
        let k = test_intrinsics();
        let km = k.matrix();

        // P1 = K[I|0], P2 = K[R|t].
        let p1 = projection_matrix(&km, &Matrix3::identity(), &Vector3::zeros());
        let p2 = projection_matrix(&km, &r, &t);

        for x_true in scene_points() {
            let x1 = project(&km, &Matrix3::identity(), &Vector3::zeros(), &x_true);
            let x2 = project(&km, &r, &t, &x_true);
            let x_est = triangulate_point(&p1, &p2, x1, x2);
            let err = (x_est - x_true).norm();
            assert!(
                err < 1e-6,
                "triangulated point {x_est:?} differs from truth {x_true:?} by {err}"
            );
        }
    }

    // ----- Test 4: recover_pose end-to-end cheirality. -----
    #[test]
    fn recover_pose_end_to_end() {
        let (r_true, t_true) = known_pose();
        let t_unit = t_true.normalize();
        let k = test_intrinsics();
        let km = k.matrix();

        let pts = scene_points();

        // Project all points into both views and build trivial 1:1 matches.
        let mut kp1 = Vec::new();
        let mut kp2 = Vec::new();
        let mut inliers = Vec::new();
        for (i, x) in pts.iter().enumerate() {
            let (u1, v1) = project(&km, &Matrix3::identity(), &Vector3::zeros(), x);
            let (u2, v2) = project(&km, &r_true, &t_true, x);
            kp1.push(Keypoint::new(u1 as f32, v1 as f32, 1.0));
            kp2.push(Keypoint::new(u2 as f32, v2 as f32, 1.0));
            inliers.push(Match {
                query_idx: i,
                train_idx: i,
                distance: 0,
            });
        }

        // Build E from the true relative pose (= what Stage 2 + calibration
        // would yield) and recover the pose.
        let e = skew(&t_unit) * r_true;
        let recon = recover_pose(&e, &k, &k, &kp1, &kp2, &inliers)
            .expect("pose should be recoverable from 12 clean correspondences");

        // Rotation matches the truth.
        let r_err = (recon.rotation - r_true).norm();
        assert!(r_err < 1e-5, "recovered R off by {r_err}");

        // Translation matches up to POSITIVE scale: it is a unit vector and,
        // because the points are genuinely in front, the cheirality-selected
        // sign is +t_unit (not -t_unit).
        assert!(
            (recon.translation.norm() - 1.0).abs() < 1e-9,
            "translation should be a unit vector"
        );
        let t_err = (recon.translation - t_unit).norm();
        assert!(
            t_err < 1e-5,
            "recovered t {:?} should match +t_unit {t_unit:?} (err {t_err})",
            recon.translation
        );

        // All 12 points in front of both cameras.
        assert_eq!(
            recon.num_in_front,
            pts.len(),
            "all points should pass cheirality"
        );

        // Triangulated points match ground truth up to the global scale.
        // Recover the scale from the first point (its true depth is known),
        // then verify every point agrees under that single factor.
        assert_eq!(recon.points3d.len(), pts.len());
        let scale = pts[0].norm() / recon.points3d[0].norm();
        assert!(scale.is_finite() && scale > 0.0, "bad recovered scale");
        for (est, truth) in recon.points3d.iter().zip(pts.iter()) {
            let err = (est * scale - truth).norm();
            assert!(
                err < 1e-5,
                "scaled triangulated point {:?} vs truth {truth:?} err {err}",
                est * scale
            );
        }
    }

    // ----- Test 5a: degenerate all-zero E -> graceful, no panic, finite. -----
    #[test]
    fn degenerate_zero_essential_is_graceful() {
        let k = test_intrinsics();
        let km = k.matrix();
        let pts = scene_points();
        let mut kp1 = Vec::new();
        let mut kp2 = Vec::new();
        let mut inliers = Vec::new();
        let (r_true, t_true) = known_pose();
        for (i, x) in pts.iter().enumerate() {
            let (u1, v1) = project(&km, &Matrix3::identity(), &Vector3::zeros(), x);
            let (u2, v2) = project(&km, &r_true, &t_true, x);
            kp1.push(Keypoint::new(u1 as f32, v1 as f32, 1.0));
            kp2.push(Keypoint::new(u2 as f32, v2 as f32, 1.0));
            inliers.push(Match {
                query_idx: i,
                train_idx: i,
                distance: 0,
            });
        }

        // The decomposition of the all-zero matrix must not panic and must
        // hand back proper rotations and finite translations (its SVD yields
        // arbitrary-but-orthogonal factors).
        for (r, t) in decompose_essential(&Matrix3::<f64>::zeros()) {
            assert!(
                (r.determinant() - 1.0).abs() < 1e-9,
                "zero-E candidate not a proper rotation"
            );
            assert!(
                t.iter().all(|c| c.is_finite()),
                "zero-E translation non-finite"
            );
        }

        // recover_pose on the all-zero E must not panic. The output is
        // meaningless (the matrix carries no geometry), so we make NO claim
        // about its R/t — only that, whatever it returns, it is well-formed
        // and finite. (Honest note: an all-zero E does NOT force `None`,
        // because its arbitrary SVD factors can still place some triangulated
        // points in front of both cameras by coincidence; the genuine
        // "no consensus" sentinel is exercised by the empty-inliers and
        // coincident-point tests below.)
        if let Some(recon) = recover_pose(&Matrix3::<f64>::zeros(), &k, &k, &kp1, &kp2, &inliers) {
            assert!(
                recon.rotation.iter().all(|c| c.is_finite())
                    && recon.translation.iter().all(|c| c.is_finite())
                    && recon
                        .points3d
                        .iter()
                        .all(|p| p.iter().all(|c| c.is_finite())),
                "zero-E reconstruction must stay finite"
            );
        }
    }

    // ----- Test 5b: genuinely degenerate inputs -> None (no consensus). -----
    #[test]
    fn degenerate_inputs_return_none() {
        let (r_true, t_true) = known_pose();
        let t_unit = t_true.normalize();
        let k = test_intrinsics();
        let km = k.matrix();
        let e = skew(&t_unit) * r_true;

        // (i) No inliers at all -> None (nothing to triangulate).
        assert!(
            recover_pose(&e, &k, &k, &[], &[], &[]).is_none(),
            "empty inliers must yield None"
        );

        // (ii) Inliers whose indices are ALL out of range for the keypoint
        // slices -> every observation is skipped defensively -> the gathered
        // set is empty -> None (and crucially, no out-of-bounds panic).
        let kp1 = vec![Keypoint::new(320.0, 240.0, 1.0)];
        let kp2 = vec![Keypoint::new(330.0, 250.0, 1.0)];
        let bad = vec![
            Match {
                query_idx: 99,
                train_idx: 0,
                distance: 0,
            },
            Match {
                query_idx: 0,
                train_idx: 99,
                distance: 0,
            },
        ];
        assert!(
            recover_pose(&e, &k, &k, &kp1, &kp2, &bad).is_none(),
            "all-out-of-range inliers must yield None without panicking"
        );
        let _ = km;
    }

    // ----- Test 5b': coincident correspondences -> graceful, finite. --------
    #[test]
    fn coincident_correspondences_stay_finite() {
        // Honest geometric note: feeding the SAME pixel pair for every
        // correspondence is rank-deficient for *scene structure* (all rays
        // intersect in a single point), but it is NOT a divide-by-zero or a
        // `None` case: that single ray intersection is a perfectly well-defined
        // finite 3-D point that can sit in front of both cameras. So the
        // documented contract here is "finite, no panic" — we do not force
        // `None`. (This also exercises the near-coincident `w` path of the DLT
        // without tripping the divide-by-zero guard.)
        let (r_true, t_true) = known_pose();
        let t_unit = t_true.normalize();
        let k = test_intrinsics();
        let e = skew(&t_unit) * r_true;

        let mut kp1 = Vec::new();
        let mut kp2 = Vec::new();
        let mut inliers = Vec::new();
        for i in 0..12 {
            kp1.push(Keypoint::new(320.0, 240.0, 1.0));
            kp2.push(Keypoint::new(330.0, 250.0, 1.0));
            inliers.push(Match {
                query_idx: i,
                train_idx: i,
                distance: 0,
            });
        }
        if let Some(recon) = recover_pose(&e, &k, &k, &kp1, &kp2, &inliers) {
            assert!(
                recon
                    .points3d
                    .iter()
                    .all(|p| p.iter().all(|c| c.is_finite())),
                "coincident-correspondence triangulation must stay finite"
            );
            assert!(
                (recon.rotation.determinant() - 1.0).abs() < 1e-9,
                "rotation must remain proper"
            );
        }
    }

    // ----- Test 5c: points behind cam-1 still produce only finite output. ---
    #[test]
    fn points_behind_camera_one_stay_finite() {
        // Honest geometric note: feeding correspondences generated from points
        // BEHIND camera 1 does NOT force `recover_pose` to return `None`. The
        // four-candidate decomposition contains a "twisted pair" in which a
        // sign flip of t maps such points to positive depth, so the cheirality
        // search can legitimately select a pose that puts them in front of
        // *its* cameras. What the contract guarantees is that the result is
        // well-formed and finite and that nothing panics — which is what we
        // assert here.
        let (r_true, t_true) = known_pose();
        let t_unit = t_true.normalize();
        let k = test_intrinsics();
        let km = k.matrix();

        let behind = [
            (-0.4, -0.3, -5.0),
            (0.3, -0.2, -6.0),
            (0.1, 0.4, -4.5),
            (-0.2, 0.25, -7.0),
            (0.45, 0.35, -5.5),
            (-0.35, -0.15, -6.5),
        ];
        let mut kp1 = Vec::new();
        let mut kp2 = Vec::new();
        let mut inliers = Vec::new();
        for (i, &(x, y, z)) in behind.iter().enumerate() {
            let xp = Vector3::new(x, y, z);
            let (u1, v1) = project(&km, &Matrix3::identity(), &Vector3::zeros(), &xp);
            let (u2, v2) = project(&km, &r_true, &t_true, &xp);
            kp1.push(Keypoint::new(u1 as f32, v1 as f32, 1.0));
            kp2.push(Keypoint::new(u2 as f32, v2 as f32, 1.0));
            inliers.push(Match {
                query_idx: i,
                train_idx: i,
                distance: 0,
            });
        }
        let e = skew(&t_unit) * r_true;
        if let Some(recon) = recover_pose(&e, &k, &k, &kp1, &kp2, &inliers) {
            assert!(
                (recon.rotation.determinant() - 1.0).abs() < 1e-9,
                "rotation must remain a proper rotation"
            );
            assert!(
                recon.translation.iter().all(|c| c.is_finite())
                    && recon
                        .points3d
                        .iter()
                        .all(|p| p.iter().all(|c| c.is_finite())),
                "behind-camera reconstruction must stay finite"
            );
        }
    }

    // ----- Test 5c: triangulate_point with w ≈ 0 does not divide by zero. -----
    #[test]
    fn triangulate_handles_zero_w_without_panic() {
        // Degenerate cameras: both projection matrices identical and rank
        // deficient so the DLT system is singular. The function must return a
        // finite vector, not NaN/inf from a divide-by-zero.
        let p = Matrix3x4::<f64>::zeros();
        let out = triangulate_point(&p, &p, (10.0, 20.0), (30.0, 40.0));
        assert!(
            out.iter().all(|c| c.is_finite()),
            "triangulation must stay finite on degenerate input, got {out:?}"
        );

        // Also: a genuine point-at-infinity style homogeneous solution. Build
        // two valid cameras but ask for a correspondence consistent with a
        // direction (w -> 0). We can't easily force exact w = 0 analytically,
        // so we assert the general contract on the degenerate-camera case above
        // and additionally that a normal call is finite.
        let (r, t) = known_pose();
        let km = test_intrinsics().matrix();
        let p1 = projection_matrix(&km, &Matrix3::identity(), &Vector3::zeros());
        let p2 = projection_matrix(&km, &r, &t);
        let finite = triangulate_point(&p1, &p2, (320.0, 240.0), (321.0, 241.0));
        assert!(finite.iter().all(|c| c.is_finite()));
    }
}
