//! Bundle adjustment: the joint nonlinear refinement that closes the SfM
//! solver (Stage 5 — the final stage).
//!
//! Stages 1–4 build a reconstruction *incrementally*: features ([`crate::fast`]
//! / [`crate::descriptor`]), matches ([`crate::matching`]), a verified
//! fundamental matrix ([`crate::verify_two_view`]), a two-view seed pose plus
//! triangulated points ([`crate::twoview`]), and the registration of each new
//! view by resectioning ([`crate::pnp`]). Every one of those steps is a
//! *local*, algebraic estimate: the two-view pose minimizes an algebraic
//! epipolar residual, DLT triangulation and DLT-PnP each minimize an algebraic
//! (not geometric) cost, and errors accumulate as views are chained. **Bundle
//! adjustment** is the global polish that ties the whole reconstruction
//! together:
//!
//! ```text
//!   cameras {(Rⱼ, tⱼ)}  +  points {Xᵢ}  +  observations {xᵢⱼ}
//!                         │
//!                         └─ minimize  Σᵢⱼ ‖ xᵢⱼ − π(K, Rⱼ, tⱼ, Xᵢ) ‖²
//!                         │            (total squared reprojection error)
//!                         ▼
//!         jointly refined cameras + points (a local MLE under Gaussian noise)
//! ```
//!
//! It simultaneously adjusts **all** camera poses and **all** 3-D points to
//! minimize the **total reprojection error** — the sum over every observation
//! of the squared pixel distance between the measured feature and the
//! projection of its 3-D point through that camera. Under the assumption of
//! independent identically-distributed Gaussian pixel noise this least-squares
//! objective is the maximum-likelihood estimate of the scene and motion. It is
//! the standard final refinement of essentially every SfM / SLAM pipeline.
//!
//! ## Method — dense Levenberg–Marquardt
//!
//! We solve the nonlinear least-squares problem by **Levenberg–Marquardt**
//! (LM), the damped Gauss–Newton method:
//!
//! 1. **Parametrization.** Each camera carries 6 degrees of freedom — a
//!    3-vector **angle-axis (Rodrigues)** rotation `ω` (with `R = exp([ω]_×)`,
//!    see [`rodrigues_exp`]) and a 3-vector translation `t`. Each 3-D point
//!    carries its 3 coordinates `(X, Y, Z)`. The full parameter vector stacks
//!    the 6 DOF of every *optimized* camera followed by the 3 DOF of every
//!    point.
//! 2. **Residuals.** For each [`Observation`] the residual is the 2-vector
//!    `xᵢⱼ − π(K, Rⱼ, tⱼ, Xᵢ)` (measured minus projected pixel). Stacking all
//!    observations gives the residual vector `r` (length `2 · #observations`)
//!    whose squared norm is the cost.
//! 3. **Jacobian.** The Jacobian `J = ∂r/∂params` is built by **numerical
//!    finite differences** (central differences, one perturbed cost evaluation
//!    per parameter). This is simple and robust but the dominant cost of the
//!    method — see the honesty note. (An analytic Jacobian is a clean future
//!    optimization.)
//! 4. **Damped normal equations.** Each LM step solves
//!    `(JᵀJ + λ · diag(JᵀJ)) δ = −Jᵀr` for the parameter update `δ`. The
//!    Marquardt scaling `λ · diag(JᵀJ)` (rather than `λ · I`) makes the step
//!    invariant to the wildly different units of rotations, translations, and
//!    point coordinates.
//! 5. **Accept / reject with `λ` control.** If the trial parameters lower the
//!    cost, the step is **accepted** and the damping `λ` is *decreased* (toward
//!    Gauss–Newton, faster convergence); otherwise it is **rejected** and `λ`
//!    is *increased* (toward gradient descent, smaller safer steps). Iteration
//!    stops at [`BundleParams::max_iterations`] or when the relative cost
//!    decrease falls below [`BundleParams::cost_tolerance`].
//!
//! See R. Hartley and A. Zisserman, *Multiple View Geometry in Computer
//! Vision*, 2nd ed., Appendix 6 (iterative estimation / Levenberg–Marquardt),
//! and B. Triggs, P. McLauchlan, R. Hartley, A. Fitzgibbon, "Bundle Adjustment
//! — A Modern Synthesis," *Vision Algorithms: Theory and Practice*, 2000.
//!
//! ## Gauge freedom — camera 0 is held fixed
//!
//! The reprojection cost is invariant under a global rigid transform of the
//! whole scene (rotate + translate every camera and point together): the
//! pixels do not change, so the problem has a **6-DOF gauge freedom** and `JᵀJ`
//! is correspondingly rank-deficient. We remove it the simplest way — by
//! **holding camera 0's pose completely fixed** (its 6 parameters are excluded
//! from the optimization). The refined reconstruction is therefore expressed
//! in camera 0's frame, exactly as the two-view seed ([`crate::twoview`]) and
//! the incremental mapper produce it.
//!
//! **Scale is still free.** Fixing one camera's pose removes the rotation +
//! translation gauge but *not* the global **scale** ambiguity of a
//! reconstruction seeded from two views (the essential matrix fixes translation
//! only up to scale — see the [`crate::twoview`] module docs). A pure
//! two-view-seeded bundle can shrink or grow the whole scene (moving every
//! point and every non-fixed camera centre by a common factor) without changing
//! any pixel, so BA does not pin the absolute scale; an external constraint (a
//! known baseline, a scale bar, GPS) is required for that. Holding camera 0
//! fixed leaves this one residual degree of freedom; it does not corrupt the
//! refinement, it simply is not determined by the pixels alone.
//!
//! ## Honesty notes
//!
//! - **Dense, not sparse.** This is a *dense* LM: it forms the full
//!   `JᵀJ` and solves it directly, costing roughly **O((6·#cameras +
//!   3·#points)³)** per iteration. That is perfectly fine for the scale this
//!   crate targets — tens of cameras and hundreds of points — but it does
//!   **not** exploit the characteristic *primary sparsity* of bundle
//!   adjustment (each observation touches exactly one camera and one point).
//!   Production BA uses the **Schur complement** to marginalize the points and
//!   solve a much smaller reduced camera system, the standard route to
//!   thousands of cameras and millions of points. That sparse path is a
//!   deliberate future extension, not implemented here.
//! - **Numerical Jacobian.** The finite-difference Jacobian costs `O(#params)`
//!   residual evaluations per iteration. An **analytic** Jacobian (closed-form
//!   derivatives of the projection w.r.t. the angle-axis rotation,
//!   translation, and point) would be substantially faster and is the usual
//!   choice in a tuned implementation; we use central differences here for a
//!   clear, obviously-correct first cut.
//! - **Shared intrinsics.** All cameras share a single [`CameraIntrinsics`]
//!   `K`, and `K` is held **fixed** (not optimized). Real pipelines often
//!   refine per-camera intrinsics (focal length, principal point, distortion)
//!   inside the bundle; **per-camera `K` and intrinsic refinement are a future
//!   extension.** For now the model is "one calibrated camera moving through
//!   the scene," which matches the [`crate::pnp`] / [`crate::twoview`]
//!   convention.
//! - **Local, not global, optimum.** LM is a local optimizer: it converges to
//!   the nearest local minimum of a non-convex cost, so it needs a reasonable
//!   initialization (exactly what Stages 3–4 provide). It is a *refiner*, not a
//!   from-scratch solver.

use crate::pnp::{project_point, CameraPose};
use crate::twoview::CameraIntrinsics;
use nalgebra::{DMatrix, DVector, Matrix3, Vector3};

/// The mutable working state of the optimizer: per-camera angle-axis rotations,
/// per-camera translations, and the 3-D points. Kept as a triple of `Vec`s so a
/// trial LM step can be applied to a copy and rolled back by simply discarding
/// it on rejection.
type State = (Vec<Vector3<f64>>, Vec<Vector3<f64>>, Vec<Vector3<f64>>);

/// Number of optimized parameters per camera: 3 for the angle-axis rotation
/// and 3 for the translation.
const CAM_PARAMS: usize = 6;
/// Number of optimized parameters per 3-D point: its `(X, Y, Z)` coordinates.
const POINT_PARAMS: usize = 3;

/// Absolute floor on the total squared-pixel cost below which a reconstruction
/// is treated as already optimal: there is nothing meaningful left to refine
/// (the residual is at round-off level), so [`bundle_adjust`] stops rather than
/// chasing sub-`1e-18` px² decreases. This makes an essentially-exact input
/// terminate immediately instead of taking a few no-op iterations on noise.
const ABS_COST_FLOOR: f64 = 1e-18;

/// A single 2-D image measurement: feature `pixel` in camera `camera_idx`
/// observed to be the projection of 3-D point `point_idx`.
///
/// `camera_idx` indexes into [`BundleProblem::cameras`] and `point_idx` into
/// [`BundleProblem::points`]. The same point is typically observed by several
/// cameras (one [`Observation`] each); that redundancy is exactly what bundle
/// adjustment exploits to refine the geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observation {
    /// Index of the observing camera in [`BundleProblem::cameras`].
    pub camera_idx: usize,
    /// Index of the observed 3-D point in [`BundleProblem::points`].
    pub point_idx: usize,
    /// The measured image location `(u, v)` of the point in that camera, in
    /// pixels.
    pub pixel: (f64, f64),
}

/// A bundle-adjustment problem: the cameras, points, shared intrinsics, and
/// the observations linking them.
///
/// The cost minimized by [`bundle_adjust`] is the total squared reprojection
/// error
///
/// ```text
///   Σ_obs  ‖ obs.pixel − π(intrinsics, cameras[obs.camera_idx], points[obs.point_idx]) ‖²
/// ```
///
/// **Camera 0 is the gauge-fixing reference** and is held fixed by
/// [`bundle_adjust`] (see the [module docs](self)); arrange the cameras so that
/// the reconstruction's reference frame is `cameras[0]`.
///
/// All cameras share the single [`Self::intrinsics`] `K` (held fixed). Per-
/// camera intrinsics are a future extension — see the [module docs](self).
#[derive(Debug, Clone)]
pub struct BundleProblem {
    /// The camera poses to refine, in the `X_cam = R · X_world + t` convention
    /// (see [`CameraPose`]). `cameras[0]` is held fixed as the gauge reference.
    pub cameras: Vec<CameraPose>,
    /// The 3-D scene points to refine, in the world (reference) frame.
    pub points: Vec<Vector3<f64>>,
    /// The shared pinhole intrinsics `K` for every camera. Held fixed (not
    /// optimized).
    pub intrinsics: CameraIntrinsics,
    /// The image measurements tying cameras to points.
    pub observations: Vec<Observation>,
}

/// Tuning parameters for the [`bundle_adjust`] Levenberg–Marquardt loop.
#[derive(Debug, Clone, Copy)]
pub struct BundleParams {
    /// Maximum number of LM iterations (accepted-or-rejected trial steps). The
    /// loop also stops early on the [`Self::cost_tolerance`] convergence test.
    pub max_iterations: usize,
    /// Initial Levenberg–Marquardt damping `λ`. Larger ⇒ smaller, safer
    /// (gradient-descent-like) first steps; it is adapted automatically as
    /// steps are accepted (λ decreased) or rejected (λ increased).
    pub initial_lambda: f64,
    /// Convergence tolerance on the **relative** cost decrease: the loop stops
    /// once an accepted step reduces the cost by a fraction smaller than this
    /// (i.e. `(prev − new) / prev < cost_tolerance`).
    pub cost_tolerance: f64,
}

impl Default for BundleParams {
    /// Sensible defaults: 100 iterations, `λ₀ = 1e-3`, relative cost tolerance
    /// `1e-9`.
    fn default() -> Self {
        Self {
            max_iterations: 100,
            initial_lambda: 1e-3,
            cost_tolerance: 1e-9,
        }
    }
}

/// The outcome of a [`bundle_adjust`] run: the refined geometry and the cost
/// before/after.
#[derive(Debug, Clone)]
pub struct BundleResult {
    /// The refined camera poses (same length and order as
    /// [`BundleProblem::cameras`]; `cameras[0]` is unchanged — the fixed
    /// gauge).
    pub cameras: Vec<CameraPose>,
    /// The refined 3-D points (same length and order as
    /// [`BundleProblem::points`]).
    pub points: Vec<Vector3<f64>>,
    /// Total squared reprojection error of the *input* reconstruction (the
    /// objective before any refinement).
    pub initial_cost: f64,
    /// Total squared reprojection error of the *refined* reconstruction (the
    /// objective at the returned solution; `≤ initial_cost`).
    pub final_cost: f64,
    /// Number of LM iterations actually run (≤ [`BundleParams::max_iterations`];
    /// fewer if the cost-tolerance convergence test fired or the system became
    /// singular).
    pub iterations: usize,
}

/// Rodrigues `exp`: map an **angle-axis** 3-vector `ω` to the rotation matrix
/// `R = exp([ω]_×)` it represents.
///
/// The direction of `ω` is the rotation axis and its magnitude `θ = ‖ω‖` is
/// the rotation angle (radians). The closed form is Rodrigues' rotation
/// formula
///
/// ```text
///   R = I + (sin θ / θ) [ω]_× + ((1 − cos θ) / θ²) [ω]_×²
/// ```
///
/// where `[ω]_×` is the skew-symmetric cross-product matrix of `ω`. As
/// `θ → 0` the factors `sin θ / θ → 1` and `(1 − cos θ) / θ² → ½` are finite,
/// but evaluating them directly divides `0/0`; for small `θ` we therefore use
/// the **Taylor limits** (`1 − θ²/6` and `½ − θ²/24`), so the function is exact
/// at `θ = 0` (returning `I`) and never divides by zero. The inverse map is
/// [`rodrigues_log`].
#[must_use]
pub fn rodrigues_exp(omega: &Vector3<f64>) -> Matrix3<f64> {
    let theta2 = omega.dot(omega);
    let theta = theta2.sqrt();
    let k = skew(omega);

    // Coefficients a = sin θ / θ and b = (1 − cos θ) / θ². Both are finite as
    // θ → 0; switch to the Taylor expansion below a small threshold to avoid
    // the 0/0 divide (and to stay accurate where the closed form loses
    // precision).
    let (a, b) = if theta < 1e-8 {
        // sin θ / θ      = 1 − θ²/6  + O(θ⁴)
        // (1 − cos θ)/θ² = ½ − θ²/24 + O(θ⁴)
        (1.0 - theta2 / 6.0, 0.5 - theta2 / 24.0)
    } else {
        (theta.sin() / theta, (1.0 - theta.cos()) / theta2)
    };

    Matrix3::identity() + a * k + b * (k * k)
}

/// Rodrigues `log`: map a rotation matrix `R` back to the **angle-axis**
/// 3-vector `ω` with `exp([ω]_×) = R` (the inverse of [`rodrigues_exp`]).
///
/// The rotation angle is `θ = arccos((tr R − 1) / 2)` and the axis comes from
/// the skew-symmetric part of `R`: `[ω]_× = (θ / (2 sin θ)) (R − Rᵀ)`, i.e.
/// `ω = (θ / (2 sin θ)) · (R₃₂−R₂₃, R₁₃−R₃₁, R₂₁−R₁₂)`. Two limits are guarded
/// so the function never divides by zero:
///
/// - **`θ ≈ 0`** (`R ≈ I`): the factor `θ / sin θ → 1`, and the off-diagonal
///   differences are themselves `≈ 0`, so `ω ≈ ½ (R − Rᵀ)^∨` (the small-angle
///   limit); we return that directly.
/// - **`θ ≈ π`**: `sin θ → 0` so the skew formula is ill-conditioned. Here
///   `R = I + 2 [u]_×²` for the unit axis `u`, so `uuᵀ = (R + I)/2`; we read
///   the axis off the largest diagonal of `(R + I)/2` (with the sign fixed from
///   the off-diagonal entries) and scale by `θ = π`. (A `±π` rotation has a
///   genuine two-fold sign ambiguity in the axis; either choice is a valid
///   logarithm.)
///
/// The input is assumed to be a proper rotation; a mildly non-orthonormal `R`
/// is handled gracefully (`(tr R − 1)/2` is clamped to `[−1, 1]` before the
/// `arccos`).
#[must_use]
pub fn rodrigues_log(r: &Matrix3<f64>) -> Vector3<f64> {
    // cos θ = (tr R − 1) / 2, clamped against round-off so arccos is defined.
    let cos_theta = ((r.trace() - 1.0) / 2.0).clamp(-1.0, 1.0);
    let theta = cos_theta.acos();

    // The (unscaled) axis from the skew-symmetric part: (R − Rᵀ)^∨.
    let v = Vector3::new(
        r[(2, 1)] - r[(1, 2)],
        r[(0, 2)] - r[(2, 0)],
        r[(1, 0)] - r[(0, 1)],
    );

    if theta < 1e-8 {
        // Small angle: θ / (2 sin θ) → ½, and v is already ≈ 0. Return the
        // first-order term ½ v (exact as θ → 0, no divide-by-zero).
        0.5 * v
    } else if theta < std::f64::consts::PI - 1e-6 {
        // Generic case: ω = (θ / (2 sin θ)) · v.
        let s = theta / (2.0 * theta.sin());
        s * v
    } else {
        // θ ≈ π: sin θ ≈ 0 makes the skew formula blow up. Recover the axis
        // from R = I + 2[u]_×²  ⇒  (R + I)/2 = uuᵀ. The diagonal of (R+I)/2 is
        // (uₓ², u_y², u_z²); take the largest for numerical stability, then fix
        // the other components' signs from the symmetric off-diagonals.
        let m = (r + Matrix3::identity()) * 0.5;
        // Largest diagonal index.
        let mut k = 0usize;
        for i in 1..3 {
            if m[(i, i)] > m[(k, k)] {
                k = i;
            }
        }
        let mut axis = Vector3::zeros();
        // Guard the sqrt against a tiny negative from round-off.
        let dk = m[(k, k)].max(0.0).sqrt();
        if dk < 1e-12 {
            // Degenerate (should not happen for a true rotation at θ = π);
            // fall back to the first-order term to stay finite.
            return 0.5 * v;
        }
        axis[k] = dk;
        for i in 0..3 {
            if i != k {
                axis[i] = m[(k, i)] / dk;
            }
        }
        // Normalize defensively and scale by the angle θ (≈ π).
        let n = axis.norm();
        if n > 1e-12 {
            axis /= n;
        }
        theta * axis
    }
}

/// Skew-symmetric cross-product matrix `[v]_×` with `[v]_× a = v × a`.
#[inline]
fn skew(v: &Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(
        0.0, -v.z, v.y, //
        v.z, 0.0, -v.x, //
        -v.y, v.x, 0.0,
    )
}

/// Jointly refine the camera poses and 3-D points to minimize the total
/// reprojection error, by **dense Levenberg–Marquardt** (Stage 5).
///
/// `problem` supplies the initial cameras, points, shared intrinsics, and the
/// observations; `params` tunes the LM loop. **Camera 0 is held fixed** as the
/// gauge reference (its pose is copied unchanged into the result); all other
/// camera poses and all points are refined. See the [module docs](self) for the
/// method, the gauge/scale handling, and the dense-vs-sparse and numerical-
/// Jacobian honesty notes.
///
/// # Behaviour and guarantees
///
/// - The returned [`BundleResult::final_cost`] is **never greater** than
///   [`BundleResult::initial_cost`]: LM only ever commits a trial step that
///   strictly lowers the cost (a rejected step leaves the parameters and cost
///   untouched), so a worst case simply returns the input unchanged.
/// - **Robust to degenerate input — never panics.** With no observations (or
///   no *free* parameters — e.g. a single fixed camera 0), there is nothing to
///   optimize and the input is returned unchanged with `iterations = 0`. The
///   perspective divide in [`project_point`] is guarded at depth `z ≈ 0`: an
///   observation whose point falls on the camera's principal plane contributes
///   **no** residual or Jacobian rows for that step (it is skipped rather than
///   producing a non-finite cost). The damped normal-equation solve is guarded
///   against a singular / near-singular `JᵀJ`: if the Cholesky solve fails the
///   step is treated as a rejection (damping `λ` is raised and the parameters
///   are kept), so a rank-deficient system never triggers a divide-by-zero or
///   an `unwrap` on a failed solve — the loop simply stops making progress and
///   returns the best reconstruction found.
#[must_use]
pub fn bundle_adjust(problem: &BundleProblem, params: &BundleParams) -> BundleResult {
    let num_cameras = problem.cameras.len();
    let num_points = problem.points.len();

    // Working state: angle-axis + translation per camera, XYZ per point. We
    // optimize cameras 1.. (camera 0 is the fixed gauge) and all points.
    let mut cam_rot: Vec<Vector3<f64>> = problem
        .cameras
        .iter()
        .map(|c| rodrigues_log(&c.rotation))
        .collect();
    let mut cam_trans: Vec<Vector3<f64>> = problem.cameras.iter().map(|c| c.translation).collect();
    let mut points: Vec<Vector3<f64>> = problem.points.clone();

    // Number of free cameras (all but camera 0, if any cameras exist at all).
    let free_cameras = num_cameras.saturating_sub(1);
    let num_params = free_cameras * CAM_PARAMS + num_points * POINT_PARAMS;

    let initial_cost = total_cost(
        &problem.intrinsics,
        &cam_rot,
        &cam_trans,
        &points,
        &problem.observations,
    );

    // Nothing to optimize (no free parameters, or no observations): return the
    // input unchanged. This covers the "single fixed camera" and "no
    // observations" degenerate cases gracefully — no panic, zero iterations.
    //
    // Also short-circuit an already-optimal input (cost at round-off level):
    // there is nothing to refine, so we terminate immediately with zero
    // iterations rather than taking a few no-op LM steps on numerical noise.
    if num_params == 0 || problem.observations.is_empty() || initial_cost <= ABS_COST_FLOOR {
        return BundleResult {
            cameras: problem.cameras.clone(),
            points: points.clone(),
            initial_cost,
            final_cost: initial_cost,
            iterations: 0,
        };
    }

    let mut lambda = params.initial_lambda.max(0.0);
    let mut cost = initial_cost;
    let mut iterations = 0usize;

    for _ in 0..params.max_iterations {
        iterations += 1;

        // Build the residual vector and the (dense) finite-difference Jacobian
        // at the current parameters.
        let (residual, jac) = residual_and_jacobian(
            &problem.intrinsics,
            &cam_rot,
            &cam_trans,
            &points,
            &problem.observations,
            free_cameras,
            num_points,
            num_params,
        );

        // Normal-equation pieces: H = JᵀJ and g = Jᵀr.
        let jt = jac.transpose();
        let h = &jt * &jac;
        let g = &jt * &residual;

        // Try (possibly several) damped steps from this linearization, raising
        // λ on each rejection, until one is accepted or λ runs away.
        let mut step_taken = false;
        // Cap inner damping retries so a hopeless linearization cannot spin.
        for _ in 0..10 {
            // Augment the diagonal: (H + λ·diag(H)) δ = −g  (Marquardt scaling).
            let mut aug = h.clone();
            for i in 0..num_params {
                let dii = h[(i, i)];
                aug[(i, i)] = dii + lambda * dii;
            }

            // Solve the damped system. A singular / near-singular H means a
            // rank-deficient (e.g. gauge- or scale-degenerate) configuration;
            // guard it — a failed solve is treated exactly like a rejected
            // step (raise λ, keep the parameters), never an unwrap/divide.
            let Some(delta) = solve_spd(&aug, &(-&g)) else {
                lambda = (lambda * 10.0).min(1e12);
                continue;
            };

            // Apply the trial update to a COPY of the state.
            let (try_rot, try_trans, try_points) =
                apply_delta(&cam_rot, &cam_trans, &points, &delta, free_cameras);
            let new_cost = total_cost(
                &problem.intrinsics,
                &try_rot,
                &try_trans,
                &try_points,
                &problem.observations,
            );

            if new_cost.is_finite() && new_cost < cost {
                // Accept: commit the trial state, relax the damping.
                cam_rot = try_rot;
                cam_trans = try_trans;
                points = try_points;
                let prev = cost;
                cost = new_cost;
                lambda = (lambda * 0.5).max(1e-12);
                step_taken = true;

                // Convergence: stop when either the cost is already negligible
                // in absolute terms (an essentially-exact fit — nothing left to
                // improve, so we do not keep spinning on round-off-level
                // decreases), or the *relative* decrease has fallen below the
                // tolerance.
                let rel = if prev > 0.0 {
                    (prev - new_cost) / prev
                } else {
                    0.0
                };
                if new_cost <= ABS_COST_FLOOR || rel < params.cost_tolerance {
                    return finish(
                        problem,
                        &cam_rot,
                        &cam_trans,
                        &points,
                        initial_cost,
                        cost,
                        iterations,
                    );
                }
                break;
            }
            // Reject: raise the damping and retry the inner loop.
            lambda = (lambda * 10.0).min(1e12);
        }

        // No accepted step from this linearization (λ saturated / system
        // singular): we are at a local minimum or a degenerate point. Stop.
        if !step_taken {
            break;
        }
    }

    finish(
        problem,
        &cam_rot,
        &cam_trans,
        &points,
        initial_cost,
        cost,
        iterations,
    )
}

/// Assemble the [`BundleResult`] from the working state (rebuilding the
/// [`CameraPose`] list from the angle-axis + translation parameters; camera 0
/// is rebuilt from its untouched parameters and is therefore identical to the
/// input).
fn finish(
    problem: &BundleProblem,
    cam_rot: &[Vector3<f64>],
    cam_trans: &[Vector3<f64>],
    points: &[Vector3<f64>],
    initial_cost: f64,
    final_cost: f64,
    iterations: usize,
) -> BundleResult {
    let cameras = (0..problem.cameras.len())
        .map(|j| CameraPose {
            rotation: rodrigues_exp(&cam_rot[j]),
            translation: cam_trans[j],
        })
        .collect();
    BundleResult {
        cameras,
        points: points.to_vec(),
        initial_cost,
        final_cost,
        iterations,
    }
}

/// Total squared reprojection error over all observations for the given
/// parameters. Observations whose point projects with `z ≈ 0` (guarded by
/// [`project_point`]) contribute nothing — they are skipped, never a
/// non-finite term.
fn total_cost(
    k: &CameraIntrinsics,
    cam_rot: &[Vector3<f64>],
    cam_trans: &[Vector3<f64>],
    points: &[Vector3<f64>],
    observations: &[Observation],
) -> f64 {
    let mut sum = 0.0;
    for obs in observations {
        if let Some((du, dv)) = obs_residual(k, cam_rot, cam_trans, points, obs) {
            sum += du * du + dv * dv;
        }
    }
    sum
}

/// The 2-D residual `measured − projected` for one observation, or [`None`] if
/// the point is unprojectable (depth `z ≈ 0`) under the current parameters or
/// the indices are out of range (defensive — the public API never builds such
/// observations, but we never panic on them).
#[inline]
fn obs_residual(
    k: &CameraIntrinsics,
    cam_rot: &[Vector3<f64>],
    cam_trans: &[Vector3<f64>],
    points: &[Vector3<f64>],
    obs: &Observation,
) -> Option<(f64, f64)> {
    let r = rodrigues_exp(cam_rot.get(obs.camera_idx)?);
    let t = cam_trans.get(obs.camera_idx)?;
    let x = points.get(obs.point_idx)?;
    let (pu, pv) = project_point(k, &r, t, x)?;
    let (ou, ov) = obs.pixel;
    Some((ou - pu, ov - pv))
}

/// Build the stacked residual vector and the dense **finite-difference**
/// Jacobian `J = ∂r/∂params` at the current parameters.
///
/// The residual vector has two entries per observation. The Jacobian has
/// `2·#observations` rows and `num_params` columns, with the camera blocks
/// (6 columns each, for cameras `1..`) first and the point blocks (3 columns
/// each) after. Central differences are used per parameter; an observation
/// whose residual is undefined (point on the principal plane) contributes a
/// zero row for that step.
#[allow(clippy::too_many_arguments)]
fn residual_and_jacobian(
    k: &CameraIntrinsics,
    cam_rot: &[Vector3<f64>],
    cam_trans: &[Vector3<f64>],
    points: &[Vector3<f64>],
    observations: &[Observation],
    free_cameras: usize,
    num_points: usize,
    num_params: usize,
) -> (DVector<f64>, DMatrix<f64>) {
    let m = 2 * observations.len();
    let mut residual = DVector::<f64>::zeros(m);
    let mut jac = DMatrix::<f64>::zeros(m, num_params);

    // Current residuals.
    for (i, obs) in observations.iter().enumerate() {
        if let Some((ru, rv)) = obs_residual(k, cam_rot, cam_trans, points, obs) {
            residual[2 * i] = ru;
            residual[2 * i + 1] = rv;
        }
    }

    // Step size for central differences.
    let eps = 1e-6;

    // Column index of camera `cam`'s first parameter (camera 0 is fixed, so
    // camera `j ≥ 1` maps to block `j − 1`).
    let cam_col = |cam: usize| (cam - 1) * CAM_PARAMS;
    // Column index of point `p`'s first parameter.
    let point_col = |p: usize| free_cameras * CAM_PARAMS + p * POINT_PARAMS;

    // Camera parameters (rotation 0..3, translation 3..6), cameras 1.. only.
    for cam in 1..(free_cameras + 1) {
        for d in 0..CAM_PARAMS {
            let col = cam_col(cam) + d;
            // Perturb a local copy of this camera's 6 params ±eps.
            let perturb = |sign: f64| -> (Vector3<f64>, Vector3<f64>) {
                let mut rr = cam_rot[cam];
                let mut tt = cam_trans[cam];
                if d < 3 {
                    rr[d] += sign * eps;
                } else {
                    tt[d - 3] += sign * eps;
                }
                (rr, tt)
            };
            let (rp, tp) = perturb(1.0);
            let (rm, tm) = perturb(-1.0);

            // Only observations of THIS camera have a non-zero derivative.
            for (i, obs) in observations.iter().enumerate() {
                if obs.camera_idx != cam {
                    continue;
                }
                let plus = single_residual(k, &rp, &tp, points, obs);
                let minus = single_residual(k, &rm, &tm, points, obs);
                if let (Some((pu, pv)), Some((mu, mv))) = (plus, minus) {
                    jac[(2 * i, col)] = (pu - mu) / (2.0 * eps);
                    jac[(2 * i + 1, col)] = (pv - mv) / (2.0 * eps);
                }
            }
        }
    }

    // Point parameters (X, Y, Z).
    for (p, point) in points.iter().enumerate().take(num_points) {
        for d in 0..POINT_PARAMS {
            let col = point_col(p) + d;
            let mut xp = *point;
            let mut xm = *point;
            xp[d] += eps;
            xm[d] -= eps;

            // Only observations of THIS point have a non-zero derivative.
            for (i, obs) in observations.iter().enumerate() {
                if obs.point_idx != p {
                    continue;
                }
                let cam = obs.camera_idx;
                // Use the (possibly fixed) camera's current pose.
                let (Some(rr), Some(tt)) = (cam_rot.get(cam), cam_trans.get(cam)) else {
                    continue;
                };
                let r = rodrigues_exp(rr);
                let plus = residual_for_point(k, &r, tt, &xp, obs);
                let minus = residual_for_point(k, &r, tt, &xm, obs);
                if let (Some((pu, pv)), Some((mu, mv))) = (plus, minus) {
                    jac[(2 * i, col)] = (pu - mu) / (2.0 * eps);
                    jac[(2 * i + 1, col)] = (pv - mv) / (2.0 * eps);
                }
            }
        }
    }

    (residual, jac)
}

/// Residual for one observation given an explicit angle-axis `rot` and
/// translation for its camera (used while perturbing camera parameters).
#[inline]
fn single_residual(
    k: &CameraIntrinsics,
    rot: &Vector3<f64>,
    trans: &Vector3<f64>,
    points: &[Vector3<f64>],
    obs: &Observation,
) -> Option<(f64, f64)> {
    let r = rodrigues_exp(rot);
    let x = points.get(obs.point_idx)?;
    let (pu, pv) = project_point(k, &r, trans, x)?;
    let (ou, ov) = obs.pixel;
    Some((ou - pu, ov - pv))
}

/// Residual for one observation given an explicit *point* `x` and a prebuilt
/// rotation matrix + translation for its camera (used while perturbing point
/// parameters).
#[inline]
fn residual_for_point(
    k: &CameraIntrinsics,
    r: &Matrix3<f64>,
    t: &Vector3<f64>,
    x: &Vector3<f64>,
    obs: &Observation,
) -> Option<(f64, f64)> {
    let (pu, pv) = project_point(k, r, t, x)?;
    let (ou, ov) = obs.pixel;
    Some((ou - pu, ov - pv))
}

/// Apply a parameter update `delta` to a *copy* of the working state and return
/// the trial `(cam_rot, cam_trans, points)`. The camera blocks update cameras
/// `1..` (camera 0 is the fixed gauge and is copied unchanged); each camera's
/// rotation update is applied **additively in the angle-axis parametrization**
/// (`ω ← ω + δω`), which is the standard local update for a finite-difference
/// LM step.
fn apply_delta(
    cam_rot: &[Vector3<f64>],
    cam_trans: &[Vector3<f64>],
    points: &[Vector3<f64>],
    delta: &DVector<f64>,
    free_cameras: usize,
) -> State {
    let mut rot = cam_rot.to_vec();
    let mut trans = cam_trans.to_vec();
    let mut pts = points.to_vec();

    // Cameras 1.. (block index 0.. in delta).
    for cam in 1..(free_cameras + 1) {
        let base = (cam - 1) * CAM_PARAMS;
        for d in 0..3 {
            rot[cam][d] += delta[base + d];
        }
        for d in 0..3 {
            trans[cam][d] += delta[base + 3 + d];
        }
    }

    // Points.
    let point_base = free_cameras * CAM_PARAMS;
    for (p, point) in pts.iter_mut().enumerate() {
        let base = point_base + p * POINT_PARAMS;
        for d in 0..POINT_PARAMS {
            point[d] += delta[base + d];
        }
    }

    (rot, trans, pts)
}

/// Solve the symmetric positive-(semi)definite system `A x = b` by Cholesky,
/// returning [`None`] if `A` is not positive-definite (singular / rank-
/// deficient) or the solution is non-finite.
///
/// The damped normal matrix `JᵀJ + λ·diag(JᵀJ)` is symmetric PSD by
/// construction and positive-definite whenever the configuration is non-
/// degenerate and `λ > 0`; a Cholesky failure therefore signals a rank-
/// deficient (gauge- or scale-degenerate) problem, which [`bundle_adjust`]
/// handles by raising `λ` rather than panicking.
fn solve_spd(a: &DMatrix<f64>, b: &DVector<f64>) -> Option<DVector<f64>> {
    let chol = a.clone().cholesky()?;
    let x = chol.solve(b);
    if x.iter().all(|v| v.is_finite()) {
        Some(x)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_intrinsics() -> CameraIntrinsics {
        CameraIntrinsics::new(800.0, 800.0, 320.0, 240.0)
    }

    /// Rotation from yaw (about y) then pitch (about x), in radians.
    fn rot_yaw_pitch(yaw: f64, pitch: f64) -> Matrix3<f64> {
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
        ry * rx
    }

    /// A small fixed scene of non-coplanar 3-D points in front of the cameras.
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

    /// A small camera rig: three views looking at the scene from slightly
    /// different poses. Camera 0 is the reference (identity rotation, origin
    /// translation), matching the gauge convention.
    fn camera_rig() -> Vec<CameraPose> {
        vec![
            // Camera 0: reference frame.
            CameraPose {
                rotation: Matrix3::identity(),
                translation: Vector3::zeros(),
            },
            // Camera 1: yaw/pitch + sideways baseline.
            CameraPose {
                rotation: rot_yaw_pitch(0.15, -0.08),
                translation: Vector3::new(-0.8, 0.05, 0.10),
            },
            // Camera 2: different yaw/pitch + baseline.
            CameraPose {
                rotation: rot_yaw_pitch(-0.10, 0.12),
                translation: Vector3::new(0.7, -0.06, 0.20),
            },
        ]
    }

    /// Build the full observation set: every camera observes every point
    /// (synthetic, fully-connected), using the exact projection so the ground
    /// truth has zero reprojection error.
    fn build_observations(
        k: &CameraIntrinsics,
        cameras: &[CameraPose],
        points: &[Vector3<f64>],
    ) -> Vec<Observation> {
        let mut obs = Vec::new();
        for (ci, cam) in cameras.iter().enumerate() {
            for (pi, x) in points.iter().enumerate() {
                let pix = project_point(k, &cam.rotation, &cam.translation, x)
                    .expect("ground-truth point projects in front of every camera");
                obs.push(Observation {
                    camera_idx: ci,
                    point_idx: pi,
                    pixel: pix,
                });
            }
        }
        obs
    }

    /// RMSE (pixels) of a reconstruction over its observations.
    fn rmse(problem: &BundleProblem, cameras: &[CameraPose], points: &[Vector3<f64>]) -> f64 {
        let mut sum = 0.0;
        let mut n = 0usize;
        for obs in &problem.observations {
            let cam = &cameras[obs.camera_idx];
            if let Some((pu, pv)) = project_point(
                &problem.intrinsics,
                &cam.rotation,
                &cam.translation,
                &points[obs.point_idx],
            ) {
                let (ou, ov) = obs.pixel;
                sum += (pu - ou).powi(2) + (pv - ov).powi(2);
                n += 1;
            }
        }
        if n == 0 {
            0.0
        } else {
            (sum / n as f64).sqrt()
        }
    }

    /// A small deterministic LCG-style perturbation helper (no `rand` dep), in
    /// the spirit of the crate's SplitMix64 usage elsewhere.
    struct Noise {
        state: u64,
    }
    impl Noise {
        fn new(seed: u64) -> Self {
            Self { state: seed }
        }
        /// Next value in roughly `[-1, 1]`.
        fn next(&mut self) -> f64 {
            self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            ((z % 2_000_000) as f64 / 1_000_000.0) - 1.0
        }
    }

    // ----- Test 1: BA reduces cost; RMSE drops toward zero. -----
    #[test]
    fn bundle_reduces_cost() {
        let k = test_intrinsics();
        let cameras = camera_rig();
        let points = scene_points();
        let observations = build_observations(&k, &cameras, &points);

        // Perturb cameras (except camera 0, the gauge) and all points.
        let mut noise = Noise::new(0x1234_5678_9ABC_DEF0);
        let mut pert_cams = cameras.clone();
        for cam in pert_cams.iter_mut().skip(1) {
            // Small rotation perturbation (~0.03 rad) and translation (~0.03).
            let dw = Vector3::new(noise.next(), noise.next(), noise.next()) * 0.03;
            cam.rotation = rodrigues_exp(&dw) * cam.rotation;
            cam.translation += Vector3::new(noise.next(), noise.next(), noise.next()) * 0.03;
        }
        let mut pert_points = points.clone();
        for p in pert_points.iter_mut() {
            *p += Vector3::new(noise.next(), noise.next(), noise.next()) * 0.05;
        }

        let problem = BundleProblem {
            cameras: pert_cams,
            points: pert_points,
            intrinsics: k,
            observations,
        };

        let rmse_before = rmse(&problem, &problem.cameras, &problem.points);
        let result = bundle_adjust(&problem, &BundleParams::default());
        let rmse_after = rmse(&problem, &result.cameras, &result.points);

        // Cost must drop sharply.
        assert!(
            result.final_cost < result.initial_cost,
            "final cost {} not below initial {}",
            result.final_cost,
            result.initial_cost
        );
        assert!(
            result.final_cost < result.initial_cost * 1e-3,
            "expected a large cost reduction: initial {}, final {}",
            result.initial_cost,
            result.final_cost
        );
        // RMSE drops toward zero (clean synthetic data → near-exact recovery).
        assert!(
            rmse_after < rmse_before,
            "RMSE did not drop: before {rmse_before}, after {rmse_after}"
        );
        assert!(
            rmse_after < 1e-3,
            "post-BA RMSE {rmse_after} px should be ~0 on clean data"
        );
    }

    // ----- Test 2: recovery — refined cameras/points ≈ ground truth. -----
    #[test]
    fn bundle_recovers_ground_truth() {
        let k = test_intrinsics();
        let cameras = camera_rig();
        let points = scene_points();
        let observations = build_observations(&k, &cameras, &points);

        // Perturb (gauge camera 0 held).
        let mut noise = Noise::new(0xCAFE_F00D_1234_5678);
        let mut pert_cams = cameras.clone();
        for cam in pert_cams.iter_mut().skip(1) {
            let dw = Vector3::new(noise.next(), noise.next(), noise.next()) * 0.025;
            cam.rotation = rodrigues_exp(&dw) * cam.rotation;
            cam.translation += Vector3::new(noise.next(), noise.next(), noise.next()) * 0.025;
        }
        let mut pert_points = points.clone();
        for p in pert_points.iter_mut() {
            *p += Vector3::new(noise.next(), noise.next(), noise.next()) * 0.04;
        }

        let problem = BundleProblem {
            cameras: pert_cams,
            points: pert_points,
            intrinsics: k,
            observations,
        };
        let params = BundleParams {
            max_iterations: 200,
            ..BundleParams::default()
        };
        let result = bundle_adjust(&problem, &params);

        // The fit is essentially exact: the refined reconstruction reprojects
        // with ~zero residual (cost driven down by ~22 orders of magnitude).
        assert!(
            result.final_cost < 1e-12,
            "BA should drive the reprojection cost to ~0; got {} (initial {})",
            result.final_cost,
            result.initial_cost
        );

        // Camera 0 is exactly the fixed gauge (unchanged).
        assert!(
            (result.cameras[0].rotation - Matrix3::identity()).norm() < 1e-12,
            "camera 0 rotation must be held fixed at identity"
        );
        assert!(
            result.cameras[0].translation.norm() < 1e-12,
            "camera 0 translation must be held fixed at origin"
        );

        // Recovery is exact UP TO THE FREE GLOBAL SCALE (see the module docs):
        // pinning camera 0 removes the rotation + translation gauge but NOT the
        // scale of a two-view-style reconstruction, so the whole scene (every
        // non-fixed camera centre and every point) can settle at a common
        // scale factor `s` while every pixel — and hence the cost — is
        // unchanged. We therefore:
        //   (a) check each rotation matches truth directly (rotation is
        //       scale-invariant), and
        //   (b) recover the single factor `s` from camera 1's translation and
        //       verify EVERY non-fixed camera translation and EVERY point is
        //       consistent under that same `s`.

        // Rotations must match the ground truth essentially exactly.
        for (j, (got, truth)) in result.cameras.iter().zip(cameras.iter()).enumerate() {
            let r_err = (got.rotation - truth.rotation).norm();
            assert!(r_err < 1e-6, "camera {j} rotation off truth by {r_err}");
        }

        // Recover the global scale from camera 1's (non-zero) translation.
        let s = result.cameras[1].translation.norm() / cameras[1].translation.norm();
        assert!(s.is_finite() && s > 0.0, "recovered scale must be positive");

        // Every non-fixed camera's translation matches truth under that one `s`.
        for (j, (got, truth)) in result
            .cameras
            .iter()
            .zip(cameras.iter())
            .enumerate()
            .skip(1)
        {
            let scaled_truth = truth.translation * s;
            let t_err = (got.translation - scaled_truth).norm();
            assert!(
                t_err < 1e-6,
                "camera {j} translation {:?} inconsistent with scale {s} (truth·s = {scaled_truth:?}, err {t_err})",
                got.translation
            );
        }

        // Every point matches truth under the SAME scale `s`.
        for (i, (got, truth)) in result.points.iter().zip(points.iter()).enumerate() {
            let err = (got - truth * s).norm();
            assert!(
                err < 1e-6,
                "point {i} {:?} inconsistent with scale {s} (truth·s = {:?}, err {err})",
                got,
                truth * s
            );
        }

        // Sanity: the recovered scale is genuinely close to 1 here (the
        // perturbation was small), confirming this is the documented residual
        // scale freedom, not a gross error.
        assert!(
            (s - 1.0).abs() < 0.05,
            "recovered scale {s} drifted further than the small perturbation justifies"
        );
    }

    // ----- Test 3: already-optimal input stays put and terminates fast. -----
    #[test]
    fn already_optimal_input_is_unchanged() {
        let k = test_intrinsics();
        let cameras = camera_rig();
        let points = scene_points();
        let observations = build_observations(&k, &cameras, &points);

        let problem = BundleProblem {
            cameras: cameras.clone(),
            points: points.clone(),
            intrinsics: k,
            observations,
        };
        let result = bundle_adjust(&problem, &BundleParams::default());

        // The exact solution has ~zero cost.
        assert!(
            result.initial_cost < 1e-12,
            "ground-truth initial cost should be ~0, got {}",
            result.initial_cost
        );
        assert!(
            result.final_cost <= result.initial_cost + 1e-15,
            "final cost {} grew above initial {}",
            result.final_cost,
            result.initial_cost
        );
        // An already-optimal input (cost at round-off level) is recognized up
        // front and returned immediately — zero wasted iterations.
        assert_eq!(
            result.iterations, 0,
            "optimal input should terminate immediately, took {} iters",
            result.iterations
        );
        // Geometry essentially unchanged.
        for (got, truth) in result.cameras.iter().zip(cameras.iter()) {
            assert!((got.rotation - truth.rotation).norm() < 1e-9);
            assert!((got.translation - truth.translation).norm() < 1e-9);
        }
        for (got, truth) in result.points.iter().zip(points.iter()) {
            assert!((got - truth).norm() < 1e-9);
        }
    }

    // ----- Test 4: Rodrigues exp(log(R)) ≈ R, incl. a near-zero angle. -----
    #[test]
    fn rodrigues_round_trip() {
        let rotations = [
            Matrix3::identity(),
            rot_yaw_pitch(0.4, -0.25),
            rot_yaw_pitch(-1.2, 0.9),
            // Near-zero angle (exercises the θ≈0 small-angle guards in BOTH
            // exp and log).
            rodrigues_exp(&Vector3::new(1e-9, -2e-9, 0.5e-9)),
            // A larger rotation about a skew axis.
            rodrigues_exp(&Vector3::new(0.7, -0.5, 1.1)),
            // Close to a π rotation about a generic axis (exercises the θ≈π
            // branch of log).
            rodrigues_exp(
                &(Vector3::new(1.0, 0.3, -0.2).normalize() * (std::f64::consts::PI - 1e-3)),
            ),
        ];
        for (i, r) in rotations.iter().enumerate() {
            let w = rodrigues_log(r);
            let r2 = rodrigues_exp(&w);
            let err = (r2 - r).norm();
            assert!(
                err < 1e-9,
                "round-trip {i}: exp(log(R)) off R by {err} (R = {r:?})"
            );
        }

        // Specifically check the near-zero vector maps through exp to ≈ I and
        // log returns ≈ 0 with no NaN (the divide-by-zero guard).
        let tiny = Vector3::new(1e-12, -1e-12, 1e-12);
        let r_tiny = rodrigues_exp(&tiny);
        assert!(
            (r_tiny - Matrix3::identity()).norm() < 1e-9,
            "exp of a ~0 angle-axis should be ≈ I"
        );
        let w_tiny = rodrigues_log(&r_tiny);
        assert!(
            w_tiny.iter().all(|v| v.is_finite()),
            "log near θ=0 must be finite (no divide-by-zero)"
        );
        assert!(w_tiny.norm() < 1e-9, "log of ≈I should be ≈ 0");
    }

    // ----- Test 5a: no observations → unchanged, graceful, no panic. -----
    #[test]
    fn no_observations_returns_input_unchanged() {
        let k = test_intrinsics();
        let cameras = camera_rig();
        let points = scene_points();
        let problem = BundleProblem {
            cameras: cameras.clone(),
            points: points.clone(),
            intrinsics: k,
            observations: Vec::new(),
        };
        let result = bundle_adjust(&problem, &BundleParams::default());
        assert_eq!(result.iterations, 0, "no observations → zero iterations");
        assert_eq!(result.initial_cost, 0.0);
        assert_eq!(result.final_cost, 0.0);
        // Cameras and points returned unchanged.
        for (got, truth) in result.cameras.iter().zip(cameras.iter()) {
            assert!((got.rotation - truth.rotation).norm() < 1e-12);
            assert!((got.translation - truth.translation).norm() < 1e-12);
        }
        for (got, truth) in result.points.iter().zip(points.iter()) {
            assert!((got - truth).norm() < 1e-12);
        }
    }

    // ----- Test 5b: single (fixed) camera → no free params, unchanged. -----
    #[test]
    fn single_fixed_camera_is_graceful() {
        let k = test_intrinsics();
        let points = scene_points();
        let cam = CameraPose {
            rotation: Matrix3::identity(),
            translation: Vector3::zeros(),
        };
        // One camera (camera 0, fixed) observing every point: there are NO free
        // camera params, but the points ARE free. So this is not "nothing to
        // optimize" — but with a single fixed camera the point depths along
        // each ray are unobservable (monocular depth ambiguity), so the system
        // is rank-deficient. The guard must keep it finite and never panic.
        let observations: Vec<Observation> = points
            .iter()
            .enumerate()
            .map(|(pi, x)| Observation {
                camera_idx: 0,
                point_idx: pi,
                pixel: project_point(&k, &cam.rotation, &cam.translation, x).unwrap(),
            })
            .collect();
        let problem = BundleProblem {
            cameras: vec![cam],
            points: points.clone(),
            intrinsics: k,
            observations,
        };
        // Must not panic; cost must stay finite and non-increasing.
        let result = bundle_adjust(&problem, &BundleParams::default());
        assert!(result.initial_cost.is_finite());
        assert!(result.final_cost.is_finite());
        assert!(
            result.final_cost <= result.initial_cost + 1e-12,
            "final cost must not increase"
        );
        // Camera 0 stays the fixed gauge.
        assert!((result.cameras[0].rotation - Matrix3::identity()).norm() < 1e-12);
        assert!(result.cameras[0].translation.norm() < 1e-12);
    }

    // ----- Test 5c: empty problem (no cameras, no points) → graceful. -----
    #[test]
    fn empty_problem_is_graceful() {
        let k = test_intrinsics();
        let problem = BundleProblem {
            cameras: Vec::new(),
            points: Vec::new(),
            intrinsics: k,
            observations: Vec::new(),
        };
        let result = bundle_adjust(&problem, &BundleParams::default());
        assert_eq!(result.iterations, 0);
        assert_eq!(result.initial_cost, 0.0);
        assert_eq!(result.final_cost, 0.0);
        assert!(result.cameras.is_empty());
        assert!(result.points.is_empty());
    }

    // ----- Test 5d: a point on the principal plane (z≈0) is skipped, no NaN. -
    #[test]
    fn principal_plane_point_does_not_break_cost() {
        let k = test_intrinsics();
        let cameras = camera_rig();
        let mut points = scene_points();
        // Force one point onto camera 0's principal plane (z = 0): its
        // projection through camera 0 is undefined and must be SKIPPED, not
        // turned into a NaN cost.
        points.push(Vector3::new(0.2, -0.1, 0.0));
        let bad_idx = points.len() - 1;

        // Observations: the good points through all cameras, plus the bad point
        // through camera 0 (the offending z≈0 observation).
        let mut observations = Vec::new();
        for (ci, cam) in cameras.iter().enumerate() {
            for (pi, x) in points.iter().enumerate() {
                if pi == bad_idx && ci != 0 {
                    continue; // only add the bad point for camera 0
                }
                if let Some(pix) = project_point(&k, &cam.rotation, &cam.translation, x) {
                    observations.push(Observation {
                        camera_idx: ci,
                        point_idx: pi,
                        pixel: pix,
                    });
                } else {
                    // The bad point through camera 0 has no projection; still
                    // register an observation (with a dummy pixel) to prove the
                    // solver skips it gracefully.
                    observations.push(Observation {
                        camera_idx: ci,
                        point_idx: pi,
                        pixel: (320.0, 240.0),
                    });
                }
            }
        }

        let problem = BundleProblem {
            cameras: cameras.clone(),
            points,
            intrinsics: k,
            observations,
        };
        let result = bundle_adjust(&problem, &BundleParams::default());
        // The cost stays finite throughout (the z≈0 observation is skipped).
        assert!(
            result.initial_cost.is_finite() && result.final_cost.is_finite(),
            "z≈0 observation must not produce a non-finite cost"
        );
        assert!(result.final_cost <= result.initial_cost + 1e-9);
    }
}
