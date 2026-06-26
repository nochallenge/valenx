//! Lens-distortion (Brown–Conrady) camera model — the calibration primitive
//! the pinhole pipeline was missing.
//!
//! Stages 3–5 of this crate model every camera as an **ideal pinhole**:
//! [`crate::project_point`] projects a world point through `K, R, t` with a
//! straight perspective divide and *no lens distortion*. Real lenses bend
//! straight lines — radially (barrel / pincushion) and tangentially (from a
//! sensor not perfectly parallel to the lens) — so before any calibrated
//! geometry step the observed pixels must be related to the ideal pinhole ones
//! through a distortion model. This module adds that missing piece: the
//! standard **Brown–Conrady "plumb-bob"** radial-tangential model used by
//! OpenCV's `cv::calibrateCamera` / `cv::undistortPoints` and by COLMAP's
//! `OPENCV` camera, with coefficients `(k1, k2, p1, p2, k3)`.
//!
//! It is a clean-room implementation of the published model (Brown 1966;
//! Heikkilä & Silvén 1997), expressed directly in the crate's existing
//! [`CameraIntrinsics`] + `nalgebra` types — **no** OpenCV / COLMAP source and
//! no new heavy dependency. The two directions are:
//!
//! 1. **Distort (forward)** — analytic. Given a normalized pinhole point
//!    `(x, y) = (X/Z, Y/Z)` on the `Z = 1` plane, apply
//!    ```text
//!    r² = x² + y²
//!    radial    = 1 + k1·r² + k2·r⁴ + k3·r⁶
//!    x_d = x·radial + 2·p1·x·y + p2·(r² + 2·x²)
//!    y_d = y·radial + p1·(r² + 2·y²) + 2·p2·x·y
//!    ```
//!    then map to pixels with `K`. This is [`distort_normalized`] /
//!    [`project_distorted`].
//!
//! 2. **Undistort (inverse)** — there is no closed form, so it is solved by
//!    **fixed-point iteration** (the algorithm in OpenCV's `undistortPoints`):
//!    start from the distorted point as the guess for the undistorted one and
//!    repeatedly subtract the tangential terms and divide out the radial
//!    factor until it converges. This is [`undistort_normalized`]; mapping a
//!    *pixel* straight to a calibrated, undistorted ray is
//!    [`unproject_to_ray`].
//!
//! ## Honesty / scope
//!
//! - Research/educational grade, matching the rest of the crate. The model is
//!   the textbook 5-coefficient plumb-bob model; it does **not** include the
//!   thin-prism or rational (`k4..k6`) terms of OpenCV's extended model, and
//!   the coefficients themselves are assumed *given* (from a prior
//!   calibration) — this module models distortion, it does not *estimate* it.
//! - The inverse is iterative; for extreme distortion at the far image corners
//!   the fixed-point map can converge slowly. The iteration count and tolerance
//!   are caller-tunable via [`undistort_normalized`], and the round-trip error
//!   is what the tests pin down.

use crate::twoview::CameraIntrinsics;
use nalgebra::Vector3;

/// Brown–Conrady radial-tangential lens-distortion coefficients.
///
/// These are the five coefficients of the standard "plumb-bob" model
/// (OpenCV's `(k1, k2, p1, p2, k3)`, COLMAP's `OPENCV` camera). `k1, k2, k3`
/// are the **radial** terms (the dominant barrel/pincushion bending, even in
/// `r`); `p1, p2` are the **tangential** terms (lens/sensor de-centering).
/// All zero is the ideal, distortion-free pinhole — see [`Distortion::none`],
/// which makes [`project_distorted`] reduce exactly to [`crate::project_point`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Distortion {
    /// First radial coefficient `k1` (multiplies `r²`). The dominant term.
    pub k1: f64,
    /// Second radial coefficient `k2` (multiplies `r⁴`).
    pub k2: f64,
    /// First tangential coefficient `p1` (de-centering, couples `x` and `y`).
    pub p1: f64,
    /// Second tangential coefficient `p2` (de-centering, couples `x` and `y`).
    pub p2: f64,
    /// Third radial coefficient `k3` (multiplies `r⁶`). Often `0` for mild
    /// lenses; significant for wide-angle / fisheye-ish optics.
    pub k3: f64,
}

impl Distortion {
    /// Construct a distortion model from its five coefficients
    /// `(k1, k2, p1, p2, k3)`, in the OpenCV / COLMAP `OPENCV` ordering.
    #[inline]
    #[must_use]
    pub fn new(k1: f64, k2: f64, p1: f64, p2: f64, k3: f64) -> Self {
        Self { k1, k2, p1, p2, k3 }
    }

    /// The zero (ideal pinhole) distortion. With this model every operation in
    /// this module reduces to the undistorted pinhole one: [`distort_normalized`]
    /// is the identity and [`project_distorted`] equals [`crate::project_point`].
    #[inline]
    #[must_use]
    pub fn none() -> Self {
        Self {
            k1: 0.0,
            k2: 0.0,
            p1: 0.0,
            p2: 0.0,
            k3: 0.0,
        }
    }

    /// `true` if every coefficient is exactly zero (the ideal pinhole). Lets
    /// callers skip the distortion arithmetic on calibrated/rectified images.
    #[inline]
    #[must_use]
    pub fn is_none(&self) -> bool {
        self.k1 == 0.0 && self.k2 == 0.0 && self.p1 == 0.0 && self.p2 == 0.0 && self.k3 == 0.0
    }
}

/// Apply the Brown–Conrady model to a **normalized** (ideal pinhole) point.
///
/// `(x, y)` is a point on the `Z = 1` normalized image plane — i.e. the camera
/// ray `(X, Y, Z)` after the perspective divide `(X/Z, Y/Z)`, **before** the
/// intrinsics `K` are applied. The returned `(x_d, y_d)` is the distorted
/// normalized point; feed it through `K` to get the observed pixel (this is
/// exactly what [`project_distorted`] does).
///
/// This is the forward (analytic) direction; the inverse is
/// [`undistort_normalized`].
#[must_use]
pub fn distort_normalized(d: &Distortion, x: f64, y: f64) -> (f64, f64) {
    let r2 = x * x + y * y;
    let r4 = r2 * r2;
    let r6 = r4 * r2;
    let radial = 1.0 + d.k1 * r2 + d.k2 * r4 + d.k3 * r6;
    let xy = x * y;
    let x_d = x * radial + 2.0 * d.p1 * xy + d.p2 * (r2 + 2.0 * x * x);
    let y_d = y * radial + d.p1 * (r2 + 2.0 * y * y) + 2.0 * d.p2 * xy;
    (x_d, y_d)
}

/// Invert the Brown–Conrady model: recover the **undistorted** normalized
/// point from a distorted one by fixed-point iteration.
///
/// Given a distorted normalized point `(x_d, y_d)` (a pixel mapped back through
/// `K⁻¹`), find the ideal-pinhole `(x, y)` such that
/// `distort_normalized(d, x, y) == (x_d, y_d)`. There is no closed form, so
/// this runs the standard OpenCV `undistortPoints` fixed-point loop: seed
/// `(x, y) = (x_d, y_d)`, then repeatedly remove the tangential terms and
/// divide out the radial factor.
///
/// `max_iters` bounds the iteration count and `tol` is the early-exit
/// convergence threshold on the normalized-plane step size (`0.0` to always run
/// the full `max_iters`). For mild-to-moderate distortion this converges in a
/// handful of iterations; OpenCV's default is 5. Returns the converged
/// `(x, y)`.
#[must_use]
pub fn undistort_normalized(
    d: &Distortion,
    x_d: f64,
    y_d: f64,
    max_iters: usize,
    tol: f64,
) -> (f64, f64) {
    // Fast path: nothing to undo.
    if d.is_none() {
        return (x_d, y_d);
    }
    let (mut x, mut y) = (x_d, y_d);
    for _ in 0..max_iters {
        let r2 = x * x + y * y;
        let r4 = r2 * r2;
        let r6 = r4 * r2;
        let radial = 1.0 + d.k1 * r2 + d.k2 * r4 + d.k3 * r6;
        // Guard a (pathological, strongly negative-distortion) zero radial.
        if radial.abs() < 1e-12 {
            break;
        }
        let xy = x * y;
        let dx_tan = 2.0 * d.p1 * xy + d.p2 * (r2 + 2.0 * x * x);
        let dy_tan = d.p1 * (r2 + 2.0 * y * y) + 2.0 * d.p2 * xy;
        let x_new = (x_d - dx_tan) / radial;
        let y_new = (y_d - dy_tan) / radial;
        let step = (x_new - x).hypot(y_new - y);
        x = x_new;
        y = y_new;
        if step <= tol {
            break;
        }
    }
    (x, y)
}

/// Project a world point through `K, R, t` **with** lens distortion to a pixel.
///
/// The distorted counterpart of [`crate::project_point`]: it computes the
/// camera-frame point `x_cam = R·X + t`, perspective-divides to the normalized
/// point `(x_cam.x/z, x_cam.y/z)`, applies the Brown–Conrady model
/// ([`distort_normalized`]), and finally maps through the intrinsics `K` to the
/// observed pixel `(fx·x_d + cx, fy·y_d + cy)`.
///
/// Returns [`None`] when the camera-frame depth `z` is at (or near) zero — the
/// point lies on the principal plane and has no finite image — matching
/// [`crate::project_point`]'s `z ≈ 0` guard. With [`Distortion::none`] the
/// result is bit-for-bit [`crate::project_point`].
#[must_use]
pub fn project_distorted(
    k: &CameraIntrinsics,
    d: &Distortion,
    r: &nalgebra::Matrix3<f64>,
    t: &Vector3<f64>,
    x: &Vector3<f64>,
) -> Option<(f64, f64)> {
    let cam = r * x + t;
    if cam.z.abs() < 1e-12 {
        return None;
    }
    let (x_d, y_d) = distort_normalized(d, cam.x / cam.z, cam.y / cam.z);
    Some((k.fx * x_d + k.cx, k.fy * y_d + k.cy))
}

/// Back-project an observed (distorted) **pixel** to a unit-`z` calibrated
/// camera ray, undistorting in the process.
///
/// The inverse of [`project_distorted`]'s last two steps: map the pixel
/// `(u, v)` back through `K⁻¹` to the distorted normalized point
/// `((u−cx)/fx, (v−cy)/fy)`, undistort it ([`undistort_normalized`]), and
/// return the homogeneous ray `(x, y, 1)`. Multiplying this ray by any positive
/// depth gives a camera-frame point that re-projects to `(u, v)`.
///
/// This is the calibrated, distortion-aware unprojection the matcher/PnP front
/// end needs to turn raw image observations into the ideal rays that the
/// pinhole geometry stages assume. `max_iters` / `tol` tune the inverse
/// iteration exactly as in [`undistort_normalized`].
#[must_use]
pub fn unproject_to_ray(
    k: &CameraIntrinsics,
    d: &Distortion,
    u: f64,
    v: f64,
    max_iters: usize,
    tol: f64,
) -> Vector3<f64> {
    let x_d = (u - k.cx) / k.fx;
    let y_d = (v - k.cy) / k.fy;
    let (x, y) = undistort_normalized(d, x_d, y_d, max_iters, tol);
    Vector3::new(x, y, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pnp::project_point;
    use nalgebra::Matrix3;

    fn k() -> CameraIntrinsics {
        CameraIntrinsics::new(800.0, 800.0, 320.0, 240.0)
    }

    /// With zero distortion the distorted projection must equal the plain
    /// pinhole projection exactly (the model degenerates correctly).
    #[test]
    fn zero_distortion_equals_pinhole() {
        let k = k();
        let d = Distortion::none();
        let r = Matrix3::identity();
        let t = Vector3::new(0.1, -0.2, 5.0);
        for x in [
            Vector3::new(0.0, 0.0, 4.0),
            Vector3::new(1.0, 0.5, 6.0),
            Vector3::new(-0.8, 1.2, 3.0),
        ] {
            let a = project_distorted(&k, &d, &r, &t, &x).unwrap();
            let b = project_point(&k, &r, &t, &x).unwrap();
            assert!((a.0 - b.0).abs() < 1e-12 && (a.1 - b.1).abs() < 1e-12);
        }
    }

    /// Known-point projection: a point on the optical axis (normalized
    /// `(0, 0)`) is a distortion fixed point and must land exactly on the
    /// principal point `(cx, cy)` regardless of the coefficients.
    #[test]
    fn principal_axis_is_distortion_fixed_point() {
        let k = k();
        let d = Distortion::new(-0.28, 0.10, 0.001, -0.002, 0.05);
        let r = Matrix3::identity();
        let t = Vector3::zeros();
        // World point straight down the +Z axis projects to (cx, cy).
        let p = project_distorted(&k, &d, &r, &t, &Vector3::new(0.0, 0.0, 2.0)).unwrap();
        assert!((p.0 - k.cx).abs() < 1e-9 && (p.1 - k.cy).abs() < 1e-9);
    }

    /// Hand-checked forward value: with only `k1`, a normalized point distorts
    /// by exactly the radial factor `1 + k1·r²`, then maps through `K`.
    #[test]
    fn forward_radial_matches_hand_value() {
        let k = k();
        let d = Distortion::new(0.1, 0.0, 0.0, 0.0, 0.0);
        let (x, y) = (0.2, 0.1); // normalized point
        let r2 = x * x + y * y; // 0.05
        let radial = 1.0 + 0.1 * r2; // 1.005
        let (xd, yd) = distort_normalized(&d, x, y);
        assert!((xd - x * radial).abs() < 1e-15);
        assert!((yd - y * radial).abs() < 1e-15);
    }

    /// Round-trip: distort a grid of normalized points, then undistort, and
    /// recover the originals to high precision. This is the core inverse test.
    #[test]
    fn distort_undistort_round_trip() {
        let d = Distortion::new(-0.25, 0.08, 0.0008, -0.0011, 0.03);
        let mut max_err: f64 = 0.0;
        for i in -4..=4 {
            for j in -4..=4 {
                let x = i as f64 * 0.12;
                let y = j as f64 * 0.12;
                let (xd, yd) = distort_normalized(&d, x, y);
                let (xu, yu) = undistort_normalized(&d, xd, yd, 20, 1e-12);
                max_err = max_err.max((xu - x).hypot(yu - y));
            }
        }
        assert!(max_err < 1e-9, "round-trip error {max_err} too large");
    }

    /// Full pixel round-trip through the camera model: project a 3-D point with
    /// distortion to a pixel, unproject that pixel back to a ray, and confirm
    /// the ray's direction matches the original point's direction.
    #[test]
    fn pixel_unproject_round_trip() {
        let k = k();
        let d = Distortion::new(-0.20, 0.05, 0.001, 0.0007, 0.0);
        let r = Matrix3::identity();
        let t = Vector3::zeros();
        let world = Vector3::new(0.7, -0.4, 3.0);
        let (u, v) = project_distorted(&k, &d, &r, &t, &world).unwrap();
        let ray = unproject_to_ray(&k, &d, u, v, 20, 1e-12);
        // The ray (x, y, 1) should be parallel to (world.x, world.y, world.z),
        // i.e. equal after scaling both to z = 1.
        let expected = Vector3::new(world.x / world.z, world.y / world.z, 1.0);
        assert!(
            (ray - expected).norm() < 1e-7,
            "ray {ray:?} vs {expected:?}"
        );
    }
}
