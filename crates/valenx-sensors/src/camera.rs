//! Simulated **pinhole camera** with Brown–Conrady lens distortion.
//!
//! A [`Camera`] projects a 3-D point given in the **camera frame** to a 2-D
//! pixel, following the standard computer-vision pipeline:
//!
//! 1. **Perspective divide** — a camera-frame point `(X, Y, Z)` with `Z > 0`
//!    (in front of the camera, looking down +Z) maps to the *normalised image
//!    coordinates* `(x, y) = (X/Z, Y/Z)`.
//! 2. **Brown–Conrady distortion** — radial (`k1, k2, k3`) and tangential
//!    (`p1, p2`) lens distortion warps `(x, y)` to `(x_d, y_d)`.
//! 3. **Intrinsics `K`** — focal lengths `(fx, fy)` and principal point
//!    `(cx, cy)` (plus an optional skew) map the distorted normalised
//!    coordinates to pixels: `u = fx·x_d + skew·y_d + cx`, `v = fy·y_d + cy`.
//!
//! [`Camera::undistort`] inverts step 2 (the radial-tangential model has no
//! closed-form inverse) by fixed-point iteration, so a pixel can be lifted back
//! to a normalised ray; [`Camera::project`] and [`Camera::undistort`] round-trip
//! to a tight tolerance, which is what the tests pin.
//!
//! Frame convention: **+Z forward (optical axis), +X right, +Y down** — the
//! usual OpenCV image convention, so `u` grows right and `v` grows down. This is
//! distinct from the +x-forward body frame used elsewhere in the crate; the
//! autonomy harness applies the camera's mount rotation to bridge the two.

use crate::error::SensorError;

/// Pinhole intrinsic matrix `K` parameters (pixels).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Intrinsics {
    /// Focal length in the x direction (pixels, > 0).
    pub fx: f64,
    /// Focal length in the y direction (pixels, > 0).
    pub fy: f64,
    /// Principal point x (pixels).
    pub cx: f64,
    /// Principal point y (pixels).
    pub cy: f64,
    /// Axis skew (pixels); `0` for square pixels with no shear (the common case).
    pub skew: f64,
}

impl Intrinsics {
    /// A skew-free intrinsic set.
    #[must_use]
    pub fn new(fx: f64, fy: f64, cx: f64, cy: f64) -> Self {
        Self {
            fx,
            fy,
            cx,
            cy,
            skew: 0.0,
        }
    }
}

/// Brown–Conrady radial-tangential distortion coefficients.
///
/// The distorted normalised coordinates are
/// ```text
/// r²   = x² + y²
/// rad  = 1 + k1·r² + k2·r⁴ + k3·r⁶
/// x_d  = x·rad + 2·p1·x·y + p2·(r² + 2x²)
/// y_d  = y·rad + p1·(r² + 2y²) + 2·p2·x·y
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Distortion {
    /// 2nd-order radial coefficient.
    pub k1: f64,
    /// 4th-order radial coefficient.
    pub k2: f64,
    /// 6th-order radial coefficient.
    pub k3: f64,
    /// First tangential coefficient.
    pub p1: f64,
    /// Second tangential coefficient.
    pub p2: f64,
}

impl Distortion {
    /// Coefficients for an ideal (distortion-free) lens — all zero.
    #[must_use]
    pub fn none() -> Self {
        Self::default()
    }

    /// Radial-only distortion `(k1, k2, k3)` with no tangential term.
    #[must_use]
    pub fn radial(k1: f64, k2: f64, k3: f64) -> Self {
        Self {
            k1,
            k2,
            k3,
            p1: 0.0,
            p2: 0.0,
        }
    }

    /// Apply the Brown–Conrady model to a normalised image point `(x, y)`,
    /// returning the distorted `(x_d, y_d)`.
    #[must_use]
    pub fn apply(&self, x: f64, y: f64) -> (f64, f64) {
        let r2 = x * x + y * y;
        let radial = 1.0 + self.k1 * r2 + self.k2 * r2 * r2 + self.k3 * r2 * r2 * r2;
        let x_d = x * radial + 2.0 * self.p1 * x * y + self.p2 * (r2 + 2.0 * x * x);
        let y_d = y * radial + self.p1 * (r2 + 2.0 * y * y) + 2.0 * self.p2 * x * y;
        (x_d, y_d)
    }

    /// Invert the distortion: given distorted `(x_d, y_d)`, recover the
    /// undistorted normalised `(x, y)` by fixed-point iteration (the OpenCV
    /// scheme). Converges for the small-to-moderate distortion typical of real
    /// lenses; returns the best iterate after `max_iter` steps.
    #[must_use]
    pub fn invert(&self, x_d: f64, y_d: f64, max_iter: usize, tol: f64) -> (f64, f64) {
        // Seed with the distorted point itself.
        let (mut x, mut y) = (x_d, y_d);
        for _ in 0..max_iter {
            let r2 = x * x + y * y;
            let radial = 1.0 + self.k1 * r2 + self.k2 * r2 * r2 + self.k3 * r2 * r2 * r2;
            // Tangential components at the current estimate.
            let dx = 2.0 * self.p1 * x * y + self.p2 * (r2 + 2.0 * x * x);
            let dy = self.p1 * (r2 + 2.0 * y * y) + 2.0 * self.p2 * x * y;
            // Solve x·radial + dx = x_d ⇒ x = (x_d − dx) / radial.
            let x_new = (x_d - dx) / radial;
            let y_new = (y_d - dy) / radial;
            let delta = (x_new - x).hypot(y_new - y);
            x = x_new;
            y = y_new;
            if delta < tol {
                break;
            }
        }
        (x, y)
    }
}

/// A pixel coordinate (sub-pixel precision).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pixel {
    /// Horizontal pixel coordinate (grows right).
    pub u: f64,
    /// Vertical pixel coordinate (grows down).
    pub v: f64,
}

/// A simulated pinhole camera with lens distortion and a finite sensor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    /// Intrinsic matrix parameters.
    pub intrinsics: Intrinsics,
    /// Lens distortion coefficients.
    pub distortion: Distortion,
    /// Image width (pixels).
    pub width: usize,
    /// Image height (pixels).
    pub height: usize,
}

impl Camera {
    /// Build a camera, validating intrinsics and image size.
    ///
    /// # Errors
    /// Returns [`SensorError::InvalidConfig`] if either focal length is not
    /// strictly positive / finite, any intrinsic is non-finite, or the image has
    /// zero width or height.
    pub fn new(
        intrinsics: Intrinsics,
        distortion: Distortion,
        width: usize,
        height: usize,
    ) -> Result<Self, SensorError> {
        if !(intrinsics.fx.is_finite() && intrinsics.fx > 0.0) {
            return Err(SensorError::InvalidConfig(format!(
                "fx must be finite and > 0, got {}",
                intrinsics.fx
            )));
        }
        if !(intrinsics.fy.is_finite() && intrinsics.fy > 0.0) {
            return Err(SensorError::InvalidConfig(format!(
                "fy must be finite and > 0, got {}",
                intrinsics.fy
            )));
        }
        if !(intrinsics.cx.is_finite() && intrinsics.cy.is_finite() && intrinsics.skew.is_finite())
        {
            return Err(SensorError::InvalidConfig(
                "principal point / skew must be finite".into(),
            ));
        }
        if width == 0 || height == 0 {
            return Err(SensorError::InvalidConfig(
                "image width and height must be ≥ 1".into(),
            ));
        }
        Ok(Self {
            intrinsics,
            distortion,
            width,
            height,
        })
    }

    /// Project a **camera-frame** 3-D point `(x, y, z)` (z forward) to a pixel.
    ///
    /// Returns `None` if the point is at or behind the image plane (`z ≤ 0`), or
    /// if any coordinate is non-finite. The returned pixel may fall outside the
    /// `[0, width) × [0, height)` image bounds (the point is in front of the
    /// camera but outside the field of view); use [`Camera::in_bounds`] to test
    /// that.
    #[must_use]
    pub fn project(&self, x: f64, y: f64, z: f64) -> Option<Pixel> {
        if !(x.is_finite() && y.is_finite() && z.is_finite()) || z <= 0.0 {
            return None;
        }
        let xn = x / z;
        let yn = y / z;
        let (xd, yd) = self.distortion.apply(xn, yn);
        let k = &self.intrinsics;
        Some(Pixel {
            u: k.fx * xd + k.skew * yd + k.cx,
            v: k.fy * yd + k.cy,
        })
    }

    /// Whether a pixel lies inside the image bounds `[0, width) × [0, height)`.
    #[must_use]
    pub fn in_bounds(&self, p: &Pixel) -> bool {
        p.u >= 0.0 && p.u < self.width as f64 && p.v >= 0.0 && p.v < self.height as f64
    }

    /// Lift a pixel back to an **undistorted normalised** image point
    /// `(x, y) = (X/Z, Y/Z)` — the inverse of intrinsics + distortion.
    ///
    /// Applies `K⁻¹` then inverts the Brown–Conrady distortion by iteration.
    /// `project(undistort(p))·Z` recovers `p` to within the iteration tolerance.
    #[must_use]
    pub fn undistort(&self, p: &Pixel) -> (f64, f64) {
        let k = &self.intrinsics;
        // Invert K (upper-triangular with skew): solve for distorted normalised.
        let yd = (p.v - k.cy) / k.fy;
        let xd = (p.u - k.cx - k.skew * yd) / k.fx;
        self.distortion.invert(xd, yd, 50, 1e-12)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vga() -> Camera {
        // 640×480, fx=fy=500, principal point at the image centre.
        Camera::new(
            Intrinsics::new(500.0, 500.0, 320.0, 240.0),
            Distortion::none(),
            640,
            480,
        )
        .unwrap()
    }

    #[test]
    fn principal_axis_point_maps_to_principal_point() {
        // A point on the optical axis projects to the principal point regardless
        // of depth.
        let cam = vga();
        let p = cam.project(0.0, 0.0, 5.0).unwrap();
        assert!((p.u - 320.0).abs() < 1e-12);
        assert!((p.v - 240.0).abs() < 1e-12);
    }

    #[test]
    fn known_point_projects_to_exact_pixel() {
        // (X,Y,Z) = (1,2,10), fx=fy=500: x_n=0.1, y_n=0.2 ⇒
        // u = 500·0.1 + 320 = 370; v = 500·0.2 + 240 = 340.
        let cam = vga();
        let p = cam.project(1.0, 2.0, 10.0).unwrap();
        assert!((p.u - 370.0).abs() < 1e-12, "u = {}", p.u);
        assert!((p.v - 340.0).abs() < 1e-12, "v = {}", p.v);
    }

    #[test]
    fn behind_camera_and_non_finite_do_not_project() {
        let cam = vga();
        assert!(cam.project(1.0, 1.0, 0.0).is_none());
        assert!(cam.project(1.0, 1.0, -2.0).is_none());
        assert!(cam.project(f64::NAN, 0.0, 5.0).is_none());
    }

    #[test]
    fn undistort_inverts_intrinsics_with_no_distortion() {
        // With no distortion, undistort is just K⁻¹; project∘undistort = identity
        // up to the perspective scale.
        let cam = vga();
        let p = Pixel { u: 400.0, v: 300.0 };
        let (x, y) = cam.undistort(&p);
        // Re-project at depth 1: pixel should match p.
        let back = cam.project(x, y, 1.0).unwrap();
        assert!((back.u - p.u).abs() < 1e-9, "u {} vs {}", back.u, p.u);
        assert!((back.v - p.v).abs() < 1e-9, "v {} vs {}", back.v, p.v);
    }

    #[test]
    fn project_undistort_roundtrip_with_radial_tangential_distortion() {
        // A realistic wide-lens distortion; round-trip several off-centre points.
        let cam = Camera::new(
            Intrinsics::new(480.0, 480.0, 320.0, 240.0),
            Distortion {
                k1: -0.28,
                k2: 0.10,
                k3: -0.02,
                p1: 0.001,
                p2: -0.0015,
            },
            640,
            480,
        )
        .unwrap();

        for &(x, y, z) in &[
            (0.30, 0.20, 1.0),
            (-0.25, 0.15, 1.0),
            (0.40, -0.35, 1.0),
            (-0.10, -0.45, 1.0),
        ] {
            let p = cam.project(x, y, z).unwrap();
            let (xu, yu) = cam.undistort(&p);
            // Undistort recovers the normalised point (since z = 1, x_n = x).
            assert!((xu - x).abs() < 1e-7, "x {xu} vs {x}");
            assert!((yu - y).abs() < 1e-7, "y {yu} vs {y}");
        }
    }

    #[test]
    fn distortion_actually_shifts_off_axis_points() {
        // With negative k1 (barrel), an off-axis point moves toward the centre
        // relative to the ideal pinhole pixel.
        let distorted = Camera::new(
            Intrinsics::new(500.0, 500.0, 320.0, 240.0),
            Distortion::radial(-0.2, 0.0, 0.0),
            640,
            480,
        )
        .unwrap();
        let ideal = vga();
        let pd = distorted.project(0.4, 0.0, 1.0).unwrap();
        let pi = ideal.project(0.4, 0.0, 1.0).unwrap();
        assert!(pd.u != pi.u, "distortion should move the pixel");
        // Barrel distortion pulls the point inward (toward cx = 320).
        assert!(pd.u < pi.u, "barrel distortion should pull inward");
    }

    #[test]
    fn in_bounds_flags_outside_pixels() {
        let cam = vga();
        assert!(cam.in_bounds(&cam.project(0.0, 0.0, 5.0).unwrap()));
        // A point far off to the side projects outside the 640-wide image.
        let off = cam.project(5.0, 0.0, 1.0).unwrap();
        assert!(
            !cam.in_bounds(&off),
            "u = {} should be out of bounds",
            off.u
        );
    }

    #[test]
    fn invalid_intrinsics_rejected() {
        assert!(Camera::new(
            Intrinsics::new(0.0, 500.0, 0.0, 0.0),
            Distortion::none(),
            64,
            48
        )
        .is_err());
        assert!(Camera::new(
            Intrinsics::new(500.0, -1.0, 0.0, 0.0),
            Distortion::none(),
            64,
            48
        )
        .is_err());
        assert!(Camera::new(
            Intrinsics::new(500.0, 500.0, 0.0, 0.0),
            Distortion::none(),
            0,
            48
        )
        .is_err());
    }
}
