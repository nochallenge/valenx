//! Paraxial ray-transfer (ABCD) matrices and sequential ray tracing.
//!
//! ## What
//!
//! The ray-transfer-matrix (a.k.a. ABCD-matrix) formalism is the standard
//! first-order tool of geometric optics. A paraxial ray is described at a
//! reference plane by a two-component state — its transverse height `y`
//! from the optical axis and its (paraxial) slope `theta` — and every
//! linear optical element acts on that state as a 2x2 matrix:
//!
//! ```text
//! [ y_out     ]   [ A  B ] [ y_in     ]
//! [ theta_out ] = [ C  D ] [ theta_in ]
//! ```
//!
//! Cascading elements is matrix multiplication (in **reverse** order of
//! travel, like all operator products), so an entire optical train
//! collapses to a single [`RayTransfer`]. This module builds the three
//! elementary matrices — a [free-space gap](RayTransfer::propagation), a
//! [thin lens](RayTransfer::thin_lens), and a [flat or spherical
//! refracting surface](RayTransfer::refraction) — lets you
//! [compose](RayTransfer::then) and [apply](RayTransfer::apply) them, and
//! extracts the system's cardinal quantities: its
//! [effective focal length](RayTransfer::effective_focal_length) and the
//! [object→image relationship](RayTransfer::image_distance) (image
//! distance and transverse magnification) for a given object distance.
//!
//! ## Model
//!
//! The elementary matrices (lens / surface in air unless an index is
//! given) are:
//!
//! ```text
//! propagation over distance d : [ 1  d ]
//!                               [ 0  1 ]
//!
//! thin lens of focal length f : [  1    0 ]
//!                               [ -1/f  1 ]
//!
//! refraction at a spherical    [    1        0    ]
//! surface, radius R, n1 -> n2: [ (n1-n2)/(n2 R)  n1/n2 ]
//! ```
//!
//! The thin-lens matrix is exactly the `R -> infinity, R -> -infinity`
//! composition of two refractions wrapping zero thickness, and reproduces
//! `1/f = 1/do + 1/di` (see the tests). A single refracting surface has
//! power `P = (n2 - n1) / R`, i.e. front focal length `f = n1 / P` and
//! back focal length `n2 / P`; this module reports the back focal length,
//! which equals `1/f` of the thin-lens convention when `n1 = n2 = 1`.
//!
//! ## Sign convention
//!
//! Consistent with [`crate::thin_lens`] and [`crate::lensmaker`]: a
//! positive focal length / positive surface power converges; distances are
//! measured positive in the direction of propagation; a surface radius `R`
//! is positive when its centre of curvature lies downstream (on the
//! transmitted side). Slopes are paraxial (`tan theta ~= theta`); the
//! reduced angle `n * theta` is *not* used here — plain geometric slope is
//! carried, which is the common single-medium textbook convention and
//! makes [`Ray`] directly comparable across air gaps.
//!
//! ## Honest scope
//!
//! Strictly first-order (paraxial): every element is linearised about the
//! axis, so this captures focal lengths, image locations, magnifications
//! and system layout — but **no** aberrations (spherical, coma,
//! astigmatism, field curvature, distortion), no chromatic effects (`n` is
//! a single number, not a dispersion curve), no apertures / vignetting /
//! finite stops, and no diffraction. It is the matrix companion to the
//! closed-form thin-lens and lensmaker models, not a real-ray sequential
//! tracer. Research / educational grade — not an optical-design tool.

use crate::error::{validate_index, OpticsError};

/// A paraxial ray at a reference plane: transverse height and slope.
///
/// `height` is the perpendicular distance from the optical axis (same
/// length unit as the system) and `angle` is the paraxial slope
/// `dy/dz ~= tan(theta) ~= theta` (radians, small). Construct directly —
/// both fields are public — or with [`Ray::new`], which rejects
/// non-finite components.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Ray {
    /// Transverse height from the optical axis.
    pub height: f64,
    /// Paraxial slope `dy/dz` (radians; small-angle).
    pub angle: f64,
}

impl Ray {
    /// Construct a ray from a height and a paraxial slope.
    ///
    /// # Errors
    ///
    /// Returns [`OpticsError::InvalidParameter`] if either component is
    /// not finite.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::abcd::Ray;
    /// // A ray parallel to the axis, 5 mm above it.
    /// let r = Ray::new(5.0, 0.0).unwrap();
    /// assert_eq!(r.angle, 0.0);
    /// assert!(Ray::new(f64::NAN, 0.0).is_err());
    /// ```
    pub fn new(height: f64, angle: f64) -> Result<Self, OpticsError> {
        if !height.is_finite() {
            return Err(OpticsError::InvalidParameter {
                name: "height",
                value: height,
                reason: "ray height must be finite",
            });
        }
        if !angle.is_finite() {
            return Err(OpticsError::InvalidParameter {
                name: "angle",
                value: angle,
                reason: "ray slope must be finite",
            });
        }
        Ok(Self { height, angle })
    }
}

/// A 2x2 paraxial ray-transfer (ABCD) matrix.
///
/// Stored row-major as `[[a, b], [c, d]]`. Build the elementary elements
/// with [`RayTransfer::propagation`], [`RayTransfer::thin_lens`] and
/// [`RayTransfer::refraction`]; cascade them with [`RayTransfer::then`];
/// push a [`Ray`] through with [`RayTransfer::apply`]; and read off the
/// cardinal quantities with [`RayTransfer::effective_focal_length`] and
/// [`RayTransfer::image_distance`].
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RayTransfer {
    /// `A` — height response to input height.
    pub a: f64,
    /// `B` — height response to input slope (units of length).
    pub b: f64,
    /// `C` — slope response to input height (units of 1/length).
    pub c: f64,
    /// `D` — slope response to input slope.
    pub d: f64,
}

impl RayTransfer {
    /// The 2x2 identity (a zero-length, no-op element).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::abcd::{Ray, RayTransfer};
    /// let r = Ray::new(2.0, 0.1).unwrap();
    /// assert_eq!(RayTransfer::identity().apply(r), r);
    /// ```
    pub fn identity() -> Self {
        Self {
            a: 1.0,
            b: 0.0,
            c: 0.0,
            d: 1.0,
        }
    }

    /// Free-space (or uniform-medium) propagation over a distance `d`.
    ///
    /// The matrix is `[[1, d], [0, 1]]`: the slope is unchanged and the
    /// height grows by `d * angle`.
    ///
    /// # Errors
    ///
    /// Returns [`OpticsError::InvalidParameter`] if `d` is not finite.
    /// Negative distances are permitted (they refer a plane backwards,
    /// e.g. to place an object in front of the first element).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::abcd::{Ray, RayTransfer};
    /// // A ray sloping up at 0.1 rad climbs 0.1*10 = 1.0 over 10 units.
    /// let out = RayTransfer::propagation(10.0).unwrap().apply(Ray::new(0.0, 0.1).unwrap());
    /// assert!((out.height - 1.0).abs() < 1e-12);
    /// assert!((out.angle - 0.1).abs() < 1e-12);
    /// ```
    pub fn propagation(d: f64) -> Result<Self, OpticsError> {
        if !d.is_finite() {
            return Err(OpticsError::InvalidParameter {
                name: "d",
                value: d,
                reason: "propagation distance must be finite",
            });
        }
        Ok(Self {
            a: 1.0,
            b: d,
            c: 0.0,
            d: 1.0,
        })
    }

    /// A thin lens of focal length `f` (matrix `[[1, 0], [-1/f, 1]]`).
    ///
    /// The height is unchanged at the lens plane; the slope is bent by
    /// `-y / f`, so a ray parallel to the axis (`angle = 0`) leaves with
    /// slope `-y / f` and therefore crosses the axis one focal length
    /// downstream — the defining property of the focal point.
    ///
    /// # Errors
    ///
    /// Returns [`OpticsError::InvalidParameter`] if `f` is not finite, or
    /// [`OpticsError::SingularGeometry`] if `f` is exactly zero (infinite
    /// power).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::abcd::{Ray, RayTransfer};
    /// // A parallel input ray leaves a +f lens aimed to cross the axis at f.
    /// let lens = RayTransfer::thin_lens(50.0).unwrap();
    /// let out = lens.apply(Ray::new(5.0, 0.0).unwrap());
    /// assert!((out.angle + 5.0 / 50.0).abs() < 1e-12); // slope = -y/f
    /// ```
    pub fn thin_lens(f: f64) -> Result<Self, OpticsError> {
        if !f.is_finite() {
            return Err(OpticsError::InvalidParameter {
                name: "f",
                value: f,
                reason: "focal length must be finite",
            });
        }
        if f == 0.0 {
            return Err(OpticsError::SingularGeometry {
                reason: "focal length of zero has undefined (infinite) optical power",
            });
        }
        Ok(Self {
            a: 1.0,
            b: 0.0,
            c: -1.0 / f,
            d: 1.0,
        })
    }

    /// Refraction at a single spherical surface of radius `r` separating a
    /// medium of index `n1` (incident) from `n2` (transmitted).
    ///
    /// The matrix is `[[1, 0], [(n1 - n2) / (n2 r), n1 / n2]]`, the
    /// surface power being `P = (n2 - n1) / r`. Pass [`f64::INFINITY`] for
    /// `r` to model a **flat** interface (a plane refraction: zero power,
    /// pure slope rescaling by `n1 / n2`).
    ///
    /// # Errors
    ///
    /// - [`OpticsError::InvalidParameter`] if either index is non-finite
    ///   or non-positive, or if `r` is `NaN` or exactly zero (a degenerate
    ///   point surface).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::abcd::RayTransfer;
    /// // A flat air->glass interface only rescales the slope by n1/n2.
    /// let flat = RayTransfer::refraction(1.0, 1.5, f64::INFINITY).unwrap();
    /// assert!((flat.d - 1.0 / 1.5).abs() < 1e-12);
    /// assert!(flat.c.abs() < 1e-15); // zero power
    /// ```
    pub fn refraction(n1: f64, n2: f64, r: f64) -> Result<Self, OpticsError> {
        let n1 = validate_index("n1", n1)?;
        let n2 = validate_index("n2", n2)?;
        if r.is_nan() {
            return Err(OpticsError::InvalidParameter {
                name: "r",
                value: r,
                reason: "surface radius must not be NaN",
            });
        }
        if r == 0.0 {
            return Err(OpticsError::InvalidParameter {
                name: "r",
                value: r,
                reason: "surface radius of zero is a degenerate (point) surface",
            });
        }
        // 1/R is zero for a flat surface, killing the power term.
        let inv_r = if r.is_infinite() { 0.0 } else { 1.0 / r };
        Ok(Self {
            a: 1.0,
            b: 0.0,
            c: (n1 - n2) / n2 * inv_r,
            d: n1 / n2,
        })
    }

    /// Compose with a **subsequent** element: the ray passes through
    /// `self` first, then `next`.
    ///
    /// This returns the matrix product `next * self` (operators compose
    /// right-to-left, so the element the ray meets first sits on the
    /// right). Building a system left-to-right in travel order is then
    /// simply `first.then(second).then(third)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::abcd::RayTransfer;
    /// // Gap, then lens, then gap — a single combined system matrix.
    /// let sys = RayTransfer::propagation(20.0).unwrap()
    ///     .then(&RayTransfer::thin_lens(50.0).unwrap())
    ///     .then(&RayTransfer::propagation(30.0).unwrap());
    /// assert!(sys.determinant().is_finite());
    /// ```
    pub fn then(&self, next: &RayTransfer) -> RayTransfer {
        // M_total = next * self.
        RayTransfer {
            a: next.a * self.a + next.b * self.c,
            b: next.a * self.b + next.b * self.d,
            c: next.c * self.a + next.d * self.c,
            d: next.c * self.b + next.d * self.d,
        }
    }

    /// The determinant `A*D - B*C`.
    ///
    /// For a system embedded in a single medium (or any all-in-air train)
    /// this is exactly `1`; for a system that starts in `n_in` and ends in
    /// `n_out` it equals `n_in / n_out`. It is a cheap self-consistency
    /// check on a composed matrix.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::abcd::RayTransfer;
    /// // Any all-in-air system has unit determinant.
    /// let sys = RayTransfer::propagation(10.0).unwrap()
    ///     .then(&RayTransfer::thin_lens(25.0).unwrap());
    /// assert!((sys.determinant() - 1.0).abs() < 1e-12);
    /// ```
    pub fn determinant(&self) -> f64 {
        self.a * self.d - self.b * self.c
    }

    /// Push a [`Ray`] through this element / system.
    ///
    /// Returns the transformed ray `[A B; C D] * [height; angle]`.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::abcd::{Ray, RayTransfer};
    /// let out = RayTransfer::propagation(5.0).unwrap().apply(Ray::new(1.0, 0.2).unwrap());
    /// assert!((out.height - (1.0 + 5.0 * 0.2)).abs() < 1e-12);
    /// ```
    pub fn apply(&self, ray: Ray) -> Ray {
        Ray {
            height: self.a * ray.height + self.b * ray.angle,
            angle: self.c * ray.height + self.d * ray.angle,
        }
    }

    /// The effective (back) focal length of this system, `f = -1 / C`.
    ///
    /// The `C` element is `-1/f` (the negative of the system's optical
    /// power), so the effective focal length is its negative reciprocal.
    /// A converging system has `C < 0` and hence `f > 0`.
    ///
    /// # Errors
    ///
    /// Returns [`OpticsError::SingularGeometry`] if `C == 0` — an
    /// afocal (telescopic / zero-power) system, whose focal length is
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::abcd::RayTransfer;
    /// // A bare thin lens reports its own focal length back.
    /// let f = RayTransfer::thin_lens(75.0).unwrap().effective_focal_length().unwrap();
    /// assert!((f - 75.0).abs() < 1e-12);
    /// ```
    pub fn effective_focal_length(&self) -> Result<f64, OpticsError> {
        if self.c == 0.0 {
            return Err(OpticsError::SingularGeometry {
                reason: "afocal system (C = 0): the effective focal length is infinite",
            });
        }
        Ok(-1.0 / self.c)
    }

    /// Locate the image of an on-axis object a distance `object_distance`
    /// upstream of this system's input plane.
    ///
    /// Pre-pends the object-space gap, then finds the downstream distance
    /// `v` that drives the total `B` element to zero — the conjugate-image
    /// condition `B_total = 0`, at which output height depends only on
    /// input height. Returns the image distance `v` (measured from the
    /// system's output plane) and the transverse magnification, which at
    /// the image plane equals the total-system `A` element.
    ///
    /// For a single thin lens this reproduces the Gaussian result of
    /// [`crate::thin_lens::image`] exactly: `v = di` and `m = -di/do`.
    ///
    /// # Errors
    ///
    /// - [`OpticsError::InvalidParameter`] if `object_distance` is not
    ///   finite.
    /// - [`OpticsError::SingularGeometry`] if no finite image plane exists
    ///   — the emerging rays are collimated (the object sits at the front
    ///   focal plane), so the image is at infinity.
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::abcd::RayTransfer;
    /// // Object 20 units before a +10 lens images at 20, m = -1.
    /// let lens = RayTransfer::thin_lens(10.0).unwrap();
    /// let img = lens.image_distance(20.0).unwrap();
    /// assert!((img.distance - 20.0).abs() < 1e-9);
    /// assert!((img.magnification + 1.0).abs() < 1e-9);
    /// ```
    pub fn image_distance(&self, object_distance: f64) -> Result<ImagePlane, OpticsError> {
        if !object_distance.is_finite() {
            return Err(OpticsError::InvalidParameter {
                name: "object_distance",
                value: object_distance,
                reason: "object distance must be finite",
            });
        }
        // System as seen from the object plane: gap of `object_distance`
        // followed by this system.
        let from_object = RayTransfer::propagation(object_distance)?.then(self);

        // The image forms at the downstream distance v where the combined
        // B element vanishes:  B(v) = B0 + v * D0 = 0  =>  v = -B0 / D0,
        // using one more propagation matrix [[1, v], [0, 1]] on the left.
        // (B_total of  prop(v) * from_object  is  from_object.b + v * from_object.d.)
        if from_object.d == 0.0 {
            return Err(OpticsError::SingularGeometry {
                reason: "image at infinity: emerging rays are collimated (object at front focus)",
            });
        }
        let v = -from_object.b / from_object.d;
        // Magnification is the A element of the full object->image matrix,
        // i.e. A_total = from_object.a + v * from_object.c.
        let magnification = from_object.a + v * from_object.c;
        Ok(ImagePlane {
            distance: v,
            magnification,
        })
    }
}

/// The conjugate image found by [`RayTransfer::image_distance`].
///
/// `distance` is measured downstream from the system's output plane
/// (positive = real image beyond the system); `magnification` is the
/// transverse magnification at that plane.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ImagePlane {
    /// Image distance from the system's output plane (positive downstream).
    pub distance: f64,
    /// Transverse magnification at the image plane.
    pub magnification: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thin_lens::image as gaussian_image;

    const EPS: f64 = 1e-9;

    #[test]
    fn identity_is_a_noop() {
        let r = Ray::new(3.0, 0.05).unwrap();
        assert_eq!(RayTransfer::identity().apply(r), r);
        assert!((RayTransfer::identity().determinant() - 1.0).abs() < EPS);
    }

    #[test]
    fn propagation_grows_height_by_slope_times_distance() {
        let out = RayTransfer::propagation(8.0)
            .unwrap()
            .apply(Ray::new(1.0, 0.25).unwrap());
        assert!((out.height - (1.0 + 8.0 * 0.25)).abs() < EPS);
        assert!((out.angle - 0.25).abs() < EPS); // slope preserved
    }

    #[test]
    fn parallel_ray_crosses_axis_at_focal_point() {
        // Ground truth: a ray parallel to the axis, refracted by a +f thin
        // lens, must cross the axis exactly one focal length downstream.
        let f = 50.0;
        let lens = RayTransfer::thin_lens(f).unwrap();
        // Parallel input ray at height h.
        let h = 4.0;
        // Lens then propagate by f.
        let sys = lens.then(&RayTransfer::propagation(f).unwrap());
        let out = sys.apply(Ray::new(h, 0.0).unwrap());
        assert!(out.height.abs() < EPS, "height at focus = {}", out.height);
        // And the slope after the lens is exactly -h/f.
        let at_lens = lens.apply(Ray::new(h, 0.0).unwrap());
        assert!(
            (at_lens.angle + h / f).abs() < EPS,
            "angle = {}",
            at_lens.angle
        );
    }

    #[test]
    fn thin_lens_matrix_effective_focal_length() {
        // The bare thin-lens matrix reports f via -1/C.
        for &f in &[10.0_f64, 25.0, -40.0, 100.0] {
            let efl = RayTransfer::thin_lens(f)
                .unwrap()
                .effective_focal_length()
                .unwrap();
            assert!((efl - f).abs() < EPS, "f = {f}, efl = {efl}");
        }
    }

    #[test]
    fn thin_lens_image_matches_gaussian_equation() {
        // Ground truth: the ABCD image distance & magnification must equal
        // the closed-form thin-lens equation 1/f = 1/do + 1/di, m = -di/do.
        let f = 10.0;
        for &do_ in &[15.0_f64, 20.0, 30.0, 5.0, 50.0] {
            let abcd = RayTransfer::thin_lens(f)
                .unwrap()
                .image_distance(do_)
                .unwrap();
            let gauss = gaussian_image(f, do_).unwrap();
            assert!(
                (abcd.distance - gauss.distance).abs() < 1e-7,
                "do = {do_}: abcd di = {}, gauss di = {}",
                abcd.distance,
                gauss.distance
            );
            assert!(
                (abcd.magnification - gauss.magnification).abs() < 1e-7,
                "do = {do_}: abcd m = {}, gauss m = {}",
                abcd.magnification,
                gauss.magnification
            );
        }
    }

    #[test]
    fn thin_lens_object_at_two_f_images_at_two_f_unit_inverted() {
        // The canonical conjugate pair, via the matrix path.
        let f = 10.0;
        let img = RayTransfer::thin_lens(f)
            .unwrap()
            .image_distance(2.0 * f)
            .unwrap();
        assert!(
            (img.distance - 2.0 * f).abs() < EPS,
            "di = {}",
            img.distance
        );
        assert!(
            (img.magnification + 1.0).abs() < EPS,
            "m = {}",
            img.magnification
        );
    }

    #[test]
    fn magnification_is_minus_v_over_u_for_a_thin_lens() {
        // Ground truth: transverse magnification m = -v/u (= -di/do).
        let f = 12.0;
        for &u in &[18.0_f64, 24.0, 40.0] {
            let img = RayTransfer::thin_lens(f)
                .unwrap()
                .image_distance(u)
                .unwrap();
            assert!(
                (img.magnification - (-img.distance / u)).abs() < 1e-9,
                "u = {u}: m = {}, -v/u = {}",
                img.magnification,
                -img.distance / u
            );
        }
    }

    #[test]
    fn refraction_thin_lens_consistency_via_two_surfaces() {
        // Ground truth: a thin lens is two refractions back-to-back at zero
        // separation. A symmetric biconvex glass lens (n=1.5, R1=+0.1,
        // R2=-0.1) should give effective focal length 0.10 (lensmaker:
        // 1/f = 0.5*(10+10) = 10), matching RayTransfer::thin_lens(0.10).
        let n = 1.5;
        let surf1 = RayTransfer::refraction(1.0, n, 0.10).unwrap();
        let surf2 = RayTransfer::refraction(n, 1.0, -0.10).unwrap();
        let lens = surf1.then(&surf2); // ray hits surf1 first
        let efl = lens.effective_focal_length().unwrap();
        assert!((efl - 0.10).abs() < 1e-12, "efl = {efl}");
        // Determinant back in air on both ends is 1.
        assert!((lens.determinant() - 1.0).abs() < 1e-12);
        // And it images like the equivalent thin lens.
        let a = lens.image_distance(0.30).unwrap();
        let b = RayTransfer::thin_lens(0.10)
            .unwrap()
            .image_distance(0.30)
            .unwrap();
        assert!((a.distance - b.distance).abs() < 1e-12);
        assert!((a.magnification - b.magnification).abs() < 1e-12);
    }

    #[test]
    fn single_surface_back_focal_length_matches_power() {
        // Single spherical surface, air (1.0) into glass (1.5), R = 0.05.
        // Power P = (n2 - n1)/R = 0.5/0.05 = 10; back focal length = n2/P
        // = 1.5/10 = 0.15.  effective_focal_length reports -1/C.
        let surf = RayTransfer::refraction(1.0, 1.5, 0.05).unwrap();
        let bfl = surf.effective_focal_length().unwrap();
        assert!((bfl - 0.15).abs() < 1e-12, "bfl = {bfl}");
    }

    #[test]
    fn flat_surface_has_zero_power_and_rescales_slope() {
        let flat = RayTransfer::refraction(1.0, 1.5, f64::INFINITY).unwrap();
        assert!(flat.c.abs() < 1e-15, "C = {}", flat.c);
        assert!((flat.d - 1.0 / 1.5).abs() < EPS);
        // No power -> no finite focal length.
        assert!(flat.effective_focal_length().is_err());
    }

    #[test]
    fn determinant_is_n_in_over_n_out() {
        // A single surface from n1 to n2 has det = n1/n2.
        let surf = RayTransfer::refraction(1.0, 1.5, 0.05).unwrap();
        assert!((surf.determinant() - 1.0 / 1.5).abs() < EPS);
        // Back-to-back into and out of glass returns to air: det = 1.
        let round = surf.then(&RayTransfer::refraction(1.5, 1.0, -0.05).unwrap());
        assert!((round.determinant() - 1.0).abs() < EPS);
    }

    #[test]
    fn composition_order_matches_sequential_apply() {
        // Pushing a ray through a composed matrix equals pushing it
        // through each element in travel order.
        let g1 = RayTransfer::propagation(7.0).unwrap();
        let l = RayTransfer::thin_lens(20.0).unwrap();
        let g2 = RayTransfer::propagation(13.0).unwrap();
        let sys = g1.then(&l).then(&g2);

        let r0 = Ray::new(2.0, -0.03).unwrap();
        let stepwise = g2.apply(l.apply(g1.apply(r0)));
        let combined = sys.apply(r0);
        assert!((stepwise.height - combined.height).abs() < EPS);
        assert!((stepwise.angle - combined.angle).abs() < EPS);
    }

    #[test]
    fn afocal_system_has_no_focal_length_or_image() {
        // Two +f lenses separated by 2f form an afocal (4f) relay: C = 0.
        let f = 25.0;
        let relay = RayTransfer::thin_lens(f)
            .unwrap()
            .then(&RayTransfer::propagation(2.0 * f).unwrap())
            .then(&RayTransfer::thin_lens(f).unwrap());
        assert!(relay.c.abs() < 1e-12, "C = {}", relay.c);
        assert!(relay.effective_focal_length().is_err());
    }

    #[test]
    fn object_at_front_focus_images_at_infinity() {
        // Object exactly one focal length before a thin lens -> collimated
        // output -> no finite image plane.
        let f = 30.0;
        let err = RayTransfer::thin_lens(f)
            .unwrap()
            .image_distance(f)
            .unwrap_err();
        assert_eq!(err.code(), "optics.singular_geometry");
    }

    #[test]
    fn invalid_inputs_fail_loud() {
        // Non-finite ray components.
        assert!(Ray::new(f64::NAN, 0.0).is_err());
        assert!(Ray::new(0.0, f64::INFINITY).is_err());
        // Zero-focal-length lens.
        assert!(RayTransfer::thin_lens(0.0).is_err());
        assert!(RayTransfer::thin_lens(f64::NAN).is_err());
        // Bad propagation distance.
        assert!(RayTransfer::propagation(f64::NAN).is_err());
        // Bad refraction inputs.
        assert!(RayTransfer::refraction(0.0, 1.5, 0.1).is_err()); // n1 <= 0
        assert!(RayTransfer::refraction(1.0, -1.0, 0.1).is_err()); // n2 < 0
        assert!(RayTransfer::refraction(1.0, 1.5, 0.0).is_err()); // r == 0
        assert!(RayTransfer::refraction(1.0, 1.5, f64::NAN).is_err()); // r NaN
                                                                       // Non-finite object distance.
        assert!(RayTransfer::thin_lens(10.0)
            .unwrap()
            .image_distance(f64::NAN)
            .is_err());
    }

    #[test]
    fn diverging_lens_images_virtual_via_matrix() {
        // f < 0: a real object gives a virtual image (v < 0) that is
        // upright and reduced (0 < m < 1) — matches the Gaussian path.
        let img = RayTransfer::thin_lens(-10.0)
            .unwrap()
            .image_distance(10.0)
            .unwrap();
        assert!(img.distance < 0.0, "di = {}", img.distance);
        assert!(
            img.magnification > 0.0 && img.magnification < 1.0,
            "m = {}",
            img.magnification
        );
        let gauss = gaussian_image(-10.0, 10.0).unwrap();
        assert!((img.distance - gauss.distance).abs() < 1e-9);
        assert!((img.magnification - gauss.magnification).abs() < 1e-9);
    }
}
