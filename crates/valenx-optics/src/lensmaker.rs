//! The lensmaker's equation: focal length from surface geometry.
//!
//! ## What
//!
//! Computes the focal length of a thin lens from the refractive index of
//! its material and the radii of curvature of its two surfaces. The
//! result feeds directly into [`crate::thin_lens`].
//!
//! ## Model
//!
//! Thin-lens lensmaker's equation, lens in air:
//!
//! ```text
//! 1/f = (n - 1) * (1/R1 - 1/R2)
//! ```
//!
//! ## Sign convention
//!
//! Radii follow the standard convention: a surface is positive if its
//! centre of curvature lies on the outgoing (far) side of the lens, and
//! negative if it lies on the incoming side. For a biconvex lens this
//! gives `R1 > 0` and `R2 < 0`. A flat surface has an infinite radius;
//! pass [`f64::INFINITY`] (its `1/R` term is zero).
//!
//! ## Honest scope
//!
//! The thin-lens form only: the lens is assumed infinitesimally thin, so
//! the `(n - 1)^2 * d / (n * R1 * R2)` thickness correction of the full
//! lensmaker's equation is dropped. The surrounding medium is air
//! (`n_medium = 1`). A single refractive index is used (no dispersion),
//! and surfaces are assumed spherical. Research / educational grade — not
//! an optical-design tool.

use crate::error::{validate_index, OpticsError};
use crate::thin_lens::ThinLens;

/// A thin lens described by its material index and two surface radii.
///
/// Use [`focal_length`] (or [`Lens::focal_length`]) to obtain the focal
/// length, or [`Lens::to_thin_lens`] to convert directly into a
/// [`ThinLens`] ready for imaging.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Lens {
    /// Refractive index of the lens material (`> 0`, and `!= 1` so the
    /// lens actually refracts).
    pub n: f64,
    /// Radius of curvature of the first (front) surface. Use
    /// [`f64::INFINITY`] for a flat surface.
    pub r1: f64,
    /// Radius of curvature of the second (back) surface. Use
    /// [`f64::INFINITY`] for a flat surface.
    pub r2: f64,
}

impl Lens {
    /// Construct a lens from its index and the two surface radii.
    ///
    /// # Errors
    ///
    /// Returns [`OpticsError::InvalidParameter`] if `n` is not finite and
    /// positive, if `n == 1` (no refraction — the focal length would be
    /// infinite), if either radius is `NaN`, or if either radius is
    /// exactly zero (a zero radius is a degenerate point surface).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::lensmaker::Lens;
    /// // Symmetric biconvex crown-glass lens.
    /// let lens = Lens::new(1.5, 0.10, -0.10).unwrap();
    /// assert!(lens.focal_length().unwrap() > 0.0);
    /// ```
    pub fn new(n: f64, r1: f64, r2: f64) -> Result<Self, OpticsError> {
        let n = validate_index("n", n)?;
        if (n - 1.0).abs() < f64::EPSILON {
            return Err(OpticsError::InvalidParameter {
                name: "n",
                value: n,
                reason: "index of 1 (matching the surrounding air) gives no refraction",
            });
        }
        validate_radius("r1", r1)?;
        validate_radius("r2", r2)?;
        Ok(Self { n, r1, r2 })
    }

    /// The focal length of this lens via the lensmaker's equation.
    ///
    /// See [`focal_length`] for the full contract.
    ///
    /// # Errors
    ///
    /// Propagates the errors of [`focal_length`].
    pub fn focal_length(&self) -> Result<f64, OpticsError> {
        focal_length(self.n, self.r1, self.r2)
    }

    /// Optical power `P = 1 / f` of this lens (dioptres if radii are in
    /// metres).
    ///
    /// # Errors
    ///
    /// Propagates the errors of [`focal_length`]; the power is its
    /// reciprocal.
    pub fn power(&self) -> Result<f64, OpticsError> {
        Ok(1.0 / self.focal_length()?)
    }

    /// Convert this lens directly into a [`ThinLens`] for imaging.
    ///
    /// # Errors
    ///
    /// Propagates the errors of [`focal_length`], plus any error from
    /// [`ThinLens::new`] (e.g. a plano-on-both-sides flat plate, whose
    /// focal length is infinite).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::lensmaker::Lens;
    /// let lens = Lens::new(1.5, 0.10, -0.10).unwrap();
    /// let tl = lens.to_thin_lens().unwrap();
    /// assert!(tl.is_converging());
    /// ```
    pub fn to_thin_lens(&self) -> Result<ThinLens, OpticsError> {
        ThinLens::new(self.focal_length()?)
    }
}

/// Validate a radius of curvature: finite-or-infinite, but never `NaN`
/// and never exactly zero.
fn validate_radius(name: &'static str, r: f64) -> Result<f64, OpticsError> {
    if r.is_nan() {
        return Err(OpticsError::InvalidParameter {
            name,
            value: r,
            reason: "radius of curvature must not be NaN",
        });
    }
    if r == 0.0 {
        return Err(OpticsError::InvalidParameter {
            name,
            value: r,
            reason: "radius of curvature of zero is a degenerate (point) surface",
        });
    }
    Ok(r)
}

/// `1/R`, defined as `0` for an infinite (flat) radius.
fn curvature(r: f64) -> f64 {
    if r.is_infinite() {
        0.0
    } else {
        1.0 / r
    }
}

/// Focal length from the lensmaker's equation.
///
/// Computes `f` from `1/f = (n - 1) * (1/R1 - 1/R2)` for a thin lens in
/// air. Flat surfaces are signalled by passing [`f64::INFINITY`] for the
/// corresponding radius (their `1/R` term vanishes).
///
/// # Errors
///
/// - [`OpticsError::InvalidParameter`] if `n` is non-finite, non-positive,
///   or equal to 1, or if either radius is `NaN` or exactly zero.
/// - [`OpticsError::SingularGeometry`] if the two surfaces have equal
///   curvature (`1/R1 == 1/R2`, e.g. a flat plate or a concentric shell),
///   so the optical power is zero and the focal length is infinite.
///
/// # Examples
///
/// ```
/// use valenx_optics::lensmaker::focal_length;
/// // Symmetric biconvex, n = 1.5, R1 = +0.10 m, R2 = -0.10 m:
/// // 1/f = 0.5 * (10 - (-10)) = 10 dioptres  =>  f = 0.10 m.
/// let f = focal_length(1.5, 0.10, -0.10).unwrap();
/// assert!((f - 0.10).abs() < 1e-12);
/// ```
pub fn focal_length(n: f64, r1: f64, r2: f64) -> Result<f64, OpticsError> {
    let n = validate_index("n", n)?;
    if (n - 1.0).abs() < f64::EPSILON {
        return Err(OpticsError::InvalidParameter {
            name: "n",
            value: n,
            reason: "index of 1 (matching the surrounding air) gives no refraction",
        });
    }
    let r1 = validate_radius("r1", r1)?;
    let r2 = validate_radius("r2", r2)?;

    let inv_f = (n - 1.0) * (curvature(r1) - curvature(r2));
    if inv_f == 0.0 {
        return Err(OpticsError::SingularGeometry {
            reason: "equal surface curvatures give zero power: focal length is infinite",
        });
    }
    Ok(1.0 / inv_f)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thin_lens::ImageKind;

    const EPS: f64 = 1e-12;

    #[test]
    fn symmetric_biconvex_is_converging() {
        // n = 1.5, R1 = +0.1, R2 = -0.1.
        // 1/f = 0.5 * (1/0.1 - 1/(-0.1)) = 0.5 * (10 + 10) = 10.
        let f = focal_length(1.5, 0.10, -0.10).unwrap();
        assert!((f - 0.10).abs() < EPS, "f = {f}");
        assert!(f > 0.0);
    }

    #[test]
    fn symmetric_biconcave_is_diverging() {
        // Flip both radii: R1 = -0.1, R2 = +0.1.
        // 1/f = 0.5 * (-10 - 10) = -10 => f = -0.1.
        let f = focal_length(1.5, -0.10, 0.10).unwrap();
        assert!((f - (-0.10)).abs() < EPS, "f = {f}");
        assert!(f < 0.0);
    }

    #[test]
    fn plano_convex_uses_only_the_curved_surface() {
        // Flat front (R1 = inf), convex back (R2 = -0.05), n = 1.5.
        // 1/f = 0.5 * (0 - 1/(-0.05)) = 0.5 * 20 = 10 => f = 0.1.
        let f = focal_length(1.5, f64::INFINITY, -0.05).unwrap();
        assert!((f - 0.10).abs() < EPS, "f = {f}");
    }

    #[test]
    fn plano_convex_orientation_independent() {
        // The two ways to orient a plano-convex lens give the same |f|.
        let a = focal_length(1.5, 0.05, f64::INFINITY).unwrap();
        let b = focal_length(1.5, f64::INFINITY, -0.05).unwrap();
        assert!((a - b).abs() < EPS, "a = {a}, b = {b}");
    }

    #[test]
    fn higher_index_shortens_focal_length() {
        // Same geometry, larger n => more power => shorter f.
        let f_low = focal_length(1.5, 0.10, -0.10).unwrap();
        let f_high = focal_length(1.8, 0.10, -0.10).unwrap();
        assert!(f_high < f_low, "f_high = {f_high}, f_low = {f_low}");
        // Quantitative: 1/f scales as (n - 1).
        let ratio = (1.0 / f_high) / (1.0 / f_low);
        assert!((ratio - (0.8 / 0.5)).abs() < 1e-9, "ratio = {ratio}");
    }

    #[test]
    fn flat_plate_has_infinite_focal_length() {
        // Both surfaces flat: zero power.
        let err = focal_length(1.5, f64::INFINITY, f64::INFINITY).unwrap_err();
        assert_eq!(err.code(), "optics.singular_geometry");
    }

    #[test]
    fn equal_curvatures_are_singular() {
        // R1 == R2 (a concentric meniscus shell): 1/R1 - 1/R2 = 0.
        let err = focal_length(1.5, 0.10, 0.10).unwrap_err();
        assert_eq!(err.code(), "optics.singular_geometry");
    }

    #[test]
    fn index_one_is_rejected() {
        let err = focal_length(1.0, 0.10, -0.10).unwrap_err();
        assert_eq!(err.code(), "optics.invalid_parameter");
    }

    #[test]
    fn invalid_radii_rejected() {
        assert!(focal_length(1.5, 0.0, -0.10).is_err()); // zero radius
        assert!(focal_length(1.5, f64::NAN, -0.10).is_err()); // NaN
        assert!(Lens::new(1.5, 0.10, 0.0).is_err());
    }

    #[test]
    fn lens_struct_matches_free_function() {
        let lens = Lens::new(1.6, 0.12, -0.18).unwrap();
        let a = lens.focal_length().unwrap();
        let b = focal_length(1.6, 0.12, -0.18).unwrap();
        assert!((a - b).abs() < EPS, "a = {a}, b = {b}");
        let p = lens.power().unwrap();
        assert!((p - 1.0 / b).abs() < EPS, "p = {p}");
    }

    #[test]
    fn lensmaker_feeds_thin_lens_imaging() {
        // End-to-end: build a converging lens, image an object at 2f,
        // expect a real inverted unit-magnification image.
        let lens = Lens::new(1.5, 0.10, -0.10).unwrap();
        let f = lens.focal_length().unwrap();
        assert!((f - 0.10).abs() < EPS);
        let tl = lens.to_thin_lens().unwrap();
        let img = tl.image(2.0 * f).unwrap();
        assert!(
            (img.distance - 2.0 * f).abs() < 1e-9,
            "di = {}",
            img.distance
        );
        assert!(
            (img.magnification + 1.0).abs() < 1e-9,
            "m = {}",
            img.magnification
        );
        assert_eq!(img.kind, ImageKind::Real);
    }
}
