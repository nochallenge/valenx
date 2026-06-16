//! The thin-lens equation and transverse magnification.
//!
//! ## What
//!
//! A thin lens of focal length `f` forms an image of an object placed a
//! distance `do` in front of it. This module solves for the image
//! distance `di`, the transverse magnification `m`, and classifies the
//! image (real vs. virtual, upright vs. inverted, magnified vs. reduced).
//! It also inverts the magnification relation: given a target `m`, it
//! returns the object distance `do = f (m - 1) / m` that produces it.
//!
//! ## Model
//!
//! Thin-lens (Gaussian) equation and the magnification:
//!
//! ```text
//! 1/f = 1/do + 1/di          =>   di = (f * do) / (do - f)
//! m   = -di / do
//! ```
//!
//! ## Sign convention
//!
//! This module uses the standard Cartesian-style convention common in
//! introductory texts:
//!
//! - `f > 0` for a converging (convex) lens, `f < 0` for a diverging
//!   (concave) lens.
//! - `do > 0` for a real object in front of the lens (the only case the
//!   solver accepts).
//! - `di > 0` is a real image on the far side; `di < 0` is a virtual
//!   image on the same side as the object.
//! - `m > 0` is upright, `m < 0` is inverted; `abs(m) > 1` is magnified,
//!   `abs(m) < 1` is reduced.
//!
//! ## Honest scope
//!
//! Paraxial (first-order, "thin-lens") approximation only: rays are
//! assumed to make small angles with the axis, the lens has negligible
//! thickness, and there is a single focal length. No spherical /
//! chromatic aberration, no thick-lens principal planes, no off-axis
//! field curvature or distortion, no finite-aperture / depth-of-field
//! effects. Research / educational grade — not an optical-design tool.

use crate::error::OpticsError;

/// Whether an image is real (rays physically converge) or virtual (rays
/// only appear to diverge from it).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum ImageKind {
    /// Rays actually cross; the image can be projected onto a screen.
    Real,
    /// Rays only appear to come from the image point; it cannot be
    /// projected.
    Virtual,
}

/// Whether an image has the same orientation as the object (upright) or
/// is flipped (inverted).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum Orientation {
    /// Same way up as the object.
    Upright,
    /// Flipped relative to the object.
    Inverted,
}

/// The full paraxial image of an object through a thin lens.
///
/// Produced by [`image`] (and the convenience constructors on
/// [`ThinLens`]). Distances are in the same length unit as the inputs;
/// the magnification is dimensionless.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Image {
    /// Image distance `di`. Positive: real image on the far side.
    /// Negative: virtual image on the object's side.
    pub distance: f64,
    /// Transverse (lateral) magnification `m = -di / do`.
    pub magnification: f64,
    /// Real or virtual.
    pub kind: ImageKind,
    /// Upright or inverted.
    pub orientation: Orientation,
}

impl Image {
    /// Whether this image is larger than the object (`abs(m) > 1`).
    ///
    /// Uses a tiny tolerance so the unit-magnification case (`abs(m)`
    /// exactly 1, e.g. an object at `2f`) is reported as neither
    /// magnified nor reduced.
    pub fn is_magnified(&self) -> bool {
        self.magnification.abs() > 1.0 + 1e-12
    }

    /// Whether this image is smaller than the object (`abs(m) < 1`).
    pub fn is_reduced(&self) -> bool {
        self.magnification.abs() < 1.0 - 1e-12
    }
}

/// A thin lens, identified solely by its focal length.
///
/// `f > 0` converging, `f < 0` diverging. A focal length of exactly zero
/// is rejected at construction because the thin-lens equation degenerates
/// there.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ThinLens {
    /// Focal length. Positive for converging, negative for diverging.
    pub focal_length: f64,
}

impl ThinLens {
    /// Construct a thin lens from its focal length.
    ///
    /// # Errors
    ///
    /// Returns [`OpticsError::InvalidParameter`] if `focal_length` is not
    /// finite, or [`OpticsError::SingularGeometry`] if it is exactly
    /// zero (an undefined / infinite-power lens).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::thin_lens::ThinLens;
    /// let converging = ThinLens::new(10.0).unwrap();
    /// assert!(converging.is_converging());
    /// assert!(ThinLens::new(0.0).is_err());
    /// ```
    pub fn new(focal_length: f64) -> Result<Self, OpticsError> {
        if !focal_length.is_finite() {
            return Err(OpticsError::InvalidParameter {
                name: "focal_length",
                value: focal_length,
                reason: "focal length must be finite",
            });
        }
        if focal_length == 0.0 {
            return Err(OpticsError::SingularGeometry {
                reason: "focal length of zero has undefined (infinite) optical power",
            });
        }
        Ok(Self { focal_length })
    }

    /// Whether this is a converging (positive-focal-length) lens.
    pub fn is_converging(&self) -> bool {
        self.focal_length > 0.0
    }

    /// Whether this is a diverging (negative-focal-length) lens.
    pub fn is_diverging(&self) -> bool {
        self.focal_length < 0.0
    }

    /// Optical power `P = 1 / f` (dioptres if `f` is in metres).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_optics::thin_lens::ThinLens;
    /// // A 0.5 m focal length is a +2 dioptre lens.
    /// let p = ThinLens::new(0.5).unwrap().power();
    /// assert!((p - 2.0).abs() < 1e-12);
    /// ```
    pub fn power(&self) -> f64 {
        1.0 / self.focal_length
    }

    /// Form an image of an object a distance `object_distance` in front
    /// of this lens.
    ///
    /// See [`image`] for the full contract.
    ///
    /// # Errors
    ///
    /// Propagates the errors of [`image`].
    pub fn image(&self, object_distance: f64) -> Result<Image, OpticsError> {
        image(self.focal_length, object_distance)
    }

    /// The object distance that yields a target transverse magnification
    /// `m` through this lens. See
    /// [`object_distance_for_magnification`].
    ///
    /// # Errors
    ///
    /// Propagates the errors of
    /// [`object_distance_for_magnification`].
    pub fn object_distance_for_magnification(&self, m: f64) -> Result<f64, OpticsError> {
        object_distance_for_magnification(self.focal_length, m)
    }
}

/// Solve the thin-lens equation for the image of a real object.
///
/// Given a focal length `f` and a (positive) object distance `do`,
/// returns the [`Image`]: its distance `di = (f * do) / (do - f)`, its
/// magnification `m = -di / do`, and its real/virtual and
/// upright/inverted classification.
///
/// # Errors
///
/// - [`OpticsError::InvalidParameter`] if `f` is non-finite, if `f` is
///   exactly zero, or if `object_distance` is non-finite or not strictly
///   positive (this solver models real objects only).
/// - [`OpticsError::SingularGeometry`] if the object sits exactly at the
///   front focal point of a converging lens (`do == f`), so the rays
///   emerge collimated and the image is at infinity.
///
/// # Examples
///
/// ```
/// use valenx_optics::thin_lens::{image, ImageKind, Orientation};
/// // Object at twice the focal length images at 2f, inverted, same size.
/// let img = image(10.0, 20.0).unwrap();
/// assert!((img.distance - 20.0).abs() < 1e-9);
/// assert!((img.magnification + 1.0).abs() < 1e-9); // m = -1
/// assert_eq!(img.kind, ImageKind::Real);
/// assert_eq!(img.orientation, Orientation::Inverted);
/// ```
pub fn image(f: f64, object_distance: f64) -> Result<Image, OpticsError> {
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
    if !object_distance.is_finite() {
        return Err(OpticsError::InvalidParameter {
            name: "object_distance",
            value: object_distance,
            reason: "object distance must be finite",
        });
    }
    if object_distance <= 0.0 {
        return Err(OpticsError::InvalidParameter {
            name: "object_distance",
            value: object_distance,
            reason: "object distance must be strictly positive (real object)",
        });
    }

    let denom = object_distance - f;
    if denom == 0.0 {
        return Err(OpticsError::SingularGeometry {
            reason: "object at the front focal point (do == f): image is at infinity",
        });
    }

    let di = (f * object_distance) / denom;
    let m = -di / object_distance;

    let kind = if di > 0.0 {
        ImageKind::Real
    } else {
        ImageKind::Virtual
    };
    let orientation = if m >= 0.0 {
        Orientation::Upright
    } else {
        Orientation::Inverted
    };

    Ok(Image {
        distance: di,
        magnification: m,
        kind,
        orientation,
    })
}

/// The object distance that yields a target transverse magnification `m`
/// through a thin lens of focal length `f`, inverting the magnification
/// relation of [`image`]:
///
/// ```text
/// do = f * (m - 1) / m
/// ```
///
/// derived from `m = -di/do` and `1/f = 1/do + 1/di`. This answers the
/// standard "where do I place the object for this magnification?"
/// question — `m = -1` for a 1:1 real (inverted) image, `m = -0.5` for a
/// reduced real image, `m = 2` for a 2x upright magnifying-glass view.
///
/// # Errors
///
/// - [`OpticsError::InvalidParameter`] if `f` is non-finite, `f` is
///   exactly zero, or `m` is non-finite or zero (`m -> 0` needs an object
///   at infinity).
/// - [`OpticsError::SingularGeometry`] if the required object distance is
///   not strictly positive — the target magnification is unreachable with
///   a real object for this lens (e.g. `m` in `(0, 1)` for a converging
///   lens, or `m = 1`, which needs `do = 0`).
///
/// # Examples
///
/// ```
/// use valenx_optics::thin_lens::{image, object_distance_for_magnification};
/// // For a 1:1 macro (m = -1) on a +10 lens, place the object at 2f.
/// let d = object_distance_for_magnification(10.0, -1.0).unwrap();
/// assert!((d - 20.0).abs() < 1e-9);
/// // ...and imaging from there indeed gives m = -1.
/// assert!((image(10.0, d).unwrap().magnification + 1.0).abs() < 1e-9);
/// ```
pub fn object_distance_for_magnification(f: f64, m: f64) -> Result<f64, OpticsError> {
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
    if !m.is_finite() || m == 0.0 {
        return Err(OpticsError::InvalidParameter {
            name: "m",
            value: m,
            reason:
                "magnification must be finite and non-zero (m -> 0 needs an object at infinity)",
        });
    }

    let object_distance = f * (m - 1.0) / m;
    if object_distance <= 0.0 {
        return Err(OpticsError::SingularGeometry {
            reason:
                "target magnification is unreachable with a real object (do <= 0) for this lens",
        });
    }
    Ok(object_distance)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn object_at_two_f_images_at_two_f_unit_inverted() {
        // The canonical "object at 2f -> image at 2f, m = -1" case.
        let f = 10.0;
        let img = image(f, 2.0 * f).unwrap();
        assert!(
            (img.distance - 2.0 * f).abs() < EPS,
            "di = {}",
            img.distance
        );
        assert!(
            (img.magnification - (-1.0)).abs() < EPS,
            "m = {}",
            img.magnification
        );
        assert_eq!(img.kind, ImageKind::Real);
        assert_eq!(img.orientation, Orientation::Inverted);
        // Exactly unit magnitude -> neither magnified nor reduced.
        assert!(!img.is_magnified());
        assert!(!img.is_reduced());
    }

    #[test]
    fn object_at_infinity_images_at_focal_point() {
        // As do -> infinity, di -> f and m -> 0. Use a very large do.
        let f = 10.0;
        let do_far = 1.0e9;
        let img = image(f, do_far).unwrap();
        assert!((img.distance - f).abs() < 1e-6, "di = {}", img.distance);
        assert!(img.magnification.abs() < 1e-6, "m = {}", img.magnification);
        assert_eq!(img.kind, ImageKind::Real);
    }

    #[test]
    fn converging_object_inside_f_gives_virtual_upright_magnified() {
        // Magnifying-glass regime: object closer than f.
        let f = 10.0;
        let img = image(f, 5.0).unwrap();
        // di = (10*5)/(5-10) = 50 / -5 = -10 (virtual, same side).
        assert!(
            (img.distance - (-10.0)).abs() < EPS,
            "di = {}",
            img.distance
        );
        // m = -di/do = -(-10)/5 = 2 (upright, 2x).
        assert!(
            (img.magnification - 2.0).abs() < EPS,
            "m = {}",
            img.magnification
        );
        assert_eq!(img.kind, ImageKind::Virtual);
        assert_eq!(img.orientation, Orientation::Upright);
        assert!(img.is_magnified());
    }

    #[test]
    fn converging_object_beyond_two_f_gives_real_inverted_reduced() {
        // Camera regime: distant object, small real inverted image.
        let f = 10.0;
        let img = image(f, 30.0).unwrap();
        // di = (10*30)/(30-10) = 300/20 = 15.
        assert!((img.distance - 15.0).abs() < EPS, "di = {}", img.distance);
        // m = -15/30 = -0.5.
        assert!(
            (img.magnification - (-0.5)).abs() < EPS,
            "m = {}",
            img.magnification
        );
        assert_eq!(img.kind, ImageKind::Real);
        assert_eq!(img.orientation, Orientation::Inverted);
        assert!(img.is_reduced());
    }

    #[test]
    fn diverging_lens_always_virtual_upright_reduced() {
        // A diverging lens (f < 0) forms a virtual, upright, reduced
        // image for ANY real object distance. Check several.
        let f = -10.0;
        for &d in &[1.0_f64, 5.0, 10.0, 25.0, 100.0] {
            let img = image(f, d).unwrap();
            assert_eq!(img.kind, ImageKind::Virtual, "do = {d}");
            assert_eq!(img.orientation, Orientation::Upright, "do = {d}");
            assert!(img.is_reduced(), "do = {d}, m = {}", img.magnification);
            // Virtual image is on the same side -> di < 0.
            assert!(img.distance < 0.0, "do = {d}, di = {}", img.distance);
            // Upright -> m > 0; reduced -> m < 1.
            assert!(
                img.magnification > 0.0 && img.magnification < 1.0,
                "do = {d}"
            );
        }
    }

    #[test]
    fn diverging_lens_exact_value() {
        // f = -10, do = 10: di = (-10*10)/(10-(-10)) = -100/20 = -5.
        let img = image(-10.0, 10.0).unwrap();
        assert!((img.distance - (-5.0)).abs() < EPS, "di = {}", img.distance);
        // m = -(-5)/10 = 0.5.
        assert!(
            (img.magnification - 0.5).abs() < EPS,
            "m = {}",
            img.magnification
        );
    }

    #[test]
    fn thin_lens_equation_is_consistent() {
        // Directly verify 1/f = 1/do + 1/di for the computed di.
        let f = 7.5;
        for &d in &[3.0_f64, 7.0, 7.6, 15.0, 50.0] {
            let img = image(f, d).unwrap();
            let lhs = 1.0 / f;
            let rhs = 1.0 / d + 1.0 / img.distance;
            assert!(
                (lhs - rhs).abs() < 1e-9,
                "do = {d}: lhs = {lhs}, rhs = {rhs}"
            );
        }
    }

    #[test]
    fn object_at_focal_point_is_singular() {
        let err = image(10.0, 10.0).unwrap_err();
        assert_eq!(err.code(), "optics.singular_geometry");
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        assert!(image(0.0, 10.0).is_err()); // f == 0
        assert!(image(f64::NAN, 10.0).is_err());
        assert!(image(10.0, 0.0).is_err()); // do == 0
        assert!(image(10.0, -5.0).is_err()); // do < 0
        assert!(image(10.0, f64::INFINITY).is_err()); // do non-finite
    }

    #[test]
    fn thin_lens_power_and_kind() {
        let l = ThinLens::new(0.25).unwrap();
        assert!((l.power() - 4.0).abs() < EPS);
        assert!(l.is_converging());
        assert!(!l.is_diverging());

        let d = ThinLens::new(-0.5).unwrap();
        assert!((d.power() - (-2.0)).abs() < EPS);
        assert!(d.is_diverging());
    }

    #[test]
    fn thin_lens_image_matches_free_function() {
        let l = ThinLens::new(12.0).unwrap();
        let a = l.image(40.0).unwrap();
        let b = image(12.0, 40.0).unwrap();
        assert!((a.distance - b.distance).abs() < EPS);
        assert!((a.magnification - b.magnification).abs() < EPS);
    }

    #[test]
    fn zero_focal_length_lens_rejected() {
        let err = ThinLens::new(0.0).unwrap_err();
        assert_eq!(err.code(), "optics.singular_geometry");
        assert!(ThinLens::new(f64::INFINITY).is_err());
    }

    #[test]
    fn object_distance_for_magnification_round_trips() {
        // Solve for do at a target m, then image() must reproduce m.
        let f = 10.0;
        for &m in &[-2.0_f64, -1.0, -0.5, 3.0] {
            let do_ = object_distance_for_magnification(f, m).unwrap();
            assert!(do_ > 0.0, "do should be a real object, got {do_}");
            let img = image(f, do_).unwrap();
            assert!(
                (img.magnification - m).abs() < 1e-9,
                "m {m} got {}",
                img.magnification
            );
        }
    }

    #[test]
    fn object_distance_for_magnification_known_values() {
        // m = -1 (1:1) -> 2f; m = -0.5 (camera) -> 3f; m = 2 (loupe) -> f/2.
        assert!((object_distance_for_magnification(10.0, -1.0).unwrap() - 20.0).abs() < EPS);
        assert!((object_distance_for_magnification(10.0, -0.5).unwrap() - 30.0).abs() < EPS);
        assert!((object_distance_for_magnification(10.0, 2.0).unwrap() - 5.0).abs() < EPS);
    }

    #[test]
    fn diverging_magnification_inverse() {
        // Diverging lens: only 0 < m < 1 is reachable. f=-10, m=0.5 -> do=10.
        let do_ = object_distance_for_magnification(-10.0, 0.5).unwrap();
        assert!((do_ - 10.0).abs() < EPS, "do = {do_}");
        assert!((image(-10.0, do_).unwrap().magnification - 0.5).abs() < EPS);
    }

    #[test]
    fn unreachable_magnification_is_rejected() {
        // Converging + m in (0,1) needs do <= 0 (no real object); m = 1
        // needs do = 0.
        assert!(object_distance_for_magnification(10.0, 0.5).is_err());
        assert!(object_distance_for_magnification(10.0, 1.0).is_err());
    }

    #[test]
    fn object_distance_for_magnification_rejects_bad_inputs() {
        assert!(object_distance_for_magnification(0.0, -1.0).is_err()); // f == 0
        assert!(object_distance_for_magnification(f64::NAN, -1.0).is_err());
        assert!(object_distance_for_magnification(10.0, 0.0).is_err()); // m == 0
        assert!(object_distance_for_magnification(10.0, f64::INFINITY).is_err());
    }

    #[test]
    fn thin_lens_magnification_inverse_matches_free_function() {
        let l = ThinLens::new(12.0).unwrap();
        let a = l.object_distance_for_magnification(-0.5).unwrap();
        let b = object_distance_for_magnification(12.0, -0.5).unwrap();
        assert!((a - b).abs() < EPS);
    }
}
