//! # valenx-optics
//!
//! Closed-form geometric (ray) optics: refraction at planar interfaces,
//! the thin-lens equation, and the lensmaker's equation.
//!
//! ## What
//!
//! Three textbook models, each with a small free-function API plus an
//! ergonomic struct wrapper:
//!
//! - **Refraction** ([`refraction`]) — Snell's law
//!   `n1 sin(t1) = n2 sin(t2)`, the critical angle `asin(n2 / n1)`, and
//!   total internal reflection. See [`refraction::refract_angle`],
//!   [`refraction::critical_angle_deg`], [`refraction::classify_ray`] and
//!   [`refraction::Interface`].
//! - **Thin lens** ([`thin_lens`]) — the Gaussian equation
//!   `1/f = 1/do + 1/di` with magnification `m = -di / do`, image
//!   classification (real/virtual, upright/inverted, magnified/reduced),
//!   and the inverse `do = f(m - 1)/m` for the object distance that gives
//!   a target magnification. See [`thin_lens::image`],
//!   [`thin_lens::object_distance_for_magnification`] and
//!   [`thin_lens::ThinLens`].
//! - **Lensmaker** ([`lensmaker`]) — focal length from material and
//!   surface geometry, `1/f = (n - 1)(1/R1 - 1/R2)`, plus its design
//!   inverse [`lensmaker::radius_for_focal_length`] (the second radius
//!   that hits a target focal length). See [`lensmaker::focal_length`]
//!   and [`lensmaker::Lens`], which can hand a [`thin_lens::ThinLens`]
//!   straight to the imaging code.
//! - **Ray-transfer (ABCD) matrices** ([`abcd`]) — the paraxial
//!   ray-transfer-matrix formalism: a [`abcd::Ray`] (height, slope) pushed
//!   through composable [`abcd::RayTransfer`] elements (free-space gaps,
//!   thin lenses, and flat / spherical refracting surfaces via Snell in
//!   the paraxial limit). Cascade with [`abcd::RayTransfer::then`], read
//!   the [effective focal length](abcd::RayTransfer::effective_focal_length),
//!   and locate conjugate images
//!   ([`abcd::RayTransfer::image_distance`]) of multi-element systems —
//!   the matrix companion to the closed-form thin-lens / lensmaker models.
//!
//! ## Model
//!
//! Everything is first-order (paraxial) ray geometry expressed as
//! elementary algebra and trigonometry. Angles are in degrees at the API
//! boundary (converted to radians internally); the focal length, object,
//! image and radius distances share whatever length unit the caller
//! supplies, and the magnification and refractive index are
//! dimensionless. A consistent sign convention is documented on each
//! module: positive focal length / radii for converging geometry,
//! positive image distance for a real image, positive magnification for
//! an upright image.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are the standard introductory
//! closed-form / first-order models — **not** a clinical, medical, or
//! production optical-engineering tool. The crate models neither
//! diffraction nor interference nor polarisation; it does not compute
//! Fresnel reflection / transmission coefficients; it assumes thin,
//! aberration-free, spherical, single-wavelength (non-dispersive) optics;
//! and it does not ray-trace through real thick-lens systems. For lens
//! design, tolerancing, or any safety-relevant optical decision, use a
//! validated optical-design package, not this crate.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod abcd;
pub mod error;
pub mod lensmaker;
pub mod refraction;
pub mod thin_lens;

pub use abcd::{ImagePlane, Ray, RayTransfer};
pub use error::{ErrorCategory, OpticsError};
pub use lensmaker::{focal_length, radius_for_focal_length, Lens};
pub use refraction::{classify_ray, critical_angle_deg, refract_angle, Interface, RayOutcome};
pub use thin_lens::{
    image, object_distance_for_magnification, Image, ImageKind, Orientation, ThinLens,
};

#[cfg(test)]
mod integration_tests {
    //! Cross-module ground-truth checks that exercise the public
    //! re-exports together.

    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn end_to_end_camera_lens() {
        // Build a +50 mm-ish converging lens from glass, then image a
        // distant object and confirm it lands near the focal plane,
        // small and inverted (the camera regime).
        let lens = Lens::new(1.5, 0.05, -0.05).unwrap();
        let f = lens.focal_length().unwrap();
        // 1/f = 0.5*(20 + 20) = 20 => f = 0.05 m.
        assert!((f - 0.05).abs() < 1e-12, "f = {f}");

        let img = image(f, 5.0).unwrap(); // object 5 m away
        assert!(img.distance > f, "di = {}", img.distance);
        assert!(img.distance < 1.02 * f, "di = {}", img.distance);
        assert_eq!(img.kind, ImageKind::Real);
        assert_eq!(img.orientation, Orientation::Inverted);
        assert!(img.is_reduced());
    }

    #[test]
    fn snell_and_critical_angle_consistency() {
        // The refraction angle approaches 90 deg exactly as the
        // incidence approaches the critical angle from below.
        let n1 = 1.5;
        let n2 = 1.0;
        let c = critical_angle_deg(n1, n2).unwrap();
        let t2 = refract_angle(n1, n2, c - 1e-7).unwrap();
        assert!((90.0 - t2).abs() < 1e-2, "t2 = {t2}");

        // And just beyond it, TIR.
        match classify_ray(n1, n2, c + 1e-3).unwrap() {
            RayOutcome::TotallyReflected { critical_deg } => {
                assert!((critical_deg - c).abs() < EPS);
            }
            other => panic!("expected TIR, got {other:?}"),
        }
    }

    #[test]
    fn diverging_lens_from_lensmaker_is_virtual_upright_reduced() {
        // A biconcave lens (built via the lensmaker) must image any real
        // object as virtual, upright, reduced.
        let lens = Lens::new(1.5, -0.10, 0.10).unwrap();
        let tl = lens.to_thin_lens().unwrap();
        assert!(tl.is_diverging());
        let img = tl.image(0.20).unwrap();
        assert_eq!(img.kind, ImageKind::Virtual);
        assert_eq!(img.orientation, Orientation::Upright);
        assert!(img.is_reduced());
    }

    #[test]
    fn interface_round_trip_via_reexport() {
        let iface = Interface::new(1.0, 1.5).unwrap();
        let into = iface.refract(40.0).unwrap();
        let back = iface.reversed().refract(into).unwrap();
        assert!((back - 40.0).abs() < 1e-7, "back = {back}");
    }

    #[test]
    fn lensmaker_abcd_and_thin_lens_all_agree() {
        // Build a converging lens with the lensmaker, then check the ABCD
        // matrix path and the Gaussian thin-lens path produce the same
        // image for the same object — three independent code paths, one
        // ground truth.
        let lens = Lens::new(1.5, 0.10, -0.10).unwrap();
        let f = lens.focal_length().unwrap(); // 0.10 m
        let object = 0.30;

        // ABCD: bare thin-lens matrix of focal length f.
        let abcd = RayTransfer::thin_lens(f)
            .unwrap()
            .image_distance(object)
            .unwrap();
        // Gaussian closed form.
        let gauss = image(f, object).unwrap();

        assert!((abcd.distance - gauss.distance).abs() < 1e-12);
        assert!((abcd.magnification - gauss.magnification).abs() < 1e-12);
        // And m = -v/u holds on the matrix result.
        assert!((abcd.magnification - (-abcd.distance / object)).abs() < 1e-12);
    }

    #[test]
    fn abcd_two_surface_lens_matches_lensmaker_focal_length() {
        // Compose two spherical refractions (air|glass and glass|air) at
        // zero thickness and confirm the effective focal length equals the
        // lensmaker's value for the same surfaces.
        let lens = Lens::new(1.5, 0.10, -0.10).unwrap();
        let f_lm = lens.focal_length().unwrap();
        let sys = RayTransfer::refraction(1.0, 1.5, 0.10)
            .unwrap()
            .then(&RayTransfer::refraction(1.5, 1.0, -0.10).unwrap());
        let f_abcd = sys.effective_focal_length().unwrap();
        assert!(
            (f_abcd - f_lm).abs() < 1e-12,
            "abcd {f_abcd} vs lensmaker {f_lm}"
        );
    }
}
