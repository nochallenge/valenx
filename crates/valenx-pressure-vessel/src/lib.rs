//! # valenx-pressure-vessel — pressure-vessel stress analysis
//!
//! Closed-form elastic stresses for internally (and externally)
//! pressurised cylindrical and spherical vessels: the thin-wall membrane
//! formulae and the exact thick-wall Lame solution they are the limit of.
//! Pure analytic functions — no meshing, no iteration, no external
//! processes.
//!
//! ## What
//!
//! - **Thin-wall cylinder** ([`cylinder`], [`ThinCylinder`]) — the
//!   membrane hoop stress `sigma_h = p r / t`, the longitudinal stress
//!   `sigma_l = p r / (2 t)`, the maximum in-plane and 3-D shear
//!   stresses, and the von Mises equivalent stress.
//! - **Thin-wall sphere** ([`sphere`], [`ThinSphere`]) — the single
//!   membrane stress `sigma = p r / (2 t)`, equal in all directions,
//!   plus its shear and von Mises values.
//! - **Thick-wall cylinder** ([`thick`], [`ThickCylinder`]) — the Lame
//!   radial and circumferential (hoop) stress at any wall radius, the
//!   peak (bore) hoop stress, the closed-ended axial stress, and the
//!   through-wall maximum shear, for internal and/or external pressure.
//!
//! All stresses come out in the same unit as the supplied pressure (the
//! formulae are dimensionless in the radius ratio), so any consistent
//! unit system works.
//!
//! ## Model
//!
//! Thin-wall (membrane) theory treats a vessel with `r / t > 10` as
//! carrying pressure in near-uniform direct stresses, ignoring the
//! through-thickness gradient. A longitudinal force balance gives the
//! cylinder hoop stress and a transverse one gives the longitudinal
//! stress, so for a cylinder `sigma_h = 2 sigma_l`; the spherical shell
//! carries one isotropic membrane stress equal to the cylinder's
//! longitudinal term, i.e. half its hoop stress.
//!
//! The thick-wall module integrates the exact axisymmetric elasticity
//! equations to the Lame solution `sigma_r = A - B/r^2`,
//! `sigma_t = A + B/r^2`, whose constants are fixed by the surface
//! pressures. The hoop stress peaks at the bore, the radial-plus-hoop sum
//! is constant through the wall, and as the wall thins the bore hoop
//! stress converges to the membrane value — the thin-wall formulae are
//! the asymptotic limit of the Lame field, and the test suite checks each
//! of these relationships against analytic ground truth.
//!
//! See Timoshenko & Goodier (*Theory of Elasticity*), Roark's *Formulas
//! for Stress and Strain* (ch. 13), and Boresi & Schmidt (*Advanced
//! Mechanics of Materials*) for the derivations.
//!
//! ## Honest scope
//!
//! Research / educational grade: textbook closed-form membrane-shell and
//! Lame elasticity for a homogeneous, isotropic, linear-elastic vessel.
//! It is **not** a clinical / medical or production engineering
//! certification tool. There are no code-of-record allowables and no
//! checks for yielding, autofrettage or plastic redistribution, fatigue,
//! weld / nozzle / discontinuity / bending stresses, buckling under
//! external pressure, thermal / residual / shrink-fit stress, or
//! ASME-BPVC / PED / EN-13445 compliance. Membrane results degrade as
//! `r / t` approaches and drops below 10 — use the thick-wall Lame
//! solution there. Do not size real pressure equipment from these numbers
//! without an accredited code calculation.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod cylinder;
pub mod error;
pub mod sphere;
pub mod thick;

pub use cylinder::ThinCylinder;
pub use error::{ErrorCategory, Result, VesselError};
pub use sphere::ThinSphere;
pub use thick::{StressState, ThickCylinder};

#[cfg(test)]
mod tests {
    use super::*;

    /// Loose absolute tolerance for f64 stress comparisons.
    const EPS: f64 = 1e-9;

    /// Cross-module ground truth tying all three vessel types together
    /// for one shared `(p, r, t)`:
    ///
    /// - cylinder hoop = 2 * cylinder longitudinal,
    /// - sphere membrane = half cylinder hoop,
    /// - thick-wall bore hoop ~ cylinder hoop when the wall is thin.
    #[test]
    fn three_vessel_consistency() {
        let p = 6.0;
        let r = 200.0;
        let t = 1.0; // r/t = 200 -> firmly thin-wall

        let cyl = ThinCylinder::new(p, r, t).unwrap();
        let sph = ThinSphere::new(p, r, t).unwrap();
        let thick = ThickCylinder::new(r, r + t, p).unwrap();

        // Cylinder hoop is twice its longitudinal stress.
        assert!((cyl.hoop_stress() - 2.0 * cyl.longitudinal_stress()).abs() < EPS);

        // Sphere membrane is half the cylinder hoop.
        assert!((sph.membrane_stress() - 0.5 * cyl.hoop_stress()).abs() < EPS);

        // Thick-wall bore hoop matches the thin-wall hoop within O(t/r).
        let rel = (thick.max_hoop_stress() - cyl.hoop_stress()).abs() / cyl.hoop_stress();
        assert!(rel < 1e-2, "thin/thick mismatch {rel}");
    }

    #[test]
    fn re_exports_are_constructible() {
        // The public surface is reachable through the crate root.
        let _: VesselError = VesselError::geometry("x");
        let _: ErrorCategory = ErrorCategory::Input;
        let _: Result<ThinCylinder> = ThinCylinder::new(1.0, 1.0, 0.05);
        let _: Result<ThinSphere> = ThinSphere::new(1.0, 1.0, 0.05);
        let s: StressState = ThickCylinder::new(1.0, 2.0, 10.0)
            .unwrap()
            .stress_at(1.5)
            .unwrap();
        assert!((s.radius - 1.5).abs() < EPS);
    }
}
