//! # valenx-bonemech — bone mechanics and damage resistance
//!
//! Closed-form structural mechanics for **bone**, modelled as a
//! homogeneous linear-elastic solid. The crate turns a cross-section and
//! a tissue material into the textbook quantities a load-vs-strength
//! check needs: axial stress, elastic strain, the fracture load, the
//! second moment of area of a hollow shaft, the bending stress it
//! carries, and an apparent-density strength law.
//!
//! ## What
//!
//! - [`Bone`] — a validated cross-section + material
//!   ([`bone`]) with:
//!   - [`Bone::axial_stress_mpa`] — axial stress `sigma = F / A` (MPa),
//!   - [`Bone::strain`] — Hooke's-law strain
//!     `epsilon = sigma / (E * 1000)` (dimensionless),
//!   - [`Bone::fracture_load_n`] — ultimate axial force
//!     `F_ult = sigma_ult * A` (N),
//!   - [`Bone::utilisation`] / [`Bone::safety_factor`] — load-vs-fracture
//!     read-outs.
//! - Cross-section geometry ([`section`]):
//!   - [`second_moment_hollow_circle_mm4`] —
//!     `I = (pi / 64)(D_o^4 - D_i^4)` (mm^4),
//!   - [`section_modulus_mm3`] — `S = I / c` (mm^3),
//!   - [`bending_stress_mpa`] — Euler-Bernoulli flexure
//!     `sigma = M c / I` (MPa),
//!   - [`bending_moment_for_stress`] — its inverse `M = sigma I / c`
//!     (N mm); feeding the ultimate stress gives the bending fracture
//!     moment, the flexural counterpart of the axial fracture load.
//! - Density scaling ([`density`]):
//!   - [`PowerLaw`] / [`strength_from_density`] — the Carter-Hayes
//!     apparent-density power law `sigma_ult ~ rho^exponent`
//!     (default exponent 2).
//!
//! ## Model
//!
//! Bone tissue is treated as **homogeneous, isotropic and
//! linear-elastic** up to a single ultimate stress. A long-bone shaft is
//! idealised as a **hollow circular tube** (cortical wall around a hollow
//! canal), and bending follows **Euler-Bernoulli** beam theory. The
//! cortical defaults ([`Bone::cortical`]) use representative round
//! textbook values — Young's modulus `17.0 GPa`, ultimate stress
//! `150.0 MPa` — and strength rises with apparent density as a power law
//! (`~rho^2`). Three unit facts make the formulae drop out with no
//! conversion factors: `1 MPa = 1 N/mm^2` (so `N / mm^2` is MPa),
//! `1 GPa = 1000 MPa` (so Hooke's law multiplies the modulus by
//! `1000.0`), and `(N mm)(mm) / mm^4 = N/mm^2` (so the flexure formula is
//! MPa). Every public entry point with a physically constrained input
//! validates it and returns a [`BoneError`] rather than producing a
//! nonsensical result.
//!
//! ## Honest scope
//!
//! **Research / educational grade only.** These are idealised, closed-form
//! textbook equations — they are **not** a clinical, medical, diagnostic,
//! or production-engineering tool and **must not** be used to assess real
//! fracture risk in a person or to certify a real structure. Real bone is
//! anisotropic, viscoelastic, tension/compression-asymmetric, remodels
//! over time, and fails by progressive micro-damage and crack growth —
//! none of which is captured by a single modulus, a single ultimate
//! stress, and a power-law density fit. The default material constants
//! are representative round figures, not measurements of any specific
//! specimen, and the density power law is unreliable outside its
//! calibrated range.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod bone;
pub mod density;
pub mod error;
pub mod section;

pub use bone::{Bone, CORTICAL_MODULUS_GPA, CORTICAL_ULTIMATE_STRESS_MPA};
pub use density::{strength_from_density, PowerLaw, DEFAULT_DENSITY_EXPONENT};
pub use error::{BoneError, ErrorCategory, Result};
pub use section::{
    bending_moment_for_stress, bending_stress_mpa, second_moment_hollow_circle_mm4,
    section_modulus_mm3,
};

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end: a cortical femoral-shaft-like tube under a bending
    /// moment. Build the section, compute its second moment of area, the
    /// peak outer-fibre bending stress, the resulting elastic strain, and
    /// confirm the stress sits below the ultimate stress (so the bone is
    /// intact in this linear model).
    #[test]
    fn bending_pipeline_end_to_end() {
        // Tube: outer 27 mm, inner 17 mm (a stand-in long-bone shaft).
        let outer_d = 27.0;
        let inner_d = 17.0;
        let i = second_moment_hollow_circle_mm4(outer_d, inner_d).unwrap();
        assert!(i > 0.0);

        // Peak stress is at the outer fibre, c = D_o / 2.
        let c = outer_d / 2.0;
        let moment_nmm = 60_000.0; // 60 N·m expressed in N·mm.
        let sigma = bending_stress_mpa(moment_nmm, c, i).unwrap();

        // Section-modulus route must agree exactly.
        let s = section_modulus_mm3(i, c).unwrap();
        assert!((sigma - moment_nmm / s).abs() < 1e-9);

        // Strain from that stress via Hooke's law on cortical bone.
        let bone = Bone::cortical(1.0).unwrap();
        let eps = bone.strain(sigma).unwrap();
        // sigma = E * eps must round-trip.
        assert!((eps * bone.elastic_modulus_gpa * 1000.0 - sigma).abs() < 1e-9);

        // The peak stress is below the cortical ultimate stress here, so
        // the section is intact (sigma < sigma_ult).
        assert!(
            sigma < bone.ultimate_stress_mpa,
            "expected intact: sigma {sigma} >= sigma_ult {}",
            bone.ultimate_stress_mpa
        );
    }

    /// End-to-end: a denser cortical cross-section (scaled via the
    /// Carter-Hayes power law) carries a strictly larger ultimate axial
    /// force than a less dense one of the same area.
    #[test]
    fn density_to_fracture_load_end_to_end() {
        let area = 320.0;
        // Anchor the squared law at the cortical defaults.
        let law = PowerLaw::carter_hayes(1.9, CORTICAL_ULTIMATE_STRESS_MPA).unwrap();

        let dense = Bone::new(area, CORTICAL_MODULUS_GPA, law.strength(2.0).unwrap()).unwrap();
        let sparse = Bone::new(area, CORTICAL_MODULUS_GPA, law.strength(1.2).unwrap()).unwrap();

        assert!(
            dense.fracture_load_n() > sparse.fracture_load_n(),
            "denser bone should carry more load: {} vs {}",
            dense.fracture_load_n(),
            sparse.fracture_load_n()
        );
    }
}
