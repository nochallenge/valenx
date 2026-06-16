//! # valenx-weld — weld static strength
//!
//! Closed-form static strength of the two everyday weld types: the
//! equal-leg fillet weld and the complete-joint-penetration (CJP) butt
//! weld. Pure analytic textbook formulas — no external processes, no
//! solver, no I/O.
//!
//! ## What
//!
//! - [`fillet`] — equal-leg fillet welds. The effective throat
//!   `0.707 * leg` ([`fillet::FILLET_THROAT_FACTOR`]), the average
//!   throat shear stress `tau = F / (throat * length)`, the
//!   allowable load `F_allow = tau_allow * throat * length`, and the
//!   weld length required to carry a load at an allowable stress
//!   ([`fillet::required_length`]) — both as a [`fillet::FilletWeld`]
//!   value and as free functions.
//! - [`butt`] — complete-joint-penetration butt welds. The direct
//!   normal stress `sigma = F / (thickness * length)`, the matching
//!   allowable load, and the required weld length
//!   ([`butt::required_length`]), as a [`butt::ButtWeld`] value and as
//!   free functions.
//! - [`error`] — the [`WeldError`] taxonomy and the shared
//!   strictly-positive-finite input gate.
//!
//! ## Model
//!
//! A fillet weld is taken to fail in shear on its effective throat. For
//! an equal-leg, 45-degree fillet the throat is the altitude of an
//! isosceles right triangle, `leg * cos(45 deg) = leg / sqrt(2)`, which
//! the AWS D1.1 / AISC tables round to `0.707 * leg`; this crate uses
//! that tabulated `0.707` so a result matches a by-hand code-table
//! calculation. The force is carried as a uniform average shear stress
//! over `throat * length`.
//!
//! A CJP butt weld fuses the full plate thickness and is treated as
//! continuous with the base metal, so a force normal to the weld
//! produces a uniform direct stress over `thickness * length` (use the
//! thinner connected part's thickness).
//!
//! Every quantity is unit-agnostic but the inputs must be consistent —
//! e.g. lengths in millimetres, force in newtons and stress in
//! `N/mm^2 = MPa`.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are textbook closed-form models
//! carrying the standard first-pass simplifying assumptions: equal-leg
//! 45-degree fillets, full-penetration (effective-throat-equals-
//! thickness) butt welds, uniform stress over the throat / section,
//! concentric in-plane (or normal) loading, the base metal not
//! governing, and no fatigue, residual stress, eccentricity, weld-group
//! geometry, or the AISC directional-strength increase for transverse
//! fillets. This is **NOT** a clinical/medical/production engineering
//! tool and is no substitute for a qualified welding engineer or a
//! governing structural-code check.
//!
//! ## Example
//!
//! ```
//! use valenx_weld::fillet::FilletWeld;
//!
//! // A 6 mm equal-leg fillet, 100 mm long, carrying 50 kN.
//! let weld = FilletWeld::new(6.0, 100.0).unwrap();
//! assert!((weld.throat() - 4.242).abs() < 1e-9);
//! let tau = weld.shear_stress(50_000.0).unwrap(); // N/mm^2 = MPa
//! assert!((tau - 117.868_929_75).abs() < 1e-6);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod butt;
pub mod error;
pub mod fillet;

pub use butt::ButtWeld;
pub use error::{ErrorCategory, Result, WeldError};
pub use fillet::{FilletWeld, FILLET_THROAT_FACTOR};

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // Cross-module sanity: a fillet and a butt weld of the SAME effective
    // load-carrying area carry the SAME force at the same stress, because
    // both reduce to stress = force / area.
    #[test]
    fn fillet_and_butt_agree_for_equal_area() {
        // Fillet: throat 0.707*10 = 7.07, length 100 -> area 707.
        let fillet = FilletWeld::new(10.0, 100.0).unwrap();
        // Butt: thickness 7.07, length 100 -> area 707.
        let butt = ButtWeld::new(fillet.throat(), 100.0).unwrap();
        assert!((fillet.throat_area() - butt.area()).abs() < EPS);

        let force = 70_700.0;
        let tau = fillet.shear_stress(force).unwrap();
        let sigma = butt.normal_stress(force).unwrap();
        assert!((tau - sigma).abs() < EPS, "tau {tau} sigma {sigma}");
        assert!((tau - 100.0).abs() < 1e-9, "got: {tau}");
    }

    // The re-exported constant is exactly the tabulated 0.707.
    #[test]
    fn reexported_throat_factor_is_tabulated_value() {
        assert!((FILLET_THROAT_FACTOR - 0.707).abs() < EPS);
    }

    // The re-exported error surface round-trips category/code.
    #[test]
    fn reexported_error_surface() {
        let err: WeldError = WeldError::degenerate("x");
        assert_eq!(err.category(), ErrorCategory::Geometry);
        let _r: Result<()> = Err(err);
    }
}
