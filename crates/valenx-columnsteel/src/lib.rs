//! # valenx-columnsteel — Euler-Johnson steel-column buckling
//!
//! ## What
//!
//! A small, dependency-light calculator for the compressive strength of
//! axially-loaded steel columns by the classic **Euler-Johnson** (AISC
//! Allowable Stress Design) method. Given a column's Young's modulus,
//! yield stress, and slenderness ratio (or its raw geometry) it returns:
//!
//! 1. the **slenderness ratio** `λ = K L / r`
//!    ([`slenderness_ratio`], [`Column::slenderness`]);
//! 2. the **column-slenderness transition** `Cc = sqrt(2 π² E / Fy)`
//!    ([`transition_slenderness`], [`Column::cc`]);
//! 3. which **regime** governs — elastic Euler (long) or inelastic
//!    Johnson (short) — ([`Column::regime`], [`Regime`]);
//! 4. the **critical buckling stress** `Fcr`
//!    ([`Column::critical_stress`]); and
//! 5. the **allowable compressive stress** `Fa = Fcr / FS` and allowable
//!    axial load ([`Column::allowable_stress`], [`Column::allowable_load`]).
//!
//! ## Model
//!
//! The strength curve has two branches that meet tangentially at the
//! transition slenderness `Cc`:
//!
//! 1. **Long columns** (`λ ≥ Cc`) buckle elastically; the **Euler**
//!    critical stress governs, `Fcr = π² E / λ²`. As `λ → ∞`, `Fcr → 0`.
//! 2. **Short / intermediate columns** (`λ < Cc`) yield before they reach
//!    the Euler load; the **Johnson parabola** governs,
//!    `Fcr = Fy [1 − λ² / (2 Cc²)]`. As `λ → 0`, `Fcr → Fy` (the squash
//!    load).
//!
//! The transition `Cc = sqrt(2 π² E / Fy)` is exactly the slenderness at
//! which the Euler stress equals `Fy / 2`; both branches pass through
//! `(Cc, Fy/2)` with matching slope, so the assembled `Fcr(λ)` curve is
//! continuous and monotonically decreasing. The allowable stress divides
//! `Fcr` by a factor of safety: the AISC-ASD **variable** safety factor
//! ramps from `5/3` at `λ = 0` to `23/12` at `λ = Cc` and stays at
//! `23/12` in the Euler regime ([`Column::factor_of_safety_aisc`]), or a
//! caller may supply a constant factor ([`Column::allowable_stress_with`]).
//!
//! Stresses are returned in whatever unit the modulus and yield stress
//! are given in (psi in, psi out; MPa in, MPa out); lengths enter only
//! through the dimensionless ratio `K L / r`, so any consistent length
//! unit works.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are the **textbook closed-form**
//! Euler-Johnson column equations from the AISC ASD tradition; they are
//! **NOT a clinical/medical or production engineering tool** and must
//! not be used as the basis for a real structural design. The model is
//! deliberately idealized: a prismatic, doubly-symmetric, concentrically
//! loaded member that buckles in pure flexure about a single principal
//! axis. It ignores local / flange-and-web plate buckling, flexural-
//! torsional and lateral-torsional buckling, residual stresses and
//! initial out-of-straightness beyond what the empirical safety factor
//! implicitly covers, combined axial-plus-bending interaction, slender
//! (non-compact) elements, and any specific edition's detailed code
//! provisions. The caller supplies the effective-length factor `K`; this
//! crate does not derive it from end conditions or perform a frame
//! buckling (alignment-chart) analysis. Real column design must follow
//! the governing building code and the judgment of a licensed engineer.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod column;
pub mod error;

pub use column::{slenderness_ratio, transition_slenderness, Column, Regime};
pub use error::{ColumnError, ErrorCategory, Result};

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end smoke test exercising the public re-exports: build a
    /// column from geometry, classify it, and derive an allowable load,
    /// confirming the crate hangs together through its public API.
    #[test]
    fn end_to_end_through_public_api() {
        // A36 steel, K=1, L=120 in, r=1.5 in -> lambda = 80 (short).
        let col = Column::from_geometry(29_000_000.0, 36_000.0, 1.0, 120.0, 1.5).unwrap();
        assert_eq!(col.regime(), Regime::Johnson);

        // Critical stress is below yield and above the allowable.
        let fcr = col.critical_stress();
        assert!(fcr <= col.yield_stress());
        assert!(col.allowable_stress() < fcr);

        // Allowable load for a 12 in^2 section is positive and equals
        // Fa * A.
        let p = col.allowable_load(12.0).unwrap();
        assert!((p - col.allowable_stress() * 12.0).abs() < 1e-6);
        assert!(p > 0.0);

        // The free functions agree with the accessors.
        let lambda = slenderness_ratio(1.0, 120.0, 1.5).unwrap();
        assert!((lambda - col.slenderness()).abs() < 1e-9);
        let cc = transition_slenderness(29_000_000.0, 36_000.0).unwrap();
        assert!((cc - col.cc()).abs() < 1e-9);
    }

    /// A genuinely slender column trips the Euler branch through the
    /// public surface.
    #[test]
    fn slender_column_is_euler_through_public_api() {
        let col = Column::new(29_000_000.0, 36_000.0, 180.0).unwrap();
        assert_eq!(col.regime(), Regime::Euler);
        assert_eq!(col.regime().as_str(), "euler");
        let fa = col.allowable_stress();
        assert!(fa > 0.0 && fa < col.critical_stress());
        // Surface an error variant through the re-exported type.
        let err: ColumnError = col.johnson_stress().unwrap_err();
        assert_eq!(err.category(), ErrorCategory::Domain);
    }
}
