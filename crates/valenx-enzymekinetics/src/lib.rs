//! # valenx-enzymekinetics — closed-form enzyme-kinetics rate laws
//!
//! A small, dependency-light library of the standard steady-state
//! enzyme-kinetics models, evaluated as exact algebraic formulas on
//! validated parameters.
//!
//! ## What
//!
//! - **Michaelis-Menten** ([`MichaelisMenten`], module
//!   [`michaelis_menten`]) — the single-substrate rate law
//!   `v = Vmax * S / (Km + S)`, the dimensionless saturation
//!   `S / (Km + S)`, and the closed-form **integrated** equation
//!   `Vmax*t = Km*ln(s0/s) + (s0 - s)`
//!   ([`MichaelisMenten::time_to_deplete`]) for the substrate-depletion
//!   progress curve.
//! - **Reversible inhibition** ([`inhibition`]) — [`Competitive`]
//!   (apparent `Km` scaled by `1 + I/Ki`), [`Noncompetitive`] (apparent
//!   `Vmax` divided by `1 + I/Ki`), [`Uncompetitive`] (both scaled by
//!   `1 + I/Ki'`), and the general [`Mixed`] case with independent
//!   free-enzyme and `ES`-complex inhibition constants.
//! - **Cooperativity** ([`Hill`], module [`hill`]) — the Hill equation
//!   `v = Vmax * S^n / (K^n + S^n)`, which collapses to Michaelis-Menten
//!   when `n = 1`.
//! - **Errors** ([`error`]) — a [`KineticsError`] type (via
//!   [`thiserror`]) returned by every validated constructor and rate
//!   function.
//!
//! Every parameter struct validates its inputs eagerly in its `new`
//! constructor, so once a model object exists, evaluating a rate on a
//! valid concentration is total. Concentrations and velocities are unit-
//! agnostic: pick a consistent set (e.g. µM for `S`, `Km`, `K`, `Ki`;
//! µM·s⁻¹ for `v` and `Vmax`) and the formulas hold.
//!
//! ## Model
//!
//! All four laws share one algebraic backbone — a rectangular hyperbola
//! (or its `n`-th-power Hill generalisation) in the substrate
//! concentration:
//!
//! ```text
//! Michaelis-Menten   v = Vmax · S / (Km + S)
//! competitive        v = Vmax · S / (Km·(1 + I/Ki) + S)
//! noncompetitive     v = (Vmax / (1 + I/Ki)) · S / (Km + S)
//! uncompetitive      v = (Vmax / (1 + I/Ki')) · S / (Km/(1 + I/Ki') + S)
//! mixed              v = Vmax · S / ((1 + I/Ki)·Km + (1 + I/Ki')·S)
//! Hill               v = Vmax · S^n / (K^n + S^n)
//! ```
//!
//! Each inhibition mode is exactly Michaelis-Menten re-evaluated with
//! *apparent* parameters; the diagnostic fingerprints (competitive moves
//! `Km` only, noncompetitive moves `Vmax` only, uncompetitive moves both
//! by the same factor) follow directly and are pinned by the unit tests,
//! along with `v(Km) = Vmax/2`, saturation `v -> Vmax` as `S -> ∞`,
//! strict monotonicity in `S`, and the `n = 1` Hill-to-Michaelis-Menten
//! reduction.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form /
//! well-established numerical models — the initial-rate Briggs-Haldane
//! steady-state laws found in any biochemistry text — and nothing more.
//! In particular this crate is **not** a clinical, medical, diagnostic,
//! or production engineering tool: do not use it for drug dosing, assay
//! certification, or any decision that affects health or safety.
//!
//! Deliberate limitations:
//!
//! - **Initial rates, plus the closed-form integrated MM curve.** The
//!   rate laws are steady-state *initial* velocities; the one progress
//!   curve provided is the exact analytic integrated Michaelis-Menten
//!   equation ([`MichaelisMenten::time_to_deplete`]). There is no
//!   *numerical* time integration, no product accumulation, and no
//!   progress-curve *fitting*.
//! - **Idealised mechanism.** Single substrate, irreversible product
//!   release, free-ligand ≈ total-ligand (the standard `[S], [I] ≫ [E]`
//!   assumption), and a phenomenological Hill coefficient rather than an
//!   explicit multi-site binding partition function.
//! - **No parameter estimation.** The crate *evaluates* rate laws from
//!   known parameters; it does not fit `Vmax`, `Km`, `Ki`, `K`, or `n`
//!   from measured data (no Lineweaver-Burk / Eadie-Hofstee /
//!   nonlinear-regression machinery).
//! - **No allosteric / multi-substrate models.** No Monod-Wyman-Changeux
//!   or Koshland-Némethy-Filmer formalisms, no ping-pong or
//!   ordered/random bi-substrate kinetics, no pH / temperature
//!   dependence.
//!
//! The crate writes no files and runs no external processes; callers
//! persist any results themselves (the parameter structs derive
//! [`serde`] for that).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod hill;
pub mod inhibition;
pub mod michaelis_menten;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{KineticsError, Result};
pub use hill::Hill;
pub use inhibition::{Competitive, Mixed, Noncompetitive, Uncompetitive};
pub use michaelis_menten::MichaelisMenten;

#[cfg(test)]
mod tests {
    //! Cross-module integration checks that exercise the re-exported API
    //! surface together — the per-model analytic properties live in each
    //! module's own `tests`.

    use super::*;

    /// Absolute-difference tolerance for the float comparisons below.
    const EPS: f64 = 1e-9;

    /// The three single-inhibitor modes and the Hill law all agree with
    /// bare Michaelis-Menten in their respective no-effect limits, so a
    /// caller can reach the baseline through any of them.
    #[test]
    fn all_models_share_the_michaelis_menten_baseline() {
        let mm = MichaelisMenten::new(4.0, 2.0).expect("valid");
        let comp = Competitive::new(mm, 1.0).expect("valid");
        let noncomp = Noncompetitive::new(mm, 1.0).expect("valid");
        let uncomp = Uncompetitive::new(mm, 1.0).expect("valid");
        let mixed = Mixed::new(mm, 1.0, 1.0).expect("valid");
        let hill = Hill::new(mm.vmax(), mm.km(), 1.0).expect("valid");

        for &s in &[0.0, 0.5, 2.0, 10.0, 100.0] {
            let base = mm.velocity(s).expect("valid");
            // I = 0 → every inhibition mode is the bare law.
            assert!((comp.velocity(s, 0.0).expect("v") - base).abs() < EPS);
            assert!((noncomp.velocity(s, 0.0).expect("v") - base).abs() < EPS);
            assert!((uncomp.velocity(s, 0.0).expect("v") - base).abs() < EPS);
            assert!((mixed.velocity(s, 0.0).expect("v") - base).abs() < EPS);
            // n = 1 → Hill is the bare law (K = Km).
            assert!((hill.velocity(s).expect("v") - base).abs() < EPS);
        }
    }

    /// The error type is a real `std::error::Error` and surfaces the
    /// offending parameter name, so callers can log it directly.
    #[test]
    fn errors_propagate_as_std_error() {
        let err = MichaelisMenten::new(1.0, -1.0).unwrap_err();
        let boxed: Box<dyn std::error::Error> = Box::new(err);
        assert!(boxed.to_string().contains("km"), "got: {boxed}");
    }
}
