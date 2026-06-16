//! # valenx-popdynamics — population & epidemic dynamics
//!
//! Textbook continuous-time models of how populations and epidemics
//! evolve, each integrated with a hand-rolled fixed-step classical
//! fourth-order Runge-Kutta (RK4) solver. Pure numerics — no external
//! processes, no neural-network weights, two lightweight dependencies
//! (`thiserror` for the error enum, `serde` for the data-model
//! derives).
//!
//! ## What
//!
//! - **Logistic growth** ([`logistic`]) — the single-species Verhulst
//!   model [`Logistic`] with `dN/dt = r N (1 - N/K)`, a closed-form
//!   analytic solution, and an RK4 [`simulate`](Logistic::simulate).
//! - **Predator-prey** ([`lotka_volterra`]) — the two-species
//!   [`LotkaVolterra`] system, its coexistence equilibrium, and the
//!   conserved quantity that pins each trajectory to a closed orbit.
//! - **SIR epidemic** ([`sir`]) — the Kermack-McKendrick [`Sir`]
//!   compartmental model with transmission `beta`, recovery `gamma`,
//!   the basic reproduction number [`R0 = beta/gamma`](Sir::r0), and the
//!   herd-immunity threshold
//!   [`HIT = 1 - 1/R0`](Sir::herd_immunity_threshold).
//! - **RK4 integrator** ([`rk4`]) — the shared fixed-step
//!   [`integrate`] routine (and a single-step [`rk4_step`]) on a
//!   fixed-length state vector.
//! - **Errors** ([`error`]) — a small validated [`PopError`] enum with
//!   stable [`code`](PopError::code) / [`category`](PopError::category)
//!   accessors.
//!
//! ## Model
//!
//! All three systems are autonomous first-order ODEs `dy/dt = f(y)`
//! that share the [`rk4`] integrator:
//!
//! ```text
//! Logistic        dN/dt = r N (1 - N/K)
//!
//! Lotka-Volterra  dx/dt =  alpha x - beta  x y
//!                 dy/dt = -gamma y + delta x y
//!
//! SIR             dS/dt = -beta S I / N
//!                 dI/dt =  beta S I / N - gamma I
//!                 dR/dt =  gamma I,        N = S + I + R.
//! ```
//!
//! Each model's documented invariants are checked against ground truth
//! in its unit tests: logistic converges to `K` (and matches its
//! closed form); `R0 = beta/gamma` and the SIR epidemic grows iff
//! `R0 > 1` while `S + I + R` is conserved; the Lotka-Volterra orbit
//! oscillates with the prey peak preceding the predator peak. The RK4
//! integrator's own fourth-order global accuracy is verified against
//! the analytic solution of `y' = y`.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are the *textbook* deterministic,
//! well-mixed, closed-population models — the ones in an introductory
//! mathematical-biology course. They deliberately omit everything that
//! makes a real model: no demography (births/deaths) in SIR, no latent
//! `E` compartment, no age/spatial/network structure, no stochasticity
//! or observation noise, no parameter inference from data, no
//! environmental carrying-capacity dynamics, and a fixed-step
//! non-adaptive integrator with no stiffness handling. This crate is
//! **NOT a clinical/medical/production engineering tool** — do not use
//! it for patient care, public-health policy, ecological management, or
//! outbreak forecasting.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod logistic;
pub mod lotka_volterra;
pub mod rk4;
pub mod sir;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{ErrorCategory, PopError, Result};
pub use logistic::Logistic;
pub use lotka_volterra::LotkaVolterra;
pub use rk4::{integrate, rk4_step, Sample, State, MAX_STEPS};
pub use sir::{Sir, SirState};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reexports_are_usable() {
        // Smoke test that the headline re-exports compose end to end.
        let logi = Logistic::new(0.5, 100.0).unwrap();
        let lt = logi.simulate(1.0, 30.0, 0.05).unwrap();
        assert!((lt.last().unwrap().y[0] - 100.0).abs() < 1.0);

        let sir = Sir::new(0.6, 0.2).unwrap();
        assert!((sir.r0() - 3.0).abs() < 1e-12);

        let lv = LotkaVolterra::new(1.0, 0.1, 1.5, 0.075).unwrap();
        assert_eq!(lv.equilibrium().len(), 2);

        // The error enum is reachable and typed.
        let e: PopError = Logistic::new(0.1, -1.0).unwrap_err();
        assert_eq!(e.category(), ErrorCategory::Input);
    }

    #[test]
    fn rk4_reexport_runs() {
        // y' = -y from the top-level re-export; exact y(1) = e^{-1}.
        let traj = integrate(|_t, s: &State<1>| [-s[0]], [1.0], 0.0, 1.0, 0.01).unwrap();
        let expected = (-1.0_f64).exp();
        assert!((traj.last().unwrap().y[0] - expected).abs() < 1e-6);
        let one = rk4_step(0.0, [1.0], 0.0, |_t, s: &State<1>| [-s[0]]);
        // A zero-size step is a no-op.
        assert!((one[0] - 1.0).abs() < 1e-12);
        // The re-exported MAX_STEPS ceiling is the one the integrator
        // enforces: a window needing more than MAX_STEPS steps is rejected.
        let over = integrate(
            |_t, s: &State<1>| [-s[0]],
            [1.0],
            0.0,
            (MAX_STEPS as f64) + 10.0,
            1.0,
        );
        assert_eq!(over.unwrap_err().code(), "popdynamics.too_many_steps");
    }
}
