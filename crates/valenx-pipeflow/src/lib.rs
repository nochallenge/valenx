//! # valenx-pipeflow
//!
//! Internal (pipe / duct) flow hydraulics — the textbook chain that
//! turns a fluid, a pipe and a velocity into a friction factor and the
//! pressure it costs to push the flow down the line.
//!
//! ## What
//!
//! Three small, composable layers:
//!
//! - [`reynolds`] — the Reynolds number `Re = rho v D / mu` and the
//!   classification of the flow into [`reynolds::FlowRegime::Laminar`] /
//!   [`reynolds::FlowRegime::Transitional`] / [`reynolds::FlowRegime::Turbulent`]
//!   at the standard `2300` / `4000` thresholds, plus the inverse
//!   [`reynolds::velocity_for_reynolds`] — the velocity `v = Re mu / (rho D)`
//!   for a target Reynolds number (e.g. the critical transition velocity).
//! - [`friction`] — the Darcy friction factor: the exact laminar
//!   `f = 64/Re`, and the explicit **Haaland** correlation for turbulent
//!   flow (an algebraic stand-in for the implicit Colebrook-White
//!   equation). [`friction::friction_factor`] dispatches on the regime.
//! - [`headloss`] — the Darcy-Weisbach head loss
//!   `h_f = f (L/D) v^2 / (2 g)`, the matching pressure drop
//!   `dP = rho g h_f`, the **wall shear stress** `tau_w = f rho v^2 / 8`
//!   and the **friction velocity** `u* = sqrt(tau_w/rho)`
//!   ([`headloss::wall_shear_stress`], [`headloss::friction_velocity`]),
//!   and an end-to-end [`headloss::solve_pipe`] that runs the whole chain
//!   from physical inputs.
//!
//! ```
//! use valenx_pipeflow::headloss::solve_pipe;
//! use valenx_pipeflow::reynolds::FlowRegime;
//!
//! // Water at 20 C through 100 m of 100 mm commercial-steel pipe at 2 m/s.
//! // rho = 998 kg/m^3, mu = 1.002e-3 Pa.s, roughness eps = 0.046 mm.
//! let r = solve_pipe(998.0, 1.002e-3, 0.1, 100.0, 4.6e-4 / 0.1, 2.0)
//!     .expect("valid pipe-flow inputs");
//!
//! assert_eq!(r.friction.regime, FlowRegime::Turbulent);
//! println!(
//!     "Re = {:.0}, f = {:.4}, head loss = {:.2} m, dP = {:.0} Pa",
//!     r.friction.reynolds,
//!     r.friction.friction_factor,
//!     r.head_loss_m,
//!     r.pressure_drop_pa,
//! );
//! ```
//!
//! ## Model
//!
//! All quantities are SI (m, m/s, kg/m^3, Pa.s, Pa) and the relations
//! are the standard incompressible single-phase pipe-flow correlations:
//!
//! ```text
//!   Re   = rho * v * D / mu                                  (Reynolds)
//!   v    = Re * mu / (rho * D)                  (velocity for a target Re)
//!   f    = 64 / Re                              (laminar, Re < 2300)
//!   1/sqrt(f) = -1.8 log10[ (eps/D / 3.7)^1.11 + 6.9/Re ]    (Haaland)
//!   h_f  = f * (L / D) * v^2 / (2 g)                  (Darcy-Weisbach)
//!   dP   = rho * g * h_f = f * (L/D) * rho * v^2 / 2  (pressure drop)
//!   tau_w = f * rho * v^2 / 8 = dP * D / (4 L)        (wall shear)
//!   u*    = sqrt(tau_w / rho) = v * sqrt(f / 8)       (friction velocity)
//! ```
//!
//! The Haaland correlation (Haaland, 1983) is the explicit algebraic
//! approximation to the implicit Colebrook-White friction law; the crate
//! tests pin it to within a couple of percent of an iteratively-solved
//! Colebrook reference across the turbulent range. Laminar friction is
//! the exact Hagen-Poiseuille result and carries no roughness term.
//!
//! Every fallible function validates its inputs and returns an
//! [`error::PipeFlowError`] (non-positive / non-finite / out-of-range)
//! rather than emitting a silent `NaN` or a negative friction factor.
//!
//! ## Honest scope
//!
//! Research / educational grade. These are textbook closed-form and
//! well-established numerical models — the Reynolds number, the
//! `64/Re` laminar law, the Haaland approximation and the
//! Darcy-Weisbach equation — exactly as they appear in an introductory
//! fluid-mechanics course (White, *Fluid Mechanics*; Munson et al.).
//! This is **not** a clinical/medical tool and **not** a production
//! engineering or piping-design tool. In particular it deliberately does
//! **not** model:
//!
//! - **Minor losses** — fittings, bends, valves, entrances/exits
//!   (the `K`-factor or equivalent-length terms). Only straight-run
//!   skin friction is computed.
//! - **Compressible, two-phase, or non-Newtonian flow** — single-phase
//!   incompressible Newtonian fluid only; the density and viscosity are
//!   constants supplied by the caller.
//! - **The transitional band** `2300 <= Re < 4000` rigorously — no
//!   correlation is reliable there, so [`friction::friction_factor`]
//!   reports the regime and returns the Haaland estimate for the caller
//!   to treat with appropriate caution.
//! - **Developing-flow / entrance-length effects, heat transfer, or
//!   elevation/pump head** — the head loss is the fully-developed
//!   frictional term in isolation, not a complete system energy balance.
//!
//! None of those omissions makes the output meaningless: for fully
//! developed single-phase flow in a straight pipe the Reynolds number,
//! the friction factor and the Darcy-Weisbach loss are the real
//! engineering numbers, validated here against analytic and
//! Colebrook-reference values.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod friction;
pub mod headloss;
pub mod reynolds;

pub use error::{ErrorCategory, PipeFlowError};
pub use friction::{
    friction_factor, haaland_friction_factor, laminar_friction_factor, FrictionResult,
};
pub use headloss::{
    friction_velocity, head_loss, head_loss_g, pressure_drop, solve_pipe, wall_shear_stress,
    PipeFlowResult, G_STANDARD,
};
pub use reynolds::{
    classify_flow, reynolds_number, reynolds_number_kinematic, velocity_for_reynolds,
    velocity_for_reynolds_kinematic, FlowRegime, RE_LAMINAR_UPPER, RE_TURBULENT_LOWER,
};

#[cfg(test)]
mod tests {
    //! Crate-level integration tests exercising the public re-exports as
    //! a user would call them.
    use super::*;

    /// The re-exported top-level functions compose into a coherent
    /// laminar solution.
    #[test]
    fn top_level_laminar_chain_is_consistent() {
        let rho = 900.0;
        let mu = 0.1;
        let d = 0.05;
        let v = 1.0;
        let re = reynolds_number(rho, v, d, mu).unwrap();
        assert!((re - 450.0).abs() < 1e-9);
        assert_eq!(FlowRegime::classify(re), FlowRegime::Laminar);
        let f = laminar_friction_factor(re).unwrap();
        assert!((f - 64.0 / 450.0).abs() < 1e-12);
        let hf = head_loss(f, 10.0, d, v).unwrap();
        assert!(hf > 0.0);
        let dp = pressure_drop(f, 10.0, d, v, rho).unwrap();
        assert!((dp - rho * G_STANDARD * hf).abs() < 1e-6);
    }

    /// `solve_pipe` agrees with the manual step-by-step composition for
    /// a turbulent case.
    #[test]
    fn solve_pipe_matches_manual_turbulent_steps() {
        let rho = 998.0;
        let mu = 1.002e-3;
        let d = 0.1;
        let l = 100.0;
        let eps_rel = 1.0e-4;
        let v = 3.0;

        let bundled = solve_pipe(rho, mu, d, l, eps_rel, v).unwrap();

        let re = reynolds_number(rho, v, d, mu).unwrap();
        let f = friction_factor(re, eps_rel).unwrap();
        let hf = head_loss(f.friction_factor, l, d, v).unwrap();
        let dp = pressure_drop(f.friction_factor, l, d, v, rho).unwrap();

        assert!((bundled.friction.reynolds - re).abs() < 1e-6);
        assert!((bundled.friction.friction_factor - f.friction_factor).abs() < 1e-12);
        assert!((bundled.head_loss_m - hf).abs() < 1e-9);
        assert!((bundled.pressure_drop_pa - dp).abs() < 1e-6);
        assert_eq!(bundled.friction.regime, FlowRegime::Turbulent);
    }

    /// The constants are the expected textbook values.
    #[test]
    fn exported_constants_have_expected_values() {
        assert!((RE_LAMINAR_UPPER - 2300.0).abs() < 1e-9);
        assert!((RE_TURBULENT_LOWER - 4000.0).abs() < 1e-9);
        assert!((G_STANDARD - 9.806_65).abs() < 1e-9);
    }
}
