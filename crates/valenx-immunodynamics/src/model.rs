//! The target-cell-limited within-host viral-dynamics model.
//!
//! ## State
//!
//! Three populations evolve in continuous time:
//!
//! - `T` — uninfected **target cells** (cells the virus can infect).
//! - `I` — productively **infected cells** (releasing virions).
//! - `V` — free **virions** (the viral load).
//!
//! ## Equations
//!
//! The standard three-compartment model (Nowak & May 2000; Perelson
//! 2002) is
//!
//! ```text
//! dT/dt = -beta * T * V
//! dI/dt =  beta * T * V - delta * I
//! dV/dt =  p * I        - c * V
//! ```
//!
//! with non-negative rate constants
//!
//! - `beta` — mass-action **infection rate** (per virion per target
//!   cell per unit time).
//! - `delta` — per-capita **death / clearance rate** of infected cells.
//! - `p` — per-infected-cell virion **production rate** (burst rate).
//! - `c` — per-virion **clearance rate** of free virus.
//!
//! This is the *target-cell-limited* form: there is no target-cell
//! replenishment (no `lambda - d*T` source term), so the infection is
//! ultimately bounded by the finite initial pool `T0` — target cells
//! deplete, the effective reproductive number falls below one, and the
//! viral load peaks and then declines. That depletion is exactly the
//! mechanism this crate is built to demonstrate.
//!
//! ## Basic reproductive number
//!
//! Linearising about the infection-free equilibrium `(T0, 0, 0)` gives
//! the basic reproductive number
//!
//! ```text
//! R0 = beta * T0 * p / (delta * c)
//! ```
//!
//! the expected number of secondary infected cells produced by one
//! infected cell when target cells are still abundant. `R0 > 1` is the
//! threshold for an initial exponential outbreak.

use serde::{Deserialize, Serialize};

use crate::error::{ImmunoError, Result};

/// Rate constants of the target-cell-limited model.
///
/// All four are non-negative; construct via [`Parameters::new`] (which
/// validates) rather than the struct literal so out-of-domain values
/// are rejected at the boundary.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Parameters {
    /// Mass-action infection rate `beta` (per virion per target cell
    /// per unit time). Non-negative.
    pub beta: f64,
    /// Per-capita death / clearance rate `delta` of infected cells.
    /// Non-negative.
    pub delta: f64,
    /// Per-infected-cell virion production (burst) rate `p`.
    /// Non-negative.
    pub p: f64,
    /// Per-virion clearance rate `c` of free virus. Non-negative.
    pub c: f64,
}

impl Parameters {
    /// Validated constructor. Every rate must be a finite, non-negative
    /// number; otherwise an [`ImmunoError::Invalid`] naming the offending
    /// parameter is returned.
    pub fn new(beta: f64, delta: f64, p: f64, c: f64) -> Result<Self> {
        check_rate("beta", beta)?;
        check_rate("delta", delta)?;
        check_rate("p", p)?;
        check_rate("c", c)?;
        Ok(Parameters { beta, delta, p, c })
    }

    /// Basic reproductive number `R0 = beta * T0 * p / (delta * c)`.
    ///
    /// The expected number of secondary infected cells produced by a
    /// single infected cell while target cells are still at their
    /// initial abundance `t0`. `R0 > 1` is the threshold for an initial
    /// exponential outbreak.
    ///
    /// Errors if `t0` is negative or non-finite, or if the denominator
    /// `delta * c` is zero (an infected cell or a virion that is never
    /// cleared makes `R0` undefined / unbounded in this formula).
    pub fn r0(&self, t0: f64) -> Result<f64> {
        if !t0.is_finite() || t0 < 0.0 {
            return Err(ImmunoError::invalid(
                "t0",
                "initial target-cell count must be finite and non-negative",
            ));
        }
        let denom = self.delta * self.c;
        if denom <= 0.0 {
            return Err(ImmunoError::invalid(
                "delta*c",
                "R0 is undefined when the infected-cell or virion clearance rate is zero",
            ));
        }
        Ok(self.beta * t0 * self.p / denom)
    }
}

/// Validate a single rate constant: finite and non-negative.
fn check_rate(what: &'static str, value: f64) -> Result<()> {
    if !value.is_finite() {
        return Err(ImmunoError::invalid(what, "must be a finite number"));
    }
    if value < 0.0 {
        return Err(ImmunoError::invalid(what, "must be non-negative"));
    }
    Ok(())
}

/// A model state: the three populations `(T, I, V)` at one instant.
///
/// In a physical infection all three are non-negative; [`State::new`]
/// enforces that at construction, and the integrator
/// ([`crate::integrate`]) keeps the trajectory non-negative.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct State {
    /// Uninfected target cells `T`.
    pub target: f64,
    /// Productively infected cells `I`.
    pub infected: f64,
    /// Free virions `V` (the viral load).
    pub virus: f64,
}

impl State {
    /// Validated constructor. Each population must be finite and
    /// non-negative.
    pub fn new(target: f64, infected: f64, virus: f64) -> Result<Self> {
        check_pop("target", target)?;
        check_pop("infected", infected)?;
        check_pop("virus", virus)?;
        Ok(State {
            target,
            infected,
            virus,
        })
    }

    /// Whether every population is non-negative (within a small
    /// tolerance to absorb floating-point round-off).
    pub fn is_non_negative(&self, tol: f64) -> bool {
        self.target >= -tol && self.infected >= -tol && self.virus >= -tol
    }

    /// Whether every component is finite.
    pub fn is_finite(&self) -> bool {
        self.target.is_finite() && self.infected.is_finite() && self.virus.is_finite()
    }

    /// The time derivative `(dT/dt, dI/dt, dV/dt)` of this state under
    /// `params`, i.e. the right-hand side of the model ODE:
    ///
    /// ```text
    /// dT/dt = -beta * T * V
    /// dI/dt =  beta * T * V - delta * I
    /// dV/dt =  p * I        - c * V
    /// ```
    pub fn derivative(&self, params: &Parameters) -> State {
        let infection = params.beta * self.target * self.virus;
        State {
            target: -infection,
            infected: infection - params.delta * self.infected,
            virus: params.p * self.infected - params.c * self.virus,
        }
    }
}

/// Validate a single population value: finite and non-negative.
fn check_pop(what: &'static str, value: f64) -> Result<()> {
    if !value.is_finite() {
        return Err(ImmunoError::invalid(what, "must be a finite number"));
    }
    if value < 0.0 {
        return Err(ImmunoError::invalid(what, "must be non-negative"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parameters_reject_negative_rates() {
        assert!(Parameters::new(-1.0, 1.0, 1.0, 1.0).is_err());
        assert!(Parameters::new(1.0, -1.0, 1.0, 1.0).is_err());
        assert!(Parameters::new(1.0, 1.0, -1.0, 1.0).is_err());
        assert!(Parameters::new(1.0, 1.0, 1.0, -1.0).is_err());
        assert!(Parameters::new(f64::NAN, 1.0, 1.0, 1.0).is_err());
        assert!(Parameters::new(1.0, 1.0, 1.0, 1.0).is_ok());
    }

    #[test]
    fn state_rejects_negative_populations() {
        assert!(State::new(-1.0, 0.0, 0.0).is_err());
        assert!(State::new(0.0, -1.0, 0.0).is_err());
        assert!(State::new(0.0, 0.0, -1.0).is_err());
        assert!(State::new(f64::INFINITY, 0.0, 0.0).is_err());
        assert!(State::new(1.0, 0.0, 0.0).is_ok());
    }

    #[test]
    fn r0_matches_closed_form() {
        // beta=2e-7, delta=1, p=100, c=23, T0=4e8 (representative
        // influenza-like numbers). R0 = beta*T0*p/(delta*c).
        let beta = 2e-7;
        let delta = 1.0;
        let p = 100.0;
        let c = 23.0;
        let t0 = 4e8;
        let params = Parameters::new(beta, delta, p, c).unwrap();
        let expect = beta * t0 * p / (delta * c);
        let got = params.r0(t0).unwrap();
        assert!(
            (got - expect).abs() < 1e-9 * expect.abs(),
            "got {got}, want {expect}"
        );
    }

    #[test]
    fn r0_handles_simple_round_numbers() {
        // Choose values that make R0 exactly 4: beta*T0*p = 4 * delta*c.
        // beta=1, T0=2, p=2, delta=1, c=1 -> 1*2*2/(1*1) = 4.
        let params = Parameters::new(1.0, 1.0, 2.0, 1.0).unwrap();
        let r0 = params.r0(2.0).unwrap();
        assert!((r0 - 4.0).abs() < 1e-12, "got {r0}");
    }

    #[test]
    fn r0_rejects_zero_clearance_and_bad_t0() {
        let params = Parameters::new(1.0, 0.0, 1.0, 1.0).unwrap();
        assert!(params.r0(1.0).is_err(), "delta=0 makes R0 undefined");
        let params = Parameters::new(1.0, 1.0, 1.0, 0.0).unwrap();
        assert!(params.r0(1.0).is_err(), "c=0 makes R0 undefined");
        let params = Parameters::new(1.0, 1.0, 1.0, 1.0).unwrap();
        assert!(params.r0(-1.0).is_err(), "negative T0 rejected");
    }

    #[test]
    fn derivative_matches_hand_computation() {
        // beta=0.1, delta=0.5, p=10, c=2; state T=3, I=4, V=5.
        let params = Parameters::new(0.1, 0.5, 10.0, 2.0).unwrap();
        let s = State::new(3.0, 4.0, 5.0).unwrap();
        let d = s.derivative(&params);
        // infection = 0.1*3*5 = 1.5
        // dT = -1.5
        // dI = 1.5 - 0.5*4 = 1.5 - 2.0 = -0.5
        // dV = 10*4 - 2*5 = 40 - 10 = 30
        assert!((d.target - (-1.5)).abs() < 1e-12, "dT = {}", d.target);
        assert!((d.infected - (-0.5)).abs() < 1e-12, "dI = {}", d.infected);
        assert!((d.virus - 30.0).abs() < 1e-12, "dV = {}", d.virus);
    }

    #[test]
    fn target_cells_never_increase() {
        // dT/dt = -beta*T*V <= 0 for non-negative populations, so the
        // target-cell derivative is never positive.
        let params = Parameters::new(0.3, 1.0, 50.0, 5.0).unwrap();
        let s = State::new(1000.0, 10.0, 20.0).unwrap();
        assert!(s.derivative(&params).target <= 0.0);
    }
}
