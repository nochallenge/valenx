//! The classic SIR compartmental epidemic model.
//!
//! ## Model
//!
//! The Kermack-McKendrick SIR model splits a fixed, closed population
//! into three compartments — Susceptible `S`, Infectious `I`,
//! Recovered (and immune) `R` — and moves people `S -> I -> R`:
//!
//! ```text
//! dS/dt = -beta S I / N
//! dI/dt =  beta S I / N - gamma I
//! dR/dt =  gamma I
//! ```
//!
//! with `N = S + I + R` the (conserved) total population, `beta` the
//! transmission rate, and `gamma` the recovery rate (`1/gamma` is the
//! mean infectious period). This is the *frequency-dependent*
//! (proportionate-mixing) form, in which the force of infection scales
//! with the *fraction* infectious `I/N`.
//!
//! ## Key quantities
//!
//! - The **basic reproduction number** `R0 = beta / gamma`: the
//!   expected number of secondary infections one infectious individual
//!   causes in a fully susceptible population.
//! - Since `dI/dt = I (beta S/N - gamma)`, an epidemic *grows* at the
//!   start (`S ~ N`) iff `beta - gamma > 0`, i.e. iff `R0 > 1`. The
//!   [`tests`](self) verify both the `R0` identity and this
//!   grows-iff-`R0 > 1` threshold, and that `S + I + R` is conserved
//!   by the integrator.

use crate::error::{PopError, Result};
use crate::rk4::{integrate, Sample};
use serde::{Deserialize, Serialize};

/// Parameters of the SIR model.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Sir {
    /// Transmission rate `beta` (per unit time). The rate at which an
    /// infectious-susceptible contact transmits, times the contact
    /// rate.
    pub beta: f64,
    /// Recovery rate `gamma` (per unit time). The reciprocal `1/gamma`
    /// is the mean infectious period.
    pub gamma: f64,
}

/// An SIR state: the three compartment sizes.
///
/// Stored as a struct for ergonomics; [`SirState::as_vec`] /
/// [`SirState::from_vec`] convert to and from the `[S, I, R]` array the
/// RK4 integrator operates on.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SirState {
    /// Susceptible count.
    pub s: f64,
    /// Infectious count.
    pub i: f64,
    /// Recovered (immune) count.
    pub r: f64,
}

impl SirState {
    /// Build a state from explicit compartment counts.
    pub fn new(s: f64, i: f64, r: f64) -> Self {
        SirState { s, i, r }
    }

    /// Total population `N = S + I + R`.
    pub fn total(&self) -> f64 {
        self.s + self.i + self.r
    }

    /// Convert to the `[S, I, R]` array form the integrator uses.
    pub fn as_vec(&self) -> [f64; 3] {
        [self.s, self.i, self.r]
    }

    /// Build from an `[S, I, R]` array.
    pub fn from_vec(v: [f64; 3]) -> Self {
        SirState {
            s: v[0],
            i: v[1],
            r: v[2],
        }
    }
}

impl Sir {
    /// Construct a validated SIR model.
    ///
    /// # Errors
    ///
    /// [`PopError::Invalid`] if `beta < 0`, `gamma <= 0`, or either is
    /// not finite. `gamma` must be strictly positive (it is the
    /// denominator of `R0` and of the mean infectious period).
    pub fn new(beta: f64, gamma: f64) -> Result<Self> {
        if !beta.is_finite() || beta < 0.0 {
            return Err(PopError::invalid(
                "beta",
                "transmission rate must be finite and non-negative",
            ));
        }
        if !gamma.is_finite() || gamma <= 0.0 {
            return Err(PopError::invalid(
                "gamma",
                "recovery rate must be finite and strictly positive",
            ));
        }
        Ok(Sir { beta, gamma })
    }

    /// The basic reproduction number `R0 = beta / gamma`.
    pub fn r0(&self) -> f64 {
        self.beta / self.gamma
    }

    /// The three time-derivatives `[dS/dt, dI/dt, dR/dt]` at a given
    /// state, using the conserved total `N = S + I + R`.
    ///
    /// If the supplied state has zero total population the force of
    /// infection is taken as zero (avoids `0/0`).
    pub fn rate(&self, state: &SirState) -> [f64; 3] {
        let n = state.total();
        let force = if n > 0.0 {
            self.beta * state.s * state.i / n
        } else {
            0.0
        };
        let recover = self.gamma * state.i;
        [-force, force - recover, recover]
    }

    /// Integrate the epidemic from an initial state to `t_end` with RK4
    /// step `dt`, returning the `[S, I, R]` trajectory.
    ///
    /// # Errors
    ///
    /// [`PopError::Invalid`] if any initial compartment is negative or
    /// non-finite, plus any error from the underlying [`integrate`].
    pub fn simulate(&self, initial: SirState, t_end: f64, dt: f64) -> Result<Vec<Sample<3>>> {
        for (name, v) in [("s0", initial.s), ("i0", initial.i), ("r0_init", initial.r)] {
            if !v.is_finite() || v < 0.0 {
                return Err(PopError::invalid(
                    name,
                    "initial compartment must be finite and non-negative",
                ));
            }
        }
        let model = *self;
        integrate(
            move |_t, y| model.rate(&SirState::from_vec(*y)),
            initial.as_vec(),
            0.0,
            t_end,
            dt,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r0_is_beta_over_gamma() {
        // VALIDATE: R0 = beta/gamma.
        let m = Sir::new(0.6, 0.2).unwrap();
        assert!((m.r0() - 3.0).abs() < 1e-12, "R0={}", m.r0());
        let m = Sir::new(0.25, 0.5).unwrap();
        assert!((m.r0() - 0.5).abs() < 1e-12, "R0={}", m.r0());
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(Sir::new(-0.1, 0.2).is_err());
        assert!(Sir::new(0.6, 0.0).is_err());
        assert!(Sir::new(0.6, -0.2).is_err());
        assert!(Sir::new(f64::INFINITY, 0.2).is_err());
    }

    #[test]
    fn conserves_total_population() {
        // VALIDATE: SIR conserves S + I + R.
        let m = Sir::new(0.5, 0.2).unwrap();
        let init = SirState::new(990.0, 10.0, 0.0);
        let n0 = init.total();
        let traj = m.simulate(init, 60.0, 0.01).unwrap();
        for s in &traj {
            let st = SirState::from_vec(s.y);
            assert!(
                (st.total() - n0).abs() < 1e-6,
                "at t={t}: total={tot} drifted from N0={n0}",
                t = s.t,
                tot = st.total()
            );
        }
    }

    #[test]
    fn epidemic_grows_when_r0_above_one() {
        // VALIDATE: epidemic grows iff R0 > 1 (R0 = 0.5/0.2 = 2.5 here).
        let m = Sir::new(0.5, 0.2).unwrap();
        assert!(m.r0() > 1.0);
        let init = SirState::new(999.0, 1.0, 0.0);
        let traj = m.simulate(init, 60.0, 0.05).unwrap();
        let peak_i = traj.iter().map(|s| s.y[1]).fold(f64::MIN, f64::max);
        // Infectious count must rise well above its seed of 1.
        assert!(peak_i > 1.0, "epidemic did not take off: peak I = {peak_i}");
    }

    #[test]
    fn epidemic_dies_out_when_r0_below_one() {
        // The complement: R0 = 0.1/0.2 = 0.5 < 1 => I decays monotonically.
        let m = Sir::new(0.1, 0.2).unwrap();
        assert!(m.r0() < 1.0);
        let init = SirState::new(999.0, 1.0, 0.0);
        let traj = m.simulate(init, 60.0, 0.05).unwrap();
        // I should never exceed its initial value.
        let peak_i = traj.iter().map(|s| s.y[1]).fold(f64::MIN, f64::max);
        assert!(
            peak_i <= 1.0 + 1e-9,
            "I rose above seed despite R0<1: peak {peak_i}"
        );
        // And it should have decayed substantially by the end.
        assert!(traj.last().unwrap().y[1] < 0.5);
    }

    #[test]
    fn susceptibles_only_decrease_recovered_only_increase() {
        // Structural sanity: S is non-increasing, R is non-decreasing.
        let m = Sir::new(0.4, 0.15).unwrap();
        let init = SirState::new(995.0, 5.0, 0.0);
        let traj = m.simulate(init, 50.0, 0.05).unwrap();
        for w in traj.windows(2) {
            assert!(w[1].y[0] <= w[0].y[0] + 1e-9, "S increased");
            assert!(w[1].y[2] >= w[0].y[2] - 1e-9, "R decreased");
        }
    }

    #[test]
    fn final_size_satisfies_implicit_relation() {
        // Cross-check against the closed-form final-size relation:
        //   s_inf = s0 * exp(-R0 * (1 - s_inf))   (frequency-dependent,
        // negligible initial infectious fraction), where s = S/N.
        let m = Sir::new(0.6, 0.2).unwrap(); // R0 = 3
        let n = 1_000_000.0;
        let init = SirState::new(n - 1.0, 1.0, 0.0);
        let traj = m.simulate(init, 120.0, 0.02).unwrap();
        let s_inf = traj.last().unwrap().y[0] / n;

        // Solve the transcendental relation by fixed-point iteration.
        let r0 = m.r0();
        let s0 = (n - 1.0) / n;
        let mut x = 0.05;
        for _ in 0..1000 {
            x = s0 * (-r0 * (1.0 - x)).exp();
        }
        assert!(
            (s_inf - x).abs() < 2e-3,
            "final S fraction {s_inf} vs analytic final-size {x}"
        );
    }
}
