//! The Lotka-Volterra predator-prey system.
//!
//! ## Model
//!
//! The classic two-species Lotka-Volterra equations couple a prey
//! population `x` and a predator population `y`:
//!
//! ```text
//! dx/dt =  alpha x - beta  x y      (prey: growth - predation)
//! dy/dt = -gamma y + delta x y      (predator: death + growth from prey)
//! ```
//!
//! with all four parameters positive: `alpha` the prey birth rate,
//! `beta` the predation rate, `gamma` the predator death rate, and
//! `delta` the predator reproduction-per-prey-eaten rate.
//!
//! ## Behaviour
//!
//! The system has a non-trivial equilibrium at
//!
//! ```text
//! x* = gamma / delta,   y* = alpha / beta,
//! ```
//!
//! around which trajectories are *closed orbits* — the populations
//! oscillate periodically and out of phase, with the **prey peak
//! preceding the predator peak** (predators rise in response to
//! abundant prey, then crash once the prey are depleted). The motion
//! conserves the quantity
//!
//! ```text
//! H = delta x - gamma ln x + beta y - alpha ln y,
//! ```
//!
//! a constant of motion. The [`tests`](self) verify the equilibrium is
//! a fixed point, that `H` is (approximately) conserved by the RK4
//! integrator, that the populations oscillate, and that the prey peak
//! comes before the predator peak.

use crate::error::{PopError, Result};
use crate::rk4::{integrate, Sample};
use serde::{Deserialize, Serialize};

/// Parameters of the Lotka-Volterra predator-prey model.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LotkaVolterra {
    /// Prey intrinsic birth rate `alpha` (per unit time), `> 0`.
    pub alpha: f64,
    /// Predation rate `beta` (per predator per unit time), `> 0`.
    pub beta: f64,
    /// Predator death rate `gamma` (per unit time), `> 0`.
    pub gamma: f64,
    /// Predator reproduction rate per prey consumed `delta`, `> 0`.
    pub delta: f64,
}

impl LotkaVolterra {
    /// Construct a validated Lotka-Volterra model.
    ///
    /// # Errors
    ///
    /// [`PopError::Invalid`] if any of the four parameters is not
    /// strictly positive or not finite. `gamma` and `delta` must be
    /// positive so the equilibrium `(gamma/delta, alpha/beta)` is well
    /// defined; `alpha` and `beta` likewise.
    pub fn new(alpha: f64, beta: f64, gamma: f64, delta: f64) -> Result<Self> {
        for (name, v) in [
            ("alpha", alpha),
            ("beta", beta),
            ("gamma", gamma),
            ("delta", delta),
        ] {
            if !v.is_finite() || v <= 0.0 {
                return Err(PopError::invalid(
                    name,
                    "rate must be finite and strictly positive",
                ));
            }
        }
        Ok(LotkaVolterra {
            alpha,
            beta,
            gamma,
            delta,
        })
    }

    /// The non-trivial coexistence equilibrium `(x*, y*) =
    /// (gamma/delta, alpha/beta)`, returned as `[prey, predator]`.
    pub fn equilibrium(&self) -> [f64; 2] {
        [self.gamma / self.delta, self.alpha / self.beta]
    }

    /// The two time-derivatives `[dx/dt, dy/dt]` at state
    /// `[prey, predator]`.
    pub fn rate(&self, state: &[f64; 2]) -> [f64; 2] {
        let (x, y) = (state[0], state[1]);
        let dx = self.alpha * x - self.beta * x * y;
        let dy = -self.gamma * y + self.delta * x * y;
        [dx, dy]
    }

    /// The conserved quantity
    /// `H = delta x - gamma ln x + beta y - alpha ln y`.
    ///
    /// Constant along exact trajectories; used to check how well the
    /// integrator preserves the closed orbit. Defined only for strictly
    /// positive populations (the logs); returns `NaN` otherwise.
    pub fn conserved_quantity(&self, state: &[f64; 2]) -> f64 {
        let (x, y) = (state[0], state[1]);
        self.delta * x - self.gamma * x.ln() + self.beta * y - self.alpha * y.ln()
    }

    /// Integrate the predator-prey system from an initial
    /// `[prey, predator]` state to `t_end` with RK4 step `dt`.
    ///
    /// # Errors
    ///
    /// [`PopError::Invalid`] if either initial population is negative or
    /// non-finite, plus any error from the underlying [`integrate`].
    pub fn simulate(&self, initial: [f64; 2], t_end: f64, dt: f64) -> Result<Vec<Sample<2>>> {
        for (name, v) in [("prey0", initial[0]), ("predator0", initial[1])] {
            if !v.is_finite() || v < 0.0 {
                return Err(PopError::invalid(
                    name,
                    "initial population must be finite and non-negative",
                ));
            }
        }
        let model = *self;
        integrate(move |_t, y| model.rate(y), initial, 0.0, t_end, dt)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Time of the first interior *local* maximum of one component of a
    /// trajectory (a sample strictly greater than its left neighbour and
    /// `>=` its right neighbour). Returns `None` if the component is
    /// monotone over the window. Used to compare the phase of the prey
    /// and predator peaks within a single oscillation cycle — the
    /// meaningful "who peaks first" question for Lotka-Volterra, which a
    /// global maximum over a multi-cycle window would answer ambiguously.
    fn first_local_max_time(traj: &[Sample<2>], comp: usize) -> Option<f64> {
        traj.windows(3).find_map(|w| {
            let rising = w[1].y[comp] > w[0].y[comp];
            let cresting = w[1].y[comp] >= w[2].y[comp];
            if rising && cresting {
                Some(w[1].t)
            } else {
                None
            }
        })
    }

    #[test]
    fn rejects_nonpositive_parameters() {
        assert!(LotkaVolterra::new(0.0, 1.0, 1.0, 1.0).is_err());
        assert!(LotkaVolterra::new(1.0, -1.0, 1.0, 1.0).is_err());
        assert!(LotkaVolterra::new(1.0, 1.0, 0.0, 1.0).is_err());
        assert!(LotkaVolterra::new(1.0, 1.0, 1.0, f64::NAN).is_err());
    }

    #[test]
    fn equilibrium_is_a_fixed_point() {
        let m = LotkaVolterra::new(1.1, 0.4, 0.4, 0.1).unwrap();
        let eq = m.equilibrium();
        // x* = gamma/delta = 0.4/0.1 = 4, y* = alpha/beta = 1.1/0.4 = 2.75.
        assert!((eq[0] - 4.0).abs() < 1e-12, "x*={}", eq[0]);
        assert!((eq[1] - 2.75).abs() < 1e-12, "y*={}", eq[1]);
        // At the equilibrium both derivatives vanish.
        let d = m.rate(&eq);
        assert!(d[0].abs() < 1e-12 && d[1].abs() < 1e-12, "rate={d:?}");
    }

    #[test]
    fn equilibrium_start_stays_put() {
        // Starting exactly at equilibrium the populations must not move.
        let m = LotkaVolterra::new(1.0, 0.5, 0.75, 0.25).unwrap();
        let eq = m.equilibrium();
        let traj = m.simulate(eq, 30.0, 0.01).unwrap();
        for s in &traj {
            assert!((s.y[0] - eq[0]).abs() < 1e-6, "prey drifted: {:?}", s.y);
            assert!((s.y[1] - eq[1]).abs() < 1e-6, "predator drifted: {:?}", s.y);
        }
    }

    #[test]
    fn conserved_quantity_is_preserved() {
        // VALIDATE: closed orbit — H is (approximately) conserved by RK4.
        let m = LotkaVolterra::new(1.0, 0.1, 1.5, 0.075).unwrap();
        let init = [10.0, 5.0];
        let h0 = m.conserved_quantity(&init);
        let traj = m.simulate(init, 30.0, 0.005).unwrap();
        for s in &traj {
            let h = m.conserved_quantity(&s.y);
            assert!(
                (h - h0).abs() < 1e-2 * h0.abs().max(1.0),
                "H drifted at t={t}: {h} vs {h0}",
                t = s.t
            );
        }
    }

    #[test]
    fn populations_oscillate() {
        // VALIDATE: Lotka-Volterra oscillates — prey both rises above
        // and falls below its starting value over a long enough run.
        let m = LotkaVolterra::new(1.0, 0.1, 1.5, 0.075).unwrap();
        let init = [10.0, 5.0];
        let traj = m.simulate(init, 30.0, 0.005).unwrap();
        let prey_max = traj.iter().map(|s| s.y[0]).fold(f64::MIN, f64::max);
        let prey_min = traj.iter().map(|s| s.y[0]).fold(f64::MAX, f64::min);
        assert!(
            prey_max > init[0] + 1.0 && prey_min < init[0] - 1.0,
            "prey did not oscillate: min={prey_min} max={prey_max}"
        );
        // Populations stay positive throughout a closed orbit.
        assert!(prey_min > 0.0);
    }

    #[test]
    fn prey_peak_precedes_predator_peak() {
        // VALIDATE: prey peak precedes predator peak.
        //
        // The phase ordering is a property of one oscillation cycle:
        // prey bloom -> prey peak -> predators (fed) rise -> predator
        // peak -> prey crash. So we compare the *first local maximum* of
        // each species, not a global maximum over a multi-cycle window
        // (which could pick a later, taller prey peak). Starting both
        // populations below their equilibrium, the prey crest leads.
        let m = LotkaVolterra::new(1.0, 0.1, 1.5, 0.075).unwrap();
        // Equilibrium: x* = 1.5/0.075 = 20, y* = 1.0/0.1 = 10.
        let init = [10.0, 5.0];
        let traj = m.simulate(init, 12.0, 0.002).unwrap();
        let prey_peak_t = first_local_max_time(&traj, 0).expect("prey oscillates");
        let pred_peak_t = first_local_max_time(&traj, 1).expect("predator oscillates");
        assert!(
            prey_peak_t < pred_peak_t,
            "prey peak (t={prey_peak_t:.3}) should precede predator peak (t={pred_peak_t:.3})"
        );
        // The lag should be a meaningful fraction of the cycle (a quarter
        // period for small oscillations), not a numerical hair.
        assert!(
            pred_peak_t - prey_peak_t > 0.1,
            "predator lag too small: prey={prey_peak_t:.3} pred={pred_peak_t:.3}"
        );
    }
}
