//! Single-species logistic (Verhulst) growth.
//!
//! ## Model
//!
//! The logistic equation models a single population `N(t)` with an
//! intrinsic per-capita growth rate `r` that is throttled as the
//! population approaches a carrying capacity `K`:
//!
//! ```text
//! dN/dt = r N (1 - N/K).
//! ```
//!
//! It has the closed-form solution
//!
//! ```text
//! N(t) = K N0 e^{rt} / (K + N0 (e^{rt} - 1)),
//! ```
//!
//! which the [`tests`](self) use as ground truth for the RK4
//! trajectory. For `r > 0` and `0 < N0 < K` the population rises
//! sigmoidally and converges monotonically to `K`; for `N0 > K` it
//! falls to `K`; `N = 0` and `N = K` are the two equilibria (`K` is
//! the stable one).

use crate::error::{PopError, Result};
use crate::rk4::{integrate, Sample};
use serde::{Deserialize, Serialize};

/// Parameters of the logistic growth model.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Logistic {
    /// Intrinsic per-capita growth rate `r` (per unit time). Positive
    /// for a growing population.
    pub r: f64,
    /// Carrying capacity `K` — the population the environment supports.
    /// Must be strictly positive (it appears in the denominator).
    pub k: f64,
}

impl Logistic {
    /// Construct a validated logistic model.
    ///
    /// # Errors
    ///
    /// [`PopError::Invalid`] if `k <= 0` or if either parameter is not
    /// finite. A non-positive or zero `K` is rejected because it
    /// appears in the `N/K` denominator.
    pub fn new(r: f64, k: f64) -> Result<Self> {
        if !r.is_finite() {
            return Err(PopError::invalid("r", "growth rate must be finite"));
        }
        if !k.is_finite() || k <= 0.0 {
            return Err(PopError::invalid(
                "k",
                "carrying capacity must be finite and strictly positive",
            ));
        }
        Ok(Logistic { r, k })
    }

    /// The instantaneous rate `dN/dt = r N (1 - N/K)`.
    pub fn rate(&self, n: f64) -> f64 {
        self.r * n * (1.0 - n / self.k)
    }

    /// Closed-form analytic solution `N(t)` given the initial
    /// population `n0` at `t = 0`.
    ///
    /// Returned for validation and for callers who want the exact curve
    /// without integrating. Valid for any `n0 >= 0`.
    pub fn analytic(&self, n0: f64, t: f64) -> f64 {
        let e = (self.r * t).exp();
        self.k * n0 * e / (self.k + n0 * (e - 1.0))
    }

    /// Integrate the population from `n0` at `t = 0` to `t_end` with RK4
    /// step `dt`, returning the `[N]` trajectory.
    ///
    /// # Errors
    ///
    /// [`PopError::Invalid`] if `n0 < 0`, plus any error from the
    /// underlying [`integrate`] (non-positive `dt`, bad window, or a
    /// step-count overflow).
    pub fn simulate(&self, n0: f64, t_end: f64, dt: f64) -> Result<Vec<Sample<1>>> {
        if !n0.is_finite() || n0 < 0.0 {
            return Err(PopError::invalid(
                "n0",
                "initial population must be finite and non-negative",
            ));
        }
        let model = *self;
        integrate(move |_t, y| [model.rate(y[0])], [n0], 0.0, t_end, dt)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_nonpositive_capacity() {
        assert!(Logistic::new(0.5, 0.0).is_err());
        assert!(Logistic::new(0.5, -10.0).is_err());
        assert!(Logistic::new(f64::NAN, 10.0).is_err());
    }

    #[test]
    fn equilibria_have_zero_rate() {
        let m = Logistic::new(0.8, 100.0).unwrap();
        // Both N = 0 and N = K are fixed points.
        assert!((m.rate(0.0)).abs() < 1e-12);
        assert!((m.rate(100.0)).abs() < 1e-12);
    }

    #[test]
    fn analytic_matches_rk4() {
        // Ground-truth check: RK4 trajectory vs the closed-form curve.
        let m = Logistic::new(0.7, 1000.0).unwrap();
        let n0 = 10.0;
        let traj = m.simulate(n0, 20.0, 0.01).unwrap();
        for s in &traj {
            let exact = m.analytic(n0, s.t);
            assert!(
                (s.y[0] - exact).abs() < 1e-4,
                "at t={t}: rk4={got} analytic={exact}",
                t = s.t,
                got = s.y[0]
            );
        }
    }

    #[test]
    fn converges_to_carrying_capacity_from_below() {
        // VALIDATE: logistic converges to K.
        let m = Logistic::new(1.0, 500.0).unwrap();
        let traj = m.simulate(1.0, 40.0, 0.01).unwrap();
        let final_n = traj.last().unwrap().y[0];
        assert!(
            (final_n - 500.0).abs() < 1e-2,
            "final N {final_n} did not converge to K=500"
        );
    }

    #[test]
    fn converges_to_carrying_capacity_from_above() {
        // Starting above K, the population must fall to K.
        let m = Logistic::new(1.0, 500.0).unwrap();
        let traj = m.simulate(900.0, 40.0, 0.01).unwrap();
        let final_n = traj.last().unwrap().y[0];
        assert!(
            (final_n - 500.0).abs() < 1e-2,
            "final N {final_n} did not settle to K=500 from above"
        );
        // And it must be monotonically decreasing toward K.
        assert!(traj[1].y[0] < traj[0].y[0]);
    }

    #[test]
    fn growth_is_monotonic_and_sigmoid_below_k() {
        // From below K with r>0 the trajectory never decreases.
        let m = Logistic::new(0.9, 200.0).unwrap();
        let traj = m.simulate(5.0, 30.0, 0.05).unwrap();
        for w in traj.windows(2) {
            assert!(
                w[1].y[0] >= w[0].y[0] - 1e-9,
                "population dipped: {} -> {}",
                w[0].y[0],
                w[1].y[0]
            );
            // Stays bounded by the capacity.
            assert!(w[1].y[0] <= 200.0 + 1e-6);
        }
    }
}
