//! ODE integrators — features 9, 10 and 11.
//!
//! Three integrators, each consuming an [`OdeSystem`]:
//!
//! - [`rk4_step`] / [`integrate_rk4`] — the classic explicit
//!   fourth-order Runge-Kutta scheme, fixed step (feature 9). Cheap
//!   and accurate for non-stiff problems.
//! - [`Rk45`] — the Dormand-Prince embedded RK4(5) pair with PI step
//!   control (feature 10). Adapts the step to a user error tolerance;
//!   the workhorse for general non-stiff kinetics.
//! - [`Bdf`] — a variable-order (1-2) backward-differentiation-formula
//!   integrator with a Newton inner solve (feature 11). Stable on
//!   stiff systems where the explicit methods would need an
//!   impractically small step.
//!
//! Each integrator returns a [`Trajectory`] — paired time and state
//! samples — or a [`SysbioError::NotConverged`] if it cannot advance.
//!
//! ## v1 caveats
//!
//! The [`Bdf`] integrator is order 1-2 only (BDF1 = implicit Euler,
//! BDF2 once enough history exists). Production stiff solvers (CVODE,
//! LSODA) go to order 5 with sophisticated order/step heuristics; this
//! v1 deliberately stops at order 2 — order 2 is A-stable and captures
//! the defining benefit (unconditional stability) without the
//! Nordsieck-history bookkeeping of the higher orders. The Newton
//! solve uses a finite-difference Jacobian refactored every step.

use crate::error::{Result, SysbioError};
use crate::ode::linalg::solve_linear;
use crate::ode::OdeSystem;

/// A solved trajectory: parallel time and state samples.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Trajectory {
    /// Sample times, strictly increasing.
    pub times: Vec<f64>,
    /// State vectors, one per time (`states[k].len() == dim`).
    pub states: Vec<Vec<f64>>,
}

impl Trajectory {
    /// Number of stored samples.
    pub fn len(&self) -> usize {
        self.times.len()
    }

    /// Whether the trajectory is empty.
    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }

    /// The final state, if any.
    pub fn final_state(&self) -> Option<&[f64]> {
        self.states.last().map(|v| v.as_slice())
    }

    /// The time series of species `i` across every sample.
    pub fn series(&self, i: usize) -> Vec<f64> {
        self.states.iter().map(|s| s[i]).collect()
    }
}

/// One explicit RK4 step of size `h` from `(t, y)`. Returns `y(t+h)`.
pub fn rk4_step(sys: &OdeSystem, t: f64, y: &[f64], h: f64) -> Vec<f64> {
    let n = y.len();
    let k1 = sys.rhs(t, y);
    let mut tmp = vec![0.0; n];
    for i in 0..n {
        tmp[i] = y[i] + 0.5 * h * k1[i];
    }
    let k2 = sys.rhs(t + 0.5 * h, &tmp);
    for i in 0..n {
        tmp[i] = y[i] + 0.5 * h * k2[i];
    }
    let k3 = sys.rhs(t + 0.5 * h, &tmp);
    for i in 0..n {
        tmp[i] = y[i] + h * k3[i];
    }
    let k4 = sys.rhs(t + h, &tmp);
    let mut out = vec![0.0; n];
    for i in 0..n {
        out[i] = y[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    out
}

/// Fixed-step RK4 integration from `t0` to `t_end` (feature 9).
///
/// Stores a sample every `n_out`-th step (and always the endpoints).
/// Errors on a non-positive step or a zero-length interval.
pub fn integrate_rk4(
    sys: &OdeSystem,
    y0: &[f64],
    t0: f64,
    t_end: f64,
    dt: f64,
    n_out: usize,
) -> Result<Trajectory> {
    if dt <= 0.0 {
        return Err(SysbioError::invalid("dt", "step must be positive"));
    }
    if t_end <= t0 {
        return Err(SysbioError::invalid("t_end", "t_end must exceed t0"));
    }
    let n_out = n_out.max(1);
    let steps = ((t_end - t0) / dt).ceil() as usize;
    let mut traj = Trajectory {
        times: vec![t0],
        states: vec![y0.to_vec()],
    };
    let mut t = t0;
    let mut y = y0.to_vec();
    for step in 0..steps {
        let h = dt.min(t_end - t);
        y = rk4_step(sys, t, &y, h);
        t += h;
        if (step + 1) % n_out == 0 || step + 1 == steps {
            traj.times.push(t);
            traj.states.push(y.clone());
        }
    }
    Ok(traj)
}

/// Dormand-Prince adaptive RK4(5) integrator (feature 10).
#[derive(Debug, Clone)]
pub struct Rk45 {
    /// Absolute error tolerance.
    pub atol: f64,
    /// Relative error tolerance.
    pub rtol: f64,
    /// Initial step-size guess.
    pub h0: f64,
    /// Smallest step the controller will accept before declaring
    /// failure.
    pub h_min: f64,
    /// Maximum number of accepted + rejected steps.
    pub max_steps: usize,
}

impl Default for Rk45 {
    fn default() -> Self {
        Rk45 {
            atol: 1e-8,
            rtol: 1e-6,
            h0: 1e-3,
            h_min: 1e-12,
            max_steps: 1_000_000,
        }
    }
}

impl Rk45 {
    /// Integrate `sys` from `t0` to `t_end`. Every accepted step
    /// appends a sample (the step size is chosen by the controller, so
    /// samples are non-uniform in time).
    pub fn integrate(
        &self,
        sys: &OdeSystem,
        y0: &[f64],
        t0: f64,
        t_end: f64,
    ) -> Result<Trajectory> {
        if t_end <= t0 {
            return Err(SysbioError::invalid("t_end", "t_end must exceed t0"));
        }
        if self.h0 <= 0.0 {
            return Err(SysbioError::invalid("h0", "initial step must be positive"));
        }
        let n = y0.len();
        // Dormand-Prince Butcher tableau.
        const C: [f64; 7] = [0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0];
        const A: [[f64; 6]; 7] = [
            [0.0; 6],
            [0.2, 0.0, 0.0, 0.0, 0.0, 0.0],
            [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
            [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
            [
                19372.0 / 6561.0,
                -25360.0 / 2187.0,
                64448.0 / 6561.0,
                -212.0 / 729.0,
                0.0,
                0.0,
            ],
            [
                9017.0 / 3168.0,
                -355.0 / 33.0,
                46732.0 / 5247.0,
                49.0 / 176.0,
                -5103.0 / 18656.0,
                0.0,
            ],
            [
                35.0 / 384.0,
                0.0,
                500.0 / 1113.0,
                125.0 / 192.0,
                -2187.0 / 6784.0,
                11.0 / 84.0,
            ],
        ];
        // 5th-order solution weights.
        const B5: [f64; 7] = [
            35.0 / 384.0,
            0.0,
            500.0 / 1113.0,
            125.0 / 192.0,
            -2187.0 / 6784.0,
            11.0 / 84.0,
            0.0,
        ];
        // 4th-order embedded weights (for error estimation).
        const B4: [f64; 7] = [
            5179.0 / 57600.0,
            0.0,
            7571.0 / 16695.0,
            393.0 / 640.0,
            -92097.0 / 339200.0,
            187.0 / 2100.0,
            1.0 / 40.0,
        ];

        let mut traj = Trajectory {
            times: vec![t0],
            states: vec![y0.to_vec()],
        };
        let mut t = t0;
        let mut y = y0.to_vec();
        let mut h = self.h0.min(t_end - t0);
        let mut count = 0usize;

        while t < t_end - 1e-15 {
            count += 1;
            if count > self.max_steps {
                return Err(SysbioError::not_converged(
                    "rk45",
                    "exceeded the maximum step count",
                ));
            }
            if h < self.h_min {
                return Err(SysbioError::not_converged(
                    "rk45",
                    "step size fell below h_min",
                ));
            }
            h = h.min(t_end - t);

            // Seven stage derivatives.
            let mut k: Vec<Vec<f64>> = Vec::with_capacity(7);
            for s in 0..7 {
                let mut ys = y.clone();
                for (j, kj) in k.iter().enumerate() {
                    let a = A[s][j];
                    if a != 0.0 {
                        for i in 0..n {
                            ys[i] += h * a * kj[i];
                        }
                    }
                }
                k.push(sys.rhs(t + C[s] * h, &ys));
            }

            // 5th- and 4th-order increments.
            let mut y5 = y.clone();
            let mut err = 0.0;
            for i in 0..n {
                let mut inc5 = 0.0;
                let mut inc4 = 0.0;
                for (s, ks) in k.iter().enumerate() {
                    inc5 += B5[s] * ks[i];
                    inc4 += B4[s] * ks[i];
                }
                y5[i] += h * inc5;
                let e = h * (inc5 - inc4);
                let scale = self.atol + self.rtol * y[i].abs().max(y5[i].abs());
                err += (e / scale).powi(2);
            }
            err = (err / n as f64).sqrt();

            if err <= 1.0 {
                // Accept.
                t += h;
                y = y5;
                traj.times.push(t);
                traj.states.push(y.clone());
            }
            // PI-ish step update with the standard safety factor.
            let fac = if err == 0.0 {
                5.0
            } else {
                (0.9 * err.powf(-0.2)).clamp(0.2, 5.0)
            };
            h *= fac;
        }
        Ok(traj)
    }
}

/// Implicit BDF integrator for stiff systems (feature 11).
#[derive(Debug, Clone)]
pub struct Bdf {
    /// Fixed step size. (A v1 simplification — see the module docs.)
    pub h: f64,
    /// Newton convergence tolerance.
    pub newton_tol: f64,
    /// Maximum Newton iterations per step.
    pub newton_max: usize,
    /// Store a sample every `n_out`-th step.
    pub n_out: usize,
}

impl Default for Bdf {
    fn default() -> Self {
        Bdf {
            h: 1e-2,
            newton_tol: 1e-9,
            newton_max: 50,
            n_out: 1,
        }
    }
}

impl Bdf {
    /// Integrate `sys` from `t0` to `t_end` with a fixed step.
    ///
    /// The first step uses BDF1 (implicit Euler); subsequent steps use
    /// BDF2 once two history points exist. Each step solves the
    /// implicit equation with a damped Newton iteration on a
    /// finite-difference Jacobian.
    pub fn integrate(
        &self,
        sys: &OdeSystem,
        y0: &[f64],
        t0: f64,
        t_end: f64,
    ) -> Result<Trajectory> {
        if self.h <= 0.0 {
            return Err(SysbioError::invalid("h", "step must be positive"));
        }
        if t_end <= t0 {
            return Err(SysbioError::invalid("t_end", "t_end must exceed t0"));
        }
        let n = y0.len();
        let n_out = self.n_out.max(1);
        let steps = ((t_end - t0) / self.h).ceil() as usize;

        let mut traj = Trajectory {
            times: vec![t0],
            states: vec![y0.to_vec()],
        };
        let mut t = t0;
        let mut y_prev2: Option<Vec<f64>> = None;
        let mut y_prev = y0.to_vec();

        for step in 0..steps {
            let h = self.h.min(t_end - t);
            // BDF coefficients: solve  a0·y - rhs_factor·h·f(y) = b·hist.
            // BDF1:  y - y_prev = h f(y)
            // BDF2:  (3/2)y - 2 y_prev + (1/2) y_prev2 = h f(y)
            let (a0, rhs_const): (f64, Vec<f64>) = match &y_prev2 {
                Some(p2) => {
                    let mut c = vec![0.0; n];
                    for i in 0..n {
                        c[i] = 2.0 * y_prev[i] - 0.5 * p2[i];
                    }
                    (1.5, c)
                }
                None => (1.0, y_prev.clone()),
            };
            let t_new = t + h;

            // Newton solve for y satisfying  G(y) = a0 y - h f(y) - rhs_const = 0.
            let mut y = y_prev.clone(); // initial guess
            let mut converged = false;
            for _ in 0..self.newton_max {
                let f = sys.rhs(t_new, &y);
                let mut g = vec![0.0; n];
                let mut gnorm = 0.0;
                for i in 0..n {
                    g[i] = a0 * y[i] - h * f[i] - rhs_const[i];
                    gnorm += g[i] * g[i];
                }
                if gnorm.sqrt() <= self.newton_tol {
                    converged = true;
                    break;
                }
                // Jacobian of G: a0 I - h J(f).
                let jf = sys.jacobian(t_new, &y);
                let mut jg = vec![vec![0.0; n]; n];
                for i in 0..n {
                    for j in 0..n {
                        jg[i][j] = -h * jf[i][j];
                    }
                    jg[i][i] += a0;
                }
                // Solve  jg · delta = -g.
                let neg_g: Vec<f64> = g.iter().map(|x| -x).collect();
                let delta = solve_linear(&jg, &neg_g).ok_or_else(|| {
                    SysbioError::not_converged("bdf", "singular Newton Jacobian")
                })?;
                // Damped update — half-step if the full step grows G.
                let mut lambda = 1.0;
                for _ in 0..8 {
                    let trial: Vec<f64> =
                        y.iter().zip(&delta).map(|(a, d)| a + lambda * d).collect();
                    let ft = sys.rhs(t_new, &trial);
                    let mut tnorm = 0.0;
                    for i in 0..n {
                        let gi = a0 * trial[i] - h * ft[i] - rhs_const[i];
                        tnorm += gi * gi;
                    }
                    if tnorm.sqrt() < gnorm.sqrt() || lambda < 1e-3 {
                        y = trial;
                        break;
                    }
                    lambda *= 0.5;
                }
            }
            if !converged {
                return Err(SysbioError::not_converged(
                    "bdf",
                    "Newton iteration did not reach the tolerance",
                ));
            }

            // Project the step onto the non-negative orthant. The state
            // of a reaction-network ODE is a vector of species
            // *amounts*, which are physically non-negative — and the
            // rate laws enforce that by clamping negative amounts to
            // zero ([`RateLaw::rate`]). On a stiff transient the BDF2
            // history term can extrapolate a species a hair below zero;
            // left unclamped, the rate (and its Jacobian) flatten in
            // that region, so the implicit equation `a0·y − h·f − c = 0`
            // degenerates to `a0·y = c` and the Newton iteration locks
            // onto a large spurious negative root. Clamping the
            // converged step to `≥ 0` keeps the integrator inside the
            // model's valid domain; the clamped excursion is within the
            // method's truncation error, so accuracy is preserved.
            for v in y.iter_mut() {
                if *v < 0.0 {
                    *v = 0.0;
                }
            }

            y_prev2 = Some(y_prev);
            y_prev = y;
            t = t_new;
            if (step + 1) % n_out == 0 || step + 1 == steps {
                traj.times.push(t);
                traj.states.push(y_prev.clone());
            }
        }
        Ok(traj)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Model, RateLaw, Reaction, Species};

    /// Pure decay A -> 0. Analytic: A(t) = A0 exp(-k t).
    fn decay(k: f64, a0: f64) -> OdeSystem {
        let mut m = Model::new("decay");
        let a = m.add_species(Species::new("A", a0));
        m.add_reaction(Reaction {
            id: "d".into(),
            reactants: vec![(a, 1.0)],
            products: vec![],
            rate_law: RateLaw::MassAction {
                k,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        OdeSystem::from_model(&m)
    }

    #[test]
    fn rk4_matches_analytic_decay() {
        let sys = decay(1.0, 10.0);
        let traj = integrate_rk4(&sys, &[10.0], 0.0, 2.0, 0.01, 10).unwrap();
        let final_a = traj.final_state().unwrap()[0];
        let expect = 10.0 * (-2.0_f64).exp();
        assert!((final_a - expect).abs() < 1e-4, "got {final_a}, want {expect}");
    }

    #[test]
    fn rk45_matches_analytic_decay() {
        let sys = decay(0.5, 4.0);
        let traj = Rk45::default().integrate(&sys, &[4.0], 0.0, 5.0).unwrap();
        let final_a = traj.final_state().unwrap()[0];
        let expect = 4.0 * (-2.5_f64).exp();
        assert!((final_a - expect).abs() < 1e-5, "got {final_a}, want {expect}");
    }

    #[test]
    fn bdf_matches_analytic_decay() {
        let sys = decay(1.0, 10.0);
        let traj = Bdf {
            h: 1e-3,
            ..Bdf::default()
        }
        .integrate(&sys, &[10.0], 0.0, 2.0)
        .unwrap();
        let final_a = traj.final_state().unwrap()[0];
        let expect = 10.0 * (-2.0_f64).exp();
        assert!((final_a - expect).abs() < 1e-2, "got {final_a}, want {expect}");
    }

    #[test]
    fn bdf_stays_stable_on_stiff_system() {
        // dy/dt = -1000 (y - cos t)-ish stiffness via a fast decay.
        // Explicit Euler at h=0.1 would blow up; BDF must not.
        let sys = decay(1000.0, 1.0);
        let traj = Bdf {
            h: 0.1,
            ..Bdf::default()
        }
        .integrate(&sys, &[1.0], 0.0, 5.0)
        .unwrap();
        let final_a = traj.final_state().unwrap()[0];
        // Should have decayed to ~0, and crucially be finite & bounded.
        assert!(final_a.is_finite());
        assert!(final_a.abs() < 1e-3, "stiff BDF unstable: {final_a}");
    }

    #[test]
    fn rejects_bad_interval() {
        let sys = decay(1.0, 1.0);
        assert!(integrate_rk4(&sys, &[1.0], 0.0, -1.0, 0.1, 1).is_err());
        assert!(Rk45::default().integrate(&sys, &[1.0], 5.0, 1.0).is_err());
    }
}
