//! Steady-state solver — feature 12.
//!
//! A steady state of `dy/dt = f(y)` is a root of `f`. [`steady_state`]
//! finds one by a **damped Newton iteration**: at each step it solves
//! the linear system `J·Δ = −f` for the Newton direction, then takes
//! the largest fraction `λ ∈ (0, 1]` of that step which actually
//! reduces `‖f‖` (a backtracking line search). Damping is what makes
//! the method robust far from the root, where an undamped Newton step
//! routinely overshoots into negative concentrations.
//!
//! When the Jacobian is rank-deficient — common for a network with a
//! conserved moiety, whose steady state is not isolated — the solver
//! falls back to the SVD pseudo-inverse direction
//! ([`solve_least_squares`]),
//! which returns the minimum-norm step and lands on *a* point of the
//! steady-state manifold.

use crate::error::{Result, SysbioError};
use crate::ode::linalg::{solve_least_squares, solve_linear};
use crate::ode::OdeSystem;

/// Outcome of a steady-state solve.
#[derive(Debug, Clone, PartialEq)]
pub struct SteadyState {
    /// The state vector at the located steady state.
    pub state: Vec<f64>,
    /// `‖f(state)‖₂` — the residual; below the tolerance on success.
    pub residual: f64,
    /// Number of Newton iterations taken.
    pub iterations: usize,
}

/// Newton-solve for a steady state starting from `y0` (feature 12).
///
/// `tol` is the residual `‖f‖₂` target; `max_iter` caps the iteration
/// count. Negative components are clamped to zero after each step
/// (amounts cannot be negative). Returns
/// [`SysbioError::NotConverged`] if the residual target is not met.
pub fn steady_state(
    sys: &OdeSystem,
    y0: &[f64],
    tol: f64,
    max_iter: usize,
) -> Result<SteadyState> {
    if tol <= 0.0 {
        return Err(SysbioError::invalid("tol", "tolerance must be positive"));
    }
    let n = y0.len();
    let mut y = y0.to_vec();

    let norm = |v: &[f64]| -> f64 { v.iter().map(|x| x * x).sum::<f64>().sqrt() };

    for iter in 0..max_iter {
        let f = sys.rhs(0.0, &y);
        let r = norm(&f);
        if r <= tol {
            return Ok(SteadyState {
                state: y,
                residual: r,
                iterations: iter,
            });
        }
        let jac = sys.jacobian(0.0, &y);
        let neg_f: Vec<f64> = f.iter().map(|x| -x).collect();
        // Prefer the exact solve; fall back to least-squares for a
        // rank-deficient Jacobian (conserved-moiety networks).
        let delta = solve_linear(&jac, &neg_f)
            .or_else(|| solve_least_squares(&jac, &neg_f))
            .ok_or_else(|| {
                SysbioError::not_converged("newton", "could not solve the Newton system")
            })?;

        // Backtracking line search on lambda.
        let mut lambda = 1.0;
        let mut accepted = false;
        for _ in 0..20 {
            let mut trial: Vec<f64> =
                y.iter().zip(&delta).map(|(a, d)| a + lambda * d).collect();
            for v in trial.iter_mut() {
                if *v < 0.0 {
                    *v = 0.0;
                }
            }
            let rt = norm(&sys.rhs(0.0, &trial));
            if rt < r {
                y = trial;
                accepted = true;
                break;
            }
            lambda *= 0.5;
        }
        if !accepted {
            // Tiny step still did not reduce the residual: take the
            // damped step anyway so a flat region does not deadlock.
            for i in 0..n {
                y[i] = (y[i] + 1e-3 * delta[i]).max(0.0);
            }
        }
    }

    let final_r = norm(&sys.rhs(0.0, &y));
    if final_r <= tol {
        Ok(SteadyState {
            state: y,
            residual: final_r,
            iterations: max_iter,
        })
    } else {
        Err(SysbioError::not_converged(
            "newton",
            format!("residual {final_r:.3e} above tolerance after {max_iter} iterations"),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Model, RateLaw, Reaction, Species};

    /// A constitutive source then first-order decay:
    /// 0 -> A (rate s), A -> 0 (rate k A). Steady state A* = s/k.
    fn source_decay(s: f64, k: f64) -> OdeSystem {
        let mut m = Model::new("srcdecay");
        let a = m.add_species(Species::new("A", 0.0));
        m.add_reaction(Reaction {
            id: "src".into(),
            reactants: vec![],
            products: vec![(a, 1.0)],
            rate_law: RateLaw::Constant { rate: s },
            reversible: false,
        });
        m.add_reaction(Reaction {
            id: "dec".into(),
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
    fn finds_source_decay_fixed_point() {
        let sys = source_decay(6.0, 2.0); // A* = 3
        let ss = steady_state(&sys, &[0.0], 1e-9, 100).unwrap();
        assert!((ss.state[0] - 3.0).abs() < 1e-6, "got {}", ss.state[0]);
        assert!(ss.residual < 1e-9);
    }

    #[test]
    fn converges_from_a_far_start() {
        let sys = source_decay(1.0, 0.25); // A* = 4
        let ss = steady_state(&sys, &[1000.0], 1e-8, 200).unwrap();
        assert!((ss.state[0] - 4.0).abs() < 1e-5);
    }

    #[test]
    fn rejects_nonpositive_tolerance() {
        let sys = source_decay(1.0, 1.0);
        assert!(steady_state(&sys, &[0.0], 0.0, 10).is_err());
    }
}
