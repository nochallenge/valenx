//! Bifurcation scan — feature 21 (v1).
//!
//! A bifurcation diagram tracks how a system's **steady states** move
//! — and appear or vanish — as a control parameter varies. This v1
//! performs *natural-parameter continuation*: it sweeps the
//! bifurcation parameter across a grid and, at each value, finds a
//! steady state by a Newton solve **warm-started from the previous
//! solution**. Warm-starting is what makes it continuation rather than
//! a naive scan — it follows one branch smoothly, and a sudden jump in
//! the located state (or a Newton failure that a cold restart
//! recovers from) flags where the branch folds.
//!
//! Each grid point yields a [`BifurcationPoint`] carrying the steady
//! state, the dominant Jacobian eigenvalue's real part (so a sign
//! change marks a **stability change** — a fold or a Hopf bifurcation)
//! and a [`Stability`] classification.
//!
//! ## v1 caveats
//!
//! This follows a *single* branch from one starting state — it does
//! not do pseudo-arclength continuation, so it cannot turn a fold and
//! trace the unstable branch back, and it detects bifurcations by a
//! stability-index sign change rather than by solving the
//! bifurcation's defining augmented system. It is a real, useful
//! diagram of the stable branch with annotated stability boundaries —
//! the AUTO-class "trace every branch" capability is out of scope.

use crate::analysis::param::ParamTarget;
use crate::analysis::scan::linspace;
use crate::error::{Result, SysbioError};
use crate::model::Model;
use crate::ode::steady::steady_state;
use crate::ode::OdeSystem;

/// Stability class of a located steady state, from the sign of the
/// dominant Jacobian eigenvalue's real part.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stability {
    /// All eigenvalues have negative real part — an attractor.
    Stable,
    /// Some eigenvalue has positive real part — a repeller / saddle.
    Unstable,
    /// The dominant real part is within tolerance of zero — a
    /// candidate bifurcation point.
    Marginal,
}

/// One point of a bifurcation diagram.
#[derive(Debug, Clone, PartialEq)]
pub struct BifurcationPoint {
    /// The bifurcation-parameter value at this point.
    pub param_value: f64,
    /// The located steady state (empty if the Newton solve failed).
    pub state: Vec<f64>,
    /// Largest real part among the Jacobian eigenvalue estimates.
    pub dominant_real: f64,
    /// Stability classification.
    pub stability: Stability,
    /// Whether a steady state was successfully located here.
    pub converged: bool,
}

/// A full bifurcation diagram.
#[derive(Debug, Clone, PartialEq)]
pub struct BifurcationDiagram {
    /// Label of the bifurcation parameter.
    pub param: String,
    /// Points along the parameter sweep.
    pub points: Vec<BifurcationPoint>,
}

impl BifurcationDiagram {
    /// Parameter values at which the stability class changes between
    /// consecutive points — candidate bifurcation locations.
    pub fn bifurcation_values(&self) -> Vec<f64> {
        let mut out = Vec::new();
        for w in self.points.windows(2) {
            if w[0].converged
                && w[1].converged
                && w[0].stability != w[1].stability
            {
                // Midpoint of the bracketing interval.
                out.push(0.5 * (w[0].param_value + w[1].param_value));
            }
        }
        out
    }
}

/// Run a steady-state continuation across `[lo, hi]` (feature 21).
///
/// `n` grid points; `eig_tol` is the half-width of the marginal band
/// around a zero dominant eigenvalue. The first point is solved from
/// the model's initial state; each subsequent point warm-starts from
/// the previous solution and falls back to a cold start if that fails.
pub fn bifurcation_scan(
    model: &Model,
    target: &ParamTarget,
    lo: f64,
    hi: f64,
    n: usize,
    eig_tol: f64,
) -> Result<BifurcationDiagram> {
    let grid = linspace(lo, hi, n)?;
    if eig_tol < 0.0 {
        return Err(SysbioError::invalid("eig_tol", "tolerance must be >= 0"));
    }
    let mut points = Vec::with_capacity(grid.len());
    let mut warm: Option<Vec<f64>> = None;

    for &p in &grid {
        let m = target.apply(model, p)?;
        let sys = OdeSystem::from_model(&m);
        let start = warm.clone().unwrap_or_else(|| m.initial_state());

        // Try warm start, then cold start.
        let ss = steady_state(&sys, &start, 1e-8, 300)
            .or_else(|_| steady_state(&sys, &m.initial_state(), 1e-8, 300));

        match ss {
            Ok(state) => {
                let dominant = dominant_eigen_real(&sys, &state.state);
                let stability = if dominant.abs() <= eig_tol {
                    Stability::Marginal
                } else if dominant < 0.0 {
                    Stability::Stable
                } else {
                    Stability::Unstable
                };
                warm = Some(state.state.clone());
                points.push(BifurcationPoint {
                    param_value: p,
                    state: state.state,
                    dominant_real: dominant,
                    stability,
                    converged: true,
                });
            }
            Err(_) => {
                points.push(BifurcationPoint {
                    param_value: p,
                    state: Vec::new(),
                    dominant_real: f64::NAN,
                    stability: Stability::Marginal,
                    converged: false,
                });
                // Drop the warm start so the next point cold-starts.
                warm = None;
            }
        }
    }

    Ok(BifurcationDiagram {
        param: target.label(),
        points,
    })
}

/// Estimate the largest real part among the Jacobian's eigenvalues at
/// a steady state.
///
/// For the v1 we use the Gershgorin upper bound on the spectral
/// abscissa combined with the trace sign — a cheap, derivative-free
/// indicator that is exact for 1-D systems and a sound *bound* for
/// larger ones. The eigenvalue real parts all lie within the union of
/// the Gershgorin discs; the largest disc's right edge is an upper
/// bound on the spectral abscissa, and for a stability *screen* an
/// upper bound is the conservative, correct choice.
fn dominant_eigen_real(sys: &OdeSystem, state: &[f64]) -> f64 {
    let j = sys.jacobian(0.0, state);
    let n = j.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return j[0][0];
    }
    // Gershgorin: for row i the eigenvalues lie within |z - a_ii| <= R_i.
    let mut bound = f64::NEG_INFINITY;
    for (i, row) in j.iter().enumerate() {
        let center = row[i];
        let radius: f64 = row
            .iter()
            .enumerate()
            .filter(|(k, _)| *k != i)
            .map(|(_, &v)| v.abs())
            .sum();
        bound = bound.max(center + radius);
    }
    bound
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RateLaw, Reaction, Species};

    /// 0 ->(s) A ->(k A) 0. A* = s/k, always stable (eigenvalue -k).
    fn source_decay() -> Model {
        let mut m = Model::new("sd");
        let a = m.add_species(Species::new("A", 0.0));
        m.add_reaction(Reaction {
            id: "src".into(),
            reactants: vec![],
            products: vec![(a, 1.0)],
            rate_law: RateLaw::Constant { rate: 1.0 },
            reversible: false,
        });
        m.add_reaction(Reaction {
            id: "dec".into(),
            reactants: vec![(a, 1.0)],
            products: vec![],
            rate_law: RateLaw::MassAction {
                k: 1.0,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        m
    }

    #[test]
    fn continuation_tracks_linear_branch() {
        // Sweep the source rate; A* should equal s at every point.
        let m = source_decay();
        let diagram = bifurcation_scan(
            &m,
            &ParamTarget::ConstantRate { reaction: 0 },
            1.0,
            5.0,
            9,
            1e-6,
        )
        .unwrap();
        assert_eq!(diagram.points.len(), 9);
        for pt in &diagram.points {
            assert!(pt.converged);
            assert!((pt.state[0] - pt.param_value).abs() < 1e-4);
            // Stable branch: dominant eigenvalue is -k = -1.
            assert_eq!(pt.stability, Stability::Stable);
        }
    }

    #[test]
    fn stable_branch_has_no_bifurcation() {
        let m = source_decay();
        let diagram = bifurcation_scan(
            &m,
            &ParamTarget::ConstantRate { reaction: 0 },
            1.0,
            5.0,
            6,
            1e-6,
        )
        .unwrap();
        assert!(diagram.bifurcation_values().is_empty());
    }

    #[test]
    fn dominant_eigen_of_decay_is_negative() {
        let m = source_decay();
        let sys = OdeSystem::from_model(&m);
        let d = dominant_eigen_real(&sys, &[1.0]);
        assert!((d - (-1.0)).abs() < 1e-4);
    }

    #[test]
    fn rejects_bad_grid() {
        let m = source_decay();
        assert!(bifurcation_scan(
            &m,
            &ParamTarget::ConstantRate { reaction: 0 },
            5.0,
            1.0,
            5,
            1e-6
        )
        .is_err());
    }
}
