//! Sensitivity analysis — feature 19.
//!
//! Sensitivity analysis asks *how much a model output moves when a
//! parameter moves*. Two complementary methods are provided:
//!
//! - [`local_sensitivity`] — **local**, finite-difference. For each
//!   targeted parameter it perturbs the value by a small relative
//!   step and measures the resulting change in a scalar readout,
//!   reporting both the raw derivative `∂y/∂p` and the dimensionless
//!   **scaled sensitivity** `(p/y)·(∂y/∂p)` (a control coefficient —
//!   the standard metabolic-control-analysis quantity). Local
//!   sensitivities are exact derivatives but valid only near the
//!   nominal point.
//!
//! - [`global_sensitivity`] — **global**, a Morris-style elementary-
//!   effects screening. It samples many random points across each
//!   parameter's full range, computes one one-at-a-time elementary
//!   effect per point, and reports the mean absolute effect `μ*`
//!   (overall influence) and the standard deviation `σ` (degree of
//!   non-linearity / interaction). Global sensitivity is robust to a
//!   non-linear, multi-modal response surface that would mislead the
//!   local method.
//!
//! ## v1 caveats
//!
//! The global method is the Morris *screening* design, not a full
//! variance-based Sobol decomposition — Morris is the standard,
//! cheap "which parameters matter" screen and is the right v1; Sobol
//! indices (with their much larger sample budget) are out of scope.

use crate::analysis::param::ParamTarget;
use crate::error::{Result, SysbioError};
use crate::model::Model;
use crate::stochastic::rng::Rng;

/// One row of a local sensitivity table.
#[derive(Debug, Clone, PartialEq)]
pub struct LocalSensitivity {
    /// Label of the parameter.
    pub param: String,
    /// Nominal parameter value.
    pub nominal: f64,
    /// Raw finite-difference derivative `∂readout/∂param`.
    pub derivative: f64,
    /// Dimensionless scaled sensitivity `(p/y)·(∂y/∂p)`.
    pub scaled: f64,
}

/// Local finite-difference sensitivity of `readout` to each `target`
/// (feature 19, local).
///
/// `rel_step` is the relative perturbation (e.g. `1e-3`). A central
/// difference is used. The nominal readout `y0` is evaluated once and
/// reused for every scaled sensitivity.
pub fn local_sensitivity<F>(
    model: &Model,
    targets: &[ParamTarget],
    rel_step: f64,
    mut readout: F,
) -> Result<Vec<LocalSensitivity>>
where
    F: FnMut(&Model) -> Result<f64>,
{
    if rel_step <= 0.0 {
        return Err(SysbioError::invalid("rel_step", "step must be positive"));
    }
    let y0 = readout(model)?;
    let mut out = Vec::with_capacity(targets.len());
    for t in targets {
        let p = t.read(model)?;
        // Absolute step scaled to the parameter magnitude.
        let h = rel_step * p.abs().max(1e-12);
        let m_plus = t.apply(model, p + h)?;
        let m_minus = t.apply(model, p - h)?;
        let y_plus = readout(&m_plus)?;
        let y_minus = readout(&m_minus)?;
        let derivative = (y_plus - y_minus) / (2.0 * h);
        let scaled = if y0.abs() > 1e-300 {
            (p / y0) * derivative
        } else {
            0.0
        };
        out.push(LocalSensitivity {
            param: t.label(),
            nominal: p,
            derivative,
            scaled,
        });
    }
    Ok(out)
}

/// One row of a Morris global-sensitivity table.
#[derive(Debug, Clone, PartialEq)]
pub struct GlobalSensitivity {
    /// Label of the parameter.
    pub param: String,
    /// `μ*` — mean of the absolute elementary effects (overall
    /// influence).
    pub mu_star: f64,
    /// `μ` — mean of the signed elementary effects (effect direction).
    pub mu: f64,
    /// `σ` — standard deviation of the elementary effects (a proxy for
    /// non-linearity and parameter interaction).
    pub sigma: f64,
}

/// Morris elementary-effects global sensitivity (feature 19, global).
///
/// For each of `n_trajectories` random base points, every parameter
/// is perturbed once by a step of `delta` of its range; the resulting
/// change in `readout`, divided by the step, is one *elementary
/// effect*. The per-parameter statistics summarise the distribution
/// of those effects over the whole input space.
///
/// `ranges` gives the `(lo, hi)` span of each `target`.
pub fn global_sensitivity<F>(
    model: &Model,
    targets: &[ParamTarget],
    ranges: &[(f64, f64)],
    n_trajectories: usize,
    delta: f64,
    seed: u64,
    mut readout: F,
) -> Result<Vec<GlobalSensitivity>>
where
    F: FnMut(&Model) -> Result<f64>,
{
    if targets.len() != ranges.len() {
        return Err(SysbioError::invalid(
            "ranges",
            "one range required per target",
        ));
    }
    if targets.is_empty() {
        return Err(SysbioError::invalid("targets", "need at least one target"));
    }
    if n_trajectories == 0 {
        return Err(SysbioError::invalid(
            "n_trajectories",
            "need at least one trajectory",
        ));
    }
    if !(0.0..1.0).contains(&delta) || delta <= 0.0 {
        return Err(SysbioError::invalid(
            "delta",
            "step fraction must lie in (0, 1)",
        ));
    }
    for &(lo, hi) in ranges {
        if lo >= hi {
            return Err(SysbioError::invalid("range", "each range needs lo < hi"));
        }
    }
    let mut rng = Rng::new(seed);
    let k = targets.len();
    // effects[param] collects one elementary effect per trajectory.
    let mut effects: Vec<Vec<f64>> = vec![Vec::with_capacity(n_trajectories); k];

    for _ in 0..n_trajectories {
        // Random base point in the unit hypercube.
        let base: Vec<f64> = (0..k).map(|_| rng.uniform()).collect();
        // Build the base model.
        let mut base_model = model.clone();
        for (i, t) in targets.iter().enumerate() {
            let (lo, hi) = ranges[i];
            base_model = t.apply(&base_model, lo + base[i] * (hi - lo))?;
        }
        let y_base = readout(&base_model)?;
        // Perturb each parameter once.
        for (i, t) in targets.iter().enumerate() {
            let (lo, hi) = ranges[i];
            // Step in the same direction unless it would leave [0,1].
            let dir = if base[i] + delta <= 1.0 { 1.0 } else { -1.0 };
            let stepped_unit = base[i] + dir * delta;
            let m_step = t.apply(&base_model, lo + stepped_unit * (hi - lo))?;
            let y_step = readout(&m_step)?;
            // Elementary effect normalised by the unit-cube step.
            let ee = (y_step - y_base) / (dir * delta);
            effects[i].push(ee);
        }
    }

    let mut out = Vec::with_capacity(k);
    for (i, t) in targets.iter().enumerate() {
        let es = &effects[i];
        let n = es.len() as f64;
        let mu = es.iter().sum::<f64>() / n;
        let mu_star = es.iter().map(|e| e.abs()).sum::<f64>() / n;
        let sigma = (es.iter().map(|e| (e - mu).powi(2)).sum::<f64>() / n).sqrt();
        out.push(GlobalSensitivity {
            param: t.label(),
            mu_star,
            mu,
            sigma,
        });
    }
    Ok(out)
}

/// Convenience: equal-width `(lo, hi)` ranges centred on each target's
/// nominal value at `±frac` (e.g. `frac = 0.5` → `[0.5·p, 1.5·p]`).
pub fn relative_ranges(
    model: &Model,
    targets: &[ParamTarget],
    frac: f64,
) -> Result<Vec<(f64, f64)>> {
    let mut out = Vec::with_capacity(targets.len());
    for t in targets {
        let p = t.read(model)?;
        let half = frac * p.abs().max(1e-9);
        out.push((p - half, p + half));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RateLaw, Reaction, Species};
    use crate::ode::{steady_state, OdeSystem};

    /// 0 ->(s) A ->(k A) 0. Steady A* = s/k.
    fn source_decay() -> Model {
        let mut m = Model::new("sd");
        let a = m.add_species(Species::new("A", 0.0));
        m.add_reaction(Reaction {
            id: "src".into(),
            reactants: vec![],
            products: vec![(a, 1.0)],
            rate_law: RateLaw::Constant { rate: 4.0 },
            reversible: false,
        });
        m.add_reaction(Reaction {
            id: "dec".into(),
            reactants: vec![(a, 1.0)],
            products: vec![],
            rate_law: RateLaw::MassAction {
                k: 2.0,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        m
    }

    fn steady_a(m: &Model) -> Result<f64> {
        let sys = OdeSystem::from_model(m);
        Ok(steady_state(&sys, &[1.0], 1e-10, 300)?.state[0])
    }

    #[test]
    fn local_sensitivity_of_source_rate() {
        // A* = s/k. ∂A*/∂s = 1/k = 0.5 at k=2. Scaled = (s/A*)(1/k) = 1.
        let m = source_decay();
        let sens = local_sensitivity(
            &m,
            &[ParamTarget::ConstantRate { reaction: 0 }],
            1e-3,
            steady_a,
        )
        .unwrap();
        assert!(
            (sens[0].derivative - 0.5).abs() < 1e-4,
            "{}",
            sens[0].derivative
        );
        assert!((sens[0].scaled - 1.0).abs() < 1e-3, "{}", sens[0].scaled);
    }

    #[test]
    fn local_sensitivity_of_decay_constant_is_negative() {
        // A* = s/k. ∂A*/∂k = -s/k^2 < 0. Scaled control coefficient -1.
        let m = source_decay();
        let sens = local_sensitivity(
            &m,
            &[ParamTarget::MassActionK { reaction: 1 }],
            1e-3,
            steady_a,
        )
        .unwrap();
        assert!(sens[0].derivative < 0.0);
        assert!((sens[0].scaled - (-1.0)).abs() < 1e-3, "{}", sens[0].scaled);
    }

    #[test]
    fn global_sensitivity_ranks_influential_parameter_higher() {
        // Readout = A* depends on both s and k; both should register
        // non-zero mu_star.
        let m = source_decay();
        let targets = vec![
            ParamTarget::ConstantRate { reaction: 0 },
            ParamTarget::MassActionK { reaction: 1 },
        ];
        let ranges = vec![(2.0, 8.0), (1.0, 4.0)];
        let g = global_sensitivity(&m, &targets, &ranges, 60, 0.1, 123, steady_a).unwrap();
        assert_eq!(g.len(), 2);
        assert!(g[0].mu_star > 0.0);
        assert!(g[1].mu_star > 0.0);
        // The decay constant enters non-linearly (1/k) -> larger sigma
        // than the linear source term.
        assert!(g[1].sigma > 0.0);
    }

    #[test]
    fn relative_ranges_centre_on_nominal() {
        let m = source_decay();
        let r = relative_ranges(&m, &[ParamTarget::ConstantRate { reaction: 0 }], 0.5).unwrap();
        // Nominal s = 4 -> [2, 6].
        assert!((r[0].0 - 2.0).abs() < 1e-9);
        assert!((r[0].1 - 6.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_arguments() {
        let m = source_decay();
        let t = vec![ParamTarget::ConstantRate { reaction: 0 }];
        assert!(local_sensitivity(&m, &t, 0.0, steady_a).is_err());
        // mismatched ranges length
        assert!(global_sensitivity(&m, &t, &[], 5, 0.1, 0, steady_a).is_err());
        // delta out of (0,1)
        assert!(global_sensitivity(&m, &t, &[(1.0, 2.0)], 5, 1.5, 0, steady_a).is_err());
    }
}
