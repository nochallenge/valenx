//! Ensemble stochastic runs and statistics — feature 17.
//!
//! A single stochastic trajectory is one sample path of a random
//! process; to characterise the *distribution* of behaviours you run
//! many trajectories and aggregate. [`run_ensemble`] launches `n_runs`
//! independent simulations — each with a distinct, deterministically
//! derived seed — samples every trajectory onto a shared uniform time
//! grid, and returns an [`EnsembleStats`] holding, per species and per
//! grid point, the **mean**, **variance**, **standard deviation** and
//! a configurable set of **percentiles** across the ensemble.
//!
//! Because the per-run seeds are derived from one master seed by a
//! fixed mixing rule, an entire ensemble is reproducible from that
//! single seed — the same guarantee the [`Rng`] gives a single run.

use crate::error::{Result, SysbioError};
use crate::stochastic::rng::Rng;
use crate::stochastic::ssa::{StochasticModel, StochasticTrace};

/// Which exact / approximate algorithm the ensemble should use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SsaMethod {
    /// Gillespie direct method (exact).
    Direct,
    /// Gibson-Bruck next-reaction method (exact).
    NextReaction,
    /// Explicit tau-leaping (approximate) with the given leap.
    TauLeap,
}

/// Aggregated statistics over an ensemble of trajectories.
///
/// Every field is indexed `[species][grid_point]`.
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleStats {
    /// The shared uniform sampling grid.
    pub grid: Vec<f64>,
    /// Ensemble mean of each species at each grid point.
    pub mean: Vec<Vec<f64>>,
    /// Ensemble variance (population variance, ÷N).
    pub variance: Vec<Vec<f64>>,
    /// Ensemble standard deviation.
    pub std_dev: Vec<Vec<f64>>,
    /// Requested percentile levels (e.g. `[5.0, 50.0, 95.0]`).
    pub percentile_levels: Vec<f64>,
    /// `percentiles[level][species][grid_point]`.
    pub percentiles: Vec<Vec<Vec<f64>>>,
    /// Number of trajectories aggregated.
    pub n_runs: usize,
}

impl EnsembleStats {
    /// The mean trajectory of species `i` as a `(time, mean)` series.
    pub fn mean_series(&self, i: usize) -> Vec<(f64, f64)> {
        self.grid
            .iter()
            .copied()
            .zip(self.mean[i].iter().copied())
            .collect()
    }
}

/// Run `n_runs` stochastic simulations and aggregate them onto a
/// uniform grid of `n_grid + 1` points spanning `[0, t_end]`.
///
/// `percentile_levels` are percentile ranks in `[0, 100]`. `tau` is
/// used only when `method` is [`SsaMethod::TauLeap`].
#[allow(clippy::too_many_arguments)]
pub fn run_ensemble(
    model: &StochasticModel,
    t_end: f64,
    n_runs: usize,
    n_grid: usize,
    method: SsaMethod,
    tau: f64,
    master_seed: u64,
    percentile_levels: &[f64],
) -> Result<EnsembleStats> {
    if t_end <= 0.0 {
        return Err(SysbioError::invalid("t_end", "t_end must be positive"));
    }
    if n_runs == 0 {
        return Err(SysbioError::invalid("n_runs", "need at least one run"));
    }
    if n_grid == 0 {
        return Err(SysbioError::invalid("n_grid", "need at least one interval"));
    }
    for &p in percentile_levels {
        if !(0.0..=100.0).contains(&p) {
            return Err(SysbioError::invalid(
                "percentile",
                "percentile levels must lie in [0, 100]",
            ));
        }
    }
    let n_species = model.initial_counts.len();
    let grid: Vec<f64> = (0..=n_grid)
        .map(|k| t_end * k as f64 / n_grid as f64)
        .collect();

    // Per-run seeds derived deterministically from the master seed.
    let mut seed_gen = Rng::new(master_seed);
    let max_steps = 5_000_000;

    // samples[species][grid_point] accumulates one value per run.
    let mut samples: Vec<Vec<Vec<f64>>> =
        vec![vec![Vec::with_capacity(n_runs); grid.len()]; n_species];

    for _ in 0..n_runs {
        let seed = seed_gen.next_u64();
        let trace: StochasticTrace = match method {
            SsaMethod::Direct => model.gillespie(t_end, seed, max_steps)?,
            SsaMethod::NextReaction => model.next_reaction(t_end, seed, max_steps)?,
            SsaMethod::TauLeap => {
                if tau <= 0.0 {
                    return Err(SysbioError::invalid("tau", "leap must be positive"));
                }
                model.tau_leap(t_end, tau, seed, max_steps)?
            }
        };
        for (sp, sp_samples) in samples.iter_mut().enumerate() {
            for (gp, &t) in grid.iter().enumerate() {
                sp_samples[gp].push(trace.count_at(sp, t) as f64);
            }
        }
    }

    // Reduce.
    let mut mean = vec![vec![0.0; grid.len()]; n_species];
    let mut variance = vec![vec![0.0; grid.len()]; n_species];
    let mut std_dev = vec![vec![0.0; grid.len()]; n_species];
    let mut percentiles = vec![vec![vec![0.0; grid.len()]; n_species]; percentile_levels.len()];

    for sp in 0..n_species {
        for gp in 0..grid.len() {
            let xs = &mut samples[sp][gp];
            let n = xs.len() as f64;
            let m = xs.iter().sum::<f64>() / n;
            let var = xs.iter().map(|x| (x - m).powi(2)).sum::<f64>() / n;
            mean[sp][gp] = m;
            variance[sp][gp] = var;
            std_dev[sp][gp] = var.sqrt();
            xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for (li, &lvl) in percentile_levels.iter().enumerate() {
                percentiles[li][sp][gp] = percentile_of_sorted(xs, lvl);
            }
        }
    }

    Ok(EnsembleStats {
        grid,
        mean,
        variance,
        std_dev,
        percentile_levels: percentile_levels.to_vec(),
        percentiles,
        n_runs,
    })
}

/// The `p`-th percentile (`p` in `[0, 100]`) of an already-sorted
/// slice, by linear interpolation between order statistics.
fn percentile_of_sorted(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let rank = (p / 100.0) * (sorted.len() - 1) as f64;
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        let w = rank - lo as f64;
        sorted[lo] * (1.0 - w) + sorted[hi] * w
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Model, RateLaw, Reaction, Species};

    fn decay(n0: i64, k: f64) -> StochasticModel {
        let mut m = Model::new("decay");
        let a = m.add_species(Species::new("A", n0 as f64));
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
        StochasticModel::from_model(&m).unwrap()
    }

    #[test]
    fn ensemble_mean_tracks_exponential_decay() {
        let m = decay(500, 1.0);
        let stats = run_ensemble(
            &m,
            3.0,
            80,
            30,
            SsaMethod::Direct,
            0.0,
            42,
            &[5.0, 50.0, 95.0],
        )
        .unwrap();
        // First grid point is t=0 -> exactly N0.
        assert!((stats.mean[0][0] - 500.0).abs() < 1e-9);
        // Mean is non-increasing along the grid (decay).
        for w in stats.mean[0].windows(2) {
            assert!(w[1] <= w[0] + 5.0);
        }
        // Mid-grid mean tracks the deterministic curve.
        let t_mid = stats.grid[15];
        let expect = 500.0 * (-t_mid).exp();
        assert!((stats.mean[0][15] - expect).abs() < 40.0);
    }

    #[test]
    fn percentiles_are_ordered() {
        let m = decay(300, 1.0);
        let stats = run_ensemble(
            &m,
            2.0,
            100,
            10,
            SsaMethod::Direct,
            0.0,
            7,
            &[10.0, 50.0, 90.0],
        )
        .unwrap();
        for gp in 0..stats.grid.len() {
            let p10 = stats.percentiles[0][0][gp];
            let p50 = stats.percentiles[1][0][gp];
            let p90 = stats.percentiles[2][0][gp];
            assert!(p10 <= p50 + 1e-9 && p50 <= p90 + 1e-9);
        }
    }

    #[test]
    fn variance_is_zero_at_deterministic_start() {
        let m = decay(100, 1.0);
        let stats = run_ensemble(&m, 1.0, 20, 5, SsaMethod::NextReaction, 0.0, 1, &[]).unwrap();
        // Every run starts at exactly N0 -> zero variance there.
        assert!(stats.variance[0][0].abs() < 1e-12);
        // Later, variance grows.
        assert!(stats.variance[0][5] > 0.0);
    }

    #[test]
    fn tau_leap_ensemble_runs() {
        let m = decay(400, 1.0);
        let stats = run_ensemble(&m, 2.0, 30, 8, SsaMethod::TauLeap, 0.02, 9, &[50.0]).unwrap();
        assert_eq!(stats.n_runs, 30);
        assert!((stats.mean[0][0] - 400.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_arguments() {
        let m = decay(10, 1.0);
        assert!(run_ensemble(&m, 1.0, 0, 5, SsaMethod::Direct, 0.0, 0, &[]).is_err());
        assert!(run_ensemble(&m, 1.0, 5, 0, SsaMethod::Direct, 0.0, 0, &[]).is_err());
        assert!(run_ensemble(&m, 1.0, 5, 5, SsaMethod::Direct, 0.0, 0, &[150.0]).is_err());
    }

    #[test]
    fn percentile_of_sorted_endpoints() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(percentile_of_sorted(&xs, 0.0), 1.0);
        assert_eq!(percentile_of_sorted(&xs, 100.0), 5.0);
        assert_eq!(percentile_of_sorted(&xs, 50.0), 3.0);
    }
}
