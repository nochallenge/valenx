//! Stochastic simulation algorithms — features 14, 15 and 16.
//!
//! The discrete-stochastic half of the simulation layer. Where the ODE
//! layer treats amounts as continuous concentrations, these
//! algorithms treat them as integer **molecule counts** and simulate
//! the Markov jump process exactly (SSA, next-reaction) or
//! approximately (tau-leaping).
//!
//! All three consume a [`StochasticModel`] — a reaction network whose
//! reactions have integer stoichiometry and a *propensity* function.
//! The propensity of a reaction is its stochastic rate: for a
//! mass-action reaction `A + B → …` with rate constant `c` it is
//! `c · #A · #B` (combinatorial — `c · #A·(#A−1)/2` for `2A → …`).
//! A non-mass-action [`RateLaw`] is evaluated directly on the counts.
//!
//! - [`StochasticModel::gillespie`] — the Gillespie **direct method**
//!   SSA (feature 14): exact, one reaction per step.
//! - [`StochasticModel::tau_leap`] — explicit **tau-leaping**
//!   (feature 15): fires many reactions per step from Poisson draws,
//!   with a negative-population guard.
//! - [`StochasticModel::next_reaction`] — the Gibson-Bruck
//!   **next-reaction method** (feature 16): exact, with an indexed
//!   priority queue of putative reaction times and a dependency graph
//!   so only affected propensities are recomputed.

use crate::error::{Result, SysbioError};
use crate::model::{Model, RateLaw};
use crate::stochastic::rng::Rng;

/// A reaction-network model prepared for discrete-stochastic
/// simulation.
///
/// Built from a [`Model`] via [`StochasticModel::from_model`]; the
/// initial amounts are rounded to integer molecule counts.
#[derive(Debug, Clone)]
pub struct StochasticModel {
    /// Integer initial molecule counts, one per species.
    pub initial_counts: Vec<i64>,
    /// Per-reaction net stoichiometry as `(species, delta)` pairs.
    pub stoich: Vec<Vec<(usize, i64)>>,
    /// Per-reaction rate law (used to derive the propensity).
    pub rate_laws: Vec<RateLaw>,
    /// Per-species `constant` flag — a constant species' count is held
    /// fixed (it still contributes to propensities).
    pub constant: Vec<bool>,
}

impl StochasticModel {
    /// Build a stochastic model from a reaction network. Initial
    /// amounts are rounded to the nearest non-negative integer.
    pub fn from_model(model: &Model) -> Result<Self> {
        model.validate()?;
        let initial_counts = model
            .species
            .iter()
            .map(|s| s.initial.round().max(0.0) as i64)
            .collect();
        let stoich = model
            .reactions
            .iter()
            .map(|r| {
                r.net_stoichiometry()
                    .into_iter()
                    .map(|(i, c)| (i, c.round() as i64))
                    .filter(|&(_, c)| c != 0)
                    .collect()
            })
            .collect();
        Ok(StochasticModel {
            initial_counts,
            stoich,
            rate_laws: model.reactions.iter().map(|r| r.rate_law.clone()).collect(),
            constant: model.species.iter().map(|s| s.constant).collect(),
        })
    }

    /// The propensity of reaction `j` at integer counts `x`.
    ///
    /// Mass-action laws get the *combinatorial* count term
    /// (`#A·(#A−1)·…` for higher orders) — the correct stochastic
    /// propensity. Other laws are evaluated on the counts as a
    /// floating-point amount, the standard hybrid treatment.
    fn propensity(&self, j: usize, x: &[i64]) -> f64 {
        match &self.rate_laws[j] {
            RateLaw::Constant { rate } => rate.max(0.0),
            RateLaw::MassAction { k, reactants } => {
                let mut a = *k;
                for &(idx, order) in reactants {
                    let n = *x.get(idx).unwrap_or(&0);
                    if n <= 0 {
                        return 0.0;
                    }
                    let ord = order.round() as i64;
                    if ord <= 1 {
                        a *= n as f64;
                    } else {
                        // Falling factorial n·(n-1)·…·(n-ord+1).
                        let mut term = 1.0;
                        for d in 0..ord {
                            term *= (n - d) as f64;
                        }
                        a *= term;
                    }
                }
                a.max(0.0)
            }
            other => {
                let amounts: Vec<f64> = x.iter().map(|&v| v as f64).collect();
                other.rate(&amounts).max(0.0)
            }
        }
    }

    /// Apply reaction `j`'s stoichiometry to the count vector in place,
    /// respecting `constant` species.
    fn fire(&self, j: usize, x: &mut [i64]) {
        for &(i, d) in &self.stoich[j] {
            if !self.constant[i] {
                x[i] += d;
            }
        }
    }

    /// The Gillespie direct-method SSA (feature 14).
    ///
    /// Exact simulation of the chemical master equation. Each step:
    /// (1) compute all propensities and their sum `a0`; (2) draw the
    /// time to the next reaction as `Exp(a0)`; (3) pick which reaction
    /// fires with probability `a_j/a0`; (4) update counts. Stops at
    /// `t_end` or when every propensity is zero (the system is
    /// absorbed). `max_steps` bounds the run.
    pub fn gillespie(
        &self,
        t_end: f64,
        seed: u64,
        max_steps: usize,
    ) -> Result<StochasticTrace> {
        if t_end <= 0.0 {
            return Err(SysbioError::invalid("t_end", "t_end must be positive"));
        }
        let mut rng = Rng::new(seed);
        let mut t = 0.0;
        let mut x = self.initial_counts.clone();
        let mut trace = StochasticTrace::seeded(&x);
        let nr = self.rate_laws.len();

        for _ in 0..max_steps {
            let mut props = vec![0.0; nr];
            let mut a0 = 0.0;
            for (j, p) in props.iter_mut().enumerate() {
                *p = self.propensity(j, &x);
                a0 += *p;
            }
            if a0 <= 0.0 {
                break; // absorbed — no reaction can fire
            }
            let dt = rng.exponential(a0);
            if t + dt > t_end {
                break;
            }
            t += dt;
            // Select the firing reaction.
            let threshold = rng.uniform() * a0;
            let mut cum = 0.0;
            let mut chosen = nr - 1;
            for (j, &p) in props.iter().enumerate() {
                cum += p;
                if cum >= threshold {
                    chosen = j;
                    break;
                }
            }
            self.fire(chosen, &mut x);
            trace.push(t, &x);
        }
        // Pin the final state at t_end so trajectories line up.
        trace.push(t_end, &x);
        Ok(trace)
    }

    /// Explicit tau-leaping (feature 15).
    ///
    /// Advances time by a fixed leap `tau`; in each leap every
    /// reaction `j` fires `Poisson(a_j · tau)` times. This is much
    /// faster than the SSA when propensities are large, at the cost of
    /// approximation. The classic failure mode — a reaction firing so
    /// many times a species count goes negative — is handled by a
    /// **negative-population guard**: a leap that would drive any
    /// count below zero is rejected, `tau` is halved, and the leap is
    /// retried; after several rejections the leap is taken with the
    /// firing counts clamped so no count crosses zero.
    pub fn tau_leap(
        &self,
        t_end: f64,
        tau: f64,
        seed: u64,
        max_steps: usize,
    ) -> Result<StochasticTrace> {
        if t_end <= 0.0 {
            return Err(SysbioError::invalid("t_end", "t_end must be positive"));
        }
        if tau <= 0.0 {
            return Err(SysbioError::invalid("tau", "leap must be positive"));
        }
        let mut rng = Rng::new(seed);
        let mut t = 0.0;
        let mut x = self.initial_counts.clone();
        let mut trace = StochasticTrace::seeded(&x);
        let nr = self.rate_laws.len();

        for _ in 0..max_steps {
            if t >= t_end {
                break;
            }
            let mut step = tau.min(t_end - t);
            // Try the leap, shrinking on a negative-population hit.
            let mut applied: Option<Vec<i64>> = None;
            for attempt in 0..10 {
                let mut firings = vec![0u64; nr];
                for (j, f) in firings.iter_mut().enumerate() {
                    let a = self.propensity(j, &x);
                    *f = rng.poisson(a * step);
                }
                let mut trial = x.clone();
                for (j, &count) in firings.iter().enumerate() {
                    for &(i, d) in &self.stoich[j] {
                        if !self.constant[i] {
                            trial[i] += d * count as i64;
                        }
                    }
                }
                if trial.iter().all(|&v| v >= 0) {
                    applied = Some(trial);
                    break;
                }
                // Negative-population hit: shrink the leap and retry,
                // or — on the final attempt — clamp at zero.
                if attempt < 9 {
                    step *= 0.5;
                } else {
                    for v in trial.iter_mut() {
                        if *v < 0 {
                            *v = 0;
                        }
                    }
                    applied = Some(trial);
                }
            }
            x = applied.expect("leap always resolves");
            t += step;
            trace.push(t, &x);
        }
        Ok(trace)
    }

    /// The Gibson-Bruck next-reaction method (feature 16).
    ///
    /// An exact algorithm equivalent to the SSA but asymptotically
    /// faster on large sparse networks. It keeps an *absolute*
    /// putative firing time for every reaction in an indexed priority
    /// queue; each step fires the soonest reaction and — using a
    /// precomputed **dependency graph** — recomputes only the
    /// propensities of reactions actually affected by the change,
    /// rescaling their putative times rather than redrawing them.
    pub fn next_reaction(
        &self,
        t_end: f64,
        seed: u64,
        max_steps: usize,
    ) -> Result<StochasticTrace> {
        if t_end <= 0.0 {
            return Err(SysbioError::invalid("t_end", "t_end must be positive"));
        }
        let mut rng = Rng::new(seed);
        let nr = self.rate_laws.len();
        let mut x = self.initial_counts.clone();
        let mut trace = StochasticTrace::seeded(&x);

        // Dependency graph: firing reaction j affects reaction k if k's
        // propensity reads a species j changes.
        let dep = self.dependency_graph();

        // Putative times and current propensities.
        let mut props = vec![0.0; nr];
        let mut times = vec![f64::INFINITY; nr];
        for j in 0..nr {
            props[j] = self.propensity(j, &x);
            times[j] = rng.exponential(props[j]);
        }

        for _ in 0..max_steps {
            // Soonest reaction (indexed min scan — the "priority
            // queue"; for the v1 sizes a linear scan is the queue).
            let mut mu = 0;
            let mut tmin = times[0];
            for (j, &tj) in times.iter().enumerate() {
                if tj < tmin {
                    tmin = tj;
                    mu = j;
                }
            }
            if !tmin.is_finite() || tmin > t_end {
                break;
            }
            let t = tmin;
            self.fire(mu, &mut x);
            trace.push(t, &x);

            // Refresh affected reactions.
            for &k in &dep[mu] {
                let old_a = props[k];
                let new_a = self.propensity(k, &x);
                props[k] = new_a;
                if k == mu {
                    // The reaction that just fired always redraws.
                    times[k] = t + rng.exponential(new_a);
                } else if new_a > 0.0 {
                    if old_a > 0.0 && times[k].is_finite() {
                        // Gibson-Bruck rescaling of the remaining time.
                        times[k] = t + (old_a / new_a) * (times[k] - t);
                    } else {
                        times[k] = t + rng.exponential(new_a);
                    }
                } else {
                    times[k] = f64::INFINITY;
                }
            }
        }
        trace.push(t_end, &x);
        Ok(trace)
    }

    /// Reaction-to-reaction dependency graph: `graph[j]` lists every
    /// reaction whose propensity must be recomputed after `j` fires
    /// (always including `j` itself).
    fn dependency_graph(&self) -> Vec<Vec<usize>> {
        let nr = self.rate_laws.len();
        // Species changed by each reaction.
        let changed: Vec<Vec<usize>> = self
            .stoich
            .iter()
            .map(|s| s.iter().map(|&(i, _)| i).collect())
            .collect();
        // Species read by each reaction's propensity.
        let read: Vec<Vec<usize>> = self
            .rate_laws
            .iter()
            .map(|l| l.dependencies())
            .collect();
        let mut graph = vec![Vec::new(); nr];
        for (j, gj) in graph.iter_mut().enumerate() {
            for (k, rk) in read.iter().enumerate() {
                if k == j || changed[j].iter().any(|s| rk.contains(s)) {
                    gj.push(k);
                }
            }
        }
        graph
    }
}

/// A single stochastic trajectory: jump times and integer states.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StochasticTrace {
    /// Event times (non-decreasing).
    pub times: Vec<f64>,
    /// Integer count vectors, one per time.
    pub states: Vec<Vec<i64>>,
}

impl StochasticTrace {
    fn seeded(x0: &[i64]) -> Self {
        StochasticTrace {
            times: vec![0.0],
            states: vec![x0.to_vec()],
        }
    }

    fn push(&mut self, t: f64, x: &[i64]) {
        self.times.push(t);
        self.states.push(x.to_vec());
    }

    /// Number of recorded events.
    pub fn len(&self) -> usize {
        self.times.len()
    }

    /// Whether the trace is empty.
    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }

    /// The final count vector, if any.
    pub fn final_state(&self) -> Option<&[i64]> {
        self.states.last().map(|v| v.as_slice())
    }

    /// The count of species `i` interpolated (left-continuous step) at
    /// time `t` — a stochastic trajectory is piecewise-constant, so the
    /// value is that of the last event at or before `t`.
    pub fn count_at(&self, i: usize, t: f64) -> i64 {
        if self.is_empty() {
            return 0;
        }
        let mut val = self.states[0][i];
        for (k, &tk) in self.times.iter().enumerate() {
            if tk <= t {
                val = self.states[k][i];
            } else {
                break;
            }
        }
        val
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Reaction, Species};

    /// Pure decay A -> 0 with N0 molecules. Every algorithm must
    /// monotonically deplete A.
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

    /// A <-> B isomerisation: total A+B is conserved.
    fn isomerise() -> StochasticModel {
        let mut m = Model::new("iso");
        let a = m.add_species(Species::new("A", 100.0));
        let b = m.add_species(Species::new("B", 0.0));
        m.add_reaction(Reaction {
            id: "fwd".into(),
            reactants: vec![(a, 1.0)],
            products: vec![(b, 1.0)],
            rate_law: RateLaw::MassAction {
                k: 1.0,
                reactants: vec![(a, 1.0)],
            },
            reversible: false,
        });
        m.add_reaction(Reaction {
            id: "rev".into(),
            reactants: vec![(b, 1.0)],
            products: vec![(a, 1.0)],
            rate_law: RateLaw::MassAction {
                k: 1.0,
                reactants: vec![(b, 1.0)],
            },
            reversible: false,
        });
        StochasticModel::from_model(&m).unwrap()
    }

    #[test]
    fn gillespie_decay_is_monotone_and_bounded() {
        let m = decay(200, 1.0);
        let tr = m.gillespie(5.0, 1, 100_000).unwrap();
        // A only ever decreases.
        for w in tr.states.windows(2) {
            assert!(w[1][0] <= w[0][0]);
        }
        // After 5 mean-lifetimes very little is left.
        let final_a = tr.final_state().unwrap()[0];
        assert!(final_a < 20, "expected near-empty, got {final_a}");
    }

    #[test]
    fn gillespie_decay_mean_matches_deterministic() {
        // Ensemble of decay runs; mean at t should track N0 exp(-k t).
        let m = decay(1000, 1.0);
        let mut sum = 0i64;
        let runs = 60;
        for s in 0..runs {
            let tr = m.gillespie(1.0, 1000 + s, 200_000).unwrap();
            sum += tr.count_at(0, 1.0);
        }
        let mean = sum as f64 / runs as f64;
        let expect = 1000.0 * (-1.0_f64).exp();
        assert!((mean - expect).abs() < 60.0, "mean {mean}, expect {expect}");
    }

    #[test]
    fn tau_leap_conserves_total_mass() {
        let m = isomerise();
        let tr = m.tau_leap(3.0, 0.01, 7, 100_000).unwrap();
        for st in &tr.states {
            assert_eq!(st[0] + st[1], 100, "A+B not conserved: {st:?}");
        }
    }

    #[test]
    fn tau_leap_never_goes_negative() {
        // Aggressive leap on a small population — the guard must hold.
        let m = decay(30, 5.0);
        let tr = m.tau_leap(2.0, 0.5, 3, 10_000).unwrap();
        for st in &tr.states {
            assert!(st[0] >= 0, "negative count: {st:?}");
        }
    }

    #[test]
    fn next_reaction_decay_matches_gillespie_distribution() {
        let m = decay(200, 1.0);
        // Both exact methods: ensemble means should agree closely.
        let gill: f64 = (0..40)
            .map(|s| m.gillespie(1.0, s, 100_000).unwrap().count_at(0, 1.0) as f64)
            .sum::<f64>()
            / 40.0;
        let nrm: f64 = (0..40)
            .map(|s| {
                m.next_reaction(1.0, 5000 + s, 100_000)
                    .unwrap()
                    .count_at(0, 1.0) as f64
            })
            .sum::<f64>()
            / 40.0;
        assert!((gill - nrm).abs() < 25.0, "gill {gill} vs nrm {nrm}");
    }

    #[test]
    fn next_reaction_conserves_isomerisation_mass() {
        let m = isomerise();
        let tr = m.next_reaction(3.0, 11, 100_000).unwrap();
        for st in &tr.states {
            assert_eq!(st[0] + st[1], 100);
        }
    }

    #[test]
    fn dependency_graph_links_coupled_reactions() {
        // In A<->B both reactions read species the other changes.
        let m = isomerise();
        let g = m.dependency_graph();
        assert!(g[0].contains(&1));
        assert!(g[1].contains(&0));
    }

    #[test]
    fn rejects_bad_horizon() {
        let m = decay(10, 1.0);
        assert!(m.gillespie(-1.0, 0, 10).is_err());
        assert!(m.tau_leap(1.0, 0.0, 0, 10).is_err());
    }
}
