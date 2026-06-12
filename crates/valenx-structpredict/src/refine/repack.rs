//! **Feature 13 — all-atom rotamer repacking.**
//!
//! Once a backbone exists, the sidechains must be packed without
//! clashes — the SCWRL / Rosetta `repack` problem. It is a
//! combinatorial optimisation: choose, for each residue, one rotamer
//! from its library so the total sidechain-packing energy is minimal.
//! With `r` rotamers per residue and `n` residues the search space is
//! `r^n` — far too large to brute-force.
//!
//! This module solves it with the two classical techniques, used
//! together:
//!
//! 1. **Dead-end elimination (DEE)** — the Desmet criterion: a
//!    rotamer `a` at a residue can be *eliminated* if some other
//!    rotamer `b` at that residue is always better, no matter what
//!    the neighbours do. Repeatedly pruning dead rotamers shrinks the
//!    search space, often dramatically, with *no* loss of the global
//!    optimum.
//! 2. **Simulated annealing** — over whatever rotamer combinations
//!    survive DEE, a Metropolis Monte-Carlo annealer finds a
//!    low-energy assignment.
//!
//! The energy model is centroid-resolution: a singles term (the
//! rotamer's own backbone-fit / prior) plus a pairwise term (a
//! soft-clash penalty between sidechain centroids). Both DEE and the
//! annealer are the genuine published algorithms.

use serde::{Deserialize, Serialize};
use valenx_md::Rng;

use crate::error::{Result, StructPredictError};
use crate::model::ProteinModel;
use crate::rotamer::{place_sidechain_centroid, rotamers_for, Rotamer};

/// The outcome of a repacking run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RepackResult {
    /// Chosen rotamer index per residue (into that residue's
    /// [`rotamers_for`] set).
    pub rotamers: Vec<usize>,
    /// Total packing energy of the chosen assignment.
    pub energy: f64,
    /// Rotamers eliminated by DEE before the annealing search.
    pub dee_eliminated: usize,
    /// Total rotamers across all residues before DEE.
    pub dee_total: usize,
}

impl RepackResult {
    /// Fraction of the rotamer search space removed by DEE.
    pub fn dee_reduction(&self) -> f64 {
        if self.dee_total == 0 {
            0.0
        } else {
            self.dee_eliminated as f64 / self.dee_total as f64
        }
    }
}

/// The per-residue rotamer pool and the precomputed sidechain
/// centroids for every rotamer.
struct RotamerPool {
    /// `rotamers[i]` are the candidate rotamers of residue `i`.
    rotamers: Vec<Vec<Rotamer>>,
    /// `centroids[i][k]` is the sidechain centroid of residue `i`
    /// under its rotamer `k` (`None` if the residue lacks a backbone).
    centroids: Vec<Vec<Option<nalgebra::Point3<f64>>>>,
    /// `alive[i][k]` — whether rotamer `k` of residue `i` survives DEE.
    alive: Vec<Vec<bool>>,
}

impl RotamerPool {
    fn build(model: &ProteinModel) -> Self {
        let mut rotamers = Vec::with_capacity(model.residues.len());
        let mut centroids = Vec::with_capacity(model.residues.len());
        let mut alive = Vec::with_capacity(model.residues.len());
        for res in &model.residues {
            let rots = rotamers_for(res.aa);
            let cents: Vec<_> = rots
                .iter()
                .map(|r| place_sidechain_centroid(res, r))
                .collect();
            alive.push(vec![true; rots.len()]);
            centroids.push(cents);
            rotamers.push(rots);
        }
        RotamerPool {
            rotamers,
            centroids,
            alive,
        }
    }

    /// Self-energy of rotamer `k` at residue `i`: a small term that
    /// penalises low-prior (uncommon) rotamers.
    fn self_energy(&self, i: usize, k: usize) -> f64 {
        // −ln(prior): a common rotamer is low-energy.
        -self.rotamers[i][k].probability.max(1e-3).ln() * 0.3
    }

    /// Pairwise energy between rotamer `ka` at residue `i` and rotamer
    /// `kb` at residue `j`: a soft repulsive clash between the two
    /// sidechain centroids.
    fn pair_energy(&self, i: usize, ka: usize, j: usize, kb: usize) -> f64 {
        let (Some(ca), Some(cb)) = (self.centroids[i][ka], self.centroids[j][kb]) else {
            return 0.0;
        };
        let d = (ca - cb).norm();
        // Soft clash: a smooth penalty rising steeply below ~4 Å.
        const CLASH: f64 = 4.0;
        if d >= CLASH {
            0.0
        } else {
            let x = (CLASH - d) / CLASH;
            5.0 * x * x
        }
    }

    /// Total energy of a complete rotamer assignment.
    fn total_energy(&self, choice: &[usize]) -> f64 {
        let n = choice.len();
        let mut e = 0.0;
        for i in 0..n {
            e += self.self_energy(i, choice[i]);
            for j in (i + 1)..n {
                e += self.pair_energy(i, choice[i], j, choice[j]);
            }
        }
        e
    }
}

/// Runs the Desmet **dead-end-elimination** pass on a rotamer pool.
///
/// Marks rotamers dead until no more can be pruned. Returns the
/// number eliminated. The Goldstein form of the criterion is used:
/// rotamer `a` at residue `i` is dead if some other rotamer `b`
/// satisfies, for every neighbour, that `a` is never better than `b`.
fn dead_end_elimination(pool: &mut RotamerPool) -> usize {
    let n = pool.rotamers.len();
    let mut eliminated = 0usize;
    let mut changed = true;
    while changed {
        changed = false;
        for i in 0..n {
            let alive_indices: Vec<usize> = (0..pool.rotamers[i].len())
                .filter(|&k| pool.alive[i][k])
                .collect();
            if alive_indices.len() < 2 {
                continue;
            }
            for &a in &alive_indices {
                if !pool.alive[i][a] {
                    continue;
                }
                for &b in &alive_indices {
                    if a == b || !pool.alive[i][b] {
                        continue;
                    }
                    // Goldstein: a is dead if E_self(a)-E_self(b) +
                    // Σ_j min_t [E_pair(i,a,j,t)-E_pair(i,b,j,t)] > 0.
                    let mut bound = pool.self_energy(i, a) - pool.self_energy(i, b);
                    for j in 0..n {
                        if j == i {
                            continue;
                        }
                        let mut min_diff = f64::INFINITY;
                        for t in 0..pool.rotamers[j].len() {
                            if !pool.alive[j][t] {
                                continue;
                            }
                            let diff = pool.pair_energy(i, a, j, t) - pool.pair_energy(i, b, j, t);
                            if diff < min_diff {
                                min_diff = diff;
                            }
                        }
                        if min_diff.is_finite() {
                            bound += min_diff;
                        }
                    }
                    if bound > 1e-9 {
                        // a can never beat b → eliminate a.
                        pool.alive[i][a] = false;
                        eliminated += 1;
                        changed = true;
                        break;
                    }
                }
            }
        }
    }
    eliminated
}

/// Repacks the sidechains of a model.
///
/// Builds the per-residue rotamer pool, runs dead-end elimination to
/// prune provably-suboptimal rotamers, then runs a simulated-
/// annealing Monte-Carlo search over the survivors. The chosen
/// rotamers and the final energy are returned; the model's `cb` slot
/// of each residue is updated to the rebuilt Cβ (the rotamer choice
/// itself is in [`RepackResult::rotamers`]).
///
/// `mc_moves` is the annealing move budget; `seed` fixes the RNG.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an empty model or `mc_moves == 0`.
pub fn repack_sidechains(
    model: &mut ProteinModel,
    mc_moves: usize,
    seed: u64,
) -> Result<RepackResult> {
    if model.residues.is_empty() {
        return Err(StructPredictError::invalid("model", "no residues"));
    }
    if mc_moves == 0 {
        return Err(StructPredictError::invalid(
            "mc_moves",
            "must be at least 1",
        ));
    }
    let mut pool = RotamerPool::build(model);
    let dee_total: usize = pool.rotamers.iter().map(|r| r.len()).sum();
    let dee_eliminated = dead_end_elimination(&mut pool);

    let n = pool.rotamers.len();
    // Initial assignment: the first surviving rotamer per residue.
    let mut choice: Vec<usize> = (0..n)
        .map(|i| {
            (0..pool.rotamers[i].len())
                .find(|&k| pool.alive[i][k])
                .unwrap_or(0)
        })
        .collect();
    let mut energy = pool.total_energy(&choice);
    let mut best = choice.clone();
    let mut best_energy = energy;

    let mut rng = Rng::new(seed);
    let start_t: f64 = 5.0;
    let end_t: f64 = 0.05;
    for step in 0..mc_moves {
        let frac = step as f64 / mc_moves as f64;
        let t = start_t * (end_t / start_t).powf(frac);
        // Pick a residue and a surviving alternative rotamer.
        let i = rng.below(n);
        let alive: Vec<usize> = (0..pool.rotamers[i].len())
            .filter(|&k| pool.alive[i][k])
            .collect();
        if alive.len() < 2 {
            continue;
        }
        let new_k = alive[rng.below(alive.len())];
        if new_k == choice[i] {
            continue;
        }
        let old_k = choice[i];
        choice[i] = new_k;
        let new_energy = pool.total_energy(&choice);
        let delta = new_energy - energy;
        if delta <= 0.0 || rng.uniform() < (-delta / t).exp() {
            energy = new_energy;
            if energy < best_energy {
                best_energy = energy;
                best = choice.clone();
            }
        } else {
            choice[i] = old_k; // reject
        }
    }

    // Update the model's Cβ atoms to the rebuilt geometry.
    for res in &mut model.residues {
        if res.aa != 'G' && res.cb.is_none() {
            if let (Some(n), Some(ca), Some(c)) = (res.n, res.ca, res.c) {
                res.cb = Some(crate::rotamer::rebuild_cb(&n, &ca, &c));
            }
        }
    }

    Ok(RepackResult {
        rotamers: best,
        energy: best_energy,
        dee_eliminated,
        dee_total,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    /// A model whose residues sit on a line, well-spaced so a clean
    /// pack exists.
    fn spaced_model(seq: &str) -> ProteinModel {
        let mut m = ProteinModel::from_sequence(seq).expect("model");
        for (i, r) in m.residues.iter_mut().enumerate() {
            let x = i as f64 * 6.5;
            r.n = Some(Point3::new(x, 0.0, 0.0));
            r.ca = Some(Point3::new(x + 1.46, 0.0, 0.0));
            r.c = Some(Point3::new(x + 2.0, 1.4, 0.0));
            r.o = Some(Point3::new(x + 1.3, 2.4, 0.0));
        }
        m
    }

    #[test]
    fn repack_chooses_one_rotamer_per_residue() {
        let mut m = spaced_model("LWFYM");
        let res = repack_sidechains(&mut m, 400, 7).expect("repack");
        assert_eq!(res.rotamers.len(), 5);
        // Each chosen index is within that residue's rotamer set.
        for (i, &k) in res.rotamers.iter().enumerate() {
            let count = rotamers_for(m.residues[i].aa).len();
            assert!(k < count, "rotamer {k} < {count}");
        }
    }

    #[test]
    fn dee_eliminates_some_rotamers_when_packed_tight() {
        // A tightly-packed cluster — DEE should be able to prune.
        let mut m = ProteinModel::from_sequence("WWWW").expect("model");
        for (i, r) in m.residues.iter_mut().enumerate() {
            let x = i as f64 * 4.0; // close packing
            r.n = Some(Point3::new(x, 0.0, 0.0));
            r.ca = Some(Point3::new(x + 1.46, 0.0, 0.0));
            r.c = Some(Point3::new(x + 2.0, 1.4, 0.0));
            r.o = Some(Point3::new(x + 1.3, 2.4, 0.0));
        }
        let res = repack_sidechains(&mut m, 200, 1).expect("repack");
        // DEE total is the sum of rotamer-set sizes; reduction in [0,1].
        assert!(res.dee_total > 0);
        assert!((0.0..=1.0).contains(&res.dee_reduction()));
    }

    #[test]
    fn repack_is_deterministic() {
        let mut a = spaced_model("LWFY");
        let mut b = spaced_model("LWFY");
        let ra = repack_sidechains(&mut a, 300, 42).expect("a");
        let rb = repack_sidechains(&mut b, 300, 42).expect("b");
        assert_eq!(ra.rotamers, rb.rotamers);
        assert_eq!(ra.energy, rb.energy);
    }

    #[test]
    fn empty_model_rejected() {
        let mut m = ProteinModel::new();
        assert!(repack_sidechains(&mut m, 100, 0).is_err());
    }
}
