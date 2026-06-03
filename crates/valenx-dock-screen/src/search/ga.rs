//! Feature 9 — AutoDock 4-class Lamarckian genetic algorithm.
//!
//! AutoDock 4's flagship search (Morris et al. 1998, *J. Comput. Chem.*
//! **19**, 1639) is a *Lamarckian* genetic algorithm (LGA). A
//! population of ligand poses is evolved over generations by
//!
//! 1. **selection** — the fittest poses are more likely to reproduce
//!    (tournament selection here, AutoDock 4's default);
//! 2. **crossover** — two parents exchange genome segments to make a
//!    child (one-point or uniform);
//! 3. **mutation** — Cauchy-distributed perturbations of translation,
//!    rotation and torsion genes;
//! 4. **local search** — a fraction of the population is refined by a
//!    local optimiser, and crucially the *refined* genome is written
//!    back into the individual (the Lamarckian step — acquired
//!    improvements are inherited).
//!
//! A pose's genome is the flat real vector `[tx ty tz | rx ry rz | t₀
//! … t_N]` — translation, axis-angle rotation, and one entry per
//! rotatable bond — exactly the encoding [`valenx_dock`]'s BFGS
//! minimiser uses.
//!
//! ## Local search
//!
//! AutoDock 4's published default local search is **Solis-Wets**
//! (Solis & Wets 1981, *Math. Oper. Res.* **6**, 19) — a
//! gradient-free random adaptive direction search with a BFGS-style
//! step-size adaptation (the consecutive-success / consecutive-failure
//! expand / contract schedule). [`SolisWetsParams`] from
//! `crate::search::solis_wets` is the default; the
//! [`GaParams::local_search`] enum still lets callers pick the legacy
//! coordinate-descent path for stability under exotic objectives.
//!
//! ## Operator schedule
//!
//! - **Selection** — tournament selection with size `tournament_size`
//!   (AutoDock 4 default = 4; matches our default).
//! - **Mutation** — Cauchy mutation, the AutoDock 4 default. Each gene
//!   is, with probability `mutation_rate`, displaced by a Cauchy-
//!   distributed kick (heavy-tailed, allowing rare long jumps while
//!   most steps stay small) whose scale is per-DOF (translation /
//!   rotation / torsion).
//! - **Crossover** — arithmetic (uniform) crossover, gene-wise random
//!   pick from one parent or the other.
//! - **Lamarckian step** — for every `local_search_fraction` of the
//!   population per generation, run a Solis-Wets local search starting
//!   from that individual's genome and *replace* its genes with the
//!   refined result. This is the "Lamarckian" inheritance of acquired
//!   improvements that gave the LGA its name.
//! - **Elitism** — the current best individual always survives into
//!   the next generation unchanged.

use nalgebra::Vector3;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use valenx_dock::pose::Pose;
use valenx_dock::search::bfgs::{pose_to_vec, vec_to_pose};

use crate::error::{DockScreenError, Result};
use crate::search::objective::PoseObjective;
use crate::search::solis_wets::{solis_wets, SolisWetsParams};

/// Which local-search operator the LGA uses for the Lamarckian step.
#[derive(Clone, Copy, Debug)]
pub enum LocalSearchKind {
    /// AutoDock 4's published default — Solis-Wets random adaptive
    /// direction search with the published bias / expansion /
    /// contraction schedule. See `crate::search::solis_wets`.
    SolisWets(SolisWetsParams),
    /// The legacy coordinate-descent path — cheap, gradient-free, and
    /// numerically very forgiving. Useful when the objective is
    /// non-smooth in ways Solis-Wets's continuous adaptation does not
    /// like (e.g. heavy clipped clashes).
    CoordinateDescent,
}

impl Default for LocalSearchKind {
    fn default() -> Self {
        LocalSearchKind::SolisWets(SolisWetsParams::lamarckian())
    }
}

/// Tunable parameters for the Lamarckian GA.
#[derive(Clone, Copy, Debug)]
pub struct GaParams {
    /// Number of individuals in the population.
    pub population: usize,
    /// Number of generations to evolve.
    pub generations: usize,
    /// Probability `[0,1]` that a child gene is mutated.
    pub mutation_rate: f64,
    /// Probability `[0,1]` that two selected parents are crossed
    /// (otherwise the fitter parent is copied through).
    pub crossover_rate: f64,
    /// Fraction `[0,1]` of the population that gets a local-search
    /// refinement each generation (the Lamarckian step).
    pub local_search_fraction: f64,
    /// Tournament size for selection.
    pub tournament_size: usize,
    /// Which local-search operator the Lamarckian step uses.
    pub local_search: LocalSearchKind,
}

impl Default for GaParams {
    fn default() -> Self {
        GaParams {
            population: 150,           // AutoDock 4 default GA pop is 150
            generations: 27,           // AutoDock's historical default GA run count
            mutation_rate: 0.02,       // AutoDock 4 default
            crossover_rate: 0.8,       // AutoDock 4 default
            local_search_fraction: 0.06, // AutoDock 4 default (~6%)
            tournament_size: 4,        // AutoDock 4 default
            local_search: LocalSearchKind::SolisWets(SolisWetsParams::lamarckian()),
        }
    }
}

impl GaParams {
    /// A small, fast configuration for quick previews and unit tests.
    pub fn fast() -> Self {
        GaParams {
            population: 12,
            generations: 6,
            mutation_rate: 0.15,
            crossover_rate: 0.8,
            local_search_fraction: 0.1,
            tournament_size: 3,
            // A short Solis-Wets budget keeps tests snappy.
            local_search: LocalSearchKind::SolisWets(SolisWetsParams::fast()),
        }
    }

    /// Validate the parameter ranges.
    fn validate(&self) -> Result<()> {
        if self.population < 2 {
            return Err(DockScreenError::invalid(
                "population",
                "GA population must be at least 2",
            ));
        }
        if self.tournament_size < 1 || self.tournament_size > self.population {
            return Err(DockScreenError::invalid(
                "tournament_size",
                "tournament size must be in 1..=population",
            ));
        }
        for (name, v) in [
            ("mutation_rate", self.mutation_rate),
            ("crossover_rate", self.crossover_rate),
            ("local_search_fraction", self.local_search_fraction),
        ] {
            if !(0.0..=1.0).contains(&v) {
                return Err(DockScreenError::invalid(name, "must be in 0..=1"));
            }
        }
        Ok(())
    }
}

/// The Lamarckian genetic-algorithm pose search.
pub struct LamarckianGa {
    params: GaParams,
}

impl LamarckianGa {
    /// Build a GA search with the given parameters.
    pub fn new(params: GaParams) -> Self {
        LamarckianGa { params }
    }

    /// Run the genetic algorithm and return the best `(Pose, score)`
    /// found, plus the score of every generation's best individual
    /// (the convergence history).
    ///
    /// Starting poses are sampled uniformly within the box defined by
    /// `center` ± `size/2`.
    pub fn run(
        &self,
        objective: &PoseObjective,
        center: Vector3<f64>,
        size: Vector3<f64>,
        seed: u64,
    ) -> Result<GaResult> {
        self.params.validate()?;
        let n_tor = objective.n_torsions();
        let mut rng = StdRng::seed_from_u64(seed);

        // --- initial population ------------------------------------
        let mut genomes: Vec<Vec<f64>> = (0..self.params.population)
            .map(|_| random_genome(&mut rng, center, size, n_tor))
            .collect();
        let mut fitness: Vec<f64> = genomes
            .iter()
            .map(|g| objective.score(&vec_to_pose(g, n_tor)))
            .collect();

        let mut best_idx = argmin(&fitness);
        let mut best_genome = genomes[best_idx].clone();
        let mut best_score = fitness[best_idx];
        let mut history = vec![best_score];

        // --- generations -------------------------------------------
        for _ in 0..self.params.generations {
            let mut next: Vec<Vec<f64>> = Vec::with_capacity(self.params.population);
            // Elitism: the current best always survives.
            next.push(best_genome.clone());
            while next.len() < self.params.population {
                let p1 = self.tournament(&fitness, &mut rng);
                let p2 = self.tournament(&fitness, &mut rng);
                let mut child = if rng.gen::<f64>() < self.params.crossover_rate {
                    crossover(&genomes[p1], &genomes[p2], &mut rng)
                } else {
                    // No crossover — copy the fitter parent.
                    if fitness[p1] <= fitness[p2] {
                        genomes[p1].clone()
                    } else {
                        genomes[p2].clone()
                    }
                };
                mutate(&mut child, self.params.mutation_rate, &mut rng);
                next.push(child);
            }
            genomes = next;
            fitness = genomes
                .iter()
                .map(|g| objective.score(&vec_to_pose(g, n_tor)))
                .collect();

            // --- Lamarckian local search ---------------------------
            let n_ls =
                (self.params.population as f64 * self.params.local_search_fraction).ceil() as usize;
            for k in 0..n_ls {
                let i = rng.gen_range(0..self.params.population);
                let pose = vec_to_pose(&genomes[i], n_tor);
                let ls_seed = seed
                    .wrapping_add(0xCAFE_F00D_DEAD_BEEF)
                    .wrapping_add((i as u64) << 16)
                    .wrapping_add(k as u64);
                let (refined, refined_score) =
                    local_search(objective, &pose, &self.params.local_search, ls_seed);
                // Write the refined genome back — the Lamarckian step.
                genomes[i] = pose_to_vec(&refined);
                fitness[i] = refined_score;
            }

            best_idx = argmin(&fitness);
            if fitness[best_idx] < best_score {
                best_score = fitness[best_idx];
                best_genome = genomes[best_idx].clone();
            }
            history.push(best_score);
        }

        Ok(GaResult {
            best_pose: vec_to_pose(&best_genome, n_tor),
            best_score,
            history,
        })
    }

    /// Tournament selection — pick `tournament_size` random
    /// individuals, return the index of the fittest.
    fn tournament(&self, fitness: &[f64], rng: &mut StdRng) -> usize {
        let mut best = rng.gen_range(0..fitness.len());
        for _ in 1..self.params.tournament_size {
            let c = rng.gen_range(0..fitness.len());
            if fitness[c] < fitness[best] {
                best = c;
            }
        }
        best
    }
}

/// The outcome of a Lamarckian-GA run.
#[derive(Clone, Debug)]
pub struct GaResult {
    /// The best pose found.
    pub best_pose: Pose,
    /// Score of [`GaResult::best_pose`].
    pub best_score: f64,
    /// Best-of-generation score history (index 0 is the initial
    /// population's best; one entry per generation thereafter).
    pub history: Vec<f64>,
}

/// A random genome: translation uniform in the box, a random
/// orientation, random torsions.
fn random_genome(
    rng: &mut StdRng,
    center: Vector3<f64>,
    size: Vector3<f64>,
    n_tor: usize,
) -> Vec<f64> {
    let half = size / 2.0;
    let mut g = Vec::with_capacity(6 + n_tor);
    g.push(center.x + rng.gen_range(-half.x..half.x));
    g.push(center.y + rng.gen_range(-half.y..half.y));
    g.push(center.z + rng.gen_range(-half.z..half.z));
    // Axis-angle rotation: a Rodrigues vector with magnitude in [-π,π].
    let axis = Vector3::new(
        rng.gen_range(-1.0..1.0),
        rng.gen_range(-1.0..1.0),
        rng.gen_range(-1.0..1.0),
    );
    let axis = if axis.norm() > 1e-9 {
        axis.normalize()
    } else {
        Vector3::x()
    };
    let angle = rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI);
    g.push(axis.x * angle);
    g.push(axis.y * angle);
    g.push(axis.z * angle);
    for _ in 0..n_tor {
        g.push(rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI));
    }
    g
}

/// Uniform crossover — each gene is taken from one parent or the other
/// with equal probability.
fn crossover(a: &[f64], b: &[f64], rng: &mut StdRng) -> Vec<f64> {
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| if rng.gen::<bool>() { x } else { y })
        .collect()
}

/// Mutate a genome — each gene, with probability `rate`, gets a
/// **Cauchy**-distributed kick (AutoDock 4's published default).
/// Translation / rotation / torsion genes use the per-DOF scale
/// AutoDock 4 uses for its mutation operator.
///
/// The Cauchy distribution is heavy-tailed, so most mutations are
/// small but rare ones can jump across the box — exactly the
/// exploration / exploitation balance the LGA needs.
fn mutate(genome: &mut [f64], rate: f64, rng: &mut StdRng) {
    for (i, gene) in genome.iter_mut().enumerate() {
        if rng.gen::<f64>() >= rate {
            continue;
        }
        // Per-DOF Cauchy scale (γ): AutoDock 4 ga.cc defaults.
        // Translation 1 Å, rotation 0.2 rad, torsion 0.2 rad.
        let gamma = if i < 3 { 1.0_f64 } else { 0.2_f64 };
        *gene += cauchy_sample(rng, gamma);
    }
}

/// Sample from a centred Cauchy distribution of scale `gamma`.
/// Inverse-CDF method: `γ · tan(π (u - 0.5))` with `u ∈ (0,1)`. The
/// `clamp` keeps the rare ultra-long tail from blowing the search out
/// of the box — AutoDock 4 does the same (`min(±10ρ, sample)`).
fn cauchy_sample(rng: &mut StdRng, gamma: f64) -> f64 {
    let u: f64 = rng.gen_range(0.001..0.999);
    let s = gamma * (std::f64::consts::PI * (u - 0.5)).tan();
    s.clamp(-10.0 * gamma, 10.0 * gamma)
}

/// Local search for the Lamarckian step.
///
/// Dispatches on [`LocalSearchKind`]: the AutoDock 4 default is
/// **Solis-Wets** (real, with the published bias / expansion /
/// contraction schedule from [`crate::search::solis_wets`]). The
/// fallback is a self-contained coordinate descent, useful when the
/// objective is non-smooth in a way Solis-Wets's continuous step
/// adaptation cannot navigate.
///
/// `valenx-dock`'s BFGS minimiser ([`valenx_dock::search::bfgs::minimize_bfgs`])
/// only optimises against the dock crate's own `GridBundle` type —
/// this crate's [`PoseObjective`] is backed by a `valenx-dock-screen`
/// affinity-map set (which can be an AutoDock 4 grid the dock crate
/// does not understand). The Vina-class ILS driver still reuses the
/// dock crate's BFGS verbatim — see [`crate::search::ils`].
fn local_search(
    objective: &PoseObjective,
    pose: &Pose,
    kind: &LocalSearchKind,
    seed: u64,
) -> (Pose, f64) {
    match kind {
        LocalSearchKind::SolisWets(p) => solis_wets(objective, pose, p, seed),
        LocalSearchKind::CoordinateDescent => coordinate_descent(objective, pose, 25),
    }
}

/// A compact local optimiser: probe each genome coordinate ± a
/// shrinking step, keep any move that lowers the score. Cheap, robust,
/// and works directly on this crate's [`PoseObjective`].
pub(crate) fn coordinate_descent(
    objective: &PoseObjective,
    start: &Pose,
    max_iter: usize,
) -> (Pose, f64) {
    let n_tor = objective.n_torsions();
    let mut x = pose_to_vec(start);
    let mut best = objective.score(start);
    let n = x.len();
    // Per-coordinate initial step: a larger one for the 3 translation
    // genes, a smaller one for the 3 rotation genes and the torsions.
    let mut step = vec![0.0; n];
    for (i, s) in step.iter_mut().enumerate() {
        *s = if i < 3 { 0.5 } else { 0.2 };
    }
    for _ in 0..max_iter {
        let mut improved = false;
        for i in 0..n {
            for &dir in &[1.0_f64, -1.0] {
                let saved = x[i];
                x[i] = saved + dir * step[i];
                let s = objective.score(&vec_to_pose(&x, n_tor));
                if s < best {
                    best = s;
                    improved = true;
                } else {
                    x[i] = saved;
                }
            }
        }
        if !improved {
            // Shrink steps; stop once they are negligible.
            for s in step.iter_mut() {
                *s *= 0.5;
            }
            if step.iter().all(|s| *s < 1e-3) {
                break;
            }
        }
    }
    (vec_to_pose(&x, n_tor), best)
}

/// Index of the smallest element. `0` for an empty slice (callers
/// guarantee non-empty).
fn argmin(v: &[f64]) -> usize {
    let mut best = 0;
    for (i, &x) in v.iter().enumerate() {
        if x < v[best] {
            best = i;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prep::gridbox::GridBox;
    use crate::score::gridmap::{AffinityMapSet, MapKind};
    use valenx_dock::atom_type::Ad4AtomType;
    use valenx_dock::ligand::Ligand;
    use valenx_dock::receptor::{Receptor, ReceptorAtom};

    fn setup() -> (Ligand, AffinityMapSet) {
        let receptor = Receptor {
            atoms: vec![ReceptorAtom {
                position: Vector3::zeros(),
                ad4_type: Ad4AtomType::C,
                partial_charge: 0.0,
            }],
        };
        let pdbqt = "ROOT
ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ENDROOT
TORSDOF 0
";
        let lig = Ligand::from_pdbqt(pdbqt).unwrap();
        let grid = GridBox::with_spacing([0.0; 3], [16.0; 3], 0.375).unwrap();
        let maps =
            AffinityMapSet::precompute(&receptor, &[Ad4AtomType::C], &grid, MapKind::Vina).unwrap();
        (lig, maps)
    }

    #[test]
    fn ga_rejects_bad_params() {
        let bad = GaParams {
            population: 1,
            ..GaParams::fast()
        };
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let res = LamarckianGa::new(bad).run(&obj, Vector3::zeros(), Vector3::new(8.0, 8.0, 8.0), 1);
        assert!(res.is_err());
    }

    #[test]
    fn ga_finds_a_favourable_pose() {
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let ga = LamarckianGa::new(GaParams::fast());
        let result = ga
            .run(&obj, Vector3::zeros(), Vector3::new(10.0, 10.0, 10.0), 42)
            .unwrap();
        // A single C-C pair: the global minimum is clearly negative.
        assert!(
            result.best_score < -0.005,
            "GA best score = {}",
            result.best_score
        );
    }

    #[test]
    fn ga_history_is_monotone_non_increasing() {
        // Best-of-generation can only stay equal or improve (elitism).
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let result = LamarckianGa::new(GaParams::fast())
            .run(&obj, Vector3::zeros(), Vector3::new(10.0, 10.0, 10.0), 7)
            .unwrap();
        for w in result.history.windows(2) {
            assert!(
                w[1] <= w[0] + 1e-9,
                "history regressed: {:?}",
                result.history
            );
        }
    }

    #[test]
    fn ga_is_deterministic_for_a_fixed_seed() {
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let ga = LamarckianGa::new(GaParams::fast());
        let a = ga
            .run(&obj, Vector3::zeros(), Vector3::new(10.0, 10.0, 10.0), 99)
            .unwrap();
        let b = ga
            .run(&obj, Vector3::zeros(), Vector3::new(10.0, 10.0, 10.0), 99)
            .unwrap();
        assert!((a.best_score - b.best_score).abs() < 1e-12);
    }

    #[test]
    fn coordinate_descent_improves_a_poor_start() {
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let mut start = Pose::identity(0);
        start.translation = Vector3::new(6.0, 0.0, 0.0);
        let before = obj.score(&start);
        let (_p, after) = coordinate_descent(&obj, &start, 30);
        assert!(after <= before, "local search must not worsen the score");
    }

    #[test]
    fn crossover_child_genes_come_from_parents() {
        let mut rng = StdRng::seed_from_u64(1);
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let child = crossover(&a, &b, &mut rng);
        for (i, &g) in child.iter().enumerate() {
            assert!(g == a[i] || g == b[i], "gene {i} not from a parent");
        }
    }

    #[test]
    fn cauchy_sample_has_zero_median_and_heavy_tails() {
        // Sanity check: out of many Cauchy samples roughly half should
        // be negative and we should see at least one large excursion.
        let mut rng = StdRng::seed_from_u64(33);
        let mut neg = 0;
        let mut max_abs: f64 = 0.0;
        let n = 1_000;
        for _ in 0..n {
            let s = cauchy_sample(&mut rng, 1.0);
            if s < 0.0 {
                neg += 1;
            }
            if s.abs() > max_abs {
                max_abs = s.abs();
            }
        }
        // Roughly symmetric — within 10% of half.
        assert!(
            (neg as f64 / n as f64 - 0.5).abs() < 0.1,
            "Cauchy negatives = {neg}/{n}, expected ~50%"
        );
        // Heavy tail: at least one sample over the median scale.
        assert!(
            max_abs > 2.0,
            "Cauchy should show heavy tails, max |s| = {max_abs}"
        );
    }

    #[test]
    fn solis_wets_local_search_kind_is_used_by_default() {
        // A regression: a GA built with `GaParams::fast()` (the test
        // default) uses the Solis-Wets local search.
        match GaParams::fast().local_search {
            LocalSearchKind::SolisWets(_) => {}
            other => panic!("expected SolisWets by default, got {other:?}"),
        }
    }

    #[test]
    fn ga_with_coordinate_descent_local_search_still_runs() {
        let (lig, maps) = setup();
        let obj = PoseObjective::uncharged(&lig, &maps);
        let params = GaParams {
            local_search: LocalSearchKind::CoordinateDescent,
            ..GaParams::fast()
        };
        let result = LamarckianGa::new(params)
            .run(&obj, Vector3::zeros(), Vector3::new(10.0, 10.0, 10.0), 13)
            .unwrap();
        assert!(result.best_score < 0.0);
    }

    #[test]
    fn ga_converges_on_a_quadratic_test_surface() {
        // Carve a deep quadratic basin into the grid map and verify the
        // LGA finds it. This is a direct convergence test on a known
        // surface — the validation the output spec asks for.
        let (lig, mut maps) = setup();
        let well;
        {
            let m = &mut maps.maps[0].1;
            let (nx, ny, nz) = m.dims;
            let cx = (nx / 2) as f64;
            let cy = (ny / 2) as f64;
            let cz = (nz / 2) as f64;
            for iz in 0..nz {
                for iy in 0..ny {
                    for ix in 0..nx {
                        let dx = ix as f64 - cx;
                        let dy = iy as f64 - cy;
                        let dz = iz as f64 - cz;
                        let idx = m.index(ix, iy, iz);
                        m.data[idx] = 0.02 * (dx * dx + dy * dy + dz * dz) - 6.0;
                    }
                }
            }
            well = m.origin + Vector3::new(cx, cy, cz) * m.spacing;
        }
        let obj = PoseObjective::uncharged(&lig, &maps);
        let res = LamarckianGa::new(GaParams::fast())
            .run(&obj, Vector3::zeros(), Vector3::new(10.0, 10.0, 10.0), 31)
            .unwrap();
        // Should land near the floor of the basin.
        assert!(res.best_score < -3.0, "got {}", res.best_score);
        let d = (res.best_pose.translation - well).norm();
        assert!(d < 3.0, "expected pose near well, got d = {d}");
    }
}
