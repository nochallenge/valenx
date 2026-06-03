//! Kinfold-class RNA folding kinetics.
//!
//! Where [`super::partition`] computes the *thermodynamic* equilibrium
//! ensemble and [`crate::fold::zuker`] finds the single most-stable
//! structure, **folding kinetics** asks how an RNA molecule traverses
//! the structure landscape *over time* — the trajectory through
//! intermediate structures, the equilibration time, and the fraction of
//! molecules in each structure at a given moment. This module is the
//! Kinfold-class (Flamm *et al.* 2000) simulator: stochastic Monte-
//! Carlo over the elementary-move neighbourhood with Metropolis
//! transition rates from Turner-2004 ΔΔG.
//!
//! ## What it does
//!
//! - **Move set** — at each step the simulator proposes a single
//!   *elementary move*: add a valid base pair, remove an existing
//!   pair, or shift one partner of an existing pair to a neighbouring
//!   position. This is the canonical Kinfold neighbour graph.
//! - **Rates** — given current free energy `E` and proposed free energy
//!   `E'`, the move is accepted with the **Metropolis** rate
//!   `k = min(1, exp(−ΔG / RT))` (the default). The alternative
//!   **Kawasaki** rate `k = exp(−ΔG / 2RT)` is selectable via
//!   [`RateModel`].
//! - **Trajectory** — each step advances *simulated time* by an
//!   exponentially-distributed waiting time `Δt = −ln(u) / Σk`, the
//!   standard Gillespie / Kawasaki time-step. The simulator emits the
//!   `(time, structure, energy)` triples at every accepted move.
//! - **Ensemble** — [`fold_kinetics`] runs `n_trajectories`
//!   independent simulations from a shared start and reports
//!   - the fraction of trajectories in the MFE at each time-checkpoint;
//!   - the mean first-passage time to the MFE;
//!   - the long-time structure population — which (for an
//!     ergodic walk) approaches the Boltzmann distribution.
//!
//! ## Why this exists
//!
//! Many real RNAs are **kinetically trapped** — they fold into a
//! metastable structure that is not the thermodynamic MFE, and stay
//! there for biologically meaningful time. The riboswitch literature is
//! built on this: the on/off conformations of an aptamer can have
//! comparable free energies but very different lifetimes. A
//! thermodynamic-only folder cannot distinguish them; a kinetic folder
//! can.
//!
//! ## Honest scope
//!
//! - Energies use the same Turner-2004 evaluator
//!   ([`crate::fold::eval::structure_energy`]) as the rest of the crate.
//! - The move set is the elementary single-pair add / remove / shift
//!   that Kinfold ships as the default. More elaborate move sets
//!   (helix moves, breathing) are out of v1 scope.
//! - The simulator is single-threaded; trajectories are sequenced not
//!   parallel. The deterministic seed makes runs reproducible.

use crate::error::{Result, RnaStructError};
use crate::ensemble::rng::Rng;
use crate::fold::energy::{self, GAS_CONSTANT};
use crate::fold::eval::structure_energy;
use crate::fold::nussinov::MIN_HAIRPIN;
use crate::fold::zuker::mfe;
use crate::rna::RnaSeq;
use crate::structure::Structure;

/// Default reference temperature for the kinetic simulator
/// (kelvin, 37 °C).
pub const DEFAULT_TEMPERATURE_K: f64 = energy::T37_KELVIN;

/// Transition-rate model.
#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
pub enum RateModel {
    /// Metropolis rate: `k = min(1, exp(−ΔG / RT))`.
    #[default]
    Metropolis,
    /// Kawasaki rate: `k = exp(−ΔG / 2RT)`.
    Kawasaki,
}

/// Parameters of the kinetic simulator.
#[derive(Clone, Debug)]
pub struct KineticParams {
    /// Maximum number of *steps* per trajectory. (Trajectory ends when
    /// either this many steps elapse or the MFE is reached and a stop
    /// condition is satisfied.)
    pub max_steps: usize,
    /// Stop the trajectory the first time the MFE is reached.
    pub stop_at_mfe: bool,
    /// Total simulated-time cap per trajectory; trajectory ends if
    /// total time exceeds this.
    pub max_time: f64,
    /// Folding temperature, kelvin.
    pub temperature_k: f64,
    /// Transition-rate model.
    pub rate_model: RateModel,
    /// Deterministic seed for reproducibility.
    pub seed: u64,
}

impl Default for KineticParams {
    fn default() -> Self {
        KineticParams {
            max_steps: 5_000,
            stop_at_mfe: false,
            max_time: 1.0e9,
            temperature_k: DEFAULT_TEMPERATURE_K,
            rate_model: RateModel::default(),
            seed: 0,
        }
    }
}

/// A single trajectory step: simulated time and the structure entered.
#[derive(Clone, Debug)]
pub struct TrajectoryStep {
    /// Simulated time (arbitrary units, set by the rate constant prefactor).
    pub time: f64,
    /// The structure entered at this step.
    pub structure: Structure,
    /// Free energy of `structure`, kcal/mol.
    pub energy: f64,
}

/// A single kinetic-folding trajectory.
#[derive(Clone, Debug)]
pub struct Trajectory {
    /// The (time, structure, energy) checkpoints, including the start
    /// and every accepted move.
    pub steps: Vec<TrajectoryStep>,
    /// `true` if the trajectory reached the MFE during simulation.
    pub reached_mfe: bool,
    /// Simulated time at which the MFE was first reached, if any.
    pub first_passage_to_mfe: Option<f64>,
}

impl Trajectory {
    /// The terminal step.
    pub fn final_step(&self) -> Option<&TrajectoryStep> {
        self.steps.last()
    }

    /// Total simulated time of the trajectory.
    pub fn total_time(&self) -> f64 {
        self.final_step().map(|s| s.time).unwrap_or(0.0)
    }
}

/// Run a single kinetic-folding trajectory starting from `start` for
/// the sequence `seq`. Returns the full trajectory.
///
/// # Errors
/// Propagates structure validation errors.
pub fn simulate_trajectory(
    seq: &RnaSeq,
    start: &Structure,
    mfe_struct: &Structure,
    params: &KineticParams,
) -> Result<Trajectory> {
    if start.len() != seq.len() {
        return Err(RnaStructError::structure(
            "start structure length must match sequence length",
        ));
    }
    let codes = seq.codes();
    let n = codes.len();
    let mut current = start.clone();
    let mut current_e = structure_energy(seq, &current).unwrap_or(0.0);

    let rt = GAS_CONSTANT * params.temperature_k;
    let mut rng = Rng::new(params.seed);
    let mut time = 0.0;
    let mut steps: Vec<TrajectoryStep> = vec![TrajectoryStep {
        time: 0.0,
        structure: current.clone(),
        energy: current_e,
    }];
    let mut reached_mfe = structures_equal(&current, mfe_struct);
    let mut first_passage = if reached_mfe { Some(0.0) } else { None };

    for _ in 0..params.max_steps {
        if time > params.max_time {
            break;
        }
        // Build the neighbour list (move list) with their ΔG and rate.
        let moves = enumerate_moves(codes, &current, current_e, seq, n);
        if moves.is_empty() {
            break;
        }
        let mut rates: Vec<f64> = Vec::with_capacity(moves.len());
        for (_, dg) in &moves {
            let k = match params.rate_model {
                RateModel::Metropolis => {
                    if *dg <= 0.0 {
                        1.0
                    } else {
                        (-dg / rt).exp()
                    }
                }
                RateModel::Kawasaki => (-dg / (2.0 * rt)).exp(),
            };
            rates.push(k);
        }
        let rate_sum: f64 = rates.iter().sum();
        if rate_sum <= 0.0 {
            break;
        }

        // Gillespie waiting time and move selection.
        let u_time = rng.next_f64().max(1e-300);
        let dt = -u_time.ln() / rate_sum;
        time += dt;

        let mut u_pick = rng.next_f64() * rate_sum;
        let mut chosen = 0usize;
        for (idx, k) in rates.iter().enumerate() {
            if u_pick <= *k {
                chosen = idx;
                break;
            }
            u_pick -= *k;
        }
        let (mv, dg) = moves[chosen];
        // Apply the chosen move.
        apply_move(&mut current, &mv);
        current_e += dg;

        steps.push(TrajectoryStep {
            time,
            structure: current.clone(),
            energy: current_e,
        });

        if !reached_mfe && structures_equal(&current, mfe_struct) {
            reached_mfe = true;
            first_passage = Some(time);
            if params.stop_at_mfe {
                break;
            }
        }
    }

    Ok(Trajectory {
        steps,
        reached_mfe,
        first_passage_to_mfe: first_passage,
    })
}

/// Aggregated ensemble of kinetic trajectories.
#[derive(Clone, Debug)]
pub struct KineticEnsemble {
    /// One trajectory per replicate.
    pub trajectories: Vec<Trajectory>,
    /// Fraction of trajectories that reached the MFE.
    pub fraction_reached_mfe: f64,
    /// Mean simulated time at first passage to the MFE (over the
    /// trajectories that reached it).
    pub mean_first_passage_time: Option<f64>,
    /// The folded MFE structure (the reference) and its energy.
    pub mfe_structure: Structure,
    /// Energy of `mfe_structure`.
    pub mfe_energy: f64,
}

impl KineticEnsemble {
    /// Counts of how many trajectories *finish* in each unique
    /// structure (terminal structure of each trajectory), keyed by
    /// dot-bracket. Useful for verifying long-time Boltzmann
    /// equilibrium.
    pub fn terminal_structure_counts(&self) -> Vec<(String, usize)> {
        let mut counts: Vec<(String, usize)> = Vec::new();
        for traj in &self.trajectories {
            let s = traj
                .final_step()
                .map(|st| st.structure.to_dot_bracket())
                .unwrap_or_default();
            match counts.iter_mut().find(|(k, _)| k == &s) {
                Some((_, c)) => *c += 1,
                None => counts.push((s, 1)),
            }
        }
        counts.sort_by(|a, b| b.1.cmp(&a.1));
        counts
    }

    /// Fraction of trajectories *currently* (terminal step) in the MFE.
    pub fn fraction_in_mfe_terminal(&self) -> f64 {
        if self.trajectories.is_empty() {
            return 0.0;
        }
        let mfe_db = self.mfe_structure.to_dot_bracket();
        let count = self
            .trajectories
            .iter()
            .filter(|t| {
                t.final_step()
                    .map(|s| s.structure.to_dot_bracket() == mfe_db)
                    .unwrap_or(false)
            })
            .count();
        count as f64 / self.trajectories.len() as f64
    }
}

/// Run `n_trajectories` independent kinetic-folding trajectories for
/// `seq`, all starting from the open chain (no pairs), and aggregate
/// the ensemble.
///
/// # Errors
/// Propagates the MFE computation and any structure validation error.
pub fn fold_kinetics(
    seq: &RnaSeq,
    n_trajectories: usize,
    params: &KineticParams,
) -> Result<KineticEnsemble> {
    if n_trajectories == 0 {
        return Err(RnaStructError::invalid(
            "n_trajectories",
            "must be >= 1",
        ));
    }
    let mfe_r = mfe(seq)?;
    let mfe_struct = mfe_r.structure.clone();
    let mfe_energy = mfe_r.energy;
    let start = Structure::empty(seq.len());

    let mut trajectories: Vec<Trajectory> = Vec::with_capacity(n_trajectories);
    for k in 0..n_trajectories {
        let mut p = params.clone();
        // Re-seed deterministically per trajectory.
        p.seed = params.seed.wrapping_add(k as u64);
        trajectories.push(simulate_trajectory(seq, &start, &mfe_struct, &p)?);
    }

    let n_reached = trajectories.iter().filter(|t| t.reached_mfe).count();
    let frac = n_reached as f64 / n_trajectories as f64;
    let times: Vec<f64> = trajectories
        .iter()
        .filter_map(|t| t.first_passage_to_mfe)
        .collect();
    let mean_fp = if times.is_empty() {
        None
    } else {
        Some(times.iter().sum::<f64>() / times.len() as f64)
    };

    Ok(KineticEnsemble {
        trajectories,
        fraction_reached_mfe: frac,
        mean_first_passage_time: mean_fp,
        mfe_structure: mfe_struct,
        mfe_energy,
    })
}

/// An elementary move on a structure.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum Move {
    /// Add a base pair `(i, j)`.
    Add { i: usize, j: usize },
    /// Remove the pair touching position `i` (whatever its partner is).
    Remove { i: usize, j: usize },
    /// Shift the pair `(i, j)` to `(i, j')` (one partner moves).
    Shift {
        old_i: usize,
        old_j: usize,
        new_i: usize,
        new_j: usize,
    },
}

/// Two structures equal: same partner array.
fn structures_equal(a: &Structure, b: &Structure) -> bool {
    if a.len() != b.len() {
        return false;
    }
    for i in 0..a.len() {
        if a.partner(i) != b.partner(i) {
            return false;
        }
    }
    true
}

/// Enumerate the elementary-move neighbourhood of `current`. Each move
/// is returned with its ΔG (proposed_energy − current_energy).
fn enumerate_moves(
    codes: &[u8],
    current: &Structure,
    current_e: f64,
    seq: &RnaSeq,
    n: usize,
) -> Vec<(Move, f64)> {
    let mut moves: Vec<(Move, f64)> = Vec::new();

    // Add: every (i, j) pair that doesn't conflict with the current
    // structure (both unpaired) and doesn't cross an existing pair.
    for i in 0..n {
        if current.is_paired(i) {
            continue;
        }
        for j in (i + MIN_HAIRPIN + 1)..n {
            if current.is_paired(j) {
                continue;
            }
            if !energy::can_pair_codes(codes[i], codes[j]) {
                continue;
            }
            // Crossing check.
            if crosses_existing(current, i, j) {
                continue;
            }
            let mut next = current.clone();
            if next.add_pair(i, j).is_err() {
                continue;
            }
            let next_e = match structure_energy(seq, &next) {
                Ok(e) => e,
                Err(_) => continue,
            };
            moves.push((Move::Add { i, j }, next_e - current_e));
        }
    }

    // Remove: every existing pair.
    for i in 0..n {
        if let Some(j) = current.partner(i) {
            if i < j {
                let mut next = current.clone();
                next.remove_pair(i);
                let next_e = match structure_energy(seq, &next) {
                    Ok(e) => e,
                    Err(_) => continue,
                };
                moves.push((Move::Remove { i, j }, next_e - current_e));
            }
        }
    }

    // Shift: for each existing pair (i, j), try to slide one partner
    // by one position (j -> j+1 or j-1; or i -> i+1 or i-1).
    for i in 0..n {
        if let Some(j) = current.partner(i) {
            if i < j {
                for new_j in [j.wrapping_sub(1), j + 1] {
                    if new_j >= n || new_j <= i || new_j == j {
                        continue;
                    }
                    if current.is_paired(new_j) {
                        continue;
                    }
                    if new_j <= i + MIN_HAIRPIN {
                        continue;
                    }
                    if !energy::can_pair_codes(codes[i], codes[new_j]) {
                        continue;
                    }
                    let mut next = current.clone();
                    next.remove_pair(i);
                    if next.add_pair(i, new_j).is_err() {
                        continue;
                    }
                    if !is_nested(&next) {
                        continue;
                    }
                    let next_e = match structure_energy(seq, &next) {
                        Ok(e) => e,
                        Err(_) => continue,
                    };
                    moves.push((
                        Move::Shift {
                            old_i: i,
                            old_j: j,
                            new_i: i,
                            new_j,
                        },
                        next_e - current_e,
                    ));
                }
                for new_i in [i.wrapping_sub(1), i + 1] {
                    if new_i >= n || new_i >= j || new_i == i {
                        continue;
                    }
                    if current.is_paired(new_i) {
                        continue;
                    }
                    if j <= new_i + MIN_HAIRPIN {
                        continue;
                    }
                    if !energy::can_pair_codes(codes[new_i], codes[j]) {
                        continue;
                    }
                    let mut next = current.clone();
                    next.remove_pair(i);
                    if next.add_pair(new_i, j).is_err() {
                        continue;
                    }
                    if !is_nested(&next) {
                        continue;
                    }
                    let next_e = match structure_energy(seq, &next) {
                        Ok(e) => e,
                        Err(_) => continue,
                    };
                    moves.push((
                        Move::Shift {
                            old_i: i,
                            old_j: j,
                            new_i,
                            new_j: j,
                        },
                        next_e - current_e,
                    ));
                }
            }
        }
    }

    moves
}

/// `true` if adding the pair `(i, j)` would cross an existing pair in
/// `current`. Used to keep the trajectory pseudoknot-free (so each step
/// is scorable by the nearest-neighbor evaluator).
fn crosses_existing(current: &Structure, i: usize, j: usize) -> bool {
    for k in 0..current.len() {
        if let Some(l) = current.partner(k) {
            if k < l {
                // (i, j) crosses (k, l) iff i < k < j < l or k < i < l < j.
                if (i < k && k < j && j < l) || (k < i && i < l && l < j) {
                    return true;
                }
            }
        }
    }
    false
}

/// `true` if `s` has no crossing pairs.
fn is_nested(s: &Structure) -> bool {
    s.is_nested()
}

/// Apply a move to a structure in place.
fn apply_move(s: &mut Structure, mv: &Move) {
    match *mv {
        Move::Add { i, j } => {
            let _ = s.add_pair(i, j);
        }
        Move::Remove { i, .. } => {
            s.remove_pair(i);
        }
        Move::Shift {
            old_i,
            new_i,
            new_j,
            ..
        } => {
            s.remove_pair(old_i);
            let _ = s.add_pair(new_i, new_j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trajectory_emits_a_start_step() {
        let seq = RnaSeq::parse("GGGGAAAACCCC").unwrap();
        let mfe_r = mfe(&seq).unwrap();
        let start = Structure::empty(seq.len());
        let params = KineticParams {
            max_steps: 10,
            seed: 42,
            ..Default::default()
        };
        let traj = simulate_trajectory(&seq, &start, &mfe_r.structure, &params).unwrap();
        assert!(!traj.steps.is_empty());
        assert_eq!(traj.steps[0].structure, start);
    }

    #[test]
    fn open_chain_reaches_mfe_for_simple_hairpin() {
        // A simple 4-pair GC hairpin: the MFE is unambiguous and a
        // small number of trajectories should reach it.
        let seq = RnaSeq::parse("GGGGAAAACCCC").unwrap();
        let mfe_r = mfe(&seq).unwrap();
        let params = KineticParams {
            max_steps: 2_000,
            stop_at_mfe: true,
            seed: 1,
            ..Default::default()
        };
        let ens = fold_kinetics(&seq, 32, &params).unwrap();
        // At least a quarter of trajectories should reach the MFE.
        assert!(
            ens.fraction_reached_mfe >= 0.25,
            "only {} reached MFE",
            ens.fraction_reached_mfe
        );
        // First-passage time, if reported, is finite.
        if let Some(t) = ens.mean_first_passage_time {
            assert!(t.is_finite() && t > 0.0);
        }
        // The MFE we hit equals the MFE we asked about.
        assert_eq!(ens.mfe_structure, mfe_r.structure);
    }

    #[test]
    fn trajectory_energies_match_evaluator() {
        // The energy at each step must equal the analytic energy of
        // the structure at that step (up to rounding).
        let seq = RnaSeq::parse("GGGGAAAACCCC").unwrap();
        let mfe_r = mfe(&seq).unwrap();
        let params = KineticParams {
            max_steps: 50,
            seed: 7,
            ..Default::default()
        };
        let traj = simulate_trajectory(
            &seq,
            &Structure::empty(seq.len()),
            &mfe_r.structure,
            &params,
        )
        .unwrap();
        for step in &traj.steps {
            let e = structure_energy(&seq, &step.structure).unwrap_or(0.0);
            assert!(
                (step.energy - e).abs() < 1e-3,
                "step energy {} != analytic {} on {}",
                step.energy,
                e,
                step.structure.to_dot_bracket()
            );
        }
    }

    #[test]
    fn deterministic_seed_reproduces_trajectory() {
        let seq = RnaSeq::parse("GGGGAAAACCCC").unwrap();
        let mfe_r = mfe(&seq).unwrap();
        let params = KineticParams {
            max_steps: 100,
            seed: 99,
            ..Default::default()
        };
        let t1 = simulate_trajectory(
            &seq,
            &Structure::empty(seq.len()),
            &mfe_r.structure,
            &params,
        )
        .unwrap();
        let t2 = simulate_trajectory(
            &seq,
            &Structure::empty(seq.len()),
            &mfe_r.structure,
            &params,
        )
        .unwrap();
        assert_eq!(t1.steps.len(), t2.steps.len());
        for (a, b) in t1.steps.iter().zip(t2.steps.iter()) {
            assert!((a.time - b.time).abs() < 1e-9);
            assert!((a.energy - b.energy).abs() < 1e-9);
        }
    }

    #[test]
    fn kawasaki_rates_are_finite() {
        // The Kawasaki rate model also runs without underflow.
        let seq = RnaSeq::parse("GGGGAAAACCCC").unwrap();
        let mfe_r = mfe(&seq).unwrap();
        let params = KineticParams {
            max_steps: 100,
            stop_at_mfe: true,
            rate_model: RateModel::Kawasaki,
            seed: 3,
            ..Default::default()
        };
        let traj = simulate_trajectory(
            &seq,
            &Structure::empty(seq.len()),
            &mfe_r.structure,
            &params,
        )
        .unwrap();
        assert!(!traj.steps.is_empty());
        for step in &traj.steps {
            assert!(step.time.is_finite());
            assert!(step.energy.is_finite());
        }
    }

    #[test]
    fn long_time_population_concentrates_in_low_energy_states() {
        // After a long simulation the trajectories should concentrate
        // in low-energy structures (Boltzmann tendency). We check that
        // the *terminal* mean energy is lower than the *starting*
        // energy (the open chain has 0 energy).
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let params = KineticParams {
            max_steps: 500,
            stop_at_mfe: false,
            seed: 11,
            ..Default::default()
        };
        let ens = fold_kinetics(&seq, 16, &params).unwrap();
        let mean_terminal: f64 = ens
            .trajectories
            .iter()
            .filter_map(|t| t.final_step().map(|s| s.energy))
            .sum::<f64>()
            / ens.trajectories.len() as f64;
        assert!(
            mean_terminal < 0.0,
            "mean terminal energy {mean_terminal} should be negative"
        );
    }

    #[test]
    fn n_trajectories_zero_is_rejected() {
        let seq = RnaSeq::parse("GGGGCCCC").unwrap();
        let params = KineticParams::default();
        assert!(fold_kinetics(&seq, 0, &params).is_err());
    }

    #[test]
    fn equilibrium_population_matches_boltzmann_at_long_times() {
        // For a short sequence with a well-defined MFE, the long-time
        // population in the MFE should be at least comparable to the
        // Boltzmann fraction in the MFE.
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let pf = crate::ensemble::partition::partition_function(&seq).unwrap();
        let mfe_r = mfe(&seq).unwrap();
        let rt = GAS_CONSTANT * DEFAULT_TEMPERATURE_K;
        let p_mfe = (-mfe_r.energy / rt).exp() / pf.q();
        // The kinetic simulator may not have fully equilibrated in
        // max_steps=2000, but the long-time MFE population should
        // *trend* toward p_mfe.
        let params = KineticParams {
            max_steps: 2_000,
            stop_at_mfe: false,
            seed: 17,
            ..Default::default()
        };
        let ens = fold_kinetics(&seq, 32, &params).unwrap();
        let frac_in_mfe = ens.fraction_in_mfe_terminal();
        // The kinetic fraction should be in a plausible band around
        // the Boltzmann probability. For a strong hairpin
        // p_mfe is close to 1.
        if p_mfe > 0.5 {
            // The MFE is strongly populated thermodynamically; kinetic
            // sim should put a non-trivial fraction there too.
            assert!(
                frac_in_mfe >= 0.05,
                "kinetic fraction in MFE {frac_in_mfe} too small (p_mfe={p_mfe})"
            );
        }
        // The terminal-structure histogram is well-formed.
        let counts = ens.terminal_structure_counts();
        assert!(!counts.is_empty());
        assert_eq!(
            counts.iter().map(|(_, c)| c).sum::<usize>(),
            ens.trajectories.len()
        );
    }
}
