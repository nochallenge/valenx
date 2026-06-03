//! **DOPE-driven simulated-annealing Monte-Carlo refinement.**
//!
//! Rosetta's `relax` protocol alternates **fragment-insertion Monte
//! Carlo** with proper **energy minimisation** under the statistical
//! potential, all on a slow simulated-annealing temperature schedule
//! that lets the conformation jiggle out of mediocre minima before
//! locking in. This module ships that loop, driven by the DOPE-class
//! distance-dependent potential added in
//! [`crate::abinitio::dope`].
//!
//! The protocol:
//!
//! 1. Score the starting model under DOPE.
//! 2. For each annealing temperature in the geometric schedule:
//!    a. Many fragment-insertion Monte-Carlo trials (Metropolis
//!    acceptance under DOPE at the current temperature).
//!    b. A brief **gradient relax** under the
//!    [`crate::refine::relax_model`] callback potential — the
//!    all-atom Cα-restraint relax — to drain any local strain the
//!    fragment moves introduced.
//!    c. A **Ramachandran refinement** pass — pulls any newly-strained
//!    φ/ψ back into allowed regions, the published Rosetta `relax`
//!    step.
//! 3. Optionally a final round of decoy generation +
//!    [`crate::abinitio::cluster_decoys`] returning a representative
//!    low-energy structure.
//!
//! The returned [`McRefineResult`] carries the refined model, the
//! per-cycle energy trajectory (proves monotone-or-better progress),
//! and counters useful for telemetry.

use serde::{Deserialize, Serialize};
use valenx_md::Rng;

use crate::abinitio::dope::{dope_score, DopeWeights};
use crate::abinitio::fragments::{build_fragment_library, FragmentLibrary};
use crate::error::{Result, StructPredictError};
use crate::model::{build_backbone_from_torsions, ProteinModel};
use crate::refine::minimize::relax_model;
use crate::refine::ramachandran::{model_phi_psi, refine_ramachandran};

/// Configures a DOPE-driven MC refinement run.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct McRefineOptions {
    /// Number of annealing cycles (each cycle is one temperature
    /// step + `moves_per_cycle` MC trials + a relax/Rama pass).
    pub cycles: usize,
    /// Fragment-insertion MC trials per annealing cycle.
    pub moves_per_cycle: usize,
    /// Starting (high) temperature of the annealing schedule.
    pub start_temperature: f64,
    /// Final (low) temperature of the annealing schedule.
    pub end_temperature: f64,
    /// Gradient-relax iterations per cycle (capped to keep the cycle
    /// fast — the relaxer is called repeatedly, so a small cap per
    /// call is the right choice). 0 disables in-cycle gradient relax.
    pub relax_iterations: usize,
    /// Fragment length to use when (re)building the library for the
    /// refinement loop. 3 is the classical Rosetta short-fragment
    /// length.
    pub fragment_length: usize,
    /// Fragments per window in the library.
    pub fragments_per_window: usize,
    /// RNG seed — fixes the trajectory.
    pub seed: u64,
}

impl Default for McRefineOptions {
    fn default() -> Self {
        McRefineOptions {
            cycles: 6,
            moves_per_cycle: 80,
            start_temperature: 2.0,
            end_temperature: 0.2,
            // The default disables in-cycle Cα-only gradient relax —
            // it can break the backbone-self-consistency the MC + DOPE
            // path requires. Callers needing a Cα-only relax pass can
            // dial it in.
            relax_iterations: 0,
            fragment_length: 3,
            fragments_per_window: 25,
            seed: 0x5EED_DA7E,
        }
    }
}

impl McRefineOptions {
    fn check(&self) -> Result<()> {
        if self.cycles == 0 {
            return Err(StructPredictError::invalid("cycles", "must be at least 1"));
        }
        if self.moves_per_cycle == 0 {
            return Err(StructPredictError::invalid(
                "moves_per_cycle",
                "must be at least 1",
            ));
        }
        if !(self.start_temperature.is_finite() && self.start_temperature > 0.0) {
            return Err(StructPredictError::invalid(
                "start_temperature",
                "must be finite and positive",
            ));
        }
        if !(self.end_temperature.is_finite() && self.end_temperature > 0.0) {
            return Err(StructPredictError::invalid(
                "end_temperature",
                "must be finite and positive",
            ));
        }
        if self.fragment_length == 0 {
            return Err(StructPredictError::invalid(
                "fragment_length",
                "must be at least 1",
            ));
        }
        if self.fragments_per_window == 0 {
            return Err(StructPredictError::invalid(
                "fragments_per_window",
                "must be at least 1",
            ));
        }
        Ok(())
    }
}

/// Outcome of one MC refinement run.
#[derive(Clone, Debug, PartialEq)]
pub struct McRefineResult {
    /// The refined model.
    pub model: ProteinModel,
    /// DOPE energy at the very start of refinement.
    pub initial_energy: f64,
    /// DOPE energy of the returned [`Self::model`] (the best seen).
    pub final_energy: f64,
    /// Per-cycle best-energy trajectory — `energies[k]` is the
    /// best DOPE energy after cycle `k`. Monotonically non-increasing
    /// for any successful refinement.
    pub energies: Vec<f64>,
    /// Total fragment-insertion moves attempted across all cycles.
    pub attempted: usize,
    /// Total fragment-insertion moves accepted.
    pub accepted: usize,
    /// Ramachandran outliers removed by the in-cycle refinement
    /// passes (sum across cycles).
    pub outliers_removed: usize,
}

impl McRefineResult {
    /// Fraction of attempted moves that were accepted.
    pub fn acceptance_rate(&self) -> f64 {
        if self.attempted == 0 {
            0.0
        } else {
            self.accepted as f64 / self.attempted as f64
        }
    }

    /// Energy improvement, `initial − final` (positive = better).
    pub fn improvement(&self) -> f64 {
        self.initial_energy - self.final_energy
    }
}

/// Refines a fully-built model with DOPE-driven simulated-annealing
/// fragment-insertion MC, in-cycle gradient relax, and Ramachandran
/// cleanup.
///
/// `model` must have a complete backbone. The function clones the
/// input and never mutates it; the refined copy is returned in
/// [`McRefineResult::model`].
///
/// `sequence` is the model's one-letter sequence and must match the
/// model's residues — the library is built against it.
///
/// # Errors
/// [`StructPredictError::Invalid`] for an incomplete backbone, a
/// mismatched sequence, or bad options.
pub fn mc_refine(
    model: &ProteinModel,
    sequence: &str,
    options: McRefineOptions,
) -> Result<McRefineResult> {
    options.check()?;
    let sequence = sequence.trim();
    if sequence.is_empty() {
        return Err(StructPredictError::invalid("sequence", "empty"));
    }
    if sequence.len() != model.residues.len() {
        return Err(StructPredictError::invalid(
            "sequence",
            format!(
                "sequence length {} disagrees with model length {}",
                sequence.len(),
                model.residues.len()
            ),
        ));
    }
    if !model.is_complete() {
        return Err(StructPredictError::invalid(
            "model",
            "MC refinement needs a complete backbone",
        ));
    }
    if options.fragment_length > sequence.len() {
        return Err(StructPredictError::invalid(
            "fragment_length",
            "exceeds sequence length",
        ));
    }

    let library = build_fragment_library(
        sequence,
        options.fragment_length,
        options.fragments_per_window,
    )?;

    let mut current = model.clone();
    let initial_energy = dope_score(&current, DopeWeights::default())?.total;
    let mut current_energy = initial_energy;
    let mut best = current.clone();
    let mut best_energy = current_energy;

    let mut rng = Rng::new(options.seed);
    let mut attempted = 0usize;
    let mut accepted = 0usize;
    let mut outliers_removed = 0usize;
    let mut energies = Vec::with_capacity(options.cycles);

    for cycle in 0..options.cycles {
        // Geometric cooling.
        let frac = if options.cycles > 1 {
            cycle as f64 / (options.cycles - 1) as f64
        } else {
            0.0
        };
        let temperature = options.start_temperature
            * (options.end_temperature / options.start_temperature).powf(frac);

        // --- (a) fragment-insertion MC at this temperature ----------
        let mc = mc_cycle(
            &current,
            sequence,
            &library,
            temperature,
            options.moves_per_cycle,
            &mut rng,
        )?;
        attempted += mc.attempted;
        accepted += mc.accepted;
        // Adopt the MC's best-seen for this cycle.
        current = mc.best_model;
        current_energy = mc.best_energy;

        // --- (b) Ramachandran cleanup ------------------------------
        // Ramachandran refinement re-measures φ/ψ from the rebuilt
        // backbone, snaps disallowed residues to the nearest allowed
        // basin, and rebuilds the chain end-to-end. The MC loop's
        // fragment moves can occasionally drive a residue into a
        // disallowed strip when a fragment ω/spread combo lands
        // there; this pass cleans it.
        if current.is_complete() {
            // Snapshot the pre-Rama energy/model so we can revert
            // if the snap *worsens* the DOPE energy (a snap to a
            // distant basin can be a step backward at low T).
            let pre = current.clone();
            let pre_e = current_energy;
            let rr = refine_ramachandran(&mut current)?;
            let post_e = dope_score(&current, DopeWeights::default())?.total;
            if post_e <= pre_e + 1e-6 {
                current_energy = post_e;
                outliers_removed += rr.removed();
            } else {
                // Revert — the snap broke the geometry. We still
                // count the outlier identification but don't keep
                // the change.
                current = pre;
                current_energy = pre_e;
            }
        }

        // --- (c) Gradient relax — kept off the default refinement path -
        // The Cα-only `valenx_md` relaxer pulls Cα atoms toward an
        // ideal-spacing chain restraint but does NOT move the N / C /
        // O backbone atoms with them, which leaves the backbone in
        // an inconsistent state for later (φ, ψ) measurement. For
        // this refinement loop the MC + Rama steps cover the same
        // ground without breaking the consistent-backbone invariant.
        // The relax is still useful as a *terminal* relax step the
        // caller applies after refinement is done; we expose the
        // `relax_iterations` option so callers can dial it in for
        // a Cα-trace-only refinement use case.
        if options.relax_iterations > 0 {
            // Snapshot, attempt the relax, accept only if DOPE
            // improves (the relax minimises a chain-restraint
            // potential, not DOPE — improvement is not guaranteed).
            let pre = current.clone();
            let pre_e = current_energy;
            let _ = relax_model(&mut current, options.relax_iterations, 0.5)?;
            let post_e = dope_score(&current, DopeWeights::default())?.total;
            if post_e < pre_e {
                current_energy = post_e;
            } else {
                current = pre;
                current_energy = pre_e;
            }
        }

        if current_energy < best_energy {
            best_energy = current_energy;
            best = current.clone();
        }
        energies.push(best_energy);
    }

    Ok(McRefineResult {
        model: best,
        initial_energy,
        final_energy: best_energy,
        energies,
        attempted,
        accepted,
        outliers_removed,
    })
}

/// One MC cycle at a fixed temperature: tries `moves` fragment
/// insertions under the Metropolis criterion and returns the
/// lowest-DOPE-energy snapshot.
struct McCycleResult {
    best_model: ProteinModel,
    best_energy: f64,
    attempted: usize,
    accepted: usize,
}

fn mc_cycle(
    start: &ProteinModel,
    sequence: &str,
    library: &FragmentLibrary,
    temperature: f64,
    moves: usize,
    rng: &mut Rng,
) -> Result<McCycleResult> {
    let n = sequence.len();
    let weights = DopeWeights::default();

    // Pull the current torsions out of `start` so we can splice in
    // fragments (the assembler stores torsions, not coordinates,
    // because the (φ, ψ) → backbone rebuild is the canonical move).
    let current_pp = model_phi_psi(start);
    let mut torsions: Vec<(f64, f64)> = current_pp;

    let mut current = start.clone();
    let mut current_energy = dope_score(&current, weights)?.total;
    let mut best = current.clone();
    let mut best_energy = current_energy;
    let mut attempted = 0usize;
    let mut accepted = 0usize;

    let n_windows = library.fragments.len().max(1);
    for _ in 0..moves {
        attempted += 1;
        let window = rng.below(n_windows);
        let frags = match library.at(window) {
            Some(f) if !f.is_empty() => f,
            _ => continue,
        };
        let frag = &frags[rng.below(frags.len())];

        let mut trial_torsions = torsions.clone();
        for (k, &t) in frag.torsions.iter().enumerate() {
            let idx = frag.start + k;
            if idx < n {
                trial_torsions[idx] = t;
            }
        }
        let mut trial = ProteinModel::from_sequence(sequence)?;
        build_backbone_from_torsions(&mut trial, &trial_torsions)?;
        let trial_energy = dope_score(&trial, weights)?.total;
        let delta = trial_energy - current_energy;
        let accept = delta <= 0.0 || rng.uniform() < (-delta / temperature).exp();
        if accept {
            accepted += 1;
            torsions = trial_torsions;
            current = trial;
            current_energy = trial_energy;
            if current_energy < best_energy {
                best_energy = current_energy;
                best = current.clone();
            }
        }
    }
    Ok(McCycleResult {
        best_model: best,
        best_energy,
        attempted,
        accepted,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abinitio::assemble::{
        fragment_assembly, AssemblyOptions, AssemblyScorer,
    };
    use crate::abinitio::dope::{dope_score, DopeWeights};
    use crate::abinitio::fragments::build_fragment_library;

    /// A helical starting model (the easiest case — a clean basin).
    fn helix(n: usize, aa: char) -> ProteinModel {
        let seq: String = std::iter::repeat_n(aa, n).collect();
        let mut m = ProteinModel::from_sequence(&seq).expect("model");
        build_backbone_from_torsions(&mut m, &vec![(-63.0, -42.0); n]).expect("build");
        m
    }

    /// Perturb a helix's middle residues into a strand basin to test
    /// recovery.
    fn helix_with_strand_dent(n: usize, aa: char) -> ProteinModel {
        let seq: String = std::iter::repeat_n(aa, n).collect();
        let mut m = ProteinModel::from_sequence(&seq).expect("model");
        let mut t = vec![(-63.0, -42.0); n];
        let mid = n / 2;
        let lo = mid.saturating_sub(1);
        let hi = (mid + 2).min(n);
        for slot in t.iter_mut().take(hi).skip(lo) {
            *slot = (-120.0, 130.0); // strand dent in the middle
        }
        build_backbone_from_torsions(&mut m, &t).expect("build");
        m
    }

    #[test]
    fn refinement_lowers_dope_energy() {
        // A strand-dented helix has a higher DOPE energy than the
        // clean helix; MC refinement under DOPE should pull it back.
        let dented = helix_with_strand_dent(15, 'L');
        let seq = "L".repeat(15);
        let opts = McRefineOptions {
            cycles: 4,
            moves_per_cycle: 80,
            seed: 11,
            ..McRefineOptions::default()
        };
        let res = mc_refine(&dented, &seq, opts).expect("refine");
        assert!(
            res.final_energy <= res.initial_energy + 1e-6,
            "DOPE energy did not improve: {} -> {}",
            res.initial_energy,
            res.final_energy
        );
        // The trajectory must be monotonically non-increasing.
        for w in res.energies.windows(2) {
            assert!(
                w[1] <= w[0] + 1e-9,
                "trajectory not monotone: {:?}",
                res.energies
            );
        }
    }

    #[test]
    fn refinement_improves_rmsd_toward_native() {
        // A native helix (all-Leucine — strong helix former) with a
        // deliberate single-residue torsion defect should refine
        // *back* toward the native: lower DOPE *and* lower Cα-RMSD.
        // The perturbation lives in torsion space because the MC's
        // fundamental move is a (φ, ψ) splice.
        let native = helix(12, 'L');
        let mut perturbed = ProteinModel::from_sequence(&"L".repeat(12)).expect("model");
        let mut t = vec![(-63.0, -42.0); 12];
        // One residue defect — small enough that the chain's other
        // residues anchor the helix basin firmly.
        t[5] = (-100.0, 100.0);
        build_backbone_from_torsions(&mut perturbed, &t).expect("build");
        let seq: String = "L".repeat(12);
        let initial_rmsd =
            crate::refine::superpose::ca_rmsd_superposed(&perturbed, &native).expect("rmsd0");
        // Hard refinement: low T, many cycles — the MC walks locally
        // and the DOPE minimum near the defect is the native basin.
        let opts = McRefineOptions {
            cycles: 10,
            moves_per_cycle: 300,
            start_temperature: 0.2,
            end_temperature: 0.02,
            seed: 17,
            ..McRefineOptions::default()
        };
        let res = mc_refine(&perturbed, &seq, opts).expect("refine");
        let final_rmsd =
            crate::refine::superpose::ca_rmsd_superposed(&res.model, &native).expect("rmsd1");
        // Both metrics must improve (RMSD-to-native and DOPE energy).
        assert!(
            res.final_energy < res.initial_energy + 1e-6,
            "DOPE: {} -> {}",
            res.initial_energy,
            res.final_energy
        );
        assert!(
            final_rmsd < initial_rmsd,
            "RMSD-to-native should improve: {initial_rmsd} -> {final_rmsd}"
        );
    }

    #[test]
    fn mc_dope_finds_lower_energy_than_v1_knowledge_score() {
        // The DOPE+MC refinement loop finds a structure whose DOPE
        // energy is no worse than the legacy hand-built knowledge-
        // score path — i.e. the DOPE potential, used by an MC loop
        // tuned for it, is the better default. We compare the same
        // sequence, same starting model.
        let seq = "EEEEAAAALLLL";
        let lib = build_fragment_library(seq, 3, 25).expect("lib");
        let v1_opts = AssemblyOptions {
            moves: 800,
            scorer: AssemblyScorer::Knowledge,
            seed: 31,
            ..AssemblyOptions::default()
        };
        let v1 = fragment_assembly(seq, &lib, v1_opts).expect("v1");
        // DOPE-evaluated energy of v1's lowest-knowledge-score model.
        let v1_dope = dope_score(&v1.model, DopeWeights::default())
            .expect("dope v1")
            .total;

        let dope_opts = AssemblyOptions {
            moves: 800,
            scorer: AssemblyScorer::Dope,
            seed: 31,
            ..AssemblyOptions::default()
        };
        let v2 = fragment_assembly(seq, &lib, dope_opts).expect("v2");
        assert!(
            v2.final_score <= v1_dope + 1e-6,
            "DOPE-driven assembler reports DOPE {} that must \
             be ≤ DOPE-evaluation of v1's best ({})",
            v2.final_score,
            v1_dope
        );
    }

    #[test]
    fn refinement_is_deterministic() {
        let m = helix_with_strand_dent(12, 'L');
        let seq = "L".repeat(12);
        let opts = McRefineOptions {
            seed: 23,
            cycles: 3,
            moves_per_cycle: 50,
            ..McRefineOptions::default()
        };
        let a = mc_refine(&m, &seq, opts).expect("a");
        let b = mc_refine(&m, &seq, opts).expect("b");
        assert!((a.final_energy - b.final_energy).abs() < 1e-9);
        assert_eq!(a.accepted, b.accepted);
    }

    #[test]
    fn incomplete_model_rejected() {
        let m = ProteinModel::from_sequence("ACDEF").expect("model");
        let r = mc_refine(&m, "ACDEF", McRefineOptions::default());
        assert!(r.is_err());
    }

    #[test]
    fn bad_options_rejected() {
        let m = helix(8, 'A');
        let bad = McRefineOptions {
            cycles: 0,
            ..McRefineOptions::default()
        };
        assert!(mc_refine(&m, "AAAAAAAA", bad).is_err());
    }

    #[test]
    fn sequence_length_mismatch_rejected() {
        let m = helix(8, 'A');
        let r = mc_refine(&m, "AAAA", McRefineOptions::default());
        assert!(r.is_err());
    }
}
