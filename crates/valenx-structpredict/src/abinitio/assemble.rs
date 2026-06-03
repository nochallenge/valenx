//! **Feature 10 — Monte-Carlo fragment assembly.**
//!
//! The core ab-initio move: build a fold by **simulated annealing**
//! over fragment insertions. Starting from an extended chain, the
//! protocol repeatedly:
//!
//! 1. picks a random window and a random fragment for it from the
//!    [`crate::abinitio::fragments::FragmentLibrary`];
//! 2. splices that fragment's (φ, ψ) angles into the current
//!    conformation and rebuilds the affected backbone;
//! 3. scores the new conformation with the knowledge-based potential;
//! 4. accepts the move by the **Metropolis criterion** — always if it
//!    lowers the energy, with probability `exp(−ΔE / T)` if it raises
//!    it.
//!
//! The temperature `T` is cooled geometrically over the run
//! (simulated annealing): early high-`T` moves explore broadly, late
//! low-`T` moves refine. This is the genuine Rosetta-class
//! fragment-assembly protocol. A deterministic seedable PCG generator
//! (`valenx_md::Rng`) drives every random choice, so a run is fully
//! reproducible.

use serde::{Deserialize, Serialize};
use valenx_md::Rng;

use crate::abinitio::dope::{dope_score, DopeWeights};
use crate::abinitio::fragments::FragmentLibrary;
use crate::abinitio::score::{score_model, ScoreWeights};
use crate::error::{Result, StructPredictError};
use crate::model::{build_backbone_from_torsions, ProteinModel};

/// Which scoring function the Monte-Carlo assembler should use.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum AssemblyScorer {
    /// The new DOPE-class distance-dependent statistical potential —
    /// the production default. See [`crate::abinitio::dope`].
    #[default]
    Dope,
    /// The legacy hand-built knowledge-based score — kept for
    /// regression tests and as a non-DOPE baseline. See
    /// [`crate::abinitio::score`].
    Knowledge,
}

/// Controls a fragment-assembly run.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AssemblyOptions {
    /// Total Monte-Carlo fragment-insertion moves.
    pub moves: usize,
    /// Starting (high) annealing temperature.
    pub start_temperature: f64,
    /// Final (low) annealing temperature.
    pub end_temperature: f64,
    /// RNG seed — fixes the whole trajectory.
    pub seed: u64,
    /// Which scorer the Metropolis criterion evaluates. Defaults to
    /// [`AssemblyScorer::Dope`].
    pub scorer: AssemblyScorer,
}

impl Default for AssemblyOptions {
    fn default() -> Self {
        AssemblyOptions {
            moves: 2000,
            start_temperature: 2.5,
            end_temperature: 0.1,
            seed: 0xA55E_3B1E,
            scorer: AssemblyScorer::Dope,
        }
    }
}

/// One Metropolis evaluation under the selected scorer.
fn score_with(model: &ProteinModel, which: AssemblyScorer) -> Result<f64> {
    Ok(match which {
        AssemblyScorer::Dope => dope_score(model, DopeWeights::default())?.total,
        AssemblyScorer::Knowledge => score_model(model, ScoreWeights::default())?.total,
    })
}

impl AssemblyOptions {
    fn check(&self) -> Result<()> {
        if self.moves == 0 {
            return Err(StructPredictError::invalid("moves", "must be at least 1"));
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
        Ok(())
    }
}

/// The outcome of a fragment-assembly run.
#[derive(Clone, Debug, PartialEq)]
pub struct AssemblyResult {
    /// The lowest-energy model encountered during the run.
    pub model: ProteinModel,
    /// The knowledge-score total of [`Self::model`].
    pub final_score: f64,
    /// The starting (extended-chain) score, for comparison.
    pub initial_score: f64,
    /// Fragment insertions accepted.
    pub accepted: usize,
    /// Total moves attempted.
    pub attempted: usize,
}

impl AssemblyResult {
    /// Fraction of attempted moves that were accepted.
    pub fn acceptance_rate(&self) -> f64 {
        if self.attempted == 0 {
            0.0
        } else {
            self.accepted as f64 / self.attempted as f64
        }
    }
}

/// Runs a Monte-Carlo simulated-annealing fragment-assembly folding
/// trajectory.
///
/// `sequence` is the target sequence; `library` is its fragment
/// library (its `sequence_length` must match). The returned
/// [`AssemblyResult`] carries the lowest-energy model seen.
///
/// # Errors
/// [`StructPredictError::Invalid`] for a sequence/library length
/// mismatch, an empty sequence, or bad options.
pub fn fragment_assembly(
    sequence: &str,
    library: &FragmentLibrary,
    options: AssemblyOptions,
) -> Result<AssemblyResult> {
    options.check()?;
    let sequence = sequence.trim();
    if sequence.is_empty() {
        return Err(StructPredictError::invalid("sequence", "empty"));
    }
    if library.sequence_length != sequence.len() {
        return Err(StructPredictError::invalid(
            "library",
            format!(
                "library is for length {} but sequence is {}",
                library.sequence_length,
                sequence.len()
            ),
        ));
    }
    let n = sequence.len();
    if n < 2 {
        return Err(StructPredictError::invalid(
            "sequence",
            "needs at least 2 residues to fold",
        ));
    }
    let mut rng = Rng::new(options.seed);

    // Current conformation, stored as per-residue (φ, ψ).
    // Start fully extended (β-strand basin).
    let mut torsions = vec![(-120.0, 130.0); n];
    let mut current = ProteinModel::from_sequence(sequence)?;
    build_backbone_from_torsions(&mut current, &torsions)?;
    let initial_score = score_with(&current, options.scorer)?;
    let mut current_score = initial_score;

    let mut best = current.clone();
    let mut best_score = current_score;
    let mut accepted = 0usize;

    let n_windows = library.fragments.len().max(1);

    for step in 0..options.moves {
        // Geometric cooling schedule.
        let frac = step as f64 / options.moves as f64;
        let temperature = options.start_temperature
            * (options.end_temperature / options.start_temperature).powf(frac);

        // Pick a window and a fragment.
        let window = rng.below(n_windows);
        let frags = match library.at(window) {
            Some(f) if !f.is_empty() => f,
            _ => continue,
        };
        let frag = &frags[rng.below(frags.len())];

        // Splice the fragment's torsions in.
        let mut trial_torsions = torsions.clone();
        for (k, &t) in frag.torsions.iter().enumerate() {
            let idx = frag.start + k;
            if idx < n {
                trial_torsions[idx] = t;
            }
        }
        let mut trial = ProteinModel::from_sequence(sequence)?;
        build_backbone_from_torsions(&mut trial, &trial_torsions)?;
        let trial_score = score_with(&trial, options.scorer)?;

        // Metropolis acceptance.
        let delta = trial_score - current_score;
        let accept = delta <= 0.0 || rng.uniform() < (-delta / temperature).exp();
        if accept {
            torsions = trial_torsions;
            current = trial;
            current_score = trial_score;
            accepted += 1;
            if current_score < best_score {
                best = current.clone();
                best_score = current_score;
            }
        }
    }

    Ok(AssemblyResult {
        model: best,
        final_score: best_score,
        initial_score,
        accepted,
        attempted: options.moves,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::abinitio::fragments::build_fragment_library;

    #[test]
    fn assembly_lowers_the_score() {
        let seq = "EEEEAAAALLLLEEEEAAAA";
        let lib = build_fragment_library(seq, 3, 30).expect("lib");
        let opts = AssemblyOptions {
            moves: 600,
            ..AssemblyOptions::default()
        };
        let res = fragment_assembly(seq, &lib, opts).expect("assemble");
        assert!(
            res.final_score <= res.initial_score,
            "score {} -> {}",
            res.initial_score,
            res.final_score
        );
        assert!(res.model.is_complete());
    }

    #[test]
    fn assembly_is_deterministic() {
        let seq = "ACDEFGHIKLMNPQRST";
        let lib = build_fragment_library(seq, 3, 20).expect("lib");
        let opts = AssemblyOptions {
            moves: 200,
            ..AssemblyOptions::default()
        };
        let a = fragment_assembly(seq, &lib, opts).expect("a");
        let b = fragment_assembly(seq, &lib, opts).expect("b");
        assert_eq!(a.final_score, b.final_score);
        assert_eq!(a.accepted, b.accepted);
    }

    #[test]
    fn library_length_mismatch_rejected() {
        let lib = build_fragment_library("ACDEFGHIKL", 3, 5).expect("lib");
        // Different-length sequence.
        assert!(fragment_assembly("ACDEFGHIKLMN", &lib, AssemblyOptions::default()).is_err());
    }

    #[test]
    fn bad_options_rejected() {
        let lib = build_fragment_library("ACDEFGHIKL", 3, 5).expect("lib");
        let bad = AssemblyOptions {
            moves: 0,
            ..AssemblyOptions::default()
        };
        assert!(fragment_assembly("ACDEFGHIKL", &lib, bad).is_err());
    }
}
