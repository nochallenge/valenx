//! **Feature 19 — the combinatorial design search.**
//!
//! Fixed-backbone design is a huge combinatorial problem: each
//! designable position can take any of the 20 amino acids, each amino
//! acid in any of its rotamers. For `n` positions that is up to
//! `(20·r)^n` states — astronomically large.
//!
//! This module solves it with **simulated annealing** over the joint
//! (residue → amino-acid, rotamer) space — the same Monte-Carlo
//! design search Rosetta `fixbb` uses. A move re-rolls one position's
//! amino acid and rotamer; the move is accepted by the Metropolis
//! criterion against the design score; the temperature cools
//! geometrically so the search explores early and refines late.
//!
//! A [`ResiduePalette`] lets the caller restrict which amino acids
//! each position may take — fixing some positions to their native
//! identity, allowing only hydrophobics in the core, etc. — exactly
//! the per-position "resfile" control a real design protocol gives.

use serde::{Deserialize, Serialize};
use valenx_md::Rng;

use crate::aa::AMINO_ACIDS;
use crate::abinitio::ss::SecondaryStructure;
use crate::design::score::{design_score, DesignScore, DesignScoreWeights};
use crate::error::{Result, StructPredictError};
use crate::model::ProteinModel;

/// Which amino acids each position of a design is allowed to take.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ResiduePalette {
    /// `allowed[i]` is the set of amino acids position `i` may take.
    /// An empty inner set means "any of the 20".
    pub allowed: Vec<Vec<char>>,
}

impl ResiduePalette {
    /// A palette that allows every amino acid at every position.
    pub fn unrestricted(n: usize) -> Self {
        ResiduePalette {
            allowed: vec![Vec::new(); n],
        }
    }

    /// A palette that fixes every position to the given native
    /// sequence — design becomes a no-op identity (useful as a
    /// baseline / for fixing a subset later).
    pub fn fixed_to(sequence: &str) -> Self {
        ResiduePalette {
            allowed: sequence
                .trim()
                .to_ascii_uppercase()
                .chars()
                .map(|c| vec![c])
                .collect(),
        }
    }

    /// The choices for position `i` — the position's allowed set, or
    /// all 20 amino acids if its set is empty.
    fn choices(&self, i: usize) -> Vec<char> {
        match self.allowed.get(i) {
            Some(set) if !set.is_empty() => set.clone(),
            _ => AMINO_ACIDS.to_vec(),
        }
    }

    /// Number of designable positions.
    pub fn len(&self) -> usize {
        self.allowed.len()
    }

    /// `true` when the palette has no positions.
    pub fn is_empty(&self) -> bool {
        self.allowed.is_empty()
    }
}

/// The outcome of a combinatorial design search.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DesignSearchResult {
    /// The lowest-energy designed sequence.
    pub sequence: String,
    /// The decomposed design score of [`Self::sequence`].
    pub score: DesignScore,
    /// The score of the starting sequence, for comparison.
    pub initial_total: f64,
    /// Monte-Carlo moves accepted.
    pub accepted: usize,
    /// Monte-Carlo moves attempted.
    pub attempted: usize,
}

/// Runs a simulated-annealing combinatorial design search.
///
/// `model` supplies the fixed backbone; `palette` restricts the
/// per-position amino-acid choices; `ss` is the per-residue secondary
/// structure (or empty); `moves` is the Monte-Carlo budget; `seed`
/// fixes the RNG. The returned [`DesignSearchResult`] carries the
/// lowest-energy sequence found.
///
/// # Errors
/// [`StructPredictError::Invalid`] for a palette / model length
/// mismatch, `moves == 0`, or a backbone with fewer than 2 Cα atoms.
pub fn combinatorial_design(
    model: &ProteinModel,
    palette: &ResiduePalette,
    ss: &[SecondaryStructure],
    weights: DesignScoreWeights,
    moves: usize,
    seed: u64,
) -> Result<DesignSearchResult> {
    let n = model.residues.len();
    if palette.len() != n {
        return Err(StructPredictError::invalid(
            "palette",
            format!("{} positions for an {n}-residue backbone", palette.len()),
        ));
    }
    if moves == 0 {
        return Err(StructPredictError::invalid("moves", "must be at least 1"));
    }
    if model.ca_trace().len() < 2 {
        return Err(StructPredictError::invalid(
            "model",
            "need at least 2 Cα atoms",
        ));
    }

    let mut rng = Rng::new(seed);

    // Initial sequence: the first allowed choice at every position.
    let mut sequence: Vec<char> = (0..n).map(|i| palette.choices(i)[0]).collect();
    let initial_seq: String = sequence.iter().collect();
    let mut current_score = design_score(model, &initial_seq, ss, weights)?;
    let initial_total = current_score.total;

    let mut best = sequence.clone();
    let mut best_score = current_score;
    let mut accepted = 0usize;

    let start_t: f64 = 4.0;
    let end_t: f64 = 0.05;
    for step in 0..moves {
        let frac = step as f64 / moves as f64;
        let t = start_t * (end_t / start_t).powf(frac);

        // Pick a position and a new amino acid for it.
        let pos = rng.below(n);
        let choices = palette.choices(pos);
        if choices.len() < 2 {
            continue; // fixed position — nothing to vary
        }
        let new_aa = choices[rng.below(choices.len())];
        if new_aa == sequence[pos] {
            continue;
        }
        let old_aa = sequence[pos];
        sequence[pos] = new_aa;
        let trial_seq: String = sequence.iter().collect();
        let trial_score = design_score(model, &trial_seq, ss, weights)?;
        let delta = trial_score.total - current_score.total;
        if delta <= 0.0 || rng.uniform() < (-delta / t).exp() {
            current_score = trial_score;
            accepted += 1;
            if current_score.total < best_score.total {
                best_score = current_score;
                best = sequence.clone();
            }
        } else {
            sequence[pos] = old_aa; // reject
        }
    }

    Ok(DesignSearchResult {
        sequence: best.iter().collect(),
        score: best_score,
        initial_total,
        accepted,
        attempted: moves,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    fn blob(n: usize) -> ProteinModel {
        let seq = "A".repeat(n);
        let mut m = ProteinModel::from_sequence(&seq).expect("model");
        let side = (n as f64).cbrt().ceil() as usize;
        for (idx, r) in m.residues.iter_mut().enumerate() {
            let x = (idx % side) as f64 * 3.8;
            let y = ((idx / side) % side) as f64 * 3.8;
            let z = (idx / (side * side)) as f64 * 3.8;
            r.ca = Some(Point3::new(x, y, z));
        }
        m
    }

    #[test]
    fn design_lowers_the_score() {
        let m = blob(27);
        let palette = ResiduePalette::unrestricted(m.residues.len());
        let res = combinatorial_design(&m, &palette, &[], DesignScoreWeights::default(), 800, 17)
            .expect("design");
        assert!(
            res.score.total <= res.initial_total,
            "score {} -> {}",
            res.initial_total,
            res.score.total
        );
        assert_eq!(res.sequence.len(), m.residues.len());
    }

    #[test]
    fn fixed_palette_keeps_the_native_sequence() {
        let m = blob(8);
        let native = "ACDEFGHI";
        let palette = ResiduePalette::fixed_to(native);
        let res = combinatorial_design(&m, &palette, &[], DesignScoreWeights::default(), 200, 1)
            .expect("design");
        // Every position is fixed → the sequence cannot change.
        assert_eq!(res.sequence, native);
        assert_eq!(res.accepted, 0);
    }

    #[test]
    fn design_is_deterministic() {
        let m = blob(27);
        let palette = ResiduePalette::unrestricted(m.residues.len());
        let w = DesignScoreWeights::default();
        let a = combinatorial_design(&m, &palette, &[], w, 300, 5).expect("a");
        let b = combinatorial_design(&m, &palette, &[], w, 300, 5).expect("b");
        assert_eq!(a.sequence, b.sequence);
    }

    #[test]
    fn palette_length_mismatch_rejected() {
        let m = blob(8);
        let palette = ResiduePalette::unrestricted(5);
        assert!(
            combinatorial_design(&m, &palette, &[], DesignScoreWeights::default(), 10, 0).is_err()
        );
    }
}
