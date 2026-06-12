//! RNA design — inverse folding and G-quadruplex prediction.
//!
//! Two "design / motif" tools:
//!
//! - **Inverse folding** ([`inverse_fold`]) — the Eterna problem:
//!   given a *target* secondary structure, find an RNA *sequence*
//!   that folds to it. This is the inverse of [`crate::fold::zuker`].
//!   The v1 is an adaptive walk: start from a random sequence
//!   compatible with the target, repeatedly mutate the worst-fitting
//!   position and keep the change if it lowers the base-pair distance
//!   to the target.
//! - **G-quadruplex prediction** ([`predict_gquadruplex`]) — G4s are
//!   four-stranded motifs formed by four runs of consecutive Gs. They
//!   are detected by the canonical G-run regular pattern and scored
//!   by run length and loop lengths (a G4Hunter / QGRS-class rule).

use crate::ensemble::rng::Rng;
use crate::error::{Result, RnaStructError};
use crate::fold::zuker::mfe;
use crate::rna::RnaSeq;
use crate::structure::Structure;

// ---------------------------------------------------------------------
// Inverse folding
// ---------------------------------------------------------------------

/// The bases the design walk may place, encoded as ASCII.
const BASES: [u8; 4] = [b'A', b'C', b'G', b'U'];
/// Canonical pairing partners the walk uses when seating a pair
/// (Watson-Crick + the G-U wobble), as ASCII pairs.
const PAIRS: [(u8, u8); 6] = [
    (b'A', b'U'),
    (b'U', b'A'),
    (b'G', b'C'),
    (b'C', b'G'),
    (b'G', b'U'),
    (b'U', b'G'),
];

/// The outcome of an inverse-folding run.
#[derive(Clone, Debug)]
pub struct InverseFoldResult {
    /// The designed RNA sequence.
    pub sequence: RnaSeq,
    /// The structure that sequence actually folds into (its MFE).
    pub achieved: Structure,
    /// The base-pair distance between [`achieved`](Self::achieved) and
    /// the target — `0` means the design folds exactly to the target.
    pub distance: usize,
    /// `true` if the design folds *exactly* to the target.
    pub solved: bool,
    /// Number of adaptive-walk iterations performed.
    pub iterations: usize,
}

/// Designs a sequence that folds to `target`, with default effort.
///
/// # Errors
/// [`RnaStructError::Structure`] if the target is empty or
/// pseudoknotted (the MFE folder is pseudoknot-free, so a
/// pseudoknotted target is unreachable).
pub fn inverse_fold(target: &Structure, seed: u64) -> Result<InverseFoldResult> {
    inverse_fold_with(target, seed, 2000)
}

/// [`inverse_fold`] with an explicit iteration budget.
///
/// # Errors
/// As [`inverse_fold`].
pub fn inverse_fold_with(
    target: &Structure,
    seed: u64,
    max_iterations: usize,
) -> Result<InverseFoldResult> {
    let n = target.len();
    if n == 0 {
        return Err(RnaStructError::structure(
            "cannot inverse-fold an empty target",
        ));
    }
    if target.has_pseudoknot() {
        return Err(RnaStructError::structure(
            "target is pseudoknotted — the MFE folder cannot reach it",
        ));
    }

    let mut rng = Rng::new(seed);
    // Initial sequence: seat every pair with a random canonical pair,
    // fill unpaired positions randomly.
    let mut ascii = vec![b'A'; n];
    let mut assigned = vec![false; n];
    for bp in target.pairs() {
        let (a, b) = PAIRS[(rng.next_u64() % 6) as usize];
        ascii[bp.i] = a;
        ascii[bp.j] = b;
        assigned[bp.i] = true;
        assigned[bp.j] = true;
    }
    for (i, a) in ascii.iter_mut().enumerate() {
        if !assigned[i] {
            *a = BASES[(rng.next_u64() % 4) as usize];
        }
    }

    let mut current = RnaSeq::parse(&ascii)?;
    let mut best_dist = fold_distance(&current, target)?;
    let mut iterations = 0;

    while best_dist > 0 && iterations < max_iterations {
        iterations += 1;
        // Pick a position that is currently mis-predicted and mutate
        // it (and its partner, if paired).
        let folded = mfe(&current)?.structure;
        let mismatched = mismatched_positions(&folded, target);
        let pos = if mismatched.is_empty() {
            (rng.next_u64() as usize) % n
        } else {
            mismatched[(rng.next_u64() as usize) % mismatched.len()]
        };

        let mut trial = ascii.clone();
        match target.partner(pos) {
            Some(partner) => {
                // Re-seat the whole pair with a fresh canonical pair.
                let (a, b) = PAIRS[(rng.next_u64() % 6) as usize];
                let (lo, hi) = if pos < partner {
                    (pos, partner)
                } else {
                    (partner, pos)
                };
                trial[lo] = a;
                trial[hi] = b;
            }
            None => {
                // Mutate a single unpaired base.
                trial[pos] = BASES[(rng.next_u64() % 4) as usize];
            }
        }

        let trial_seq = RnaSeq::parse(&trial)?;
        let trial_dist = fold_distance(&trial_seq, target)?;
        if trial_dist <= best_dist {
            ascii = trial;
            current = trial_seq;
            best_dist = trial_dist;
        }
    }

    let achieved = mfe(&current)?.structure;
    Ok(InverseFoldResult {
        sequence: current,
        achieved,
        distance: best_dist,
        solved: best_dist == 0,
        iterations,
    })
}

/// Base-pair distance between the MFE fold of `seq` and `target`.
fn fold_distance(seq: &RnaSeq, target: &Structure) -> Result<usize> {
    let folded = mfe(seq)?.structure;
    crate::compare::distance::base_pair_distance(&folded, target)
}

/// Positions where the folded structure disagrees with the target —
/// either paired-vs-unpaired or paired to a different partner.
fn mismatched_positions(folded: &Structure, target: &Structure) -> Vec<usize> {
    let n = folded.len().min(target.len());
    (0..n)
        .filter(|&i| folded.partner(i) != target.partner(i))
        .collect()
}

// ---------------------------------------------------------------------
// G-quadruplex prediction
// ---------------------------------------------------------------------

/// One predicted G-quadruplex (G4) motif.
#[derive(Clone, Debug, PartialEq)]
pub struct GQuadruplex {
    /// Start index of the motif (inclusive).
    pub start: usize,
    /// End index of the motif (exclusive).
    pub end: usize,
    /// The length of each of the four G-runs (the G-tract length).
    pub g_run_length: usize,
    /// The three loop lengths between the four G-runs.
    pub loop_lengths: [usize; 3],
    /// A heuristic stability score — higher is a more stable G4.
    pub score: f64,
}

impl GQuadruplex {
    /// The motif length in nucleotides.
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    /// `true` if the motif spans zero nucleotides (never produced).
    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }
}

/// Parameters for G-quadruplex detection.
#[derive(Copy, Clone, Debug)]
pub struct GQuadParams {
    /// Minimum length of each G-run (a "G-tract"). Canonical G4s use
    /// runs of ≥ 3 Gs.
    pub min_g_run: usize,
    /// Maximum length of a loop between two G-runs.
    pub max_loop: usize,
    /// Minimum loop length (1 — a loop cannot be empty).
    pub min_loop: usize,
}

impl Default for GQuadParams {
    fn default() -> Self {
        GQuadParams {
            min_g_run: 3,
            max_loop: 7,
            min_loop: 1,
        }
    }
}

/// Predicts G-quadruplex motifs in `seq` with default parameters.
///
/// A G4 is four G-runs of length ≥ `min_g_run` separated by three
/// loops. Overlapping candidates are returned; sort / filter by
/// [`GQuadruplex::score`] as needed.
pub fn predict_gquadruplex(seq: &RnaSeq) -> Vec<GQuadruplex> {
    predict_gquadruplex_with(seq, GQuadParams::default())
}

/// [`predict_gquadruplex`] with explicit [`GQuadParams`].
pub fn predict_gquadruplex_with(seq: &RnaSeq, params: GQuadParams) -> Vec<GQuadruplex> {
    let bytes = seq.as_bytes();
    let n = bytes.len();
    let mut out = Vec::new();
    if params.min_g_run == 0 || n < 4 * params.min_g_run + 3 * params.min_loop {
        return out;
    }

    // Find every maximal G-run.
    let mut runs: Vec<(usize, usize)> = Vec::new(); // (start, len)
    let mut i = 0;
    while i < n {
        if bytes[i] == b'G' {
            let start = i;
            while i < n && bytes[i] == b'G' {
                i += 1;
            }
            runs.push((start, i - start));
        } else {
            i += 1;
        }
    }

    // A G4 is four runs r0..r3, each contributing >= min_g_run Gs,
    // with the inter-run gaps acting as loops.
    let qualifying: Vec<&(usize, usize)> = runs
        .iter()
        .filter(|(_, l)| *l >= params.min_g_run)
        .collect();

    for w in qualifying.windows(4) {
        let r: [&(usize, usize); 4] = [w[0], w[1], w[2], w[3]];
        // The G-tract used is min_g_run Gs from the *end* of each run
        // for runs 0..2 and the *start* for the layout; for scoring we
        // treat the tract length as min_g_run (extra Gs add to score).
        let tract = params.min_g_run;
        // Loop lengths = gap between consecutive runs.
        let mut loops = [0usize; 3];
        let mut ok = true;
        for k in 0..3 {
            let end_prev = r[k].0 + r[k].1;
            let start_next = r[k + 1].0;
            if start_next < end_prev {
                ok = false;
                break;
            }
            let gap = start_next - end_prev;
            if gap < params.min_loop || gap > params.max_loop {
                ok = false;
                break;
            }
            loops[k] = gap;
        }
        if !ok {
            continue;
        }
        let start = r[0].0 + r[0].1 - tract.min(r[0].1);
        let end = r[3].0 + tract.min(r[3].1);
        // Score: longer tracts and shorter loops are more stable.
        let extra_g: usize = r.iter().map(|(_, l)| l.saturating_sub(tract)).sum();
        let loop_penalty: f64 = loops.iter().map(|&l| l as f64).sum();
        let score = 4.0 * tract as f64 + extra_g as f64 - 0.5 * loop_penalty;
        out.push(GQuadruplex {
            start,
            end,
            g_run_length: tract,
            loop_lengths: loops,
            score,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_fold_simple_hairpin() {
        // A small hairpin target should be solvable.
        let target = Structure::from_dot_bracket("((((....))))").unwrap();
        let r = inverse_fold(&target, 1).unwrap();
        // The designed sequence has the right length...
        assert_eq!(r.sequence.len(), 12);
        // ...and the walk makes progress (distance is small).
        assert!(
            r.distance <= 2,
            "inverse fold left distance {} (iters {})",
            r.distance,
            r.iterations
        );
    }

    #[test]
    fn inverse_fold_unstructured_target_is_trivial() {
        // An all-unpaired target: any sequence that folds open solves it.
        let target = Structure::empty(10);
        let r = inverse_fold(&target, 5).unwrap();
        assert!(r.solved, "all-unpaired target should always be solvable");
        assert_eq!(r.distance, 0);
    }

    #[test]
    fn inverse_fold_rejects_pseudoknot() {
        let pk = Structure::from_dot_bracket("((..[[..))..]]").unwrap();
        assert!(inverse_fold(&pk, 1).is_err());
    }

    #[test]
    fn inverse_fold_is_deterministic() {
        let target = Structure::from_dot_bracket("(((...)))").unwrap();
        let a = inverse_fold(&target, 42).unwrap();
        let b = inverse_fold(&target, 42).unwrap();
        assert_eq!(a.sequence.as_str(), b.sequence.as_str());
    }

    #[test]
    fn gquadruplex_canonical_motif_detected() {
        // four GGG runs separated by short loops
        let seq = RnaSeq::parse("GGGUUGGGUUGGGUUGGG").unwrap();
        let g4 = predict_gquadruplex(&seq);
        assert!(!g4.is_empty(), "a canonical G4 should be found");
        assert_eq!(g4[0].g_run_length, 3);
        assert_eq!(g4[0].loop_lengths, [2, 2, 2]);
    }

    #[test]
    fn gquadruplex_none_in_g_poor_sequence() {
        let seq = RnaSeq::parse("AAUUAAUUAAUUAAUU").unwrap();
        assert!(predict_gquadruplex(&seq).is_empty());
    }

    #[test]
    fn gquadruplex_long_runs_score_higher() {
        let short = RnaSeq::parse("GGGUGGGUGGGUGGG").unwrap();
        let long = RnaSeq::parse("GGGGGUGGGGGUGGGGGUGGGGG").unwrap();
        let s = predict_gquadruplex(&short);
        let l = predict_gquadruplex(&long);
        assert!(!s.is_empty() && !l.is_empty());
        assert!(l[0].score > s[0].score, "longer G-runs should score higher");
    }

    #[test]
    fn gquadruplex_respects_max_loop() {
        // loops of length 9 exceed the default max of 7
        let seq = RnaSeq::parse("GGGAAAAAAAAAGGGAAAAAAAAAGGGAAAAAAAAAGGG").unwrap();
        assert!(predict_gquadruplex(&seq).is_empty());
    }

    #[test]
    fn pairs_helper_is_canonical() {
        use crate::fold::energy::{can_pair_codes, encode_base};
        for (a, b) in PAIRS {
            assert!(can_pair_codes(
                encode_base(a).unwrap(),
                encode_base(b).unwrap()
            ));
        }
    }
}
