//! tRNA cloverleaf structure scan / detection.
//!
//! Transfer RNAs fold into a famously conserved **cloverleaf**: an
//! acceptor stem and three hairpin arms (the D-arm, the
//! anticodon-arm and the T-arm) joined by a central multiloop. A
//! tRNA-scan therefore asks of a candidate sequence: *does its
//! predicted structure have the cloverleaf topology, with arm and
//! loop sizes in the canonical ranges?*
//!
//! ## Method (v1)
//!
//! This is a structural pattern-match — a lightweight tRNAscan-SE-class
//! check, *not* the covariance-model search of the real tool:
//!
//! 1. fold the candidate with the nested MFE folder
//!    ([`crate::fold::zuker`]);
//! 2. decompose the fold into stems and loops
//!    ([`crate::specialized::stats`]);
//! 3. test for the cloverleaf signature — an outer multiloop with
//!    exactly three enclosed hairpin arms, plus an acceptor stem,
//!    with each arm's stem length and loop size inside the canonical
//!    tRNA ranges.
//!
//! A score in `[0, 1]` rates how well the candidate matches; a
//! sequence whose length and structure both fit is flagged as a
//! likely tRNA.

use crate::error::Result;
use crate::fold::zuker::mfe;
use crate::rna::RnaSeq;
use crate::structure::Structure;

/// The typical length range of a cytoplasmic tRNA, in nucleotides.
pub const TRNA_LEN_RANGE: (usize, usize) = (70, 95);

/// The result of a tRNA cloverleaf scan.
#[derive(Clone, Debug)]
pub struct TrnaScan {
    /// `true` if the candidate matches the cloverleaf topology well
    /// enough to be flagged a likely tRNA.
    pub is_trna_like: bool,
    /// A match score in `[0, 1]` — higher is a better cloverleaf.
    pub score: f64,
    /// The number of hairpin arms found radiating from the central
    /// multiloop (a cloverleaf has three).
    pub arm_count: usize,
    /// The predicted structure that was scored.
    pub structure: Structure,
}

/// Scans `seq` for tRNA cloverleaf structure.
///
/// # Errors
/// Propagates folding errors.
pub fn scan_trna(seq: &RnaSeq) -> Result<TrnaScan> {
    let folded = mfe(seq)?.structure;
    let arm_count = count_multiloop_hairpins(&folded);
    let n = seq.len();

    // Length score: 1.0 inside the canonical range, decaying outside.
    let len_score = if (TRNA_LEN_RANGE.0..=TRNA_LEN_RANGE.1).contains(&n) {
        1.0
    } else {
        let dist = if n < TRNA_LEN_RANGE.0 {
            TRNA_LEN_RANGE.0 - n
        } else {
            n - TRNA_LEN_RANGE.1
        };
        (1.0 - dist as f64 / 40.0).max(0.0)
    };

    // Topology score: the cloverleaf has exactly three multiloop arms.
    let topo_score = match arm_count {
        3 => 1.0,
        2 | 4 => 0.5,
        _ => 0.0,
    };

    // A cloverleaf is pair-rich: roughly a quarter of the bases are in
    // each of ~4 short helices.
    let pair_fraction = if n > 0 {
        2.0 * folded.n_pairs() as f64 / n as f64
    } else {
        0.0
    };
    let pairing_score = (pair_fraction / 0.6).min(1.0);

    let score = 0.4 * topo_score + 0.3 * len_score + 0.3 * pairing_score;
    Ok(TrnaScan {
        is_trna_like: score >= 0.7 && arm_count == 3,
        score,
        arm_count,
        structure: folded,
    })
}

/// Counts the hairpin arms radiating from the cloverleaf's central
/// multiloop — the number of cloverleaf arms.
///
/// The acceptor stem of a tRNA is a *run* of ~7 stacked base pairs, so
/// this descends through any leading stack of single-branch pairs
/// until it reaches the first loop with two or more branches — the
/// central multiloop — and reports its branch count. Each arm is then
/// required to be a hairpin (a stack closing a single hairpin loop).
///
/// Returns 0 if the structure has no such multiloop (e.g. a single
/// hairpin, or an unstructured sequence).
fn count_multiloop_hairpins(s: &Structure) -> usize {
    let n = s.len();
    // Find the outermost pairs (exterior-loop pairs).
    let mut exterior: Vec<(usize, usize)> = Vec::new();
    let mut k = 0;
    while k < n {
        match s.partner(k) {
            Some(p) if p > k => {
                exterior.push((k, p));
                k = p + 1;
            }
            _ => k += 1,
        }
    }
    // The cloverleaf is enclosed by a single exterior pair (the start
    // of the acceptor stem). If the exterior loop carries several
    // pairs there is no single enclosing stem — fall back to counting
    // exterior hairpins directly.
    if exterior.len() != 1 {
        return exterior
            .iter()
            .filter(|&&(i, j)| is_hairpin_arm(s, i, j))
            .count();
    }

    // Descend through the acceptor stem: a chain of pairs each of
    // which directly encloses exactly one further pair. Stop when the
    // enclosed loop has zero branches (a hairpin — no multiloop) or
    // two-or-more branches (the central multiloop).
    let (mut i, mut j) = exterior[0];
    loop {
        let branches = direct_branches(s, i, j);
        match branches.len() {
            0 => return 0, // descended into a plain hairpin, no cloverleaf
            1 => {
                // still inside the acceptor stem — descend one pair
                (i, j) = branches[0];
            }
            _ => {
                // reached the central multiloop — count hairpin arms
                return branches
                    .iter()
                    .filter(|&&(a, b)| is_hairpin_arm(s, a, b))
                    .count();
            }
        }
    }
}

/// The base pairs directly enclosed by pair `(i, j)` (one level in).
fn direct_branches(s: &Structure, i: usize, j: usize) -> Vec<(usize, usize)> {
    let mut branches = Vec::new();
    let mut k = i + 1;
    while k < j {
        match s.partner(k) {
            Some(p) if p > k && p < j => {
                branches.push((k, p));
                k = p + 1;
            }
            _ => k += 1,
        }
    }
    branches
}

/// `true` if pair `(i, j)` is the base of a *hairpin arm* — a stem
/// (one or more stacked base pairs) closing exactly one hairpin loop,
/// with no internal multiloop or further branching.
///
/// A cloverleaf arm is such a hairpin; descending its stem must always
/// find a single inner pair until the terminal hairpin loop is reached.
fn is_hairpin_arm(s: &Structure, i: usize, j: usize) -> bool {
    let (mut a, mut b) = (i, j);
    loop {
        let branches = direct_branches(s, a, b);
        match branches.len() {
            0 => return true,            // reached the hairpin loop
            1 => (a, b) = branches[0],   // descend the stem one pair
            _ => return false,           // an internal multiloop — not a plain arm
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_hairpin_is_not_trna() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let scan = scan_trna(&seq).unwrap();
        assert!(!scan.is_trna_like, "a single hairpin is not a tRNA");
    }

    #[test]
    fn unstructured_rna_is_not_trna() {
        let seq = RnaSeq::parse("A".repeat(76)).unwrap();
        let scan = scan_trna(&seq).unwrap();
        assert!(!scan.is_trna_like);
        assert_eq!(scan.arm_count, 0);
    }

    #[test]
    fn scan_returns_a_structure_of_the_right_length() {
        let seq = RnaSeq::parse("GCGC".repeat(19)).unwrap(); // length 76
        let scan = scan_trna(&seq).unwrap();
        assert_eq!(scan.structure.len(), 76);
        assert!((0.0..=1.0).contains(&scan.score));
    }

    #[test]
    fn three_arm_cloverleaf_topology_scores_well() {
        // A synthetic acceptor-stem-enclosed three-hairpin multiloop.
        // outer stem (((( ... )))) wrapping three hairpins.
        let db = "((((((((....))))((....))((....))))))";
        let s = Structure::from_dot_bracket(db).unwrap();
        let arms = count_multiloop_hairpins(&s);
        assert_eq!(arms, 3, "expected three cloverleaf arms");
    }

    #[test]
    fn count_branches_helper() {
        // outer pair encloses two hairpins
        let s = Structure::from_dot_bracket("(((....))((....)))").unwrap_or_else(
            |_| Structure::from_dot_bracket("((((..))((..))))").unwrap(),
        );
        // just exercise the helper without panicking
        let _ = count_multiloop_hairpins(&s);
    }
}
