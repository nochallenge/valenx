//! Progressive multiple-sequence alignment (Clustal-class).
//!
//! Progressive alignment is the classic MSA heuristic: build a guide
//! tree, then walk it bottom-up, at each internal node aligning the
//! two child sub-alignments with profile-profile DP. Closely related
//! sequences are merged first, so an early-introduced gap is rarely
//! wrong ("once a gap, always a gap" — the known progressive-MSA
//! limitation that the [`mod@crate::msa::refine`] module then
//! mitigates).
//!
//! [`align`] is the one-call entry point: sequences in → an
//! [`Msa`] out. It runs [`distance_matrix`] →
//! [`upgma`] → tree walk with [`align_profiles`].

use super::guidetree::{distance_matrix, upgma, GuideTree, TreeNode};
use super::profile::{align_profiles, Profile};
use crate::error::{AlignError, Result};
use crate::matrix::ScoringScheme;

/// A finished multiple-sequence alignment: equal-length gapped rows in
/// the original input order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Msa {
    /// One gapped row per input sequence, all the same length, in the
    /// order the sequences were supplied.
    pub rows: Vec<Vec<u8>>,
}

impl Msa {
    /// Wraps a set of rows, validating equal length.
    pub fn new(rows: Vec<Vec<u8>>) -> Result<Self> {
        if let Some(first) = rows.first() {
            let w = first.len();
            for r in &rows {
                if r.len() != w {
                    return Err(AlignError::dimension(format!(
                        "MSA rows differ: {} vs {w}",
                        r.len()
                    )));
                }
            }
        }
        Ok(Msa { rows })
    }

    /// Number of aligned sequences.
    pub fn depth(&self) -> usize {
        self.rows.len()
    }

    /// Number of alignment columns (`0` for an empty MSA).
    pub fn width(&self) -> usize {
        self.rows.first().map(Vec::len).unwrap_or(0)
    }

    /// `true` if the MSA has no sequences.
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Row `i` as a `&str`.
    pub fn row_str(&self, i: usize) -> Option<&str> {
        self.rows
            .get(i)
            .map(|r| std::str::from_utf8(r).unwrap_or("<non-utf8>"))
    }

    /// The sum-of-pairs score: the total pairwise score over every
    /// distinct row pair and every column, under `scheme`. The
    /// standard MSA objective — higher is better.
    pub fn sum_of_pairs(&self, scheme: &ScoringScheme) -> i32 {
        let mut total = 0i64;
        let depth = self.depth();
        for i in 0..depth {
            for j in (i + 1)..depth {
                total += pair_score(&self.rows[i], &self.rows[j], scheme) as i64;
            }
        }
        total as i32
    }

    /// Builds a [`Profile`] over the alignment.
    pub fn to_profile(&self) -> Result<Profile> {
        Profile::from_alignment(&self.rows)
    }
}

/// Score of one pair of gapped rows: substitution scores for aligned
/// residues, `gap.extend` penalty per gap column (gap/gap free).
fn pair_score(a: &[u8], b: &[u8], scheme: &ScoringScheme) -> i32 {
    let mut s = 0;
    for (&x, &y) in a.iter().zip(b) {
        match (x == b'-', y == b'-') {
            (false, false) => s += scheme.sub(x, y),
            (true, true) => {}
            _ => s -= scheme.gap.extend,
        }
    }
    s
}

/// Progressive MSA of `seqs` under `scheme`.
///
/// Builds a UPGMA guide tree from pairwise identities, then merges
/// profiles bottom-up. Returns [`AlignError::Invalid`] for an empty
/// input; a single sequence is returned unchanged.
pub fn align(seqs: &[&[u8]], scheme: &ScoringScheme) -> Result<Msa> {
    if seqs.is_empty() {
        return Err(AlignError::invalid("seqs", "MSA needs >= 1 sequence"));
    }
    if seqs.len() == 1 {
        return Msa::new(vec![seqs[0].to_vec()]);
    }
    let dm = distance_matrix(seqs, scheme)?;
    let tree = upgma(&dm)?;
    align_along_tree(seqs, &tree, scheme)
}

/// Progressive MSA using a *caller-supplied* guide tree — used by the
/// iterative-refinement pass, which re-aligns sub-trees.
pub fn align_along_tree(seqs: &[&[u8]], tree: &GuideTree, scheme: &ScoringScheme) -> Result<Msa> {
    let root = tree
        .root()
        .ok_or_else(|| AlignError::invalid("tree", "empty guide tree"))?;

    // Each node yields a (profile, original-sequence-index list). The
    // profile's rows are in the order of that index list.
    let (profile, order) = build_node(root, tree, seqs, scheme)?;

    // Reorder the profile rows back to the original input order.
    let mut rows: Vec<Vec<u8>> = vec![Vec::new(); seqs.len()];
    for (slot, &seq_idx) in order.iter().enumerate() {
        rows[seq_idx] = profile.rows[slot].clone();
    }
    Msa::new(rows)
}

/// Recursively aligns a guide-tree node, returning its profile and the
/// input-sequence indices its rows correspond to (in row order).
fn build_node(
    idx: usize,
    tree: &GuideTree,
    seqs: &[&[u8]],
    scheme: &ScoringScheme,
) -> Result<(Profile, Vec<usize>)> {
    match &tree.nodes[idx] {
        TreeNode::Leaf(s) => {
            let p = Profile::from_sequence(seqs[*s])?;
            Ok((p, vec![*s]))
        }
        TreeNode::Internal { left, right, .. } => {
            let (lp, lorder) = build_node(*left, tree, seqs, scheme)?;
            let (rp, rorder) = build_node(*right, tree, seqs, scheme)?;
            let merged = align_profiles(&lp, &rp, scheme)?;
            // align_profiles stacks A's rows then B's rows.
            let mut order = lorder;
            order.extend(rorder);
            let p = Profile::from_alignment(&merged.rows)?;
            Ok((p, order))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};

    fn dna() -> ScoringScheme {
        ScoringScheme::new(SubstitutionMatrix::dna_simple(2, -1), GapCost::new(4, 1))
    }

    #[test]
    fn single_sequence_passthrough() {
        let msa = align(&[b"ACGTACGT"], &dna()).unwrap();
        assert_eq!(msa.depth(), 1);
        assert_eq!(msa.rows[0], b"ACGTACGT");
    }

    #[test]
    fn identical_sequences_align_without_gaps() {
        let seqs: &[&[u8]] = &[b"ACGTACGT", b"ACGTACGT", b"ACGTACGT"];
        let msa = align(seqs, &dna()).unwrap();
        assert_eq!(msa.depth(), 3);
        assert_eq!(msa.width(), 8);
        for r in &msa.rows {
            assert_eq!(r, b"ACGTACGT");
            assert!(!r.contains(&b'-'));
        }
    }

    #[test]
    fn rows_stay_equal_length_and_recover_inputs() {
        let seqs: &[&[u8]] = &[b"ACGTACGT", b"ACGTCGT", b"ACGTACGTT", b"ACGACGT"];
        let msa = align(seqs, &dna()).unwrap();
        let w = msa.width();
        assert!(msa.rows.iter().all(|r| r.len() == w));
        // Stripping gaps recovers each original, in order.
        for (i, &orig) in seqs.iter().enumerate() {
            let stripped: Vec<u8> = msa.rows[i].iter().copied().filter(|&c| c != b'-').collect();
            assert_eq!(stripped, orig, "row {i} should recover its input");
        }
    }

    #[test]
    fn gap_introduced_for_indel() {
        // One sequence is missing a residue; the MSA must gap it.
        let seqs: &[&[u8]] = &[b"ACGTACGT", b"ACGTACGT", b"ACGTCGT"];
        let msa = align(seqs, &dna()).unwrap();
        // The short row carries exactly one gap.
        let gapped = msa.rows.iter().find(|r| r.contains(&b'-')).unwrap();
        assert_eq!(gapped.iter().filter(|&&c| c == b'-').count(), 1);
    }

    #[test]
    fn sum_of_pairs_positive_for_similar() {
        let seqs: &[&[u8]] = &[b"ACGTACGT", b"ACGTACGT", b"ACGTACGA"];
        let msa = align(seqs, &dna()).unwrap();
        assert!(msa.sum_of_pairs(&dna()) > 0);
    }

    #[test]
    fn protein_msa() {
        let scheme = ScoringScheme::new(SubstitutionMatrix::blosum62(), GapCost::new(11, 1));
        let seqs: &[&[u8]] = &[b"MKVLAAGG", b"MKVLAAGG", b"MKVLAGG"];
        let msa = align(seqs, &scheme).unwrap();
        assert_eq!(msa.depth(), 3);
        let w = msa.width();
        assert!(msa.rows.iter().all(|r| r.len() == w));
    }

    #[test]
    fn empty_input_rejected() {
        assert!(align(&[], &dna()).is_err());
    }
}
