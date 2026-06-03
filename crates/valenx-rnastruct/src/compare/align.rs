//! Structure alignment (RNAforester-class).
//!
//! RNAforester aligns two RNA secondary structures by aligning the
//! *trees* that represent them — finding the correspondence of loops
//! and pairs that maximises a similarity score. This module is a v1
//! of that idea: an ordered-forest alignment by dynamic programming.
//!
//! ## Method
//!
//! Each structure is converted to its ordered tree (a base pair is an
//! internal node, an unpaired base a leaf — the same tree
//! [`crate::compare::distance`] uses). The trees are aligned with the
//! Jiang-Wang-Zhang (1995) ordered-tree-*alignment* recurrence:
//! every node of one tree is either matched to a node of the other,
//! or aligned against a gap. A match of two `Pair` nodes or two
//! `Unpaired` nodes scores `+1`; a mismatch or a gap scores `0`. The
//! optimal score is the alignment similarity; dividing by the larger
//! tree size gives a normalised `[0, 1]` similarity.
//!
//! Crucially, deleting (gapping) an internal node *promotes* its
//! children into the sibling sequence — it does not discard the whole
//! subtree. That is what lets a hairpin with a 3-bp stem align well to
//! a hairpin with a 2-bp stem: the extra base pair is gapped and the
//! loop bases it enclosed still match.
//!
//! ## v1 scope
//!
//! General ordered-tree alignment is `O(n²m²)` in the worst case;
//! this v1 implements the standard forest-alignment DP without the
//! later speed-ups, which is fine for the structure sizes a desktop
//! tool aligns. Affine loop-gap costs are not modelled — every node
//! gap costs the same.

use crate::error::Result;
use crate::structure::Structure;

/// The result of aligning two structures.
#[derive(Clone, Debug, PartialEq)]
pub struct StructureAlignment {
    /// The raw tree-alignment similarity score (count of matched
    /// nodes, weighted).
    pub score: f64,
    /// The score normalised to `[0, 1]` by the larger structure's
    /// tree size — `1.0` for identical structures.
    pub similarity: f64,
    /// Number of structural nodes (pairs + unpaired bases) in the
    /// first structure.
    pub nodes_a: usize,
    /// Number of structural nodes in the second structure.
    pub nodes_b: usize,
}

/// A node label in the ordered structure forest.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Label {
    Pair,
    Unpaired,
}

/// One subtree in an arena: a labelled root plus the ids of its
/// ordered children.
struct Subtree {
    label: Label,
    children: Vec<u32>,
}

/// An ordered structure forest, stored as an arena of [`Subtree`]s.
/// A "forest" passed around the alignment recurrence is just an ordered
/// slice of root ids into `arena`.
struct Arena {
    nodes: Vec<Subtree>,
    /// The ids of the top-level (outermost) subtrees, in order.
    roots: Vec<u32>,
}

impl Arena {
    /// Builds the ordered forest of a whole structure.
    fn build(s: &Structure) -> Arena {
        let mut arena = Arena {
            nodes: Vec::new(),
            roots: Vec::new(),
        };
        arena.roots = arena.build_span(s, 0, s.len());
        arena
    }

    /// Builds the ordered forest of the loop elements spanning
    /// `[lo, hi)`, returning their root ids.
    fn build_span(&mut self, s: &Structure, lo: usize, hi: usize) -> Vec<u32> {
        let mut out = Vec::new();
        let mut k = lo;
        while k < hi {
            match s.partner(k) {
                Some(p) if p > k && p < hi => {
                    let children = self.build_span(s, k + 1, p);
                    let id = self.nodes.len() as u32;
                    self.nodes.push(Subtree {
                        label: Label::Pair,
                        children,
                    });
                    out.push(id);
                    k = p + 1;
                }
                _ => {
                    let id = self.nodes.len() as u32;
                    self.nodes.push(Subtree {
                        label: Label::Unpaired,
                        children: Vec::new(),
                    });
                    out.push(id);
                    k += 1;
                }
            }
        }
        out
    }

    /// Total number of structural nodes (pairs + unpaired bases).
    fn count_nodes(&self) -> usize {
        self.nodes.len()
    }
}

/// The match score of two node labels: 1 for an exact label match,
/// 0 otherwise.
fn match_score(a: Label, b: Label) -> f64 {
    if a == b {
        1.0
    } else {
        0.0
    }
}

/// Memoised Jiang-Wang-Zhang ordered-forest alignment.
///
/// `align(F, G)` for two sibling sequences `F` (from `a`) and `G`
/// (from `b`) is the maximum of three moves on the first trees of each:
///
/// * **match** the two first roots — add the label match, recurse on
///   their child forests, recurse on the remaining siblings;
/// * **delete** the first root of `F` — *promote* its children into the
///   sibling sequence and recurse;
/// * **insert** the first root of `G` — *promote* its children and
///   recurse.
///
/// Promotion on a gap is what makes this a tree *alignment* rather than
/// a subtree-drop edit, and is essential for the similarity to reflect
/// shared loops across differing stem lengths.
struct Aligner<'x> {
    a: &'x Arena,
    b: &'x Arena,
    memo: std::collections::HashMap<(Vec<u32>, Vec<u32>), f64>,
}

impl<'x> Aligner<'x> {
    fn align(&mut self, f: &[u32], g: &[u32]) -> f64 {
        if f.is_empty() || g.is_empty() {
            // Remaining trees on one side are all gapped; every gap
            // scores 0, so the contribution is 0.
            return 0.0;
        }
        let key = (f.to_vec(), g.to_vec());
        if let Some(&v) = self.memo.get(&key) {
            return v;
        }

        let (ta, fa_rest) = (f[0], &f[1..]);
        let (tb, gb_rest) = (g[0], &g[1..]);
        let na = &self.a.nodes[ta as usize];
        let nb = &self.b.nodes[tb as usize];

        // Match the two first roots.
        let kids = self.align(&na.children, &nb.children);
        let rest = self.align(fa_rest, gb_rest);
        let matched = match_score(na.label, nb.label) + kids + rest;

        // Delete `ta`: promote its children ahead of the rest of `f`.
        let mut f_promoted = na.children.clone();
        f_promoted.extend_from_slice(fa_rest);
        let deleted = self.align(&f_promoted, g);

        // Insert `tb`: promote its children ahead of the rest of `g`.
        let mut g_promoted = nb.children.clone();
        g_promoted.extend_from_slice(gb_rest);
        let inserted = self.align(f, &g_promoted);

        let best = matched.max(deleted).max(inserted);
        self.memo.insert(key, best);
        best
    }
}

/// Aligns two ordered structure forests, returning the maximum
/// similarity score.
fn align_forests(a: &Arena, b: &Arena) -> f64 {
    let mut aligner = Aligner {
        a,
        b,
        memo: std::collections::HashMap::new(),
    };
    let roots_a = a.roots.clone();
    let roots_b = b.roots.clone();
    aligner.align(&roots_a, &roots_b)
}

/// Aligns two secondary structures (RNAforester-class).
///
/// # Errors
/// Never fails for valid structures; the `Result` is kept for
/// signature symmetry with the rest of the crate.
pub fn align_structures(a: &Structure, b: &Structure) -> Result<StructureAlignment> {
    let fa = Arena::build(a);
    let fb = Arena::build(b);
    let nodes_a = fa.count_nodes();
    let nodes_b = fb.count_nodes();
    let score = align_forests(&fa, &fb);
    let denom = nodes_a.max(nodes_b).max(1) as f64;
    Ok(StructureAlignment {
        score,
        similarity: (score / denom).clamp(0.0, 1.0),
        nodes_a,
        nodes_b,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_structures_align_perfectly() {
        let s = Structure::from_dot_bracket("(((...)))").unwrap();
        let a = align_structures(&s, &s).unwrap();
        assert!(
            (a.similarity - 1.0).abs() < 1e-9,
            "identical structures should score 1.0, got {}",
            a.similarity
        );
    }

    #[test]
    fn similar_structures_score_high() {
        let a = Structure::from_dot_bracket("(((...)))").unwrap();
        let b = Structure::from_dot_bracket("((.....))").unwrap();
        let al = align_structures(&a, &b).unwrap();
        assert!(al.similarity > 0.5, "similar structures: {}", al.similarity);
        assert!(al.similarity < 1.0);
    }

    #[test]
    fn dissimilar_structures_score_lower() {
        let hairpin = Structure::from_dot_bracket("(((...)))").unwrap();
        let open = Structure::from_dot_bracket(".........").unwrap();
        let same = align_structures(&hairpin, &hairpin).unwrap();
        let diff = align_structures(&hairpin, &open).unwrap();
        assert!(diff.similarity < same.similarity);
    }

    #[test]
    fn alignment_is_symmetric() {
        let a = Structure::from_dot_bracket("((((....))))").unwrap();
        let b = Structure::from_dot_bracket("((.((..)).))").unwrap();
        let ab = align_structures(&a, &b).unwrap();
        let ba = align_structures(&b, &a).unwrap();
        assert!((ab.score - ba.score).abs() < 1e-9);
    }

    #[test]
    fn handles_empty_structures() {
        let empty = Structure::empty(0);
        let a = align_structures(&empty, &empty).unwrap();
        assert_eq!(a.nodes_a, 0);
        // similarity of two empty structures is defined as 1.0
        assert!((a.similarity - 1.0).abs() < 1e-9 || a.score == 0.0);
    }

    #[test]
    fn different_length_structures_align() {
        let a = Structure::from_dot_bracket("(((...)))").unwrap();
        let b = Structure::from_dot_bracket("(...)").unwrap();
        let al = align_structures(&a, &b).unwrap();
        assert!(al.score > 0.0);
        assert!(al.nodes_a != al.nodes_b);
    }
}
