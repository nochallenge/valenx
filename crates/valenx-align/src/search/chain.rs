//! Anchor chaining — colinear DP chaining of seed anchors.
//!
//! After seeding (k-mers or minimizers) a long-read mapper holds a
//! cloud of [`Anchor`]s — short exact matches, each a
//! `(query-position, target-position)` pair. The real homology is the
//! longest *colinear* run of anchors: a subset that increases in both
//! coordinates and stays close to one diagonal. [`chain_anchors`]
//! finds the maximum-score such chain with the minimap2-style
//! dynamic-programming recurrence
//!
//! ```text
//! f(i) = anchor_weight(i) + max( 0, max_{j<i, j->i colinear} f(j) − gap_cost(j, i) )
//! ```
//!
//! where the gap cost penalises the difference between the query and
//! target offsets of consecutive anchors (indels) and the raw
//! distance jumped.

/// A seed anchor — one short exact match shared by query and target.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct Anchor {
    /// 0-based start of the match in the query.
    pub query_pos: usize,
    /// 0-based start of the match in the target.
    pub target_pos: usize,
    /// Length of the exact match (e.g. the k-mer length).
    pub len: usize,
}

impl Anchor {
    /// Builds an anchor.
    pub fn new(query_pos: usize, target_pos: usize, len: usize) -> Self {
        Anchor {
            query_pos,
            target_pos,
            len,
        }
    }

    /// The anchor's diagonal `target_pos − query_pos` (may be negative).
    pub fn diagonal(&self) -> isize {
        self.target_pos as isize - self.query_pos as isize
    }

    /// Query end coordinate (exclusive).
    pub fn query_end(&self) -> usize {
        self.query_pos + self.len
    }

    /// Target end coordinate (exclusive).
    pub fn target_end(&self) -> usize {
        self.target_pos + self.len
    }
}

/// Tunable parameters for [`chain_anchors`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ChainParams {
    /// Penalty weight on the indel (diagonal difference) between two
    /// chained anchors. Higher = stricter colinearity.
    pub gap_weight: f64,
    /// Maximum allowed gap, in residues, between consecutive anchors
    /// on either axis. Anchors farther apart cannot chain.
    pub max_gap: usize,
}

impl Default for ChainParams {
    fn default() -> Self {
        ChainParams {
            gap_weight: 0.5,
            max_gap: 10_000,
        }
    }
}

/// A colinear chain of anchors and its score.
#[derive(Clone, Debug, PartialEq)]
pub struct Chain {
    /// The chained anchors, ordered by increasing coordinate.
    pub anchors: Vec<Anchor>,
    /// The chaining score (sum of anchor weights minus gap penalties).
    pub score: f64,
}

impl Chain {
    /// Half-open query span `[start, end)` covered by the chain.
    /// `(0, 0)` for an empty chain.
    pub fn query_span(&self) -> (usize, usize) {
        match (self.anchors.first(), self.anchors.last()) {
            (Some(f), Some(l)) => (f.query_pos, l.query_end()),
            _ => (0, 0),
        }
    }

    /// Half-open target span `[start, end)` covered by the chain.
    pub fn target_span(&self) -> (usize, usize) {
        match (self.anchors.first(), self.anchors.last()) {
            (Some(f), Some(l)) => (f.target_pos, l.target_end()),
            _ => (0, 0),
        }
    }

    /// Number of anchors in the chain.
    pub fn len(&self) -> usize {
        self.anchors.len()
    }

    /// `true` if the chain has no anchors.
    pub fn is_empty(&self) -> bool {
        self.anchors.is_empty()
    }
}

/// Finds the single highest-scoring colinear chain among `anchors`.
///
/// Anchors are sorted by `(query_pos, target_pos)` and an
/// O(n²) DP picks the best predecessor for each. Returns an empty
/// [`Chain`] (score `0`) when `anchors` is empty.
pub fn chain_anchors(anchors: &[Anchor], params: ChainParams) -> Chain {
    if anchors.is_empty() {
        return Chain {
            anchors: Vec::new(),
            score: 0.0,
        };
    }

    // Sort by query then target so a valid chain is a subsequence.
    let mut sorted: Vec<Anchor> = anchors.to_vec();
    sorted.sort_by(|a, b| {
        a.query_pos
            .cmp(&b.query_pos)
            .then(a.target_pos.cmp(&b.target_pos))
    });
    let n = sorted.len();

    let mut f = vec![0f64; n]; // best chain score ending at i
    let mut prev = vec![usize::MAX; n]; // predecessor index

    for i in 0..n {
        f[i] = sorted[i].len as f64;
        for j in 0..i {
            if !colinear(&sorted[j], &sorted[i], params.max_gap) {
                continue;
            }
            let gap = gap_cost(&sorted[j], &sorted[i], params.gap_weight);
            let cand = f[j] + sorted[i].len as f64 - gap;
            if cand > f[i] {
                f[i] = cand;
                prev[i] = j;
            }
        }
    }

    // Best endpoint.
    let mut best_i = 0;
    for i in 1..n {
        if f[i] > f[best_i] {
            best_i = i;
        }
    }

    // Backtrack.
    let mut chain = Vec::new();
    let mut idx = best_i;
    while idx != usize::MAX {
        chain.push(sorted[idx]);
        idx = prev[idx];
    }
    chain.reverse();

    Chain {
        anchors: chain,
        score: f[best_i],
    }
}

/// `true` if anchor `b` can follow anchor `a` in a colinear chain:
/// `b` must start strictly after `a` on *both* axes and the jump on
/// each axis must not exceed `max_gap`.
fn colinear(a: &Anchor, b: &Anchor, max_gap: usize) -> bool {
    if b.query_pos < a.query_end() || b.target_pos < a.target_end() {
        return false; // overlap or out of order
    }
    let qgap = b.query_pos - a.query_end();
    let tgap = b.target_pos - a.target_end();
    qgap <= max_gap && tgap <= max_gap
}

/// The gap penalty between two consecutive chained anchors: the indel
/// (difference of axis gaps) scaled by `gap_weight`.
fn gap_cost(a: &Anchor, b: &Anchor, gap_weight: f64) -> f64 {
    let qgap = (b.query_pos - a.query_end()) as f64;
    let tgap = (b.target_pos - a.target_end()) as f64;
    let indel = (qgap - tgap).abs();
    gap_weight * indel
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_input() {
        let c = chain_anchors(&[], ChainParams::default());
        assert!(c.is_empty());
        assert_eq!(c.score, 0.0);
    }

    #[test]
    fn perfectly_colinear_anchors_all_chain() {
        // Five anchors on one diagonal, evenly spaced — all chain.
        let anchors: Vec<Anchor> = (0..5).map(|i| Anchor::new(i * 20, i * 20, 10)).collect();
        let c = chain_anchors(&anchors, ChainParams::default());
        assert_eq!(c.len(), 5);
        assert_eq!(c.query_span(), (0, 90));
    }

    #[test]
    fn off_diagonal_outlier_excluded() {
        // Four colinear anchors plus one far off-diagonal noise anchor.
        let mut anchors: Vec<Anchor> = (0..4).map(|i| Anchor::new(i * 20, i * 20, 10)).collect();
        anchors.push(Anchor::new(35, 5000, 10)); // wild outlier
        let c = chain_anchors(&anchors, ChainParams::default());
        // The outlier should not be part of the best chain.
        assert!(!c.anchors.iter().any(|a| a.target_pos == 5000));
        assert_eq!(c.len(), 4);
    }

    #[test]
    fn anti_diagonal_anchors_do_not_chain() {
        // Anchors that decrease in target as query increases: no chain
        // longer than 1.
        let anchors = vec![
            Anchor::new(0, 100, 10),
            Anchor::new(20, 70, 10),
            Anchor::new(40, 40, 10),
        ];
        let c = chain_anchors(&anchors, ChainParams::default());
        assert_eq!(c.len(), 1, "non-colinear anchors cannot chain");
    }

    #[test]
    fn small_indel_still_chains() {
        // A 5-residue indel between two anchors: tolerated by chaining.
        let anchors = vec![
            Anchor::new(0, 0, 10),
            Anchor::new(20, 25, 10), // target 5 ahead => indel of 5
        ];
        let c = chain_anchors(&anchors, ChainParams::default());
        assert_eq!(c.len(), 2);
    }

    #[test]
    fn max_gap_breaks_chain() {
        // A jump larger than max_gap cannot be bridged.
        let params = ChainParams {
            gap_weight: 0.5,
            max_gap: 50,
        };
        let anchors = vec![
            Anchor::new(0, 0, 10),
            Anchor::new(1000, 1000, 10), // 990 gap > max_gap
        ];
        let c = chain_anchors(&anchors, params);
        assert_eq!(c.len(), 1);
    }
}
