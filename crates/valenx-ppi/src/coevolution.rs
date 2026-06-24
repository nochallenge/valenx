//! Coevolution / mutual-information contact prediction over a **paired**
//! multiple-sequence alignment.
//!
//! ## The idea
//!
//! Two residues that touch across a protein-protein interface tend to
//! *coevolve*: a substitution on one side is compensated by a
//! substitution on the other, so the two alignment columns carry
//! correlated residue identities. Mutual information (MI) between a
//! column in chain A and a column in chain B quantifies that
//! correlation; the highest-MI inter-chain column pairs are the
//! predicted interface contacts. This is the classic
//! correlated-mutation / direct-coupling family of methods, here in its
//! simplest MI form with the standard **average-product correction
//! (APC)** that removes per-column entropy and phylogenetic background.
//!
//! ## Paired MSA layout
//!
//! A [`PairedMsa`] is two [`Msa`]s — one per
//! chain — sharing the **same number of rows in the same order**: row
//! `k` of each is the orthologue pair from organism `k`. Columns of the
//! A-half are indexed `0..width_a`; columns of the B-half `0..width_b`.
//! Only **inter-chain** pairs `(i in A, j in B)` are scored — that is
//! what an interface contact is.
//!
//! ## Honest scope
//!
//! Plain MI (even APC-corrected) does **not** disentangle direct from
//! transitive (chained) couplings the way a full DCA / pseudolikelihood
//! model (plmDCA, GREMLIN, EVcomplex) does, and it needs a deep, well
//! paired alignment to be meaningful at all. Treat the ranking as a
//! research heuristic, never a validated contact.

use std::collections::HashMap;

use valenx_align::msa::Msa;

use crate::error::PpiError;

/// Smallest paired-MSA depth at which MI is computed at all. Below this,
/// column statistics are pure noise; we fail loud rather than emit a
/// meaningless score. (A genuine analysis wants hundreds of sequences;
/// this floor only rejects the degenerate case.)
pub const MIN_PAIRED_DEPTH: usize = 3;

/// A paired multiple-sequence alignment: aligned orthologues of chain A
/// and of chain B, one row per organism, **in the same row order**.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PairedMsa {
    /// Aligned rows of chain A (the host / first partner).
    pub a: Msa,
    /// Aligned rows of chain B (the pathogen / second partner).
    pub b: Msa,
}

impl PairedMsa {
    /// Build a paired MSA, validating equal depth and non-empty,
    /// equal-length rows on each side.
    ///
    /// # Errors
    /// - [`PpiError::DepthMismatch`] if the two halves differ in depth.
    /// - [`PpiError::TooFewSequences`] if depth `< `[`MIN_PAIRED_DEPTH`].
    /// - [`PpiError::EmptyAlignment`] if either side has zero columns.
    /// - [`PpiError::RaggedRows`] if a side's rows differ in length.
    pub fn new(a: Msa, b: Msa) -> Result<Self, PpiError> {
        if a.depth() != b.depth() {
            return Err(PpiError::DepthMismatch {
                a: a.depth(),
                b: b.depth(),
            });
        }
        if a.depth() < MIN_PAIRED_DEPTH {
            return Err(PpiError::TooFewSequences {
                got: a.depth(),
                need: MIN_PAIRED_DEPTH,
            });
        }
        check_rectangular(&a)?;
        check_rectangular(&b)?;
        if a.width() == 0 || b.width() == 0 {
            return Err(PpiError::EmptyAlignment);
        }
        Ok(PairedMsa { a, b })
    }

    /// Number of paired organisms (rows on each side).
    pub fn depth(&self) -> usize {
        self.a.depth()
    }
}

/// Verify every row of an MSA has the same length, returning a
/// fail-loud error otherwise. (`Msa::new` already enforces this, but a
/// caller can construct an `Msa` directly, so we re-check.)
fn check_rectangular(m: &Msa) -> Result<(), PpiError> {
    let w = m.width();
    for r in &m.rows {
        if r.len() != w {
            return Err(PpiError::RaggedRows {
                got: r.len(),
                expected: w,
            });
        }
    }
    Ok(())
}

/// One predicted inter-chain interface contact: a column of chain A, a
/// column of chain B, and both the APC-corrected and the raw
/// mutual-information scores.
///
/// `score` (APC-corrected MI, the MIp statistic) is the **primary
/// ranking key** — it removes the per-column entropy / phylogenetic
/// background that dominates raw MI on real deep alignments. `raw_mi`
/// is kept visible because (a) it is the honest absolute measure of how
/// much two columns covary, and (b) APC can fully cancel a genuine
/// signal on a degenerate, near-rank-1 MI matrix (a known property on
/// shallow alignments), in which case raw MI is the tie-breaker that
/// keeps a truly coupled pair ahead of an uncoupled one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactPrediction {
    /// Zero-based column index in chain A's alignment.
    pub col_a: usize,
    /// Zero-based column index in chain B's alignment.
    pub col_b: usize,
    /// APC-corrected mutual information (bits) — the primary ranking
    /// key. May be slightly negative after APC subtraction.
    pub score: f64,
    /// Raw (uncorrected) mutual information (bits), always `>= 0`. The
    /// tie-breaker and the basis of the aggregate [`signal`].
    ///
    /// [`signal`]: CoevolutionResult::signal
    pub raw_mi: f64,
}

/// The full inter-chain coevolution analysis: every `(A-column,
/// B-column)` pair ranked by APC-corrected MI, plus a single aggregate
/// signal.
#[derive(Clone, Debug, PartialEq)]
pub struct CoevolutionResult {
    /// All inter-chain column pairs, sorted by descending corrected MI.
    /// Length is `width_a * width_b`.
    pub ranked: Vec<ContactPrediction>,
    /// Number of columns in chain A.
    pub width_a: usize,
    /// Number of columns in chain B.
    pub width_b: usize,
}

impl CoevolutionResult {
    /// The top-`k` predicted contacts (by corrected MI). Fewer are
    /// returned if the interface has fewer than `k` column pairs.
    pub fn top(&self, k: usize) -> &[ContactPrediction] {
        let k = k.min(self.ranked.len());
        &self.ranked[..k]
    }

    /// A single `[0, 1]` coevolution signal: the mean **raw** MI of the
    /// top `L/5` inter-chain pairs (the conventional contact budget,
    /// `L = min(width_a, width_b)`), squashed through `1 - exp(-x)` so a
    /// stronger top-of-list maps monotonically toward `1`. `0` when
    /// there is no coevolution. This is the value fed to the aggregate
    /// [`PpiScore`](crate::score::PpiScore).
    ///
    /// Raw (not APC-corrected) MI is used **deliberately**: APC is a
    /// background-subtraction that sharpens the *relative ranking* of
    /// contacts, but on a shallow or low-rank alignment it can drive the
    /// corrected value of even a genuinely coupled pair to ~0 (it sums
    /// to a near-rank-1 matrix that APC cancels). The aggregate
    /// confidence must reflect *how much covariation actually exists*,
    /// which is what raw MI measures; the ranking inside [`ranked`] /
    /// [`top`] still uses the APC-corrected `score`.
    ///
    /// [`ranked`]: CoevolutionResult::ranked
    /// [`top`]: CoevolutionResult::top
    pub fn signal(&self) -> f64 {
        let l = self.width_a.min(self.width_b);
        let budget = (l / 5).max(1);
        let top = self.top(budget);
        if top.is_empty() {
            return 0.0;
        }
        let mean = top.iter().map(|c| c.raw_mi.max(0.0)).sum::<f64>() / top.len() as f64;
        // Monotone squash to [0, 1); MI in bits, so ~1 bit -> ~0.63.
        1.0 - (-mean).exp()
    }
}

/// Mutual information (bits) between two equal-length, gap-aware
/// alignment columns.
///
/// Gaps (`b'-'`) are treated as a distinct symbol so a gap that tracks a
/// gap on the other side still registers as shared information (the
/// conservative choice; an alternative is gap-masking, documented as a
/// future option). Residues are upper-cased. Returns `0.0` for a column
/// pair with no variation (MI of a constant is zero).
pub fn column_pair_mi(col_a: &[u8], col_b: &[u8]) -> f64 {
    debug_assert_eq!(col_a.len(), col_b.len());
    let n = col_a.len() as f64;
    if n == 0.0 {
        return 0.0;
    }
    let mut joint: HashMap<(u8, u8), f64> = HashMap::new();
    let mut pa: HashMap<u8, f64> = HashMap::new();
    let mut pb: HashMap<u8, f64> = HashMap::new();
    for (&x, &y) in col_a.iter().zip(col_b) {
        let x = x.to_ascii_uppercase();
        let y = y.to_ascii_uppercase();
        *joint.entry((x, y)).or_insert(0.0) += 1.0;
        *pa.entry(x).or_insert(0.0) += 1.0;
        *pb.entry(y).or_insert(0.0) += 1.0;
    }
    let mut mi = 0.0;
    for (&(x, y), &nxy) in &joint {
        let pxy = nxy / n;
        let px = pa[&x] / n;
        let py = pb[&y] / n;
        mi += pxy * (pxy / (px * py)).log2();
    }
    mi.max(0.0)
}

/// Collect the bytes of column `c` down every row of `m`.
fn column_bytes(m: &Msa, c: usize) -> Vec<u8> {
    m.rows.iter().map(|r| r[c]).collect()
}

/// Predict inter-chain interface contacts from a paired MSA by
/// APC-corrected mutual information.
///
/// For every column `i` of chain A and `j` of chain B it computes raw MI
/// `MI(i, j)`, then subtracts the **average product correction**
///
/// ```text
///   APC(i, j) = ( mean_j MI(i, ·) * mean_i MI(·, j) ) / mean MI
/// ```
///
/// over the inter-chain MI matrix. APC removes the per-column entropy /
/// phylogenetic background that otherwise dominates raw MI, and is the
/// standard correction in the MIp / coevolution literature
/// (Dunn et al. 2008). The corrected pairs are returned sorted by
/// descending score.
///
/// # Errors
/// Propagates [`PairedMsa`] validation errors via the constructor; this
/// function additionally never panics on a well-formed `PairedMsa`.
pub fn predict_contacts(paired: &PairedMsa) -> Result<CoevolutionResult, PpiError> {
    let wa = paired.a.width();
    let wb = paired.b.width();
    // PairedMsa::new already guarantees wa>0, wb>0, equal depth, but a
    // PairedMsa can be field-constructed; re-validate fail-loud.
    if wa == 0 || wb == 0 {
        return Err(PpiError::EmptyAlignment);
    }

    // Precompute columns once.
    let cols_a: Vec<Vec<u8>> = (0..wa).map(|i| column_bytes(&paired.a, i)).collect();
    let cols_b: Vec<Vec<u8>> = (0..wb).map(|j| column_bytes(&paired.b, j)).collect();

    // Raw inter-chain MI matrix.
    let mut raw = vec![vec![0.0f64; wb]; wa];
    let mut sum = 0.0;
    for (i, ca) in cols_a.iter().enumerate() {
        for (j, cb) in cols_b.iter().enumerate() {
            let mi = column_pair_mi(ca, cb);
            raw[i][j] = mi;
            sum += mi;
        }
    }

    // Row / column means and overall mean for APC.
    let total = (wa * wb) as f64;
    let overall = if total > 0.0 { sum / total } else { 0.0 };
    let row_mean: Vec<f64> = raw
        .iter()
        .map(|row| row.iter().sum::<f64>() / wb as f64)
        .collect();
    let mut col_mean = vec![0.0f64; wb];
    for row in &raw {
        for (j, &v) in row.iter().enumerate() {
            col_mean[j] += v;
        }
    }
    for v in &mut col_mean {
        *v /= wa as f64;
    }

    // Corrected score = MI - APC. If the whole matrix is zero (no
    // variation anywhere), APC is zero and every score is zero.
    let mut ranked: Vec<ContactPrediction> = Vec::with_capacity(wa * wb);
    for (i, (row, &rm)) in raw.iter().zip(&row_mean).enumerate() {
        for (j, (&mi, &cm)) in row.iter().zip(&col_mean).enumerate() {
            let apc = if overall > 0.0 {
                rm * cm / overall
            } else {
                0.0
            };
            ranked.push(ContactPrediction {
                col_a: i,
                col_b: j,
                score: mi - apc,
                raw_mi: mi,
            });
        }
    }

    // Deterministic sort: descending APC-corrected score (primary),
    // then descending raw MI (so a genuinely coupled pair wins when APC
    // ties it to zero on a degenerate matrix), then (col_a, col_b) so
    // the order never depends on Vec insertion races.
    ranked.sort_by(|p, q| {
        q.score
            .partial_cmp(&p.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(
                q.raw_mi
                    .partial_cmp(&p.raw_mi)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
            .then(p.col_a.cmp(&q.col_a))
            .then(p.col_b.cmp(&q.col_b))
    });

    Ok(CoevolutionResult {
        ranked,
        width_a: wa,
        width_b: wb,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn msa(rows: &[&[u8]]) -> Msa {
        Msa::new(rows.iter().map(|r| r.to_vec()).collect()).unwrap()
    }

    #[test]
    fn mi_of_constant_columns_is_zero() {
        assert!(column_pair_mi(b"AAAA", b"CCCC").abs() < 1e-12);
        assert!(column_pair_mi(b"AAAA", b"ACGT").abs() < 1e-12);
    }

    #[test]
    fn mi_of_perfectly_correlated_columns_is_one_bit() {
        // Two states, perfectly coupled: A<->C, T<->G. 1 bit of MI.
        let mi = column_pair_mi(b"AATT", b"CCGG");
        assert!((mi - 1.0).abs() < 1e-9, "mi = {mi}");
    }

    #[test]
    fn paired_msa_rejects_depth_mismatch() {
        let a = msa(&[b"AA", b"AA", b"AA"]);
        let b = msa(&[b"CC", b"CC"]);
        let err = PairedMsa::new(a, b).unwrap_err();
        assert_eq!(err.code(), "depth_mismatch");
    }

    #[test]
    fn paired_msa_rejects_too_few_sequences() {
        let a = msa(&[b"AA", b"AA"]);
        let b = msa(&[b"CC", b"CC"]);
        let err = PairedMsa::new(a, b).unwrap_err();
        assert_eq!(err.code(), "too_few_sequences");
    }

    #[test]
    fn coevolving_pair_ranks_first() {
        // Chain A col 1 perfectly tracks chain B col 0; all other
        // columns are conserved (zero MI). The coupled pair must top
        // the ranking.
        let a = msa(&[b"MA", b"MA", b"MT", b"MT"]); // col1 varies A/T
        let b = msa(&[b"CK", b"CK", b"GK", b"GK"]); // col0 varies C/G, tracks A-col1
        let paired = PairedMsa::new(a, b).unwrap();
        let res = predict_contacts(&paired).unwrap();
        let top = res.ranked[0];
        // The coupled pair wins the ranking. Its APC-corrected `score`
        // is ~0 here (APC fully cancels a rank-1 MI matrix — a known
        // property), but its raw MI is a full bit and breaks the tie, so
        // it still ranks first and drives a positive aggregate signal.
        assert_eq!((top.col_a, top.col_b), (1, 0));
        assert!((top.raw_mi - 1.0).abs() < 1e-9, "raw_mi = {}", top.raw_mi);
        assert!(res.signal() > 0.0);
    }
}
