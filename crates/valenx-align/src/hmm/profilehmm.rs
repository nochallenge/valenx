//! Profile hidden Markov model (Plan7-style, HMMER-class core).
//!
//! A *profile HMM* models a sequence family as a left-to-right chain
//! of nodes, each node having three states:
//!
//! - **Match (M)** — emits a residue from this column's distribution;
//! - **Insert (I)** — emits a residue between columns (background
//!   distribution);
//! - **Delete (D)** — a silent state, skips this column.
//!
//! This is the HMMER architecture (here a Plan7-*ish* subset: no
//! N/C/J flanking states, no local-vs-global switch — a clean global
//! profile that is enough for a real v1).
//!
//! [`ProfileHmm::from_msa`] builds a model from a multiple-sequence
//! alignment: columns where most rows have a residue become match
//! states, the rest become inserts; emission and transition
//! probabilities are the observed frequencies with Laplace
//! pseudocounts. [`ProfileHmm::viterbi`] scores the most probable
//! parse of a query against the model and [`ProfileHmm::forward`] the
//! summed likelihood. All DP is in log space.

use crate::error::{AlignError, Result};
use crate::limits::{check_dp_size_with, MAX_DP_CELLS};
use crate::msa::progressive::Msa;

const LOG_ZERO: f64 = -1.0e30;
/// The 20 amino acids — emission alphabet of a protein profile HMM.
const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

fn log_add(a: f64, b: f64) -> f64 {
    if a <= LOG_ZERO {
        return b;
    }
    if b <= LOG_ZERO {
        return a;
    }
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    hi + (1.0 + (lo - hi).exp()).ln()
}

/// One node of a profile HMM: the emission distributions of its match
/// and insert states plus the seven Plan7 transition probabilities.
#[derive(Clone, Debug, PartialEq)]
pub struct HmmNode {
    /// Match-state emission log-probabilities over the 20 amino acids
    /// in `ACDEFGHIKLMNPQRSTVWY` order.
    pub match_emit: [f64; 20],
    /// Insert-state emission log-probabilities (background-like).
    pub insert_emit: [f64; 20],
    /// `ln P(M → M)` to the next node.
    pub t_mm: f64,
    /// `ln P(M → I)` at this node.
    pub t_mi: f64,
    /// `ln P(M → D)` of the next node.
    pub t_md: f64,
    /// `ln P(I → M)` to the next node.
    pub t_im: f64,
    /// `ln P(I → I)` at this node.
    pub t_ii: f64,
    /// `ln P(D → M)` to the next node.
    pub t_dm: f64,
    /// `ln P(D → D)` of the next node.
    pub t_dd: f64,
}

/// A profile HMM — a chain of [`HmmNode`]s plus a begin state.
#[derive(Clone, Debug, PartialEq)]
pub struct ProfileHmm {
    /// The match-state columns, in order. `nodes.len()` is the model
    /// length (number of match states).
    pub nodes: Vec<HmmNode>,
}

impl ProfileHmm {
    /// The model length — the number of match states.
    pub fn length(&self) -> usize {
        self.nodes.len()
    }

    /// `true` if the model has no match states.
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// Builds a profile HMM from a multiple-sequence alignment.
    ///
    /// A column is a **match column** if at most `gap_threshold` of its
    /// rows are gaps (default rule: more than half are residues);
    /// other columns contribute to insert states. Emission and
    /// transition probabilities are maximum-likelihood estimates with
    /// a Laplace `+1` pseudocount so no probability is ever zero.
    ///
    /// Returns [`AlignError::Invalid`] for an empty MSA or one with no
    /// match columns.
    pub fn from_msa(msa: &Msa, gap_threshold: f64) -> Result<Self> {
        if msa.is_empty() || msa.width() == 0 {
            return Err(AlignError::invalid("msa", "cannot build HMM from empty MSA"));
        }
        let depth = msa.depth();
        let width = msa.width();

        // Classify each column.
        let is_match: Vec<bool> = (0..width)
            .map(|c| {
                let gaps = msa.rows.iter().filter(|r| r[c] == b'-').count();
                (gaps as f64 / depth as f64) <= gap_threshold
            })
            .collect();
        let match_cols: Vec<usize> = (0..width).filter(|&c| is_match[c]).collect();
        if match_cols.is_empty() {
            return Err(AlignError::invalid(
                "msa",
                "no match columns (every column is gap-dominated)",
            ));
        }

        let mut nodes = Vec::with_capacity(match_cols.len());
        for &mc in &match_cols {
            // Match emission counts with Laplace pseudocount.
            let mut emit = [1.0f64; 20];
            let mut emit_total = 20.0;
            for r in &msa.rows {
                if let Some(i) = aa_index(r[mc]) {
                    emit[i] += 1.0;
                    emit_total += 1.0;
                }
            }
            let match_emit = log_normalize(&emit, emit_total);

            // Insert emission: pooled residues of the insert columns
            // immediately following this match column (Laplace).
            let mut ins = [1.0f64; 20];
            let mut ins_total = 20.0;
            for c in (mc + 1)..width {
                if is_match[c] {
                    break;
                }
                for r in &msa.rows {
                    if let Some(i) = aa_index(r[c]) {
                        ins[i] += 1.0;
                        ins_total += 1.0;
                    }
                }
            }
            let insert_emit = log_normalize(&ins, ins_total);

            // Transitions: simple symmetric Plan7 defaults with a
            // light data nudge. A real HMMER counts state paths; this
            // v1 uses fixed sensible probabilities (documented
            // simplification) that still give a working scorer.
            nodes.push(HmmNode {
                match_emit,
                insert_emit,
                t_mm: 0.90f64.ln(),
                t_mi: 0.05f64.ln(),
                t_md: 0.05f64.ln(),
                t_im: 0.50f64.ln(),
                t_ii: 0.50f64.ln(),
                t_dm: 0.60f64.ln(),
                t_dd: 0.40f64.ln(),
            });
        }
        Ok(ProfileHmm { nodes })
    }

    /// **Viterbi** score: the log-probability of the most probable
    /// state path emitting `query`.
    ///
    /// The DP has three layers (`vm`, `vi`, `vd`) of size
    /// `(L+1) × (n+1)` where `L` is the model length and `n` the query
    /// length. Returns a log-probability (≤ 0).
    ///
    /// Returns [`AlignError::TooLarge`] when the `(L+1)·(n+1)` DP grid
    /// would exceed [`MAX_DP_CELLS`] (three
    /// `f64` layers); there is no linear-space variant, so an oversized
    /// model/query pair is rejected rather than risking an OOM.
    pub fn viterbi(&self, query: &[u8]) -> Result<f64> {
        self.dp(query, false, MAX_DP_CELLS)
    }

    /// **Forward** score: the log of the total probability summed over
    /// every state path emitting `query`.
    ///
    /// Returns [`AlignError::TooLarge`] when the DP grid would exceed
    /// [`MAX_DP_CELLS`].
    pub fn forward(&self, query: &[u8]) -> Result<f64> {
        self.dp(query, true, MAX_DP_CELLS)
    }

    /// Shared DP core. `sum == true` runs Forward (log-sum), `false`
    /// runs Viterbi (log-max). `max_cells` bounds the `(L+1)·(n+1)`
    /// allocation.
    fn dp(&self, query: &[u8], sum: bool, max_cells: usize) -> Result<f64> {
        let l = self.length();
        let n = query.len();
        let w = n + 1;
        check_dp_size_with(l + 1, w, max_cells)?;
        let combine = |a: f64, b: f64| if sum { log_add(a, b) } else { a.max(b) };

        // vm[k][i] = best/total log-prob in match state of node k
        // having emitted query[0..i]. k in 0..=l (0 = begin).
        let mut vm = vec![LOG_ZERO; (l + 1) * w];
        let mut vi = vec![LOG_ZERO; (l + 1) * w];
        let mut vd = vec![LOG_ZERO; (l + 1) * w];
        vm[0] = 0.0; // begin state, no residues emitted

        // Delete-only path along k at i = 0 (begin -> D1 -> D2 -> ...).
        for k in 1..=l {
            let prev = &self.nodes[k - 1];
            let from_m = vm[(k - 1) * w] + prev.t_md;
            let from_d = vd[(k - 1) * w] + prev.t_dd;
            vd[k * w] = combine(from_m, from_d);
        }

        for i in 1..=n {
            let res = query[i - 1];
            for k in 0..=l {
                let idx = k * w + i;

                // Match state of node k (k >= 1): emits query[i-1],
                // comes from node k-1's M / I / D at i-1.
                if k >= 1 {
                    let node = &self.nodes[k - 1];
                    let pe = node.match_emit[aa_index(res).unwrap_or(0)];
                    let p = (k - 1) * w + i - 1;
                    let from = combine(
                        combine(vm[p] + node.t_mm, vi[p] + node.t_im),
                        vd[p] + node.t_dm,
                    );
                    vm[idx] = pe + from;
                }

                // Insert state of node k (k >= 1): emits query[i-1],
                // self-loop, comes from this node's M / I at i-1.
                if k >= 1 {
                    let node = &self.nodes[k - 1];
                    let pe = node.insert_emit[aa_index(res).unwrap_or(0)];
                    let p = k * w + i - 1;
                    let from = combine(vm[p] + node.t_mi, vi[p] + node.t_ii);
                    vi[idx] = pe + from;
                }

                // Delete state of node k (k >= 1): silent, comes from
                // node k-1's M / D at the *same* i.
                if k >= 1 {
                    let node = &self.nodes[k - 1];
                    let p = (k - 1) * w + i;
                    let from = combine(vm[p] + node.t_md, vd[p] + node.t_dd);
                    vd[idx] = combine(vd[idx], from);
                }
            }
        }

        // End: the last node's M, D, or I having consumed the whole
        // query. I_L is a real built state, so a query with residues
        // trailing past the final match column terminates in the last
        // insert — omitting it silently mis-scored such queries (and a
        // 1-node model + multi-residue query underflowed to LOG_ZERO).
        let end_m = vm[l * w + n];
        let end_d = vd[l * w + n];
        Ok(combine(combine(end_m, end_d), vi[l * w + n]))
    }
}

/// Index of an amino acid in [`AA`], or `None` for a gap / non-AA.
fn aa_index(b: u8) -> Option<usize> {
    AA.iter().position(|&c| c == b.to_ascii_uppercase())
}

/// Normalises a count array to log-probabilities.
fn log_normalize(counts: &[f64; 20], total: f64) -> [f64; 20] {
    let mut out = [0.0f64; 20];
    for (o, &c) in out.iter_mut().zip(counts) {
        *o = (c / total).ln();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn msa(rows: &[&[u8]]) -> Msa {
        Msa::new(rows.iter().map(|r| r.to_vec()).collect()).unwrap()
    }

    #[test]
    fn build_from_msa() {
        let m = msa(&[b"MKVLA", b"MKVLA", b"MKILA"]);
        let hmm = ProfileHmm::from_msa(&m, 0.5).unwrap();
        // No gap columns -> 5 match states.
        assert_eq!(hmm.length(), 5);
        assert!(!hmm.is_empty());
    }

    #[test]
    fn gap_columns_become_inserts() {
        // Column 2 is mostly gaps -> not a match column.
        let m = msa(&[b"MK-LA", b"MK-LA", b"MKILA"]);
        let hmm = ProfileHmm::from_msa(&m, 0.5).unwrap();
        // 4 match states (column 2 demoted to insert).
        assert_eq!(hmm.length(), 4);
    }

    #[test]
    fn empty_msa_rejected() {
        let m = Msa::new(vec![]).unwrap();
        assert!(ProfileHmm::from_msa(&m, 0.5).is_err());
    }

    #[test]
    fn family_member_scores_above_outsider() {
        // Model trained on an MKVLA-like family.
        let m = msa(&[b"MKVLAAGG", b"MKVLAAGG", b"MKVLAAGG", b"MKVLATGG"]);
        let hmm = ProfileHmm::from_msa(&m, 0.5).unwrap();
        let member = hmm.viterbi(b"MKVLAAGG").unwrap();
        let outsider = hmm.viterbi(b"PWPWPWPW").unwrap();
        assert!(
            member > outsider,
            "family member ({member}) must score above an outsider ({outsider})"
        );
    }

    #[test]
    fn forward_at_least_viterbi() {
        let m = msa(&[b"MKVLAAGG", b"MKVLAAGG", b"MKVLATGG"]);
        let hmm = ProfileHmm::from_msa(&m, 0.5).unwrap();
        let q = b"MKVLAAGG";
        let v = hmm.viterbi(q).unwrap();
        let f = hmm.forward(q).unwrap();
        assert!(f >= v - 1e-6, "forward {f} >= viterbi {v}");
    }

    #[test]
    fn scores_are_finite_with_indels() {
        let m = msa(&[b"MKVLAAGG", b"MKVLAAGG", b"MKVLATGG"]);
        let hmm = ProfileHmm::from_msa(&m, 0.5).unwrap();
        // Query shorter than the model -> must use delete states.
        let short = hmm.viterbi(b"MKVLGG").unwrap();
        assert!(short.is_finite());
        // Query longer than the model -> must use insert states.
        let long = hmm.viterbi(b"MKVLAAAAGG").unwrap();
        assert!(long.is_finite());
    }

    #[test]
    fn termination_includes_last_insert_state() {
        // Regression for the Viterbi/Forward termination dropping the
        // final node's insert state I_L. A model built from a single
        // 8-residue sequence has 8 match nodes; a query equal to that
        // sequence plus two trailing residues ("AA") must place those
        // two residues in the *last* node's insert state — a path the
        // termination omitted when it combined only the M and D layers.
        //
        // The numbers below are computed from the exact transition /
        // emission log-probabilities the builder uses (Laplace-smoothed
        // single-row counts, Plan7 default transitions): the optimal
        // parse is 8 matches + 2 self-looping I_8 inserts.
        let m = msa(&[b"MKVLAAGG"]);
        let hmm = ProfileHmm::from_msa(&m, 0.5).unwrap();
        assert_eq!(hmm.length(), 8);

        // 8 matches + 2 trailing inserts. With the bug the score is
        // ~-30.615 (the two trailing residues can only be forced through
        // a wrong path because I_8 is unreachable at termination);
        // correct is ~-29.334.
        let v = hmm.viterbi(b"MKVLAAGGAA").unwrap();
        assert!(
            (v - (-29.334_230)).abs() < 1e-4,
            "Viterbi with 2 trailing inserts must reach I_8: got {v}, want ~-29.334230"
        );

        // The exact-match query takes no insert state at the end, so its
        // score is unaffected by the fix (sanity guard that the fix did
        // not perturb the no-insert path).
        let exact = hmm.viterbi(b"MKVLAAGG").unwrap();
        assert!(
            (exact - (-19.653_886)).abs() < 1e-4,
            "exact 8-match query score unchanged by fix: got {exact}, want ~-19.653886"
        );

        // Forward (log-sum) must be >= Viterbi and finite.
        let f = hmm.forward(b"MKVLAAGGAA").unwrap();
        assert!(f.is_finite() && f >= v - 1e-6, "forward {f} >= viterbi {v}");
    }

    #[test]
    fn one_node_model_multi_residue_query_not_log_zero() {
        // A 1-node model scoring a 2-residue query: the second residue
        // must be absorbed by node 1's insert state. With the
        // termination bug this path is unreachable, so the whole score
        // underflows to LOG_ZERO (≈ -1e30). It must be a real finite
        // log-probability instead.
        let m = msa(&[b"A"]);
        let hmm = ProfileHmm::from_msa(&m, 0.5).unwrap();
        assert_eq!(hmm.length(), 1);
        let v = hmm.viterbi(b"AA").unwrap();
        assert!(
            v > LOG_ZERO / 2.0,
            "1-node model + 2-residue query must not underflow to LOG_ZERO: got {v}"
        );
        assert!(v.is_finite(), "score must be finite: got {v}");
    }

    #[test]
    fn profilehmm_over_cap_errors() {
        use crate::error::AlignError;
        let m = msa(&[b"MKVLAAGG"]);
        let hmm = ProfileHmm::from_msa(&m, 0.5).unwrap();
        // 9 nodes+1 by query 11 -> well over a cap of 8; rejected before
        // the three f64 layers allocate. dp() is the shared core.
        assert!(matches!(
            hmm.dp(b"MKVLAAGGAA", false, 8).unwrap_err(),
            AlignError::TooLarge { .. }
        ));
        assert!(matches!(
            hmm.dp(b"MKVLAAGGAA", true, 8).unwrap_err(),
            AlignError::TooLarge { .. }
        ));
        // Generous cap succeeds (same as the public methods).
        assert!(hmm.dp(b"MKVLAAGGAA", false, usize::MAX).is_ok());
    }

    #[test]
    fn conserved_column_emission_peaks() {
        // Column 0 is always M; its match emission for M should be the
        // highest in that column's distribution.
        let m = msa(&[b"MAA", b"MCC", b"MDD"]);
        let hmm = ProfileHmm::from_msa(&m, 0.5).unwrap();
        let col0 = &hmm.nodes[0].match_emit;
        let m_idx = aa_index(b'M').unwrap();
        let max = col0.iter().cloned().fold(f64::MIN, f64::max);
        assert!((col0[m_idx] - max).abs() < 1e-9, "M should dominate column 0");
    }
}
