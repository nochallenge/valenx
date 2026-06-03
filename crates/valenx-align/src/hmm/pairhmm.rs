//! Pair hidden Markov model — forward algorithm and Viterbi.
//!
//! A *pair HMM* generates two sequences jointly. It has three emitting
//! states (Durbin et al., *Biological Sequence Analysis*, ch. 4):
//!
//! - **M** — emits an aligned pair `(x_i, y_j)`;
//! - **X** — emits `x_i` against a gap in `y`;
//! - **Y** — emits `y_j` against a gap in `x`.
//!
//! Transitions are governed by `δ` (probability of opening a gap from
//! M) and `ε` (probability of extending a gap). Two algorithms:
//!
//! - [`PairHmm::viterbi`] — the single most probable state path,
//!   i.e. the most probable *alignment* (the probabilistic analogue of
//!   Gotoh).
//! - [`PairHmm::forward`] — the total probability summed over *all*
//!   alignments, the proper likelihood that the two sequences are
//!   related.
//!
//! All DP is done in **log space** (natural log) so long sequences do
//! not underflow.

use crate::error::{AlignError, Result};
use crate::limits::{check_dp_size_with, MAX_DP_CELLS};

/// `ln(0)` sentinel — a very negative number that stays finite under
/// addition of normal log-probabilities.
const LOG_ZERO: f64 = -1.0e30;

/// Adds two probabilities held in log space: returns
/// `ln(e^a + e^b)`, computed stably.
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

/// A pair HMM with a symmetric three-state (M/X/Y) topology.
///
/// Emission probabilities are supplied as a *match* matrix `p(x, y)`
/// and a background `q(x)`; transition probabilities by the gap-open
/// `delta` and gap-extend `epsilon` parameters. All inputs are
/// probabilities in `[0, 1]`.
#[derive(Clone, Debug)]
pub struct PairHmm {
    /// `match_emit[a][b]` = probability state M emits residues `(a,b)`.
    /// Indexed by [`residue_index`].
    match_emit: [[f64; ALPHA]; ALPHA],
    /// `background[a]` = probability an X/Y state emits residue `a`.
    background: [f64; ALPHA],
    /// Gap-open probability (M → X or M → Y).
    delta: f64,
    /// Gap-extend probability (X → X or Y → Y).
    epsilon: f64,
}

/// Alphabet size: the 20 amino acids + a wildcard slot. DNA uses the
/// `ACGT` subset of the same indices.
const ALPHA: usize = 21;

/// Maps a residue byte to an index in `[0, ALPHA)`. Unknown residues
/// fold to the last (wildcard) slot.
pub fn residue_index(b: u8) -> usize {
    const ORDER: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
    ORDER
        .iter()
        .position(|&c| c == b.to_ascii_uppercase())
        .unwrap_or(ALPHA - 1)
}

impl PairHmm {
    /// Builds a pair HMM.
    ///
    /// `match_prob(a, b)` returns the M-state emission probability for
    /// residues `a, b`; `background(a)` the X/Y emission probability.
    /// `delta` and `epsilon` are the gap-open / gap-extend
    /// probabilities and must satisfy `0 < delta`, `0 <= epsilon < 1`,
    /// `2·delta < 1`.
    pub fn new(
        match_prob: impl Fn(u8, u8) -> f64,
        background: impl Fn(u8) -> f64,
        delta: f64,
        epsilon: f64,
    ) -> Result<Self> {
        if !(0.0..0.5).contains(&delta) {
            return Err(AlignError::invalid("delta", "gap-open must be in [0, 0.5)"));
        }
        if !(0.0..1.0).contains(&epsilon) {
            return Err(AlignError::invalid(
                "epsilon",
                "gap-extend must be in [0, 1)",
            ));
        }
        const ORDER: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
        let mut match_emit = [[0.0f64; ALPHA]; ALPHA];
        let mut bg = [0.0f64; ALPHA];
        for (i, &a) in ORDER.iter().enumerate() {
            bg[i] = background(a).max(1e-12);
            for (j, &b) in ORDER.iter().enumerate() {
                match_emit[i][j] = match_prob(a, b).max(1e-12);
            }
        }
        // Wildcard slot: small uniform fallback.
        bg[ALPHA - 1] = 1.0 / ALPHA as f64;
        for row in match_emit.iter_mut() {
            row[ALPHA - 1] = 1e-6;
        }
        for j in 0..ALPHA {
            match_emit[ALPHA - 1][j] = 1e-6;
        }
        Ok(PairHmm {
            match_emit,
            background: bg,
            delta,
            epsilon,
        })
    }

    /// A simple symmetric DNA pair HMM: `p_match` for identical bases,
    /// `p_mismatch` for differing bases, uniform `0.25` background.
    pub fn dna(p_match: f64, p_mismatch: f64, delta: f64, epsilon: f64) -> Result<Self> {
        Self::new(
            |a, b| {
                if a.eq_ignore_ascii_case(&b) {
                    p_match
                } else {
                    p_mismatch
                }
            },
            |_| 0.25,
            delta,
            epsilon,
        )
    }

    /// Log emission probability of the M state for residues `(a, b)`.
    fn log_match(&self, a: u8, b: u8) -> f64 {
        self.match_emit[residue_index(a)][residue_index(b)].ln()
    }

    /// Log emission probability of an X/Y state for residue `a`.
    fn log_bg(&self, a: u8) -> f64 {
        self.background[residue_index(a)].ln()
    }

    /// **Viterbi**: the log-probability of the single most probable
    /// state path generating `x` and `y` — the most probable
    /// alignment's score.
    ///
    /// Returns [`AlignError::TooLarge`] when the `(n+1)·(m+1)` DP grid
    /// would exceed [`MAX_DP_CELLS`] (three
    /// `f64` layers); a pair HMM has no linear-space variant, so an
    /// oversized input is rejected rather than risking an OOM.
    pub fn viterbi(&self, x: &[u8], y: &[u8]) -> Result<f64> {
        self.viterbi_capped(x, y, MAX_DP_CELLS)
    }

    fn viterbi_capped(&self, x: &[u8], y: &[u8], max_cells: usize) -> Result<f64> {
        let n = x.len();
        let m = y.len();
        let w = m + 1;
        check_dp_size_with(n + 1, m + 1, max_cells)?;
        let trans_mm = (1.0 - 2.0 * self.delta).ln();
        let trans_mg = self.delta.ln();
        let trans_gg = self.epsilon.ln();
        let trans_gm = (1.0 - self.epsilon).ln();

        let mut vm = vec![LOG_ZERO; (n + 1) * w];
        let mut vx = vec![LOG_ZERO; (n + 1) * w];
        let mut vy = vec![LOG_ZERO; (n + 1) * w];
        vm[0] = 0.0; // begin in M with log-prob 0

        for i in 0..=n {
            for j in 0..=m {
                if i == 0 && j == 0 {
                    continue;
                }
                let idx = i * w + j;
                if i > 0 && j > 0 {
                    let prev = (i - 1) * w + j - 1;
                    let best = vm[prev]
                        .max(vx[prev] - 0.0)
                        .max(vy[prev]);
                    let from = (vm[prev] + trans_mm)
                        .max(vx[prev] + trans_gm)
                        .max(vy[prev] + trans_gm);
                    let _ = best;
                    vm[idx] = self.log_match(x[i - 1], y[j - 1]) + from;
                }
                if i > 0 {
                    let prev = (i - 1) * w + j;
                    let from = (vm[prev] + trans_mg).max(vx[prev] + trans_gg);
                    vx[idx] = self.log_bg(x[i - 1]) + from;
                }
                if j > 0 {
                    let prev = i * w + j - 1;
                    let from = (vm[prev] + trans_mg).max(vy[prev] + trans_gg);
                    vy[idx] = self.log_bg(y[j - 1]) + from;
                }
            }
        }
        let last = n * w + m;
        Ok(vm[last].max(vx[last]).max(vy[last]))
    }

    /// **Forward**: the log of the *total* probability summed over
    /// every alignment of `x` and `y` — the likelihood the two
    /// sequences were generated by the model.
    ///
    /// Returns [`AlignError::TooLarge`] when the `(n+1)·(m+1)` DP grid
    /// would exceed [`MAX_DP_CELLS`].
    pub fn forward(&self, x: &[u8], y: &[u8]) -> Result<f64> {
        self.forward_capped(x, y, MAX_DP_CELLS)
    }

    fn forward_capped(&self, x: &[u8], y: &[u8], max_cells: usize) -> Result<f64> {
        let n = x.len();
        let m = y.len();
        let w = m + 1;
        check_dp_size_with(n + 1, m + 1, max_cells)?;
        let trans_mm = (1.0 - 2.0 * self.delta).ln();
        let trans_mg = self.delta.ln();
        let trans_gg = self.epsilon.ln();
        let trans_gm = (1.0 - self.epsilon).ln();

        let mut fm = vec![LOG_ZERO; (n + 1) * w];
        let mut fx = vec![LOG_ZERO; (n + 1) * w];
        let mut fy = vec![LOG_ZERO; (n + 1) * w];
        fm[0] = 0.0;

        for i in 0..=n {
            for j in 0..=m {
                if i == 0 && j == 0 {
                    continue;
                }
                let idx = i * w + j;
                if i > 0 && j > 0 {
                    let prev = (i - 1) * w + j - 1;
                    let from = log_add(
                        log_add(fm[prev] + trans_mm, fx[prev] + trans_gm),
                        fy[prev] + trans_gm,
                    );
                    fm[idx] = self.log_match(x[i - 1], y[j - 1]) + from;
                }
                if i > 0 {
                    let prev = (i - 1) * w + j;
                    let from = log_add(fm[prev] + trans_mg, fx[prev] + trans_gg);
                    fx[idx] = self.log_bg(x[i - 1]) + from;
                }
                if j > 0 {
                    let prev = i * w + j - 1;
                    let from = log_add(fm[prev] + trans_mg, fy[prev] + trans_gg);
                    fy[idx] = self.log_bg(y[j - 1]) + from;
                }
            }
        }
        let last = n * w + m;
        Ok(log_add(log_add(fm[last], fx[last]), fy[last]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hmm() -> PairHmm {
        // Strong match signal, modest gap penalties.
        PairHmm::dna(0.9, 0.033, 0.1, 0.2).unwrap()
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(PairHmm::dna(0.9, 0.03, 0.6, 0.2).is_err()); // delta >= 0.5
        assert!(PairHmm::dna(0.9, 0.03, 0.1, 1.5).is_err()); // epsilon >= 1
    }

    #[test]
    fn forward_at_least_viterbi() {
        // Forward sums over all paths => >= the single best path.
        let h = hmm();
        let x = b"ACGTACGT";
        let y = b"ACGTACGT";
        let v = h.viterbi(x, y).unwrap();
        let f = h.forward(x, y).unwrap();
        assert!(f >= v - 1e-6, "forward {f} should be >= viterbi {v}");
    }

    #[test]
    fn identical_scores_higher_than_divergent() {
        let h = hmm();
        let q = b"ACGTACGT";
        let same = h.viterbi(q, b"ACGTACGT").unwrap();
        let diff = h.viterbi(q, b"TGCATGCA").unwrap();
        assert!(same > diff, "identical pair must score higher");
    }

    #[test]
    fn longer_identical_still_consistent() {
        let h = hmm();
        let x = b"ACGTACGTACGTACGT";
        let f = h.forward(x, x).unwrap();
        let v = h.viterbi(x, x).unwrap();
        assert!(f.is_finite() && v.is_finite());
        assert!(f >= v - 1e-6);
    }

    #[test]
    fn gap_tolerated_by_forward() {
        // y is x with one residue deleted; forward should stay finite
        // and not catastrophically low.
        let h = hmm();
        let f = h.forward(b"ACGTACGT", b"ACGTCGT").unwrap();
        assert!(f.is_finite());
        assert!(f > -200.0);
    }

    #[test]
    fn pairhmm_over_cap_errors() {
        use crate::error::AlignError;
        let h = hmm();
        let x = b"ACGTACGT";
        let y = b"ACGTACGT";
        // 9*9 = 81 cells; a cap of 8 rejects without the three f64
        // matrices being allocated.
        assert!(matches!(
            h.viterbi_capped(x, y, 8).unwrap_err(),
            AlignError::TooLarge { .. }
        ));
        assert!(matches!(
            h.forward_capped(x, y, 8).unwrap_err(),
            AlignError::TooLarge { .. }
        ));
        // Generous cap matches the public methods.
        assert!(h.viterbi_capped(x, y, usize::MAX).is_ok());
        assert!(h.forward_capped(x, y, usize::MAX).is_ok());
    }

    #[test]
    fn log_add_is_stable() {
        // log_add(ln a, ln b) == ln(a+b).
        let a = 0.3f64;
        let b = 0.7f64;
        let got = log_add(a.ln(), b.ln());
        assert!((got - (a + b).ln()).abs() < 1e-9);
        // Identity with LOG_ZERO.
        assert!((log_add(LOG_ZERO, a.ln()) - a.ln()).abs() < 1e-9);
    }
}
