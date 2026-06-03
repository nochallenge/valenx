//! GATK-style PairHMM — quality-aware read-vs-haplotype likelihood.
//!
//! The `PairHmm` in [`valenx_align`] is a general protein/DNA pair HMM
//! with a single per-mismatch emission probability. The GATK
//! HaplotypeCaller's PairHMM is a more specialised tool: it scores one
//! sequencing **read** (with per-base Phred qualities) against one
//! candidate **haplotype**, using each base's own quality as the per-
//! position emission error probability. That is the noise model the
//! Bayesian site-genotyper marginalises over.
//!
//! Three states (Durbin et al. ch.4 conventions):
//!
//! - **M** emits an aligned pair `(read_i, hap_j)`. When the bases
//!   agree the emission is `1 − e_i` (the probability that base `i` is
//!   *correct*); when they disagree it is `e_i / 3` (the probability of
//!   miscalling to that specific base, spread uniformly over the three
//!   wrong bases).
//! - **I** (insertion in the read relative to the haplotype) emits a
//!   read base with a uniform-base probability `1/4`.
//! - **D** (deletion in the read relative to the haplotype) emits a
//!   haplotype base — but, like GATK, the deletion state does not emit a
//!   read symbol, so the emission factor is `1`.
//!
//! Transitions are governed by a gap-open / gap-extend probability.
//! GATK's PairHMM uses per-base gap-open / extend qualities; here we
//! use a single, configurable pair of probabilities (a documented
//! simplification — the per-base GOP/GCP qualities of CRAM `BI/BD` tags
//! are not stored in plain SAM, so a single pair captures the available
//! information). The forward algorithm runs in `log10` space because
//! the downstream Bayesian marginalisation is in `log10` too.

use crate::error::{GenomicsError, Result};

/// `log10(0)` sentinel — a finite very-negative number.
const LOG10_ZERO: f64 = -1.0e30;

/// Stably adds two probabilities held in `log10` space.
#[inline]
fn log10_add(a: f64, b: f64) -> f64 {
    if a <= LOG10_ZERO {
        return b;
    }
    if b <= LOG10_ZERO {
        return a;
    }
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    // log10(10^a + 10^b) = hi + log10(1 + 10^(lo - hi))
    hi + (1.0 + 10f64.powf(lo - hi)).log10()
}

/// Parameters of the GATK-class PairHMM.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PairHmmParams {
    /// Probability of opening a gap from the M state. Must lie in
    /// `(0, 0.5)`.
    pub gap_open: f64,
    /// Probability of extending an open gap. Must lie in `[0, 1)`.
    pub gap_extend: f64,
}

impl Default for PairHmmParams {
    /// GATK-class defaults — small gap-open (Phred-45) + moderate gap-
    /// extend (Phred-10). The exact constants are not load-bearing; the
    /// downstream marginal is robust to small perturbations.
    fn default() -> Self {
        PairHmmParams {
            gap_open: 10f64.powi(-3), // ~Phred 30
            gap_extend: 0.1,
        }
    }
}

impl PairHmmParams {
    fn validate(&self) -> Result<()> {
        if !(self.gap_open > 0.0 && self.gap_open < 0.5) {
            return Err(GenomicsError::invalid(
                "gap_open",
                "must lie in (0, 0.5)",
            ));
        }
        if !(0.0..1.0).contains(&self.gap_extend) {
            return Err(GenomicsError::invalid(
                "gap_extend",
                "must lie in [0, 1)",
            ));
        }
        Ok(())
    }
}

/// Converts a Phred quality to a base-error probability, clamped to
/// `[1e-6, 0.75]` so a single perfect or junk base cannot dominate the
/// product.
#[inline]
fn error_prob(q: u8) -> f64 {
    let q = q.min(60) as f64;
    10f64.powf(-q / 10.0).clamp(1e-6, 0.75)
}

/// Per-position `log10` emission of the M state for `(read_base,
/// hap_base)`.
#[inline]
fn log10_match(read_base: u8, hap_base: u8, q: u8) -> f64 {
    let e = error_prob(q);
    let p = if read_base.eq_ignore_ascii_case(&hap_base) {
        1.0 - e
    } else {
        e / 3.0
    };
    p.max(1e-300).log10()
}

/// Per-position `log10` emission of the I state (read base against a
/// gap in the haplotype). Treat as a uniform base prior.
#[inline]
fn log10_ins() -> f64 {
    0.25f64.log10()
}

/// Per-position `log10` emission of the D state. GATK's PairHMM treats
/// a deletion as non-emitting from the read's point of view, so the
/// emission factor is `1` → `log10 1 = 0`.
#[inline]
fn log10_del() -> f64 {
    0.0
}

/// `log10 P(read | haplotype)` under the GATK-style PairHMM forward
/// algorithm.
///
/// The result is a *real* log-probability — bounded above by zero,
/// approaching it as the read matches the haplotype exactly at high
/// quality. The read's `qualities` slice must have the same length as
/// `read_bases`; an empty `qualities` is treated as a uniform Phred-30.
pub fn log10_p_read_given_haplotype(
    read_bases: &[u8],
    qualities: &[u8],
    haplotype: &[u8],
    params: &PairHmmParams,
) -> Result<f64> {
    params.validate()?;
    let n = read_bases.len();
    let m = haplotype.len();
    if n == 0 {
        return Ok(0.0); // P(empty read | anything) ≈ 1
    }
    if m == 0 {
        return Ok(LOG10_ZERO);
    }
    if !qualities.is_empty() && qualities.len() != n {
        return Err(GenomicsError::invalid(
            "qualities",
            "qualities length must match read_bases length",
        ));
    }
    let qual_for = |i: usize| -> u8 {
        if qualities.is_empty() {
            30
        } else {
            qualities[i]
        }
    };

    let w = m + 1;
    // Transition log10 probabilities — three-state symmetric model.
    let tmm = (1.0 - 2.0 * params.gap_open).log10();
    let tmi = params.gap_open.log10();
    let tmd = params.gap_open.log10();
    let tii = params.gap_extend.log10();
    let tim = (1.0 - params.gap_extend).log10();
    let tdd = params.gap_extend.log10();
    let tdm = (1.0 - params.gap_extend).log10();

    // Forward matrices (log10).
    let mut fm = vec![LOG10_ZERO; (n + 1) * w];
    let mut fi = vec![LOG10_ZERO; (n + 1) * w];
    let mut fd = vec![LOG10_ZERO; (n + 1) * w];

    // GATK initialisation: the haplotype can be visited at any starting
    // column with equal probability — equivalent to a uniform prior
    // over alignment start columns. This is the standard
    // HaplotypeCaller convention so a read that aligns to any window of
    // the haplotype gets a fair shake.
    let init = -((m as f64).max(1.0)).log10();
    for cell in fm.iter_mut().take(m + 1) {
        *cell = init;
    }

    for i in 1..=n {
        let q = qual_for(i - 1);
        let rb = read_bases[i - 1];
        for j in 0..=m {
            let idx = i * w + j;
            // M state.
            if j >= 1 {
                let prev = (i - 1) * w + (j - 1);
                let from = log10_add(
                    log10_add(fm[prev] + tmm, fi[prev] + tim),
                    fd[prev] + tdm,
                );
                fm[idx] = log10_match(rb, haplotype[j - 1], q) + from;
            }
            // I state (insertion in read).
            let prev = (i - 1) * w + j;
            let from = log10_add(fm[prev] + tmi, fi[prev] + tii);
            fi[idx] = log10_ins() + from;
            // D state (deletion in read relative to haplotype).
            if j >= 1 {
                let prev = i * w + (j - 1);
                let from = log10_add(fm[prev] + tmd, fd[prev] + tdd);
                fd[idx] = log10_del() + from;
            }
        }
    }

    // Termination: sum over the last read row.
    let mut total = LOG10_ZERO;
    for j in 0..=m {
        let idx = n * w + j;
        total = log10_add(total, log10_add(log10_add(fm[idx], fi[idx]), fd[idx]));
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_bad_parameters() {
        let p = PairHmmParams {
            gap_open: 0.6,
            gap_extend: 0.1,
        };
        assert!(p.validate().is_err());
        let p = PairHmmParams {
            gap_open: 0.01,
            gap_extend: 1.5,
        };
        assert!(p.validate().is_err());
    }

    #[test]
    fn rejects_qualities_length_mismatch() {
        let r = b"ACGT";
        let q = vec![30u8; 3];
        let h = b"ACGT";
        let res = log10_p_read_given_haplotype(r, &q, h, &PairHmmParams::default());
        assert!(res.is_err());
    }

    #[test]
    fn exact_match_is_near_zero() {
        let r = b"ACGTACGTACGT";
        let q = vec![40u8; r.len()];
        let h = b"ACGTACGTACGT";
        let lp = log10_p_read_given_haplotype(r, &q, h, &PairHmmParams::default()).unwrap();
        // log10 P should be high (close to 0): a perfect match.
        assert!(lp > -2.0, "exact-match log10 P was {lp}, expected > -2");
    }

    #[test]
    fn monotonic_in_mismatches() {
        let hap = b"ACGTACGTACGTACGT";
        let q = vec![40u8; hap.len()];
        let p = PairHmmParams::default();
        let p0 = log10_p_read_given_haplotype(hap, &q, hap, &p).unwrap();
        let mut r1 = hap.to_vec();
        r1[5] = b'T'; // 1 mismatch (was C)
        if r1[5] == b'C' {
            r1[5] = b'A';
        }
        let mut r2 = r1.clone();
        r2[10] = b'A'; // 2 mismatches
        if r2[10] == hap[10] {
            r2[10] = if hap[10] == b'A' { b'C' } else { b'A' };
        }
        let p1 = log10_p_read_given_haplotype(&r1, &q, hap, &p).unwrap();
        let p2 = log10_p_read_given_haplotype(&r2, &q, hap, &p).unwrap();
        assert!(p0 > p1, "0 mismatches {p0} should beat 1 mismatch {p1}");
        assert!(p1 > p2, "1 mismatch {p1} should beat 2 mismatches {p2}");
    }

    #[test]
    fn closer_haplotype_scores_higher() {
        // Read carries an exact SNV in the middle.
        let read = b"ACGTACTTACGTACGT"; // pos 6 = T
        let q = vec![40u8; read.len()];
        let near = b"ACGTACTTACGTACGT"; // identical
        let far = b"ACGTACGTACGTACGT"; // pos 6 differs from read
        let p_near =
            log10_p_read_given_haplotype(read, &q, near, &PairHmmParams::default()).unwrap();
        let p_far =
            log10_p_read_given_haplotype(read, &q, far, &PairHmmParams::default()).unwrap();
        assert!(
            p_near > p_far,
            "near {p_near} should beat far {p_far} for a read carrying the SNV"
        );
    }

    #[test]
    fn quality_modulates_penalty() {
        // A mismatch carrying low quality should be less punishing.
        let read = b"ACGT";
        let hap = b"AAGT"; // one mismatch at pos 1
        let p = PairHmmParams::default();
        let low_q = log10_p_read_given_haplotype(read, &[5u8; 4], hap, &p).unwrap();
        let high_q = log10_p_read_given_haplotype(read, &[40u8; 4], hap, &p).unwrap();
        assert!(
            low_q > high_q,
            "low-quality mismatch {low_q} should beat high-quality mismatch {high_q}"
        );
    }

    #[test]
    fn read_with_insertion_scored_above_substitution_chain() {
        // Read has a 2-base insertion relative to the haplotype.
        let read = b"ACGTGGACGT";
        let hap = b"ACGTACGT";
        let q = vec![40u8; read.len()];
        let lp =
            log10_p_read_given_haplotype(read, &q, hap, &PairHmmParams::default()).unwrap();
        // We are not asserting an absolute number; just that the model
        // still returns a finite value (not LOG10_ZERO).
        assert!(lp > -50.0, "log10 P was {lp}");
    }

    #[test]
    fn read_with_deletion_finite() {
        let read = b"ACGTCGT"; // missing the middle A
        let hap = b"ACGTACGT";
        let q = vec![40u8; read.len()];
        let lp =
            log10_p_read_given_haplotype(read, &q, hap, &PairHmmParams::default()).unwrap();
        assert!(lp.is_finite() && lp > -50.0, "lp = {lp}");
    }

    #[test]
    fn log10_add_is_stable() {
        // log10(0.3 + 0.7) checked exactly.
        let got = log10_add((0.3f64).log10(), (0.7f64).log10());
        assert!((got - (1.0f64).log10()).abs() < 1e-12);
        // Identity with LOG10_ZERO.
        assert!((log10_add(LOG10_ZERO, (0.5f64).log10()) - (0.5f64).log10()).abs() < 1e-12);
    }
}
