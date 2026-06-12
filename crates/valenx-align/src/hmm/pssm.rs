//! Position-specific scoring matrices, position-weight matrices, and
//! motif search.
//!
//! A **PSSM** (a.k.a. PWM) is a per-position log-odds score table: for
//! each position of a fixed-length motif and each residue it stores
//! `log2( observed-frequency / background-frequency )`. Sliding it
//! across a sequence and summing the per-position scores is the
//! classic transcription-factor-binding-site / sequence-motif scanner.
//!
//! [`Pssm::from_sites`] builds a matrix from a set of equal-length
//! aligned sites (a Laplace pseudocount avoids `log 0`);
//! [`Pssm::scan`] reports every window whose score clears a threshold.

use crate::error::{AlignError, Result};

/// A position-specific scoring matrix over a fixed residue alphabet.
#[derive(Clone, Debug, PartialEq)]
pub struct Pssm {
    /// The residue alphabet, uppercased — e.g. `b"ACGT"`.
    alphabet: Vec<u8>,
    /// Dense byte → alphabet-index map (`255` = absent).
    idx: [u8; 256],
    /// `scores[pos][a]` = log2-odds score of residue index `a` at
    /// position `pos`. Outer length = motif width.
    scores: Vec<Vec<f64>>,
}

impl Pssm {
    /// Builds a PSSM from a set of *aligned, equal-length* sites.
    ///
    /// `background[i]` is the prior probability of `alphabet[i]`
    /// (uniform if you have none). A Laplace `+1` pseudocount is added
    /// to every observed count. Returns [`AlignError::Invalid`] for an
    /// empty site set / alphabet and [`AlignError::Dimension`] if the
    /// sites are not all the same length.
    pub fn from_sites(sites: &[&[u8]], alphabet: &[u8], background: &[f64]) -> Result<Self> {
        if sites.is_empty() {
            return Err(AlignError::invalid("sites", "need >= 1 site"));
        }
        if alphabet.is_empty() {
            return Err(AlignError::invalid("alphabet", "alphabet is empty"));
        }
        if background.len() != alphabet.len() {
            return Err(AlignError::dimension(format!(
                "background has {} entries, alphabet {}",
                background.len(),
                alphabet.len()
            )));
        }
        let width = sites[0].len();
        if width == 0 {
            return Err(AlignError::invalid("sites", "sites are zero-length"));
        }
        for s in sites {
            if s.len() != width {
                return Err(AlignError::dimension(format!(
                    "site length {} != motif width {width}",
                    s.len()
                )));
            }
        }

        let upper: Vec<u8> = alphabet.iter().map(|b| b.to_ascii_uppercase()).collect();
        let mut idx = [255u8; 256];
        for (i, &b) in upper.iter().enumerate() {
            idx[b as usize] = i as u8;
        }
        let a = upper.len();

        let mut scores = vec![vec![0.0f64; a]; width];
        for (pos, score_row) in scores.iter_mut().enumerate() {
            // Counts with Laplace pseudocount.
            let mut counts = vec![1.0f64; a];
            let mut total = a as f64;
            for s in sites {
                let c = idx[s[pos].to_ascii_uppercase() as usize];
                if c != 255 {
                    counts[c as usize] += 1.0;
                    total += 1.0;
                }
            }
            for (ai, slot) in score_row.iter_mut().enumerate() {
                let freq = counts[ai] / total;
                let bg = background[ai].max(1e-9);
                *slot = (freq / bg).log2();
            }
        }

        Ok(Pssm {
            alphabet: upper,
            idx,
            scores,
        })
    }

    /// Convenience: build a DNA PSSM with a uniform `0.25` background.
    pub fn dna_from_sites(sites: &[&[u8]]) -> Result<Self> {
        Self::from_sites(sites, b"ACGT", &[0.25, 0.25, 0.25, 0.25])
    }

    /// The motif width (number of positions).
    pub fn width(&self) -> usize {
        self.scores.len()
    }

    /// The residue alphabet.
    pub fn alphabet(&self) -> &[u8] {
        &self.alphabet
    }

    /// Log-odds score of residue `residue` at position `pos`. Residues
    /// outside the alphabet score `0` (neutral).
    pub fn score_at(&self, pos: usize, residue: u8) -> f64 {
        let c = self.idx[residue.to_ascii_uppercase() as usize];
        if pos >= self.width() || c == 255 {
            0.0
        } else {
            self.scores[pos][c as usize]
        }
    }

    /// Total PSSM score of one window — `window.len()` must equal the
    /// motif width. Returns [`AlignError::Dimension`] otherwise.
    pub fn score_window(&self, window: &[u8]) -> Result<f64> {
        if window.len() != self.width() {
            return Err(AlignError::dimension(format!(
                "window length {} != motif width {}",
                window.len(),
                self.width()
            )));
        }
        Ok(window
            .iter()
            .enumerate()
            .map(|(pos, &r)| self.score_at(pos, r))
            .sum())
    }

    /// The maximum score the PSSM can award — the sum of each
    /// position's best residue. Useful for normalising hit scores.
    pub fn max_score(&self) -> f64 {
        self.scores
            .iter()
            .map(|row| row.iter().cloned().fold(f64::MIN, f64::max))
            .sum()
    }

    /// The minimum possible score — the sum of each position's worst
    /// residue.
    pub fn min_score(&self) -> f64 {
        self.scores
            .iter()
            .map(|row| row.iter().cloned().fold(f64::MAX, f64::min))
            .sum()
    }

    /// Scans `seq`, returning every `(offset, score)` window whose
    /// score is `>= threshold`, sorted by descending score.
    pub fn scan(&self, seq: &[u8], threshold: f64) -> Vec<MotifHit> {
        let w = self.width();
        let mut hits = Vec::new();
        if seq.len() >= w {
            for offset in 0..=seq.len() - w {
                let score = self
                    .score_window(&seq[offset..offset + w])
                    .expect("window width matches motif");
                if score >= threshold {
                    hits.push(MotifHit { offset, score });
                }
            }
        }
        hits.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        hits
    }

    /// The single best-scoring window of `seq`, or `None` if `seq` is
    /// shorter than the motif.
    pub fn best_match(&self, seq: &[u8]) -> Option<MotifHit> {
        self.scan(seq, f64::NEG_INFINITY).into_iter().next()
    }
}

/// A motif occurrence found by [`Pssm::scan`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MotifHit {
    /// 0-based start offset of the window in the scanned sequence.
    pub offset: usize,
    /// The PSSM log-odds score of that window.
    pub score: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_bad_input() {
        assert!(Pssm::dna_from_sites(&[]).is_err());
        // Ragged sites.
        let sites: &[&[u8]] = &[b"ACGT", b"ACG"];
        assert!(Pssm::dna_from_sites(sites).is_err());
    }

    #[test]
    fn perfect_motif_scores_high() {
        // Every site is "ACGT" -> consensus ACGT scores near-maximal.
        let sites: &[&[u8]] = &[b"ACGT", b"ACGT", b"ACGT", b"ACGT"];
        let pssm = Pssm::dna_from_sites(sites).unwrap();
        let consensus = pssm.score_window(b"ACGT").unwrap();
        let anti = pssm.score_window(b"TGCA").unwrap();
        assert!(consensus > anti);
        assert!(consensus > 0.0, "consensus should be a positive log-odds");
    }

    #[test]
    fn score_window_dimension_check() {
        let sites: &[&[u8]] = &[b"ACGT", b"ACGT"];
        let pssm = Pssm::dna_from_sites(sites).unwrap();
        assert!(pssm.score_window(b"ACG").is_err());
        assert!(pssm.score_window(b"ACGT").is_ok());
    }

    #[test]
    fn scan_finds_motif_occurrence() {
        // Build a motif from "GATC" sites, then scan a sequence that
        // contains "GATC" embedded in noise.
        let sites: &[&[u8]] = &[b"GATC", b"GATC", b"GATC", b"GATC"];
        let pssm = Pssm::dna_from_sites(sites).unwrap();
        let seq = b"TTTTGATCTTTT";
        let hits = pssm.scan(seq, 0.0);
        assert!(!hits.is_empty());
        // The top hit should be at offset 4 (where GATC sits).
        assert_eq!(hits[0].offset, 4);
    }

    #[test]
    fn best_match_picks_strongest() {
        let sites: &[&[u8]] = &[b"AAAA", b"AAAA", b"AAAA"];
        let pssm = Pssm::dna_from_sites(sites).unwrap();
        let seq = b"CCAAAACC";
        let best = pssm.best_match(seq).unwrap();
        assert_eq!(best.offset, 2);
    }

    #[test]
    fn max_and_min_bounds() {
        let sites: &[&[u8]] = &[b"ACGT", b"ACGT"];
        let pssm = Pssm::dna_from_sites(sites).unwrap();
        let max = pssm.max_score();
        let min = pssm.min_score();
        assert!(max > min);
        // Any real window scores within [min, max].
        let s = pssm.score_window(b"AAAA").unwrap();
        assert!(s <= max + 1e-9 && s >= min - 1e-9);
    }

    #[test]
    fn short_sequence_no_hits() {
        let sites: &[&[u8]] = &[b"ACGTACGT", b"ACGTACGT"];
        let pssm = Pssm::dna_from_sites(sites).unwrap();
        assert!(pssm.scan(b"ACG", 0.0).is_empty());
        assert!(pssm.best_match(b"ACG").is_none());
    }

    #[test]
    fn protein_pssm() {
        let sites: &[&[u8]] = &[b"MKVL", b"MKVL", b"MRVL"];
        let bg = [0.05f64; 20];
        let pssm = Pssm::from_sites(sites, b"ACDEFGHIKLMNPQRSTVWY", &bg).unwrap();
        assert_eq!(pssm.width(), 4);
        // Position 0 is always M -> M scores positively there.
        assert!(pssm.score_at(0, b'M') > 0.0);
    }
}
