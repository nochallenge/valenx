//! Guide-RNA design — PAM scanning and on-target efficiency scoring.
//!
//! CRISPR guide design (CHOPCHOP, Benchling, GuideScan) starts by
//! finding every candidate protospacer in a target sequence — a
//! 20-mer immediately 5′ of a PAM motif — then ranking the candidates
//! by predicted cutting efficiency.
//!
//! This module implements:
//!
//! - a configurable [`PamSpec`] (`NGG` for SpCas9, `NNGRRT` for SaCas9,
//!   `TTTV` for Cas12a, …) with IUPAC-aware matching and a 5′-vs-3′
//!   PAM side;
//! - [`scan_guides`] — a both-strand protospacer scan;
//! - an on-target efficiency score in the spirit of the Doench
//!   Rule-Set-2 model — a sequence-feature linear score over
//!   position-specific nucleotides, GC content and the PAM-proximal
//!   context.
//!
//! ## v1 scope
//!
//! The on-target score is a **transparent feature-weighted linear
//! model** — position-dependent nucleotide preferences, a GC-content
//! optimum and homopolymer / poly-T penalties, all in the *spirit* of
//! Doench Rule-Set-2 and the Moreno-Mateos rules. It is not the
//! published gradient-boosted-regression weights (that is a trained
//! model the project's "no llms / no trained weights" rule keeps out);
//! the score is documented as a heuristic and lands in `[0, 1]` with
//! the right qualitative ranking.

use crate::error::{GenomicsError, Result};

/// Which side of the protospacer the PAM sits on.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum PamSide {
    /// PAM is 3′ of the protospacer (SpCas9, SaCas9).
    ThreePrime,
    /// PAM is 5′ of the protospacer (Cas12a / Cpf1).
    FivePrime,
}

/// A PAM specification: the motif (IUPAC codes allowed), the
/// protospacer length and which side the PAM sits on.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PamSpec {
    /// The PAM motif, e.g. `"NGG"`. IUPAC ambiguity codes are honoured.
    pub motif: String,
    /// The protospacer (guide) length in bases.
    pub protospacer_len: usize,
    /// Which side of the protospacer the PAM is on.
    pub side: PamSide,
}

impl PamSpec {
    /// SpCas9 — `NGG` PAM, 20 nt protospacer, PAM 3′.
    pub fn spcas9() -> Self {
        PamSpec {
            motif: "NGG".to_string(),
            protospacer_len: 20,
            side: PamSide::ThreePrime,
        }
    }

    /// SaCas9 — `NNGRRT` PAM, 21 nt protospacer, PAM 3′.
    pub fn sacas9() -> Self {
        PamSpec {
            motif: "NNGRRT".to_string(),
            protospacer_len: 21,
            side: PamSide::ThreePrime,
        }
    }

    /// Cas12a / Cpf1 — `TTTV` PAM, 23 nt protospacer, PAM 5′.
    pub fn cas12a() -> Self {
        PamSpec {
            motif: "TTTV".to_string(),
            protospacer_len: 23,
            side: PamSide::FivePrime,
        }
    }

    /// The PAM length.
    pub fn pam_len(&self) -> usize {
        self.motif.len()
    }
}

/// `true` when nucleotide `base` is matched by IUPAC `code`.
pub fn iupac_match(code: u8, base: u8) -> bool {
    let code = code.to_ascii_uppercase();
    let base = base.to_ascii_uppercase();
    let set: &[u8] = match code {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"GC",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => return false,
    };
    set.contains(&base)
}

/// `true` when `window` matches the IUPAC `motif` base-for-base.
fn motif_matches(motif: &[u8], window: &[u8]) -> bool {
    motif.len() == window.len() && motif.iter().zip(window).all(|(&c, &b)| iupac_match(c, b))
}

/// Strand a candidate guide was found on.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum GuideStrand {
    /// The forward strand of the target sequence.
    Forward,
    /// The reverse-complement strand.
    Reverse,
}

/// One candidate guide RNA.
#[derive(Clone, Debug, PartialEq)]
pub struct Guide {
    /// The 20-mer (or N-mer) protospacer sequence, 5′→3′ on its strand.
    pub protospacer: String,
    /// The PAM sequence as found.
    pub pam: String,
    /// 0-based start of the protospacer on the *forward* target
    /// sequence.
    pub start: usize,
    /// Strand the guide was found on.
    pub strand: GuideStrand,
    /// On-target efficiency score in `[0, 1]` — higher is better.
    pub on_target_score: f64,
    /// GC fraction of the protospacer.
    pub gc_content: f64,
}

fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// GC fraction of a base slice.
fn gc_fraction(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
        .count();
    gc as f64 / seq.len() as f64
}

/// A Doench-Rule-Set-2-*style* on-target efficiency score in `[0, 1]`.
///
/// The score is a transparent feature-weighted heuristic, **not** the
/// trained model:
///
/// - **GC content** — peaks around 40-60 %; the score is penalised
///   away from that band (Doench: extreme GC hurts efficiency).
/// - **Position-specific nucleotides** — a `G` at the PAM-proximal
///   position 20 and a purine just upstream are favourable; a `T` at
///   PAM-proximal positions is penalised (the well-known Rule-Set-2
///   signs).
/// - **Poly-T** — a run of 4+ `T`s is a Pol-III terminator and zeroes
///   out a usable guide; it is heavily penalised.
///
/// `protospacer` is the guide 5′→3′; an empty guide scores `0`.
pub fn on_target_score(protospacer: &[u8]) -> f64 {
    if protospacer.is_empty() {
        return 0.0;
    }
    let g: Vec<u8> = protospacer.iter().map(|b| b.to_ascii_uppercase()).collect();
    let n = g.len();

    let mut score = 0.5f64; // neutral baseline

    // GC-content term: triangular peak at 0.5, zero contribution at the
    // 0.2 / 0.8 extremes.
    let gc = gc_fraction(&g);
    let gc_term = 1.0 - ((gc - 0.5).abs() / 0.3).min(1.0);
    score += 0.25 * (gc_term - 0.5);

    // Position-specific term (PAM-proximal = last base for a 3' PAM).
    // A favourable PAM-proximal G; an unfavourable PAM-proximal T.
    if let Some(&last) = g.last() {
        match last {
            b'G' => score += 0.08,
            b'C' => score += 0.03,
            b'T' => score -= 0.10,
            _ => {}
        }
    }
    // The base just upstream of the PAM-proximal position: a purine
    // helps.
    if n >= 2 {
        match g[n - 2] {
            b'A' | b'G' => score += 0.04,
            b'T' => score -= 0.04,
            _ => {}
        }
    }
    // PAM-distal (5') position: a `G` start is mildly favourable
    // (U6-promoter transcription).
    if g[0] == b'G' {
        score += 0.03;
    }

    // Poly-T terminator penalty: 4+ consecutive T anywhere.
    let mut run = 0usize;
    let mut max_t_run = 0usize;
    for &b in &g {
        if b == b'T' {
            run += 1;
            max_t_run = max_t_run.max(run);
        } else {
            run = 0;
        }
    }
    if max_t_run >= 4 {
        score -= 0.45;
    }

    // Homopolymer penalty for any 5+ run (folding / synthesis issues).
    let mut hp = 1usize;
    let mut max_hp = 1usize;
    for w in g.windows(2) {
        if w[0] == w[1] {
            hp += 1;
            max_hp = max_hp.max(hp);
        } else {
            hp = 1;
        }
    }
    if max_hp >= 5 {
        score -= 0.15;
    }

    score.clamp(0.0, 1.0)
}

/// Scans a target sequence on both strands for candidate guides.
///
/// Every position whose flanking window matches the [`PamSpec`] motif
/// yields a [`Guide`]; the reverse strand is scanned by reverse-
/// complementing. Each guide is scored with [`on_target_score`].
/// Returns the guides sorted by descending on-target score.
pub fn scan_guides(target: &[u8], pam: &PamSpec) -> Result<Vec<Guide>> {
    if pam.protospacer_len == 0 {
        return Err(GenomicsError::invalid(
            "protospacer_len",
            "must be positive",
        ));
    }
    if pam.motif.is_empty() {
        return Err(GenomicsError::invalid("motif", "PAM motif is empty"));
    }
    let total = pam.protospacer_len + pam.pam_len();
    if target.len() < total {
        return Ok(Vec::new());
    }
    let fwd: Vec<u8> = target.iter().map(|b| b.to_ascii_uppercase()).collect();
    let rev = revcomp(&fwd);
    let motif = pam.motif.as_bytes();
    let mut guides = Vec::new();

    // Forward strand.
    collect_strand(&fwd, motif, pam, GuideStrand::Forward, &fwd, &mut guides);
    // Reverse strand — the start is re-mapped onto the forward coords.
    collect_strand(&rev, motif, pam, GuideStrand::Reverse, &fwd, &mut guides);

    guides.sort_by(|a, b| {
        b.on_target_score
            .partial_cmp(&a.on_target_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(guides)
}

/// Collects guides from one already-oriented strand sequence.
fn collect_strand(
    seq: &[u8],
    motif: &[u8],
    pam: &PamSpec,
    strand: GuideStrand,
    forward_ref: &[u8],
    out: &mut Vec<Guide>,
) {
    let plen = pam.protospacer_len;
    let pamlen = pam.pam_len();
    let total = plen + pamlen;
    for i in 0..=seq.len() - total {
        let (proto, pam_seq): (&[u8], &[u8]) = match pam.side {
            PamSide::ThreePrime => (&seq[i..i + plen], &seq[i + plen..i + total]),
            PamSide::FivePrime => (&seq[i + pamlen..i + total], &seq[i..i + pamlen]),
        };
        if !motif_matches(motif, pam_seq) {
            continue;
        }
        // Skip guides containing ambiguous bases.
        if proto
            .iter()
            .any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T'))
        {
            continue;
        }
        // Map the protospacer start back to forward-strand coordinates.
        let fwd_start = match strand {
            GuideStrand::Forward => match pam.side {
                PamSide::ThreePrime => i,
                PamSide::FivePrime => i + pamlen,
            },
            GuideStrand::Reverse => {
                // `i` is the index in the reverse-complement string.
                let rc_proto_start = match pam.side {
                    PamSide::ThreePrime => i,
                    PamSide::FivePrime => i + pamlen,
                };
                // forward coordinate of the protospacer's forward-most
                // base.
                forward_ref.len() - rc_proto_start - plen
            }
        };
        out.push(Guide {
            protospacer: String::from_utf8_lossy(proto).into_owned(),
            pam: String::from_utf8_lossy(pam_seq).into_owned(),
            start: fwd_start,
            strand,
            on_target_score: on_target_score(proto),
            gc_content: gc_fraction(proto),
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iupac_matching() {
        assert!(iupac_match(b'N', b'A'));
        assert!(iupac_match(b'R', b'G'));
        assert!(!iupac_match(b'R', b'C'));
        assert!(iupac_match(b'V', b'C'));
        assert!(!iupac_match(b'V', b'T'));
    }

    #[test]
    fn finds_ngg_guide() {
        // 20-mer protospacer + NGG PAM.
        let proto = "ACGTACGTACGTACGTACGT";
        let target = format!("{proto}AGG"); // PAM = AGG matches NGG
        let guides = scan_guides(target.as_bytes(), &PamSpec::spcas9()).unwrap();
        let fwd: Vec<_> = guides
            .iter()
            .filter(|g| g.strand == GuideStrand::Forward)
            .collect();
        assert_eq!(fwd.len(), 1);
        assert_eq!(fwd[0].protospacer, proto);
        assert_eq!(fwd[0].pam, "AGG");
        assert_eq!(fwd[0].start, 0);
    }

    #[test]
    fn rejects_non_pam() {
        // PAM = AAA does NOT match NGG.
        let target = "ACGTACGTACGTACGTACGTAAA";
        let guides = scan_guides(target.as_bytes(), &PamSpec::spcas9()).unwrap();
        assert!(guides
            .iter()
            .all(|g| g.strand != GuideStrand::Forward || g.pam != "AAA"));
    }

    #[test]
    fn scans_reverse_strand() {
        // Place a guide only on the reverse strand: the forward
        // sequence ends with CCN read 5'->3', which is NGG on the
        // reverse strand.
        let target = "CCAACGTACGTACGTACGTACGT"; // CC at the 5' end
        let guides = scan_guides(target.as_bytes(), &PamSpec::spcas9()).unwrap();
        assert!(guides.iter().any(|g| g.strand == GuideStrand::Reverse));
    }

    #[test]
    fn poly_t_guide_scores_low() {
        // A guide with a 4-T run (Pol-III terminator).
        let with_polyt = b"ACGTTTTACGTACGTACGTA";
        let clean = b"ACGTACGTACGTACGTACGT";
        assert!(on_target_score(with_polyt) < on_target_score(clean));
    }

    #[test]
    fn extreme_gc_scores_lower_than_balanced() {
        let all_gc = b"GCGCGCGCGCGCGCGCGCGC";
        let balanced = b"ACGTACGTACGTACGTACGT";
        assert!(on_target_score(all_gc) < on_target_score(balanced));
    }

    #[test]
    fn score_in_unit_range() {
        for g in [
            b"ACGTACGTACGTACGTACGT".as_slice(),
            b"TTTTTTTTTTTTTTTTTTTT".as_slice(),
            b"GGGGGGGGGGGGGGGGGGGG".as_slice(),
        ] {
            let s = on_target_score(g);
            assert!((0.0..=1.0).contains(&s), "score {s} out of range");
        }
    }

    #[test]
    fn guides_sorted_by_score() {
        // A target with several PAM sites.
        let target = "ACGTACGTACGTACGTACGTAGGACGTACGTACGTACGTACGTCGGACGTACGTACGTACGTACGTTGG";
        let guides = scan_guides(target.as_bytes(), &PamSpec::spcas9()).unwrap();
        for w in guides.windows(2) {
            assert!(w[0].on_target_score >= w[1].on_target_score);
        }
    }

    #[test]
    fn short_target_yields_nothing() {
        let guides = scan_guides(b"ACGT", &PamSpec::spcas9()).unwrap();
        assert!(guides.is_empty());
    }

    #[test]
    fn cas12a_five_prime_pam() {
        // Cas12a: TTTV PAM 5' of a 23-mer protospacer.
        let proto = "ACGTACGTACGTACGTACGTACG"; // 23
        let target = format!("TTTA{proto}"); // TTTA matches TTTV
        let guides = scan_guides(target.as_bytes(), &PamSpec::cas12a()).unwrap();
        let fwd: Vec<_> = guides
            .iter()
            .filter(|g| g.strand == GuideStrand::Forward)
            .collect();
        assert_eq!(fwd.len(), 1);
        assert_eq!(fwd[0].protospacer, proto);
        assert_eq!(fwd[0].start, 4);
    }
}
