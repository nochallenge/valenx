//! Off-target enumeration — mismatch-tolerant genome scan + CFD score.
//!
//! A guide RNA can cut at genomic sites that differ from the intended
//! target by a few mismatches. Cas-OFFinder enumerates every such site
//! by scanning the genome for protospacer-PAM windows within a
//! mismatch budget; CRISPOR / the Doench lab score each site with the
//! **Cutting Frequency Determination (CFD)** matrix — a
//! position-and-mismatch-specific activity weight.
//!
//! This module implements both:
//!
//! - [`enumerate_off_targets`] — a both-strand genome scan that finds
//!   every protospacer window within `max_mismatches` of the guide and
//!   adjacent to a valid PAM;
//! - [`cfd_score`] — a CFD-style off-target activity score, the
//!   product of per-position mismatch penalties.
//!
//! ## v1 scope
//!
//! The CFD penalty table is a **structured position-weighted model**:
//! a mismatch is penalised more the closer it sits to the PAM (the
//! universal CRISPR specificity gradient), with a milder penalty for
//! transition-like mismatches. It is the *shape* of the published CFD
//! matrix, not the exact ~400 trained coefficients (those are a
//! trained model the project rule excludes); it is documented as a
//! heuristic. The scan is a direct O(genome × guide) sweep — correct
//! and exhaustive within the mismatch budget, not a seed-indexed
//! accelerated search. Bulges (RNA / DNA insertions) are not modelled.

use crate::crispr::guide::{iupac_match, PamSide, PamSpec};
use crate::error::{GenomicsError, Result};

/// One enumerated off-target site.
#[derive(Clone, Debug, PartialEq)]
pub struct OffTarget {
    /// The genomic contig the site is on.
    pub chrom: String,
    /// 0-based start of the protospacer on the forward strand.
    pub start: usize,
    /// `true` when the match is on the reverse-complement strand.
    pub reverse: bool,
    /// The protospacer sequence as found at the site (5′→3′ on its
    /// strand).
    pub protospacer: String,
    /// The PAM sequence at the site.
    pub pam: String,
    /// Number of mismatches between the guide and this protospacer.
    pub mismatches: usize,
    /// 0-based positions (in guide coordinates, 5′→3′) of the
    /// mismatches.
    pub mismatch_positions: Vec<usize>,
    /// CFD-style off-target activity score in `[0, 1]`; `1.0` is a
    /// perfect match, lower means a less active off-target.
    pub cfd_score: f64,
}

impl OffTarget {
    /// `true` when this is a perfect (zero-mismatch) match — the
    /// on-target site itself when scanning the source genome.
    pub fn is_perfect(&self) -> bool {
        self.mismatches == 0
    }
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

/// `true` when `window` matches the IUPAC `motif` base for base.
fn pam_ok(motif: &[u8], window: &[u8]) -> bool {
    motif.len() == window.len() && motif.iter().zip(window).all(|(&c, &b)| iupac_match(c, b))
}

/// A CFD-style off-target activity score in `[0, 1]`.
///
/// The score is the product over guide positions of a per-position
/// "retained activity" factor:
///
/// - a matched position contributes `1.0`;
/// - a mismatched position contributes a penalty that **rises toward
///   the PAM** — a mismatch in the PAM-distal seed-free region barely
///   hurts, a PAM-proximal mismatch (the "seed") nearly abolishes
///   cutting;
/// - the score is also multiplied by a PAM-activity factor (`NGG`
///   active = 1.0, `NAG` ≈ 0.2, others low).
///
/// `guide` and `protospacer` must be the same length, both 5′→3′;
/// `pam` is the off-target PAM. For a 3′-PAM nuclease, guide position
/// 0 is PAM-distal and the last position is PAM-proximal.
pub fn cfd_score(guide: &[u8], protospacer: &[u8], pam: &[u8]) -> f64 {
    if guide.is_empty() || guide.len() != protospacer.len() {
        return 0.0;
    }
    let n = guide.len();
    let mut score = 1.0f64;
    for i in 0..n {
        let g = guide[i].to_ascii_uppercase();
        let p = protospacer[i].to_ascii_uppercase();
        if g == p {
            continue;
        }
        // PAM-distal fraction → [0, 1]. For a 3′-PAM nuclease the guide
        // is 5′→3′ with the PAM at the 3′ end, so position 0 is the
        // PAM-*distal* 5′ end (→ 1.0) and the last position is the
        // PAM-proximal seed (→ 0.0).
        let pam_distal = 1.0 - i as f64 / (n - 1).max(1) as f64;
        // Retained activity: a PAM-distal mismatch keeps ~0.85 of the
        // activity, a PAM-proximal (seed) one keeps ~0.05.
        let mut retained = 0.05 + 0.80 * pam_distal;
        // Milder penalty for a transition-like mismatch (rG:dT, rU:dG
        // wobble pairs tolerate better).
        if is_wobble_like(g, p) {
            retained = (retained + 0.15).min(1.0);
        }
        score *= retained;
    }
    // PAM-activity factor — applied for an NGG-class PAM.
    score *= pam_activity(pam);
    score.clamp(0.0, 1.0)
}

/// `true` for a guide:DNA mismatch that behaves like a tolerated
/// wobble (rG:dT and rU:dG — i.e. guide `G` vs target `T`, guide `T`
/// vs target `G`, since the guide RNA `U` is written `T`).
fn is_wobble_like(guide: u8, target: u8) -> bool {
    matches!((guide, target), (b'G', b'T') | (b'T', b'G'))
}

/// A PAM-activity multiplier for an NGG-class nuclease: `NGG` cuts
/// fully, `NAG` is a weak alternative PAM, anything else is near-dead.
fn pam_activity(pam: &[u8]) -> f64 {
    if pam.len() < 3 {
        return 1.0; // unknown PAM length — do not penalise
    }
    let p: Vec<u8> = pam
        .iter()
        .rev()
        .take(2)
        .map(|b| b.to_ascii_uppercase())
        .collect();
    // p[0] is the last base, p[1] the second-to-last.
    match (p[1], p[0]) {
        (b'G', b'G') => 1.0,
        (b'A', b'G') => 0.20,
        (b'G', b'A') => 0.10,
        _ => 0.03,
    }
}

/// Enumerates off-target sites for a guide across a set of named
/// contigs.
///
/// `guide` is the 5′→3′ protospacer the user wants. Every contig is
/// scanned on both strands; a window is an off-target hit when its
/// protospacer is within `max_mismatches` of the guide **and** the
/// adjacent PAM matches `pam.motif`. Each hit is CFD-scored. Results
/// are sorted by descending CFD score (the most dangerous off-targets
/// first).
pub fn enumerate_off_targets(
    guide: &[u8],
    genome: &[(String, Vec<u8>)],
    pam: &PamSpec,
    max_mismatches: usize,
) -> Result<Vec<OffTarget>> {
    if guide.is_empty() {
        return Err(GenomicsError::invalid("guide", "guide is empty"));
    }
    if guide.len() != pam.protospacer_len {
        return Err(GenomicsError::invalid(
            "guide",
            format!(
                "guide length {} != PAM-spec protospacer length {}",
                guide.len(),
                pam.protospacer_len
            ),
        ));
    }
    let guide_u: Vec<u8> = guide.iter().map(|b| b.to_ascii_uppercase()).collect();
    let motif = pam.motif.as_bytes();
    let plen = pam.protospacer_len;
    let pamlen = pam.pam_len();
    let total = plen + pamlen;
    let mut hits = Vec::new();

    for (name, seq) in genome {
        let fwd: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        if fwd.len() < total {
            continue;
        }
        let rev = revcomp(&fwd);

        // Forward strand.
        scan_one(
            &fwd,
            &guide_u,
            motif,
            plen,
            pamlen,
            total,
            pam.side,
            name,
            false,
            &fwd,
            max_mismatches,
            &mut hits,
        );
        // Reverse strand.
        scan_one(
            &rev,
            &guide_u,
            motif,
            plen,
            pamlen,
            total,
            pam.side,
            name,
            true,
            &fwd,
            max_mismatches,
            &mut hits,
        );
    }

    hits.sort_by(|a, b| {
        b.cfd_score
            .partial_cmp(&a.cfd_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(hits)
}

#[allow(clippy::too_many_arguments)]
fn scan_one(
    seq: &[u8],
    guide: &[u8],
    motif: &[u8],
    plen: usize,
    pamlen: usize,
    total: usize,
    side: PamSide,
    chrom: &str,
    reverse: bool,
    forward_ref: &[u8],
    max_mm: usize,
    out: &mut Vec<OffTarget>,
) {
    for i in 0..=seq.len() - total {
        let (proto, pam_seq): (&[u8], &[u8]) = match side {
            PamSide::ThreePrime => (&seq[i..i + plen], &seq[i + plen..i + total]),
            PamSide::FivePrime => (&seq[i + pamlen..i + total], &seq[i..i + pamlen]),
        };
        if !pam_ok(motif, pam_seq) {
            continue;
        }
        // Count mismatches, bailing once the budget is blown.
        let mut mm = 0usize;
        let mut positions = Vec::new();
        let mut over = false;
        for (k, (&g, &p)) in guide.iter().zip(proto).enumerate() {
            if g != p {
                mm += 1;
                positions.push(k);
                if mm > max_mm {
                    over = true;
                    break;
                }
            }
        }
        if over {
            continue;
        }
        // Map back to forward-strand coordinates.
        let fwd_start = if !reverse {
            match side {
                PamSide::ThreePrime => i,
                PamSide::FivePrime => i + pamlen,
            }
        } else {
            let rc_proto_start = match side {
                PamSide::ThreePrime => i,
                PamSide::FivePrime => i + pamlen,
            };
            forward_ref.len() - rc_proto_start - plen
        };
        let cfd = cfd_score(guide, proto, pam_seq);
        out.push(OffTarget {
            chrom: chrom.to_string(),
            start: fwd_start,
            reverse,
            protospacer: String::from_utf8_lossy(proto).into_owned(),
            pam: String::from_utf8_lossy(pam_seq).into_owned(),
            mismatches: mm,
            mismatch_positions: positions,
            cfd_score: cfd,
        });
    }
}

/// Aggregate specificity score for a guide given its off-target set —
/// the CRISPOR-style guide score `100 / (100 + Σ cfd_off)`, where the
/// sum runs over the **non-perfect** off-targets. A guide with no
/// off-targets scores `1.0`; many active off-targets drive it toward
/// `0`.
pub fn guide_specificity_score(off_targets: &[OffTarget]) -> f64 {
    let sum: f64 = off_targets
        .iter()
        .filter(|o| !o.is_perfect())
        .map(|o| o.cfd_score)
        .sum();
    100.0 / (100.0 + 100.0 * sum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crispr::guide::PamSpec;

    fn genome(seqs: &[(&str, &str)]) -> Vec<(String, Vec<u8>)> {
        seqs.iter()
            .map(|(n, s)| (n.to_string(), s.as_bytes().to_vec()))
            .collect()
    }

    #[test]
    fn perfect_match_scores_one() {
        let guide = b"ACGTACGTACGTACGTACGT";
        let s = cfd_score(guide, guide, b"AGG");
        assert!((s - 1.0).abs() < 1e-9);
    }

    #[test]
    fn pam_proximal_mismatch_hurts_more() {
        let guide = b"ACGTACGTACGTACGTACGT";
        // Mismatch at position 0 (PAM-distal).
        let mut distal = guide.to_vec();
        distal[0] = b'T';
        // Mismatch at position 19 (PAM-proximal).
        let mut proximal = guide.to_vec();
        proximal[19] = b'A';
        let s_distal = cfd_score(guide, &distal, b"AGG");
        let s_proximal = cfd_score(guide, &proximal, b"AGG");
        assert!(
            s_proximal < s_distal,
            "proximal {s_proximal} should be < distal {s_distal}"
        );
    }

    #[test]
    fn nag_pam_scores_lower_than_ngg() {
        let guide = b"ACGTACGTACGTACGTACGT";
        let ngg = cfd_score(guide, guide, b"AGG");
        let nag = cfd_score(guide, guide, b"AAG");
        assert!(nag < ngg);
    }

    #[test]
    fn finds_perfect_on_target() {
        let proto = "ACGTACGTACGTACGTACGT";
        let g = genome(&[("chr1", &format!("{proto}AGG"))]);
        let hits = enumerate_off_targets(proto.as_bytes(), &g, &PamSpec::spcas9(), 3).unwrap();
        let perfect: Vec<_> = hits.iter().filter(|h| h.is_perfect()).collect();
        assert_eq!(perfect.len(), 1);
        assert_eq!(perfect[0].chrom, "chr1");
        assert_eq!(perfect[0].mismatches, 0);
    }

    #[test]
    fn finds_mismatched_off_target() {
        let guide = "ACGTACGTACGTACGTACGT";
        // An off-target with 2 mismatches followed by a TGG PAM.
        let off = "TCGTACGTACGTACGTACGA"; // pos 0 and pos 19 differ
        let g = genome(&[("chr1", &format!("{off}TGG"))]);
        let hits = enumerate_off_targets(guide.as_bytes(), &g, &PamSpec::spcas9(), 3).unwrap();
        let two_mm: Vec<_> = hits.iter().filter(|h| h.mismatches == 2).collect();
        assert!(!two_mm.is_empty());
        assert_eq!(two_mm[0].mismatch_positions, vec![0, 19]);
    }

    #[test]
    fn mismatch_budget_excludes_distant_sites() {
        let guide = "ACGTACGTACGTACGTACGT";
        // A site with 5 mismatches — over a budget of 2.
        let off = "TTTTTCGTACGTACGTACGT";
        let g = genome(&[("chr1", &format!("{off}AGG"))]);
        let hits = enumerate_off_targets(guide.as_bytes(), &g, &PamSpec::spcas9(), 2).unwrap();
        assert!(hits.iter().all(|h| h.mismatches <= 2));
    }

    #[test]
    fn rejects_length_mismatch() {
        let g = genome(&[("chr1", "ACGT")]);
        // 4-base guide vs a 20-base PAM spec.
        assert!(enumerate_off_targets(b"ACGT", &g, &PamSpec::spcas9(), 3).is_err());
    }

    #[test]
    fn specificity_score_high_with_few_off_targets() {
        let perfect = OffTarget {
            chrom: "c".to_string(),
            start: 0,
            reverse: false,
            protospacer: "A".to_string(),
            pam: "AGG".to_string(),
            mismatches: 0,
            mismatch_positions: vec![],
            cfd_score: 1.0,
        };
        // Only the perfect on-target -> specificity 1.0.
        assert!((guide_specificity_score(&[perfect.clone()]) - 1.0).abs() < 1e-9);

        let weak_off = OffTarget {
            mismatches: 3,
            cfd_score: 0.01,
            ..perfect.clone()
        };
        let s = guide_specificity_score(&[perfect, weak_off]);
        assert!(s < 1.0 && s > 0.9);
    }

    #[test]
    fn hits_sorted_by_cfd() {
        let guide = "ACGTACGTACGTACGTACGT";
        // Two PAM sites: a perfect match and a 2-mismatch site.
        let perfect = format!("{guide}AGG");
        let off = format!("TCGTACGTACGTACGTACGA{}", "TGG");
        let g = genome(&[("chr1", &format!("{perfect}{off}"))]);
        let hits = enumerate_off_targets(guide.as_bytes(), &g, &PamSpec::spcas9(), 3).unwrap();
        for w in hits.windows(2) {
            assert!(w[0].cfd_score >= w[1].cfd_score);
        }
    }
}
