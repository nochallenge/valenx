//! Features 18–19 — 5′UTR and 3′UTR design.
//!
//! The untranslated regions flank the CDS and tune translation and
//! stability:
//!
//! - the **5′UTR** carries the **Kozak context** around the start
//!   codon — the consensus `gccRccATGG` (a purine at −3 and a `G` at
//!   +4 matter most) — and should be short, unstructured and free of
//!   upstream `AUG`s (uORFs);
//! - the **3′UTR** carries **stability elements** (the well-known
//!   human α-/β-globin 3′UTRs are the workhorse stabilisers) and must
//!   **avoid AU-rich elements (AREs)** — `AUUUA` pentamers that
//!   recruit the mRNA-decay machinery.
//!
//! This module **scores** a 5′UTR for Kozak quality and uORF burden
//! ([`analyze_utr5`]), **scores** a 3′UTR for ARE load
//! ([`analyze_utr3`]), and offers **reference UTRs** known to express
//! well ([`reference_utr5`], [`reference_utr3`]).
//!
//! ## v1 scope
//!
//! The Kozak score is a transparent position-weighted match to the
//! consensus, not a trained translation-initiation model. The 3′UTR
//! analysis counts AREs and known destabilising motifs; it does not
//! model the full miRNA-target / RBP-site landscape.

use crate::error::{GeneditingError, Result};
use crate::sequtil::{is_acgu, transcribe};
use serde::{Deserialize, Serialize};

/// Quality analysis of a 5′UTR (feature 18).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Utr5Analysis {
    /// Kozak-context score in `[0, 1]` — how well the bases around the
    /// start codon match the `gccRccATGG` consensus.
    pub kozak_score: f64,
    /// `true` when the strongest determinants are met — a purine
    /// (`A`/`G`) at position −3 and a `G` at position +4.
    pub strong_kozak: bool,
    /// Number of upstream `AUG`s in the 5′UTR (uORF starts — each one
    /// can sequester scanning ribosomes).
    pub uorf_count: usize,
    /// 5′UTR length in nucleotides.
    pub length: usize,
    /// `true` when the UTR is in the recommended length band
    /// (~20–100 nt).
    pub length_ok: bool,
    /// A combined design verdict in `[0, 1]` — Kozak quality, no
    /// uORFs, sensible length.
    pub design_score: f64,
}

/// The Kozak consensus around an `AUG`, positions −6..+4 (the start
/// codon's `A` is position +1). Upper-case = strongly conserved.
/// `gcc gcc A U G g` — here as the per-position preferred base, with
/// `R` (purine) at −3.
const KOZAK_CONSENSUS: &[(i32, u8)] = &[
    (-6, b'G'),
    (-5, b'C'),
    (-4, b'C'),
    (-3, b'R'), // purine — the single most important position
    (-2, b'C'),
    (-1, b'C'),
    (4, b'G'), // +4 — the second most important position
];

/// `true` when RNA base `base` is matched by IUPAC `code` (a tiny
/// subset — `R` purine and the four plain bases are all we need here).
fn kozak_match(code: u8, base: u8) -> bool {
    let b = base.to_ascii_uppercase();
    match code {
        b'R' => matches!(b, b'A' | b'G'),
        b'A' => b == b'A',
        b'C' => b == b'C',
        b'G' => b == b'G',
        b'U' => b == b'U',
        _ => false,
    }
}

/// Counts non-overlapping `AUG` triplets in a sequence (uORF starts).
fn count_aug(seq: &[u8]) -> usize {
    let mut n = 0;
    let mut i = 0;
    while i + 3 <= seq.len() {
        if seq[i..i + 3].eq_ignore_ascii_case(b"AUG") {
            n += 1;
            i += 3;
        } else {
            i += 1;
        }
    }
    n
}

/// Analyses a 5′UTR for Kozak quality and uORF burden (feature 18).
///
/// `utr5` is the 5′UTR (DNA or RNA — transcribed); `cds` is the coding
/// sequence whose start codon the Kozak context wraps. The Kozak score
/// weighs each consensus position, with the −3 and +4 positions
/// dominating.
///
/// # Errors
/// - [`GeneditingError::InvalidTarget`] for a non-ACGU UTR or a CDS too
///   short to read the `+4` position.
pub fn analyze_utr5(utr5: &[u8], cds: &[u8]) -> Result<Utr5Analysis> {
    let utr = transcribe(utr5);
    if !utr.is_empty() && !is_acgu(&utr) {
        return Err(GeneditingError::invalid_target(
            "region",
            "5'UTR must be A/C/G/U after transcription",
        ));
    }
    let cds_rna = transcribe(cds);
    if cds_rna.len() < 4 {
        return Err(GeneditingError::invalid_target(
            "cds",
            "CDS too short to evaluate the Kozak +4 position",
        ));
    }
    // Build the context window: 5'UTR bases (the last 6) + start codon
    // + the CDS base at +4.
    //   position +1 = cds_rna[0]; +4 = cds_rna[3].
    let mut total_weight = 0.0f64;
    let mut score = 0.0f64;
    let mut purine_minus3 = false;
    let mut g_plus4 = false;
    for &(pos, code) in KOZAK_CONSENSUS {
        // Weight: -3 and +4 are the dominant determinants.
        let weight = if pos == -3 || pos == 4 { 3.0 } else { 1.0 };
        total_weight += weight;
        let base: Option<u8> = if pos < 0 {
            // Position -k is the k-th base back from the start codon —
            // i.e. utr[len - k].
            let k = (-pos) as usize;
            utr.len().checked_sub(k).map(|i| utr[i])
        } else {
            // pos == 4 → cds_rna[3].
            cds_rna.get((pos - 1) as usize).copied()
        };
        if let Some(b) = base {
            if kozak_match(code, b) {
                score += weight;
                if pos == -3 {
                    purine_minus3 = true;
                }
                if pos == 4 {
                    g_plus4 = true;
                }
            }
        }
    }
    let kozak_score = if total_weight > 0.0 {
        score / total_weight
    } else {
        0.0
    };
    let uorf_count = count_aug(&utr);
    let length = utr.len();
    let length_ok = (20..=100).contains(&length);

    // Design score: Kozak quality, a uORF penalty, a length bonus.
    let uorf_penalty = 0.20 * uorf_count as f64;
    let length_bonus = if length_ok { 0.10 } else { 0.0 };
    let design_score = (kozak_score + length_bonus - uorf_penalty).clamp(0.0, 1.0);

    Ok(Utr5Analysis {
        kozak_score,
        strong_kozak: purine_minus3 && g_plus4,
        uorf_count,
        length,
        length_ok,
        design_score,
    })
}

/// One AU-rich / destabilising element found in a 3′UTR.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DestabilizingElement {
    /// 0-based start of the motif in the 3′UTR.
    pub pos: usize,
    /// The motif kind (`"ARE pentamer"`, `"ARE nonamer"`,
    /// `"GU-rich"`).
    pub kind: String,
}

/// Quality analysis of a 3′UTR (feature 19).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Utr3Analysis {
    /// AU-rich and other destabilising elements found.
    pub destabilizers: Vec<DestabilizingElement>,
    /// 3′UTR length in nucleotides.
    pub length: usize,
    /// AU fraction of the 3′UTR.
    pub au_fraction: f64,
    /// A combined stability verdict in `[0, 1]` — higher = a more
    /// stable, ARE-free 3′UTR.
    pub stability_score: f64,
}

impl Utr3Analysis {
    /// Number of destabilising elements found.
    pub fn destabilizer_count(&self) -> usize {
        self.destabilizers.len()
    }

    /// `true` when no ARE / destabilising motif was found.
    pub fn is_clean(&self) -> bool {
        self.destabilizers.is_empty()
    }
}

/// Analyses a 3′UTR for AU-rich elements and stability (feature 19).
///
/// Scans for `AUUUA` ARE pentamers, the `UUAUUUAUU`-class ARE nonamer
/// and long `GU`-rich stretches; reports the AU fraction and a
/// stability score that penalises every destabiliser and an extreme
/// AU content.
///
/// # Errors
/// [`GeneditingError::InvalidTarget`] for a non-ACGU 3′UTR.
pub fn analyze_utr3(utr3: &[u8]) -> Result<Utr3Analysis> {
    let utr = transcribe(utr3);
    if !utr.is_empty() && !is_acgu(&utr) {
        return Err(GeneditingError::invalid_target(
            "region",
            "3'UTR must be A/C/G/U after transcription",
        ));
    }
    let mut destabilizers = Vec::new();
    // ARE pentamer AUUUA.
    for i in 0..utr.len().saturating_sub(4) {
        if utr[i..i + 5].eq_ignore_ascii_case(b"AUUUA") {
            destabilizers.push(DestabilizingElement {
                pos: i,
                kind: "ARE pentamer (AUUUA)".to_string(),
            });
        }
    }
    // ARE nonamer UUAUUUAUU (a high-affinity ARE).
    for i in 0..utr.len().saturating_sub(8) {
        if utr[i..i + 9].eq_ignore_ascii_case(b"UUAUUUAUU") {
            destabilizers.push(DestabilizingElement {
                pos: i,
                kind: "ARE nonamer (UUAUUUAUU)".to_string(),
            });
        }
    }
    // A long GU-rich stretch (>= 10 nt of only G/U) — a GU-rich
    // destabilising element.
    let mut run = 0usize;
    let mut run_start = 0usize;
    for (i, &b) in utr.iter().enumerate() {
        if matches!(b.to_ascii_uppercase(), b'G' | b'U') {
            if run == 0 {
                run_start = i;
            }
            run += 1;
        } else {
            if run >= 10 {
                destabilizers.push(DestabilizingElement {
                    pos: run_start,
                    kind: "GU-rich element".to_string(),
                });
            }
            run = 0;
        }
    }
    if run >= 10 {
        destabilizers.push(DestabilizingElement {
            pos: run_start,
            kind: "GU-rich element".to_string(),
        });
    }

    let length = utr.len();
    let au = utr
        .iter()
        .filter(|&&b| matches!(b.to_ascii_uppercase(), b'A' | b'U'))
        .count();
    let au_fraction = if length == 0 {
        0.0
    } else {
        au as f64 / length as f64
    };
    // Stability: 1.0 minus a penalty per destabiliser and an extreme-AU
    // penalty (a very AT-rich UTR decays faster even without a discrete
    // ARE).
    let destab_penalty = 0.18 * destabilizers.len() as f64;
    let au_penalty = if au_fraction > 0.65 {
        0.20 * (au_fraction - 0.65) / 0.35
    } else {
        0.0
    };
    let stability_score = (1.0 - destab_penalty - au_penalty).clamp(0.0, 1.0);

    Ok(Utr3Analysis {
        destabilizers,
        length,
        au_fraction,
        stability_score,
    })
}

/// A reference 5′UTR known to express well — a short, unstructured,
/// uORF-free leader with a strong Kozak context (an HBA/HBB-style
/// minimal leader). RNA.
pub fn reference_utr5() -> Vec<u8> {
    // Short GC-balanced leader ending in a strong Kozak context
    // (...GCCACC before the ATG, purine at -3, the +4 G comes from the
    // CDS).
    b"GGGAAAUAAGAGAGAAAAGAAGAGUAAGAAGAAAUAUAAGAGCCACC".to_vec()
}

/// A reference 3′UTR known to stabilise mRNA — a human β-globin-style
/// 3′UTR, ARE-free. RNA.
pub fn reference_utr3() -> Vec<u8> {
    // Representative GC-balanced, ARE-free trailer (the human HBB 3'UTR
    // is the canonical stabiliser; this is a clean stand-in).
    b"GCUCGCUUUCUUGCUGUCCAAUUUCUAUUAAAGGUUCCUUUGUUCCCUAAGUCCAACUACUAAACUGGG".to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strong_kozak_scores_high() {
        // 5'UTR ending GCCACC, CDS = ATGG... → purine at -3 (A), G at +4.
        let a = analyze_utr5(b"AAAAAAAAAAAAAAAAAAAAGCCACC", b"ATGGCCTAA").unwrap();
        assert!(a.strong_kozak);
        assert!(a.kozak_score > 0.7);
    }

    #[test]
    fn weak_kozak_scores_low() {
        // 5'UTR ending in U at -3 (a pyrimidine), CDS = ATGA → no +4 G.
        let a = analyze_utr5(b"AAAAAAAAAAAAAAAAAAAAUUUUUU", b"ATGATATAA").unwrap();
        assert!(!a.strong_kozak);
        // Strong Kozak should outscore this.
        let strong = analyze_utr5(b"AAAAAAAAAAAAAAAAAAAAGCCACC", b"ATGGCCTAA").unwrap();
        assert!(strong.kozak_score > a.kozak_score);
    }

    #[test]
    fn detects_uorfs() {
        // A 5'UTR containing two AUGs.
        let a = analyze_utr5(b"AUGCCCAUGCCCGCCACCGCCACC", b"ATGGCCTAA").unwrap();
        assert_eq!(a.uorf_count, 2);
        // A uORF-free UTR scores higher.
        let clean = analyze_utr5(b"CCCCCCCCCCCCCCCCCCCCGCCACC", b"ATGGCCTAA").unwrap();
        assert!(clean.design_score > a.design_score);
    }

    #[test]
    fn rejects_non_acgu_utr5() {
        assert!(analyze_utr5(b"NNNN", b"ATGGCCTAA").is_err());
    }

    #[test]
    fn rejects_short_cds() {
        assert!(analyze_utr5(b"GCCACC", b"AT").is_err());
    }

    #[test]
    fn length_band_flag() {
        let short = analyze_utr5(b"GCCACC", b"ATGGCCTAA").unwrap();
        assert!(!short.length_ok); // 6 nt < 20
        let ok = analyze_utr5(b"AAAAAAAAAAAAAAAAAAAAGCCACC", b"ATGGCCTAA").unwrap();
        assert!(ok.length_ok); // 26 nt
    }

    #[test]
    fn detects_are_pentamer() {
        let a = analyze_utr3(b"GGGGAUUUAGGGG").unwrap();
        assert!(a
            .destabilizers
            .iter()
            .any(|d| d.kind.contains("pentamer")));
        assert!(!a.is_clean());
    }

    #[test]
    fn detects_are_nonamer() {
        let a = analyze_utr3(b"GGGUUAUUUAUUGGG").unwrap();
        assert!(a.destabilizers.iter().any(|d| d.kind.contains("nonamer")));
    }

    #[test]
    fn detects_gu_rich_element() {
        let a = analyze_utr3(b"CCCGUGUGUGUGUGUCCC").unwrap();
        assert!(a.destabilizers.iter().any(|d| d.kind.contains("GU-rich")));
    }

    #[test]
    fn clean_utr3_is_stable() {
        let a = analyze_utr3(b"GCGCGCGCGCGCGCGC").unwrap();
        assert!(a.is_clean());
        assert!(a.stability_score > 0.8);
    }

    #[test]
    fn are_load_lowers_stability() {
        let clean = analyze_utr3(b"GCGCGCGCGCGCGCGC").unwrap();
        let aring = analyze_utr3(b"AUUUAGGGAUUUAGGGAUUUA").unwrap();
        assert!(clean.stability_score > aring.stability_score);
    }

    #[test]
    fn reference_utr3_is_are_free() {
        let a = analyze_utr3(&reference_utr3()).unwrap();
        assert!(a.is_clean(), "the reference 3'UTR must be ARE-free");
    }

    #[test]
    fn reference_utr5_has_a_strong_kozak() {
        let a = analyze_utr5(&reference_utr5(), b"ATGGCCTAA").unwrap();
        assert!(a.strong_kozak);
    }
}
