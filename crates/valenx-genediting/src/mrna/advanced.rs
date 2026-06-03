//! Feature 23 — self-amplifying mRNA (saRNA) and circular-RNA
//! construct design (v1).
//!
//! Two newer mRNA formats trade construct complexity for performance:
//!
//! - **Self-amplifying mRNA (saRNA)** — derived from an alphavirus
//!   replicon. The construct encodes the viral **replicase (nsP1-4)**
//!   *and* the gene of interest behind a **subgenomic promoter**;
//!   once in the cell the replicase copies the RNA, so a tiny dose
//!   gives prolonged, amplified expression. The construct still needs
//!   a cap, the conserved sequence elements and a poly-A tail.
//! - **Circular RNA (circRNA)** — a covalently closed circle with no
//!   free ends, so exonucleases cannot degrade it. It has no cap and
//!   no poly-A; translation is driven by an **IRES** (internal
//!   ribosome entry site). The construct is laid out for back-splicing
//!   or a permuted-intron-exon (PIE) ligation.
//!
//! This module assembles **construct layouts** for both — the part
//! order, the sizes, and validation that the gene of interest is a
//! proper CDS. It does not simulate replicase kinetics or splicing.
//!
//! ## v1 scope
//!
//! These are v1 *construct-layout* designers — they order and validate
//! the parts and report the total size and a layout description. They
//! do not model alphavirus replication, IRES strength or
//! back-splicing efficiency; the replicase and IRES are represented as
//! sized, named elements supplied by the caller (or as documented
//! placeholders).

use crate::error::{GeneditingError, Result};
use crate::sequtil::{is_acgu, transcribe, upper};
use serde::{Deserialize, Serialize};
use valenx_bioseq::ops::translate::GeneticCode;

/// Validates that `cds` is a proper coding sequence (RNA), returning
/// the transcribed RNA. Shared by both designers here.
fn validate_cds(cds: &[u8]) -> Result<Vec<u8>> {
    let rna = transcribe(cds);
    if rna.is_empty() || rna.len() % 3 != 0 || !is_acgu(&rna) {
        return Err(GeneditingError::invalid_target(
            "cds",
            "gene of interest must be a non-empty A/C/G/U CDS of length divisible by 3",
        ));
    }
    let code = GeneticCode::standard();
    let first = dna_codon(&rna[0..3]);
    if !code.is_start_codon(&first) {
        return Err(GeneditingError::invalid_target(
            "cds",
            "gene of interest does not begin with a start codon",
        ));
    }
    let last = dna_codon(&rna[rna.len() - 3..]);
    if !code.is_stop_codon(&last) {
        return Err(GeneditingError::invalid_target(
            "cds",
            "gene of interest does not end with a stop codon",
        ));
    }
    Ok(rna)
}

/// Reverse-transcribes a 3-base RNA codon to a DNA codon array.
fn dna_codon(rna: &[u8]) -> [u8; 3] {
    let mut out = [b'N'; 3];
    for (i, &b) in rna.iter().take(3).enumerate() {
        let u = b.to_ascii_uppercase();
        out[i] = if u == b'U' { b'T' } else { u };
    }
    out
}

/// A self-amplifying-mRNA (saRNA) construct layout.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SaRnaConstruct {
    /// The alphavirus replicase (nsP1-4) coding region (RNA). When the
    /// caller has no replicase sequence, a documented placeholder of
    /// the typical ~7.5 kb length is used.
    pub replicase: Vec<u8>,
    /// The subgenomic promoter element separating the replicase from
    /// the gene of interest (RNA).
    pub subgenomic_promoter: Vec<u8>,
    /// The gene-of-interest CDS (RNA).
    pub gene_of_interest: Vec<u8>,
    /// Poly-A tail length in nucleotides.
    pub poly_a_len: usize,
    /// `true` when a real replicase sequence was supplied (vs. the
    /// sized placeholder).
    pub replicase_supplied: bool,
}

impl SaRnaConstruct {
    /// Total construct length in nucleotides (replicase + promoter +
    /// gene + poly-A).
    pub fn len(&self) -> usize {
        self.replicase.len()
            + self.subgenomic_promoter.len()
            + self.gene_of_interest.len()
            + self.poly_a_len
    }

    /// `true` when the construct is empty (never produced).
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// A one-line layout description.
    pub fn layout(&self) -> String {
        format!(
            "saRNA: cap - replicase nsP1-4 ({} nt{}) - subgenomic promoter ({} nt) - \
             gene of interest ({} nt) - poly(A) {} nt; total {} nt",
            self.replicase.len(),
            if self.replicase_supplied {
                ""
            } else {
                ", placeholder"
            },
            self.subgenomic_promoter.len(),
            self.gene_of_interest.len(),
            self.poly_a_len,
            self.len(),
        )
    }
}

/// The conserved alphavirus subgenomic-promoter core (a representative
/// minimal element — the real promoter spans more context).
const SUBGENOMIC_PROMOTER: &[u8] = b"AUAGGCGGCGCAUGAGAGAAGCCCAGACCAAUUACCUACCCAAA";

/// The typical length of an alphavirus replicase (nsP1-4) ORF, used
/// when the caller supplies no replicase sequence.
const TYPICAL_REPLICASE_LEN: usize = 7500;

/// Designs a self-amplifying-mRNA construct layout (feature 23).
///
/// `gene_of_interest` is the CDS to express; `replicase` is the
/// alphavirus nsP1-4 region (pass an empty slice to use a sized
/// placeholder). The layout is `cap - replicase - subgenomic promoter
/// - gene - poly(A)`.
///
/// # Errors
/// [`GeneditingError::InvalidTarget`] if the gene of interest is not a
/// valid CDS, or a supplied replicase is non-ACGU.
pub fn design_sarna(
    gene_of_interest: &[u8],
    replicase: &[u8],
    poly_a_len: usize,
) -> Result<SaRnaConstruct> {
    let gene = validate_cds(gene_of_interest)?;
    let (replicase_rna, supplied) = if replicase.is_empty() {
        (vec![b'N'; TYPICAL_REPLICASE_LEN], false)
    } else {
        let r = transcribe(replicase);
        if !is_acgu(&r) {
            return Err(GeneditingError::invalid_target(
                "region",
                "replicase region must be A/C/G/U",
            ));
        }
        (r, true)
    };
    if poly_a_len > 500 {
        return Err(GeneditingError::invalid_target(
            "region",
            "poly-A tail length is beyond any physiological value",
        ));
    }
    Ok(SaRnaConstruct {
        replicase: replicase_rna,
        subgenomic_promoter: SUBGENOMIC_PROMOTER.to_vec(),
        gene_of_interest: gene,
        poly_a_len,
        replicase_supplied: supplied,
    })
}

/// A circular-RNA (circRNA) construct layout.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CircRnaConstruct {
    /// The 5′ homology / back-splicing arm (RNA).
    pub five_prime_arm: Vec<u8>,
    /// The IRES driving cap-independent translation (RNA). A
    /// placeholder of the typical length is used when none is supplied.
    pub ires: Vec<u8>,
    /// The gene-of-interest CDS (RNA).
    pub gene_of_interest: Vec<u8>,
    /// The 3′ homology / back-splicing arm (RNA).
    pub three_prime_arm: Vec<u8>,
    /// `true` when a real IRES sequence was supplied (vs. placeholder).
    pub ires_supplied: bool,
}

impl CircRnaConstruct {
    /// The length of the *circular* body once ligated (arms + IRES +
    /// gene). The homology arms participate in back-splicing; the
    /// mature circle is the whole thing closed end-to-end.
    pub fn circle_len(&self) -> usize {
        self.five_prime_arm.len()
            + self.ires.len()
            + self.gene_of_interest.len()
            + self.three_prime_arm.len()
    }

    /// `true` when the construct is empty.
    pub fn is_empty(&self) -> bool {
        self.circle_len() == 0
    }

    /// A one-line layout description.
    pub fn layout(&self) -> String {
        format!(
            "circRNA (no cap, no poly-A): 5' arm ({} nt) - IRES ({} nt{}) - \
             gene of interest ({} nt) - 3' arm ({} nt), back-spliced to a \
             {} nt covalently closed circle",
            self.five_prime_arm.len(),
            self.ires.len(),
            if self.ires_supplied {
                ""
            } else {
                ", placeholder"
            },
            self.gene_of_interest.len(),
            self.three_prime_arm.len(),
            self.circle_len(),
        )
    }
}

/// The typical length of a viral IRES (e.g. an EMCV / CVB3 IRES),
/// used when the caller supplies no IRES sequence.
const TYPICAL_IRES_LEN: usize = 600;

/// The minimum homology-arm length for a circRNA back-splice / PIE
/// ligation.
const MIN_CIRC_ARM_LEN: usize = 20;

/// Designs a circular-RNA construct layout (feature 23).
///
/// `gene_of_interest` is the CDS; `ires` is the internal-ribosome-entry
/// site driving translation (pass an empty slice to use a sized
/// placeholder); `arm_len` is the length of each back-splicing
/// homology arm. The mature product is a covalently closed circle with
/// no free ends.
///
/// # Errors
/// - [`GeneditingError::InvalidTarget`] if the gene of interest is not
///   a valid CDS or a supplied IRES is non-ACGU.
/// - [`GeneditingError::Invalid`] if `arm_len` is below the minimum
///   needed for a back-splice.
pub fn design_circrna(
    gene_of_interest: &[u8],
    ires: &[u8],
    arm_len: usize,
) -> Result<CircRnaConstruct> {
    let gene = validate_cds(gene_of_interest)?;
    if arm_len < MIN_CIRC_ARM_LEN {
        return Err(GeneditingError::invalid(
            "arm_len",
            format!("back-splicing homology arms must be at least {MIN_CIRC_ARM_LEN} nt"),
        ));
    }
    let (ires_rna, supplied) = if ires.is_empty() {
        (vec![b'N'; TYPICAL_IRES_LEN], false)
    } else {
        let r = transcribe(ires);
        if !is_acgu(&r) {
            return Err(GeneditingError::invalid_target(
                "region",
                "IRES must be A/C/G/U",
            ));
        }
        (r, true)
    };
    // The two homology arms are complementary repeats that bring the
    // ends together; here they are sized, neutral placeholder repeats.
    let five = vec![b'A'; arm_len];
    let three = vec![b'U'; arm_len];
    Ok(CircRnaConstruct {
        five_prime_arm: upper(&five),
        ires: ires_rna,
        gene_of_interest: gene,
        three_prime_arm: upper(&three),
        ires_supplied: supplied,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn good_cds() -> &'static [u8] {
        b"ATGGCCGCCGCCTAA"
    }

    #[test]
    fn sarna_layout_orders_parts() {
        let c = design_sarna(good_cds(), b"", 100).unwrap();
        assert!(!c.replicase_supplied); // placeholder
        assert_eq!(c.replicase.len(), TYPICAL_REPLICASE_LEN);
        assert_eq!(c.poly_a_len, 100);
        assert!(c.layout().contains("replicase"));
        assert!(c.layout().contains("subgenomic promoter"));
    }

    #[test]
    fn sarna_uses_supplied_replicase() {
        let c = design_sarna(good_cds(), b"ACGUACGUACGU", 100).unwrap();
        assert!(c.replicase_supplied);
        assert_eq!(c.replicase.len(), 12);
    }

    #[test]
    fn sarna_total_length_sums_parts() {
        let c = design_sarna(good_cds(), b"ACGUACGU", 50).unwrap();
        assert_eq!(
            c.len(),
            8 + SUBGENOMIC_PROMOTER.len() + good_cds().len() + 50
        );
    }

    #[test]
    fn sarna_rejects_bad_cds() {
        assert!(design_sarna(b"GCCGCCTAA", b"", 100).is_err()); // no start
        assert!(design_sarna(b"ATGGCCGCC", b"", 100).is_err()); // no stop
    }

    #[test]
    fn circrna_layout_has_no_cap_or_tail() {
        let c = design_circrna(good_cds(), b"", 30).unwrap();
        assert!(c.layout().contains("no cap"));
        assert!(c.layout().contains("no poly-A"));
        assert!(!c.ires_supplied);
        assert_eq!(c.ires.len(), TYPICAL_IRES_LEN);
    }

    #[test]
    fn circrna_circle_length_sums_parts() {
        let c = design_circrna(good_cds(), b"ACGUACGU", 25).unwrap();
        assert_eq!(c.circle_len(), 25 + 8 + good_cds().len() + 25);
    }

    #[test]
    fn circrna_uses_supplied_ires() {
        let c = design_circrna(good_cds(), b"ACGUACGUACGUACGU", 30).unwrap();
        assert!(c.ires_supplied);
        assert_eq!(c.ires.len(), 16);
    }

    #[test]
    fn circrna_rejects_short_arms() {
        let err = design_circrna(good_cds(), b"", 5).unwrap_err();
        assert_eq!(err.code(), "genediting.invalid");
    }

    #[test]
    fn circrna_rejects_bad_cds() {
        assert!(design_circrna(b"ATGGCC", b"", 30).is_err()); // not /3 + no stop
    }

    #[test]
    fn arms_have_the_requested_length() {
        let c = design_circrna(good_cds(), b"", 40).unwrap();
        assert_eq!(c.five_prime_arm.len(), 40);
        assert_eq!(c.three_prime_arm.len(), 40);
    }
}
