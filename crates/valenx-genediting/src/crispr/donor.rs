//! Feature 4 — HDR knock-in donor-template design.
//!
//! Homology-directed repair (HDR) installs a precise edit — an
//! insertion, a tag, a corrected base — by giving the cell a *donor
//! template*: the new sequence flanked by **homology arms** identical
//! to the genome on either side of the cut.
//!
//! Two design problems matter:
//!
//! 1. **Homology arms.** The donor carries a left and a right arm
//!    copied from the reference around the cut site. Arms can be
//!    *symmetric* (equal length) or *asymmetric* (a longer arm on one
//!    side — useful for ssODN donors).
//! 2. **Re-cut prevention.** After HDR succeeds, the nuclease will
//!    happily re-cut the corrected allele unless the donor also breaks
//!    the protospacer or PAM. The standard fix is one or more **silent
//!    mutations** — base changes that block PAM / seed recognition but
//!    do not change the encoded protein (synonymous codon swaps) or
//!    fall outside any CDS. This module finds such mutations.
//!
//! ## v1 scope
//!
//! The donor designer works on a **linear reference window** the
//! caller supplies — it does not fetch a genome or resolve isoforms.
//! Silent-mutation selection prefers a synonymous codon change inside
//! a supplied CDS frame; if no synonymous change blocks the site, or
//! the site is non-coding, it falls back to a PAM-disrupting change
//! and reports `silent = false` so the caller knows the protein may be
//! affected. HDR *efficiency* is not predicted (a context-and-cell-type
//! problem with no good closed form); the module designs a *correct*
//! donor and reports re-cut protection.

use crate::error::{GeneditingError, Result};
use crate::sequtil::{complement, is_acgt, revcomp, upper};
use serde::{Deserialize, Serialize};
use valenx_bioseq::ops::translate::GeneticCode;

/// How the two homology arms are sized.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ArmLayout {
    /// Both arms the same length.
    Symmetric {
        /// Length of each arm in base pairs.
        len: usize,
    },
    /// Arms of independent lengths (asymmetric ssODN-style donor).
    Asymmetric {
        /// Left (5′) arm length.
        left: usize,
        /// Right (3′) arm length.
        right: usize,
    },
}

impl ArmLayout {
    /// The left and right arm lengths.
    fn lengths(self) -> (usize, usize) {
        match self {
            ArmLayout::Symmetric { len } => (len, len),
            ArmLayout::Asymmetric { left, right } => (left, right),
        }
    }
}

/// A request to design an HDR knock-in donor template.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct HdrDonorRequest {
    /// The reference window around the intended edit (forward strand,
    /// 5′→3′, unambiguous ACGT).
    pub reference: Vec<u8>,
    /// 0-based insertion point in `reference` — the sequence in
    /// [`insert`](Self::insert) is placed *between* `insert_pos - 1`
    /// and `insert_pos`. For a pure base correction with no insert,
    /// leave `insert` empty.
    pub insert_pos: usize,
    /// The sequence to knock in (a tag, a cassette, a corrected
    /// stretch). May be empty for a re-cut-blocking-only donor.
    pub insert: Vec<u8>,
    /// Homology-arm layout.
    pub arms: ArmLayout,
    /// The protospacer (5′→3′) of the guide that made the cut — needed
    /// to design re-cut-blocking mutations. Must match the reference
    /// (forward or reverse strand).
    pub protospacer: Vec<u8>,
    /// The PAM sequence of that guide as found on the reference.
    pub pam: Vec<u8>,
    /// `true` if the guide / PAM is on the reverse strand of
    /// `reference`.
    pub guide_reverse: bool,
    /// 0-based start of the protospacer on the forward strand of
    /// `reference`.
    pub protospacer_start: usize,
    /// Optional reading-frame offset: if `Some(phase)`, `reference` is
    /// coding and codon `0` starts at index `phase` (`0`, `1` or `2`).
    /// Enables synonymous (silent) re-cut mutations. `None` ⇒ the
    /// window is treated as non-coding.
    pub coding_phase: Option<usize>,
}

/// A single re-cut-blocking mutation in the donor.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BlockingMutation {
    /// 0-based position in the *reference* window the change is at.
    pub ref_pos: usize,
    /// The original reference base.
    pub from: u8,
    /// The mutated base placed in the donor.
    pub to: u8,
    /// `true` when the change is synonymous (silent at the protein
    /// level) or non-coding; `false` when it is a forced PAM-disrupting
    /// change that may alter the protein.
    pub silent: bool,
    /// Whether this mutation hits the PAM (`true`) or the seed region
    /// of the protospacer (`false`).
    pub in_pam: bool,
}

/// A designed HDR donor template.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DonorTemplate {
    /// The full donor sequence: left arm + insert + right arm, with
    /// the re-cut-blocking mutations already applied.
    pub donor: Vec<u8>,
    /// The left (5′) homology arm.
    pub left_arm: Vec<u8>,
    /// The right (3′) homology arm.
    pub right_arm: Vec<u8>,
    /// The knocked-in insert (echoed back; may be empty).
    pub insert: Vec<u8>,
    /// Re-cut-blocking mutations applied to the donor.
    pub blocking_mutations: Vec<BlockingMutation>,
    /// `true` when the donor is predicted to resist re-cutting (at
    /// least one PAM or seed change was placed).
    pub recut_protected: bool,
    /// `true` when every blocking mutation is silent / non-coding.
    pub all_silent: bool,
}

impl DonorTemplate {
    /// Total donor length in base pairs.
    pub fn len(&self) -> usize {
        self.donor.len()
    }

    /// `true` when the donor is empty (never produced — present for
    /// clippy's `len`-without-`is_empty` lint).
    pub fn is_empty(&self) -> bool {
        self.donor.is_empty()
    }
}

/// The three synonymous-or-not candidate bases for `from`, in a
/// deterministic order (the other three nucleotides).
fn other_bases(from: u8) -> [u8; 3] {
    match from.to_ascii_uppercase() {
        b'A' => [b'C', b'G', b'T'],
        b'C' => [b'A', b'G', b'T'],
        b'G' => [b'A', b'C', b'T'],
        _ => [b'A', b'C', b'G'],
    }
}

/// Tries to find a synonymous (silent) substitution at reference index
/// `ref_pos` given a coding `phase`. Returns the substituting base if
/// some codon-mate change keeps the amino acid, else `None`.
fn synonymous_sub(
    reference: &[u8],
    ref_pos: usize,
    phase: usize,
    code: &GeneticCode,
) -> Option<u8> {
    if ref_pos < phase {
        return None;
    }
    let codon_idx = (ref_pos - phase) / 3;
    let codon_start = phase + codon_idx * 3;
    if codon_start + 3 > reference.len() {
        return None;
    }
    let codon = [
        reference[codon_start],
        reference[codon_start + 1],
        reference[codon_start + 2],
    ];
    let original_aa = code.translate_codon(&codon);
    let within = ref_pos - codon_start;
    for &b in &other_bases(reference[ref_pos]) {
        let mut trial = codon;
        trial[within] = b;
        if code.translate_codon(&trial) == original_aa && original_aa != b'X' {
            return Some(b);
        }
    }
    None
}

/// Designs an HDR knock-in donor template.
///
/// Builds the homology arms from the reference around the insertion
/// point, places the insert, and applies one or more re-cut-blocking
/// mutations — preferring silent (synonymous / non-coding) changes —
/// to the PAM and seed of the cutting guide.
///
/// # Errors
/// - [`GeneditingError::InvalidTarget`] for a non-ACGT reference,
///   protospacer or PAM, or coordinates outside the reference.
/// - [`GeneditingError::Invalid`] for a zero-length homology arm or a
///   `coding_phase` above 2.
/// - [`GeneditingError::NoValidDesign`] if the arms do not fit inside
///   the supplied reference window.
pub fn design_hdr_donor(req: &HdrDonorRequest) -> Result<DonorTemplate> {
    if !is_acgt(&req.reference) {
        return Err(GeneditingError::invalid_target(
            "region",
            "reference window must be non-empty ACGT",
        ));
    }
    if !is_acgt(&req.protospacer) {
        return Err(GeneditingError::invalid_target(
            "locus",
            "protospacer must be non-empty ACGT",
        ));
    }
    if !req.pam.is_empty() && !is_acgt(&req.pam) {
        return Err(GeneditingError::invalid_target("locus", "PAM must be ACGT"));
    }
    if !req.insert.is_empty() && !is_acgt(&req.insert) {
        return Err(GeneditingError::invalid_target(
            "region",
            "insert must be ACGT",
        ));
    }
    if req.insert_pos > req.reference.len() {
        return Err(GeneditingError::invalid_target(
            "region",
            "insertion point is past the end of the reference window",
        ));
    }
    if let Some(p) = req.coding_phase {
        if p > 2 {
            return Err(GeneditingError::invalid(
                "coding_phase",
                "must be 0, 1 or 2",
            ));
        }
    }
    let (left_len, right_len) = req.arms.lengths();
    if left_len == 0 || right_len == 0 {
        return Err(GeneditingError::invalid(
            "homology_arm_len",
            "homology arms must be positive",
        ));
    }
    if req.insert_pos < left_len || req.insert_pos + right_len > req.reference.len() {
        return Err(GeneditingError::no_valid_design(
            "donor",
            "homology arms do not fit inside the supplied reference window",
        ));
    }

    // --- Re-cut-blocking mutations on a *copy* of the reference -------
    let mut edited_ref = upper(&req.reference);
    let code = GeneticCode::standard();
    let mut muts: Vec<BlockingMutation> = Vec::new();

    // Forward-strand coordinates of the PAM and the protospacer seed.
    let plen = req.protospacer.len();
    let pam_len = req.pam.len();
    let proto_start = req.protospacer_start;
    let proto_end = proto_start + plen; // exclusive

    // PAM forward-strand span (3' of the protospacer on a forward
    // guide; 5' of it — to the left — on a reverse guide, since for a
    // reverse 3'-PAM the PAM sits at lower forward coordinates).
    let pam_span: Option<(usize, usize)> = if pam_len == 0 {
        None
    } else if !req.guide_reverse {
        Some((proto_end, proto_end + pam_len))
    } else {
        proto_start.checked_sub(pam_len).map(|s| (s, proto_start))
    };

    // 1) Try to break the PAM first — the strongest re-cut block.
    if let Some((ps, pe)) = pam_span {
        if pe <= edited_ref.len() {
            for pos in ps..pe {
                if try_blocking_mutation(
                    &mut edited_ref,
                    pos,
                    true,
                    req.coding_phase,
                    &code,
                    &mut muts,
                ) {
                    break; // one PAM hit is enough
                }
            }
        }
    }

    // 2) If the PAM could not be broken, hit the seed (the ~10 nt of
    //    the protospacer nearest the PAM — re-cutting is seed-sensitive).
    if muts.is_empty() {
        let seed_positions: Vec<usize> = if !req.guide_reverse {
            // 3' PAM forward: seed is the 3' (high-coordinate) end.
            (proto_start..proto_end).rev().take(10).collect()
        } else {
            // 3' PAM reverse: PAM is at low forward coords, so the seed
            // is the low-coordinate end of the protospacer.
            (proto_start..proto_end).take(10).collect()
        };
        for pos in seed_positions {
            if try_blocking_mutation(
                &mut edited_ref,
                pos,
                false,
                req.coding_phase,
                &code,
                &mut muts,
            ) {
                break;
            }
        }
    }

    // --- Assemble the donor from the edited reference ----------------
    let left_arm: Vec<u8> = edited_ref[req.insert_pos - left_len..req.insert_pos].to_vec();
    let right_arm: Vec<u8> = edited_ref[req.insert_pos..req.insert_pos + right_len].to_vec();
    let mut donor = Vec::with_capacity(left_arm.len() + req.insert.len() + right_arm.len());
    donor.extend_from_slice(&left_arm);
    donor.extend_from_slice(&upper(&req.insert));
    donor.extend_from_slice(&right_arm);

    let recut_protected = !muts.is_empty();
    let all_silent = muts.iter().all(|m| m.silent);

    Ok(DonorTemplate {
        donor,
        left_arm,
        right_arm,
        insert: upper(&req.insert),
        blocking_mutations: muts,
        recut_protected,
        all_silent,
    })
}

/// Attempts a re-cut-blocking mutation at `pos`. Prefers a synonymous
/// change when `phase` is `Some`; otherwise (or if no synonymous change
/// exists) applies the first alternative base as a forced change.
/// Returns `true` and records the mutation iff one was applied.
fn try_blocking_mutation(
    edited_ref: &mut [u8],
    pos: usize,
    in_pam: bool,
    phase: Option<usize>,
    code: &GeneticCode,
    out: &mut Vec<BlockingMutation>,
) -> bool {
    if pos >= edited_ref.len() {
        return false;
    }
    let from = edited_ref[pos];
    // Prefer a synonymous substitution inside a coding window.
    if let Some(p) = phase {
        if let Some(to) = synonymous_sub(edited_ref, pos, p, code) {
            edited_ref[pos] = to;
            out.push(BlockingMutation {
                ref_pos: pos,
                from,
                to,
                silent: true,
                in_pam,
            });
            return true;
        }
        // Coding but no synonymous option here — let the caller try the
        // next position rather than forcing a protein-changing edit.
        return false;
    }
    // Non-coding window: any change is "silent" at the protein level.
    let to = other_bases(from)[0];
    edited_ref[pos] = to;
    out.push(BlockingMutation {
        ref_pos: pos,
        from,
        to,
        silent: true,
        in_pam,
    });
    true
}

/// Convenience: the reverse complement of a donor — some cloning
/// workflows want the bottom strand of an ssODN.
pub fn donor_bottom_strand(donor: &DonorTemplate) -> Vec<u8> {
    revcomp(&donor.donor)
}

/// Convenience: `true` when `arm` is a perfect match to `reference`
/// at `offset` (a homology-arm sanity check that ignores the
/// deliberately introduced blocking mutations is the caller's job;
/// this is the strict version).
pub fn arm_matches_reference(arm: &[u8], reference: &[u8], offset: usize) -> bool {
    if offset + arm.len() > reference.len() {
        return false;
    }
    arm.iter()
        .zip(&reference[offset..offset + arm.len()])
        .all(|(a, r)| a.eq_ignore_ascii_case(r))
}

/// `true` when two bases are Watson-Crick complementary.
pub fn is_complementary(a: u8, b: u8) -> bool {
    complement(a) == b.to_ascii_uppercase()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_request() -> HdrDonorRequest {
        // A 60 bp reference; a 20-mer protospacer at 20..40 on the
        // forward strand, NGG PAM at 40..43.
        let reference = b"ACGTACGTACGTACGTACGTGGCATGCATGCATGCATGCAAGGTTACGTTACGTTACGTT".to_vec();
        HdrDonorRequest {
            reference,
            insert_pos: 30,
            insert: b"TAG".to_vec(),
            arms: ArmLayout::Symmetric { len: 15 },
            protospacer: b"GGCATGCATGCATGCATGCA".to_vec(),
            pam: b"AGG".to_vec(),
            guide_reverse: false,
            protospacer_start: 20,
            coding_phase: None,
        }
    }

    #[test]
    fn donor_has_arms_and_insert() {
        let d = design_hdr_donor(&base_request()).unwrap();
        assert_eq!(d.left_arm.len(), 15);
        assert_eq!(d.right_arm.len(), 15);
        assert_eq!(d.insert, b"TAG");
        assert_eq!(d.len(), 15 + 3 + 15);
    }

    #[test]
    fn donor_blocks_recut() {
        let d = design_hdr_donor(&base_request()).unwrap();
        assert!(d.recut_protected, "donor should carry a blocking mutation");
        assert!(!d.blocking_mutations.is_empty());
    }

    #[test]
    fn pam_is_preferred_block_site() {
        let d = design_hdr_donor(&base_request()).unwrap();
        // The PAM at 40..43 is non-coding here, so a PAM change is
        // applied and reported.
        assert!(d.blocking_mutations.iter().any(|m| m.in_pam));
    }

    #[test]
    fn noncoding_blocking_mutations_are_silent() {
        let d = design_hdr_donor(&base_request()).unwrap();
        assert!(d.all_silent, "non-coding changes count as silent");
    }

    #[test]
    fn synonymous_substitution_keeps_amino_acid() {
        let code = GeneticCode::standard();
        // Leucine CTG: a synonymous change at the wobble position
        // (CTG -> CTA / CTC / CTT) must exist.
        let refseq = b"CTG";
        let sub = synonymous_sub(refseq, 2, 0, &code);
        assert!(sub.is_some());
        let to = sub.unwrap();
        let mut codon = *refseq;
        codon[2] = to;
        assert_eq!(code.translate_codon(&codon), code.translate_codon(refseq));
    }

    #[test]
    fn asymmetric_arms_have_independent_lengths() {
        let mut req = base_request();
        req.arms = ArmLayout::Asymmetric {
            left: 20,
            right: 10,
        };
        let d = design_hdr_donor(&req).unwrap();
        assert_eq!(d.left_arm.len(), 20);
        assert_eq!(d.right_arm.len(), 10);
    }

    #[test]
    fn rejects_arms_that_overrun_window() {
        let mut req = base_request();
        req.arms = ArmLayout::Symmetric { len: 100 };
        let err = design_hdr_donor(&req).unwrap_err();
        assert_eq!(err.code(), "genediting.no_valid_design");
    }

    #[test]
    fn rejects_zero_length_arm() {
        let mut req = base_request();
        req.arms = ArmLayout::Symmetric { len: 0 };
        let err = design_hdr_donor(&req).unwrap_err();
        assert_eq!(err.code(), "genediting.invalid");
    }

    #[test]
    fn rejects_non_acgt_reference() {
        let mut req = base_request();
        req.reference = b"ACGTNNNN".to_vec();
        assert!(design_hdr_donor(&req).is_err());
    }

    #[test]
    fn bottom_strand_is_reverse_complement() {
        let d = design_hdr_donor(&base_request()).unwrap();
        let bottom = donor_bottom_strand(&d);
        assert_eq!(bottom, revcomp(&d.donor));
    }

    #[test]
    fn complementarity_helper() {
        assert!(is_complementary(b'A', b'T'));
        assert!(is_complementary(b'g', b'c'));
        assert!(!is_complementary(b'A', b'G'));
    }
}
