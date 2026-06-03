//! `valenx-bioseq` interop — align [`Seq`] values directly.
//!
//! The DP routines in [`crate::pairwise`] / [`crate::msa`] operate on
//! raw `&[u8]` for speed. This module is the thin, ergonomic layer on
//! top: it accepts validated [`Seq`] values from Block 6.1, checks
//! that the sequences are alignment-compatible (same [`SeqKind`]),
//! picks a sensible default [`crate::matrix::ScoringScheme`] for the
//! kind, and forwards to the byte-slice routines.
//!
//! Use these when you already hold `Seq` values; drop to the
//! `&[u8]` routines when you want a custom scoring scheme or are in a
//! hot loop.

use crate::error::{AlignError, Result};
use crate::matrix::ScoringScheme;
use crate::msa::progressive::{align as align_bytes, Msa};
use crate::pairwise::global::gotoh;
use crate::pairwise::local::smith_waterman;
use crate::pairwise::result::Alignment;
use valenx_bioseq::{Seq, SeqKind};

/// The default [`ScoringScheme`] for a sequence kind: `NUC.4.4` with
/// `10/1` gaps for nucleotides, BLOSUM62 with `11/1` gaps for protein.
pub fn default_scheme(kind: SeqKind) -> ScoringScheme {
    if kind.is_nucleotide() {
        ScoringScheme::dna_default()
    } else {
        ScoringScheme::blosum62_default()
    }
}

/// Checks that every sequence shares one [`SeqKind`], returning it.
/// Returns [`AlignError::Invalid`] on an empty list or a kind mismatch.
fn unify_kind(seqs: &[&Seq]) -> Result<SeqKind> {
    let first = seqs
        .first()
        .ok_or_else(|| AlignError::invalid("seqs", "need >= 1 sequence"))?
        .kind();
    for s in seqs {
        if s.kind() != first {
            return Err(AlignError::invalid(
                "kind",
                format!(
                    "cannot align {} with {}",
                    first.name(),
                    s.kind().name()
                ),
            ));
        }
    }
    Ok(first)
}

/// Global (Needleman-Wunsch) alignment of two [`Seq`] values using the
/// kind's [`default_scheme`]. Both sequences must be the same kind.
///
/// This uses the *affine*-gap global DP ([`gotoh`]): the default
/// schemes carry an affine `open/extend` gap cost (`10/1` for DNA,
/// `11/1` for protein), and the plain linear-gap
/// [`needleman_wunsch`](crate::pairwise::global::needleman_wunsch)
/// ignores the `open` term — pairing the two would make a gap-pair
/// cheaper than a single mismatch.
pub fn global_align(a: &Seq, b: &Seq) -> Result<Alignment> {
    let kind = unify_kind(&[a, b])?;
    gotoh(a.as_bytes(), b.as_bytes(), &default_scheme(kind))
}

/// Global alignment of two [`Seq`] values with a caller-supplied
/// scoring scheme. Uses the affine-gap DP so the scheme's `open` term
/// is honoured.
pub fn global_align_with(a: &Seq, b: &Seq, scheme: &ScoringScheme) -> Result<Alignment> {
    unify_kind(&[a, b])?;
    gotoh(a.as_bytes(), b.as_bytes(), scheme)
}

/// Local (Smith-Waterman) alignment of two [`Seq`] values using the
/// kind's [`default_scheme`].
pub fn local_align(a: &Seq, b: &Seq) -> Result<Alignment> {
    let kind = unify_kind(&[a, b])?;
    smith_waterman(a.as_bytes(), b.as_bytes(), &default_scheme(kind))
}

/// Percent identity of two [`Seq`] values from their global alignment,
/// in `[0, 1]`.
pub fn identity(a: &Seq, b: &Seq) -> Result<f64> {
    Ok(global_align(a, b)?.percent_identity())
}

/// Progressive multiple-sequence alignment of a set of [`Seq`] values
/// using the kind's [`default_scheme`]. All sequences must share a
/// kind.
pub fn multiple_align(seqs: &[&Seq]) -> Result<Msa> {
    let kind = unify_kind(seqs)?;
    let byte_seqs: Vec<&[u8]> = seqs.iter().map(|s| s.as_bytes()).collect();
    align_bytes(&byte_seqs, &default_scheme(kind))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_scheme_by_kind() {
        // Nucleotide default uses NUC.4.4 (A/A == 5).
        assert_eq!(default_scheme(SeqKind::Dna).sub(b'A', b'A'), 5);
        // Protein default uses BLOSUM62 (A/A == 4).
        assert_eq!(default_scheme(SeqKind::Protein).sub(b'A', b'A'), 4);
    }

    #[test]
    fn global_align_dna_seqs() {
        let a = Seq::new(SeqKind::Dna, "ACGTACGT").unwrap();
        let b = Seq::new(SeqKind::Dna, "ACGTACGT").unwrap();
        let al = global_align(&a, &b).unwrap();
        assert_eq!(al.row1_str(), "ACGTACGT");
        assert!((al.percent_identity() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn kind_mismatch_rejected() {
        let dna = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        let protein = Seq::new(SeqKind::Protein, "MKVL").unwrap();
        assert!(global_align(&dna, &protein).is_err());
        assert!(multiple_align(&[&dna, &protein]).is_err());
    }

    #[test]
    fn local_align_finds_core() {
        let a = Seq::new(SeqKind::Dna, "TTTTGATTACATTTT").unwrap();
        let b = Seq::new(SeqKind::Dna, "CCCGATTACACCC").unwrap();
        let al = local_align(&a, &b).unwrap();
        assert_eq!(al.row1_str(), "GATTACA");
    }

    #[test]
    fn identity_of_seqs() {
        let a = Seq::new(SeqKind::Dna, "ACGTACGT").unwrap();
        let b = Seq::new(SeqKind::Dna, "ACGTACGA").unwrap();
        let id = identity(&a, &b).unwrap();
        // 7 of 8 identical.
        assert!((id - 0.875).abs() < 1e-9);
    }

    #[test]
    fn protein_global_align_with_custom_scheme() {
        let a = Seq::new(SeqKind::Protein, "MKVLAAGG").unwrap();
        let b = Seq::new(SeqKind::Protein, "MKVLAAGG").unwrap();
        let scheme = ScoringScheme::blosum62_default();
        let al = global_align_with(&a, &b, &scheme).unwrap();
        assert_eq!(al.row1_str(), "MKVLAAGG");
    }

    #[test]
    fn multiple_align_dna_seqs() {
        let s1 = Seq::new(SeqKind::Dna, "ACGTACGT").unwrap();
        let s2 = Seq::new(SeqKind::Dna, "ACGTACGT").unwrap();
        let s3 = Seq::new(SeqKind::Dna, "ACGTCGT").unwrap();
        let msa = multiple_align(&[&s1, &s2, &s3]).unwrap();
        assert_eq!(msa.depth(), 3);
        // All rows equal length.
        let w = msa.width();
        assert!(msa.rows.iter().all(|r| r.len() == w));
    }

    #[test]
    fn empty_input_rejected() {
        let empty: &[&Seq] = &[];
        assert!(multiple_align(empty).is_err());
    }
}
