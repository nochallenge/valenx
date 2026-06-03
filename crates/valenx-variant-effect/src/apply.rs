//! Applying a [`Variant`] to a reference [`Seq`].
//!
//! [`apply`] validates that a variant is consistent with a reference
//! sequence (right kind, in range, and the named wild-type residue
//! actually present) and returns the mutant sequence with the
//! substitution made. [`translate_coding`] turns a coding-sequence `Seq`
//! into its protein consequence using the standard genetic code, so a
//! caller can compare wild-type and mutant proteins.
//!
//! All functions are pure: they return new [`Seq`] values and never write
//! to disk.

use valenx_bioseq::ops::translate::{translate_default, GeneticCode};
use valenx_bioseq::seq::{Seq, SeqKind};

use crate::error::VariantError;
use crate::variant::Variant;

/// Applies `variant` to `wild_type`, returning the mutant sequence.
///
/// Validation, by variant kind:
///
/// - [`Variant::ProteinSub`]: `wild_type` must be a protein sequence; the
///   1-based `pos` must be in range; and the residue at `pos` must equal
///   the variant's wild-type amino acid. The returned sequence is the
///   protein with the mutant residue substituted at `pos`.
/// - [`Variant::CodingSub`]: `wild_type` must be DNA; the 1-based `pos`
///   must be in range; and the base at `pos` must equal the variant's
///   wild-type base. The returned sequence is the DNA with the mutant
///   base substituted.
///
/// # Errors
///
/// - [`VariantError::Bioseq`] if `wild_type` is the wrong kind for the
///   variant (e.g. a protein variant against a DNA reference).
/// - [`VariantError::PositionOutOfRange`] if `pos` exceeds the length.
/// - [`VariantError::WildTypeMismatch`] if the reference residue at `pos`
///   differs from the variant's declared wild-type residue.
pub fn apply(variant: &Variant, wild_type: &Seq) -> Result<Seq, VariantError> {
    match *variant {
        Variant::ProteinSub { wt, pos, mt } => {
            if wild_type.kind() != SeqKind::Protein {
                return Err(VariantError::Bioseq(format!(
                    "protein variant requires a protein reference, got {}",
                    wild_type.kind().name()
                )));
            }
            substitute(wild_type, pos, wt.as_byte(), mt.as_byte(), SeqKind::Protein)
        }
        Variant::CodingSub { pos, wt, mt } => {
            if wild_type.kind() != SeqKind::Dna {
                return Err(VariantError::Bioseq(format!(
                    "coding variant requires a DNA reference, got {}",
                    wild_type.kind().name()
                )));
            }
            substitute(wild_type, pos, wt.as_byte(), mt.as_byte(), SeqKind::Dna)
        }
    }
}

/// Shared substitution core: bounds-check the 1-based `pos`, verify the
/// reference byte equals `expected`, and rebuild the sequence with `mt` at
/// `pos`. `kind` is the kind the mutant `Seq` is constructed as.
fn substitute(
    wild_type: &Seq,
    pos: usize,
    expected: u8,
    mt: u8,
    kind: SeqKind,
) -> Result<Seq, VariantError> {
    let len = wild_type.len();
    // 1-based: valid positions are 1..=len.
    if pos == 0 || pos > len {
        return Err(VariantError::PositionOutOfRange { position: pos, len });
    }
    let idx = pos - 1;
    let found = wild_type
        .get(idx)
        .expect("idx < len was just checked, so get() is Some");
    if found != expected {
        return Err(VariantError::WildTypeMismatch {
            position: pos,
            expected: expected as char,
            found: found as char,
        });
    }
    let mut residues = wild_type.as_bytes().to_vec();
    residues[idx] = mt;
    // The mutant residue is a single validated AA / base, so the rebuilt
    // sequence is guaranteed valid for `kind`; `Seq::new` re-validates.
    Seq::new(kind, residues).map_err(VariantError::from)
}

/// Translates a coding-sequence DNA [`Seq`] into its protein consequence,
/// using the standard (NCBI table 1) genetic code with full read-through
/// (stop codons appear as `*`).
///
/// # Errors
///
/// - [`VariantError::Bioseq`] if `dna` is not a DNA sequence.
/// - [`VariantError::Bioseq`] if translation itself fails.
pub fn translate_coding(dna: &Seq) -> Result<Seq, VariantError> {
    if dna.kind() != SeqKind::Dna {
        return Err(VariantError::Bioseq(format!(
            "translate_coding requires a DNA sequence, got {}",
            dna.kind().name()
        )));
    }
    translate_default(dna, &GeneticCode::standard()).map_err(VariantError::from)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::parse;

    fn dna(s: &str) -> Seq {
        Seq::new(SeqKind::Dna, s).unwrap()
    }
    fn protein(s: &str) -> Seq {
        Seq::new(SeqKind::Protein, s).unwrap()
    }

    // ---- protein substitution -------------------------------------------

    #[test]
    fn protein_substitution_succeeds() {
        // Reference protein ...R at position 3... -> H.
        let wt = protein("MKRTA");
        let v = parse("p.R3H").unwrap();
        let mutant = apply(&v, &wt).unwrap();
        assert_eq!(mutant.as_str(), "MKHTA");
        assert_eq!(mutant.kind(), SeqKind::Protein);
        // The original reference is untouched.
        assert_eq!(wt.as_str(), "MKRTA");
    }

    #[test]
    fn protein_substitution_at_first_and_last_position() {
        let wt = protein("MKRTA");
        assert_eq!(apply(&parse("p.M1L").unwrap(), &wt).unwrap().as_str(), "LKRTA");
        assert_eq!(apply(&parse("p.A5G").unwrap(), &wt).unwrap().as_str(), "MKRTG");
    }

    #[test]
    fn protein_wild_type_mismatch_is_detected() {
        let wt = protein("MKRTA");
        // Position 3 is R, not Q.
        let v = parse("p.Q3H").unwrap();
        match apply(&v, &wt) {
            Err(VariantError::WildTypeMismatch { position, expected, found }) => {
                assert_eq!(position, 3);
                assert_eq!(expected, 'Q');
                assert_eq!(found, 'R');
            }
            other => panic!("expected WildTypeMismatch, got {other:?}"),
        }
    }

    #[test]
    fn protein_position_out_of_range_is_detected() {
        let wt = protein("MKRTA"); // length 5
        let v = parse("p.A6G").unwrap();
        match apply(&v, &wt) {
            Err(VariantError::PositionOutOfRange { position, len }) => {
                assert_eq!(position, 6);
                assert_eq!(len, 5);
            }
            other => panic!("expected PositionOutOfRange, got {other:?}"),
        }
    }

    #[test]
    fn protein_variant_against_dna_reference_errors() {
        let wt = dna("ATGCGT");
        let v = parse("p.R1H").unwrap();
        assert!(matches!(apply(&v, &wt), Err(VariantError::Bioseq(_))));
    }

    // ---- coding substitution --------------------------------------------

    #[test]
    fn coding_substitution_succeeds() {
        // ATGCGT, positions A1 T2 G3 C4 G5 T6. c.5G>A -> ATGCAT.
        let wt = dna("ATGCGT");
        let v = parse("c.5G>A").unwrap();
        let mutant = apply(&v, &wt).unwrap();
        assert_eq!(mutant.as_str(), "ATGCAT");
        assert_eq!(mutant.kind(), SeqKind::Dna);
    }

    #[test]
    fn coding_wild_type_mismatch_is_detected() {
        let wt = dna("ATGCGT");
        // Position 4 is C, not G.
        let v = parse("c.4G>A").unwrap();
        match apply(&v, &wt) {
            Err(VariantError::WildTypeMismatch { position, expected, found }) => {
                assert_eq!(position, 4);
                assert_eq!(expected, 'G');
                assert_eq!(found, 'C');
            }
            other => panic!("expected WildTypeMismatch, got {other:?}"),
        }
    }

    #[test]
    fn coding_position_out_of_range_is_detected() {
        let wt = dna("ATGCGT"); // length 6
        let v = parse("c.7C>T").unwrap();
        assert!(matches!(
            apply(&v, &wt),
            Err(VariantError::PositionOutOfRange { position: 7, len: 6 })
        ));
    }

    #[test]
    fn coding_variant_against_protein_reference_errors() {
        let wt = protein("MKRTA");
        let v = parse("c.1A>T").unwrap();
        assert!(matches!(apply(&v, &wt), Err(VariantError::Bioseq(_))));
    }

    // ---- translate_coding: missense vs synonymous -----------------------

    #[test]
    fn translate_coding_reference() {
        // ATG CGT = Met Arg.
        let wt = dna("ATGCGT");
        assert_eq!(translate_coding(&wt).unwrap().as_str(), "MR");
    }

    #[test]
    fn coding_missense_changes_the_protein() {
        // c.5G>A: CGT (Arg) -> CAT (His). Protein MR -> MH.
        let wt = dna("ATGCGT");
        let mutant = apply(&parse("c.5G>A").unwrap(), &wt).unwrap();
        let wt_protein = translate_coding(&wt).unwrap();
        let mt_protein = translate_coding(&mutant).unwrap();
        assert_eq!(wt_protein.as_str(), "MR");
        assert_eq!(mt_protein.as_str(), "MH");
        assert_ne!(wt_protein, mt_protein, "missense must change the protein");
    }

    #[test]
    fn coding_synonymous_leaves_the_protein_unchanged() {
        // c.6T>A: CGT (Arg) -> CGA (Arg), both encode Arg. Protein stays MR.
        let wt = dna("ATGCGT");
        let mutant = apply(&parse("c.6T>A").unwrap(), &wt).unwrap();
        assert_eq!(mutant.as_str(), "ATGCGA");
        let wt_protein = translate_coding(&wt).unwrap();
        let mt_protein = translate_coding(&mutant).unwrap();
        assert_eq!(mt_protein.as_str(), "MR");
        assert_eq!(
            wt_protein, mt_protein,
            "synonymous substitution must leave the protein unchanged"
        );
    }

    #[test]
    fn translate_coding_rejects_non_dna() {
        let p = protein("MKRTA");
        assert!(matches!(translate_coding(&p), Err(VariantError::Bioseq(_))));
    }

    #[test]
    fn coding_sub_in_dropped_partial_codon_leaves_protein_unchanged() {
        // Contract: `translate_coding` treats the reference as an in-frame
        // CDS and drops a partial trailing codon (len % 3 leftover bases).
        // A substitution that lands entirely within those dropped base(s)
        // therefore cannot change the translated protein.
        //
        // ATGCGTA is 7 bases -> 2 whole codons ATG|CGT (= Met Arg) plus a
        // leftover trailing `A` at position 7 that translation drops.
        let wt = dna("ATGCGTA");
        assert_eq!(wt.len() % 3, 1, "reference must have a partial trailing codon");

        // c.7A>C mutates only that dropped 7th base.
        let mutant = apply(&parse("c.7A>C").unwrap(), &wt).unwrap();
        assert_eq!(mutant.as_str(), "ATGCGTC");

        let wt_protein = translate_coding(&wt).unwrap();
        let mt_protein = translate_coding(&mutant).unwrap();
        assert_eq!(wt_protein.as_str(), "MR");
        assert_eq!(
            wt_protein, mt_protein,
            "a substitution in the dropped partial codon must not change the protein"
        );
    }
}
