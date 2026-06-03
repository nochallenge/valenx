use valenx_bio::{Alignment, AlignmentError, Alphabet, Sequence};

#[test]
fn build_alignment_from_aligned_sequences() {
    let s1 = Sequence::new("a", Alphabet::Protein, "ACDE").unwrap();
    let s2 = Sequence::new("b", Alphabet::Protein, "AC-E").unwrap();
    let s3 = Sequence::new("c", Alphabet::Protein, "-CDE").unwrap();
    let aln = Alignment::new(vec![s1, s2, s3]).expect("aligned");
    assert_eq!(aln.row_count(), 3);
    assert_eq!(aln.column_count(), 4);
    assert_eq!(aln.alphabet(), Alphabet::Protein);
}

#[test]
fn rejects_unequal_row_lengths() {
    let s1 = Sequence::new("a", Alphabet::Dna, "ACGT").unwrap();
    let s2 = Sequence::new("b", Alphabet::Dna, "AC").unwrap();
    let err = Alignment::new(vec![s1, s2]).unwrap_err();
    assert!(matches!(err, AlignmentError::UnequalRows { .. }));
}

#[test]
fn rejects_mixed_alphabets() {
    let s1 = Sequence::new("a", Alphabet::Dna, "ACGT").unwrap();
    let s2 = Sequence::new("b", Alphabet::Rna, "ACGU").unwrap();
    let err = Alignment::new(vec![s1, s2]).unwrap_err();
    assert!(matches!(err, AlignmentError::AlphabetMismatch { .. }));
}

#[test]
fn rejects_empty_alignment() {
    let err = Alignment::new(vec![]).unwrap_err();
    assert!(matches!(err, AlignmentError::Empty));
}

#[test]
fn gap_count_helper() {
    let s1 = Sequence::new("a", Alphabet::Protein, "AC-E").unwrap();
    let s2 = Sequence::new("b", Alphabet::Protein, "--DE").unwrap();
    let aln = Alignment::new(vec![s1, s2]).unwrap();
    assert_eq!(aln.gap_count(), 3);
}
