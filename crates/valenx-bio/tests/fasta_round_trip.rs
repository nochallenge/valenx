use valenx_bio::format::fasta;
use valenx_bio::{Alphabet, Sequence};

const SAMPLE_FASTA: &str = "\
>p53|tp53_HUMAN cellular tumor antigen
MEEPQSDPSV
EPPLSQETFS
>ubq|ubiquitin
MQIFVKTLTG
KTITLEVEPS
DTIENVKAKI
QDKEGIPPDQ
QRLIFAGKQL
EDGRTLSDYN
IQKESTLHLV
LRLRGG
";

#[test]
fn read_round_trip() {
    let seqs: Vec<Sequence> = fasta::read(SAMPLE_FASTA, Alphabet::Protein).unwrap();
    assert_eq!(seqs.len(), 2);
    assert_eq!(seqs[0].name, "p53|tp53_HUMAN cellular tumor antigen");
    // Multi-line body concatenated into one sequence.
    assert_eq!(seqs[0].as_str(), "MEEPQSDPSVEPPLSQETFS");
    assert_eq!(seqs[1].name, "ubq|ubiquitin");
    assert_eq!(seqs[1].len(), 76);
}

#[test]
fn write_round_trip() {
    let s1 = Sequence::new("a", Alphabet::Dna, "ACGT").unwrap();
    let s2 = Sequence::new("b", Alphabet::Dna, "GGGGAAAA").unwrap();
    let body = fasta::write(&[s1.clone(), s2.clone()]);
    let parsed = fasta::read(&body, Alphabet::Dna).unwrap();
    assert_eq!(parsed.len(), 2);
    assert_eq!(parsed[0], s1);
    assert_eq!(parsed[1], s2);
}

#[test]
fn empty_input_returns_empty_vec() {
    let seqs = fasta::read("", Alphabet::Dna).unwrap();
    assert!(seqs.is_empty());
}

#[test]
fn invalid_alphabet_byte_surfaces_position() {
    // Capital E (glutamate) is invalid in DNA → error.
    // Don't use `M` — it's a valid IUPAC ambiguity code for
    // `aMino` = A or C, so DNA accepts it.
    let err = fasta::read(">x\nACGTE\n", Alphabet::Dna).unwrap_err();
    assert!(format!("{err}").contains("position 4"));
}
