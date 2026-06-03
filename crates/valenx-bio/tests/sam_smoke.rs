use valenx_bio::format::sam;

const SAMPLE: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
@PG\tID:bwa\tPN:bwa\tVN:0.7.17
read1\t0\tchr1\t100\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII
read2\t4\t*\t0\t0\t*\t*\t0\t0\tNNNNNNNN\t!!!!!!!!
";

#[test]
fn parses_sam_with_header_and_records() {
    let sam = sam::read_str(SAMPLE).expect("parse");
    assert_eq!(sam.header.len(), 3);
    assert_eq!(sam.records.len(), 2);
    let r0 = &sam.records[0];
    assert_eq!(r0.qname, "read1");
    assert_eq!(r0.flag, 0);
    assert_eq!(r0.rname, "chr1");
    assert_eq!(r0.pos, 100);
    assert_eq!(r0.mapq, 60);
    assert_eq!(r0.cigar, "8M");
    assert_eq!(r0.seq, "ACGTACGT");
    assert!(r0.is_mapped());
    let r1 = &sam.records[1];
    assert!(!r1.is_mapped()); // flag 4 = unmapped
}

#[test]
fn empty_sam_parses_to_empty_struct() {
    let sam = sam::read_str("").expect("parse");
    assert!(sam.header.is_empty());
    assert!(sam.records.is_empty());
}

#[test]
fn rejects_record_with_too_few_fields() {
    let bad = "read1\t0\tchr1\n";
    assert!(sam::read_str(bad).is_err());
}

#[test]
fn rejects_seq_qual_length_mismatch() {
    // SEQ has 8 bases, QUAL has 4 — SAM v1 §1.4 says these MUST
    // match when neither is "*". Without the spec check, downstream
    // `seq.bytes().zip(qual.bytes())` silently truncates to 4 and
    // drops half the alignment.
    let bad = "read1\t0\tchr1\t100\t60\t8M\t*\t0\t0\tACGTACGT\tIIII\n";
    let err = sam::read_str(bad).unwrap_err();
    match err {
        sam::SamError::Malformed(msg) => {
            assert!(
                msg.contains("SEQ length 8 != QUAL length 4"),
                "expected length-mismatch message, got: {msg}"
            );
        }
        other => panic!("expected SamError::Malformed, got {other:?}"),
    }
}

#[test]
fn star_seq_or_qual_skips_length_check() {
    // Either side "*" means "unspecified" — must NOT trigger the
    // length-mismatch spec check.
    let star_qual = "read1\t0\tchr1\t100\t60\t8M\t*\t0\t0\tACGTACGT\t*\n";
    let star_seq = "read2\t0\tchr1\t100\t60\t8M\t*\t0\t0\t*\tIIII\n";
    let star_both = "read3\t0\tchr1\t100\t60\t8M\t*\t0\t0\t*\t*\n";
    sam::read_str(star_qual).expect("star qual ok");
    sam::read_str(star_seq).expect("star seq ok");
    sam::read_str(star_both).expect("star both ok");
}
