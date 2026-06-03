use valenx_bio::format::fastq;
use valenx_bio::FastqRecord;

const SAMPLE: &str = "\
@read1
ACGTACGT
+
IIIIIIII
@read2 desc here
NNNNAAAA
+
!!!!IIII
";

#[test]
fn parses_two_records() {
    let recs = fastq::read_str(SAMPLE).expect("parse");
    assert_eq!(recs.len(), 2);
    assert_eq!(recs[0].name, "read1");
    assert_eq!(recs[0].sequence, b"ACGTACGT");
    assert_eq!(recs[0].quality, b"IIIIIIII");
    assert_eq!(recs[1].name, "read2");
    assert_eq!(recs[1].description.as_deref(), Some("desc here"));
}

#[test]
fn round_trip_preserves_records() {
    let recs = fastq::read_str(SAMPLE).expect("parse");
    let out = fastq::write_string(&recs).expect("write");
    let recs2 = fastq::read_str(&out).expect("re-parse");
    assert_eq!(recs.len(), recs2.len());
    for (a, b) in recs.iter().zip(recs2.iter()) {
        assert_eq!(a.name, b.name);
        assert_eq!(a.sequence, b.sequence);
        assert_eq!(a.quality, b.quality);
    }
}

#[test]
fn rejects_quality_length_mismatch() {
    let bad = "@r1\nACGT\n+\nII\n";
    assert!(fastq::read_str(bad).is_err());
}

#[test]
fn rejects_missing_plus_separator() {
    let bad = "@r1\nACGT\nNOT-PLUS\nIIII\n";
    assert!(fastq::read_str(bad).is_err());
}

#[test]
fn record_phred_min_max_helpers() {
    let rec = FastqRecord::new("r".into(), None, b"ACGT".to_vec(), b"!*FI".to_vec()).unwrap();
    // phred + 33 offset: '!' = 0, 'I' = 40
    assert_eq!(rec.min_quality(), Some(0));
    assert_eq!(rec.max_quality(), Some(40));
}
