use valenx_bio::format::vcf;
use valenx_bio::VcfRecord;

const SAMPLE: &str = "\
##fileformat=VCFv4.2
##source=test
##contig=<ID=chr1,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1
chr1\t100\trs1\tA\tG\t40.0\tPASS\tDP=20\tGT:DP\t0/1:20
chr1\t200\t.\tA\tT,C\t.\t.\tDP=10\tGT:DP\t1/2:10
chr1\t300\trs3\tACGT\tA\t99.5\tLowQual;Strand\tDP=5\tGT:DP\t1/1:5
";

#[test]
fn parses_header_samples_and_records() {
    let v = vcf::read_str(SAMPLE).expect("parse");
    assert_eq!(v.header.len(), 3);
    assert_eq!(v.samples, vec!["sample1".to_string()]);
    assert_eq!(v.records.len(), 3);
    let r0: &VcfRecord = &v.records[0];
    assert_eq!(r0.chrom, "chr1");
    assert_eq!(r0.pos, 100);
    assert_eq!(r0.id.as_deref(), Some("rs1"));
    assert_eq!(r0.alt, vec!["G".to_string()]);
    assert_eq!(r0.qual, Some(40.0));
    assert_eq!(r0.filter, vec!["PASS".to_string()]);
    assert!(r0.is_pass());
    assert!(r0.has_alt());
}

#[test]
fn handles_dot_for_missing_fields() {
    let v = vcf::read_str(SAMPLE).expect("parse");
    let r1 = &v.records[1];
    assert_eq!(r1.id, None);
    assert_eq!(r1.qual, None);
    assert!(r1.filter.is_empty());
    assert!(r1.is_pass()); // empty filter == unfiltered
    assert_eq!(r1.alt, vec!["T".to_string(), "C".to_string()]); // comma split
}

#[test]
fn parses_multiple_filters() {
    let v = vcf::read_str(SAMPLE).expect("parse");
    let r2 = &v.records[2];
    assert_eq!(r2.filter, vec!["LowQual".to_string(), "Strand".to_string()]);
    assert!(!r2.is_pass());
}

#[test]
fn empty_vcf_parses_to_empty_struct() {
    let v = vcf::read_str("").expect("parse");
    assert!(v.header.is_empty());
    assert!(v.samples.is_empty());
    assert!(v.records.is_empty());
}

#[test]
fn rejects_record_with_too_few_fields() {
    let bad = "chr1\t100\t.\tA\n";
    assert!(vcf::read_str(bad).is_err());
}

#[test]
fn rejects_non_integer_pos() {
    let bad = "chr1\tnot_a_pos\t.\tA\tG\t.\t.\t.\n";
    assert!(vcf::read_str(bad).is_err());
}

#[test]
fn rejects_header_with_no_samples_but_data_has_them() {
    // #CHROM header lacks sample names (only the 8 mandatory + FORMAT,
    // so `Vcf.samples` stays empty) but the data row supplies FORMAT
    // plus a sample column. Without the spec check, the sample
    // string lands at index 0 of `record.samples` with `Vcf.samples`
    // still empty — the consumer reads "sample[0]" with no header
    // label and silently mis-assigns.
    let bad = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT
chr1\t100\trs1\tA\tG\t40.0\tPASS\tDP=20\tGT:DP\t0/1:20
";
    let err = vcf::read_str(bad).unwrap_err();
    assert!(
        matches!(
            err,
            valenx_bio::vcf::VcfError::SampleCountMismatch {
                expected: 0,
                got: 1,
                ..
            }
        ),
        "expected SampleCountMismatch, got {err:?}"
    );
}

#[test]
fn rejects_header_with_samples_but_data_missing_one() {
    // Header advertises two samples; data row only supplies one.
    let bad = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2
chr1\t100\trs1\tA\tG\t40.0\tPASS\tDP=20\tGT:DP\t0/1:20
";
    let err = vcf::read_str(bad).unwrap_err();
    assert!(
        matches!(
            err,
            valenx_bio::vcf::VcfError::SampleCountMismatch {
                expected: 2,
                got: 1,
                ..
            }
        ),
        "expected SampleCountMismatch, got {err:?}"
    );
}
