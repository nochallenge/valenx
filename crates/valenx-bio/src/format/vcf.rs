//! Minimal VCF text reader. Handles 8 mandatory columns + optional
//! FORMAT + per-sample columns. BCF (binary) and bgzf-compressed VCF
//! are out of scope — convert with `bcftools view` first.
//!
//! `"."` is the missing-value sentinel for ID / QUAL / ALT / FILTER.
//! INFO is special: `"."` there means "no annotations" and rides
//! through as the literal `"."` string.

use crate::vcf::{Vcf, VcfError, VcfRecord};

/// Parse a VCF text document.
///
/// Empty lines are skipped. `##` lines are collected verbatim into
/// [`Vcf::header`]. The `#CHROM` line populates [`Vcf::samples`] from
/// columns 9 onward (after the 8 mandatory + FORMAT columns).
///
/// Records with fewer than 8 TAB-separated fields are an error. Once
/// a `#CHROM` header has been seen, every data row's per-sample column
/// count must match the header's sample count exactly — otherwise the
/// per-sample arrays would silently shift onto the wrong subject. A
/// `#CHROM` line with fewer than 10 columns followed by data rows with
/// sample columns triggers [`VcfError::SampleCountMismatch`].
pub fn read_str(s: &str) -> Result<Vcf, VcfError> {
    let mut vcf = Vcf::default();
    // Whether we've seen the `#CHROM` line. Once true, every data
    // row's per-sample column count is checked against
    // `vcf.samples.len()`. Pre-`#CHROM` data rows (rare) are not
    // validated because we don't yet know the expected count.
    let mut chrom_seen = false;
    for (i, line) in s.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        if let Some(meta) = line.strip_prefix("##") {
            vcf.header.push(format!("##{meta}"));
            continue;
        }
        if let Some(cols) = line.strip_prefix('#') {
            // #CHROM line — column header. Pull samples from
            // index 9 onward (after the 8 mandatory + FORMAT).
            let fields: Vec<&str> = cols.split('\t').collect();
            if fields.len() >= 10 {
                vcf.samples = fields[9..].iter().map(|s| s.to_string()).collect();
            }
            chrom_seen = true;
            continue;
        }
        // Data row.
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            return Err(VcfError::Bad {
                line: i + 1,
                msg: format!("expected at least 8 fields, got {}", fields.len()),
            });
        }
        let line_no = i + 1;
        let pos: u64 = fields[1].parse().map_err(|_| VcfError::Bad {
            line: line_no,
            msg: format!("POS not an integer: `{}`", fields[1]),
        })?;
        let id = match fields[2] {
            "." | "" => None,
            s => Some(s.to_string()),
        };
        let alt: Vec<String> = match fields[4] {
            "." | "" => Vec::new(),
            s => s.split(',').map(|x| x.to_string()).collect(),
        };
        let qual = match fields[5] {
            "." | "" => None,
            s => Some(s.parse().map_err(|_| VcfError::Bad {
                line: line_no,
                msg: format!("QUAL not a number: `{s}`"),
            })?),
        };
        let filter: Vec<String> = match fields[6] {
            "." | "" => Vec::new(),
            s => s.split(';').map(|x| x.to_string()).collect(),
        };
        let info = fields[7].to_string();
        let format = fields.get(8).map(|s| s.to_string());
        let samples: Vec<String> = if fields.len() > 9 {
            fields[9..].iter().map(|s| s.to_string()).collect()
        } else {
            Vec::new()
        };
        // Header/data drift check: once `#CHROM` is parsed, every
        // data row's sample count must match `vcf.samples.len()`. The
        // canonical buggy case is a header with no sample names
        // (fields.len() < 10, so `vcf.samples` stays empty) followed
        // by data rows that DO carry sample columns; the per-sample
        // arrays would otherwise be assigned to non-existent subjects.
        if chrom_seen && samples.len() != vcf.samples.len() {
            return Err(VcfError::SampleCountMismatch {
                line: line_no,
                expected: vcf.samples.len(),
                got: samples.len(),
            });
        }
        vcf.records.push(VcfRecord {
            chrom: fields[0].to_string(),
            pos,
            id,
            ref_allele: fields[3].to_string(),
            alt,
            qual,
            filter,
            info,
            format,
            samples,
        });
    }
    Ok(vcf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn happy_path_parses_header_and_one_record() {
        let vcf_text = "##fileformat=VCFv4.2\n\
                        ##source=test\n\
                        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                        chr1\t100\t.\tA\tG\t.\tPASS\t.\n";
        let vcf = read_str(vcf_text).unwrap();
        assert_eq!(vcf.header.len(), 2);
        assert_eq!(vcf.samples.len(), 0);
        assert_eq!(vcf.records.len(), 1);
        let r = &vcf.records[0];
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.pos, 100);
        assert_eq!(r.ref_allele, "A");
        assert_eq!(r.alt, vec!["G".to_string()]);
        assert!(r.is_pass());
    }

    #[test]
    fn parses_per_sample_columns() {
        let vcf_text = "##fileformat=VCFv4.2\n\
                        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA0001\tNA0002\n\
                        chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1\n";
        let vcf = read_str(vcf_text).unwrap();
        assert_eq!(
            vcf.samples,
            vec!["NA0001".to_string(), "NA0002".to_string()]
        );
        assert_eq!(
            vcf.records[0].samples,
            vec!["0/1".to_string(), "1/1".to_string()]
        );
    }

    #[test]
    fn data_row_with_more_samples_than_header_rejected() {
        // Header has 8 cols (no FORMAT, no samples) but data row carries
        // FORMAT + GT for one phantom sample. Round-3 added the validation.
        let vcf = "##fileformat=VCFv4.2\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"; // header 0 samples, row 2 sample-area cols
        let result = read_str(vcf);
        assert!(matches!(result, Err(VcfError::SampleCountMismatch { .. })));
    }
}
