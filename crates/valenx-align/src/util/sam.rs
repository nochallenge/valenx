//! SAM CIGAR ⇄ [`Alignment`] conversion helpers.
//!
//! The SAM/BAM format encodes an alignment as a *reference start
//! position* plus a *CIGAR string*. This module bridges that compact
//! representation and the crate's explicit two-row
//! [`crate::pairwise::Alignment`]:
//!
//! - [`alignment_to_cigar`] / [`Alignment::cigar`] — gapped rows → CIGAR.
//! - [`cigar_to_rows`] — CIGAR + the two ungapped sequences →
//!   reconstructed gapped rows.
//! - [`SamRecord`] — a minimal SAM alignment record (query name,
//!   reference name, 1-based POS, CIGAR, mapping quality) with a
//!   [`SamRecord::to_sam_line`] tab-delimited serialiser.
//!
//! Convention throughout: **row1 / the query is the read**, **row2 /
//! the target is the reference**. A gap in the reference row is an
//! insertion (`I`), a gap in the read row a deletion (`D`).

use crate::error::{AlignError, Result};
use crate::pairwise::result::{Alignment, Cigar, CigarOp};

/// Builds a SAM CIGAR from an [`Alignment`], treating `row2` as the
/// reference. Equivalent to [`Alignment::cigar`]; provided as a
/// free function for symmetry with [`cigar_to_rows`].
pub fn alignment_to_cigar(alignment: &Alignment) -> Cigar {
    alignment.cigar()
}

/// Reconstructs the two gapped alignment rows from a CIGAR plus the two
/// *ungapped* sequences it was computed over.
///
/// `query` is the read, `reference` the reference; `ref_start` is the
/// 0-based offset into `reference` where the alignment begins (CIGAR
/// `D`/`M` consume reference, `I`/`M` consume query — soft clips
/// consume query only).
///
/// Returns [`AlignError::Dimension`] if the CIGAR consumes more
/// residues than the supplied sequences provide.
pub fn cigar_to_rows(
    cigar: &Cigar,
    query: &[u8],
    reference: &[u8],
    ref_start: usize,
) -> Result<(Vec<u8>, Vec<u8>)> {
    let mut row_query = Vec::new();
    let mut row_ref = Vec::new();
    let mut qi = 0usize;
    let mut ri = ref_start;

    for &(len, op) in &cigar.ops {
        match op {
            CigarOp::Match | CigarOp::Equal | CigarOp::Diff => {
                for _ in 0..len {
                    let q = *query
                        .get(qi)
                        .ok_or_else(|| AlignError::dimension("CIGAR overruns the query"))?;
                    let r = *reference
                        .get(ri)
                        .ok_or_else(|| AlignError::dimension("CIGAR overruns the reference"))?;
                    row_query.push(q);
                    row_ref.push(r);
                    qi += 1;
                    ri += 1;
                }
            }
            CigarOp::Ins => {
                for _ in 0..len {
                    let q = *query
                        .get(qi)
                        .ok_or_else(|| AlignError::dimension("CIGAR insertion overruns the query"))?;
                    row_query.push(q);
                    row_ref.push(b'-');
                    qi += 1;
                }
            }
            CigarOp::Del => {
                for _ in 0..len {
                    let r = *reference
                        .get(ri)
                        .ok_or_else(|| AlignError::dimension("CIGAR deletion overruns the reference"))?;
                    row_query.push(b'-');
                    row_ref.push(r);
                    ri += 1;
                }
            }
            CigarOp::SoftClip => {
                // Soft-clipped query residues are present but unaligned;
                // they do not appear in the alignment rows.
                qi += len;
            }
        }
    }
    Ok((row_query, row_ref))
}

/// Rebuilds a full [`Alignment`] from a CIGAR and the two ungapped
/// sequences, scoring it with a substitution-free identity count
/// (`+1` per match column, `0` otherwise) — enough for round-tripping
/// and inspection. For a properly *scored* alignment, run a pairwise
/// routine instead.
pub fn alignment_from_cigar(
    cigar: &Cigar,
    query: &[u8],
    reference: &[u8],
    ref_start: usize,
) -> Result<Alignment> {
    let (row_q, row_r) = cigar_to_rows(cigar, query, reference, ref_start)?;
    let score = row_q
        .iter()
        .zip(&row_r)
        .filter(|(&a, &b)| a != b'-' && b != b'-' && a.eq_ignore_ascii_case(&b))
        .count() as i32;
    let q_consumed = cigar.query_len() - soft_clip_len(cigar);
    let r_consumed = cigar.ref_len();
    Alignment::new(
        row_q,
        row_r,
        score,
        (0, q_consumed),
        (ref_start, ref_start + r_consumed),
    )
}

/// Total number of soft-clipped query residues in a CIGAR.
fn soft_clip_len(cigar: &Cigar) -> usize {
    cigar
        .ops
        .iter()
        .filter(|(_, op)| *op == CigarOp::SoftClip)
        .map(|(l, _)| l)
        .sum()
}

/// A minimal SAM alignment record — the fields a read mapper needs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SamRecord {
    /// Query (read) name — SAM column 1 `QNAME`.
    pub qname: String,
    /// Bitwise flag — SAM column 2 `FLAG`. `4` means unmapped.
    pub flag: u16,
    /// Reference name — SAM column 3 `RNAME` (`*` if unmapped).
    pub rname: String,
    /// 1-based leftmost reference position — SAM column 4 `POS`
    /// (`0` if unmapped).
    pub pos: usize,
    /// Mapping quality — SAM column 5 `MAPQ`.
    pub mapq: u8,
    /// CIGAR string — SAM column 6 (`*` if unmapped).
    pub cigar: Cigar,
    /// Read sequence — SAM column 10 `SEQ`.
    pub seq: Vec<u8>,
}

impl SamRecord {
    /// Builds a mapped SAM record. `pos` is **0-based** here and is
    /// converted to SAM's 1-based `POS` on output.
    pub fn mapped(
        qname: impl Into<String>,
        rname: impl Into<String>,
        pos_zero_based: usize,
        mapq: u8,
        cigar: Cigar,
        seq: Vec<u8>,
    ) -> Self {
        SamRecord {
            qname: qname.into(),
            flag: 0,
            rname: rname.into(),
            pos: pos_zero_based + 1,
            mapq,
            cigar,
            seq,
        }
    }

    /// Builds an unmapped SAM record (flag `4`, `RNAME`/`CIGAR` `*`).
    pub fn unmapped(qname: impl Into<String>, seq: Vec<u8>) -> Self {
        SamRecord {
            qname: qname.into(),
            flag: 4,
            rname: "*".into(),
            pos: 0,
            mapq: 0,
            cigar: Cigar::new(),
            seq,
        }
    }

    /// `true` if the record is unmapped (flag bit `0x4` set).
    pub fn is_unmapped(&self) -> bool {
        self.flag & 0x4 != 0
    }

    /// Serialises the record as a tab-delimited SAM alignment line
    /// (the 11 mandatory columns; `RNEXT`/`PNEXT`/`TLEN` are `*`/`0`,
    /// `QUAL` is `*`).
    pub fn to_sam_line(&self) -> String {
        let seq = if self.seq.is_empty() {
            "*".to_string()
        } else {
            String::from_utf8_lossy(&self.seq).into_owned()
        };
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t*",
            self.qname,
            self.flag,
            self.rname,
            self.pos,
            self.mapq,
            self.cigar,
            seq,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alignment_cigar_roundtrip() {
        // read AC-GT vs ref ACTGT: 2M 1D 2M.
        let al = Alignment::new(
            b"AC-GT".to_vec(),
            b"ACTGT".to_vec(),
            0,
            (0, 4),
            (0, 5),
        )
        .unwrap();
        let cigar = alignment_to_cigar(&al);
        assert_eq!(cigar.to_string(), "2M1D2M");

        // Reconstruct the rows from the CIGAR + ungapped sequences.
        let (rq, rr) = cigar_to_rows(&cigar, b"ACGT", b"ACTGT", 0).unwrap();
        assert_eq!(rq, b"AC-GT");
        assert_eq!(rr, b"ACTGT");
    }

    #[test]
    fn cigar_with_insertion() {
        // read ACTGT vs ref AC-GT: 2M 1I 2M.
        let cigar = Cigar::parse("2M1I2M").unwrap();
        let (rq, rr) = cigar_to_rows(&cigar, b"ACTGT", b"ACGT", 0).unwrap();
        assert_eq!(rq, b"ACTGT");
        assert_eq!(rr, b"AC-GT");
    }

    #[test]
    fn cigar_with_ref_offset() {
        // Alignment starts 3 residues into the reference.
        let cigar = Cigar::parse("4M").unwrap();
        let (rq, rr) = cigar_to_rows(&cigar, b"ACGT", b"TTTACGTAA", 3).unwrap();
        assert_eq!(rq, b"ACGT");
        assert_eq!(rr, b"ACGT");
    }

    #[test]
    fn soft_clip_skips_query_residues() {
        // 2S4M: first two read residues clipped, four aligned.
        let cigar = Cigar::parse("2S4M").unwrap();
        let (rq, rr) = cigar_to_rows(&cigar, b"NNACGT", b"ACGT", 0).unwrap();
        assert_eq!(rq, b"ACGT"); // the clipped NN do not appear
        assert_eq!(rr, b"ACGT");
    }

    #[test]
    fn cigar_overrun_is_error() {
        let cigar = Cigar::parse("10M").unwrap();
        // Only 4 residues available -> dimension error.
        assert!(cigar_to_rows(&cigar, b"ACGT", b"ACGT", 0).is_err());
    }

    #[test]
    fn alignment_from_cigar_scores_identities() {
        let cigar = Cigar::parse("4M").unwrap();
        let al = alignment_from_cigar(&cigar, b"ACGT", b"ACTT", 0).unwrap();
        // 3 of 4 columns identical.
        assert_eq!(al.score, 3);
        assert_eq!(al.span2, (0, 4));
    }

    #[test]
    fn sam_record_mapped_line() {
        let cigar = Cigar::parse("8M").unwrap();
        let rec = SamRecord::mapped("read1", "chr1", 41, 60, cigar, b"ACGTACGT".to_vec());
        let line = rec.to_sam_line();
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols[0], "read1");
        assert_eq!(cols[1], "0"); // mapped
        assert_eq!(cols[2], "chr1");
        assert_eq!(cols[3], "42"); // 0-based 41 -> 1-based 42
        assert_eq!(cols[4], "60");
        assert_eq!(cols[5], "8M");
        assert_eq!(cols[9], "ACGTACGT");
        assert!(!rec.is_unmapped());
    }

    #[test]
    fn sam_record_unmapped() {
        let rec = SamRecord::unmapped("read2", b"ACGT".to_vec());
        assert!(rec.is_unmapped());
        let line = rec.to_sam_line();
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols[1], "4");
        assert_eq!(cols[2], "*");
        assert_eq!(cols[3], "0");
        assert_eq!(cols[5], "*");
    }
}
