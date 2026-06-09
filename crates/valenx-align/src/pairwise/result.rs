//! The [`Alignment`] result type — aligned rows, score, CIGAR, stats.
//!
//! Every pairwise routine in [`crate::pairwise`] returns an
//! `Alignment`. It stores the two gapped rows as raw bytes (gaps are
//! the ASCII `-`), the integer score, and the half-open coordinate
//! span each original sequence contributed. From those it derives a
//! [`Cigar`] string, identity / similarity statistics
//! ([`AlignStats`]), and a human-readable [`pretty`](Alignment::pretty)
//! block.

use crate::error::{AlignError, Result};
use crate::matrix::SubstitutionMatrix;
use std::fmt;

/// One operation in a CIGAR string.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum CigarOp {
    /// Aligned column — a residue in both sequences (`M`, "alignment
    /// match", which per the SAM spec covers both identities and
    /// substitutions).
    Match,
    /// Sequence-level identity (`=`) — a residue in both, equal.
    Equal,
    /// Sequence-level mismatch (`X`) — a residue in both, different.
    Diff,
    /// Insertion relative to the reference (`I`) — a residue in the
    /// query, a gap in the reference / target.
    Ins,
    /// Deletion relative to the reference (`D`) — a gap in the query,
    /// a residue in the reference / target.
    Del,
    /// Soft clip (`S`) — query residues present but not aligned.
    SoftClip,
}

impl CigarOp {
    /// The single-character SAM code for this op.
    pub fn code(self) -> char {
        match self {
            CigarOp::Match => 'M',
            CigarOp::Equal => '=',
            CigarOp::Diff => 'X',
            CigarOp::Ins => 'I',
            CigarOp::Del => 'D',
            CigarOp::SoftClip => 'S',
        }
    }

    /// Parses a SAM CIGAR op code; `None` for an unrecognised char.
    pub fn from_code(c: char) -> Option<Self> {
        Some(match c {
            'M' => CigarOp::Match,
            '=' => CigarOp::Equal,
            'X' => CigarOp::Diff,
            'I' => CigarOp::Ins,
            'D' => CigarOp::Del,
            'S' => CigarOp::SoftClip,
            _ => return None,
        })
    }

    /// `true` if this op consumes a residue from the query sequence.
    pub fn consumes_query(self) -> bool {
        matches!(
            self,
            CigarOp::Match | CigarOp::Equal | CigarOp::Diff | CigarOp::Ins | CigarOp::SoftClip
        )
    }

    /// `true` if this op consumes a residue from the reference.
    pub fn consumes_ref(self) -> bool {
        matches!(
            self,
            CigarOp::Match | CigarOp::Equal | CigarOp::Diff | CigarOp::Del
        )
    }
}

/// A run-length-encoded CIGAR string: a list of `(length, op)` pairs.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Cigar {
    /// The `(run-length, operation)` pairs in left-to-right order.
    pub ops: Vec<(usize, CigarOp)>,
}

impl Cigar {
    /// An empty CIGAR.
    pub fn new() -> Self {
        Cigar::default()
    }

    /// Appends one op of length 1, merging with the previous run if the
    /// op matches.
    pub fn push(&mut self, op: CigarOp) {
        match self.ops.last_mut() {
            Some((len, last)) if *last == op => *len += 1,
            _ => self.ops.push((1, op)),
        }
    }

    /// Total number of reference residues consumed.
    pub fn ref_len(&self) -> usize {
        self.ops
            .iter()
            .filter(|(_, op)| op.consumes_ref())
            .map(|(l, _)| l)
            .sum()
    }

    /// Total number of query residues consumed.
    pub fn query_len(&self) -> usize {
        self.ops
            .iter()
            .filter(|(_, op)| op.consumes_query())
            .map(|(l, _)| l)
            .sum()
    }

    /// Parses a SAM CIGAR string (`"10M2I5M"`).
    pub fn parse(s: &str) -> Result<Self> {
        let mut ops = Vec::new();
        let mut num = String::new();
        for c in s.chars() {
            if c.is_ascii_digit() {
                num.push(c);
            } else {
                let len: usize = num.parse().map_err(|_| {
                    AlignError::parse("cigar", format!("bad run length before `{c}`"))
                })?;
                let op = CigarOp::from_code(c)
                    .ok_or_else(|| AlignError::parse("cigar", format!("unknown op `{c}`")))?;
                if len > 0 {
                    ops.push((len, op));
                }
                num.clear();
            }
        }
        if !num.is_empty() {
            return Err(AlignError::parse("cigar", "trailing run length"));
        }
        Ok(Cigar { ops })
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.ops.is_empty() {
            return write!(f, "*");
        }
        for (len, op) in &self.ops {
            write!(f, "{}{}", len, op.code())?;
        }
        Ok(())
    }
}

/// Identity / similarity / gap statistics over an [`Alignment`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct AlignStats {
    /// Number of aligned columns (gaps included).
    pub columns: usize,
    /// Columns where both residues are identical.
    pub identities: usize,
    /// Columns where both residues align with a positive matrix score
    /// (a superset of identities — the "positives" line in BLAST).
    pub similarities: usize,
    /// Columns containing a gap in either row.
    pub gaps: usize,
    /// Number of distinct gap runs (gap-opening events).
    pub gap_opens: usize,
}

impl AlignStats {
    /// Fraction of non-gap columns that are identical, in `[0, 1]`.
    /// Returns `0.0` when there are no aligned (non-gap) columns.
    pub fn percent_identity(&self) -> f64 {
        let aligned = self.columns - self.gaps;
        if aligned == 0 {
            0.0
        } else {
            self.identities as f64 / aligned as f64
        }
    }

    /// Fraction of non-gap columns that are positive-scoring, `[0, 1]`.
    pub fn percent_similarity(&self) -> f64 {
        let aligned = self.columns - self.gaps;
        if aligned == 0 {
            0.0
        } else {
            self.similarities as f64 / aligned as f64
        }
    }

    /// Fraction of columns containing a gap, in `[0, 1]`.
    pub fn percent_gaps(&self) -> f64 {
        if self.columns == 0 {
            0.0
        } else {
            self.gaps as f64 / self.columns as f64
        }
    }

    /// Number of aligned (non-gap) columns — the alignment length after removing gap
    /// columns (`columns − gaps`).
    pub fn aligned_length(&self) -> usize {
        self.columns.saturating_sub(self.gaps)
    }

    /// Number of aligned (non-gap) columns whose residues differ — the complement of the
    /// identities within the aligned region, so `identities + mismatch_count ==
    /// aligned_length`.
    pub fn mismatch_count(&self) -> usize {
        self.aligned_length().saturating_sub(self.identities)
    }
}

/// A pairwise alignment result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alignment {
    /// First sequence, gapped (ASCII; `-` for gaps).
    pub row1: Vec<u8>,
    /// Second sequence, gapped (same length as `row1`).
    pub row2: Vec<u8>,
    /// The optimal alignment score under the scoring scheme used.
    pub score: i32,
    /// Half-open `[start, end)` span of the *first* sequence covered.
    pub span1: (usize, usize),
    /// Half-open `[start, end)` span of the *second* sequence covered.
    pub span2: (usize, usize),
}

impl Alignment {
    /// Builds an alignment, validating that both rows are the same
    /// length. Returns [`AlignError::Dimension`] otherwise.
    pub fn new(
        row1: Vec<u8>,
        row2: Vec<u8>,
        score: i32,
        span1: (usize, usize),
        span2: (usize, usize),
    ) -> Result<Self> {
        if row1.len() != row2.len() {
            return Err(AlignError::dimension(format!(
                "alignment rows differ: {} vs {}",
                row1.len(),
                row2.len()
            )));
        }
        Ok(Alignment {
            row1,
            row2,
            score,
            span1,
            span2,
        })
    }

    /// Number of aligned columns.
    pub fn len(&self) -> usize {
        self.row1.len()
    }

    /// `true` if the alignment has no columns.
    pub fn is_empty(&self) -> bool {
        self.row1.is_empty()
    }

    /// The first row as a `&str` (rows are always ASCII).
    pub fn row1_str(&self) -> &str {
        std::str::from_utf8(&self.row1).unwrap_or("<non-utf8>")
    }

    /// The second row as a `&str`.
    pub fn row2_str(&self) -> &str {
        std::str::from_utf8(&self.row2).unwrap_or("<non-utf8>")
    }

    /// The CIGAR string with `row2` treated as the reference: a gap in
    /// `row2` is an insertion (`I`), a gap in `row1` a deletion (`D`),
    /// and an aligned column an `M`.
    pub fn cigar(&self) -> Cigar {
        let mut c = Cigar::new();
        for (&a, &b) in self.row1.iter().zip(&self.row2) {
            let op = match (a == b'-', b == b'-') {
                (false, false) => CigarOp::Match,
                (true, false) => CigarOp::Del,
                (false, true) => CigarOp::Ins,
                (true, true) => continue, // gap/gap column: skip
            };
            c.push(op);
        }
        c
    }

    /// An extended CIGAR distinguishing identities (`=`) from
    /// substitutions (`X`).
    pub fn cigar_extended(&self) -> Cigar {
        let mut c = Cigar::new();
        for (&a, &b) in self.row1.iter().zip(&self.row2) {
            let op = match (a == b'-', b == b'-') {
                (false, false) => {
                    if a.eq_ignore_ascii_case(&b) {
                        CigarOp::Equal
                    } else {
                        CigarOp::Diff
                    }
                }
                (true, false) => CigarOp::Del,
                (false, true) => CigarOp::Ins,
                (true, true) => continue,
            };
            c.push(op);
        }
        c
    }

    /// Identity / similarity / gap statistics. `matrix` defines which
    /// substitutions count as "similar" (positive score).
    pub fn stats(&self, matrix: &SubstitutionMatrix) -> AlignStats {
        let mut s = AlignStats {
            columns: self.len(),
            identities: 0,
            similarities: 0,
            gaps: 0,
            gap_opens: 0,
        };
        let mut in_gap = false;
        for (&a, &b) in self.row1.iter().zip(&self.row2) {
            let is_gap = a == b'-' || b == b'-';
            if is_gap {
                s.gaps += 1;
                if !in_gap {
                    s.gap_opens += 1;
                }
                in_gap = true;
            } else {
                in_gap = false;
                if a.eq_ignore_ascii_case(&b) {
                    s.identities += 1;
                    s.similarities += 1;
                } else if matrix.score(a, b) > 0 {
                    s.similarities += 1;
                }
            }
        }
        s
    }

    /// Fraction of non-gap columns that are identical, in `[0, 1]`.
    /// A convenience for [`AlignStats::percent_identity`] using an
    /// identity matrix (substitution scores are irrelevant to it).
    pub fn percent_identity(&self) -> f64 {
        self.stats(&SubstitutionMatrix::identity(1, -1))
            .percent_identity()
    }

    /// A human-readable multi-line rendering. `width` is the residues
    /// per block (e.g. `60`); a `|` marks identities and a `+` marks
    /// positive-scoring substitutions on the middle "match" line.
    pub fn pretty(&self, width: usize, matrix: &SubstitutionMatrix) -> String {
        let width = width.max(1);
        let mut out = String::new();
        let mut pos = 0;
        while pos < self.len() {
            let end = (pos + width).min(self.len());
            let r1 = &self.row1[pos..end];
            let r2 = &self.row2[pos..end];
            let mid: String = r1
                .iter()
                .zip(r2)
                .map(|(&a, &b)| {
                    if a == b'-' || b == b'-' {
                        ' '
                    } else if a.eq_ignore_ascii_case(&b) {
                        '|'
                    } else if matrix.score(a, b) > 0 {
                        '+'
                    } else {
                        '.'
                    }
                })
                .collect();
            out.push_str(&format!("seq1 {}\n", String::from_utf8_lossy(r1)));
            out.push_str(&format!("     {mid}\n"));
            out.push_str(&format!("seq2 {}\n", String::from_utf8_lossy(r2)));
            if end < self.len() {
                out.push('\n');
            }
            pos = end;
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigar_roundtrip() {
        let c = Cigar::parse("10M2I5M3D").unwrap();
        assert_eq!(c.ops.len(), 4);
        assert_eq!(c.to_string(), "10M2I5M3D");
        assert_eq!(c.ref_len(), 18); // 10 + 5 + 3
        assert_eq!(c.query_len(), 17); // 10 + 2 + 5
    }

    #[test]
    fn cigar_rejects_garbage() {
        assert!(Cigar::parse("10Z").is_err());
        assert!(Cigar::parse("M").is_err());
        assert!(Cigar::parse("5").is_err());
        assert_eq!(Cigar::new().to_string(), "*");
    }

    #[test]
    fn alignment_rows_must_match() {
        assert!(Alignment::new(b"ACGT".to_vec(), b"ACG".to_vec(), 0, (0, 4), (0, 3)).is_err());
        let a = Alignment::new(b"ACGT".to_vec(), b"ACGT".to_vec(), 8, (0, 4), (0, 4)).unwrap();
        assert_eq!(a.len(), 4);
    }

    #[test]
    fn cigar_from_alignment() {
        // row1 = AC-GT, row2 = ACTGT -> M M D M M  =>  2M1D2M
        let a = Alignment::new(
            b"AC-GT".to_vec(),
            b"ACTGT".to_vec(),
            0,
            (0, 4),
            (0, 5),
        )
        .unwrap();
        assert_eq!(a.cigar().to_string(), "2M1D2M");
        // row1 has a residue where row2 has a gap -> insertion
        let b = Alignment::new(
            b"ACTGT".to_vec(),
            b"AC-GT".to_vec(),
            0,
            (0, 5),
            (0, 4),
        )
        .unwrap();
        assert_eq!(b.cigar().to_string(), "2M1I2M");
    }

    #[test]
    fn extended_cigar_distinguishes_mismatch() {
        // AAGT / AACT -> = = X =
        let a = Alignment::new(b"AAGT".to_vec(), b"AACT".to_vec(), 0, (0, 4), (0, 4)).unwrap();
        assert_eq!(a.cigar_extended().to_string(), "2=1X1=");
    }

    #[test]
    fn stats_and_identity() {
        // 4 identical, 1 mismatch, 1 gap column
        let a = Alignment::new(
            b"ACGTA-".to_vec(),
            b"ACGTTC".to_vec(),
            0,
            (0, 5),
            (0, 6),
        )
        .unwrap();
        let m = SubstitutionMatrix::identity(1, -1);
        let s = a.stats(&m);
        assert_eq!(s.columns, 6);
        assert_eq!(s.identities, 4);
        assert_eq!(s.gaps, 1);
        assert_eq!(s.gap_opens, 1);
        // 4 identical of 5 aligned columns
        assert!((s.percent_identity() - 0.8).abs() < 1e-9);
        assert!((a.percent_identity() - 0.8).abs() < 1e-9);
    }

    #[test]
    fn aligned_length_and_mismatch_count() {
        // "ACGTA-" vs "ACGTTC": 6 columns, 4 identities, 1 mismatch (col 4), 1 gap (col 5).
        let a = Alignment::new(b"ACGTA-".to_vec(), b"ACGTTC".to_vec(), 0, (0, 5), (0, 6)).unwrap();
        let s = a.stats(&SubstitutionMatrix::identity(1, -1));
        // aligned_length = columns − gaps = 6 − 1 = 5.
        assert_eq!(s.aligned_length(), 5);
        // mismatch_count = aligned − identities = 5 − 4 = 1.
        assert_eq!(s.mismatch_count(), 1);
        // partition invariant: identities + mismatches = aligned length.
        assert_eq!(s.identities + s.mismatch_count(), s.aligned_length());
    }

    #[test]
    fn pretty_print_runs() {
        let a = Alignment::new(b"ACGT".to_vec(), b"ACTT".to_vec(), 0, (0, 4), (0, 4)).unwrap();
        let m = SubstitutionMatrix::identity(1, -1);
        let p = a.pretty(60, &m);
        assert!(p.contains("seq1 ACGT"));
        assert!(p.contains("seq2 ACTT"));
        assert!(p.contains("||")); // at least two identities marked
    }
}
