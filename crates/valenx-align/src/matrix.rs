//! Substitution matrices and the affine gap-cost model.
//!
//! This module supplies the *scoring* half of every alignment in the
//! crate. Two pieces:
//!
//! - [`SubstitutionMatrix`] — a residue × residue integer score table.
//!   Built-in tables: the BLOSUM family (45 / 62 / 80), the PAM family
//!   (30 / 70 / 250), a generic [`identity`](SubstitutionMatrix::identity)
//!   match/mismatch matrix, and `NUC.4.4`, the EDNAFULL DNA matrix.
//! - [`GapCost`] — the affine gap penalty (`open` for the first gap
//!   residue, `extend` for each subsequent one) and [`ScoringScheme`],
//!   which bundles a matrix with a gap cost so the DP routines take a
//!   single parameter.
//!
//! ## v1 scope
//!
//! The BLOSUM and PAM tables are the canonical NCBI values for the 20
//! standard amino acids plus `B Z X *`. The `NUC.4.4` table is the
//! EMBOSS `EDNAFULL` matrix over the 15 IUPAC nucleotide codes.
//! Scores are integers (the reference matrices are integer log-odds);
//! this keeps the DP exact and fast. Unknown residues fall back to the
//! `*` (or `X`) row when present, else to a configurable mismatch
//! score.

use crate::error::{AlignError, Result};

/// The amino-acid column order shared by all built-in protein matrices.
///
/// This is the classic BLOSUM/PAM ordering. Index `i` of every row in
/// the embedded tables corresponds to residue `AA_ORDER[i]`.
pub const AA_ORDER: &[u8; 24] = b"ARNDCQEGHILKMFPSTWYVBZX*";

/// The nucleotide column order for the `NUC.4.4` table (EDNAFULL).
pub const NUC_ORDER: &[u8; 15] = b"ATGCSWRYKMBVHDN";

/// An integer residue-substitution score table.
///
/// Constructed from one of the named factory methods. Look up a score
/// with [`score`](Self::score), which is symmetric and case-insensitive
/// and gracefully handles residues outside the table.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SubstitutionMatrix {
    /// Human-readable matrix name (`"BLOSUM62"`, `"NUC.4.4"`, …).
    name: &'static str,
    /// The residue letters indexing rows / columns, uppercased.
    order: Vec<u8>,
    /// `idx[b]` is the row/column for ASCII byte `b`, or `255` if `b`
    /// is not in the table. 256-entry dense map for O(1) lookup.
    idx: [u8; 256],
    /// Flattened `order.len() * order.len()` score grid, row-major.
    scores: Vec<i32>,
    /// Score used when *either* residue is outside the table and no
    /// fallback row exists.
    default_mismatch: i32,
}

impl SubstitutionMatrix {
    /// Builds a matrix from an explicit row order and a flattened
    /// row-major score grid. `scores.len()` must equal `order.len()²`.
    fn from_table(
        name: &'static str,
        order: &[u8],
        scores: Vec<i32>,
        default_mismatch: i32,
    ) -> Self {
        assert_eq!(
            scores.len(),
            order.len() * order.len(),
            "score grid size must be order²"
        );
        let mut idx = [255u8; 256];
        for (i, &b) in order.iter().enumerate() {
            idx[b as usize] = i as u8;
        }
        SubstitutionMatrix {
            name,
            order: order.to_vec(),
            idx,
            scores,
            default_mismatch,
        }
    }

    /// The matrix name.
    pub fn name(&self) -> &'static str {
        self.name
    }

    /// The residue letters indexing the table.
    pub fn alphabet(&self) -> &[u8] {
        &self.order
    }

    /// Score for substituting residue `a` with residue `b`.
    ///
    /// Case-insensitive. If a residue is not in the table the lookup
    /// retries with `X` (proteins) or `N` (nucleotides) when that
    /// fallback row exists; failing that it returns
    /// `default_mismatch`. Gap characters `-` always score
    /// `default_mismatch` (real gap handling is the DP's job, not the
    /// matrix's).
    pub fn score(&self, a: u8, b: u8) -> i32 {
        let au = a.to_ascii_uppercase();
        let bu = b.to_ascii_uppercase();
        let ia = self.resolve(au);
        let ib = self.resolve(bu);
        match (ia, ib) {
            (Some(i), Some(j)) => self.scores[i * self.order.len() + j],
            _ => self.default_mismatch,
        }
    }

    /// Resolves an uppercased residue to a table row, applying the
    /// `X` / `N` fallback for unknown letters.
    fn resolve(&self, u: u8) -> Option<usize> {
        if u == b'-' {
            return None;
        }
        let direct = self.idx[u as usize];
        if direct != 255 {
            return Some(direct as usize);
        }
        // Fallback to the wildcard row if present.
        for &fallback in b"XN" {
            let f = self.idx[fallback as usize];
            if f != 255 {
                return Some(f as usize);
            }
        }
        None
    }

    /// `true` if `residue` (case-insensitive) has its own row.
    pub fn contains(&self, residue: u8) -> bool {
        self.idx[residue.to_ascii_uppercase() as usize] != 255
    }

    /// A generic match/mismatch ("identity") matrix over the 20 amino
    /// acids plus `B Z X *`: `match_score` on the diagonal,
    /// `mismatch_score` everywhere else.
    pub fn identity(match_score: i32, mismatch_score: i32) -> Self {
        let n = AA_ORDER.len();
        let mut scores = vec![mismatch_score; n * n];
        for i in 0..n {
            scores[i * n + i] = match_score;
        }
        Self::from_table("IDENTITY", AA_ORDER, scores, mismatch_score)
    }

    /// A generic nucleotide match/mismatch matrix over `ACGT` only.
    /// Useful when a simple `+1 / -1` DNA score is wanted instead of
    /// the full `NUC.4.4` ambiguity table.
    pub fn dna_simple(match_score: i32, mismatch_score: i32) -> Self {
        let order = b"ACGT";
        let n = order.len();
        let mut scores = vec![mismatch_score; n * n];
        for i in 0..n {
            scores[i * n + i] = match_score;
        }
        Self::from_table("DNA-SIMPLE", order, scores, mismatch_score)
    }

    /// `NUC.4.4` — the EMBOSS `EDNAFULL` DNA matrix over the 15 IUPAC
    /// nucleotide codes. Exact matches score `+5`, transitions and
    /// transversions `-4`, and ambiguity codes get fractional-rounded
    /// partial credit.
    pub fn nuc44() -> Self {
        Self::from_table("NUC.4.4", NUC_ORDER, nuc44_table(), -4)
    }

    /// BLOSUM62 — the default protein matrix for most homology search.
    pub fn blosum62() -> Self {
        Self::from_table("BLOSUM62", AA_ORDER, blosum62_table(), -4)
    }

    /// BLOSUM45 — for more divergent protein sequences than BLOSUM62.
    pub fn blosum45() -> Self {
        Self::from_table("BLOSUM45", AA_ORDER, blosum45_table(), -5)
    }

    /// BLOSUM80 — for closely related protein sequences.
    pub fn blosum80() -> Self {
        Self::from_table("BLOSUM80", AA_ORDER, blosum80_table(), -6)
    }

    /// PAM30 — short, near-identical protein alignment.
    pub fn pam30() -> Self {
        Self::from_table("PAM30", AA_ORDER, pam30_table(), -17)
    }

    /// PAM70 — moderately related protein sequences.
    pub fn pam70() -> Self {
        Self::from_table("PAM70", AA_ORDER, pam70_table(), -11)
    }

    /// PAM250 — distantly related protein sequences.
    pub fn pam250() -> Self {
        Self::from_table("PAM250", AA_ORDER, pam250_table(), -8)
    }

    /// Looks a built-in matrix up by name (case-insensitive). Accepts
    /// `BLOSUM45/62/80`, `PAM30/70/250`, `NUC.4.4` / `NUC44`,
    /// `IDENTITY` and `DNA`.
    pub fn by_name(name: &str) -> Result<Self> {
        Ok(match name.to_ascii_uppercase().as_str() {
            "BLOSUM45" => Self::blosum45(),
            "BLOSUM62" => Self::blosum62(),
            "BLOSUM80" => Self::blosum80(),
            "PAM30" => Self::pam30(),
            "PAM70" => Self::pam70(),
            "PAM250" => Self::pam250(),
            "NUC.4.4" | "NUC44" | "EDNAFULL" => Self::nuc44(),
            "IDENTITY" => Self::identity(1, -1),
            "DNA" | "DNA-SIMPLE" => Self::dna_simple(5, -4),
            other => {
                return Err(AlignError::invalid(
                    "matrix",
                    format!("unknown substitution matrix `{other}`"),
                ))
            }
        })
    }
}

/// An affine gap-cost model.
///
/// A run of `n` consecutive gap residues costs `open + extend * n`
/// (Gotoh convention — the first residue pays `open + extend`). Both
/// fields are stored as non-negative penalties; the DP subtracts them.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct GapCost {
    /// Penalty applied once when a gap is opened.
    pub open: i32,
    /// Penalty applied for every gap residue (including the first).
    pub extend: i32,
}

impl GapCost {
    /// A new affine gap cost. `open` and `extend` are penalties and
    /// should be non-negative; negative values are clamped to `0`.
    pub fn new(open: i32, extend: i32) -> Self {
        GapCost {
            open: open.max(0),
            extend: extend.max(0),
        }
    }

    /// A *linear* gap cost — every gap residue costs `per_residue` and
    /// there is no separate open penalty. Equivalent to
    /// `GapCost::new(0, per_residue)`.
    pub fn linear(per_residue: i32) -> Self {
        GapCost::new(0, per_residue.max(0))
    }

    /// Total penalty for a gap of `len` residues (`len >= 1`); `0` for
    /// `len == 0`.
    pub fn total(&self, len: usize) -> i32 {
        if len == 0 {
            0
        } else {
            self.open + self.extend * len as i32
        }
    }

    /// `true` if this is a pure linear gap cost (no open penalty).
    pub fn is_linear(&self) -> bool {
        self.open == 0
    }
}

impl Default for GapCost {
    /// The common BLOSUM62 default: open `11`, extend `1`.
    fn default() -> Self {
        GapCost::new(11, 1)
    }
}

/// A complete scoring scheme: a substitution matrix plus a gap cost.
///
/// The pairwise / MSA / search routines all take a `&ScoringScheme` so
/// scoring is configured in exactly one place.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ScoringScheme {
    /// The residue substitution matrix.
    pub matrix: SubstitutionMatrix,
    /// The affine gap-cost model.
    pub gap: GapCost,
}

impl ScoringScheme {
    /// Bundles a matrix and a gap cost.
    pub fn new(matrix: SubstitutionMatrix, gap: GapCost) -> Self {
        ScoringScheme { matrix, gap }
    }

    /// The standard protein default: BLOSUM62 with gap open `11`,
    /// extend `1`.
    pub fn blosum62_default() -> Self {
        ScoringScheme::new(SubstitutionMatrix::blosum62(), GapCost::new(11, 1))
    }

    /// A standard DNA default: `NUC.4.4` with gap open `10`, extend `1`.
    pub fn dna_default() -> Self {
        ScoringScheme::new(SubstitutionMatrix::nuc44(), GapCost::new(10, 1))
    }

    /// Substitution score for residues `a` and `b`.
    pub fn sub(&self, a: u8, b: u8) -> i32 {
        self.matrix.score(a, b)
    }
}

// =====================================================================
// Embedded matrix tables
// =====================================================================
//
// All protein tables are 24×24 over `AA_ORDER` ("ARNDCQEGHILKMFPSTWYVBZX*").
// Values are the canonical NCBI integer log-odds scores. The `*` row /
// column carries the standard "-4 vs everything, +1 on the diagonal"
// convention (BLOSUM62) so a stop codon scores like a strong mismatch.

/// Helper: builds a 24×24 protein table from the canonical
/// upper/lower-triangular text form used by NCBI matrix files. `rows`
/// is 24 slices, each `i+1` long (the lower triangle including the
/// diagonal); the matrix is mirrored to fill the upper triangle.
fn protein_table(rows: &[&[i32]]) -> Vec<i32> {
    let n = AA_ORDER.len();
    assert_eq!(rows.len(), n);
    let mut m = vec![0i32; n * n];
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row.len(), i + 1, "row {i} must have {} entries", i + 1);
        for (j, &v) in row.iter().enumerate() {
            m[i * n + j] = v;
            m[j * n + i] = v;
        }
    }
    m
}

/// BLOSUM62 (NCBI canonical, integer log-odds).
fn blosum62_table() -> Vec<i32> {
    protein_table(&[
        &[4],
        &[-1, 5],
        &[-2, 0, 6],
        &[-2, -2, 1, 6],
        &[0, -3, -3, -3, 9],
        &[-1, 1, 0, 0, -3, 5],
        &[-1, 0, 0, 2, -4, 2, 5],
        &[0, -2, 0, -1, -3, -2, -2, 6],
        &[-2, 0, 1, -1, -3, 0, 0, -2, 8],
        &[-1, -3, -3, -3, -1, -3, -3, -4, -3, 4],
        &[-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4],
        &[-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5],
        &[-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5],
        &[-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6],
        &[-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7],
        &[1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4],
        &[0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5],
        &[
            -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11,
        ],
        &[
            -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7,
        ],
        &[
            0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4,
        ],
        &[
            -2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4,
        ],
        &[
            -1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4,
        ],
        &[
            0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1,
        ],
        &[
            -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
            -4, 1,
        ],
    ])
}

/// BLOSUM45 (NCBI canonical).
fn blosum45_table() -> Vec<i32> {
    protein_table(&[
        &[5],
        &[-2, 7],
        &[-1, 0, 6],
        &[-2, -1, 2, 7],
        &[-1, -3, -2, -3, 12],
        &[-1, 1, 0, 0, -3, 6],
        &[-1, 0, 0, 2, -3, 2, 6],
        &[0, -2, 0, -1, -3, -2, -2, 7],
        &[-2, 0, 1, 0, -3, 1, 0, -2, 10],
        &[-1, -3, -2, -4, -3, -2, -3, -4, -3, 5],
        &[-1, -2, -3, -3, -2, -2, -2, -3, -2, 2, 5],
        &[-1, 3, 0, 0, -3, 1, 1, -2, -1, -3, -3, 5],
        &[-1, -1, -2, -3, -2, 0, -2, -2, 0, 2, 2, -1, 6],
        &[-2, -2, -2, -4, -2, -4, -3, -3, -2, 0, 1, -3, 0, 8],
        &[-1, -2, -2, -1, -4, -1, 0, -2, -2, -2, -3, -1, -2, -3, 9],
        &[1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -3, -1, -2, -2, -1, 4],
        &[0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1, 2, 5],
        &[
            -2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2, 1, -3, -4, -3, 15,
        ],
        &[
            -2, -1, -2, -2, -3, -1, -2, -3, 2, 0, 0, -1, 0, 3, -3, -2, -1, 3, 8,
        ],
        &[
            0, -2, -3, -3, -1, -3, -3, -3, -3, 3, 1, -2, 1, 0, -3, -1, 0, -3, -1, 5,
        ],
        &[
            -1, -1, 4, 5, -2, 0, 1, -1, 0, -3, -3, 0, -2, -3, -2, 0, 0, -4, -2, -3, 4,
        ],
        &[
            -1, 0, 0, 1, -3, 4, 4, -2, 0, -3, -2, 1, -1, -3, -1, 0, -1, -2, -2, -3, 2, 4,
        ],
        &[
            -1, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -2, -1, -1, -1, -1, -1,
        ],
        &[
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, 1,
        ],
    ])
}

/// BLOSUM80 (NCBI canonical).
fn blosum80_table() -> Vec<i32> {
    protein_table(&[
        &[7],
        &[-3, 9],
        &[-3, -1, 9],
        &[-3, -3, 2, 10],
        &[-1, -6, -5, -7, 13],
        &[-2, 1, 0, -1, -5, 9],
        &[-2, -1, -1, 2, -7, 3, 8],
        &[0, -4, -1, -3, -6, -4, -4, 9],
        &[-3, 0, 1, -2, -7, 1, 0, -4, 12],
        &[-3, -5, -6, -7, -2, -5, -6, -7, -6, 7],
        &[-3, -4, -6, -7, -3, -4, -6, -7, -5, 2, 6],
        &[-1, 3, 0, -2, -6, 2, 1, -3, -1, -5, -4, 8],
        &[-2, -3, -4, -6, -3, -1, -4, -5, -4, 2, 3, -3, 9],
        &[-4, -5, -6, -6, -4, -5, -6, -6, -2, -1, 0, -5, 0, 10],
        &[-1, -3, -4, -3, -6, -3, -2, -5, -4, -5, -6, -2, -4, -6, 12],
        &[2, -2, 1, -1, -2, -1, -1, -1, -2, -4, -4, -1, -3, -4, -2, 7],
        &[0, -2, 0, -2, -2, -1, -2, -3, -3, -2, -3, -1, -1, -4, -3, 2, 8],
        &[
            -5, -5, -7, -8, -5, -4, -6, -6, -4, -5, -4, -6, -3, 0, -7, -6, -5, 16,
        ],
        &[
            -4, -4, -4, -6, -5, -3, -5, -6, 3, -3, -2, -4, -3, 4, -6, -3, -3, 3, 11,
        ],
        &[
            -1, -4, -5, -6, -2, -4, -4, -6, -5, 4, 1, -4, 1, -2, -4, -3, 0, -5, -3, 7,
        ],
        &[
            -3, -2, 5, 6, -4, -1, 1, -2, -1, -6, -7, -1, -5, -6, -4, 0, -1, -8, -5, -6, 6,
        ],
        &[
            -2, 0, -1, 1, -7, 5, 6, -4, 0, -6, -5, 1, -3, -6, -2, -1, -2, -5, -4, -4, 0, 6,
        ],
        &[
            -1, -2, -2, -3, -4, -2, -2, -3, -3, -2, -3, -2, -2, -3, -3, -1, -1, -5, -3, -2, -3, -1,
            -2,
        ],
        &[
            -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,
            -8, 1,
        ],
    ])
}

/// PAM30 (NCBI canonical).
fn pam30_table() -> Vec<i32> {
    protein_table(&[
        &[6],
        &[-7, 8],
        &[-4, -6, 8],
        &[-3, -10, 2, 8],
        &[-6, -8, -11, -14, 10],
        &[-4, -2, -3, -2, -14, 8],
        &[-2, -9, -2, 2, -14, 1, 8],
        &[-2, -9, -3, -3, -9, -7, -4, 6],
        &[-7, -2, 0, -4, -7, 1, -5, -9, 9],
        &[-5, -5, -5, -7, -6, -8, -5, -11, -9, 8],
        &[-6, -8, -7, -12, -15, -5, -9, -10, -6, -1, 7],
        &[-7, 0, -1, -4, -14, -3, -4, -7, -6, -6, -8, 7],
        &[-5, -4, -9, -11, -13, -4, -7, -8, -10, -1, 1, -2, 11],
        &[-8, -9, -9, -15, -13, -13, -14, -9, -6, -2, -3, -14, -4, 9],
        &[-2, -4, -6, -8, -8, -3, -5, -6, -4, -8, -7, -6, -8, -10, 8],
        &[0, -3, 0, -4, -3, -5, -4, -2, -6, -7, -8, -4, -5, -6, -2, 6],
        &[-1, -6, -2, -5, -8, -5, -6, -6, -7, -2, -7, -3, -4, -9, -4, 0, 7],
        &[
            -13, -2, -8, -15, -15, -13, -17, -15, -7, -14, -6, -12, -13, -4, -14, -5, -13, 13,
        ],
        &[
            -8, -10, -4, -11, -4, -12, -8, -14, -3, -6, -7, -9, -11, 2, -13, -7, -6, -5, 10,
        ],
        &[
            -2, -8, -8, -8, -6, -7, -6, -5, -6, 2, -2, -9, -1, -8, -6, -6, -3, -15, -7, 7,
        ],
        &[
            -3, -7, 6, 6, -12, -3, 1, -3, -1, -6, -9, -2, -10, -10, -7, -1, -3, -10, -6, -8, 6,
        ],
        &[
            -3, -4, -3, 1, -14, 6, 6, -5, -1, -6, -7, -4, -5, -13, -4, -5, -6, -14, -9, -6, 0, 6,
        ],
        &[
            -3, -6, -3, -5, -9, -5, -5, -5, -5, -6, -6, -5, -5, -8, -5, -3, -4, -11, -7, -5, -5, -5,
            -5,
        ],
        &[
            -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17,
            -17, -17, -17, -17, -17, -17, 1,
        ],
    ])
}

/// PAM70 (NCBI canonical).
fn pam70_table() -> Vec<i32> {
    protein_table(&[
        &[5],
        &[-4, 8],
        &[-2, -3, 6],
        &[-1, -6, 3, 6],
        &[-4, -5, -6, -9, 9],
        &[-2, 0, -1, 0, -9, 7],
        &[-1, -4, 0, 3, -9, 2, 6],
        &[0, -6, -1, -1, -6, -4, -2, 6],
        &[-4, 0, 1, -1, -5, 2, -1, -6, 8],
        &[-2, -3, -3, -5, -4, -5, -4, -6, -6, 7],
        &[-4, -6, -5, -8, -10, -3, -6, -7, -4, 1, 6],
        &[-4, 2, 0, -2, -9, -1, -2, -5, -3, -4, -5, 6],
        &[-3, -2, -5, -7, -9, -2, -4, -6, -6, 1, 2, 0, 10],
        &[-6, -7, -6, -10, -8, -9, -9, -7, -4, 0, -1, -9, -2, 8],
        &[0, -2, -3, -4, -5, -1, -3, -3, -2, -5, -5, -4, -5, -7, 7],
        &[1, -1, 1, -1, -1, -3, -2, 0, -3, -4, -6, -2, -3, -4, 0, 5],
        &[1, -3, 0, -2, -5, -3, -3, -3, -4, -1, -4, -1, -2, -6, -2, 2, 6],
        &[
            -9, 0, -5, -10, -11, -8, -11, -10, -5, -9, -3, -6, -8, -2, -9, -3, -8, 13,
        ],
        &[
            -5, -7, -3, -7, -1, -8, -6, -9, -1, -4, -4, -7, -7, 4, -9, -5, -4, -3, 9,
        ],
        &[
            -1, -5, -5, -5, -4, -4, -4, -3, -4, 3, 0, -6, 0, -5, -3, -3, -1, -10, -5, 6,
        ],
        &[
            -1, -4, 5, 5, -8, -1, 2, -1, 0, -4, -6, -1, -6, -7, -4, 0, -1, -7, -5, -5, 5,
        ],
        &[
            -2, -2, -1, 2, -9, 5, 5, -3, 0, -4, -4, -2, -3, -8, -2, -3, -3, -10, -7, -4, 1, 5,
        ],
        &[
            -2, -3, -2, -3, -6, -2, -3, -3, -3, -3, -4, -3, -3, -5, -3, -1, -2, -7, -5, -2, -2, -2,
            -3,
        ],
        &[
            -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11, -11,
            -11, -11, -11, -11, -11, -11, 1,
        ],
    ])
}

/// PAM250 (NCBI canonical — the classic Dayhoff matrix).
fn pam250_table() -> Vec<i32> {
    protein_table(&[
        &[2],
        &[-2, 6],
        &[0, 0, 2],
        &[0, -1, 2, 4],
        &[-2, -4, -4, -5, 12],
        &[0, 1, 1, 2, -5, 4],
        &[0, -1, 1, 3, -5, 2, 4],
        &[1, -3, 0, 1, -3, -1, 0, 5],
        &[-1, 2, 2, 1, -3, 3, 1, -2, 6],
        &[-1, -2, -2, -2, -2, -2, -2, -3, -2, 5],
        &[-2, -3, -3, -4, -6, -2, -3, -4, -2, 2, 6],
        &[-1, 3, 1, 0, -5, 1, 0, -2, 0, -2, -3, 5],
        &[-1, 0, -2, -3, -5, -1, -2, -3, -2, 2, 4, 0, 6],
        &[-3, -4, -3, -6, -4, -5, -5, -5, -2, 1, 2, -5, 0, 9],
        &[1, 0, 0, -1, -3, 0, -1, 0, 0, -2, -3, -1, -2, -5, 6],
        &[1, 0, 1, 0, 0, -1, 0, 1, -1, -1, -3, 0, -2, -3, 1, 2],
        &[1, -1, 0, 0, -2, -1, 0, 0, -1, 0, -2, 0, -1, -3, 0, 1, 3],
        &[
            -6, 2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4, 0, -6, -2, -5, 17,
        ],
        &[
            -3, -4, -2, -4, 0, -4, -4, -5, 0, -1, -1, -4, -2, 7, -5, -3, -3, 0, 10,
        ],
        &[
            0, -2, -2, -2, -2, -2, -2, -1, -2, 4, 2, -2, 2, -1, -1, -1, 0, -6, -2, 4,
        ],
        &[
            0, -1, 2, 3, -4, 1, 2, 0, 1, -2, -3, 1, -2, -4, -1, 0, 0, -5, -3, -2, 3,
        ],
        &[
            0, 0, 1, 3, -5, 3, 3, 0, 2, -2, -3, 0, -2, -5, 0, 0, -1, -6, -4, -2, 2, 3,
        ],
        &[
            0, -1, 0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, 0, 0, -4, -2, -1, -1, -1, -1,
        ],
        &[
            -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,
            -8, 1,
        ],
    ])
}

/// NUC.4.4 — the EMBOSS `EDNAFULL` 15×15 nucleotide matrix over
/// `ATGCSWRYKMBVHDN`. Exact `A/T/G/C` matches score `+5`; mismatches
/// `-4`; ambiguity codes get the canonical EDNAFULL partial scores.
fn nuc44_table() -> Vec<i32> {
    // Canonical EDNAFULL values, written row-major over NUC_ORDER.
    #[rustfmt::skip]
    let rows: [[i32; 15]; 15] = [
        // A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
        [  5, -4, -4, -4, -4,  1,  1, -4, -4,  1, -4, -1, -1, -1, -2], // A
        [ -4,  5, -4, -4, -4,  1, -4,  1,  1, -4, -1, -4, -1, -1, -2], // T
        [ -4, -4,  5, -4,  1, -4,  1, -4,  1, -4, -1, -1, -4, -1, -2], // G
        [ -4, -4, -4,  5,  1, -4, -4,  1, -4,  1, -1, -1, -1, -4, -2], // C
        [ -4, -4,  1,  1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1], // S
        [  1,  1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1], // W
        [  1, -4,  1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1], // R
        [ -4,  1, -4,  1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1], // Y
        [ -4,  1,  1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1], // K
        [  1, -4, -4,  1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1], // M
        [ -4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1], // B
        [ -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1], // V
        [ -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1], // H
        [ -1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1], // D
        [ -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1], // N
    ];
    rows.iter().flat_map(|r| r.iter().copied()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blosum62_known_values() {
        let m = SubstitutionMatrix::blosum62();
        assert_eq!(m.score(b'A', b'A'), 4);
        assert_eq!(m.score(b'W', b'W'), 11);
        assert_eq!(m.score(b'A', b'R'), -1);
        assert_eq!(m.score(b'R', b'A'), -1, "matrix must be symmetric");
        assert_eq!(m.score(b'C', b'C'), 9);
        assert_eq!(m.score(b'L', b'I'), 2); // conservative substitution
    }

    #[test]
    fn blosum_family_diagonal_positive() {
        for m in [
            SubstitutionMatrix::blosum45(),
            SubstitutionMatrix::blosum62(),
            SubstitutionMatrix::blosum80(),
        ] {
            for &r in b"ARNDCQEGHILKMFPSTWYV" {
                assert!(m.score(r, r) > 0, "{} diag {}", m.name(), r as char);
            }
        }
    }

    #[test]
    fn pam_family_loads_and_symmetric() {
        for m in [
            SubstitutionMatrix::pam30(),
            SubstitutionMatrix::pam70(),
            SubstitutionMatrix::pam250(),
        ] {
            for &a in b"ACDEFG" {
                for &b in b"ACDEFG" {
                    assert_eq!(m.score(a, b), m.score(b, a), "{} not symmetric", m.name());
                }
            }
        }
        // PAM250 W/W is the famous +17.
        assert_eq!(SubstitutionMatrix::pam250().score(b'W', b'W'), 17);
    }

    #[test]
    fn nuc44_values() {
        let m = SubstitutionMatrix::nuc44();
        assert_eq!(m.score(b'A', b'A'), 5);
        assert_eq!(m.score(b'A', b'T'), -4);
        assert_eq!(m.score(b'A', b'G'), -4);
        assert_eq!(m.score(b'N', b'A'), -2);
        // Symmetric.
        assert_eq!(m.score(b'R', b'G'), m.score(b'G', b'R'));
    }

    #[test]
    fn identity_matrix() {
        let m = SubstitutionMatrix::identity(2, -3);
        assert_eq!(m.score(b'K', b'K'), 2);
        assert_eq!(m.score(b'K', b'L'), -3);
    }

    #[test]
    fn unknown_residue_falls_back_to_wildcard() {
        let m = SubstitutionMatrix::blosum62();
        // 'J' is not a standard row; should fall back to X.
        let viaj = m.score(b'J', b'A');
        let viax = m.score(b'X', b'A');
        assert_eq!(viaj, viax);
    }

    #[test]
    fn by_name_dispatch() {
        assert_eq!(SubstitutionMatrix::by_name("blosum62").unwrap().name(), "BLOSUM62");
        assert_eq!(SubstitutionMatrix::by_name("PAM250").unwrap().name(), "PAM250");
        assert_eq!(SubstitutionMatrix::by_name("nuc44").unwrap().name(), "NUC.4.4");
        assert!(SubstitutionMatrix::by_name("nope").is_err());
    }

    #[test]
    fn gap_cost_arithmetic() {
        let g = GapCost::new(11, 1);
        assert_eq!(g.total(0), 0);
        assert_eq!(g.total(1), 12);
        assert_eq!(g.total(3), 14);
        assert!(!g.is_linear());
        let lin = GapCost::linear(2);
        assert!(lin.is_linear());
        assert_eq!(lin.total(4), 8);
        // Negative penalties clamp to zero.
        assert_eq!(GapCost::new(-5, -2), GapCost::new(0, 0));
    }

    #[test]
    fn scoring_scheme_defaults() {
        let s = ScoringScheme::blosum62_default();
        assert_eq!(s.gap.open, 11);
        assert_eq!(s.sub(b'A', b'A'), 4);
        let d = ScoringScheme::dna_default();
        assert_eq!(d.sub(b'A', b'A'), 5);
    }
}

/// Reference-value validation against published textbook results — the
/// canonical Needleman-Wunsch and Smith-Waterman worked examples and the
/// published BLOSUM62 matrix.
#[cfg(test)]
mod validation {
    use super::*;
    use crate::pairwise::global::needleman_wunsch;
    use crate::pairwise::local::smith_waterman;

    /// The canonical Needleman-Wunsch worked example (the one on the
    /// Wikipedia "Needleman-Wunsch algorithm" page): aligning
    /// `GCATGCG` with `GATTACA` under match +1, mismatch -1, gap -1 has
    /// optimal global score 0.
    #[test]
    fn needleman_wunsch_textbook_example_scores_zero() {
        let scheme = ScoringScheme::new(
            SubstitutionMatrix::dna_simple(1, -1),
            GapCost::new(0, 1),
        );
        let al = needleman_wunsch(b"GCATGCG", b"GATTACA", &scheme).unwrap();
        assert_eq!(al.score, 0, "textbook NW optimal score is 0");
        // Both rows are the same length.
        assert_eq!(al.row1.len(), al.row2.len());
        // 4 matched identities on the optimal path.
        let ident = al
            .row1
            .iter()
            .zip(&al.row2)
            .filter(|(&a, &b)| a == b && a != b'-')
            .count();
        assert_eq!(ident, 4);
    }

    /// The canonical Smith-Waterman worked example (Wikipedia
    /// "Smith-Waterman algorithm"): the best local alignment of
    /// `TGTTACGG` and `GGTTGACTA` under match +3, mismatch -3, gap -2
    /// has score 13, with the aligned core `GTT-AC` / `GTTGAC`.
    #[test]
    fn smith_waterman_textbook_example_scores_13() {
        let scheme = ScoringScheme::new(
            SubstitutionMatrix::dna_simple(3, -3),
            GapCost::new(0, 2),
        );
        let al = smith_waterman(b"TGTTACGG", b"GGTTGACTA", &scheme).unwrap();
        assert_eq!(al.score, 13, "textbook SW optimal local score is 13");
        assert_eq!(al.row1_str(), "GTT-AC");
        assert_eq!(al.row2_str(), "GTTGAC");
    }

    /// Spot-check the embedded BLOSUM62 against the published matrix —
    /// diagonal self-scores and a few well-known off-diagonal entries.
    #[test]
    fn blosum62_matches_the_published_matrix() {
        let m = SubstitutionMatrix::blosum62();
        // Published BLOSUM62 diagonal self-scores.
        let diag: &[(u8, i32)] = &[
            (b'A', 4),
            (b'R', 5),
            (b'N', 6),
            (b'D', 6),
            (b'C', 9),
            (b'Q', 5),
            (b'E', 5),
            (b'G', 6),
            (b'H', 8),
            (b'I', 4),
            (b'L', 4),
            (b'K', 5),
            (b'M', 5),
            (b'F', 6),
            (b'P', 7),
            (b'S', 4),
            (b'T', 5),
            (b'W', 11),
            (b'Y', 7),
            (b'V', 4),
        ];
        for &(r, expected) in diag {
            assert_eq!(
                m.score(r, r),
                expected,
                "BLOSUM62 diagonal for {}",
                r as char
            );
        }
        // Well-known off-diagonal entries from the published matrix.
        assert_eq!(m.score(b'W', b'C'), -2);
        assert_eq!(m.score(b'L', b'I'), 2); // conservative hydrophobic swap
        assert_eq!(m.score(b'F', b'Y'), 3); // aromatic pair
        assert_eq!(m.score(b'K', b'R'), 2); // basic pair
        assert_eq!(m.score(b'D', b'E'), 2); // acidic pair
        // The matrix is symmetric.
        assert_eq!(m.score(b'A', b'R'), m.score(b'R', b'A'));
    }
}
