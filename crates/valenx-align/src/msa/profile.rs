//! Sequence profiles and profile-profile alignment.
//!
//! A [`Profile`] is the column-wise residue-frequency model of a set
//! of aligned sequences — the data structure progressive MSA merges.
//! Each column holds the count of every residue plus a gap count; the
//! profile-vs-profile score of two columns is the
//! frequency-weighted average of the underlying substitution scores.
//!
//! [`align_profiles`] runs a Needleman-Wunsch global alignment over
//! two profiles (treating each as a single "super-sequence" of
//! columns) and returns the column-level edit script needed to merge
//! them; [`align_profile_sequence`] is the common special case of
//! adding one plain sequence to a growing profile.

use crate::error::{AlignError, Result};
use crate::limits::{check_dp_size_with, MAX_DP_CELLS};
use crate::matrix::ScoringScheme;

/// The residue alphabet a [`Profile`] tallies — the 20 amino acids
/// plus `X`. DNA columns use the same machinery (the unused rows just
/// stay zero), keeping one code path.
pub const PROFILE_ALPHABET: &[u8; 21] = b"ACDEFGHIKLMNPQRSTVWYX";

/// One column of a [`Profile`]: a residue-count histogram.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct ProfileColumn {
    /// `counts[i]` = number of sequences with `PROFILE_ALPHABET[i]` in
    /// this column.
    pub counts: [u32; 21],
    /// Number of sequences with a gap in this column.
    pub gaps: u32,
}

impl ProfileColumn {
    /// Total number of sequences contributing to this column
    /// (residues + gaps).
    pub fn depth(&self) -> u32 {
        self.counts.iter().sum::<u32>() + self.gaps
    }

    /// The most frequent residue (ties broken by alphabet order), or
    /// `-` if the column is all gaps.
    pub fn consensus(&self) -> u8 {
        let (mut best_i, mut best_c) = (usize::MAX, 0u32);
        for (i, &c) in self.counts.iter().enumerate() {
            if c > best_c {
                best_c = c;
                best_i = i;
            }
        }
        if best_i == usize::MAX || best_c == 0 {
            b'-'
        } else {
            PROFILE_ALPHABET[best_i]
        }
    }

    /// Adds one residue (or `-`) to the column's tally.
    fn add(&mut self, residue: u8) {
        if residue == b'-' {
            self.gaps += 1;
            return;
        }
        match index_of(residue) {
            Some(i) => self.counts[i] += 1,
            None => self.counts[20] += 1, // unknown -> X bucket
        }
    }
}

/// A column-wise residue-frequency profile of aligned sequences.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Profile {
    /// The aligned rows the profile was built from (kept so a merge can
    /// re-emit the full MSA). Each row is an equal-length gapped byte
    /// string.
    pub rows: Vec<Vec<u8>>,
    /// One [`ProfileColumn`] per alignment column.
    pub columns: Vec<ProfileColumn>,
}

impl Profile {
    /// Builds a profile from a set of *already aligned* equal-length
    /// rows. Returns [`AlignError::Dimension`] on a length mismatch and
    /// [`AlignError::Invalid`] for an empty set.
    pub fn from_alignment(rows: &[Vec<u8>]) -> Result<Self> {
        if rows.is_empty() {
            return Err(AlignError::invalid("rows", "profile needs >= 1 sequence"));
        }
        let width = rows[0].len();
        for r in rows {
            if r.len() != width {
                return Err(AlignError::dimension(format!(
                    "profile rows differ: {} vs {width}",
                    r.len()
                )));
            }
        }
        let mut columns = vec![ProfileColumn::default(); width];
        for r in rows {
            for (c, &residue) in r.iter().enumerate() {
                columns[c].add(residue);
            }
        }
        Ok(Profile {
            rows: rows.to_vec(),
            columns,
        })
    }

    /// Builds a single-sequence profile from one ungapped sequence.
    pub fn from_sequence(seq: &[u8]) -> Result<Self> {
        Self::from_alignment(&[seq.to_vec()])
    }

    /// Number of alignment columns.
    pub fn width(&self) -> usize {
        self.columns.len()
    }

    /// Number of sequences in the profile.
    pub fn depth(&self) -> usize {
        self.rows.len()
    }

    /// `true` if the profile has no columns.
    pub fn is_empty(&self) -> bool {
        self.columns.is_empty()
    }

    /// The consensus sequence — the per-column most-frequent residue.
    pub fn consensus(&self) -> Vec<u8> {
        self.columns.iter().map(ProfileColumn::consensus).collect()
    }

    /// The frequency-weighted substitution score of profile column `i`
    /// against profile column `j` of `other`, under `scheme`. Gaps in
    /// either column contribute zero (gap-vs-gap is free; the DP
    /// handles real gap opening).
    fn column_score(&self, i: usize, other: &Profile, j: usize, scheme: &ScoringScheme) -> i32 {
        let ca = &self.columns[i];
        let cb = &other.columns[j];
        let mut total = 0i64;
        let mut weight = 0i64;
        for (ai, &na) in ca.counts.iter().enumerate() {
            if na == 0 {
                continue;
            }
            for (bi, &nb) in cb.counts.iter().enumerate() {
                if nb == 0 {
                    continue;
                }
                let s = scheme.sub(PROFILE_ALPHABET[ai], PROFILE_ALPHABET[bi]);
                total += s as i64 * na as i64 * nb as i64;
                weight += na as i64 * nb as i64;
            }
        }
        if weight == 0 {
            0
        } else {
            (total / weight) as i32
        }
    }
}

/// The merge plan produced by [`align_profiles`]: the two input
/// profiles re-expanded so every row has the same final width.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ProfileAlignment {
    /// The merged alignment rows: every row of profile A followed by
    /// every row of profile B, all padded to the common width.
    pub rows: Vec<Vec<u8>>,
    /// The global alignment score over the profile columns.
    pub score: i32,
}

/// Aligns two profiles with **affine-gap** global (Gotoh) DP over their
/// columns and returns the merged alignment.
///
/// A "match" places column `i` of A opposite column `j` of B; a "gap"
/// inserts an all-gap column into one side. The gap penalty is the
/// scheme's **affine** [`crate::matrix::GapCost`]: opening a run of
/// gap columns costs `open + extend`, each further gap column costs
/// `extend`.
///
/// Honouring `open` is essential, not cosmetic: with a linear penalty
/// of just `extend` per column, a single cheap gap column can beat a
/// mismatch column (e.g. NUC.4.4 mismatch −4 vs `extend` 1), so the DP
/// would gap two *substituted* profiles apart instead of aligning the
/// mismatches — producing a spuriously wide, all-gaps merge. The
/// affine open penalty restores the intended behaviour: a substitution
/// aligns as a mismatch column, a true indel opens a gap.
pub fn align_profiles(
    a: &Profile,
    b: &Profile,
    scheme: &ScoringScheme,
) -> Result<ProfileAlignment> {
    align_profiles_capped(a, b, scheme, MAX_DP_CELLS)
}

/// [`align_profiles`] with an explicit cell cap (test seam). Over the
/// cap it returns
/// [`AlignError::TooLarge`](crate::error::AlignError::TooLarge):
/// profile-profile alignment holds six matrices and has no linear-space
/// variant, so an oversized merge is rejected.
fn align_profiles_capped(
    a: &Profile,
    b: &Profile,
    scheme: &ScoringScheme,
    max_cells: usize,
) -> Result<ProfileAlignment> {
    let n = a.width();
    let m = b.width();
    let open = scheme.gap.open;
    let ext = scheme.gap.extend;
    let w = m + 1;
    // A large but arithmetic-safe "negative infinity" — small enough
    // that adding gap penalties to it never overflows `i32`.
    const NEG_INF: i32 = i32::MIN / 4;

    // Six matrices of (n+1)·(m+1) — bound the allocation.
    check_dp_size_with(n + 1, m + 1, max_cells)?;

    // Three Gotoh matrices over the (n+1) x (m+1) grid:
    //   mm — best score ending with column i of A opposite column j of B
    //   gb — best score ending with a gap column in B (A advanced)
    //   ga — best score ending with a gap column in A (B advanced)
    let mut mm = vec![NEG_INF; (n + 1) * w];
    let mut gb = vec![NEG_INF; (n + 1) * w];
    let mut ga = vec![NEG_INF; (n + 1) * w];
    // Traceback: which matrix the cell's best came from. 0 diag(mm),
    // 1 up(gb), 2 left(ga); a separate "previous state" per matrix.
    let mut tb_m = vec![0u8; (n + 1) * w];
    let mut tb_gb = vec![0u8; (n + 1) * w];
    let mut tb_ga = vec![0u8; (n + 1) * w];

    mm[0] = 0;
    for j in 1..=m {
        // A leading run of gap-in-A columns. Computed in i64 and clamped
        // so a very wide profile cannot overflow i32 in the init row.
        ga[j] = crate::pairwise::affine_init(open, ext, j, NEG_INF);
        tb_ga[j] = 2; // extend an existing gap-in-A run
    }
    for i in 1..=n {
        gb[i * w] = crate::pairwise::affine_init(open, ext, i, NEG_INF);
        tb_gb[i * w] = 1; // extend an existing gap-in-B run
    }
    for i in 1..=n {
        for j in 1..=m {
            let idx = i * w + j;
            // --- match/substitution column (came from any state) ---
            let col = a.column_score(i - 1, b, j - 1, scheme);
            let dprev = (i - 1) * w + (j - 1);
            let (m_from, m_src) = max3(mm[dprev], gb[dprev], ga[dprev]);
            mm[idx] = m_from.saturating_add(col);
            tb_m[idx] = m_src;
            // --- gap column in B: A advances (open a new run from mm,
            //     or extend the existing gb run) ---
            let uprev = (i - 1) * w + j;
            let open_gb = mm[uprev].saturating_sub(open + ext);
            let ext_gb = gb[uprev].saturating_sub(ext);
            if ext_gb >= open_gb {
                gb[idx] = ext_gb;
                tb_gb[idx] = 1; // extend
            } else {
                gb[idx] = open_gb;
                tb_gb[idx] = 0; // opened from mm
            }
            // --- gap column in A: B advances ---
            let lprev = i * w + (j - 1);
            let open_ga = mm[lprev].saturating_sub(open + ext);
            let ext_ga = ga[lprev].saturating_sub(ext);
            if ext_ga >= open_ga {
                ga[idx] = ext_ga;
                tb_ga[idx] = 2; // extend
            } else {
                ga[idx] = open_ga;
                tb_ga[idx] = 0; // opened from mm
            }
        }
    }

    // Traceback from the best of the three end states.
    let end = n * w + m;
    let (_, mut state) = max3(mm[end], gb[end], ga[end]);
    let mut ops: Vec<u8> = Vec::new();
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        match state {
            0 => {
                // a match column — move diagonally; the next state is
                // whatever fed this match cell.
                ops.push(0);
                let src = tb_m[i * w + j];
                i -= 1;
                j -= 1;
                state = src;
            }
            1 => {
                // a gap-in-B column — A advances.
                ops.push(1);
                let src = tb_gb[i * w + j];
                i -= 1;
                state = src;
            }
            _ => {
                // a gap-in-A column — B advances.
                ops.push(2);
                let src = tb_ga[i * w + j];
                j -= 1;
                state = src;
            }
        }
    }
    ops.reverse();

    // Re-expand both profiles' rows following the op script.
    let na_rows = a.rows.len();
    let nb_rows = b.rows.len();
    let mut out: Vec<Vec<u8>> = vec![Vec::new(); na_rows + nb_rows];
    let (mut ca, mut cb) = (0usize, 0usize);
    for op in ops {
        match op {
            0 => {
                for (r, row) in a.rows.iter().enumerate() {
                    out[r].push(row[ca]);
                }
                for (r, row) in b.rows.iter().enumerate() {
                    out[na_rows + r].push(row[cb]);
                }
                ca += 1;
                cb += 1;
            }
            1 => {
                // gap in B: A advances, B gets a gap column.
                for (r, row) in a.rows.iter().enumerate() {
                    out[r].push(row[ca]);
                }
                for r in 0..nb_rows {
                    out[na_rows + r].push(b'-');
                }
                ca += 1;
            }
            _ => {
                // gap in A.
                for row in out.iter_mut().take(na_rows) {
                    row.push(b'-');
                }
                for (r, row) in b.rows.iter().enumerate() {
                    out[na_rows + r].push(row[cb]);
                }
                cb += 1;
            }
        }
    }

    let (final_score, _) = max3(mm[end], gb[end], ga[end]);
    Ok(ProfileAlignment {
        rows: out,
        score: final_score,
    })
}

/// Returns the maximum of three scores and a tag for which one won:
/// `0` for the first (`mm`), `1` for the second (`gb`), `2` for the
/// third (`ga`). Ties prefer `mm`, then `gb` — biasing the traceback
/// toward a match column over a gap.
fn max3(m: i32, gb: i32, ga: i32) -> (i32, u8) {
    let mut best = m;
    let mut tag = 0u8;
    if gb > best {
        best = gb;
        tag = 1;
    }
    if ga > best {
        best = ga;
        tag = 2;
    }
    (best, tag)
}

/// Aligns a single plain (ungapped) sequence to a profile — the common
/// "add one sequence to the growing MSA" step. A thin wrapper over
/// [`align_profiles`].
pub fn align_profile_sequence(
    profile: &Profile,
    seq: &[u8],
    scheme: &ScoringScheme,
) -> Result<ProfileAlignment> {
    let seq_profile = Profile::from_sequence(seq)?;
    align_profiles(profile, &seq_profile, scheme)
}

/// Alphabet index of a residue, or `None` if outside the 20 AAs + X.
fn index_of(residue: u8) -> Option<usize> {
    PROFILE_ALPHABET
        .iter()
        .position(|&c| c == residue.to_ascii_uppercase())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::{GapCost, ScoringScheme, SubstitutionMatrix};

    fn scheme() -> ScoringScheme {
        ScoringScheme::new(SubstitutionMatrix::blosum62(), GapCost::new(0, 4))
    }

    #[test]
    fn profile_from_alignment() {
        let rows = vec![b"MKV".to_vec(), b"MRV".to_vec(), b"MK-".to_vec()];
        let p = Profile::from_alignment(&rows).unwrap();
        assert_eq!(p.width(), 3);
        assert_eq!(p.depth(), 3);
        // Column 0 is all M.
        assert_eq!(p.columns[0].consensus(), b'M');
        // Column 2 has one gap.
        assert_eq!(p.columns[2].gaps, 1);
    }

    #[test]
    fn profile_rejects_ragged_rows() {
        let rows = vec![b"MKV".to_vec(), b"MK".to_vec()];
        assert!(Profile::from_alignment(&rows).is_err());
        assert!(Profile::from_alignment(&[]).is_err());
    }

    #[test]
    fn consensus_sequence() {
        let rows = vec![b"AAAA".to_vec(), b"AACA".to_vec(), b"AAGA".to_vec()];
        let p = Profile::from_alignment(&rows).unwrap();
        // Column 2 splits A/C/G evenly; first-by-alphabet wins (A).
        assert_eq!(p.consensus(), b"AAAA");
    }

    #[test]
    fn align_two_single_sequence_profiles() {
        // Profile-profile align of two singletons == pairwise NW.
        let a = Profile::from_sequence(b"MKVLA").unwrap();
        let b = Profile::from_sequence(b"MKVLA").unwrap();
        let merged = align_profiles(&a, &b, &scheme()).unwrap();
        assert_eq!(merged.rows.len(), 2);
        assert_eq!(merged.rows[0], b"MKVLA");
        assert_eq!(merged.rows[1], b"MKVLA");
    }

    #[test]
    fn align_profile_to_sequence_inserts_gap() {
        // Profile of two identical 5-mers; add a sequence missing one
        // residue -> a gap column appears.
        let p = Profile::from_alignment(&[b"MKVLA".to_vec(), b"MKVLA".to_vec()]).unwrap();
        let merged = align_profile_sequence(&p, b"MKLA", &scheme()).unwrap();
        assert_eq!(merged.rows.len(), 3);
        // All rows equal length.
        let w = merged.rows[0].len();
        assert!(merged.rows.iter().all(|r| r.len() == w));
        // The added sequence row carries a gap.
        assert!(merged.rows[2].contains(&b'-'));
    }

    #[test]
    fn merge_preserves_row_count() {
        let a = Profile::from_alignment(&[b"AAA".to_vec(), b"AAA".to_vec()]).unwrap();
        let b = Profile::from_alignment(&[b"AAA".to_vec()]).unwrap();
        let merged = align_profiles(&a, &b, &scheme()).unwrap();
        assert_eq!(merged.rows.len(), 3);
    }

    #[test]
    fn align_profiles_over_cap_errors() {
        use crate::error::AlignError;
        let a = Profile::from_sequence(b"MKVLAAGG").unwrap();
        let b = Profile::from_sequence(b"MKVLAAGG").unwrap();
        // 9*9 = 81 cells; a cap of 8 rejects without the six-matrix
        // allocation.
        let err = align_profiles_capped(&a, &b, &scheme(), 8).unwrap_err();
        assert!(matches!(err, AlignError::TooLarge { .. }), "got {err:?}");
        // Generous cap computes normally.
        assert!(align_profiles_capped(&a, &b, &scheme(), usize::MAX).is_ok());
    }
}
