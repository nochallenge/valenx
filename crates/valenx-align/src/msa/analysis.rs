//! MSA analysis — conservation, Shannon entropy, consensus.
//!
//! Once an alignment exists the natural questions are *which columns
//! are conserved* and *what is the consensus sequence*. This module
//! answers both:
//!
//! - [`column_entropy`] / [`entropy_profile`] — the Shannon entropy of
//!   each column's residue distribution (low entropy = conserved).
//! - [`column_conservation`] / [`conservation_profile`] — a `[0, 1]`
//!   conservation score (`1 −` normalised entropy).
//! - [`consensus`] — the per-column majority residue, with a
//!   configurable gap / ambiguity threshold.
//! - [`identity_columns`] — count of fully conserved (single-residue)
//!   columns.

use super::progressive::Msa;

/// Shannon entropy (in bits) of one alignment column.
///
/// `column` is the slice of residues down one MSA column. Gaps are
/// ignored (they are not a residue state). An all-gap column has
/// entropy `0`. A column of one residue has entropy `0`; a uniform
/// mix of `n` residues has entropy `log2(n)`.
pub fn column_entropy(column: &[u8]) -> f64 {
    use std::collections::HashMap;
    let mut counts: HashMap<u8, usize> = HashMap::new();
    let mut total = 0usize;
    for &c in column {
        if c == b'-' {
            continue;
        }
        *counts.entry(c.to_ascii_uppercase()).or_insert(0) += 1;
        total += 1;
    }
    if total == 0 {
        return 0.0;
    }
    let mut h = 0.0;
    for &cnt in counts.values() {
        let p = cnt as f64 / total as f64;
        h -= p * p.log2();
    }
    h
}

/// The Shannon entropy of every column of `msa`, left to right.
pub fn entropy_profile(msa: &Msa) -> Vec<f64> {
    (0..msa.width())
        .map(|c| column_entropy(&column_bytes(msa, c)))
        .collect()
}

/// Conservation score of one column in `[0, 1]`.
///
/// `1.0` means perfectly conserved (one residue); `0.0` means maximal
/// diversity. Computed as `1 − H / Hmax` where `Hmax = log2(alphabet)`
/// — `log2(4)` for nucleotides, `log2(20)` for protein. Pass the
/// expected alphabet size as `alphabet_size`.
pub fn column_conservation(column: &[u8], alphabet_size: usize) -> f64 {
    let h = column_entropy(column);
    let hmax = (alphabet_size.max(2) as f64).log2();
    (1.0 - h / hmax).clamp(0.0, 1.0)
}

/// The conservation score of every column of `msa` (see
/// [`column_conservation`]).
pub fn conservation_profile(msa: &Msa, alphabet_size: usize) -> Vec<f64> {
    (0..msa.width())
        .map(|c| column_conservation(&column_bytes(msa, c), alphabet_size))
        .collect()
}

/// Options controlling [`consensus`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ConsensusOptions {
    /// A column whose majority residue has frequency below this
    /// threshold (fraction of non-gap residues) emits the ambiguity
    /// character instead. `0.5` is the common majority rule.
    pub majority_threshold: f64,
    /// A column with a gap fraction at or above this emits a gap in
    /// the consensus.
    pub gap_threshold: f64,
    /// The character emitted for an under-threshold column.
    pub ambiguity_char: u8,
}

impl Default for ConsensusOptions {
    fn default() -> Self {
        ConsensusOptions {
            majority_threshold: 0.5,
            gap_threshold: 0.5,
            ambiguity_char: b'X',
        }
    }
}

/// The consensus sequence of an MSA under `opts`.
///
/// Each column emits its majority residue, or the ambiguity character
/// if no residue clears `majority_threshold`, or a gap if the column
/// is gap-heavy past `gap_threshold`.
pub fn consensus(msa: &Msa, opts: ConsensusOptions) -> Vec<u8> {
    let mut out = Vec::with_capacity(msa.width());
    for c in 0..msa.width() {
        let col = column_bytes(msa, c);
        let depth = col.len();
        let gap_count = col.iter().filter(|&&b| b == b'-').count();
        if depth > 0 && gap_count as f64 / depth as f64 >= opts.gap_threshold {
            out.push(b'-');
            continue;
        }
        // Tally residues into a fixed 256-slot array keyed by byte value
        // so the majority pick is deterministic: ties resolve to the
        // lowest byte (first by alphabet), matching
        // `ProfileColumn::consensus`. A `HashMap` + `max_by_key` made
        // the winner depend on iteration order — nondeterministic.
        let mut counts = [0usize; 256];
        let mut residues = 0usize;
        for &b in &col {
            if b == b'-' {
                continue;
            }
            counts[b.to_ascii_uppercase() as usize] += 1;
            residues += 1;
        }
        if residues == 0 {
            out.push(b'-');
            continue;
        }
        // Strict `>` keeps the first (lowest-byte) residue on a tie.
        let mut best_byte = 0usize;
        let mut best_count = 0usize;
        for (byte, &n) in counts.iter().enumerate() {
            if n > best_count {
                best_count = n;
                best_byte = byte;
            }
        }
        if best_count as f64 / residues as f64 >= opts.majority_threshold {
            out.push(best_byte as u8);
        } else {
            out.push(opts.ambiguity_char);
        }
    }
    out
}

/// Number of fully conserved columns — columns containing exactly one
/// distinct residue (gaps ignored).
pub fn identity_columns(msa: &Msa) -> usize {
    (0..msa.width())
        .filter(|&c| {
            let col = column_bytes(msa, c);
            let mut seen: Option<u8> = None;
            for &b in &col {
                if b == b'-' {
                    continue;
                }
                let u = b.to_ascii_uppercase();
                match seen {
                    None => seen = Some(u),
                    Some(s) if s != u => return false,
                    _ => {}
                }
            }
            seen.is_some()
        })
        .count()
}

/// Mean column conservation of an MSA in `[0, 1]` — a single-number
/// quality summary.
pub fn mean_conservation(msa: &Msa, alphabet_size: usize) -> f64 {
    let prof = conservation_profile(msa, alphabet_size);
    if prof.is_empty() {
        0.0
    } else {
        prof.iter().sum::<f64>() / prof.len() as f64
    }
}

/// Collects the bytes of column `c` down every row.
fn column_bytes(msa: &Msa, c: usize) -> Vec<u8> {
    msa.rows.iter().map(|r| r[c]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn msa(rows: &[&[u8]]) -> Msa {
        Msa::new(rows.iter().map(|r| r.to_vec()).collect()).unwrap()
    }

    #[test]
    fn entropy_of_conserved_column_is_zero() {
        assert!(column_entropy(b"AAAA").abs() < 1e-9);
        assert!(column_entropy(b"----").abs() < 1e-9);
    }

    #[test]
    fn entropy_of_uniform_column() {
        // 4 distinct residues, equal frequency -> log2(4) = 2 bits.
        assert!((column_entropy(b"ACGT") - 2.0).abs() < 1e-9);
        // 2 distinct, equal -> 1 bit.
        assert!((column_entropy(b"AATT") - 1.0).abs() < 1e-9);
    }

    #[test]
    fn entropy_ignores_gaps() {
        // Gaps dropped: "A-A-" is effectively "AA" -> entropy 0.
        assert!(column_entropy(b"A-A-").abs() < 1e-9);
    }

    #[test]
    fn conservation_endpoints() {
        // Fully conserved -> 1.0.
        assert!((column_conservation(b"AAAA", 4) - 1.0).abs() < 1e-9);
        // Maximally diverse over the alphabet -> 0.0.
        assert!(column_conservation(b"ACGT", 4).abs() < 1e-9);
    }

    #[test]
    fn entropy_and_conservation_profiles() {
        let m = msa(&[b"AAAA", b"AACA", b"AAGA"]);
        let ent = entropy_profile(&m);
        assert_eq!(ent.len(), 4);
        // Columns 0,1,3 conserved -> entropy 0; column 2 diverse.
        assert!(ent[0].abs() < 1e-9);
        assert!(ent[2] > 0.0);
        let cons = conservation_profile(&m, 4);
        assert!(cons[0] > cons[2]);
    }

    #[test]
    fn consensus_majority_rule() {
        // Column 2 is A,A,C -> 2/3 A clears the 0.5 threshold.
        let m = msa(&[b"GGAGG", b"GGAGG", b"GGCGG"]);
        let c = consensus(&m, ConsensusOptions::default());
        assert_eq!(c, b"GGAGG");
    }

    #[test]
    fn consensus_ambiguity_when_split() {
        // Column 0 splits A/C/G — none clears 0.5 -> ambiguity 'X'.
        let m = msa(&[b"AT", b"CT", b"GT"]);
        let c = consensus(&m, ConsensusOptions::default());
        assert_eq!(c[0], b'X');
        assert_eq!(c[1], b'T');
    }

    #[test]
    fn consensus_gap_column() {
        // Column 1 is gap-heavy (>= 0.5) -> consensus gap.
        let m = msa(&[b"A-C", b"A-C", b"AGC"]);
        let c = consensus(&m, ConsensusOptions::default());
        assert_eq!(c, b"A-C");
    }

    #[test]
    fn consensus_tie_is_deterministic_first_by_alphabet() {
        // A 4-way tie A/C/G/T (each 1/4). With a 0.0 majority threshold
        // every residue clears it, so the column emits its "winner".
        // The previous HashMap `max_by_key` made that winner depend on
        // iteration order — nondeterministic across runs. The rule must
        // be the same as ProfileColumn::consensus: first by byte value,
        // i.e. 'A'.
        let opts = ConsensusOptions {
            majority_threshold: 0.0,
            gap_threshold: 1.0,
            ambiguity_char: b'X',
        };
        let m = msa(&[b"A", b"C", b"G", b"T"]);
        // Run many times; the result must never vary.
        let first = consensus(&m, opts);
        for _ in 0..64 {
            assert_eq!(consensus(&m, opts), first, "consensus must be deterministic");
        }
        assert_eq!(first, b"A", "tie must resolve to the alphabetically-first residue");

        // A two-way tie that clears the default 0.5 threshold: C vs A,
        // depth 2 -> each 0.5. Must pick 'A' (lower byte), every run.
        let m2 = msa(&[b"C", b"A"]);
        let c2 = consensus(&m2, ConsensusOptions::default());
        assert_eq!(c2, b"A", "0.5/0.5 tie resolves to 'A', deterministically");
    }

    #[test]
    fn identity_column_count() {
        let m = msa(&[b"AAAA", b"AACA", b"AAGA"]);
        // Columns 0,1,3 fully conserved; column 2 not.
        assert_eq!(identity_columns(&m), 3);
    }

    #[test]
    fn mean_conservation_summary() {
        let perfect = msa(&[b"ACGT", b"ACGT"]);
        assert!((mean_conservation(&perfect, 4) - 1.0).abs() < 1e-9);
    }
}
