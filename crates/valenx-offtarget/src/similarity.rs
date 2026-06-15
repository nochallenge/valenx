//! The two sequence-similarity measures used for off-target screening.

use std::collections::HashSet;

use crate::aa::first_invalid;
use crate::error::OffTargetError;

/// The best fractional residue identity over any ungapped overlap of `a` and
/// `b`, in `[0, 1]`.
///
/// The shorter sequence is slid across the longer one; at each offset the
/// fraction of matching positions (over the shorter length) is taken, and the
/// maximum is returned. Comparison is case-insensitive. Identical sequences
/// score 1.0; sequences sharing no residue in any frame score 0.0.
///
/// Both sequences must be non-empty and contain only standard amino acids.
pub fn best_ungapped_identity(a: &str, b: &str) -> Result<f64, OffTargetError> {
    validate("a", a)?;
    validate("b", b)?;
    let ua = a.to_ascii_uppercase();
    let ub = b.to_ascii_uppercase();
    let (short, long) = if ua.len() <= ub.len() {
        (ua.as_bytes(), ub.as_bytes())
    } else {
        (ub.as_bytes(), ua.as_bytes())
    };
    let mut best = 0.0_f64;
    for start in 0..=(long.len() - short.len()) {
        let matches = short
            .iter()
            .zip(&long[start..start + short.len()])
            .filter(|(x, y)| x == y)
            .count();
        let identity = matches as f64 / short.len() as f64;
        if identity > best {
            best = identity;
        }
    }
    Ok(best)
}

/// The Jaccard overlap of the length-`k` k-mer (k-peptide) sets of `a` and `b`,
/// in `[0, 1]`.
///
/// `|A ∩ B| / |A ∪ B|`, comparing the *sets* of distinct k-mers. A sequence
/// shorter than `k` contributes no k-mers; if neither sequence has any k-mer the
/// union is empty and the result is `0.0`. Comparison is case-insensitive.
///
/// `k` must be at least 1; both sequences must be non-empty and all-standard.
pub fn kmer_jaccard(a: &str, b: &str, k: usize) -> Result<f64, OffTargetError> {
    if k == 0 {
        return Err(OffTargetError::ZeroK);
    }
    validate("a", a)?;
    validate("b", b)?;
    let ua = a.to_ascii_uppercase();
    let ub = b.to_ascii_uppercase();
    let sa = kmer_set(&ua, k);
    let sb = kmer_set(&ub, k);
    let inter = sa.intersection(&sb).count();
    let union = sa.union(&sb).count();
    if union == 0 {
        return Ok(0.0);
    }
    Ok(inter as f64 / union as f64)
}

fn kmer_set(seq: &str, k: usize) -> HashSet<&[u8]> {
    let bytes = seq.as_bytes();
    let mut set = HashSet::new();
    if bytes.len() >= k {
        for i in 0..=(bytes.len() - k) {
            set.insert(&bytes[i..i + k]);
        }
    }
    set
}

fn validate(which: &str, seq: &str) -> Result<(), OffTargetError> {
    if seq.is_empty() {
        return Err(OffTargetError::Empty { what: "sequence" });
    }
    if let Some((pos, residue)) = first_invalid(seq) {
        return Err(OffTargetError::InvalidResidue {
            which: which.to_string(),
            pos,
            residue,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_of_identical_is_one() {
        assert!((best_ungapped_identity("AAAA", "AAAA").unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn identity_of_disjoint_is_zero() {
        assert!(best_ungapped_identity("AAAA", "WWWW").unwrap().abs() < 1e-12);
    }

    #[test]
    fn identity_finds_best_frame() {
        // "AC" matches the "AC" frame inside "AACC" exactly -> 1.0.
        assert!((best_ungapped_identity("AC", "AACC").unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn identity_equal_length_partial() {
        // ACD vs AGD: A and D match, G != C -> 2/3.
        let id = best_ungapped_identity("ACD", "AGD").unwrap();
        assert!((id - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn identity_is_case_insensitive() {
        assert!((best_ungapped_identity("acd", "ACD").unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn jaccard_identical_is_one() {
        assert!((kmer_jaccard("AAAA", "AAAA", 2).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn jaccard_known_overlap() {
        // ACDE -> {AC,CD,DE}; CDEF -> {CD,DE,EF}; inter=2, union=4 -> 0.5.
        let j = kmer_jaccard("ACDE", "CDEF", 2).unwrap();
        assert!((j - 0.5).abs() < 1e-12);
    }

    #[test]
    fn jaccard_disjoint_is_zero() {
        assert!(kmer_jaccard("AAAA", "CCCC", 2).unwrap().abs() < 1e-12);
    }

    #[test]
    fn jaccard_k_larger_than_seq_is_zero() {
        // Neither 2-letter sequence has a 5-mer -> empty union -> 0.0.
        assert!(kmer_jaccard("AC", "AC", 5).unwrap().abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(kmer_jaccard("AC", "AC", 0).unwrap_err().code(), "zero_k");
        assert_eq!(
            best_ungapped_identity("", "AC").unwrap_err().code(),
            "empty"
        );
        assert_eq!(
            best_ungapped_identity("AC", "AX").unwrap_err().code(),
            "invalid_residue"
        );
    }
}
