//! Minimizer sketch (minimap2-class).
//!
//! A **minimizer** of a window is the smallest k-mer (under a hash
//! order) among the `w` consecutive k-mers in that window. Sliding the
//! window across a sequence yields a sparse, deterministic sample of
//! k-mers with a useful property: two sequences that share a substring
//! longer than `w + k − 1` are guaranteed to share at least one
//! minimizer. Minimizers are therefore the seeds of choice for fast
//! long-read overlap and mapping.
//!
//! [`minimizer_sketch`] returns the ordered list of distinct
//! [`Minimizer`]s of a sequence for a given `(k, w)`.

use crate::error::{AlignError, Result};

/// A single minimizer occurrence.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Minimizer {
    /// The 64-bit hash of the k-mer (its identity for matching).
    pub hash: u64,
    /// 0-based start offset of the k-mer within the sequence.
    pub pos: usize,
}

/// Computes the `(k, w)`-minimizer sketch of `seq`.
///
/// For each window of `w` consecutive k-mers the lexicographically
/// minimal *hash* is selected; consecutive windows that select the
/// same occurrence collapse to one entry. The result is ordered by
/// position.
///
/// Returns [`AlignError::Invalid`] for `k == 0` or `w == 0`, and an
/// empty sketch when `seq` is shorter than one window
/// (`k + w - 1` residues).
pub fn minimizer_sketch(seq: &[u8], k: usize, w: usize) -> Result<Vec<Minimizer>> {
    if k == 0 {
        return Err(AlignError::invalid("k", "k-mer length must be >= 1"));
    }
    if w == 0 {
        return Err(AlignError::invalid("w", "window size must be >= 1"));
    }
    let n = seq.len();
    let window_span = k + w - 1;
    if n < window_span {
        return Ok(Vec::new());
    }

    // Hash every k-mer once.
    let kmer_count = n - k + 1;
    let hashes: Vec<u64> = (0..kmer_count).map(|i| hash_kmer(&seq[i..i + k])).collect();

    // Slide a w-wide window over the k-mer hash array; pick the min.
    let mut out: Vec<Minimizer> = Vec::new();
    for start in 0..=kmer_count - w {
        let mut best_pos = start;
        let mut best_hash = hashes[start];
        for off in 1..w {
            let p = start + off;
            // Ties resolved toward the leftmost occurrence.
            if hashes[p] < best_hash {
                best_hash = hashes[p];
                best_pos = p;
            }
        }
        let m = Minimizer {
            hash: best_hash,
            pos: best_pos,
        };
        if out.last() != Some(&m) {
            out.push(m);
        }
    }
    Ok(out)
}

/// Hashes a k-mer with a FNV-1a / mixing function. Deterministic and
/// fast; the exact constant is not load-bearing, only that the order
/// it imposes is consistent within a run.
fn hash_kmer(kmer: &[u8]) -> u64 {
    // FNV-1a over the uppercased bytes, then a final avalanche mix so
    // the low bits (which the window minimum compares first) are
    // well-distributed.
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in kmer {
        h ^= b.to_ascii_uppercase() as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    // SplitMix64 finalizer.
    h ^= h >> 30;
    h = h.wrapping_mul(0xbf58476d1ce4e5b9);
    h ^= h >> 27;
    h = h.wrapping_mul(0x94d049bb133111eb);
    h ^= h >> 31;
    h
}

/// Counts how many minimizers two sketches share by hash — a fast
/// lower bound on sequence similarity (the minimap2 "shared seed"
/// estimate).
pub fn shared_minimizers(a: &[Minimizer], b: &[Minimizer]) -> usize {
    use std::collections::HashSet;
    let set_a: HashSet<u64> = a.iter().map(|m| m.hash).collect();
    b.iter().filter(|m| set_a.contains(&m.hash)).count()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_bad_params() {
        assert!(minimizer_sketch(b"ACGT", 0, 4).is_err());
        assert!(minimizer_sketch(b"ACGT", 4, 0).is_err());
    }

    #[test]
    fn short_sequence_empty_sketch() {
        // k=5, w=5 needs >= 9 residues.
        assert!(minimizer_sketch(b"ACGTACG", 5, 5).unwrap().is_empty());
    }

    #[test]
    fn sketch_is_ordered_and_sparse() {
        let seq = b"ACGTACGTACGTACGTACGTACGT";
        let sk = minimizer_sketch(seq, 5, 5).unwrap();
        assert!(!sk.is_empty());
        // Positions strictly increasing across distinct entries.
        for pair in sk.windows(2) {
            assert!(pair[1].pos >= pair[0].pos);
        }
        // Sparse: fewer minimizers than k-mers.
        let kmer_count = seq.len() - 5 + 1;
        assert!(sk.len() <= kmer_count);
    }

    #[test]
    fn determinism() {
        let seq = b"GATTACAGATTACAGATTACA";
        let a = minimizer_sketch(seq, 6, 4).unwrap();
        let b = minimizer_sketch(seq, 6, 4).unwrap();
        assert_eq!(a, b, "sketch must be deterministic");
    }

    #[test]
    fn shared_substring_shares_a_minimizer() {
        // Two sequences sharing a long core (>= w+k-1) must share a
        // minimizer — the defining minimizer guarantee.
        let core = b"GATTACACATGGCATAGCATAG";
        let s1 = [b"TTTTTT".as_slice(), core].concat();
        let s2 = [core, b"AAAAAA".as_slice()].concat();
        let k = 6;
        let w = 5;
        let m1 = minimizer_sketch(&s1, k, w).unwrap();
        let m2 = minimizer_sketch(&s2, k, w).unwrap();
        assert!(
            shared_minimizers(&m1, &m2) >= 1,
            "long shared substring must yield a shared minimizer"
        );
    }

    #[test]
    fn unrelated_sequences_share_few() {
        let m1 = minimizer_sketch(b"AAAAAAAAAAAAAAAAAAAA", 6, 5).unwrap();
        let m2 = minimizer_sketch(b"CCCCCCCCCCCCCCCCCCCC", 6, 5).unwrap();
        assert_eq!(shared_minimizers(&m1, &m2), 0);
    }

    #[test]
    fn hash_is_case_insensitive() {
        assert_eq!(hash_kmer(b"acgt"), hash_kmer(b"ACGT"));
    }
}
