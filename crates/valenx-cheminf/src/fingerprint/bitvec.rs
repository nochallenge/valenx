//! A fixed-length bit-vector — the common representation behind every
//! fingerprint kind in this crate.
//!
//! [`FingerprintBits`] stores `n_bits` bits packed into `u64` words. It
//! supports the population-count and bitwise-overlap operations the
//! Tanimoto / Dice similarity coefficients need, plus folding (hashing
//! a sparse feature id into a dense bit).

use serde::{Deserialize, Serialize};

/// A dense, fixed-length bit fingerprint.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FingerprintBits {
    words: Vec<u64>,
    n_bits: usize,
}

impl FingerprintBits {
    /// A zeroed fingerprint of `n_bits` bits. Panics on `n_bits == 0`.
    pub fn new(n_bits: usize) -> Self {
        assert!(n_bits > 0, "fingerprint must have at least one bit");
        FingerprintBits {
            words: vec![0; n_bits.div_ceil(64)],
            n_bits,
        }
    }

    /// Number of bits in the fingerprint.
    pub fn len(&self) -> usize {
        self.n_bits
    }

    /// `true` — a fingerprint always has ≥ 1 bit. Present so Clippy is
    /// satisfied alongside [`len`](Self::len).
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Set bit `index` (no-op if out of range).
    pub fn set(&mut self, index: usize) {
        if index < self.n_bits {
            self.words[index / 64] |= 1u64 << (index % 64);
        }
    }

    /// Fold an arbitrary feature id (a hash) into a bit and set it.
    /// This is how circular / path features map into a dense vector.
    pub fn set_hashed(&mut self, feature: u64) {
        let bit = (feature % self.n_bits as u64) as usize;
        self.set(bit);
    }

    /// Read bit `index`.
    pub fn get(&self, index: usize) -> bool {
        index < self.n_bits && (self.words[index / 64] >> (index % 64)) & 1 == 1
    }

    /// Number of set bits.
    pub fn count_ones(&self) -> usize {
        self.words.iter().map(|w| w.count_ones() as usize).sum()
    }

    /// Number of bits set in *both* fingerprints (`|A ∩ B|`). Returns
    /// `0` if the lengths differ.
    pub fn intersection_count(&self, other: &Self) -> usize {
        if self.n_bits != other.n_bits {
            return 0;
        }
        self.words
            .iter()
            .zip(&other.words)
            .map(|(a, b)| (a & b).count_ones() as usize)
            .sum()
    }

    /// Number of bits set in *either* fingerprint (`|A ∪ B|`).
    pub fn union_count(&self, other: &Self) -> usize {
        if self.n_bits != other.n_bits {
            return 0;
        }
        self.words
            .iter()
            .zip(&other.words)
            .map(|(a, b)| (a | b).count_ones() as usize)
            .sum()
    }

    /// Bitwise OR into `self` (used to merge feature sets).
    pub fn or_with(&mut self, other: &Self) {
        if self.n_bits != other.n_bits {
            return;
        }
        for (a, b) in self.words.iter_mut().zip(&other.words) {
            *a |= b;
        }
    }

    /// Fraction of bits set — fingerprint density. A cheap, sometimes
    /// useful descriptor.
    pub fn density(&self) -> f64 {
        self.count_ones() as f64 / self.n_bits as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_get_count() {
        let mut fp = FingerprintBits::new(128);
        assert_eq!(fp.len(), 128);
        assert_eq!(fp.count_ones(), 0);
        fp.set(0);
        fp.set(63);
        fp.set(64);
        fp.set(127);
        assert!(fp.get(0));
        assert!(fp.get(64));
        assert!(!fp.get(1));
        assert_eq!(fp.count_ones(), 4);
    }

    #[test]
    fn intersection_and_union() {
        let mut a = FingerprintBits::new(64);
        let mut b = FingerprintBits::new(64);
        a.set(1);
        a.set(2);
        a.set(3);
        b.set(2);
        b.set(3);
        b.set(4);
        assert_eq!(a.intersection_count(&b), 2);
        assert_eq!(a.union_count(&b), 4);
    }

    #[test]
    fn mismatched_lengths_are_safe() {
        let a = FingerprintBits::new(64);
        let b = FingerprintBits::new(128);
        assert_eq!(a.intersection_count(&b), 0);
        assert_eq!(a.union_count(&b), 0);
    }

    #[test]
    fn hashed_folding_stays_in_range() {
        let mut fp = FingerprintBits::new(100);
        fp.set_hashed(123_456_789);
        fp.set_hashed(u64::MAX);
        assert!(fp.count_ones() >= 1);
    }

    #[test]
    fn or_with_merges() {
        let mut a = FingerprintBits::new(64);
        let mut b = FingerprintBits::new(64);
        a.set(1);
        b.set(2);
        a.or_with(&b);
        assert!(a.get(1) && a.get(2));
    }
}
