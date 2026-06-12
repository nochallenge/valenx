//! K-mer index — hash every length-`k` substring to its positions.
//!
//! [`KmerIndex`] is the seeding data structure behind the
//! BLAST-class heuristic search ([`crate::search::seed`]). It can
//! index a single sequence or a whole database of sequences; a lookup
//! returns every `(sequence-id, offset)` where a given k-mer occurs.

use crate::error::{AlignError, Result};
use std::collections::HashMap;

/// A position of a k-mer occurrence: which indexed sequence, and the
/// 0-based start offset within it.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct KmerHit {
    /// Index of the sequence in the order it was added.
    pub seq_id: usize,
    /// 0-based offset of the k-mer start within that sequence.
    pub offset: usize,
}

/// An inverted index from k-mer → occurrence list.
///
/// K-mers are stored as owned `Vec<u8>` keys (uppercased). For typical
/// `k` (8–15) this is simple and fast enough; a 2-bit packed key would
/// save memory but is not needed for the v1.
#[derive(Clone, Debug, Default)]
pub struct KmerIndex {
    k: usize,
    map: HashMap<Vec<u8>, Vec<KmerHit>>,
    /// Lengths of the indexed sequences, by id.
    seq_lens: Vec<usize>,
}

impl KmerIndex {
    /// Creates an empty index for k-mers of length `k`.
    ///
    /// Returns [`AlignError::Invalid`] if `k == 0`.
    pub fn new(k: usize) -> Result<Self> {
        if k == 0 {
            return Err(AlignError::invalid("k", "k-mer length must be >= 1"));
        }
        Ok(KmerIndex {
            k,
            map: HashMap::new(),
            seq_lens: Vec::new(),
        })
    }

    /// Builds an index over a single sequence (assigned `seq_id` 0).
    pub fn build(seq: &[u8], k: usize) -> Result<Self> {
        let mut idx = Self::new(k)?;
        idx.add_sequence(seq);
        Ok(idx)
    }

    /// Builds an index over many sequences; `seq_id`s are the slice
    /// indices.
    pub fn build_many(seqs: &[&[u8]], k: usize) -> Result<Self> {
        let mut idx = Self::new(k)?;
        for s in seqs {
            idx.add_sequence(s);
        }
        Ok(idx)
    }

    /// The k-mer length.
    pub fn k(&self) -> usize {
        self.k
    }

    /// The number of sequences indexed so far.
    pub fn sequence_count(&self) -> usize {
        self.seq_lens.len()
    }

    /// Length of the indexed sequence with id `seq_id`.
    pub fn sequence_len(&self, seq_id: usize) -> Option<usize> {
        self.seq_lens.get(seq_id).copied()
    }

    /// The number of *distinct* k-mers in the index.
    pub fn distinct_kmers(&self) -> usize {
        self.map.len()
    }

    /// Adds one sequence, returning the `seq_id` it was assigned.
    ///
    /// A sequence shorter than `k` contributes no k-mers but still
    /// consumes an id (so ids stay aligned with the caller's list).
    pub fn add_sequence(&mut self, seq: &[u8]) -> usize {
        let seq_id = self.seq_lens.len();
        self.seq_lens.push(seq.len());
        if seq.len() >= self.k {
            for offset in 0..=seq.len() - self.k {
                let kmer: Vec<u8> = seq[offset..offset + self.k]
                    .iter()
                    .map(|b| b.to_ascii_uppercase())
                    .collect();
                self.map
                    .entry(kmer)
                    .or_default()
                    .push(KmerHit { seq_id, offset });
            }
        }
        seq_id
    }

    /// All occurrences of `kmer`. Returns an empty slice for an absent
    /// or wrong-length k-mer.
    pub fn lookup(&self, kmer: &[u8]) -> &[KmerHit] {
        if kmer.len() != self.k {
            return &[];
        }
        let up: Vec<u8> = kmer.iter().map(|b| b.to_ascii_uppercase()).collect();
        self.map.get(&up).map(Vec::as_slice).unwrap_or(&[])
    }

    /// `true` if `kmer` occurs at least once.
    pub fn contains(&self, kmer: &[u8]) -> bool {
        !self.lookup(kmer).is_empty()
    }

    /// Total occurrence count of `kmer` across every indexed sequence.
    pub fn count(&self, kmer: &[u8]) -> usize {
        self.lookup(kmer).len()
    }

    /// Every `(query-offset, hit)` pair where a k-mer of `query`
    /// matches the index — the raw seed list for seed-and-extend.
    ///
    /// The result is sorted by query offset then by hit.
    pub fn seed_query(&self, query: &[u8]) -> Vec<(usize, KmerHit)> {
        let mut seeds = Vec::new();
        if query.len() >= self.k {
            for q_off in 0..=query.len() - self.k {
                for &hit in self.lookup(&query[q_off..q_off + self.k]) {
                    seeds.push((q_off, hit));
                }
            }
        }
        seeds
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_zero_k() {
        assert!(KmerIndex::new(0).is_err());
    }

    #[test]
    fn single_sequence_index() {
        // "ACGTACGT", k=3 -> ACG appears twice (offsets 0, 4).
        let idx = KmerIndex::build(b"ACGTACGT", 3).unwrap();
        assert_eq!(idx.sequence_count(), 1);
        let hits = idx.lookup(b"ACG");
        assert_eq!(hits.len(), 2);
        assert_eq!(
            hits[0],
            KmerHit {
                seq_id: 0,
                offset: 0
            }
        );
        assert_eq!(
            hits[1],
            KmerHit {
                seq_id: 0,
                offset: 4
            }
        );
        assert!(idx.contains(b"CGT"));
        assert!(!idx.contains(b"TTT"));
    }

    #[test]
    fn case_insensitive() {
        // "acgtACGT" has 5 length-4 windows (acgt, cgtA, gtAC, tACG,
        // ACGT) but only TWO of them — the first and the last — spell
        // ACGT case-insensitively, so count(ACGT) == 2.
        let idx = KmerIndex::build(b"acgtACGT", 4).unwrap();
        assert!(idx.contains(b"ACGT"));
        assert!(idx.contains(b"acgt"));
        assert_eq!(idx.count(b"ACGT"), 2);
    }

    #[test]
    fn multi_sequence_ids() {
        let seqs: &[&[u8]] = &[b"AAAA", b"AACC", b"CCAA"];
        let idx = KmerIndex::build_many(seqs, 2).unwrap();
        assert_eq!(idx.sequence_count(), 3);
        // "AA" occurs in all three.
        let hits = idx.lookup(b"AA");
        let ids: Vec<usize> = hits.iter().map(|h| h.seq_id).collect();
        assert!(ids.contains(&0));
        assert!(ids.contains(&1));
        assert!(ids.contains(&2));
    }

    #[test]
    fn short_sequence_no_kmers_but_keeps_id() {
        let mut idx = KmerIndex::new(5).unwrap();
        let id0 = idx.add_sequence(b"AC"); // too short for k=5
        let id1 = idx.add_sequence(b"ACGTACGT");
        assert_eq!(id0, 0);
        assert_eq!(id1, 1);
        assert_eq!(idx.sequence_len(0), Some(2));
        assert!(idx.contains(b"ACGTA"));
    }

    #[test]
    fn seed_query_pairs() {
        let idx = KmerIndex::build(b"GATTACAGATTACA", 4).unwrap();
        let seeds = idx.seed_query(b"GATTACA");
        // "GATT" occurs at query 0; matches index offsets 0 and 7.
        assert!(seeds.iter().any(|&(q, h)| q == 0 && h.offset == 0));
        assert!(seeds.iter().any(|&(q, h)| q == 0 && h.offset == 7));
        assert!(!seeds.is_empty());
    }

    #[test]
    fn wrong_length_kmer_returns_empty() {
        let idx = KmerIndex::build(b"ACGTACGT", 3).unwrap();
        assert!(idx.lookup(b"AC").is_empty());
        assert!(idx.lookup(b"ACGT").is_empty());
    }
}
