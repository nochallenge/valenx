//! k-mer counting, k-mer spectrum, and k-mer-based sequence distance.
//!
//! A k-mer is a length-`k` substring. k-mer counts underpin genome
//! assembly, alignment-free comparison, and read error correction.
//! This module counts overlapping k-mers, builds the count histogram
//! (the "k-mer spectrum"), and computes a cosine-style distance
//! between two sequences' k-mer profiles.

use crate::error::{BioseqError, Result};
use crate::seq::Seq;
use std::collections::BTreeMap;

/// Counts every overlapping length-`k` k-mer in `seq`.
///
/// k-mers are keyed by their uppercase string. Returns
/// [`BioseqError::Invalid`] if `k` is `0` or larger than the
/// sequence.
pub fn count_kmers(seq: &Seq, k: usize) -> Result<BTreeMap<String, usize>> {
    if k == 0 {
        return Err(BioseqError::invalid("k", "k must be > 0"));
    }
    let bytes = seq.as_bytes();
    if k > bytes.len() {
        return Err(BioseqError::invalid(
            "k",
            format!("k={k} exceeds sequence length {}", bytes.len()),
        ));
    }
    let mut map: BTreeMap<String, usize> = BTreeMap::new();
    for w in bytes.windows(k) {
        let key = std::str::from_utf8(w).expect("residues are ASCII").to_string();
        *map.entry(key).or_insert(0) += 1;
    }
    Ok(map)
}

/// The number of overlapping length-`k` windows in a sequence of
/// length `n`: `n - k + 1` (or `0` if `k > n`).
pub fn kmer_window_count(n: usize, k: usize) -> usize {
    if k == 0 || k > n {
        0
    } else {
        n - k + 1
    }
}

/// The k-mer spectrum: a histogram mapping a multiplicity (how many
/// times a k-mer occurs) to how many distinct k-mers have that
/// multiplicity.
///
/// The spectrum's shape diagnoses sequencing-error and
/// repeat-content; the low-multiplicity peak is mostly errors.
pub fn kmer_spectrum(seq: &Seq, k: usize) -> Result<BTreeMap<usize, usize>> {
    let counts = count_kmers(seq, k)?;
    let mut spectrum: BTreeMap<usize, usize> = BTreeMap::new();
    for &mult in counts.values() {
        *spectrum.entry(mult).or_insert(0) += 1;
    }
    Ok(spectrum)
}

/// The `top_n` most frequent k-mers, as `(kmer, count)` pairs sorted by
/// descending count (ties broken lexicographically for determinism).
pub fn most_frequent_kmers(seq: &Seq, k: usize, top_n: usize) -> Result<Vec<(String, usize)>> {
    let counts = count_kmers(seq, k)?;
    let mut pairs: Vec<(String, usize)> = counts.into_iter().collect();
    pairs.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    pairs.truncate(top_n);
    Ok(pairs)
}

/// Number of distinct k-mers in a sequence.
pub fn distinct_kmer_count(seq: &Seq, k: usize) -> Result<usize> {
    Ok(count_kmers(seq, k)?.len())
}

/// k-mer-based distance between two sequences in `[0, 1]`.
///
/// Treats each sequence's k-mer counts as a vector and returns
/// `1 − cosine_similarity` — `0` for identical k-mer profiles, `1` for
/// profiles sharing no k-mer. This is an alignment-free distance used
/// for fast clustering. Both sequences must admit length-`k` k-mers.
pub fn kmer_distance(a: &Seq, b: &Seq, k: usize) -> Result<f64> {
    let ca = count_kmers(a, k)?;
    let cb = count_kmers(b, k)?;
    let mut dot = 0.0f64;
    let mut na = 0.0f64;
    let mut nb = 0.0f64;
    for &v in ca.values() {
        na += (v as f64) * (v as f64);
    }
    for &v in cb.values() {
        nb += (v as f64) * (v as f64);
    }
    for (kmer, &va) in &ca {
        if let Some(&vb) = cb.get(kmer) {
            dot += (va as f64) * (vb as f64);
        }
    }
    if na == 0.0 || nb == 0.0 {
        return Ok(1.0);
    }
    let cosine = dot / (na.sqrt() * nb.sqrt());
    Ok((1.0 - cosine).clamp(0.0, 1.0))
}

/// Jaccard distance on k-mer *sets* (ignoring multiplicity).
///
/// `1 − |A∩B| / |A∪B|`. A complement to [`kmer_distance`] when
/// presence/absence matters more than abundance.
pub fn kmer_jaccard_distance(a: &Seq, b: &Seq, k: usize) -> Result<f64> {
    let ca = count_kmers(a, k)?;
    let cb = count_kmers(b, k)?;
    let mut intersection = 0usize;
    for kmer in ca.keys() {
        if cb.contains_key(kmer) {
            intersection += 1;
        }
    }
    let union = ca.len() + cb.len() - intersection;
    if union == 0 {
        return Ok(0.0);
    }
    Ok(1.0 - intersection as f64 / union as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::SeqKind;

    #[test]
    fn count_overlapping_kmers() {
        let s = Seq::new(SeqKind::Dna, "AAAA").unwrap();
        let c = count_kmers(&s, 2).unwrap();
        // windows: AA, AA, AA.
        assert_eq!(c["AA"], 3);
        assert_eq!(c.len(), 1);
    }

    #[test]
    fn window_count_helper() {
        assert_eq!(kmer_window_count(10, 3), 8);
        assert_eq!(kmer_window_count(3, 3), 1);
        assert_eq!(kmer_window_count(2, 3), 0);
    }

    #[test]
    fn invalid_k_errors() {
        let s = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        assert!(count_kmers(&s, 0).is_err());
        assert!(count_kmers(&s, 99).is_err());
    }

    #[test]
    fn distinct_count() {
        let s = Seq::new(SeqKind::Dna, "ACGTACGT").unwrap();
        // 3-mers: ACG CGT GTA TAC ACG CGT -> distinct: ACG CGT GTA TAC.
        assert_eq!(distinct_kmer_count(&s, 3).unwrap(), 4);
    }

    #[test]
    fn spectrum_histogram() {
        let s = Seq::new(SeqKind::Dna, "ACGTACGT").unwrap();
        let spec = kmer_spectrum(&s, 3).unwrap();
        // ACG x2, CGT x2, GTA x1, TAC x1 -> mult 2 has 2 kmers, mult 1 has 2.
        assert_eq!(spec[&2], 2);
        assert_eq!(spec[&1], 2);
    }

    #[test]
    fn most_frequent() {
        let s = Seq::new(SeqKind::Dna, "AAAAGC").unwrap();
        let top = most_frequent_kmers(&s, 2, 2).unwrap();
        assert_eq!(top[0], ("AA".to_string(), 3));
        assert_eq!(top.len(), 2);
    }

    #[test]
    fn distance_identical_is_zero() {
        let a = Seq::new(SeqKind::Dna, "ACGTACGT").unwrap();
        let b = Seq::new(SeqKind::Dna, "ACGTACGT").unwrap();
        assert!(kmer_distance(&a, &b, 3).unwrap().abs() < 1e-9);
        assert!(kmer_jaccard_distance(&a, &b, 3).unwrap().abs() < 1e-9);
    }

    #[test]
    fn distance_disjoint_is_one() {
        let a = Seq::new(SeqKind::Dna, "AAAAAA").unwrap();
        let b = Seq::new(SeqKind::Dna, "GGGGGG").unwrap();
        assert!((kmer_distance(&a, &b, 3).unwrap() - 1.0).abs() < 1e-9);
        assert!((kmer_jaccard_distance(&a, &b, 3).unwrap() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn distance_partial_overlap_is_between() {
        let a = Seq::new(SeqKind::Dna, "ACGTACGT").unwrap();
        let b = Seq::new(SeqKind::Dna, "ACGTAAAA").unwrap();
        let d = kmer_distance(&a, &b, 3).unwrap();
        assert!(d > 0.0 && d < 1.0, "got {d}");
    }
}
