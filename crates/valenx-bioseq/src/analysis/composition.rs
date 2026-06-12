//! Sequence composition: GC content, GC skew, residue / dinucleotide
//! frequencies and Shannon entropy.

use crate::error::{BioseqError, Result};
use crate::seq::{Seq, SeqKind};
use std::collections::BTreeMap;

/// Overall GC fraction (0.0–1.0) of a nucleotide sequence.
///
/// Counts `G`, `C`, and the ambiguity code `S` (G or C). The
/// denominator is the count of unambiguous A/C/G/T(U) plus S/W — so
/// `N` runs do not deflate the result. Returns
/// [`BioseqError::Invalid`] for a protein sequence and `0.0` for an
/// all-ambiguous / empty sequence.
pub fn gc_content(seq: &Seq) -> Result<f64> {
    if seq.kind() == SeqKind::Protein {
        return Err(BioseqError::invalid(
            "kind",
            "GC content needs a nucleotide sequence",
        ));
    }
    let (gc, at) = gc_at_counts(seq.as_bytes());
    let total = gc + at;
    Ok(if total == 0 {
        0.0
    } else {
        gc as f64 / total as f64
    })
}

/// AT fraction — the complement of [`gc_content`] over the same
/// unambiguous denominator.
pub fn at_content(seq: &Seq) -> Result<f64> {
    Ok(1.0 - gc_content(seq)?)
}

/// Count of **purine** bases — adenine + guanine — in a nucleotide sequence.
/// Returns an error for a protein. A distinct grouping from GC (purines span both the
/// GC and AT pairs); for an unambiguous ACGT/ACGU sequence
/// `purine_count + pyrimidine_count == len`. Ambiguity codes are not counted.
pub fn purine_count(seq: &Seq) -> Result<usize> {
    if seq.kind() == SeqKind::Protein {
        return Err(BioseqError::invalid(
            "kind",
            "purine count needs a nucleotide sequence",
        ));
    }
    Ok(seq.count(b'A') + seq.count(b'G'))
}

/// Count of **pyrimidine** bases — cytosine + thymine (DNA) or uracil (RNA) — in a
/// nucleotide sequence. Returns an error for a protein; ambiguity codes are not
/// counted. The complement grouping to [`purine_count`].
pub fn pyrimidine_count(seq: &Seq) -> Result<usize> {
    if seq.kind() == SeqKind::Protein {
        return Err(BioseqError::invalid(
            "kind",
            "pyrimidine count needs a nucleotide sequence",
        ));
    }
    Ok(seq.count(b'C') + seq.count(b'T') + seq.count(b'U'))
}

/// `(gc_count, at_count)` over A/C/G/T(U) plus the two-base ambiguity
/// codes `S` (→GC) and `W` (→AT).
fn gc_at_counts(bytes: &[u8]) -> (usize, usize) {
    let mut gc = 0;
    let mut at = 0;
    for &b in bytes {
        match b.to_ascii_uppercase() {
            b'G' | b'C' | b'S' => gc += 1,
            b'A' | b'T' | b'U' | b'W' => at += 1,
            _ => {}
        }
    }
    (gc, at)
}

/// GC content in a sliding window — returns one value per window
/// position. `window` is the window size; `step` is the stride.
///
/// The output `Vec` has `floor((len - window) / step) + 1` entries.
/// Returns [`BioseqError::Invalid`] for a zero window/step or a
/// window larger than the sequence.
pub fn gc_sliding_window(seq: &Seq, window: usize, step: usize) -> Result<Vec<f64>> {
    if seq.kind() == SeqKind::Protein {
        return Err(BioseqError::invalid(
            "kind",
            "GC content needs a nucleotide sequence",
        ));
    }
    if window == 0 || step == 0 {
        return Err(BioseqError::invalid(
            "window",
            "window and step must be > 0",
        ));
    }
    let bytes = seq.as_bytes();
    if window > bytes.len() {
        return Err(BioseqError::invalid(
            "window",
            format!("window {window} exceeds sequence length {}", bytes.len()),
        ));
    }
    let mut out = Vec::new();
    let mut start = 0;
    while start + window <= bytes.len() {
        let (gc, at) = gc_at_counts(&bytes[start..start + window]);
        let total = gc + at;
        out.push(if total == 0 {
            0.0
        } else {
            gc as f64 / total as f64
        });
        start += step;
    }
    Ok(out)
}

/// GC skew `(G − C) / (G + C)` over the whole sequence.
///
/// GC skew is used to locate replication origins/termini in bacterial
/// genomes. Returns `0.0` when `G + C == 0`.
pub fn gc_skew(seq: &Seq) -> Result<f64> {
    if seq.kind() == SeqKind::Protein {
        return Err(BioseqError::invalid(
            "kind",
            "GC skew needs a nucleotide sequence",
        ));
    }
    let g = seq.count(b'G') as f64;
    let c = seq.count(b'C') as f64;
    Ok(if g + c == 0.0 { 0.0 } else { (g - c) / (g + c) })
}

/// Cumulative GC skew per window — the running sum of per-window
/// `(G−C)/(G+C)` values, the classic skew plot. See
/// [`gc_sliding_window`] for the windowing.
pub fn cumulative_gc_skew(seq: &Seq, window: usize, step: usize) -> Result<Vec<f64>> {
    if seq.kind() == SeqKind::Protein {
        return Err(BioseqError::invalid(
            "kind",
            "GC skew needs a nucleotide sequence",
        ));
    }
    if window == 0 || step == 0 {
        return Err(BioseqError::invalid(
            "window",
            "window and step must be > 0",
        ));
    }
    let bytes = seq.as_bytes();
    if window > bytes.len() {
        return Err(BioseqError::invalid(
            "window",
            format!("window {window} exceeds length {}", bytes.len()),
        ));
    }
    let mut out = Vec::new();
    let mut running = 0.0;
    let mut start = 0;
    while start + window <= bytes.len() {
        let slice = &bytes[start..start + window];
        let g = slice
            .iter()
            .filter(|&&b| b.eq_ignore_ascii_case(&b'G'))
            .count() as f64;
        let c = slice
            .iter()
            .filter(|&&b| b.eq_ignore_ascii_case(&b'C'))
            .count() as f64;
        let skew = if g + c == 0.0 { 0.0 } else { (g - c) / (g + c) };
        running += skew;
        out.push(running);
        start += step;
    }
    Ok(out)
}

/// Counts of every residue that appears, keyed by uppercase residue.
pub fn residue_counts(seq: &Seq) -> BTreeMap<char, usize> {
    let mut map: BTreeMap<char, usize> = BTreeMap::new();
    for b in seq.iter() {
        *map.entry(b as char).or_insert(0) += 1;
    }
    map
}

/// Residue frequencies (each count divided by the sequence length).
/// Empty for an empty sequence.
pub fn residue_frequencies(seq: &Seq) -> BTreeMap<char, f64> {
    let n = seq.len();
    let mut map: BTreeMap<char, f64> = BTreeMap::new();
    if n == 0 {
        return map;
    }
    for (k, v) in residue_counts(seq) {
        map.insert(k, v as f64 / n as f64);
    }
    map
}

/// Counts of every overlapping dinucleotide (length-2 window, step 1).
/// Keyed by the uppercase two-character string.
pub fn dinucleotide_counts(seq: &Seq) -> BTreeMap<String, usize> {
    let mut map: BTreeMap<String, usize> = BTreeMap::new();
    let bytes = seq.as_bytes();
    for w in bytes.windows(2) {
        let k = format!("{}{}", w[0] as char, w[1] as char);
        *map.entry(k).or_insert(0) += 1;
    }
    map
}

/// Dinucleotide frequencies — [`dinucleotide_counts`] normalized by the
/// number of dinucleotide windows.
pub fn dinucleotide_frequencies(seq: &Seq) -> BTreeMap<String, f64> {
    let counts = dinucleotide_counts(seq);
    let total: usize = counts.values().sum();
    let mut map: BTreeMap<String, f64> = BTreeMap::new();
    if total == 0 {
        return map;
    }
    for (k, v) in counts {
        map.insert(k, v as f64 / total as f64);
    }
    map
}

/// Shannon entropy of the residue distribution, in bits.
///
/// `H = -Σ pᵢ·log₂(pᵢ)`. The maximum is `log₂(alphabet_size)` — about
/// 2 bits for DNA, ~4.32 for protein. Returns `0.0` for an empty or
/// single-residue-type sequence.
pub fn shannon_entropy(seq: &Seq) -> f64 {
    let n = seq.len();
    if n == 0 {
        return 0.0;
    }
    let counts = residue_counts(seq);
    let mut h = 0.0;
    for &c in counts.values() {
        if c == 0 {
            continue;
        }
        let p = c as f64 / n as f64;
        h -= p * p.log2();
    }
    h
}

/// Hamming distance — the number of positions at which two **equal-length** sequences
/// differ. A foundational position-exact pairwise metric (works on any [`SeqKind`]).
/// Returns an error if the two sequences have different lengths.
pub fn hamming_distance(a: &Seq, b: &Seq) -> Result<usize> {
    if a.len() != b.len() {
        return Err(BioseqError::invalid(
            "sequences",
            format!("lengths must match: {} vs {}", a.len(), b.len()),
        ));
    }
    Ok(a.as_bytes()
        .iter()
        .zip(b.as_bytes().iter())
        .filter(|(x, y)| x != y)
        .count())
}

/// N50 — the assembly-contiguity statistic over a multiset of contig/fragment `lengths`:
/// sort descending, accumulate, and return the length at which the running total first
/// reaches at least half the grand total. Higher N50 = longer contiguity. Returns `0` for
/// an empty input. A set-level statistic, so it takes a length slice (not a [`Seq`]).
pub fn n50(lengths: &[usize]) -> usize {
    let total: usize = lengths.iter().sum();
    if total == 0 {
        return 0;
    }
    let mut sorted = lengths.to_vec();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    let mut cumulative = 0usize;
    for &len in &sorted {
        cumulative += len;
        if 2 * cumulative >= total {
            return len;
        }
    }
    0
}

/// L50 — the assembly-fragmentation companion to [`n50`]: the number of the longest
/// contigs needed for their cumulative length to reach at least half the grand total.
/// N50 is the *length* at that point; L50 is the *count* (1-indexed). Returns `0` for an
/// empty input. A count, not a transform of N50 — e.g. `n50([5,5])` = 5 but `l50([5,5])` = 1.
pub fn l50(lengths: &[usize]) -> usize {
    let total: usize = lengths.iter().sum();
    if total == 0 {
        return 0;
    }
    let mut sorted = lengths.to_vec();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    let mut cumulative = 0usize;
    let mut count = 0usize;
    for &len in &sorted {
        count += 1;
        cumulative += len;
        if 2 * cumulative >= total {
            return count;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_content_basic() {
        let s = Seq::new(SeqKind::Dna, "GGCC").unwrap();
        assert!((gc_content(&s).unwrap() - 1.0).abs() < 1e-12);
        let s = Seq::new(SeqKind::Dna, "ATAT").unwrap();
        assert!((gc_content(&s).unwrap() - 0.0).abs() < 1e-12);
        let s = Seq::new(SeqKind::Dna, "ATGC").unwrap();
        assert!((gc_content(&s).unwrap() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn gc_ignores_n_runs() {
        let s = Seq::new(SeqKind::Dna, "GCNNNN").unwrap();
        // 2 of 2 unambiguous bases are G/C.
        assert!((gc_content(&s).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn at_content_complements_gc() {
        let s = Seq::new(SeqKind::Dna, "ATGC").unwrap();
        assert!((at_content(&s).unwrap() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn purine_pyrimidine_counts() {
        // ATGC: purines A,G = 2; pyrimidines T,C = 2.
        let dna = Seq::new(SeqKind::Dna, "ATGC").unwrap();
        assert_eq!(purine_count(&dna).unwrap(), 2);
        assert_eq!(pyrimidine_count(&dna).unwrap(), 2);
        assert_eq!(
            purine_count(&Seq::new(SeqKind::Dna, "AAGG").unwrap()).unwrap(),
            4
        );
        assert_eq!(
            pyrimidine_count(&Seq::new(SeqKind::Dna, "CCTT").unwrap()).unwrap(),
            4
        );
        // RNA: U is a pyrimidine.
        let rna = Seq::new(SeqKind::Rna, "AUGC").unwrap();
        assert_eq!(purine_count(&rna).unwrap(), 2);
        assert_eq!(pyrimidine_count(&rna).unwrap(), 2);
        // Non-tautological partition over an unambiguous sequence: purine + pyrimidine = len.
        let s = Seq::new(SeqKind::Dna, "ATGCATGCAT").unwrap();
        assert_eq!(
            purine_count(&s).unwrap() + pyrimidine_count(&s).unwrap(),
            s.len()
        );
        // Proteins are rejected (A, G are also amino-acid codes).
        let prot = Seq::new(SeqKind::Protein, "MAGK").unwrap();
        assert!(purine_count(&prot).is_err());
        assert!(pyrimidine_count(&prot).is_err());
    }

    #[test]
    fn protein_rejected() {
        let p = Seq::new(SeqKind::Protein, "MKVL").unwrap();
        assert!(gc_content(&p).is_err());
        assert!(gc_skew(&p).is_err());
    }

    #[test]
    fn hamming_distance_counts_differing_positions() {
        let a = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        // one mismatch (position 1: C vs G).
        assert_eq!(
            hamming_distance(&a, &Seq::new(SeqKind::Dna, "AGGT").unwrap()).unwrap(),
            1
        );
        // all four positions differ.
        assert_eq!(
            hamming_distance(
                &Seq::new(SeqKind::Dna, "AAAA").unwrap(),
                &Seq::new(SeqKind::Dna, "TTTT").unwrap()
            )
            .unwrap(),
            4
        );
        // identical → 0.
        assert_eq!(hamming_distance(&a, &a).unwrap(), 0);
        // works on protein (any SeqKind).
        let p1 = Seq::new(SeqKind::Protein, "MVKL").unwrap();
        let p2 = Seq::new(SeqKind::Protein, "MVQL").unwrap();
        assert_eq!(hamming_distance(&p1, &p2).unwrap(), 1);
        // unequal length → error.
        assert!(hamming_distance(&a, &Seq::new(SeqKind::Dna, "AC").unwrap()).is_err());
        // non-tautological triangle inequality: d(x,z) ≤ d(x,y) + d(y,z).
        let x = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        let y = Seq::new(SeqKind::Dna, "AGGT").unwrap();
        let z = Seq::new(SeqKind::Dna, "AGGA").unwrap();
        let dxz = hamming_distance(&x, &z).unwrap();
        let dxy = hamming_distance(&x, &y).unwrap();
        let dyz = hamming_distance(&y, &z).unwrap();
        assert!(dxz <= dxy + dyz, "triangle inequality");
    }

    #[test]
    fn n50_worked_examples() {
        // [2,3,4,5,6]: total 20, half 10, sorted desc [6,5,4,3,2], cumsum 6,11(≥10) → 5.
        assert_eq!(n50(&[2, 3, 4, 5, 6]), 5);
        assert_eq!(n50(&[10]), 10);
        assert_eq!(n50(&[5, 5]), 5);
        assert_eq!(n50(&[100, 1, 1]), 100);
        assert_eq!(n50(&[]), 0);
        // Input order doesn't matter (the algorithm sorts).
        assert_eq!(
            n50(&[3, 1, 4, 1, 5, 9, 2, 6]),
            n50(&[9, 6, 5, 4, 3, 2, 1, 1])
        );
        // Non-tautological: N50 is one of the input lengths for non-empty input.
        let input = [7, 3, 15, 8];
        assert!(input.contains(&n50(&input)));
    }

    #[test]
    fn l50_worked_examples() {
        // [2,3,4,5,6]: sorted [6,5,4,3,2], cum 6 (1 contig, <10), 11 (2 contigs, ≥10) → 2.
        assert_eq!(l50(&[2, 3, 4, 5, 6]), 2);
        assert_eq!(l50(&[10]), 1);
        // n50([5,5]) = 5 but l50([5,5]) = 1 — distinct quantities (length vs count).
        assert_eq!(l50(&[5, 5]), 1);
        assert_eq!(l50(&[100, 1, 1]), 1);
        assert_eq!(l50(&[]), 0);
        // Non-tautological tie: the L50-th largest contig equals N50.
        let input = [2, 3, 4, 5, 6];
        let mut sorted = input.to_vec();
        sorted.sort_unstable_by(|a, b| b.cmp(a));
        assert_eq!(sorted[l50(&input) - 1], n50(&input));
    }

    #[test]
    fn sliding_window_shape() {
        let s = Seq::new(SeqKind::Dna, "GGGGCCCCAAAATTTT").unwrap(); // 16 nt
        let w = gc_sliding_window(&s, 4, 4).unwrap();
        assert_eq!(w.len(), 4);
        assert!((w[0] - 1.0).abs() < 1e-12); // GGGG
        assert!((w[2] - 0.0).abs() < 1e-12); // AAAA
        assert!(gc_sliding_window(&s, 0, 1).is_err());
        assert!(gc_sliding_window(&s, 99, 1).is_err());
    }

    #[test]
    fn gc_skew_sign() {
        let s = Seq::new(SeqKind::Dna, "GGGC").unwrap();
        // (3-1)/(3+1) = 0.5.
        assert!((gc_skew(&s).unwrap() - 0.5).abs() < 1e-12);
        let s = Seq::new(SeqKind::Dna, "AAAA").unwrap();
        assert_eq!(gc_skew(&s).unwrap(), 0.0);
    }

    #[test]
    fn cumulative_skew_is_running_sum() {
        let s = Seq::new(SeqKind::Dna, "GGCC").unwrap();
        let c = cumulative_gc_skew(&s, 2, 2).unwrap();
        // window GG: skew +1 ; window CC: skew -1 -> cumulative [1, 0].
        assert_eq!(c.len(), 2);
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 0.0).abs() < 1e-12);
    }

    #[test]
    fn residue_counts_and_freqs() {
        let s = Seq::new(SeqKind::Dna, "AAGC").unwrap();
        let counts = residue_counts(&s);
        assert_eq!(counts[&'A'], 2);
        assert_eq!(counts[&'G'], 1);
        let freqs = residue_frequencies(&s);
        assert!((freqs[&'A'] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn dinucleotide_counts_overlap() {
        let s = Seq::new(SeqKind::Dna, "ATAT").unwrap();
        let d = dinucleotide_counts(&s);
        // windows: AT, TA, AT.
        assert_eq!(d["AT"], 2);
        assert_eq!(d["TA"], 1);
        let f = dinucleotide_frequencies(&s);
        assert!((f["AT"] - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn entropy_extremes() {
        // Uniform 4-letter DNA -> 2 bits.
        let s = Seq::new(SeqKind::Dna, "ACGT").unwrap();
        assert!((shannon_entropy(&s) - 2.0).abs() < 1e-12);
        // Single residue type -> 0 bits.
        let s = Seq::new(SeqKind::Dna, "AAAA").unwrap();
        assert!(shannon_entropy(&s).abs() < 1e-12);
    }
}
