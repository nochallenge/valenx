//! Fingerprint similarity coefficients.
//!
//! Given two [`FingerprintBits`] of equal length, these are the
//! standard chemical-similarity measures:
//!
//! - [`tanimoto`] — `|A ∩ B| / |A ∪ B|`, the Jaccard coefficient; the
//!   default for virtual screening (range 0..1);
//! - [`dice`] — `2|A ∩ B| / (|A| + |B|)`, the Sørensen-Dice
//!   coefficient (range 0..1, slightly more lenient than Tanimoto);
//! - [`cosine`] — `|A ∩ B| / √(|A|·|B|)`, the Ochiai coefficient;
//! - [`tversky`] — the asymmetric generalisation parameterised by
//!   `(α, β)` (`α = β = 1` recovers Tanimoto, `α = β = 0.5` recovers
//!   Dice).
//!
//! All return `1.0` for two empty (all-zero) fingerprints of equal
//! length — two molecules with no features are "identical" — and
//! `0.0` when the lengths differ.

use super::bitvec::FingerprintBits;

/// Tanimoto / Jaccard similarity, `|A ∩ B| / |A ∪ B|`.
pub fn tanimoto(a: &FingerprintBits, b: &FingerprintBits) -> f64 {
    if a.len() != b.len() {
        return 0.0;
    }
    let union = a.union_count(b);
    if union == 0 {
        return 1.0;
    }
    a.intersection_count(b) as f64 / union as f64
}

/// Sørensen-Dice similarity, `2|A ∩ B| / (|A| + |B|)`.
pub fn dice(a: &FingerprintBits, b: &FingerprintBits) -> f64 {
    if a.len() != b.len() {
        return 0.0;
    }
    let sum = a.count_ones() + b.count_ones();
    if sum == 0 {
        return 1.0;
    }
    2.0 * a.intersection_count(b) as f64 / sum as f64
}

/// Cosine / Ochiai similarity, `|A ∩ B| / √(|A|·|B|)`.
pub fn cosine(a: &FingerprintBits, b: &FingerprintBits) -> f64 {
    if a.len() != b.len() {
        return 0.0;
    }
    let (na, nb) = (a.count_ones(), b.count_ones());
    if na == 0 && nb == 0 {
        return 1.0;
    }
    if na == 0 || nb == 0 {
        return 0.0;
    }
    a.intersection_count(b) as f64 / ((na as f64) * (nb as f64)).sqrt()
}

/// Tversky asymmetric similarity:
/// `c / (c + α·(a−c) + β·(b−c))` where `c = |A ∩ B|`, `a = |A|`,
/// `b = |B|`. `α` weights the features unique to `A`, `β` those unique
/// to `B`.
pub fn tversky(a: &FingerprintBits, b: &FingerprintBits, alpha: f64, beta: f64) -> f64 {
    if a.len() != b.len() {
        return 0.0;
    }
    let c = a.intersection_count(b) as f64;
    let na = a.count_ones() as f64;
    let nb = b.count_ones() as f64;
    let denom = c + alpha * (na - c) + beta * (nb - c);
    if denom <= 0.0 {
        return 1.0;
    }
    c / denom
}

/// Tanimoto *distance*, `1 − tanimoto` — a proper metric, handy for
/// clustering.
pub fn tanimoto_distance(a: &FingerprintBits, b: &FingerprintBits) -> f64 {
    1.0 - tanimoto(a, b)
}

/// For a query fingerprint and a library, return `(index, similarity)`
/// of the most similar library entry by Tanimoto. `None` for an empty
/// library.
pub fn nearest_neighbor(
    query: &FingerprintBits,
    library: &[FingerprintBits],
) -> Option<(usize, f64)> {
    library
        .iter()
        .enumerate()
        .map(|(i, fp)| (i, tanimoto(query, fp)))
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fp(bits: &[usize], len: usize) -> FingerprintBits {
        let mut f = FingerprintBits::new(len);
        for &b in bits {
            f.set(b);
        }
        f
    }

    #[test]
    fn identical_is_one() {
        let a = fp(&[1, 2, 3], 64);
        assert!((tanimoto(&a, &a) - 1.0).abs() < 1e-12);
        assert!((dice(&a, &a) - 1.0).abs() < 1e-12);
        assert!((cosine(&a, &a) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn disjoint_is_zero() {
        let a = fp(&[1, 2], 64);
        let b = fp(&[10, 11], 64);
        assert_eq!(tanimoto(&a, &b), 0.0);
        assert_eq!(dice(&a, &b), 0.0);
        assert_eq!(cosine(&a, &b), 0.0);
    }

    #[test]
    fn partial_overlap_values() {
        // A = {1,2,3}, B = {2,3,4}; intersection 2, union 4
        let a = fp(&[1, 2, 3], 64);
        let b = fp(&[2, 3, 4], 64);
        assert!((tanimoto(&a, &b) - 0.5).abs() < 1e-12);
        // dice = 2*2 / (3+3) = 2/3
        assert!((dice(&a, &b) - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn tversky_recovers_tanimoto_and_dice() {
        let a = fp(&[1, 2, 3], 64);
        let b = fp(&[2, 3, 4], 64);
        assert!((tversky(&a, &b, 1.0, 1.0) - tanimoto(&a, &b)).abs() < 1e-12);
        assert!((tversky(&a, &b, 0.5, 0.5) - dice(&a, &b)).abs() < 1e-12);
    }

    #[test]
    fn empty_fingerprints_identical() {
        let a = FingerprintBits::new(64);
        let b = FingerprintBits::new(64);
        assert_eq!(tanimoto(&a, &b), 1.0);
        assert_eq!(dice(&a, &b), 1.0);
    }

    #[test]
    fn mismatched_length_zero() {
        let a = FingerprintBits::new(64);
        let b = FingerprintBits::new(128);
        assert_eq!(tanimoto(&a, &b), 0.0);
    }

    #[test]
    fn nearest_neighbor_finds_best() {
        let q = fp(&[1, 2, 3], 64);
        let lib = vec![fp(&[10, 11], 64), fp(&[1, 2, 3], 64), fp(&[1, 2], 64)];
        let (idx, sim) = nearest_neighbor(&q, &lib).unwrap();
        assert_eq!(idx, 1);
        assert!((sim - 1.0).abs() < 1e-12);
    }
}
