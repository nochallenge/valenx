//! Descriptor matching (Stage 2, part 1).
//!
//! Given two sets of binary descriptors — typically the `[u8; 32]`
//! signatures produced for two images by [`crate::detect_and_describe`] —
//! this module finds putative feature correspondences between them.
//!
//! The pipeline is the standard one from the SfM / wide-baseline-stereo
//! literature:
//!
//! 1. **Brute-force nearest neighbour.** For every *query* descriptor, scan
//!    *all* *train* descriptors and keep the two with the smallest
//!    [`crate::hamming_distance`] (the nearest and the second-nearest).
//!
//! 2. **Lowe ratio test** (Lowe 2004). A match is only kept when the best
//!    distance is clearly better than the second-best:
//!    `best_dist < ratio * second_dist`. This rejects ambiguous matches
//!    where two train descriptors are almost equally close — the hallmark
//!    of a repeated texture or a non-discriminative point. A typical
//!    `ratio` is `0.7`–`0.8`.
//!
//! 3. **Mutual (cross-check) consistency.** The surviving query→train pair
//!    is only kept if the train descriptor's own nearest neighbour back in
//!    the query set is exactly that query descriptor. This symmetric check
//!    removes one-sided matches and is a cheap, effective outlier filter.
//!
//! ## Cost (stated honestly)
//!
//! This is an exhaustive **brute-force** matcher: with `n` query and `m`
//! train descriptors it performs `Θ(n · m)` Hamming comparisons for the
//! forward pass, plus another `Θ(m · n)` for the cross-check direction —
//! i.e. `Θ(n · m)` overall, with no spatial or LSH acceleration. Each
//! comparison is 32 byte-XOR + `count_ones` operations (the 256-bit
//! Hamming distance). For the few-thousand-keypoint image pairs typical of
//! incremental SfM this is entirely adequate and gives *exact* nearest
//! neighbours; for very large descriptor sets a multi-probe LSH or a
//! vocabulary tree would be the next optimisation, and is intentionally
//! out of scope for this stage.

use crate::descriptor::{hamming_distance, DESCRIPTOR_BYTES};

/// A putative correspondence between a query descriptor and a train
/// descriptor, with the Hamming distance between them.
///
/// `query_idx` indexes the `query` slice passed to [`match_descriptors`]
/// and `train_idx` indexes the `train` slice; `distance` is the 256-bit
/// [`crate::hamming_distance`] of the two descriptors (range `0..=256`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Match {
    /// Index of the matched descriptor in the `query` slice.
    pub query_idx: usize,
    /// Index of the matched descriptor in the `train` slice.
    pub train_idx: usize,
    /// Hamming distance between the two descriptors (`0..=256`).
    pub distance: u32,
}

/// The nearest and second-nearest train descriptors for one query
/// descriptor, by Hamming distance.
#[derive(Clone, Copy)]
struct TwoNearest {
    /// Index of the closest train descriptor, if any train descriptors
    /// exist.
    best_idx: Option<usize>,
    /// Distance to the closest train descriptor.
    best_dist: u32,
    /// Distance to the second-closest train descriptor (`u32::MAX` if there
    /// are fewer than two train descriptors).
    second_dist: u32,
}

/// Find the nearest and second-nearest descriptor in `train` to
/// `query[qi]`, by Hamming distance, via an exhaustive scan.
fn two_nearest(
    query: &[[u8; DESCRIPTOR_BYTES]],
    qi: usize,
    train: &[[u8; DESCRIPTOR_BYTES]],
) -> TwoNearest {
    let mut best_idx: Option<usize> = None;
    let mut best_dist = u32::MAX;
    let mut second_dist = u32::MAX;
    let q = &query[qi];
    for (ti, t) in train.iter().enumerate() {
        let d = hamming_distance(q, t);
        if d < best_dist {
            second_dist = best_dist;
            best_dist = d;
            best_idx = Some(ti);
        } else if d < second_dist {
            second_dist = d;
        }
    }
    TwoNearest {
        best_idx,
        best_dist,
        second_dist,
    }
}

/// Match two sets of binary descriptors with the Lowe ratio test and a
/// mutual cross-check.
///
/// For every descriptor in `query`, the nearest and second-nearest
/// descriptors in `train` are found by [`crate::hamming_distance`]
/// (brute-force). The match is kept only if **both**:
///
/// - the **Lowe ratio test** passes — `best_dist < ratio * second_dist`,
///   so the best match is unambiguously closer than the runner-up; and
/// - the **cross-check** passes — `train[train_idx]`'s own nearest
///   neighbour, searched back over `query`, is exactly `query[query_idx]`.
///
/// The returned [`Match`]es are unique in both `query_idx` and `train_idx`
/// (a consequence of the mutual nearest-neighbour requirement) and are
/// listed in increasing `query_idx` order.
///
/// `ratio` is the Lowe threshold; values `≤ 0.0` reject everything and
/// values `≥ 1.0` disable the ratio test (cross-check still applies). A
/// `ratio` between `0.7` and `0.8` is the usual choice. When a query has
/// only one candidate in `train` there is no second-nearest, so the ratio
/// test cannot discriminate and the match is rejected (conservative).
///
/// # Cost
///
/// `Θ(n · m)` Hamming comparisons for `n = query.len()`,
/// `m = train.len()` — see the [module docs](self) for the full
/// discussion. Returns an empty vector if either slice is empty.
#[must_use]
pub fn match_descriptors(
    query: &[[u8; DESCRIPTOR_BYTES]],
    train: &[[u8; DESCRIPTOR_BYTES]],
    ratio: f32,
) -> Vec<Match> {
    if query.is_empty() || train.is_empty() {
        return Vec::new();
    }

    // Forward pass: for each train descriptor, remember which query
    // descriptor is its nearest neighbour, so the cross-check is O(1) per
    // candidate rather than another full scan inside the loop.
    let mut train_best_query: Vec<Option<usize>> = vec![None; train.len()];
    {
        let mut train_best_dist: Vec<u32> = vec![u32::MAX; train.len()];
        for (qi, q) in query.iter().enumerate() {
            for (ti, t) in train.iter().enumerate() {
                let d = hamming_distance(q, t);
                if d < train_best_dist[ti] {
                    train_best_dist[ti] = d;
                    train_best_query[ti] = Some(qi);
                }
            }
        }
    }

    let mut matches = Vec::new();
    for qi in 0..query.len() {
        let nn = two_nearest(query, qi, train);
        let Some(ti) = nn.best_idx else {
            continue;
        };

        // Lowe ratio test. Requires a genuine second-nearest neighbour;
        // with only one candidate the test is undefined, so we reject.
        // best_dist < ratio * second_dist, evaluated in f32: both distances
        // are integers in 0..=256 (exactly representable), so the only
        // rounding is in the intended `ratio * second_dist` scaling itself.
        if nn.second_dist == u32::MAX {
            continue;
        }
        if (nn.best_dist as f32) >= ratio * (nn.second_dist as f32) {
            continue;
        }

        // Mutual cross-check: train[ti]'s nearest query must be qi.
        if train_best_query[ti] != Some(qi) {
            continue;
        }

        matches.push(Match {
            query_idx: qi,
            train_idx: ti,
            distance: nn.best_dist,
        });
    }
    matches
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a deterministic, well-separated set of 256-bit descriptors.
    /// Each descriptor `k` is a distinct bit pattern derived from `k`, with
    /// enough spread that the nearest neighbour of a descriptor is itself.
    fn distinct_descriptors(n: usize) -> Vec<[u8; DESCRIPTOR_BYTES]> {
        (0..n)
            .map(|k| {
                let mut d = [0u8; DESCRIPTOR_BYTES];
                // Spread bits across all 32 bytes so distinct k values are
                // far apart in Hamming space: byte j gets a function of
                // both k and j.
                for (j, b) in d.iter_mut().enumerate() {
                    let v = (k.wrapping_mul(31).wrapping_add(j.wrapping_mul(97))) as u8;
                    *b = v ^ ((k >> (j % 8)) as u8);
                }
                d
            })
            .collect()
    }

    /// Flip `bits` low bits of a descriptor to simulate a noisy but still
    /// nearest observation of the same point.
    fn perturb(mut d: [u8; DESCRIPTOR_BYTES], bits: usize) -> [u8; DESCRIPTOR_BYTES] {
        for i in 0..bits {
            d[i >> 3] ^= 1u8 << (i & 7);
        }
        d
    }

    // Test 1: matching recovers a known permutation.
    #[test]
    fn recovers_known_permutation() {
        let a = distinct_descriptors(20);
        // A fixed, non-trivial permutation of the indices.
        let perm: Vec<usize> = (0..20).map(|i| (i * 7 + 3) % 20).collect();
        // Build B so that B[i] == A[perm[i]]; then for query=A, the match
        // for A[k] is the train index i with perm[i] == k.
        let b: Vec<[u8; DESCRIPTOR_BYTES]> = perm.iter().map(|&p| a[p]).collect();

        let matches = match_descriptors(&a, &b, 0.8);
        assert_eq!(
            matches.len(),
            20,
            "every descriptor should match exactly once"
        );

        for m in &matches {
            // A[query_idx] must equal B[train_idx], i.e. perm[train_idx] == query_idx.
            assert_eq!(
                perm[m.train_idx], m.query_idx,
                "query {} should map to the train slot holding it",
                m.query_idx
            );
            assert_eq!(m.distance, 0, "identical descriptors are an exact match");
        }

        // train indices are a bijection (each used once).
        let mut seen = [false; 20];
        for m in &matches {
            assert!(!seen[m.train_idx], "train index reused");
            seen[m.train_idx] = true;
        }
    }

    // Test 2a: the ratio test rejects ambiguous matches.
    #[test]
    fn ratio_test_rejects_ambiguous() {
        // One query descriptor, two near-identical train descriptors: the
        // best and second-best distances are close, so the ratio test must
        // reject (ambiguous repeated texture).
        let base = distinct_descriptors(1)[0];
        let query = vec![base];
        let train = vec![perturb(base, 2), perturb(base, 3)];

        // Distances are 2 and 3: 2 >= 0.6 * 3 = 1.8 is true, so rejected.
        let strict = match_descriptors(&query, &train, 0.6);
        assert!(
            strict.is_empty(),
            "ambiguous pair must be rejected by a strict ratio"
        );

        // A permissive ratio (1.0 disables the ratio gate) lets the nearest
        // through (cross-check still holds: query is the train's nearest too).
        let loose = match_descriptors(&query, &train, 1.0);
        assert_eq!(loose.len(), 1);
        assert_eq!(
            loose[0].train_idx, 0,
            "the closer (2-bit) train descriptor wins"
        );
    }

    // Test 2b: the cross-check rejects non-mutual matches.
    #[test]
    fn cross_check_rejects_non_mutual() {
        // q0 and q1 are both closest to t0; t0 is closest to q0. So q1->t0
        // is a one-sided (non-mutual) match and must be dropped, while
        // q0->t0 survives.
        let descs = distinct_descriptors(3);
        let t0 = descs[0];
        let q0 = perturb(t0, 1); // 1 bit from t0
        let q1 = perturb(t0, 5); // 5 bits from t0, but still closest to t0
                                 // t0's nearest query is q0 (closer), so q1->t0 is one-sided; t1 is
                                 // an unrelated/far descriptor.
        let t1 = descs[1];
        let query = vec![q0, q1];
        let train = vec![t0, t1];

        // Use a permissive ratio so only the cross-check does the filtering.
        let matches = match_descriptors(&query, &train, 1.0);
        // q1->t0 is non-mutual (t0's nearest is q0), so it is rejected.
        assert!(
            matches
                .iter()
                .all(|m| !(m.query_idx == 1 && m.train_idx == 0)),
            "non-mutual match q1->t0 must be rejected by cross-check"
        );
        // q0->t0 is mutual and should be present.
        assert!(
            matches.iter().any(|m| m.query_idx == 0 && m.train_idx == 0),
            "mutual match q0->t0 should survive"
        );
    }

    #[test]
    fn empty_inputs_yield_no_matches() {
        let a = distinct_descriptors(3);
        assert!(match_descriptors(&[], &a, 0.8).is_empty());
        assert!(match_descriptors(&a, &[], 0.8).is_empty());
        assert!(match_descriptors(&[], &[], 0.8).is_empty());
    }
}
