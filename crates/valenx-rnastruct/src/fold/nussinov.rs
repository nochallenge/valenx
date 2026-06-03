//! Nussinov maximum-base-pairing folding.
//!
//! The Nussinov-Jacobson (1978) algorithm is the simplest secondary-
//! structure DP: it finds the nested structure that *maximises the
//! number of base pairs*, with no energy model at all. It is `O(n³)`
//! time / `O(n²)` space and is exact for that objective.
//!
//! It is included as a baseline / teaching folder and because the
//! maximum-pairing structure is a useful sanity reference. For
//! thermodynamically meaningful structures use [`crate::fold::zuker`].
//!
//! The recurrence over the score matrix `m[i][j]` (best pairs in the
//! sub-sequence `i..=j`) is:
//!
//! ```text
//! m[i][j] = max(
//!     m[i+1][j],                                  // i unpaired
//!     m[i][j-1],                                  // j unpaired
//!     m[i+1][j-1] + pair(i,j),                    // i pairs j
//!     max over i<k<j of  m[i][k] + m[k+1][j],     // bifurcation
//! )
//! ```

use crate::error::Result;
use crate::fold::energy::can_pair_codes;
use crate::rna::RnaSeq;
use crate::structure::Structure;

/// The minimum number of unpaired bases a hairpin loop must enclose.
/// A pair `(i, j)` is only allowed when `j - i - 1 >= MIN_HAIRPIN`.
pub const MIN_HAIRPIN: usize = 3;

/// Result of a Nussinov fold: the maximum-pairing [`Structure`] and the
/// pair count it achieves.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NussinovResult {
    /// The maximum-base-pairing nested structure.
    pub structure: Structure,
    /// The number of base pairs (the optimised objective).
    pub pairs: usize,
}

/// Folds `seq` to the nested structure with the most base pairs.
///
/// Only canonical / wobble pairs separated by at least
/// [`MIN_HAIRPIN`] unpaired bases are allowed.
///
/// # Errors
/// Never fails for a valid [`RnaSeq`]; the `Result` is kept for
/// signature symmetry with the energy folders.
pub fn fold(seq: &RnaSeq) -> Result<NussinovResult> {
    let codes = seq.codes();
    let n = codes.len();
    if n == 0 {
        return Ok(NussinovResult {
            structure: Structure::empty(0),
            pairs: 0,
        });
    }

    // m[i][j] = best #pairs in i..=j. Use a flat n*n buffer.
    let mut m = vec![0u32; n * n];
    let at = |i: usize, j: usize| i * n + j;

    // span = j - i, growing from MIN_HAIRPIN+1 (shortest pair-able
    // window) to n-1.
    for span in 1..n {
        for i in 0..(n - span) {
            let j = i + span;
            // i unpaired, or j unpaired
            let mut best = m[at(i + 1, j)].max(m[at(i, j - 1)]);
            // i pairs with j
            if j - i > MIN_HAIRPIN && can_pair_codes(codes[i], codes[j]) {
                let inner = if i + 2 <= j {
                    m[at(i + 1, j - 1)]
                } else {
                    0
                };
                best = best.max(inner + 1);
            }
            // bifurcation
            for k in (i + 1)..j {
                let cand = m[at(i, k)] + m[at(k + 1, j)];
                if cand > best {
                    best = cand;
                }
            }
            m[at(i, j)] = best;
        }
    }

    // Traceback.
    let mut partner: Vec<Option<usize>> = vec![None; n];
    let mut stack: Vec<(usize, usize)> = vec![(0, n - 1)];
    while let Some((i, j)) = stack.pop() {
        if i >= j {
            continue;
        }
        let here = m[at(i, j)];
        // case: i unpaired
        if m[at(i + 1, j)] == here {
            stack.push((i + 1, j));
            continue;
        }
        // case: j unpaired
        if m[at(i, j - 1)] == here {
            stack.push((i, j - 1));
            continue;
        }
        // case: i pairs j
        if j - i > MIN_HAIRPIN && can_pair_codes(codes[i], codes[j]) {
            let inner = if i + 2 <= j {
                m[at(i + 1, j - 1)]
            } else {
                0
            };
            if inner + 1 == here {
                partner[i] = Some(j);
                partner[j] = Some(i);
                if i + 2 <= j {
                    stack.push((i + 1, j - 1));
                }
                continue;
            }
        }
        // case: bifurcation
        for k in (i + 1)..j {
            if m[at(i, k)] + m[at(k + 1, j)] == here {
                stack.push((i, k));
                stack.push((k + 1, j));
                break;
            }
        }
    }

    let structure = Structure::from_partner(partner)?;
    let pairs = structure.n_pairs();
    Ok(NussinovResult { structure, pairs })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_sequence() {
        let r = fold(&RnaSeq::parse("A").unwrap()).unwrap();
        assert_eq!(r.pairs, 0);
    }

    #[test]
    fn simple_hairpin_pairs_maximally() {
        // GGGG....CCCC: 4 G-C pairs are possible
        let r = fold(&RnaSeq::parse("GGGGAAAACCCC").unwrap()).unwrap();
        assert_eq!(r.pairs, 4);
        assert!(r.structure.is_nested());
        // outermost pair joins position 0 and 11
        assert_eq!(r.structure.partner(0), Some(11));
    }

    #[test]
    fn respects_min_hairpin() {
        // GGCC has a pairable G-C at (0,3) but only 2 unpaired between
        // -> still allowed (2 >= 3 is false) so 0 pairs
        let r = fold(&RnaSeq::parse("GCGC").unwrap()).unwrap();
        // (0,3): j-i-1 = 2 < 3 -> forbidden; (1,2): j-i-1=0 -> forbidden
        assert_eq!(r.pairs, 0);
    }

    #[test]
    fn no_pairs_for_unpairable() {
        let r = fold(&RnaSeq::parse("AAAAAAAA").unwrap()).unwrap();
        assert_eq!(r.pairs, 0);
    }

    #[test]
    fn structure_is_always_valid() {
        let r = fold(&RnaSeq::parse("GGGCAUGCCCAAAGGGCAUGCCC").unwrap()).unwrap();
        // a Nussinov structure must always be nested
        assert!(r.structure.is_nested());
        assert_eq!(r.structure.n_pairs(), r.pairs);
    }
}
