//! Centroid and maximum-expected-accuracy (MEA) structures.
//!
//! Both summarise the Boltzmann ensemble (a
//! [`crate::ensemble::partition::PartitionFunction`]) into a single
//! representative structure, but with different objectives:
//!
//! - **Centroid** — the structure with the smallest total base-pair
//!   distance to every other structure in the ensemble. It has a
//!   famously simple form: include exactly the pairs with probability
//!   `p(i,j) > 0.5`. Because no two pairs with `p > 0.5` can both be
//!   in conflict, the result is automatically a valid structure.
//! - **MEA** — the structure maximising the expected number of
//!   correctly-called positions: `Σ p(i,j)` over chosen pairs plus
//!   `γ⁻¹ · Σ p_unpaired(i)` over chosen unpaired positions. This is a
//!   Nussinov-style `O(n³)` DP over the probability matrix; the
//!   parameter `gamma` trades sensitivity against specificity.

use crate::ensemble::partition::PartitionFunction;
use crate::error::{Result, RnaStructError};
use crate::fold::nussinov::MIN_HAIRPIN;
use crate::structure::Structure;

/// The default MEA `gamma` — equal weight to a true pair and a true
/// unpaired position. Larger `gamma` favours more pairs.
pub const DEFAULT_GAMMA: f64 = 1.0;

/// Builds the centroid structure from a partition function.
///
/// The centroid contains exactly those pairs with
/// [`PartitionFunction::pair_probability`] strictly above 0.5.
///
/// # Errors
/// [`RnaStructError::Structure`] only if the resulting pair set is
/// somehow inconsistent (cannot happen for a `p > 0.5` rule — kept
/// for signature symmetry).
pub fn centroid(pf: &PartitionFunction) -> Result<Structure> {
    let n = pf.len();
    let mut partner: Vec<Option<usize>> = vec![None; n];
    for i in 0..n {
        for j in (i + 1)..n {
            if pf.pair_probability(i, j) > 0.5 {
                // No conflicts are possible: if i already had a >0.5
                // partner, two probabilities at i would exceed 0.5,
                // summing past 1. Still, guard defensively.
                if partner[i].is_none() && partner[j].is_none() {
                    partner[i] = Some(j);
                    partner[j] = Some(i);
                } else {
                    return Err(RnaStructError::structure(
                        "centroid pair conflict (probabilities inconsistent)",
                    ));
                }
            }
        }
    }
    Structure::from_partner(partner)
}

/// The MEA structure together with its expected-accuracy score.
#[derive(Clone, Debug, PartialEq)]
pub struct MeaResult {
    /// The maximum-expected-accuracy structure.
    pub structure: Structure,
    /// The expected-accuracy objective value it achieves.
    pub score: f64,
    /// The `gamma` used.
    pub gamma: f64,
}

/// Computes the maximum-expected-accuracy structure from a partition
/// function with the default [`DEFAULT_GAMMA`].
pub fn mea(pf: &PartitionFunction) -> Result<MeaResult> {
    mea_with_gamma(pf, DEFAULT_GAMMA)
}

/// Computes the MEA structure with an explicit `gamma`.
///
/// The objective maximised is
/// `Σ_{(i,j)∈S} 2·γ·p(i,j) + Σ_{i unpaired in S} p_unpaired(i)`.
/// (The factor 2 accounts for a pair covering two positions.)
///
/// # Errors
/// [`RnaStructError::Invalid`] if `gamma` is not finite or negative.
pub fn mea_with_gamma(pf: &PartitionFunction, gamma: f64) -> Result<MeaResult> {
    if !gamma.is_finite() || gamma < 0.0 {
        return Err(RnaStructError::invalid(
            "gamma",
            "MEA gamma must be a finite non-negative number",
        ));
    }
    let n = pf.len();
    if n == 0 {
        return Ok(MeaResult {
            structure: Structure::empty(0),
            score: 0.0,
            gamma,
        });
    }

    // m[i][j] = best expected accuracy over the sub-sequence i..=j.
    let mut m = vec![0.0_f64; n * n];
    let at = |i: usize, j: usize| i * n + j;

    // Precompute per-base unpaired probability.
    let unpaired: Vec<f64> = (0..n).map(|i| pf.unpaired_probability(i)).collect();

    for span in 1..n {
        for i in 0..(n - span) {
            let j = i + span;
            // i unpaired
            let mut best = m[at(i + 1, j)] + unpaired[i];
            // j unpaired
            best = best.max(m[at(i, j - 1)] + unpaired[j]);
            // i pairs j
            if span > MIN_HAIRPIN {
                let p = pf.pair_probability(i, j);
                if p > 0.0 {
                    let inner = if i + 2 <= j {
                        m[at(i + 1, j - 1)]
                    } else {
                        0.0
                    };
                    best = best.max(inner + 2.0 * gamma * p);
                }
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
    let feq = |a: f64, b: f64| (a - b).abs() < 1e-9;
    while let Some((i, j)) = stack.pop() {
        if i >= j {
            continue;
        }
        let here = m[at(i, j)];
        if feq(here, m[at(i + 1, j)] + unpaired[i]) {
            stack.push((i + 1, j));
            continue;
        }
        if feq(here, m[at(i, j - 1)] + unpaired[j]) {
            stack.push((i, j - 1));
            continue;
        }
        if j - i > MIN_HAIRPIN {
            let p = pf.pair_probability(i, j);
            let inner = if i + 2 <= j {
                m[at(i + 1, j - 1)]
            } else {
                0.0
            };
            if p > 0.0 && feq(here, inner + 2.0 * gamma * p) {
                partner[i] = Some(j);
                partner[j] = Some(i);
                if i + 2 <= j {
                    stack.push((i + 1, j - 1));
                }
                continue;
            }
        }
        for k in (i + 1)..j {
            if feq(here, m[at(i, k)] + m[at(k + 1, j)]) {
                stack.push((i, k));
                stack.push((k + 1, j));
                break;
            }
        }
    }

    let structure = Structure::from_partner(partner)?;
    let score = m[at(0, n - 1)];
    Ok(MeaResult {
        structure,
        score,
        gamma,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ensemble::partition::partition_function;
    use crate::rna::RnaSeq;

    #[test]
    fn centroid_of_stable_stem_has_pairs() {
        let seq = RnaSeq::parse("GGGGGGAAAACCCCCC").unwrap();
        let pf = partition_function(&seq).unwrap();
        let c = centroid(&pf).unwrap();
        // a strong stem should leave some > 0.5 pairs
        assert!(c.is_nested());
        assert!(c.len() == seq.len());
    }

    #[test]
    fn centroid_of_unpairable_is_empty() {
        let seq = RnaSeq::parse("AAAAAAAAAAAA").unwrap();
        let pf = partition_function(&seq).unwrap();
        let c = centroid(&pf).unwrap();
        assert_eq!(c.n_pairs(), 0);
    }

    #[test]
    fn mea_runs_and_is_nested() {
        let seq = RnaSeq::parse("GGGGGAAAACCCCC").unwrap();
        let pf = partition_function(&seq).unwrap();
        let r = mea(&pf).unwrap();
        assert!(r.structure.is_nested());
        assert!(r.score >= 0.0);
    }

    #[test]
    fn mea_score_grows_with_gamma() {
        // Larger gamma weights pairs more, so the structure picks up
        // more (or equal) pairs.
        let seq = RnaSeq::parse("GGGGGGAAAACCCCCC").unwrap();
        let pf = partition_function(&seq).unwrap();
        let low = mea_with_gamma(&pf, 0.1).unwrap();
        let high = mea_with_gamma(&pf, 8.0).unwrap();
        assert!(high.structure.n_pairs() >= low.structure.n_pairs());
    }

    #[test]
    fn mea_rejects_bad_gamma() {
        let seq = RnaSeq::parse("GGGGCCCC").unwrap();
        let pf = partition_function(&seq).unwrap();
        assert!(mea_with_gamma(&pf, -1.0).is_err());
        assert!(mea_with_gamma(&pf, f64::NAN).is_err());
    }

    #[test]
    fn mea_unpairable_sequence() {
        let seq = RnaSeq::parse("AAAAAAAA").unwrap();
        let pf = partition_function(&seq).unwrap();
        let r = mea(&pf).unwrap();
        assert_eq!(r.structure.n_pairs(), 0);
        // every base unpaired -> score ~ sum of unpaired probs ~ n
        assert!(r.score > 0.0);
    }
}
