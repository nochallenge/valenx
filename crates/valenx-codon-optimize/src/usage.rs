//! Host codon-usage weights (relative adaptiveness) for CAI and optimization.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::code::{synonymous_codons, translate_codon};
use crate::error::CodonError;

/// Floor applied to a zero/near-zero relative adaptiveness so weights stay in
/// `(0, 1]` (a zero weight would make CAI's log undefined). A common convention.
const WEIGHT_FLOOR: f64 = 1e-3;

/// Per-codon relative-adaptiveness weights `w(codon) ∈ (0, 1]` for a host.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CodonUsage {
    weights: HashMap<String, f64>,
}

fn valid_sense_codon(codon: &str) -> bool {
    matches!(translate_codon(codon), Some(aa) if aa != '*')
}

impl CodonUsage {
    /// Build from precomputed relative-adaptiveness weights. Every key must be a
    /// valid sense codon and every value must lie in `(0, 1]`.
    pub fn from_weights(weights: HashMap<String, f64>) -> Result<Self, CodonError> {
        if weights.is_empty() {
            return Err(CodonError::Empty { what: "weights" });
        }
        for (codon, &w) in &weights {
            if !valid_sense_codon(codon) {
                return Err(CodonError::InvalidCodon {
                    codon: codon.clone(),
                    index: 0,
                });
            }
            if !w.is_finite() || w <= 0.0 || w > 1.0 {
                return Err(CodonError::WeightOutOfRange {
                    codon: codon.clone(),
                    value: w,
                });
            }
        }
        let weights = weights
            .into_iter()
            .map(|(c, w)| (c.to_ascii_uppercase(), w))
            .collect();
        Ok(Self { weights })
    }

    /// Build from raw codon frequencies/counts: `w(codon) = f(codon) / max f`
    /// over its synonymous group, floored to `(0, 1]`.
    pub fn from_frequencies(freqs: &HashMap<String, f64>) -> Result<Self, CodonError> {
        if freqs.is_empty() {
            return Err(CodonError::Empty {
                what: "frequencies",
            });
        }
        let mut max_for_aa: HashMap<char, f64> = HashMap::new();
        for (codon, &f) in freqs {
            let aa = translate_codon(codon).ok_or_else(|| CodonError::InvalidCodon {
                codon: codon.clone(),
                index: 0,
            })?;
            if aa == '*' {
                continue;
            }
            if !f.is_finite() || f < 0.0 {
                return Err(CodonError::WeightOutOfRange {
                    codon: codon.clone(),
                    value: f,
                });
            }
            max_for_aa
                .entry(aa)
                .and_modify(|m| {
                    if f > *m {
                        *m = f;
                    }
                })
                .or_insert(f);
        }
        let mut weights = HashMap::new();
        for (codon, &f) in freqs {
            let aa = translate_codon(codon).unwrap();
            if aa == '*' {
                continue;
            }
            let m = max_for_aa[&aa];
            let w = if m > 0.0 {
                (f / m).max(WEIGHT_FLOOR)
            } else {
                WEIGHT_FLOOR
            };
            weights.insert(codon.to_ascii_uppercase(), w.min(1.0));
        }
        Self::from_weights(weights)
    }

    /// The relative adaptiveness of `codon`, if known.
    pub fn weight(&self, codon: &str) -> Option<f64> {
        self.weights.get(&codon.to_ascii_uppercase()).copied()
    }

    /// The highest-adaptiveness synonymous codon for `aa`, or `None` if `aa` is
    /// not a sense residue or none of its codons have a weight.
    pub fn optimal_codon(&self, aa: char) -> Option<String> {
        let mut best: Option<(String, f64)> = None;
        for codon in synonymous_codons(aa) {
            if let Some(w) = self.weight(&codon) {
                if best.as_ref().is_none_or(|(_, bw)| w > *bw) {
                    best = Some((codon, w));
                }
            }
        }
        best.map(|(c, _)| c)
    }
}

/// An **illustrative** weight set: the first synonymous codon of each amino acid
/// gets weight `1.0`, the rest `0.5`. This is a teaching default, **not** a real
/// organism's measured codon-usage table — supply your own (e.g. from Kazusa)
/// for real work.
pub fn illustrative_weights() -> CodonUsage {
    let mut weights = HashMap::new();
    for aa in "ACDEFGHIKLMNPQRSTVWY".chars() {
        for (i, codon) in synonymous_codons(aa).into_iter().enumerate() {
            weights.insert(codon, if i == 0 { 1.0 } else { 0.5 });
        }
    }
    CodonUsage::from_weights(weights).expect("illustrative weights are valid")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn illustrative_has_a_unit_optimal_per_aa() {
        let u = illustrative_weights();
        for aa in "ACDEFGHIKLMNPQRSTVWY".chars() {
            let opt = u.optimal_codon(aa).unwrap();
            assert!((u.weight(&opt).unwrap() - 1.0).abs() < 1e-12, "aa {aa}");
        }
    }

    #[test]
    fn single_codon_aa_optimal_is_the_only_codon() {
        let u = illustrative_weights();
        assert_eq!(u.optimal_codon('M').unwrap(), "ATG");
        assert_eq!(u.optimal_codon('W').unwrap(), "TGG");
    }

    #[test]
    fn from_frequencies_normalizes_to_group_max() {
        let mut f = HashMap::new();
        f.insert("TTT".to_string(), 30.0); // Phe
        f.insert("TTC".to_string(), 10.0); // Phe -> w = 10/30
        let u = CodonUsage::from_frequencies(&f).unwrap();
        assert!((u.weight("TTT").unwrap() - 1.0).abs() < 1e-12);
        assert!((u.weight("TTC").unwrap() - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_weights() {
        let mut w = HashMap::new();
        w.insert("TTT".to_string(), 1.5);
        assert_eq!(
            CodonUsage::from_weights(w).unwrap_err().code(),
            "weight_out_of_range"
        );
        let mut w2 = HashMap::new();
        w2.insert("ZZZ".to_string(), 0.5);
        assert_eq!(
            CodonUsage::from_weights(w2).unwrap_err().code(),
            "invalid_codon"
        );
    }
}
