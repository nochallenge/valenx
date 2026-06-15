//! Fuse heterogeneous evidence channels into one comparable, ranked score.

use serde::{Deserialize, Serialize};

use crate::error::ScoreError;

/// Soft scale (kcal/mol) for mapping an energy / dock score to `[0, 1]`. A
/// convenience default — the right answer is to *calibrate* against data, not
/// to trust this transform as a probability.
pub const AFFINITY_SCALE: f64 = 5.0;

/// The raw evidence available for a candidate. Any channel may be absent.
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct ScoreComponents {
    /// Interface predicted TM-score, in `[0, 1]` (higher is better).
    pub iptm: Option<f64>,
    /// Mean pLDDT confidence, rescaled to `[0, 1]` (higher is better).
    pub plddt: Option<f64>,
    /// Docking score in kcal/mol (more negative is better).
    pub dock_score: Option<f64>,
    /// Endpoint binding-energy estimate in kcal/mol (more negative is better).
    pub dg_bind_kcal: Option<f64>,
}

/// Per-channel weights for the fused score (default `1.0` each). A weight of
/// `0` drops the channel.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ScoreWeights {
    /// Weight on ipTM.
    pub iptm: f64,
    /// Weight on pLDDT.
    pub plddt: f64,
    /// Weight on the docking score.
    pub dock: f64,
    /// Weight on the binding-energy estimate.
    pub dg: f64,
}

impl Default for ScoreWeights {
    fn default() -> Self {
        Self {
            iptm: 1.0,
            plddt: 1.0,
            dock: 1.0,
            dg: 1.0,
        }
    }
}

/// One normalized channel inside a [`ComparableScore`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScoreChannel {
    /// Channel name.
    pub name: String,
    /// The raw input value.
    pub raw: f64,
    /// The value mapped to `[0, 1]` (higher is better).
    pub normalized: f64,
    /// The weight applied.
    pub weight: f64,
}

/// A single comparable score in `[0, 1]` with its per-channel breakdown.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ComparableScore {
    /// The fused, weighted-mean score in `[0, 1]` (higher is better).
    pub value: f64,
    /// The contributing channels, normalized and weighted.
    pub breakdown: Vec<ScoreChannel>,
}

/// Map an energy / dock score (more-negative-is-better) to `[0, 1]` via
/// `1 / (1 + exp(x / scale))`: very negative → 1, zero → 0.5, positive → 0.
fn normalize_affinity(x: f64, scale: f64) -> f64 {
    1.0 / (1.0 + (x / scale).exp())
}

fn check_confidence(what: &'static str, v: f64) -> Result<(), ScoreError> {
    if !v.is_finite() {
        return Err(ScoreError::NonFinite { what });
    }
    if !(0.0..=1.0).contains(&v) {
        return Err(ScoreError::ConfidenceOutOfRange { what, value: v });
    }
    Ok(())
}

impl ComparableScore {
    /// Compute the comparable score from the present channels and `weights`.
    ///
    /// Confidences pass through; the docking score and ΔG are mapped through
    /// a `more-negative-is-better` transform with [`AFFINITY_SCALE`]. The value is the
    /// weighted mean of the present, positively-weighted channels. Errors if no
    /// channel contributes, or a confidence is out of `[0, 1]`.
    pub fn compute(
        components: &ScoreComponents,
        weights: &ScoreWeights,
    ) -> Result<Self, ScoreError> {
        let mut breakdown: Vec<ScoreChannel> = Vec::new();
        let mut push = |name: &str, raw: f64, normalized: f64, weight: f64| {
            if weight > 0.0 {
                breakdown.push(ScoreChannel {
                    name: name.to_string(),
                    raw,
                    normalized,
                    weight,
                });
            }
        };

        if let Some(v) = components.iptm {
            check_confidence("iptm", v)?;
            push("iptm", v, v, weights.iptm);
        }
        if let Some(v) = components.plddt {
            check_confidence("plddt", v)?;
            push("plddt", v, v, weights.plddt);
        }
        if let Some(v) = components.dock_score {
            if !v.is_finite() {
                return Err(ScoreError::NonFinite { what: "dock_score" });
            }
            push(
                "dock_score",
                v,
                normalize_affinity(v, AFFINITY_SCALE),
                weights.dock,
            );
        }
        if let Some(v) = components.dg_bind_kcal {
            if !v.is_finite() {
                return Err(ScoreError::NonFinite {
                    what: "dg_bind_kcal",
                });
            }
            push(
                "dg_bind",
                v,
                normalize_affinity(v, AFFINITY_SCALE),
                weights.dg,
            );
        }

        let wsum: f64 = breakdown.iter().map(|c| c.weight).sum();
        if breakdown.is_empty() || wsum <= 0.0 {
            return Err(ScoreError::NoComponents);
        }
        let value = breakdown
            .iter()
            .map(|c| c.normalized * c.weight)
            .sum::<f64>()
            / wsum;
        Ok(Self { value, breakdown })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_best_channels_score_near_one() {
        let c = ScoreComponents {
            iptm: Some(1.0),
            plddt: Some(1.0),
            dock_score: Some(-25.0),
            dg_bind_kcal: Some(-25.0),
        };
        let s = ComparableScore::compute(&c, &ScoreWeights::default()).unwrap();
        assert!(s.value > 0.99, "value = {}", s.value);
        assert_eq!(s.breakdown.len(), 4);
    }

    #[test]
    fn more_negative_dg_scores_higher() {
        let weak = ScoreComponents {
            dg_bind_kcal: Some(-1.0),
            ..Default::default()
        };
        let strong = ScoreComponents {
            dg_bind_kcal: Some(-15.0),
            ..Default::default()
        };
        let w = ScoreWeights::default();
        let a = ComparableScore::compute(&weak, &w).unwrap().value;
        let b = ComparableScore::compute(&strong, &w).unwrap().value;
        assert!(b > a);
    }

    #[test]
    fn absent_channels_are_excluded() {
        let c = ScoreComponents {
            iptm: Some(0.7),
            ..Default::default()
        };
        let s = ComparableScore::compute(&c, &ScoreWeights::default()).unwrap();
        assert_eq!(s.breakdown.len(), 1);
        assert!((s.value - 0.7).abs() < 1e-12);
    }

    #[test]
    fn zero_weight_drops_a_channel() {
        let c = ScoreComponents {
            iptm: Some(0.9),
            plddt: Some(0.1),
            ..Default::default()
        };
        let w = ScoreWeights {
            plddt: 0.0,
            ..Default::default()
        };
        let s = ComparableScore::compute(&c, &w).unwrap();
        assert_eq!(s.breakdown.len(), 1);
        assert!((s.value - 0.9).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_confidence_and_empty() {
        let bad = ScoreComponents {
            iptm: Some(1.5),
            ..Default::default()
        };
        assert_eq!(
            ComparableScore::compute(&bad, &ScoreWeights::default())
                .unwrap_err()
                .code(),
            "confidence_out_of_range"
        );
        assert_eq!(
            ComparableScore::compute(&ScoreComponents::default(), &ScoreWeights::default())
                .unwrap_err()
                .code(),
            "no_components"
        );
    }

    #[test]
    fn serde_round_trips() {
        let c = ScoreComponents {
            iptm: Some(0.8),
            dg_bind_kcal: Some(-10.0),
            ..Default::default()
        };
        let s = ComparableScore::compute(&c, &ScoreWeights::default()).unwrap();
        let j = serde_json::to_string(&s).unwrap();
        let back: ComparableScore = serde_json::from_str(&j).unwrap();
        assert_eq!(s, back);
    }
}
