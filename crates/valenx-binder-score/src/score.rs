//! The binder-quality scoring function.

use serde::{Deserialize, Serialize};

use crate::error::BinderError;

/// Soft scale (kcal/mol) for the binding-ΔG → `[0, 1]` transform. A convenience
/// default; calibrate against data rather than trusting it as a probability.
pub const AFFINITY_SCALE: f64 = 5.0;

/// The evidence for one candidate binder. Any channel may be absent.
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct BinderInputs {
    /// Binding free-energy estimate in kcal/mol (more negative is better).
    pub dg_bind_kcal: Option<f64>,
    /// Developability summary in `[0, 1]` (higher is better).
    pub developability: Option<f64>,
    /// Worst safety severity, `0..=4` (0 = no flags, 4 = critical).
    pub safety_severity: Option<u8>,
}

/// Per-channel weights (default `1.0`; a weight of `0` drops the channel).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BinderWeights {
    /// Weight on binding ΔG.
    pub dg: f64,
    /// Weight on developability.
    pub developability: f64,
    /// Weight on safety.
    pub safety: f64,
}

impl Default for BinderWeights {
    fn default() -> Self {
        Self {
            dg: 1.0,
            developability: 1.0,
            safety: 1.0,
        }
    }
}

/// A fused binder-quality score with its per-channel breakdown.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BinderScore {
    /// The fused score in `[0, 1]` (higher is better).
    pub value: f64,
    /// `(channel, normalized)` for each contributing channel.
    pub components: Vec<(String, f64)>,
    /// Always `true`: a score ranks, it never clears a candidate for the lab.
    pub requires_review: bool,
}

fn normalize_affinity(dg: f64) -> f64 {
    1.0 / (1.0 + (dg / AFFINITY_SCALE).exp())
}

fn check_weight(what: &'static str, w: f64) -> Result<(), BinderError> {
    if !w.is_finite() || w < 0.0 {
        return Err(BinderError::BadWeight { what, value: w });
    }
    Ok(())
}

/// Fold the present channels into a single `[0, 1]` binder-quality score.
pub fn score(inputs: &BinderInputs, weights: &BinderWeights) -> Result<BinderScore, BinderError> {
    check_weight("dg", weights.dg)?;
    check_weight("developability", weights.developability)?;
    check_weight("safety", weights.safety)?;

    let mut components: Vec<(String, f64)> = Vec::new();
    let mut wsum = 0.0;
    let mut acc = 0.0;
    let mut add = |name: &str, normalized: f64, weight: f64| {
        if weight > 0.0 {
            components.push((name.to_string(), normalized));
            acc += normalized * weight;
            wsum += weight;
        }
    };

    if let Some(dg) = inputs.dg_bind_kcal {
        if !dg.is_finite() {
            return Err(BinderError::NonFinite {
                what: "dg_bind_kcal",
            });
        }
        add("binding", normalize_affinity(dg), weights.dg);
    }
    if let Some(d) = inputs.developability {
        if !d.is_finite() {
            return Err(BinderError::NonFinite {
                what: "developability",
            });
        }
        if !(0.0..=1.0).contains(&d) {
            return Err(BinderError::DevelopabilityOutOfRange { value: d });
        }
        add("developability", d, weights.developability);
    }
    if let Some(sev) = inputs.safety_severity {
        if sev > 4 {
            return Err(BinderError::SeverityOutOfRange { value: sev });
        }
        add("safety", 1.0 - f64::from(sev) / 4.0, weights.safety);
    }

    if components.is_empty() || wsum <= 0.0 {
        return Err(BinderError::NoComponents);
    }
    Ok(BinderScore {
        value: acc / wsum,
        components,
        requires_review: true,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strong_candidate_scores_high_and_flags_review() {
        let s = score(
            &BinderInputs {
                dg_bind_kcal: Some(-20.0),
                developability: Some(1.0),
                safety_severity: Some(0),
            },
            &BinderWeights::default(),
        )
        .unwrap();
        assert!(s.value > 0.95, "value = {}", s.value);
        assert_eq!(s.components.len(), 3);
        assert!(s.requires_review); // always
    }

    #[test]
    fn monotone_in_binding_and_severity() {
        let w = BinderWeights::default();
        let weak = BinderInputs {
            dg_bind_kcal: Some(-1.0),
            ..Default::default()
        };
        let strong = BinderInputs {
            dg_bind_kcal: Some(-15.0),
            ..Default::default()
        };
        assert!(score(&strong, &w).unwrap().value > score(&weak, &w).unwrap().value);

        let safe = BinderInputs {
            safety_severity: Some(0),
            ..Default::default()
        };
        let risky = BinderInputs {
            safety_severity: Some(4),
            ..Default::default()
        };
        assert!(score(&safe, &w).unwrap().value > score(&risky, &w).unwrap().value);
    }

    #[test]
    fn critical_severity_zeroes_that_channel() {
        let s = score(
            &BinderInputs {
                safety_severity: Some(4),
                ..Default::default()
            },
            &BinderWeights::default(),
        )
        .unwrap();
        assert!(s.value.abs() < 1e-12); // 1 - 4/4 = 0
    }

    #[test]
    fn absent_and_zero_weight_channels_excluded() {
        let only_dev = BinderInputs {
            developability: Some(0.7),
            ..Default::default()
        };
        let s = score(&only_dev, &BinderWeights::default()).unwrap();
        assert_eq!(s.components.len(), 1);
        assert!((s.value - 0.7).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(
            score(&BinderInputs::default(), &BinderWeights::default())
                .unwrap_err()
                .code(),
            "no_components"
        );
        let bad_dev = BinderInputs {
            developability: Some(1.5),
            ..Default::default()
        };
        assert_eq!(
            score(&bad_dev, &BinderWeights::default())
                .unwrap_err()
                .code(),
            "developability_out_of_range"
        );
        let bad_sev = BinderInputs {
            safety_severity: Some(9),
            ..Default::default()
        };
        assert_eq!(
            score(&bad_sev, &BinderWeights::default())
                .unwrap_err()
                .code(),
            "severity_out_of_range"
        );
    }

    #[test]
    fn serde_round_trips() {
        let s = score(
            &BinderInputs {
                dg_bind_kcal: Some(-8.0),
                developability: Some(0.6),
                safety_severity: Some(1),
            },
            &BinderWeights::default(),
        )
        .unwrap();
        let j = serde_json::to_string(&s).unwrap();
        let back: BinderScore = serde_json::from_str(&j).unwrap();
        assert_eq!(s, back);
    }
}
