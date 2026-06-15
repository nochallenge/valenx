//! A scored candidate row in a dossier.

use serde::{Deserialize, Serialize};

use crate::error::DossierError;

/// One ranked candidate: its unified comparable score, the broken-out score
/// components, an optional calibrated confidence, and any safety flags.
///
/// An **empty `safety_flags` list is not a safety assertion** — it means no
/// screen raised a flag, which is different from "safe".
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScoredCandidate {
    /// Candidate identifier.
    pub id: String,
    /// Unified comparable score (higher is better).
    pub comparable_score: f64,
    /// Broken-out score components `(name, value)` (e.g. dock, MM-GBSA, ipTM).
    pub components: Vec<(String, f64)>,
    /// Calibrated confidence in `[0, 1]`, if calibration was available.
    pub calibrated_confidence: Option<f64>,
    /// Safety flags raised by upstream screens (off-target, immunogenicity …).
    pub safety_flags: Vec<String>,
}

impl ScoredCandidate {
    /// A new candidate with a finite comparable score.
    pub fn new(id: impl Into<String>, comparable_score: f64) -> Result<Self, DossierError> {
        let id = id.into();
        if id.trim().is_empty() {
            return Err(DossierError::Empty {
                what: "candidate id",
            });
        }
        if !comparable_score.is_finite() {
            return Err(DossierError::NonFinite {
                what: "comparable_score",
            });
        }
        Ok(Self {
            id,
            comparable_score,
            components: Vec::new(),
            calibrated_confidence: None,
            safety_flags: Vec::new(),
        })
    }

    /// Add a broken-out score component (value must be finite).
    pub fn with_component(
        mut self,
        name: impl Into<String>,
        value: f64,
    ) -> Result<Self, DossierError> {
        if !value.is_finite() {
            return Err(DossierError::NonFinite { what: "component" });
        }
        self.components.push((name.into(), value));
        Ok(self)
    }

    /// Attach a calibrated confidence (must be in `[0, 1]`).
    pub fn with_confidence(mut self, confidence: f64) -> Result<Self, DossierError> {
        if !confidence.is_finite() || !(0.0..=1.0).contains(&confidence) {
            return Err(DossierError::ConfidenceOutOfRange { value: confidence });
        }
        self.calibrated_confidence = Some(confidence);
        Ok(self)
    }

    /// Add a safety flag.
    pub fn with_flag(mut self, flag: impl Into<String>) -> Self {
        self.safety_flags.push(flag.into());
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_and_validates() {
        let c = ScoredCandidate::new("design_A", 0.8)
            .unwrap()
            .with_component("dock", -9.1)
            .unwrap()
            .with_confidence(0.7)
            .unwrap()
            .with_flag("GDF11_crossreactivity");
        assert_eq!(c.id, "design_A");
        assert_eq!(c.components.len(), 1);
        assert_eq!(c.calibrated_confidence, Some(0.7));
        assert_eq!(c.safety_flags, vec!["GDF11_crossreactivity"]);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(ScoredCandidate::new("", 1.0).unwrap_err().code(), "empty");
        assert_eq!(
            ScoredCandidate::new("x", f64::NAN).unwrap_err().code(),
            "non_finite"
        );
        assert_eq!(
            ScoredCandidate::new("x", 1.0)
                .unwrap()
                .with_confidence(1.5)
                .unwrap_err()
                .code(),
            "confidence_out_of_range"
        );
    }
}
