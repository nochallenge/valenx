//! Risk flags and the per-screen constructors that produce them.

use serde::{Deserialize, Serialize};

use crate::error::{check_finite, SafetyError};
use crate::severity::Severity;

/// One risk flag raised by a screen: its `source`, `severity`, a human `detail`,
/// and optional numeric `evidence` (e.g. the identity or density that triggered
/// it).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RiskFlag {
    /// The screen that raised the flag (e.g. `"off-target"`).
    pub source: String,
    /// How serious the flag is.
    pub severity: Severity,
    /// Human-readable description.
    pub detail: String,
    /// The numeric evidence that triggered the flag, if any.
    pub evidence: Option<f64>,
}

impl RiskFlag {
    /// A new flag.
    pub fn new(source: impl Into<String>, severity: Severity, detail: impl Into<String>) -> Self {
        Self {
            source: source.into(),
            severity,
            detail: detail.into(),
            evidence: None,
        }
    }

    /// Attach numeric evidence.
    pub fn with_evidence(mut self, value: f64) -> Self {
        self.evidence = Some(value);
        self
    }
}

/// Build an off-target / cross-reactivity flag from a sequence-identity hit.
/// Returns `None` when `identity` is below `threshold`. Severity bands
/// (illustrative): `>= 0.95` critical, `>= 0.85` high, `>= 0.70` moderate, else
/// low.
pub fn offtarget_flag(
    reference: &str,
    identity: f64,
    threshold: f64,
) -> Result<Option<RiskFlag>, SafetyError> {
    check_finite("identity", identity)?;
    check_finite("threshold", threshold)?;
    if identity < threshold {
        return Ok(None);
    }
    let severity = if identity >= 0.95 {
        Severity::Critical
    } else if identity >= 0.85 {
        Severity::High
    } else if identity >= 0.70 {
        Severity::Moderate
    } else {
        Severity::Low
    };
    Ok(Some(
        RiskFlag::new(
            "off-target",
            severity,
            format!("{:.1}% sequence identity to {reference}", identity * 100.0),
        )
        .with_evidence(identity),
    ))
}

/// Build an immunogenicity flag from a predicted T-cell epitope density. Returns
/// `None` below `threshold`. Severity bands (illustrative): `>= 0.50` high,
/// `>= 0.25` moderate, else low.
pub fn immunogenicity_flag(
    epitope_density: f64,
    threshold: f64,
) -> Result<Option<RiskFlag>, SafetyError> {
    check_finite("epitope_density", epitope_density)?;
    check_finite("threshold", threshold)?;
    if epitope_density < threshold {
        return Ok(None);
    }
    let severity = if epitope_density >= 0.50 {
        Severity::High
    } else if epitope_density >= 0.25 {
        Severity::Moderate
    } else {
        Severity::Low
    };
    Ok(Some(
        RiskFlag::new(
            "immunogenicity",
            severity,
            format!("predicted T-cell epitope density {epitope_density:.2}"),
        )
        .with_evidence(epitope_density),
    ))
}

/// Build a CRISPR off-target flag from a predicted off-target edit-site count.
/// Returns `None` below `threshold`. Severity bands (illustrative): `>= 10`
/// high, `>= 3` moderate, else low.
pub fn crispr_offtarget_flag(
    predicted_sites: u32,
    threshold: u32,
) -> Result<Option<RiskFlag>, SafetyError> {
    if predicted_sites < threshold {
        return Ok(None);
    }
    let severity = if predicted_sites >= 10 {
        Severity::High
    } else if predicted_sites >= 3 {
        Severity::Moderate
    } else {
        Severity::Low
    };
    Ok(Some(
        RiskFlag::new(
            "crispr-off-target",
            severity,
            format!("{predicted_sites} predicted off-target edit site(s)"),
        )
        .with_evidence(f64::from(predicted_sites)),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn offtarget_gdf11_is_high() {
        // GDF-11 vs an anti-myostatin design: ~90% identity -> High.
        let f = offtarget_flag("GDF-11", 0.899, 0.80).unwrap().unwrap();
        assert_eq!(f.severity, Severity::High);
        assert_eq!(f.source, "off-target");
        assert_eq!(f.evidence, Some(0.899));
    }

    #[test]
    fn offtarget_below_threshold_is_none() {
        assert!(offtarget_flag("BMP-7", 0.23, 0.80).unwrap().is_none());
    }

    #[test]
    fn immunogenicity_bands() {
        assert_eq!(
            immunogenicity_flag(0.6, 0.1).unwrap().unwrap().severity,
            Severity::High
        );
        assert_eq!(
            immunogenicity_flag(0.3, 0.1).unwrap().unwrap().severity,
            Severity::Moderate
        );
        assert!(immunogenicity_flag(0.05, 0.1).unwrap().is_none());
    }

    #[test]
    fn crispr_bands() {
        assert_eq!(
            crispr_offtarget_flag(12, 1).unwrap().unwrap().severity,
            Severity::High
        );
        assert_eq!(
            crispr_offtarget_flag(4, 1).unwrap().unwrap().severity,
            Severity::Moderate
        );
        assert!(crispr_offtarget_flag(0, 1).unwrap().is_none());
    }

    #[test]
    fn rejects_non_finite() {
        assert_eq!(
            offtarget_flag("x", f64::NAN, 0.8).unwrap_err().code(),
            "non_finite"
        );
    }
}
