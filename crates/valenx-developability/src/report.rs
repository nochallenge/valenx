//! The consolidated developability report.

use serde::{Deserialize, Serialize};

use crate::charge::{isoelectric_point, net_charge_at_ph};
use crate::error::DevelopabilityError;
use crate::hydrophobicity::{aggregation_prone_regions, gravy};

/// Default sliding window for aggregation-prone-region detection.
pub const DEFAULT_APR_WINDOW: usize = 7;
/// Default mean-hydropathy threshold for an aggregation-prone window.
pub const DEFAULT_APR_THRESHOLD: f64 = 1.5;
/// Physiological pH used for the reported net charge.
pub const PHYSIOLOGICAL_PH: f64 = 7.4;

/// A sequence's developability summary plus heuristic flags.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DevelopabilityReport {
    /// Grand average of hydropathy (mean Kyte-Doolittle value).
    pub gravy: f64,
    /// Isoelectric point (pH of zero net charge).
    pub isoelectric_point: f64,
    /// Net charge at physiological pH (7.4).
    pub net_charge_physiological: f64,
    /// Aggregation-prone-region spans `(start, end)` (end exclusive).
    pub aggregation_prone_regions: Vec<(usize, usize)>,
    /// Heuristic developability flags for human follow-up.
    pub flags: Vec<String>,
}

/// Assess a sequence with default settings (`DEFAULT_APR_WINDOW`,
/// `DEFAULT_APR_THRESHOLD`). The window is clamped to the sequence length so
/// short sequences are still assessed.
pub fn assess(seq: &str) -> Result<DevelopabilityReport, DevelopabilityError> {
    let g = gravy(seq)?;
    let pi = isoelectric_point(seq)?;
    let net = net_charge_at_ph(seq, PHYSIOLOGICAL_PH)?;
    let window = DEFAULT_APR_WINDOW.min(seq.len());
    let aprs = aggregation_prone_regions(seq, window, DEFAULT_APR_THRESHOLD)?;

    let mut flags = Vec::new();
    if g > 0.0 {
        flags.push(format!(
            "net hydrophobic (GRAVY {g:.2}); review solubility / aggregation"
        ));
    }
    if !aprs.is_empty() {
        flags.push(format!("{} aggregation-prone region(s)", aprs.len()));
    }
    if !(5.0..=9.0).contains(&pi) {
        flags.push(format!(
            "extreme isoelectric point ({pi:.1}); low solubility near physiological pH is possible"
        ));
    }

    Ok(DevelopabilityReport {
        gravy: g,
        isoelectric_point: pi,
        net_charge_physiological: net,
        aggregation_prone_regions: aprs,
        flags,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hydrophobic_sequence_is_flagged() {
        let r = assess("IIIIVVVVLLLL").unwrap();
        assert!(r.gravy > 0.0);
        assert!(!r.aggregation_prone_regions.is_empty());
        assert!(r.flags.iter().any(|f| f.contains("hydrophobic")));
        assert!(r.flags.iter().any(|f| f.contains("aggregation-prone")));
    }

    #[test]
    fn hydrophilic_balanced_sequence_is_calmer() {
        // a mixed, hydrophilic-leaning sequence: negative GRAVY, no APRs
        let r = assess("STNQGDEKRH").unwrap();
        assert!(r.gravy < 0.0);
        assert!(r.aggregation_prone_regions.is_empty());
        assert!(!r.flags.iter().any(|f| f.contains("aggregation-prone")));
    }

    #[test]
    fn report_fields_are_populated() {
        let r = assess("ACDEFGHIKLMNPQRSTVWY").unwrap();
        assert!(r.isoelectric_point > 0.0 && r.isoelectric_point < 14.0);
        assert!(r.gravy.is_finite());
        assert!(r.net_charge_physiological.is_finite());
    }

    #[test]
    fn serde_round_trips() {
        let r = assess("IIIIVVVVLLLL").unwrap();
        let j = serde_json::to_string(&r).unwrap();
        let back: DevelopabilityReport = serde_json::from_str(&j).unwrap();
        assert_eq!(r, back);
    }

    #[test]
    fn rejects_empty() {
        assert_eq!(assess("").unwrap_err().code(), "empty_sequence");
    }
}
