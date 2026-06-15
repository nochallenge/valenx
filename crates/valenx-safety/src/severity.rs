//! The ordered severity scale for risk flags.

use serde::{Deserialize, Serialize};

/// Risk severity, ordered `Info < Low < Moderate < High < Critical`.
///
/// The aggregate severity of a [`crate::RiskReport`] is the worst flag's
/// severity — never a "safe" verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Severity {
    /// Informational; no concern, but recorded.
    Info,
    /// Low concern.
    Low,
    /// Moderate concern; worth a look.
    Moderate,
    /// High concern; a reviewer should scrutinise this.
    High,
    /// Critical concern; a likely blocker pending review.
    Critical,
}

impl Severity {
    /// A lower-case label for the severity.
    pub fn as_str(&self) -> &'static str {
        match self {
            Severity::Info => "info",
            Severity::Low => "low",
            Severity::Moderate => "moderate",
            Severity::High => "high",
            Severity::Critical => "critical",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ordering_is_increasing() {
        assert!(Severity::Info < Severity::Low);
        assert!(Severity::Low < Severity::Moderate);
        assert!(Severity::Moderate < Severity::High);
        assert!(Severity::High < Severity::Critical);
        // max picks the worst
        assert_eq!(
            [Severity::Low, Severity::Critical, Severity::Moderate]
                .into_iter()
                .max()
                .unwrap(),
            Severity::Critical
        );
    }
}
