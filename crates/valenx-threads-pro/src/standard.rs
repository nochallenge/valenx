//! Extended thread-standard family — adds Acme / Trapezoidal /
//! Whitworth on top of the Phase 13 [`valenx_feature_tree::threads::ThreadStandard`].
//!
//! The legacy four (`IsoMetric`, `UnifiedNational`, `BSPP`, `NPT`)
//! remain in valenx-feature-tree because the bolt / nut / hole
//! features take a [`valenx_feature_tree::threads::ThreadSpec`]. The
//! [`ThreadStandardPro`] enum here is a *superset* tag that the new
//! tables (full metric / UN / BSP / NPT / Acme / ...) return so
//! callers can reason about thread shape (`V`, `Acme`, `Trapezoidal`,
//! `Buttress`).

use serde::{Deserialize, Serialize};

/// Extended thread-standard family.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum ThreadStandardPro {
    /// ISO 261 metric coarse-pitch series (M0.5 … M300).
    IsoMetric,
    /// ISO metric fine (M0.5x0.3 … M300x6).
    MetricFine,
    /// ISO metric extra-fine.
    MetricExtraFine,
    /// ANSI/ASME B1.1 Unified National Coarse (UNC).
    UnifiedNationalCoarse,
    /// Unified National Fine (UNF).
    UnifiedNationalFine,
    /// Unified National Extra-Fine (UNEF).
    UnifiedNationalExtraFine,
    /// BS EN ISO 228 British Standard Pipe Parallel (BSPP, G1/16 … G6).
    BSPP,
    /// British Standard Pipe Taper (BSPT, R1/16 … R6).
    BSPT,
    /// ANSI/ASME B1.20.1 American National Standard Pipe Taper (NPT).
    NPT,
    /// Dryseal Pipe Taper (NPTF).
    NPTF,
    /// American National Standard Straight Pipe Thread (NPS).
    NPS,
    /// ASME B1.5 Acme thread (29° trapezoidal, general-purpose).
    Acme,
    /// ISO 2904 metric trapezoidal (30°, ISO 2904 / DIN 103).
    TrapezoidalIso,
    /// BS 84 / DIN 11 Whitworth British Standard Whitworth (BSW).
    WhitworthBSW,
    /// ASME B1.9 Buttress thread (push-direction-asymmetric).
    Buttress,
}

impl ThreadStandardPro {
    /// Short label for UI display.
    pub fn short_label(self) -> &'static str {
        match self {
            ThreadStandardPro::IsoMetric => "ISO Metric",
            ThreadStandardPro::MetricFine => "ISO Metric Fine",
            ThreadStandardPro::MetricExtraFine => "ISO Metric XFine",
            ThreadStandardPro::UnifiedNationalCoarse => "UNC",
            ThreadStandardPro::UnifiedNationalFine => "UNF",
            ThreadStandardPro::UnifiedNationalExtraFine => "UNEF",
            ThreadStandardPro::BSPP => "BSPP (G)",
            ThreadStandardPro::BSPT => "BSPT (R)",
            ThreadStandardPro::NPT => "NPT",
            ThreadStandardPro::NPTF => "NPTF",
            ThreadStandardPro::NPS => "NPS",
            ThreadStandardPro::Acme => "Acme",
            ThreadStandardPro::TrapezoidalIso => "Trapezoidal (ISO 2904)",
            ThreadStandardPro::WhitworthBSW => "Whitworth BSW",
            ThreadStandardPro::Buttress => "Buttress",
        }
    }

    /// True for tapered-pipe families.
    pub fn is_tapered(self) -> bool {
        matches!(
            self,
            ThreadStandardPro::NPT | ThreadStandardPro::NPTF | ThreadStandardPro::BSPT
        )
    }
}

/// The geometric cross-section of one thread tooth.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum ProfileShape {
    /// 60° symmetric V (ISO metric / UN / BSPP / NPT / Whitworth-55°
    /// variant rounding is treated as V here).
    V,
    /// 29° flat-flank Acme.
    Acme,
    /// 30° flat-flank ISO trapezoidal.
    Trapezoidal,
    /// 7°/45° asymmetric Buttress.
    Buttress,
}

impl ThreadStandardPro {
    /// Map family → profile shape.
    pub fn profile_shape(self) -> ProfileShape {
        match self {
            ThreadStandardPro::Acme => ProfileShape::Acme,
            ThreadStandardPro::TrapezoidalIso => ProfileShape::Trapezoidal,
            ThreadStandardPro::Buttress => ProfileShape::Buttress,
            _ => ProfileShape::V,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn short_label_matches_family() {
        assert_eq!(ThreadStandardPro::Acme.short_label(), "Acme");
        assert_eq!(ThreadStandardPro::BSPT.short_label(), "BSPT (R)");
    }

    #[test]
    fn is_tapered_recognises_pipe_taper_families() {
        assert!(ThreadStandardPro::NPT.is_tapered());
        assert!(ThreadStandardPro::BSPT.is_tapered());
        assert!(!ThreadStandardPro::BSPP.is_tapered());
        assert!(!ThreadStandardPro::Acme.is_tapered());
    }

    #[test]
    fn profile_shape_matches_family() {
        assert_eq!(ThreadStandardPro::Acme.profile_shape(), ProfileShape::Acme);
        assert_eq!(ThreadStandardPro::Buttress.profile_shape(), ProfileShape::Buttress);
        assert_eq!(ThreadStandardPro::IsoMetric.profile_shape(), ProfileShape::V);
    }
}
