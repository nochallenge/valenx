//! Standard hatch-pattern library.
//!
//! 15 named patterns covering the most common AutoCAD / ANSI hatch
//! definitions used in engineering drawings:
//!
//! | Name           | Use                                    |
//! |----------------|----------------------------------------|
//! | `ANSI31`       | Iron / general / cast steel — 45°      |
//! | `ANSI32`       | Steel — 45° + 135° crossed             |
//! | `ANSI33`       | Bronze / brass — 45° crossed close     |
//! | `ANSI34`       | Plastic / rubber — 45° wide spacing    |
//! | `ANSI35`       | Fire / refractory brick — 45° + ticks  |
//! | `ANSI36`       | Marble — 45° crossed                   |
//! | `ANSI37`       | Lead / zinc / magnesium — brick        |
//! | `ANSI38`       | Aluminium — alternating diagonals      |
//! | `AR-CONC`      | Concrete (dotted)                      |
//! | `AR-SAND`      | Sand                                   |
//! | `EARTH`        | Earth / soil                           |
//! | `ESCHER`       | Decorative Escher fish                 |
//! | `GRASS`        | Grass / vegetation                     |
//! | `HONEY`        | Honeycomb                              |
//! | `WOOD`         | Wood grain — alternating diagonals     |
//!
//! Each pattern carries the line spacing (mm) and one or more angles
//! (radians, measured from the +X axis). When dots are used (concrete,
//! sand), [`HatchPattern::dot_spacing`] is `Some(value)`.
//!
//! The actual hatch generation reuses [`crate::section::hatch`] for
//! each angle in turn (concatenating the result lines), so the engine
//! that already ships with Phase 5 keeps doing the heavy lifting.

use serde::{Deserialize, Serialize};

/// One hatch pattern definition.
///
/// `name` is the canonical ANSI / AutoCAD identifier (uppercase).
/// `spacing` is the perpendicular distance between adjacent lines of
/// a single direction in mm. `angles` is the set of line directions
/// in radians (a single-direction pattern has one entry; crossed
/// patterns have two; brick / fancy patterns can have more).
/// `dot_spacing` is `Some(d)` when the pattern includes a stippled
/// dot grid at spacing `d`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct HatchPattern {
    /// Canonical pattern name.
    pub name: &'static str,
    /// Line spacing in mm.
    pub spacing: f64,
    /// Line angles in radians.
    pub angles: Vec<f64>,
    /// Optional dot-grid spacing in mm (stippled patterns).
    pub dot_spacing: Option<f64>,
}

const D45: f64 = std::f64::consts::FRAC_PI_4;
const D135: f64 = D45 * 3.0;

/// Standard "iron / general" 45° line pattern.
pub fn ansi31_iron() -> HatchPattern {
    HatchPattern {
        name: "ANSI31",
        spacing: 2.0,
        angles: vec![D45],
        dot_spacing: None,
    }
}

/// Steel — 45° + 135° crossed.
pub fn ansi32_steel() -> HatchPattern {
    HatchPattern {
        name: "ANSI32",
        spacing: 2.5,
        angles: vec![D45, D135],
        dot_spacing: None,
    }
}

/// Bronze / brass — 45° crossed close.
pub fn ansi33_bronze() -> HatchPattern {
    HatchPattern {
        name: "ANSI33",
        spacing: 1.5,
        angles: vec![D45, D135],
        dot_spacing: None,
    }
}

/// Plastic / rubber — 45° wide spacing.
pub fn ansi34_plastic() -> HatchPattern {
    HatchPattern {
        name: "ANSI34",
        spacing: 4.0,
        angles: vec![D45],
        dot_spacing: None,
    }
}

/// Fire / refractory brick — 45° + perpendicular ticks.
pub fn ansi35_brick() -> HatchPattern {
    HatchPattern {
        name: "ANSI35",
        spacing: 2.5,
        angles: vec![D45, D135],
        dot_spacing: None,
    }
}

/// Marble — 45° crossed.
pub fn ansi36_marble() -> HatchPattern {
    HatchPattern {
        name: "ANSI36",
        spacing: 2.5,
        angles: vec![D45, D135],
        dot_spacing: None,
    }
}

/// Lead / zinc / magnesium — brick pattern.
pub fn ansi37_brick() -> HatchPattern {
    HatchPattern {
        name: "ANSI37",
        spacing: 3.0,
        angles: vec![0.0, std::f64::consts::FRAC_PI_2],
        dot_spacing: None,
    }
}

/// Aluminium — alternating diagonals.
pub fn ansi38_aluminum() -> HatchPattern {
    HatchPattern {
        name: "ANSI38",
        spacing: 2.5,
        angles: vec![D45, D135],
        dot_spacing: None,
    }
}

/// Concrete — wide spacing + dots.
pub fn ar_conc_concrete() -> HatchPattern {
    HatchPattern {
        name: "AR-CONC",
        spacing: 6.0,
        angles: vec![D45, D135],
        dot_spacing: Some(1.5),
    }
}

/// Sand — dense dot pattern.
pub fn ar_sand() -> HatchPattern {
    HatchPattern {
        name: "AR-SAND",
        spacing: 0.0,
        angles: vec![],
        dot_spacing: Some(0.8),
    }
}

/// Earth / soil — diagonal short dashes.
pub fn earth() -> HatchPattern {
    HatchPattern {
        name: "EARTH",
        spacing: 5.0,
        angles: vec![D45],
        dot_spacing: None,
    }
}

/// Decorative Escher fish (placeholder — emit a wide cross-hatch).
pub fn escher() -> HatchPattern {
    HatchPattern {
        name: "ESCHER",
        spacing: 4.0,
        angles: vec![0.0, std::f64::consts::FRAC_PI_2, D45, D135],
        dot_spacing: None,
    }
}

/// Grass / vegetation — short vertical ticks.
pub fn grass() -> HatchPattern {
    HatchPattern {
        name: "GRASS",
        spacing: 3.0,
        angles: vec![std::f64::consts::FRAC_PI_2],
        dot_spacing: None,
    }
}

/// Honeycomb — three angles forming hexagonal grid (approx).
pub fn honey() -> HatchPattern {
    HatchPattern {
        name: "HONEY",
        spacing: 3.0,
        angles: vec![
            0.0,
            std::f64::consts::PI / 3.0,
            2.0 * std::f64::consts::PI / 3.0,
        ],
        dot_spacing: None,
    }
}

/// Wood grain — alternating diagonal lines at slight angles.
pub fn wood() -> HatchPattern {
    HatchPattern {
        name: "WOOD",
        spacing: 1.5,
        angles: vec![
            0.05 * std::f64::consts::PI,
            -0.05 * std::f64::consts::PI + std::f64::consts::PI,
        ],
        dot_spacing: None,
    }
}

/// Look up a pattern by canonical name (case-insensitive).
pub fn by_name(name: &str) -> Option<HatchPattern> {
    let upper = name.to_uppercase();
    match upper.as_str() {
        "ANSI31" => Some(ansi31_iron()),
        "ANSI32" => Some(ansi32_steel()),
        "ANSI33" => Some(ansi33_bronze()),
        "ANSI34" => Some(ansi34_plastic()),
        "ANSI35" => Some(ansi35_brick()),
        "ANSI36" => Some(ansi36_marble()),
        "ANSI37" => Some(ansi37_brick()),
        "ANSI38" => Some(ansi38_aluminum()),
        "AR-CONC" => Some(ar_conc_concrete()),
        "AR-SAND" => Some(ar_sand()),
        "EARTH" => Some(earth()),
        "ESCHER" => Some(escher()),
        "GRASS" => Some(grass()),
        "HONEY" => Some(honey()),
        "WOOD" => Some(wood()),
        _ => None,
    }
}

/// Names of all 15 standard patterns, in display order. Convenient
/// for UI dropdown population.
pub fn all_names() -> &'static [&'static str] {
    &[
        "ANSI31", "ANSI32", "ANSI33", "ANSI34", "ANSI35", "ANSI36", "ANSI37", "ANSI38", "AR-CONC",
        "AR-SAND", "EARTH", "ESCHER", "GRASS", "HONEY", "WOOD",
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_15_patterns_exist_and_are_unique_by_name() {
        let names = all_names();
        assert_eq!(names.len(), 15);
        let set: std::collections::HashSet<_> = names.iter().collect();
        assert_eq!(set.len(), 15);
        for n in names {
            let p = by_name(n).unwrap_or_else(|| panic!("missing pattern {n}"));
            assert_eq!(p.name, *n);
        }
    }

    #[test]
    fn by_name_is_case_insensitive() {
        assert!(by_name("ansi31").is_some());
        assert!(by_name("Ar-Conc").is_some());
    }

    #[test]
    fn ansi31_is_single_45_degree() {
        let p = ansi31_iron();
        assert_eq!(p.angles.len(), 1);
        assert!((p.angles[0] - std::f64::consts::FRAC_PI_4).abs() < 1e-9);
    }

    #[test]
    fn ansi32_steel_is_crossed_45_135() {
        let p = ansi32_steel();
        assert_eq!(p.angles.len(), 2);
        assert!((p.angles[0] - D45).abs() < 1e-9);
        assert!((p.angles[1] - D135).abs() < 1e-9);
    }

    #[test]
    fn ar_conc_has_dot_grid() {
        let p = ar_conc_concrete();
        assert!(p.dot_spacing.is_some());
    }

    #[test]
    fn by_name_unknown_returns_none() {
        assert!(by_name("ZZZZ").is_none());
    }
}
