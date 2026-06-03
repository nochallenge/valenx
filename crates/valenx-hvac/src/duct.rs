//! Ducting types + CAD emit.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use valenx_cad::primitives::{box_solid, cylinder};
use valenx_cad::Solid;

use crate::error::HvacError;

/// Cross-section of a duct.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum CrossSection {
    /// Rectangular duct (width × height, mm).
    Rect {
        /// Width in mm.
        w: f64,
        /// Height in mm.
        h: f64,
    },
    /// Round duct (diameter in mm).
    Round {
        /// Diameter in mm.
        d: f64,
    },
}

impl CrossSection {
    /// Hydraulic diameter `D_h = 4A / P`.
    pub fn hydraulic_diameter_mm(self) -> f64 {
        match self {
            CrossSection::Round { d } => d,
            CrossSection::Rect { w, h } => 2.0 * w * h / (w + h),
        }
    }

    /// Flow area in mm².
    pub fn area_mm2(self) -> f64 {
        match self {
            CrossSection::Round { d } => std::f64::consts::PI * (d * 0.5).powi(2),
            CrossSection::Rect { w, h } => w * h,
        }
    }
}

/// A length of ducting following a polyline path with optional
/// insulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Duct {
    /// Cross-section.
    pub cross_section: CrossSection,
    /// Polyline path in world coordinates (mm).
    pub path: Vec<Vector3<f64>>,
    /// Insulation thickness in mm (0.0 = none).
    pub insulation_thickness: f64,
}

impl Duct {
    /// Construct directly.
    pub fn new(cross_section: CrossSection, path: Vec<Vector3<f64>>, insulation_thickness: f64) -> Self {
        Self {
            cross_section,
            path,
            insulation_thickness,
        }
    }

    /// Total polyline length in mm.
    pub fn path_length_mm(&self) -> f64 {
        self.path
            .windows(2)
            .map(|s| (s[1] - s[0]).norm())
            .sum::<f64>()
    }
}

/// Emit a CAD [`Solid`] for the duct. v1: builds a single straight
/// extrusion of the cross-section along +Z with length =
/// `path_length_mm`. v2 will sweep the cross-section along the actual
/// path polyline.
pub fn to_solid(d: &Duct) -> Result<Solid, HvacError> {
    let len = d.path_length_mm();
    if len <= 0.0 {
        return Err(HvacError::BadParameter {
            name: "duct.path",
            reason: "path length must be > 0".into(),
        });
    }
    match d.cross_section {
        CrossSection::Round { d: dia } => {
            cylinder(dia / 2.0, len).map_err(|e| HvacError::Cad(e.to_string()))
        }
        CrossSection::Rect { w, h } => {
            box_solid(w, h, len).map_err(|e| HvacError::Cad(e.to_string()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hydraulic_diameter_round_equals_diameter() {
        let s = CrossSection::Round { d: 200.0 };
        assert!((s.hydraulic_diameter_mm() - 200.0).abs() < 1e-9);
    }

    #[test]
    fn hydraulic_diameter_rect_matches_formula() {
        let s = CrossSection::Rect { w: 200.0, h: 100.0 };
        // 4·A/P = 4·200·100 / (2·(200+100)) = 80000 / 600 = 133.333
        assert!((s.hydraulic_diameter_mm() - 133.3333).abs() < 1e-3);
    }

    #[test]
    fn area_round_matches_pi_r_squared() {
        let s = CrossSection::Round { d: 100.0 };
        // π·(50)^2 = 7853.98
        assert!((s.area_mm2() - 7853.98).abs() < 1e-1);
    }

    #[test]
    fn path_length_sums_segments() {
        let d = Duct::new(
            CrossSection::Round { d: 100.0 },
            vec![
                Vector3::zeros(),
                Vector3::new(100.0, 0.0, 0.0),
                Vector3::new(100.0, 50.0, 0.0),
            ],
            0.0,
        );
        assert!((d.path_length_mm() - 150.0).abs() < 1e-9);
    }

    #[test]
    fn to_solid_rejects_empty_path() {
        let d = Duct::new(
            CrossSection::Round { d: 100.0 },
            vec![Vector3::zeros()],
            0.0,
        );
        let err = to_solid(&d).unwrap_err();
        assert!(matches!(err, HvacError::BadParameter { .. }));
    }
}
