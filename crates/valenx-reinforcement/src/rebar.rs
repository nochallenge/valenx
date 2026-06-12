//! Rebar — one reinforcing bar.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Steel grade (yield-stress class, US standard).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum RebarGrade {
    /// Grade 40 (40 ksi yield).
    G40,
    /// Grade 60 (60 ksi yield).
    G60,
    /// Grade 75 (75 ksi yield).
    G75,
}

impl RebarGrade {
    /// Yield stress in ksi.
    pub fn yield_ksi(self) -> u32 {
        match self {
            Self::G40 => 40,
            Self::G60 => 60,
            Self::G75 => 75,
        }
    }

    /// UI dropdown label.
    pub fn label(self) -> &'static str {
        match self {
            Self::G40 => "Grade 40",
            Self::G60 => "Grade 60",
            Self::G75 => "Grade 75",
        }
    }
}

impl Default for RebarGrade {
    fn default() -> Self {
        Self::G60
    }
}

/// Rebar shape recipe — the centreline geometry of the bar.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum RebarShape {
    /// Straight bar of `length` along +X starting at origin.
    Straight {
        /// Bar length (m).
        length: f64,
    },
    /// L-shape: a `leg_a`-long horizontal leg along +X then bends
    /// 90° downward into a `leg_b`-long leg along -Z.
    L {
        /// First leg length (m).
        leg_a: f64,
        /// Second leg length (m).
        leg_b: f64,
    },
    /// U-shape (stirrup): rectangle with `width` × `height`, open at
    /// the top. The base centreline runs along +X at z = -height.
    U {
        /// Base width (m).
        width: f64,
        /// Side height (m).
        height: f64,
    },
    /// 180° hook at the end of a straight bar — a small semicircle of
    /// radius `bend_radius` after a `straight` length along +X.
    Hook {
        /// Straight tail length (m).
        straight: f64,
        /// Hook bend radius (m).
        bend_radius: f64,
    },
    /// Spiral hoop — helix of `pitch`-per-turn, `radius`, `turns`
    /// turns, centred on the Z axis.
    Spiral {
        /// Helix radius (m).
        radius: f64,
        /// Vertical rise per full turn (m).
        pitch: f64,
        /// Number of full turns.
        turns: f64,
    },
}

impl Default for RebarShape {
    fn default() -> Self {
        Self::Straight { length: 1.0 }
    }
}

impl RebarShape {
    /// Sample the bar's centreline into a 3D polyline. Vertex count
    /// scales with the visual complexity of each shape variant.
    pub fn to_polyline(&self) -> Vec<Vector3<f64>> {
        match self {
            Self::Straight { length } => vec![Vector3::zeros(), Vector3::new(*length, 0.0, 0.0)],
            Self::L { leg_a, leg_b } => vec![
                Vector3::zeros(),
                Vector3::new(*leg_a, 0.0, 0.0),
                Vector3::new(*leg_a, 0.0, -*leg_b),
            ],
            Self::U { width, height } => vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(0.0, 0.0, -*height),
                Vector3::new(*width, 0.0, -*height),
                Vector3::new(*width, 0.0, 0.0),
            ],
            Self::Hook {
                straight,
                bend_radius,
            } => {
                // Straight tail then a 180° arc.
                let mut pts = vec![Vector3::zeros(), Vector3::new(*straight, 0.0, 0.0)];
                let r = *bend_radius;
                let cx = *straight;
                let cz = -r;
                let segs = 12;
                for i in 1..=segs {
                    let theta = std::f64::consts::FRAC_PI_2
                        - (i as f64 / segs as f64) * std::f64::consts::PI;
                    pts.push(Vector3::new(
                        cx + r * theta.sin(),
                        0.0,
                        cz + r * theta.cos(),
                    ));
                }
                pts
            }
            Self::Spiral {
                radius,
                pitch,
                turns,
            } => {
                let total = (*turns) * 16.0;
                // A spiral with no (or negative) turns has no geometry; bail
                // before the loop so `t = i / total` can't divide by zero and
                // emit NaN vertices.
                if total <= 0.0 {
                    return Vec::new();
                }
                let n = total.ceil() as usize + 1;
                let mut pts = Vec::with_capacity(n);
                for i in 0..n {
                    let t = (i as f64) / total;
                    let theta = t * turns * std::f64::consts::TAU;
                    pts.push(Vector3::new(
                        radius * theta.cos(),
                        radius * theta.sin(),
                        t * turns * pitch,
                    ));
                }
                pts
            }
        }
    }

    /// UI label.
    pub fn label(&self) -> &'static str {
        match self {
            Self::Straight { .. } => "Straight",
            Self::L { .. } => "L",
            Self::U { .. } => "U",
            Self::Hook { .. } => "Hook",
            Self::Spiral { .. } => "Spiral",
        }
    }
}

/// One rebar — geometry + steel class.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Rebar {
    /// Bar nominal diameter in mm.
    pub diameter_mm: f64,
    /// Centreline length in metres (set by callers; cage generators
    /// fill this from the shape).
    pub length_m: f64,
    /// Shape recipe.
    pub shape: RebarShape,
    /// Steel grade.
    pub grade: RebarGrade,
}

impl Default for Rebar {
    fn default() -> Self {
        Self {
            diameter_mm: 16.0,
            length_m: 1.0,
            shape: RebarShape::default(),
            grade: RebarGrade::default(),
        }
    }
}

/// Cross-sectional area of a circular reinforcing bar in mm², `A = π·d²/4` for nominal diameter
/// `diameter_mm`. Feeds reinforcement ratios (ρ = As/Ag) and steel-area takedowns. Returns `0.0`
/// for a non-positive or non-finite diameter.
pub fn rebar_cross_section_area_mm2(diameter_mm: f64) -> f64 {
    if !diameter_mm.is_finite() || diameter_mm <= 0.0 {
        return 0.0;
    }
    let r = diameter_mm * 0.5;
    std::f64::consts::PI * r * r
}

/// Longitudinal reinforcement ratio `ρ = As / Ag` (ACI 318) — total steel area over the gross
/// concrete section area. Returns `0.0` if either argument is non-finite or the gross area is
/// non-positive.
pub fn reinforcement_ratio(steel_area_mm2: f64, gross_area_mm2: f64) -> f64 {
    if !steel_area_mm2.is_finite() || !gross_area_mm2.is_finite() || gross_area_mm2 <= 0.0 {
        return 0.0;
    }
    steel_area_mm2 / gross_area_mm2
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grade_yields() {
        assert_eq!(RebarGrade::G60.yield_ksi(), 60);
    }

    #[test]
    fn straight_polyline_has_2_pts() {
        let s = RebarShape::Straight { length: 3.0 };
        assert_eq!(s.to_polyline().len(), 2);
    }

    #[test]
    fn l_polyline_has_3_pts() {
        let s = RebarShape::L {
            leg_a: 1.0,
            leg_b: 0.5,
        };
        assert_eq!(s.to_polyline().len(), 3);
    }

    #[test]
    fn hook_polyline_has_2_plus_arc() {
        let s = RebarShape::Hook {
            straight: 0.5,
            bend_radius: 0.1,
        };
        let p = s.to_polyline();
        assert!(p.len() >= 14);
    }

    #[test]
    fn spiral_polyline_scales_with_turns() {
        let s = RebarShape::Spiral {
            radius: 0.3,
            pitch: 0.1,
            turns: 5.0,
        };
        assert!(s.to_polyline().len() > 40);
    }

    #[test]
    fn spiral_with_zero_turns_is_empty_not_nan() {
        // turns = 0 previously made `t = i / (turns*16)` compute 0/0 = NaN,
        // emitting (NaN, NaN, NaN) vertices. It must yield no geometry.
        let zero = RebarShape::Spiral {
            radius: 0.3,
            pitch: 0.1,
            turns: 0.0,
        };
        assert!(
            zero.to_polyline().is_empty(),
            "a 0-turn spiral has no geometry"
        );
        // A normal spiral's vertices are all finite (no NaN leaks).
        let good = RebarShape::Spiral {
            radius: 0.3,
            pitch: 0.1,
            turns: 3.0,
        };
        assert!(good
            .to_polyline()
            .iter()
            .all(|p| p.iter().all(|c| c.is_finite())));
    }

    #[test]
    fn rebar_cross_section_area_is_pi_r_squared() {
        // Ø16 → π·8² ≈ 201.06 mm²; Ø20 → π·10² ≈ 314.16 mm².
        assert!((rebar_cross_section_area_mm2(16.0) - 201.062).abs() < 0.01);
        assert!((rebar_cross_section_area_mm2(20.0) - 314.159).abs() < 0.01);
        // Quadruples when the diameter doubles.
        assert!(
            (rebar_cross_section_area_mm2(20.0) - 4.0 * rebar_cross_section_area_mm2(10.0)).abs()
                < 1e-9
        );
        // Guards: non-positive / non-finite → 0.
        assert_eq!(rebar_cross_section_area_mm2(0.0), 0.0);
        assert_eq!(rebar_cross_section_area_mm2(-5.0), 0.0);
        assert_eq!(rebar_cross_section_area_mm2(f64::NAN), 0.0);
    }

    #[test]
    fn reinforcement_ratio_is_steel_over_gross() {
        // Canonical: 2000 mm² steel in a 100000 mm² section → ρ = 0.02.
        assert!((reinforcement_ratio(2000.0, 100000.0) - 0.02).abs() < 1e-12);
        // Compose with the bar-area fn: 4·Ø16 bars in a 300×300 column.
        let steel = 4.0 * rebar_cross_section_area_mm2(16.0);
        let rho = reinforcement_ratio(steel, 300.0 * 300.0);
        assert!((rho - 0.008936).abs() < 1e-5);
        // Guards: gross ≤ 0 / non-finite → 0.
        assert_eq!(reinforcement_ratio(2000.0, 0.0), 0.0);
        assert_eq!(reinforcement_ratio(2000.0, -1.0), 0.0);
        assert_eq!(reinforcement_ratio(f64::NAN, 100000.0), 0.0);
    }
}
