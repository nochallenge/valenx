//! Cross-section profile catalogue. Each variant is a parametric
//! 2D shape; [`cross_section_polygon`] returns the section's outline
//! polygon as a closed CCW polygon in the local (u, v) frame.

use serde::{Deserialize, Serialize};

/// Standard structural profiles.
///
/// All dimensions in millimetres. The profile is positioned with its
/// centroid at the origin and the strong-axis ("h") aligned with the
/// local +v axis; the local +u axis points to the right.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Profile {
    /// I-beam (wide-flange) — `h` overall height, `b` flange width,
    /// `tw` web thickness, `tf` flange thickness.
    IBeam {
        /// Overall height (mm).
        h: f64,
        /// Flange width (mm).
        b: f64,
        /// Web thickness (mm).
        tw: f64,
        /// Flange thickness (mm).
        tf: f64,
    },
    /// C-channel — `h` overall height, `b` flange width, `tw` web
    /// thickness (we approximate flange thickness == tw).
    CChannel {
        /// Overall height (mm).
        h: f64,
        /// Flange width (mm).
        b: f64,
        /// Web thickness (mm).
        tw: f64,
    },
    /// L-angle (equal or unequal) — `h` vertical leg height, `b`
    /// horizontal leg width, `t` leg thickness.
    LAngle {
        /// Vertical leg length (mm).
        h: f64,
        /// Horizontal leg length (mm).
        b: f64,
        /// Leg thickness (mm).
        t: f64,
    },
    /// Rectangular hollow section (HSS rect) — `h` outer height, `b`
    /// outer width, `t` wall thickness.
    RhsRect {
        /// Outer height (mm).
        h: f64,
        /// Outer width (mm).
        b: f64,
        /// Wall thickness (mm).
        t: f64,
    },
    /// Round hollow section (CHS) — `d` outer diameter, `t` wall
    /// thickness.
    ChsRound {
        /// Outer diameter (mm).
        d: f64,
        /// Wall thickness (mm).
        t: f64,
    },
    /// T-beam — `h` overall height, `b` flange width, `tw` web
    /// thickness, `tf` flange thickness.
    TBeam {
        /// Overall height (mm).
        h: f64,
        /// Flange width (mm).
        b: f64,
        /// Web thickness (mm).
        tw: f64,
        /// Flange thickness (mm).
        tf: f64,
    },
}

impl Profile {
    /// Short UI label.
    pub fn label(self) -> &'static str {
        match self {
            Self::IBeam { .. } => "I-beam",
            Self::CChannel { .. } => "C-channel",
            Self::LAngle { .. } => "L-angle",
            Self::RhsRect { .. } => "Rect HSS",
            Self::ChsRound { .. } => "Round CHS",
            Self::TBeam { .. } => "T-beam",
        }
    }

    /// Convenient default for the UI — an IPE 200 wide-flange.
    pub fn default_ipe200() -> Self {
        Self::IBeam {
            h: 200.0,
            b: 100.0,
            tw: 5.6,
            tf: 8.5,
        }
    }

    /// Exact cross-sectional area (mm²) of the idealized (sharp-cornered, no-fillet) profile —
    /// a sum of rectangles, or an annulus for [`ChsRound`](Profile::ChsRound). Matches the
    /// geometry of [`cross_section_polygon`]; differs from a real rolled section that has root
    /// fillets. Non-finite or degenerate dimensions clamp the result to `0.0`.
    pub fn cross_section_area_mm2(self) -> f64 {
        // Inner/remaining dimensions are clamped to ≥0 so an over-thick wall collapses to a solid
        // section rather than producing a spurious negative-squared hole.
        let area = match self {
            Self::IBeam { h, b, tw, tf } => 2.0 * b * tf + tw * (h - 2.0 * tf).max(0.0),
            Self::CChannel { h, b, tw } => 2.0 * b * tw + tw * (h - 2.0 * tw).max(0.0),
            Self::LAngle { h, b, t } => t * h + t * (b - t).max(0.0),
            Self::RhsRect { h, b, t } => b * h - (b - 2.0 * t).max(0.0) * (h - 2.0 * t).max(0.0),
            Self::ChsRound { d, t } => {
                let inner = (d - 2.0 * t).max(0.0);
                std::f64::consts::PI / 4.0 * (d * d - inner * inner)
            }
            Self::TBeam { h, b, tw, tf } => b * tf + tw * (h - tf).max(0.0),
        };
        if area.is_finite() && area > 0.0 {
            area
        } else {
            0.0
        }
    }

    /// Exact outer-boundary perimeter (mm) of the idealized (sharp-cornered, no-fillet) profile —
    /// the sum of the outer outline's edge lengths; for [`ChsRound`](Profile::ChsRound) it is the
    /// circle perimeter `π·d` (not a polygon approximation). The wall-thickness terms cancel in the
    /// open profiles, so this depends only on the outer dimensions. Non-finite/degenerate → `0.0`.
    pub fn cross_section_perimeter_mm(self) -> f64 {
        let perimeter = match self {
            Self::IBeam { h, b, tw, .. } => 2.0 * (2.0 * b + h - tw),
            Self::CChannel { h, b, tw } => 2.0 * (2.0 * b + h - tw),
            Self::LAngle { h, b, .. } => 2.0 * (b + h),
            Self::RhsRect { h, b, .. } => 2.0 * (b + h),
            Self::ChsRound { d, .. } => std::f64::consts::PI * d,
            Self::TBeam { h, b, .. } => 2.0 * (b + h),
        };
        if perimeter.is_finite() && perimeter > 0.0 {
            perimeter
        } else {
            0.0
        }
    }
}

impl Default for Profile {
    fn default() -> Self {
        Self::default_ipe200()
    }
}

/// Build the closed CCW outline polygon for `p` in local (u, v).
/// The polygon is **not** explicitly closed (last vertex != first);
/// the consumer must close it if it needs an explicit duplicate.
pub fn cross_section_polygon(p: Profile) -> Vec<[f64; 2]> {
    match p {
        Profile::IBeam { h, b, tw, tf } => {
            // 12-vertex I-beam, centroid at origin.
            let h2 = h * 0.5;
            let b2 = b * 0.5;
            let tw2 = tw * 0.5;
            let yfb = h2 - tf;
            vec![
                [-b2, -h2],
                [b2, -h2],
                [b2, -yfb],
                [tw2, -yfb],
                [tw2, yfb],
                [b2, yfb],
                [b2, h2],
                [-b2, h2],
                [-b2, yfb],
                [-tw2, yfb],
                [-tw2, -yfb],
                [-b2, -yfb],
            ]
        }
        Profile::CChannel { h, b, tw } => {
            // 8-vertex C, centroid roughly at origin (-b/4 offset
            // ignored for v1 — kept axis-aligned for predictability).
            let h2 = h * 0.5;
            let b2 = b * 0.5;
            vec![
                [-b2, -h2],
                [b2, -h2],
                [b2, -h2 + tw],
                [-b2 + tw, -h2 + tw],
                [-b2 + tw, h2 - tw],
                [b2, h2 - tw],
                [b2, h2],
                [-b2, h2],
            ]
        }
        Profile::LAngle { h, b, t } => {
            // 6-vertex L with vertical leg height `h` and horizontal
            // leg width `b`; origin at the lower-left corner of the
            // bounding box, then re-centred on the leg-thickness
            // intersection.
            vec![
                [0.0, 0.0],
                [b, 0.0],
                [b, t],
                [t, t],
                [t, h],
                [0.0, h],
            ]
        }
        Profile::RhsRect { h, b, t } => {
            // Outer rectangle CCW. The hollow interior is not yet
            // represented — v1 emits the solid outer outline only
            // (the visual proxy under a sweep is still recognisable
            // as a rect HSS section).
            let _ = t;
            let h2 = h * 0.5;
            let b2 = b * 0.5;
            vec![[-b2, -h2], [b2, -h2], [b2, h2], [-b2, h2]]
        }
        Profile::ChsRound { d, t } => {
            let _ = t;
            let r = d * 0.5;
            let n = 24;
            (0..n)
                .map(|i| {
                    let th = (i as f64 / n as f64) * std::f64::consts::TAU;
                    [r * th.cos(), r * th.sin()]
                })
                .collect()
        }
        Profile::TBeam { h, b, tw, tf } => {
            // 8-vertex T, centroid at origin.
            let b2 = b * 0.5;
            let tw2 = tw * 0.5;
            let v_top = h - tf;
            vec![
                [-tw2, 0.0],
                [tw2, 0.0],
                [tw2, v_top],
                [b2, v_top],
                [b2, h],
                [-b2, h],
                [-b2, v_top],
                [-tw2, v_top],
            ]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ibeam_polygon_has_12_verts() {
        let p = Profile::IBeam {
            h: 200.0,
            b: 100.0,
            tw: 6.0,
            tf: 10.0,
        };
        assert_eq!(cross_section_polygon(p).len(), 12);
    }

    #[test]
    fn round_chs_polygon_has_24_verts() {
        let p = Profile::ChsRound { d: 100.0, t: 4.0 };
        assert_eq!(cross_section_polygon(p).len(), 24);
    }

    #[test]
    fn rect_hss_polygon_is_a_quad() {
        let p = Profile::RhsRect {
            h: 80.0,
            b: 40.0,
            t: 3.0,
        };
        assert_eq!(cross_section_polygon(p).len(), 4);
    }

    #[test]
    fn cross_section_area_exact_per_variant() {
        // I-beam {200,100,5.6,8.5}: 2·100·8.5 + 5.6·(200−17) = 2724.8.
        assert!(
            (Profile::IBeam { h: 200.0, b: 100.0, tw: 5.6, tf: 8.5 }.cross_section_area_mm2()
                - 2724.8)
                .abs()
                < 1e-9
        );
        // RHS rect {100,50,5}: 50·100 − 40·90 = 1400.
        assert!(
            (Profile::RhsRect { h: 100.0, b: 50.0, t: 5.0 }.cross_section_area_mm2() - 1400.0)
                .abs()
                < 1e-9
        );
        // CHS round {100,4}: π/4·(100²−92²) ≈ 1206.37.
        assert!(
            (Profile::ChsRound { d: 100.0, t: 4.0 }.cross_section_area_mm2()
                - std::f64::consts::PI / 4.0 * (100.0_f64 * 100.0 - 92.0 * 92.0))
                .abs()
                < 1e-9
        );
        // L-angle {100,80,5}: 5·100 + 5·75 = 875. T-beam {200,100,8,12}: 1200 + 1504 = 2704.
        assert!(
            (Profile::LAngle { h: 100.0, b: 80.0, t: 5.0 }.cross_section_area_mm2() - 875.0).abs()
                < 1e-9
        );
        assert!(
            (Profile::TBeam { h: 200.0, b: 100.0, tw: 8.0, tf: 12.0 }.cross_section_area_mm2()
                - 2704.0)
                .abs()
                < 1e-9
        );
        // An over-thick wall collapses to a solid section (no spurious negative hole).
        assert_eq!(
            Profile::RhsRect { h: 100.0, b: 100.0, t: 60.0 }.cross_section_area_mm2(),
            10000.0
        );
    }

    #[test]
    fn cross_section_perimeter_exact_per_variant() {
        // I-beam {200,100,5.6,8.5}: 2·(2·100+200−5.6) = 788.8 (tf cancels).
        assert!(
            (Profile::IBeam { h: 200.0, b: 100.0, tw: 5.6, tf: 8.5 }.cross_section_perimeter_mm()
                - 788.8)
                .abs()
                < 1e-9
        );
        // C-channel {200,100,5.6}: 2·(2·100+200−5.6) = 788.8.
        assert!(
            (Profile::CChannel { h: 200.0, b: 100.0, tw: 5.6 }.cross_section_perimeter_mm() - 788.8)
                .abs()
                < 1e-9
        );
        // L-angle {100,80,5}: 2·(80+100) = 360; RHS rect {100,50,5}: 2·(50+100) = 300.
        assert!(
            (Profile::LAngle { h: 100.0, b: 80.0, t: 5.0 }.cross_section_perimeter_mm() - 360.0)
                .abs()
                < 1e-9
        );
        assert!(
            (Profile::RhsRect { h: 100.0, b: 50.0, t: 5.0 }.cross_section_perimeter_mm() - 300.0)
                .abs()
                < 1e-9
        );
        // CHS round {100,4}: the exact circle perimeter π·d ≈ 314.16 (not the 24-gon).
        assert!(
            (Profile::ChsRound { d: 100.0, t: 4.0 }.cross_section_perimeter_mm()
                - std::f64::consts::PI * 100.0)
                .abs()
                < 1e-9
        );
        // T-beam {200,100,8,12}: 2·(100+200) = 600.
        assert!(
            (Profile::TBeam { h: 200.0, b: 100.0, tw: 8.0, tf: 12.0 }.cross_section_perimeter_mm()
                - 600.0)
                .abs()
                < 1e-9
        );
    }
}
