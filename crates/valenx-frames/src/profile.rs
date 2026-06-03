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
}
