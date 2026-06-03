//! Surface-finish symbols per ISO 1302.
//!
//! The classic check-mark triangle with optional branches:
//!
//! ```text
//!       process
//!        ___
//!        \ /  Ra
//!         V
//! ```
//!
//! The triangle itself has three forms:
//!
//! - **Required** — plain `V` (machining or other process is required
//!   but unspecified).
//! - **Machined** — `V` with a horizontal bar across the top (material
//!   removal required).
//! - **AsCast** — `V` with a small circle in the apex (material
//!   removal *not* permitted).
//! - **Removed** — `V` with an "X" through it (specific process
//!   required, defined in the tail).
//!
//! Above the triangle: roughness value (Ra in µm) + optional process /
//! lay symbol. Below: extra notes (Phase 18 v1: empty).

use serde::{Deserialize, Serialize};

/// Material-removal requirement.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SurfaceProcess {
    /// Plain V — process unspecified.
    Required,
    /// V + horizontal bar — machining (material removal) required.
    Machined,
    /// V + apex circle — material removal not permitted.
    AsCast,
    /// V + X — specific process (defined in tail).
    Removed,
}

impl SurfaceProcess {
    /// UI label.
    pub fn label(self) -> &'static str {
        match self {
            SurfaceProcess::Required => "Required",
            SurfaceProcess::Machined => "Machined",
            SurfaceProcess::AsCast => "AsCast",
            SurfaceProcess::Removed => "Removed",
        }
    }
}

/// ISO 1302 lay pattern (direction of the tool marks).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum LayPattern {
    /// `=` — parallel to the view plane.
    Parallel,
    /// `⊥` — perpendicular to the view plane.
    Perpendicular,
    /// `X` — crossed in two directions.
    Crossed,
    /// `M` — multi-directional / random.
    Multi,
    /// `R` — radial.
    Radial,
    /// `C` — circular.
    Circular,
}

impl LayPattern {
    /// Single-character lay symbol per ISO 1302.
    pub fn glyph(self) -> &'static str {
        match self {
            LayPattern::Parallel => "=",
            LayPattern::Perpendicular => "⊥",
            LayPattern::Crossed => "X",
            LayPattern::Multi => "M",
            LayPattern::Radial => "R",
            LayPattern::Circular => "C",
        }
    }
    /// UI label.
    pub fn label(self) -> &'static str {
        match self {
            LayPattern::Parallel => "Parallel",
            LayPattern::Perpendicular => "Perpendicular",
            LayPattern::Crossed => "Crossed",
            LayPattern::Multi => "Multi",
            LayPattern::Radial => "Radial",
            LayPattern::Circular => "Circular",
        }
    }
}

/// One ISO 1302 surface-finish callout.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SurfaceFinish {
    /// Apex of the triangle in sheet mm.
    pub position: [f64; 2],
    /// Roughness value (Ra) in µm. Set to 0 to omit the label.
    pub roughness_value: f64,
    /// Process requirement.
    pub process: SurfaceProcess,
    /// Lay pattern (tool-mark direction).
    pub lay_pattern: LayPattern,
}

impl SurfaceFinish {
    /// Default Ra=1.6 µm machined finish with parallel lay.
    pub fn new(position: [f64; 2], roughness_value: f64) -> Self {
        Self {
            position,
            roughness_value,
            process: SurfaceProcess::Machined,
            lay_pattern: LayPattern::Parallel,
        }
    }
}

/// Render the surface-finish symbol to an SVG fragment.
pub fn render_svg(s: &SurfaceFinish, sheet_height: f64) -> String {
    let (px, py) = (s.position[0], sheet_height - s.position[1]);
    let tri_w = 4.0;
    let tri_h = 5.0;
    let mut out = String::new();
    out.push_str("    <g class=\"surface-finish\">\n");
    // V triangle: apex at (px, py), opens upward (in flipped-Y SVG, points down).
    let left = (px - tri_w * 0.5, py - tri_h);
    let right = (px + tri_w * 0.5, py - tri_h);
    out.push_str(&format!(
        "      <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.4\"/>\n\
         <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.4\"/>\n",
        px, py, left.0, left.1, px, py, right.0, right.1,
    ));
    match s.process {
        SurfaceProcess::Required => {}
        SurfaceProcess::Machined => {
            // Horizontal bar across top of V.
            out.push_str(&format!(
                "      <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.4\"/>\n",
                left.0, left.1, right.0, right.1
            ));
        }
        SurfaceProcess::AsCast => {
            // Small circle at apex.
            out.push_str(&format!(
                "      <circle cx=\"{px}\" cy=\"{}\" r=\"0.8\" fill=\"none\" stroke=\"black\" stroke-width=\"0.3\"/>\n",
                py - tri_h * 0.5
            ));
        }
        SurfaceProcess::Removed => {
            // X across the V.
            out.push_str(&format!(
                "      <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.4\"/>\n\
                 <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.4\"/>\n",
                left.0,
                left.1,
                right.0,
                py,
                right.0,
                left.1,
                left.0,
                py,
            ));
        }
    }
    // Ra value above the V (between the branches).
    if s.roughness_value > 0.0 {
        out.push_str(&format!(
            "      <text x=\"{}\" y=\"{}\" font-size=\"2.5\" text-anchor=\"middle\">Ra {:.2}</text>\n",
            px,
            py - tri_h - 1.0,
            s.roughness_value
        ));
    }
    // Lay glyph at right of triangle.
    out.push_str(&format!(
        "      <text x=\"{}\" y=\"{}\" font-size=\"2.5\">{}</text>\n",
        right.0 + 1.0,
        left.1 + 2.0,
        s.lay_pattern.glyph()
    ));
    out.push_str("    </g>\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_to_machined_parallel_ra() {
        let s = SurfaceFinish::new([10.0, 20.0], 1.6);
        assert_eq!(s.process, SurfaceProcess::Machined);
        assert_eq!(s.lay_pattern, LayPattern::Parallel);
        assert!((s.roughness_value - 1.6).abs() < 1e-9);
    }

    #[test]
    fn render_includes_ra_value() {
        let s = SurfaceFinish::new([10.0, 20.0], 3.2);
        let svg = render_svg(&s, 100.0);
        assert!(svg.contains("Ra 3.20"));
        assert!(svg.contains("class=\"surface-finish\""));
    }

    #[test]
    fn machined_adds_top_bar() {
        let mut s = SurfaceFinish::new([0.0, 0.0], 1.6);
        s.process = SurfaceProcess::Machined;
        let svg = render_svg(&s, 100.0);
        // 2 V lines + 1 top bar + lay glyph (no circle/X).
        let n_lines = svg.matches("<line").count();
        assert_eq!(n_lines, 3);
    }

    #[test]
    fn ascast_adds_circle() {
        let mut s = SurfaceFinish::new([0.0, 0.0], 0.0);
        s.process = SurfaceProcess::AsCast;
        let svg = render_svg(&s, 100.0);
        assert!(svg.contains("<circle"));
    }

    #[test]
    fn removed_adds_x() {
        let mut s = SurfaceFinish::new([0.0, 0.0], 0.0);
        s.process = SurfaceProcess::Removed;
        let svg = render_svg(&s, 100.0);
        // 2 V lines + 2 X lines = 4 lines.
        let n_lines = svg.matches("<line").count();
        assert_eq!(n_lines, 4);
    }

    #[test]
    fn lay_pattern_glyphs_distinct() {
        let glyphs: Vec<_> = [
            LayPattern::Parallel,
            LayPattern::Perpendicular,
            LayPattern::Crossed,
            LayPattern::Multi,
            LayPattern::Radial,
            LayPattern::Circular,
        ]
        .iter()
        .map(|l| l.glyph())
        .collect();
        let set: std::collections::HashSet<_> = glyphs.iter().collect();
        assert_eq!(set.len(), glyphs.len());
    }

    #[test]
    fn ra_zero_suppresses_label() {
        let s = SurfaceFinish::new([0.0, 0.0], 0.0);
        let svg = render_svg(&s, 100.0);
        assert!(!svg.contains("Ra "));
    }

    #[test]
    fn process_labels_distinct() {
        let labels = [
            SurfaceProcess::Required.label(),
            SurfaceProcess::Machined.label(),
            SurfaceProcess::AsCast.label(),
            SurfaceProcess::Removed.label(),
        ];
        let set: std::collections::HashSet<_> = labels.iter().collect();
        assert_eq!(set.len(), labels.len());
    }
}
