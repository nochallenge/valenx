//! Leader lines — a line + arrowhead pointing from a label to a part
//! feature.
//!
//! A [`Leader`] carries `start` (where the text label sits) and `end`
//! (the arrow tip on the part), an optional `text` label, and an
//! [`ArrowKind`] for the arrow style. [`polyline`] and [`with_jog`]
//! (Phase 18G) build multi-segment leaders for tight layouts.
//!
//! Coordinates are sheet millimeters. Renderers consume [`render_svg`]
//! for the SVG fragment.

use serde::{Deserialize, Serialize};

/// Arrowhead style at the end of a leader.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ArrowKind {
    /// Filled solid triangle (default for dimensions / leaders).
    Closed,
    /// Two-line open V.
    Open,
    /// Small filled dot (used for "this point" callouts).
    Dot,
    /// Short perpendicular tick (architectural convention).
    Tick,
}

impl ArrowKind {
    /// UI label.
    pub fn label(self) -> &'static str {
        match self {
            ArrowKind::Closed => "Closed",
            ArrowKind::Open => "Open",
            ArrowKind::Dot => "Dot",
            ArrowKind::Tick => "Tick",
        }
    }
}

/// A leader line with an arrowhead at `end`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Leader {
    /// Where the text label / start of the leader sits.
    pub start: [f64; 2],
    /// Arrow tip on the part feature.
    pub end: [f64; 2],
    /// Optional text label drawn near `start`.
    pub text: String,
    /// Arrowhead style at `end`.
    pub arrow_kind: ArrowKind,
}

impl Leader {
    /// Standard closed-arrow leader with `text`.
    pub fn new(start: [f64; 2], end: [f64; 2], text: impl Into<String>) -> Self {
        Self {
            start,
            end,
            text: text.into(),
            arrow_kind: ArrowKind::Closed,
        }
    }
}

/// SVG `<g>` fragment for a leader line + arrow + label.
///
/// `sheet_height` flips Y (paper origin bottom-left → SVG origin
/// top-left).
pub fn render_svg(l: &Leader, sheet_height: f64) -> String {
    let mut out = String::new();
    out.push_str("    <g class=\"leader\">\n");
    out.push_str(&line_svg(l.start, l.end, sheet_height));
    out.push_str(&arrowhead_svg(l.start, l.end, l.arrow_kind, sheet_height));
    if !l.text.is_empty() {
        let (sx, sy) = (l.start[0], sheet_height - l.start[1]);
        out.push_str(&format!(
            "      <text x=\"{}\" y=\"{}\" font-size=\"3\">{}</text>\n",
            sx + 1.0,
            sy - 1.0,
            escape_xml(&l.text)
        ));
    }
    out.push_str("    </g>\n");
    out
}

fn line_svg(a: [f64; 2], b: [f64; 2], h: f64) -> String {
    format!(
        "      <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.3\"/>\n",
        a[0],
        h - a[1],
        b[0],
        h - b[1]
    )
}

/// Arrow-tip SVG fragment. Shared by [`render_svg`] +
/// [`polyline`]/[`with_jog`] so the arrow style stays consistent.
fn arrowhead_svg(start: [f64; 2], end: [f64; 2], kind: ArrowKind, h: f64) -> String {
    let dx = end[0] - start[0];
    let dy = end[1] - start[1];
    let len = (dx * dx + dy * dy).sqrt().max(1e-9);
    let ux = dx / len;
    let uy = dy / len;
    let nx = -uy;
    let ny = ux;
    let tip = (end[0], h - end[1]);
    let arrow_len = 2.0;
    let arrow_w = 0.7;
    match kind {
        ArrowKind::Closed => {
            let base_x = end[0] - ux * arrow_len;
            let base_y = end[1] - uy * arrow_len;
            let l = (base_x - nx * arrow_w, h - (base_y - ny * arrow_w));
            let r = (base_x + nx * arrow_w, h - (base_y + ny * arrow_w));
            format!(
                "      <polygon points=\"{} {} {} {} {} {}\" fill=\"black\"/>\n",
                tip.0, tip.1, l.0, l.1, r.0, r.1
            )
        }
        ArrowKind::Open => {
            let base_x = end[0] - ux * arrow_len;
            let base_y = end[1] - uy * arrow_len;
            let l = (base_x - nx * arrow_w, h - (base_y - ny * arrow_w));
            let r = (base_x + nx * arrow_w, h - (base_y + ny * arrow_w));
            format!(
                "      <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.3\"/>\n\
                 <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.3\"/>\n",
                tip.0, tip.1, l.0, l.1, tip.0, tip.1, r.0, r.1
            )
        }
        ArrowKind::Dot => {
            format!(
                "      <circle cx=\"{}\" cy=\"{}\" r=\"0.6\" fill=\"black\"/>\n",
                tip.0, tip.1
            )
        }
        ArrowKind::Tick => {
            let tw = 1.2;
            let l = (end[0] - nx * tw, h - (end[1] - ny * tw));
            let r = (end[0] + nx * tw, h - (end[1] + ny * tw));
            format!(
                "      <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.4\"/>\n",
                l.0, l.1, r.0, r.1
            )
        }
    }
}

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// Multi-segment leader (Phase 18G). Polyline through `points` with
/// `arrow_kind` at the last vertex. Returns the SVG fragment.
pub fn polyline(points: &[[f64; 2]], arrow_kind: ArrowKind, sheet_height: f64) -> String {
    if points.len() < 2 {
        return String::new();
    }
    let mut out = String::new();
    out.push_str("    <g class=\"leader polyline\">\n");
    for w in points.windows(2) {
        out.push_str(&line_svg(w[0], w[1], sheet_height));
    }
    let last_two = (points[points.len() - 2], points[points.len() - 1]);
    out.push_str(&arrowhead_svg(
        last_two.0,
        last_two.1,
        arrow_kind,
        sheet_height,
    ));
    out.push_str("    </g>\n");
    out
}

/// Leader with a horizontal jog (Phase 18G). Starts at `start`, runs
/// to `(end.x, start.y + jog_offset)`, then drops down (or up) to
/// `end`. The arrow tip is at `end`. Returns the SVG fragment.
pub fn with_jog(
    start: [f64; 2],
    end: [f64; 2],
    jog_offset: f64,
    arrow_kind: ArrowKind,
    sheet_height: f64,
) -> String {
    let knee = [end[0], start[1] + jog_offset];
    polyline(&[start, knee, end], arrow_kind, sheet_height)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_leader_uses_closed_arrow() {
        let l = Leader::new([0.0, 0.0], [10.0, 0.0], "A");
        assert_eq!(l.arrow_kind, ArrowKind::Closed);
        assert_eq!(l.text, "A");
    }

    #[test]
    fn render_closed_uses_filled_polygon() {
        let l = Leader::new([0.0, 0.0], [10.0, 0.0], "ABC");
        let svg = render_svg(&l, 100.0);
        assert!(svg.contains("<polygon"));
        assert!(svg.contains("fill=\"black\""));
        assert!(svg.contains(">ABC<"));
    }

    #[test]
    fn render_dot_uses_circle() {
        let mut l = Leader::new([0.0, 0.0], [10.0, 0.0], "");
        l.arrow_kind = ArrowKind::Dot;
        let svg = render_svg(&l, 100.0);
        assert!(svg.contains("<circle"));
    }

    #[test]
    fn polyline_three_points_emits_two_segments() {
        let svg = polyline(
            &[[0.0, 0.0], [5.0, 0.0], [10.0, 0.0]],
            ArrowKind::Open,
            100.0,
        );
        // two segments + one arrowhead which is also two lines (open).
        // Count <line> occurrences in SVG fragment.
        let lines = svg.matches("<line").count();
        // 2 polyline segments + 2 arrowhead lines = 4.
        assert_eq!(lines, 4);
    }

    #[test]
    fn polyline_too_short_returns_empty() {
        let svg = polyline(&[[0.0, 0.0]], ArrowKind::Closed, 100.0);
        assert!(svg.is_empty());
    }

    #[test]
    fn with_jog_yields_three_points() {
        let svg = with_jog([0.0, 0.0], [20.0, 10.0], 5.0, ArrowKind::Closed, 100.0);
        // three points → two segments → two <line>.
        let lines = svg.matches("<line").count();
        assert_eq!(lines, 2);
        assert!(svg.contains("<polygon"));
    }

    #[test]
    fn arrow_kind_labels_distinct() {
        let labels = [
            ArrowKind::Closed.label(),
            ArrowKind::Open.label(),
            ArrowKind::Dot.label(),
            ArrowKind::Tick.label(),
        ];
        let set: std::collections::HashSet<_> = labels.iter().collect();
        assert_eq!(set.len(), labels.len());
    }
}
