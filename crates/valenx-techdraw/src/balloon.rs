//! Balloons — circular / polygonal item callouts pointing at a part.
//!
//! A [`Balloon`] is placed at `position` (where the number text sits)
//! and points back to `target_point` via a thin leader. The `number`
//! is rendered inside the bubble; `style` controls the bubble shape
//! per ASME Y14.34 (Circle is the universal default, the others show
//! up in process-flow drawings).
//!
//! The geometry is pure data; [`render_svg`] turns a balloon into an
//! SVG `<g>` fragment the export pipeline can splice into the
//! drawing. The fragment uses sheet-relative coordinates so the
//! exporter doesn't have to apply any further transform.

use serde::{Deserialize, Serialize};

/// Bubble shape for a balloon.
///
/// Circle is the universal default. Square/Triangle/Hexagon are used
/// in some process-engineering and aerospace drawing styles.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BalloonStyle {
    /// Round bubble (most common, ASME Y14.34 default).
    Circle,
    /// Square bubble.
    Square,
    /// Six-sided bubble.
    Hexagon,
    /// Triangle bubble.
    Triangle,
}

impl BalloonStyle {
    /// UI label.
    pub fn label(self) -> &'static str {
        match self {
            BalloonStyle::Circle => "Circle",
            BalloonStyle::Square => "Square",
            BalloonStyle::Hexagon => "Hexagon",
            BalloonStyle::Triangle => "Triangle",
        }
    }
}

/// One balloon callout.
///
/// All coordinates are in **sheet millimeters**. `position` is where
/// the bubble's center lands; the renderer draws the leader from the
/// edge of the bubble to `target_point`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Balloon {
    /// Bubble center (sheet mm).
    pub position: [f64; 2],
    /// Item number / letter inside the bubble.
    pub number: String,
    /// Point on the part the leader points at.
    pub target_point: [f64; 2],
    /// Bubble shape.
    pub style: BalloonStyle,
    /// Radius of the bubble outline in mm. Defaults to ~4.0 (sized to
    /// fit a 1–2-digit number at 3 mm font).
    pub radius: f64,
}

impl Balloon {
    /// Construct a standard 4 mm Circle balloon.
    pub fn new(position: [f64; 2], number: impl Into<String>, target_point: [f64; 2]) -> Self {
        Self {
            position,
            number: number.into(),
            target_point,
            style: BalloonStyle::Circle,
            radius: 4.0,
        }
    }
}

/// Render a balloon to an SVG fragment using sheet-mm coordinates.
///
/// `sheet_height` is needed because SVG's Y axis points down while
/// our drawing coordinates assume Y-up (bottom-left origin, matching
/// paper conventions). Pass `drawing.sheet.dimensions_mm().1`.
pub fn render_svg(b: &Balloon, sheet_height: f64) -> String {
    let cx = b.position[0];
    let cy = sheet_height - b.position[1];
    let tx = b.target_point[0];
    let ty = sheet_height - b.target_point[1];
    let r = b.radius.max(0.5);

    // Leader: from balloon edge (toward target) to target_point.
    let dx = tx - cx;
    let dy = ty - cy;
    let len = (dx * dx + dy * dy).sqrt().max(1e-9);
    let ex = cx + dx / len * r;
    let ey = cy + dy / len * r;

    let mut out = String::new();
    out.push_str("    <g class=\"balloon\">\n");
    out.push_str(&format!(
        "      <line x1=\"{ex}\" y1=\"{ey}\" x2=\"{tx}\" y2=\"{ty}\" stroke=\"black\" stroke-width=\"0.3\"/>\n"
    ));
    match b.style {
        BalloonStyle::Circle => out.push_str(&format!(
            "      <circle cx=\"{cx}\" cy=\"{cy}\" r=\"{r}\" fill=\"white\" stroke=\"black\" stroke-width=\"0.3\"/>\n"
        )),
        BalloonStyle::Square => out.push_str(&format!(
            "      <rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"white\" stroke=\"black\" stroke-width=\"0.3\"/>\n",
            cx - r,
            cy - r,
            r * 2.0,
            r * 2.0
        )),
        BalloonStyle::Triangle => {
            let p1 = (cx, cy - r);
            let p2 = (cx - r, cy + r * 0.866);
            let p3 = (cx + r, cy + r * 0.866);
            out.push_str(&format!(
                "      <polygon points=\"{} {} {} {} {} {}\" fill=\"white\" stroke=\"black\" stroke-width=\"0.3\"/>\n",
                p1.0, p1.1, p2.0, p2.1, p3.0, p3.1
            ));
        }
        BalloonStyle::Hexagon => {
            let mut pts = String::new();
            for i in 0..6 {
                let a = (i as f64) * std::f64::consts::PI / 3.0;
                pts.push_str(&format!("{} {} ", cx + r * a.cos(), cy + r * a.sin()));
            }
            out.push_str(&format!(
                "      <polygon points=\"{}\" fill=\"white\" stroke=\"black\" stroke-width=\"0.3\"/>\n",
                pts.trim()
            ));
        }
    }
    out.push_str(&format!(
        "      <text x=\"{cx}\" y=\"{}\" font-size=\"3\" text-anchor=\"middle\">{}</text>\n",
        cy + 1.0,
        escape_xml(&b.number)
    ));
    out.push_str("    </g>\n");
    out
}

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn balloon_defaults_to_circle_4mm() {
        let b = Balloon::new([10.0, 20.0], "5", [30.0, 40.0]);
        assert_eq!(b.style, BalloonStyle::Circle);
        assert!((b.radius - 4.0).abs() < 1e-9);
        assert_eq!(b.number, "5");
    }

    #[test]
    fn render_circle_includes_circle_and_text() {
        let b = Balloon::new([10.0, 20.0], "7", [30.0, 40.0]);
        let svg = render_svg(&b, 200.0);
        assert!(svg.contains("<circle"));
        assert!(svg.contains(">7<"));
        assert!(svg.contains("class=\"balloon\""));
    }

    #[test]
    fn render_square_uses_rect_element() {
        let mut b = Balloon::new([0.0, 0.0], "1", [10.0, 0.0]);
        b.style = BalloonStyle::Square;
        let svg = render_svg(&b, 100.0);
        assert!(svg.contains("<rect"));
    }

    #[test]
    fn render_hexagon_uses_polygon_with_six_points() {
        let mut b = Balloon::new([0.0, 0.0], "1", [10.0, 0.0]);
        b.style = BalloonStyle::Hexagon;
        let svg = render_svg(&b, 100.0);
        assert!(svg.contains("<polygon"));
        // 6 "x y " pairs.
        let pts_attr = svg
            .split("points=\"")
            .nth(1)
            .unwrap_or("")
            .split('"')
            .next()
            .unwrap_or("");
        assert_eq!(pts_attr.split_whitespace().count(), 12);
    }

    #[test]
    fn style_labels_distinct() {
        let labels = [
            BalloonStyle::Circle.label(),
            BalloonStyle::Square.label(),
            BalloonStyle::Hexagon.label(),
            BalloonStyle::Triangle.label(),
        ];
        let set: std::collections::HashSet<_> = labels.iter().collect();
        assert_eq!(set.len(), labels.len());
    }

    #[test]
    fn xml_escape_handles_special_chars() {
        let b = Balloon::new([0.0, 0.0], "<&>", [10.0, 0.0]);
        let svg = render_svg(&b, 100.0);
        assert!(svg.contains("&lt;&amp;&gt;"));
    }
}
