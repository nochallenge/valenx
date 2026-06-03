//! Weld symbols per ISO 2553.
//!
//! A welding symbol has three layout components:
//!
//! 1. **Reference line** — a horizontal line that carries the symbol.
//! 2. **Arrow** — a leader pointing from the reference line at the
//!    joint. Symbols above the line apply to the "other side"
//!    (opposite the arrow); symbols below apply to the "arrow side."
//! 3. **Tail** — optional flag at the far end of the reference line
//!    carrying process notes (E70xx, GMAW, etc.). Phase 18 v1 keeps
//!    the tail empty.
//!
//! The symbol itself is a small graphic per [`WeldType`]:
//!
//! - **Fillet** — right-triangle.
//! - **Square** — vertical bar.
//! - **V** — two slanted lines forming a "V".
//! - **U** — half-circle dished upward.
//! - **Bevel** — single slanted line.
//! - **J** — vertical line + quarter-circle.
//! - **Flare** — curved bevel.
//! - **Plug** — small filled square.
//! - **Seam** — long horizontal rectangle.
//!
//! See the canonical ISO 2553 weld-symbol reference (or Wikipedia's
//! "ISO 2553" article) for the exact pictograms.

use serde::{Deserialize, Serialize};

/// Which side of the reference line the symbol belongs on.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum WeldPosition {
    /// Below the reference line — applies to the side the arrow points at.
    Arrow,
    /// Above the reference line — applies to the opposite side.
    Other,
}

impl WeldPosition {
    /// UI label.
    pub fn label(self) -> &'static str {
        match self {
            WeldPosition::Arrow => "Arrow",
            WeldPosition::Other => "Other",
        }
    }
}

/// One of the ISO 2553 weld pictograms.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum WeldType {
    /// Right-triangle fillet weld (the most common pictogram).
    Fillet,
    /// Square-groove weld (vertical bar).
    Square,
    /// V-groove weld (two slanted lines).
    V,
    /// U-groove weld (half-circle).
    U,
    /// Bevel weld (single slanted line).
    Bevel,
    /// J-groove weld (vertical line + quarter-circle).
    J,
    /// Flare-V weld (curved).
    Flare,
    /// Plug weld (small filled square).
    Plug,
    /// Seam weld (long rectangle).
    Seam,
}

impl WeldType {
    /// UI label.
    pub fn label(self) -> &'static str {
        match self {
            WeldType::Fillet => "Fillet",
            WeldType::Square => "Square",
            WeldType::V => "V",
            WeldType::U => "U",
            WeldType::Bevel => "Bevel",
            WeldType::J => "J",
            WeldType::Flare => "Flare",
            WeldType::Plug => "Plug",
            WeldType::Seam => "Seam",
        }
    }
}

/// Pictogram SVG path (without leading `<path>` wrapper) drawn around
/// a unit square at the origin, scaled to ~3 mm tall.
pub fn weld_type_glyph(t: WeldType) -> &'static str {
    match t {
        WeldType::Fillet => "M0 0 L3 0 L3 3 Z",
        WeldType::Square => "M1.5 0 L1.5 3",
        WeldType::V => "M0 3 L1.5 0 L3 3",
        WeldType::U => "M0 3 A1.5 1.5 0 0 1 3 3",
        WeldType::Bevel => "M0 3 L1.5 0",
        WeldType::J => "M0 3 L0 1 A1 1 0 0 1 1 0",
        WeldType::Flare => "M0 3 Q1.5 -1 3 3",
        WeldType::Plug => "M0.5 0.5 L2.5 0.5 L2.5 2.5 L0.5 2.5 Z",
        WeldType::Seam => "M0 1 L4 1 L4 2 L0 2 Z",
    }
}

/// A single weld-symbol annotation.
///
/// `position_2d` is the sheet-mm origin of the reference line's left
/// end. `arrow_target` is where the leader arrowhead lands. `size`,
/// `length`, and `pitch` annotate the weld dimensions; an empty
/// string (default) suppresses that part of the label.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct WeldSymbol {
    /// Left end of the reference line (sheet mm).
    pub position: [f64; 2],
    /// Where the leader arrow points (sheet mm).
    pub arrow_target: [f64; 2],
    /// Whether the pictogram sits above (Other) or below (Arrow) the
    /// reference line.
    pub weld_position: WeldPosition,
    /// Pictogram type.
    pub weld_type: WeldType,
    /// Weld throat / size label.
    pub size: String,
    /// Weld length label.
    pub length: String,
    /// Intermittent-weld pitch.
    pub pitch: String,
    /// "Weld in the field" flag — adds a black filled flag at the
    /// arrow / reference junction.
    pub field_weld: bool,
    /// "Weld all around" flag — adds an open circle at the junction.
    pub all_around: bool,
}

impl WeldSymbol {
    /// Standard 8 mm fillet on the arrow side.
    pub fn new_fillet(position: [f64; 2], arrow_target: [f64; 2], size: impl Into<String>) -> Self {
        Self {
            position,
            arrow_target,
            weld_position: WeldPosition::Arrow,
            weld_type: WeldType::Fillet,
            size: size.into(),
            length: String::new(),
            pitch: String::new(),
            field_weld: false,
            all_around: false,
        }
    }
}

/// Render the weld symbol as an SVG fragment.
///
/// Layout:
/// ```text
///    ┌─symbol (above when Other)
///    ├─────────────── reference line ──────────
///    └─symbol (below when Arrow)
///    │
///    arrow ↓ to arrow_target
/// ```
///
/// `sheet_height` flips Y for SVG.
pub fn render_svg(s: &WeldSymbol, sheet_height: f64) -> String {
    let (px, py) = (s.position[0], sheet_height - s.position[1]);
    let ref_len = 14.0;
    let mut out = String::new();
    out.push_str("    <g class=\"weld-symbol\">\n");
    // Reference line.
    out.push_str(&format!(
        "      <line x1=\"{px}\" y1=\"{py}\" x2=\"{}\" y2=\"{py}\" stroke=\"black\" stroke-width=\"0.4\"/>\n",
        px + ref_len
    ));
    // Tail (Phase 18 v1: simple two short lines forming an open "<").
    let tx = px + ref_len;
    out.push_str(&format!(
        "      <line x1=\"{tx}\" y1=\"{py}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.3\"/>\n\
         <line x1=\"{tx}\" y1=\"{py}\" x2=\"{}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.3\"/>\n",
        tx + 3.0,
        py - 2.0,
        tx + 3.0,
        py + 2.0
    ));
    // Arrow leader from position to arrow_target.
    let (ax, ay) = (s.arrow_target[0], sheet_height - s.arrow_target[1]);
    out.push_str(&format!(
        "      <line x1=\"{px}\" y1=\"{py}\" x2=\"{ax}\" y2=\"{ay}\" stroke=\"black\" stroke-width=\"0.3\"/>\n"
    ));
    // Arrowhead (filled triangle, same convention as Leader::Closed).
    let dx = ax - px;
    let dy = ay - py;
    let len = (dx * dx + dy * dy).sqrt().max(1e-9);
    let ux = dx / len;
    let uy = dy / len;
    let nx = -uy;
    let ny = ux;
    let alen = 2.0;
    let aw = 0.7;
    let bx = ax - ux * alen;
    let by = ay - uy * alen;
    out.push_str(&format!(
        "      <polygon points=\"{} {} {} {} {} {}\" fill=\"black\"/>\n",
        ax,
        ay,
        bx - nx * aw,
        by - ny * aw,
        bx + nx * aw,
        by + ny * aw,
    ));
    // Pictogram. Above reference line for Other, below for Arrow.
    let glyph_offset_y = match s.weld_position {
        WeldPosition::Other => -4.0,
        WeldPosition::Arrow => 1.0,
    };
    let glyph_tx = px + 4.0;
    let glyph_ty = py + glyph_offset_y;
    out.push_str(&format!(
        "      <g transform=\"translate({glyph_tx} {glyph_ty})\">\n\
         <path d=\"{}\" fill=\"none\" stroke=\"black\" stroke-width=\"0.4\"/>\n\
         </g>\n",
        weld_type_glyph(s.weld_type)
    ));
    // Optional flags.
    if s.all_around {
        out.push_str(&format!(
            "      <circle cx=\"{px}\" cy=\"{py}\" r=\"1.5\" fill=\"none\" stroke=\"black\" stroke-width=\"0.3\"/>\n"
        ));
    }
    if s.field_weld {
        // Black flag attached at junction.
        out.push_str(&format!(
            "      <polygon points=\"{} {} {} {} {} {}\" fill=\"black\"/>\n",
            px,
            py - 1.0,
            px + 3.0,
            py - 1.0,
            px,
            py - 3.5,
        ));
    }
    // Dimension labels: size before reference, length after, pitch as suffix.
    if !s.size.is_empty() {
        out.push_str(&format!(
            "      <text x=\"{}\" y=\"{}\" font-size=\"2.5\">{}</text>\n",
            px,
            py + glyph_offset_y + 4.5,
            escape_xml(&s.size)
        ));
    }
    if !s.length.is_empty() {
        let label = if s.pitch.is_empty() {
            s.length.clone()
        } else {
            format!("{}-{}", s.length, s.pitch)
        };
        out.push_str(&format!(
            "      <text x=\"{}\" y=\"{}\" font-size=\"2.5\">{}</text>\n",
            px + 8.0,
            py + glyph_offset_y + 4.5,
            escape_xml(&label)
        ));
    }
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
    fn fillet_constructor_sets_defaults() {
        let w = WeldSymbol::new_fillet([0.0, 0.0], [10.0, 10.0], "8");
        assert_eq!(w.weld_type, WeldType::Fillet);
        assert_eq!(w.weld_position, WeldPosition::Arrow);
        assert_eq!(w.size, "8");
        assert!(!w.field_weld);
        assert!(!w.all_around);
    }

    #[test]
    fn render_fillet_includes_reference_line_arrow_and_glyph() {
        let w = WeldSymbol::new_fillet([10.0, 50.0], [30.0, 70.0], "6");
        let svg = render_svg(&w, 200.0);
        assert!(svg.contains("class=\"weld-symbol\""));
        assert!(svg.contains("<polygon")); // arrowhead
        assert!(svg.contains("<path")); // glyph
        assert!(svg.contains(">6<")); // size label
    }

    #[test]
    fn all_around_adds_circle() {
        let mut w = WeldSymbol::new_fillet([0.0, 0.0], [10.0, 10.0], "5");
        w.all_around = true;
        let svg = render_svg(&w, 100.0);
        assert!(svg.contains("<circle"));
    }

    #[test]
    fn field_weld_adds_filled_flag() {
        let mut w = WeldSymbol::new_fillet([0.0, 0.0], [10.0, 10.0], "5");
        w.field_weld = true;
        let svg = render_svg(&w, 100.0);
        // The flag is a black triangle — we should see at least 2 fill="black" polygons (arrow + flag).
        assert!(svg.matches("fill=\"black\"").count() >= 2);
    }

    #[test]
    fn glyph_changes_per_weld_type() {
        let glyph_a = weld_type_glyph(WeldType::Fillet);
        let glyph_b = weld_type_glyph(WeldType::V);
        let glyph_c = weld_type_glyph(WeldType::U);
        assert_ne!(glyph_a, glyph_b);
        assert_ne!(glyph_b, glyph_c);
    }

    #[test]
    fn position_above_vs_below_changes_glyph_y() {
        let mut w_above = WeldSymbol::new_fillet([0.0, 0.0], [10.0, 10.0], "5");
        w_above.weld_position = WeldPosition::Other;
        let mut w_below = WeldSymbol::new_fillet([0.0, 0.0], [10.0, 10.0], "5");
        w_below.weld_position = WeldPosition::Arrow;
        let svg_a = render_svg(&w_above, 100.0);
        let svg_b = render_svg(&w_below, 100.0);
        // glyph transform line differs.
        let g_a = svg_a
            .split("translate(")
            .nth(1)
            .unwrap_or("")
            .split(')')
            .next()
            .unwrap_or("");
        let g_b = svg_b
            .split("translate(")
            .nth(1)
            .unwrap_or("")
            .split(')')
            .next()
            .unwrap_or("");
        assert_ne!(g_a, g_b);
    }

    #[test]
    fn weld_type_labels_unique() {
        let labels = [
            WeldType::Fillet.label(),
            WeldType::Square.label(),
            WeldType::V.label(),
            WeldType::U.label(),
            WeldType::Bevel.label(),
            WeldType::J.label(),
            WeldType::Flare.label(),
            WeldType::Plug.label(),
            WeldType::Seam.label(),
        ];
        let set: std::collections::HashSet<_> = labels.iter().collect();
        assert_eq!(set.len(), labels.len());
    }
}
