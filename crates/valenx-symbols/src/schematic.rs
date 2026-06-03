//! Schematic = placed symbols + wires + a renderer to SVG.

use serde::{Deserialize, Serialize};

use crate::error::SymbolError;
use crate::symbol::SymbolKind;

/// A symbol placed at a 2D position with a rotation in degrees.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlacedSymbol {
    /// Which glyph.
    pub kind: SymbolKind,
    /// Centre point in schematic units (mm).
    pub position: [f64; 2],
    /// Rotation about the centre, in degrees.
    pub rotation_deg: f64,
    /// Optional reference designator displayed below the glyph
    /// (e.g. `"R1"`).
    pub designator: String,
}

impl PlacedSymbol {
    /// Convenience constructor.
    pub fn new(kind: SymbolKind, position: [f64; 2]) -> Self {
        Self {
            kind,
            position,
            rotation_deg: 0.0,
            designator: String::new(),
        }
    }
}

/// A schematic wire — a polyline plus a net-name label.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Wire {
    /// Vertices in schematic units (mm).
    pub polyline: Vec<[f64; 2]>,
    /// Net / signal name (e.g. `"VCC"`, `"GND"`, `"SIG1"`).
    pub label: String,
}

impl Wire {
    /// Construct a wire and reject degenerate (< 2 vertex) polylines.
    pub fn new(polyline: Vec<[f64; 2]>, label: impl Into<String>) -> Result<Self, SymbolError> {
        if polyline.len() < 2 {
            return Err(SymbolError::DegenerateWire(format!(
                "wire needs >= 2 vertices, got {}",
                polyline.len()
            )));
        }
        Ok(Self {
            polyline,
            label: label.into(),
        })
    }
}

/// A schematic = ordered collections of placed symbols + wires.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Schematic {
    /// All symbol placements.
    pub symbols: Vec<PlacedSymbol>,
    /// All wires.
    pub wires: Vec<Wire>,
}

impl Schematic {
    /// Empty schematic.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a placed symbol, returning the new index.
    pub fn push_symbol(&mut self, s: PlacedSymbol) -> usize {
        self.symbols.push(s);
        self.symbols.len() - 1
    }

    /// Add a wire, returning the new index.
    pub fn push_wire(&mut self, w: Wire) -> usize {
        self.wires.push(w);
        self.wires.len() - 1
    }

    /// Total entity count (symbols + wires).
    pub fn entity_count(&self) -> usize {
        self.symbols.len() + self.wires.len()
    }
}

/// Render the schematic as a standalone SVG document. The viewBox is
/// auto-sized to cover all symbol centres + wire endpoints with a
/// 40-unit margin. Returns the SVG source text.
pub fn to_svg(s: &Schematic) -> String {
    // Bounding box of all entity points.
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for ps in &s.symbols {
        min_x = min_x.min(ps.position[0]);
        min_y = min_y.min(ps.position[1]);
        max_x = max_x.max(ps.position[0]);
        max_y = max_y.max(ps.position[1]);
    }
    for w in &s.wires {
        for p in &w.polyline {
            min_x = min_x.min(p[0]);
            min_y = min_y.min(p[1]);
            max_x = max_x.max(p[0]);
            max_y = max_y.max(p[1]);
        }
    }
    if !min_x.is_finite() {
        // Empty schematic.
        min_x = -50.0;
        min_y = -50.0;
        max_x = 50.0;
        max_y = 50.0;
    }
    let margin = 40.0;
    let vb_x = min_x - margin;
    let vb_y = min_y - margin;
    let vb_w = (max_x - min_x).max(1.0) + 2.0 * margin;
    let vb_h = (max_y - min_y).max(1.0) + 2.0 * margin;

    let mut out = String::new();
    out.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" \
         viewBox=\"{vb_x:.2} {vb_y:.2} {vb_w:.2} {vb_h:.2}\" \
         width=\"800\" height=\"600\">\n"
    ));
    out.push_str(
        "  <style>\
         .sym{fill:none;stroke:#000;stroke-width:1.2;stroke-linecap:round;stroke-linejoin:round}\
         .wire{fill:none;stroke:#1a5fb4;stroke-width:1.2}\
         .label{font:10px sans-serif;fill:#333}\
         .desig{font:bold 10px sans-serif;fill:#000}\
         </style>\n",
    );

    // Wires first (so symbols overlay).
    for w in &s.wires {
        out.push_str("  <polyline class=\"wire\" points=\"");
        for (i, p) in w.polyline.iter().enumerate() {
            if i > 0 {
                out.push(' ');
            }
            out.push_str(&format!("{:.2},{:.2}", p[0], p[1]));
        }
        out.push_str("\"/>\n");
        if !w.label.is_empty() {
            // Label at the midpoint of the first segment.
            let a = w.polyline[0];
            let b = w.polyline[1];
            let mx = (a[0] + b[0]) * 0.5;
            let my = (a[1] + b[1]) * 0.5 - 4.0;
            out.push_str(&format!(
                "  <text class=\"label\" x=\"{mx:.2}\" y=\"{my:.2}\">{}</text>\n",
                escape_xml(&w.label)
            ));
        }
    }

    // Symbols.
    for ps in &s.symbols {
        out.push_str(&format!(
            "  <g class=\"sym\" transform=\"translate({:.2} {:.2}) rotate({:.2})\">\n    \
             <path d=\"{}\"/>\n",
            ps.position[0],
            ps.position[1],
            ps.rotation_deg,
            ps.kind.to_svg_path()
        ));
        if !ps.designator.is_empty() {
            out.push_str(&format!(
                "    <text class=\"desig\" x=\"0\" y=\"35\" text-anchor=\"middle\">{}</text>\n",
                escape_xml(&ps.designator)
            ));
        }
        out.push_str("  </g>\n");
    }

    out.push_str("</svg>\n");
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
    fn wire_rejects_singleton() {
        assert!(matches!(
            Wire::new(vec![[0.0, 0.0]], "a"),
            Err(SymbolError::DegenerateWire(_))
        ));
    }

    #[test]
    fn empty_schematic_renders_a_valid_root_svg() {
        let s = Schematic::new();
        let svg = to_svg(&s);
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("viewBox="));
        assert!(svg.trim_end().ends_with("</svg>"));
    }

    #[test]
    fn small_schematic_renders_all_glyphs_and_wire() {
        let mut s = Schematic::new();
        s.push_symbol(PlacedSymbol::new(SymbolKind::Resistor, [0.0, 0.0]));
        s.push_symbol(PlacedSymbol::new(SymbolKind::Capacitor, [100.0, 0.0]));
        s.push_wire(Wire::new(vec![[30.0, 0.0], [70.0, 0.0]], "SIG").unwrap());
        let svg = to_svg(&s);
        assert!(svg.contains("class=\"sym\""));
        assert!(svg.contains("class=\"wire\""));
        assert!(svg.contains(">SIG<"));
        assert!(svg.contains("viewBox="));
    }
}
