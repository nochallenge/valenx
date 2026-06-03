//! Geometric Dimensioning and Tolerancing per ASME Y14.5.
//!
//! A [`GdtSymbol`] renders to a "feature control frame": a rectangle
//! split into compartments, the first holding the geometric
//! characteristic glyph (flatness ⏥, position ⌖, etc.), the second
//! the tolerance value (with optional material-condition modifier
//! ⓂⒻⓁ), and the trailing compartments holding the datum-feature
//! references (A, B, C ...).
//!
//! A [`Datum`] is the matching datum-feature symbol — a square box
//! around a letter, attached to a feature by a perpendicular leader.

use serde::{Deserialize, Serialize};

/// One geometric-tolerance characteristic per ASME Y14.5 / ISO 1101.
///
/// The 14 standard characteristics are split across form, profile,
/// orientation, location, and runout. Each maps to a single Unicode
/// glyph the renderer drops into the leading compartment of the
/// feature control frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GeometricCharacteristic {
    /// Form — straightness ⏤.
    Straightness,
    /// Form — flatness ⏥.
    Flatness,
    /// Form — circularity ○.
    Circularity,
    /// Form — cylindricity ⌭.
    Cylindricity,
    /// Profile — line profile ⌒.
    ProfileLine,
    /// Profile — surface profile ⌓.
    ProfileSurface,
    /// Orientation — perpendicularity ⊥.
    Perpendicularity,
    /// Orientation — angularity ∠.
    Angularity,
    /// Orientation — parallelism ∥.
    Parallelism,
    /// Location — position ⌖.
    Position,
    /// Location — concentricity ◎.
    Concentricity,
    /// Location — symmetry ⌯.
    Symmetry,
    /// Runout — circular runout ↗.
    CircularRunout,
    /// Runout — total runout ⌰.
    TotalRunout,
}

impl GeometricCharacteristic {
    /// Unicode glyph for the characteristic — drops directly into an
    /// SVG `<text>` element.
    pub fn glyph(self) -> &'static str {
        match self {
            GeometricCharacteristic::Straightness => "⏤",
            GeometricCharacteristic::Flatness => "⏥",
            GeometricCharacteristic::Circularity => "○",
            GeometricCharacteristic::Cylindricity => "⌭",
            GeometricCharacteristic::ProfileLine => "⌒",
            GeometricCharacteristic::ProfileSurface => "⌓",
            GeometricCharacteristic::Perpendicularity => "⊥",
            GeometricCharacteristic::Angularity => "∠",
            GeometricCharacteristic::Parallelism => "∥",
            GeometricCharacteristic::Position => "⌖",
            GeometricCharacteristic::Concentricity => "◎",
            GeometricCharacteristic::Symmetry => "⌯",
            GeometricCharacteristic::CircularRunout => "↗",
            GeometricCharacteristic::TotalRunout => "⌰",
        }
    }

    /// UI dropdown label.
    pub fn label(self) -> &'static str {
        match self {
            GeometricCharacteristic::Straightness => "Straightness",
            GeometricCharacteristic::Flatness => "Flatness",
            GeometricCharacteristic::Circularity => "Circularity",
            GeometricCharacteristic::Cylindricity => "Cylindricity",
            GeometricCharacteristic::ProfileLine => "Profile (line)",
            GeometricCharacteristic::ProfileSurface => "Profile (surface)",
            GeometricCharacteristic::Perpendicularity => "Perpendicularity",
            GeometricCharacteristic::Angularity => "Angularity",
            GeometricCharacteristic::Parallelism => "Parallelism",
            GeometricCharacteristic::Position => "Position",
            GeometricCharacteristic::Concentricity => "Concentricity",
            GeometricCharacteristic::Symmetry => "Symmetry",
            GeometricCharacteristic::CircularRunout => "Circular runout",
            GeometricCharacteristic::TotalRunout => "Total runout",
        }
    }
}

/// Material condition modifier appended to the tolerance value.
///
/// `RFS` is the default ("regardless of feature size"); MMC ⓂMMR and
/// LMC ⓁLMC tighten or loosen the tolerance based on the feature's
/// produced size.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MaterialCondition {
    /// Regardless of feature size (default, no modifier).
    Rfs,
    /// Maximum material condition Ⓜ.
    Mmc,
    /// Least material condition Ⓛ.
    Lmc,
}

impl MaterialCondition {
    /// Unicode modifier glyph (empty for RFS).
    pub fn glyph(self) -> &'static str {
        match self {
            MaterialCondition::Rfs => "",
            MaterialCondition::Mmc => "Ⓜ",
            MaterialCondition::Lmc => "Ⓛ",
        }
    }
    /// UI label.
    pub fn label(self) -> &'static str {
        match self {
            MaterialCondition::Rfs => "RFS",
            MaterialCondition::Mmc => "MMC",
            MaterialCondition::Lmc => "LMC",
        }
    }
}

/// A datum-feature reference inside a feature control frame.
///
/// `letter` is the datum letter (A, B, C, ...). `modifier` is the
/// material-condition glyph appended to the datum letter (rare —
/// usually RFS).
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct DatumRef {
    /// Datum letter (A, B, C, ...). Typically one uppercase letter.
    pub letter: String,
    /// Material-condition modifier applied to the datum.
    pub modifier: MaterialCondition,
}

impl DatumRef {
    /// Plain datum reference at RFS.
    pub fn new(letter: impl Into<String>) -> Self {
        Self {
            letter: letter.into(),
            modifier: MaterialCondition::Rfs,
        }
    }
}

/// One ASME Y14.5 feature control frame.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GdtSymbol {
    /// Top-left corner of the frame in sheet mm.
    pub position: [f64; 2],
    /// Geometric characteristic.
    pub geometric_characteristic: GeometricCharacteristic,
    /// Tolerance value (e.g. "0.1" or "⌀0.1" — caller-provided).
    pub tolerance_value: String,
    /// Material-condition modifier on the tolerance value.
    pub material_condition: MaterialCondition,
    /// Datum references (primary, secondary, tertiary).
    pub datums: Vec<DatumRef>,
}

impl GdtSymbol {
    /// Construct a position-tolerance feature control frame at RFS
    /// with no datums.
    pub fn new(
        position: [f64; 2],
        characteristic: GeometricCharacteristic,
        tolerance_value: impl Into<String>,
    ) -> Self {
        Self {
            position,
            geometric_characteristic: characteristic,
            tolerance_value: tolerance_value.into(),
            material_condition: MaterialCondition::Rfs,
            datums: Vec::new(),
        }
    }
}

/// Render the feature control frame to an SVG fragment.
pub fn render_frame_svg(g: &GdtSymbol, sheet_height: f64) -> String {
    let (px, py) = (g.position[0], sheet_height - g.position[1]);
    // Compartments: characteristic + tol+modifier + each datum.
    let cell_h = 5.0;
    let cell_w_char = 5.0;
    let cell_w_tol = (g.tolerance_value.chars().count().max(2) as f64) * 1.8
        + 2.0
        + if g.material_condition == MaterialCondition::Rfs {
            0.0
        } else {
            2.5
        };
    let cell_w_datum = 3.5;
    let total_w = cell_w_char + cell_w_tol + cell_w_datum * (g.datums.len() as f64);
    let mut out = String::new();
    out.push_str("    <g class=\"gdt-feature-control-frame\">\n");
    // Outer rectangle.
    out.push_str(&format!(
        "      <rect x=\"{px}\" y=\"{py}\" width=\"{total_w}\" height=\"{cell_h}\" fill=\"white\" stroke=\"black\" stroke-width=\"0.4\"/>\n"
    ));
    // Vertical dividers.
    let mut x = px + cell_w_char;
    out.push_str(&format!(
        "      <line x1=\"{x}\" y1=\"{py}\" x2=\"{x}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.4\"/>\n",
        py + cell_h
    ));
    x += cell_w_tol;
    for _ in &g.datums {
        out.push_str(&format!(
            "      <line x1=\"{x}\" y1=\"{py}\" x2=\"{x}\" y2=\"{}\" stroke=\"black\" stroke-width=\"0.4\"/>\n",
            py + cell_h
        ));
        x += cell_w_datum;
    }
    // Characteristic glyph.
    out.push_str(&format!(
        "      <text x=\"{}\" y=\"{}\" font-size=\"3.5\" text-anchor=\"middle\">{}</text>\n",
        px + cell_w_char * 0.5,
        py + cell_h * 0.75,
        g.geometric_characteristic.glyph()
    ));
    // Tolerance + material condition.
    let tol_label = format!("{}{}", g.tolerance_value, g.material_condition.glyph());
    out.push_str(&format!(
        "      <text x=\"{}\" y=\"{}\" font-size=\"3\" text-anchor=\"middle\">{}</text>\n",
        px + cell_w_char + cell_w_tol * 0.5,
        py + cell_h * 0.75,
        escape_xml(&tol_label),
    ));
    // Each datum.
    let mut dx = px + cell_w_char + cell_w_tol;
    for d in &g.datums {
        let label = format!("{}{}", d.letter, d.modifier.glyph());
        out.push_str(&format!(
            "      <text x=\"{}\" y=\"{}\" font-size=\"3\" text-anchor=\"middle\">{}</text>\n",
            dx + cell_w_datum * 0.5,
            py + cell_h * 0.75,
            escape_xml(&label),
        ));
        dx += cell_w_datum;
    }
    out.push_str("    </g>\n");
    out
}

/// A datum-feature symbol — square box around a letter, attached to a
/// feature by a perpendicular leader.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Datum {
    /// Center of the datum box in sheet mm.
    pub position: [f64; 2],
    /// Datum letter (A, B, C, ...).
    pub letter: String,
    /// Point on the part the leader perpendicularly attaches to.
    pub leader_target: [f64; 2],
}

impl Datum {
    /// New datum at `position` labelled `letter`.
    pub fn new(position: [f64; 2], letter: impl Into<String>, leader_target: [f64; 2]) -> Self {
        Self {
            position,
            letter: letter.into(),
            leader_target,
        }
    }
}

/// Render a datum-feature symbol to an SVG fragment.
pub fn render_datum_svg(d: &Datum, sheet_height: f64) -> String {
    let (px, py) = (d.position[0], sheet_height - d.position[1]);
    let (tx, ty) = (d.leader_target[0], sheet_height - d.leader_target[1]);
    let box_w = 5.0;
    let mut out = String::new();
    out.push_str("    <g class=\"gdt-datum\">\n");
    // Square box around letter.
    out.push_str(&format!(
        "      <rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"white\" stroke=\"black\" stroke-width=\"0.4\"/>\n",
        px - box_w * 0.5,
        py - box_w * 0.5,
        box_w,
        box_w
    ));
    // Letter.
    out.push_str(&format!(
        "      <text x=\"{px}\" y=\"{}\" font-size=\"3.5\" text-anchor=\"middle\">{}</text>\n",
        py + 1.0,
        escape_xml(&d.letter)
    ));
    // Leader: line from box edge to target, terminating in a filled
    // datum triangle. The Y14.5 spec requires a filled equilateral
    // triangle at the feature.
    out.push_str(&format!(
        "      <line x1=\"{px}\" y1=\"{py}\" x2=\"{tx}\" y2=\"{ty}\" stroke=\"black\" stroke-width=\"0.3\"/>\n"
    ));
    let dx = tx - px;
    let dy = ty - py;
    let len = (dx * dx + dy * dy).sqrt().max(1e-9);
    let ux = dx / len;
    let uy = dy / len;
    let nx = -uy;
    let ny = ux;
    let alen = 2.5;
    let aw = 1.0;
    let bx = tx - ux * alen;
    let by = ty - uy * alen;
    out.push_str(&format!(
        "      <polygon points=\"{} {} {} {} {} {}\" fill=\"black\"/>\n",
        tx,
        ty,
        bx - nx * aw,
        by - ny * aw,
        bx + nx * aw,
        by + ny * aw,
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
    fn gdt_new_defaults_to_rfs_no_datums() {
        let g = GdtSymbol::new([10.0, 20.0], GeometricCharacteristic::Position, "0.1");
        assert_eq!(g.material_condition, MaterialCondition::Rfs);
        assert!(g.datums.is_empty());
        assert_eq!(g.tolerance_value, "0.1");
    }

    #[test]
    fn render_frame_includes_rect_and_glyphs() {
        let mut g = GdtSymbol::new([10.0, 20.0], GeometricCharacteristic::Position, "0.1");
        g.datums.push(DatumRef::new("A"));
        g.datums.push(DatumRef::new("B"));
        let svg = render_frame_svg(&g, 200.0);
        assert!(svg.contains("class=\"gdt-feature-control-frame\""));
        assert!(svg.contains("<rect"));
        assert!(svg.contains("⌖"));
        assert!(svg.contains(">0.1<"));
        assert!(svg.contains(">A<"));
        assert!(svg.contains(">B<"));
    }

    #[test]
    fn frame_mmc_modifier_appears_in_tolerance_compartment() {
        let mut g = GdtSymbol::new([0.0, 0.0], GeometricCharacteristic::Position, "0.5");
        g.material_condition = MaterialCondition::Mmc;
        let svg = render_frame_svg(&g, 100.0);
        assert!(svg.contains("0.5Ⓜ"));
    }

    #[test]
    fn datum_render_includes_box_letter_and_triangle() {
        let d = Datum::new([10.0, 20.0], "A", [30.0, 40.0]);
        let svg = render_datum_svg(&d, 100.0);
        assert!(svg.contains("class=\"gdt-datum\""));
        assert!(svg.contains("<rect"));
        assert!(svg.contains(">A<"));
        assert!(svg.contains("<polygon"));
        assert!(svg.contains("fill=\"black\""));
    }

    #[test]
    fn characteristic_glyphs_distinct() {
        let glyphs: Vec<_> = [
            GeometricCharacteristic::Straightness,
            GeometricCharacteristic::Flatness,
            GeometricCharacteristic::Circularity,
            GeometricCharacteristic::Cylindricity,
            GeometricCharacteristic::ProfileLine,
            GeometricCharacteristic::ProfileSurface,
            GeometricCharacteristic::Perpendicularity,
            GeometricCharacteristic::Angularity,
            GeometricCharacteristic::Parallelism,
            GeometricCharacteristic::Position,
            GeometricCharacteristic::Concentricity,
            GeometricCharacteristic::Symmetry,
            GeometricCharacteristic::CircularRunout,
            GeometricCharacteristic::TotalRunout,
        ]
        .iter()
        .map(|c| c.glyph())
        .collect();
        let set: std::collections::HashSet<_> = glyphs.iter().collect();
        assert_eq!(set.len(), glyphs.len(), "all 14 glyphs must be unique");
    }

    #[test]
    fn material_condition_labels_distinct() {
        let labels = [
            MaterialCondition::Rfs.label(),
            MaterialCondition::Mmc.label(),
            MaterialCondition::Lmc.label(),
        ];
        let set: std::collections::HashSet<_> = labels.iter().collect();
        assert_eq!(set.len(), labels.len());
    }

    #[test]
    fn datumref_new_defaults_to_rfs() {
        let d = DatumRef::new("C");
        assert_eq!(d.letter, "C");
        assert_eq!(d.modifier, MaterialCondition::Rfs);
    }
}
