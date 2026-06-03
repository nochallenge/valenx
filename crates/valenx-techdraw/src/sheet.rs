//! Paper-sized drawing sheet + title-block template.
//!
//! A [`Sheet`] is the canvas a [`crate::Drawing`] renders onto: a
//! paper size (A4 through A0, or a custom mm × mm) and the
//! title-block fields (title, author, revision) that get stamped in
//! the bottom-right corner at export time.
//!
//! [`SheetTemplate`] is the renderer-agnostic helper that turns a
//! [`Sheet`] into geometric primitives (rectangle outlines, text
//! positions). The export crates (SVG / PDF / DXF) consume those
//! primitives — keeping all the layout math in one place means the
//! three exporters stay in lockstep.

use serde::{Deserialize, Serialize};

/// One of the ISO 216 A-series sheet sizes (landscape), or a custom
/// (width, height) pair in millimeters.
///
/// `dimensions_mm` returns `(width, height)` — width is always the
/// longer edge for the named A-series variants since they're
/// landscape-oriented.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum SheetSize {
    /// 297 × 210 mm.
    A4,
    /// 420 × 297 mm.
    A3,
    /// 594 × 420 mm.
    A2,
    /// 841 × 594 mm.
    A1,
    /// 1189 × 841 mm.
    A0,
    /// Free-form width × height in millimeters.
    Custom {
        /// Width in mm. Must be > 0.
        width: f64,
        /// Height in mm. Must be > 0.
        height: f64,
    },
}

impl SheetSize {
    /// Return `(width_mm, height_mm)` for this sheet size.
    pub fn dimensions_mm(&self) -> (f64, f64) {
        match *self {
            SheetSize::A4 => (297.0, 210.0),
            SheetSize::A3 => (420.0, 297.0),
            SheetSize::A2 => (594.0, 420.0),
            SheetSize::A1 => (841.0, 594.0),
            SheetSize::A0 => (1189.0, 841.0),
            SheetSize::Custom { width, height } => (width, height),
        }
    }

    /// Short human label for UI display.
    pub fn label(&self) -> String {
        match self {
            SheetSize::A4 => "A4".into(),
            SheetSize::A3 => "A3".into(),
            SheetSize::A2 => "A2".into(),
            SheetSize::A1 => "A1".into(),
            SheetSize::A0 => "A0".into(),
            SheetSize::Custom { width, height } => format!("Custom {width:.0}×{height:.0}"),
        }
    }
}

/// A paper sheet with title-block fields.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Sheet {
    /// Paper size.
    pub size: SheetSize,
    /// Title of the drawing — stamped in the title block.
    pub title: String,
    /// Author / drafter name.
    pub author: String,
    /// Revision label (e.g. `"A"`, `"v1.2"`).
    pub revision: String,
}

impl Sheet {
    /// Construct an A4 landscape sheet with the given title-block
    /// fields. Mirrors the FreeCAD-style "default drawing template".
    pub fn a4_landscape(title: &str, author: &str, revision: &str) -> Self {
        Self::with_size(SheetSize::A4, title, author, revision)
    }

    /// A3 landscape constructor.
    pub fn a3_landscape(title: &str, author: &str, revision: &str) -> Self {
        Self::with_size(SheetSize::A3, title, author, revision)
    }

    /// A2 landscape constructor.
    pub fn a2_landscape(title: &str, author: &str, revision: &str) -> Self {
        Self::with_size(SheetSize::A2, title, author, revision)
    }

    /// A1 landscape constructor.
    pub fn a1_landscape(title: &str, author: &str, revision: &str) -> Self {
        Self::with_size(SheetSize::A1, title, author, revision)
    }

    /// A0 landscape constructor.
    pub fn a0_landscape(title: &str, author: &str, revision: &str) -> Self {
        Self::with_size(SheetSize::A0, title, author, revision)
    }

    /// Generic constructor — pick any `SheetSize` plus the title-block
    /// fields.
    pub fn with_size(size: SheetSize, title: &str, author: &str, revision: &str) -> Self {
        Self {
            size,
            title: title.into(),
            author: author.into(),
            revision: revision.into(),
        }
    }

    /// Convenience: `(width_mm, height_mm)` of this sheet's paper.
    pub fn dimensions_mm(&self) -> (f64, f64) {
        self.size.dimensions_mm()
    }
}

/// Geometric / textual data for the title block (Task 24-25).
///
/// The standard title block is a 180 × 60 mm rectangle anchored to the
/// bottom-right corner of the sheet, divided horizontally for the
/// title / author / revision / date fields.
pub struct SheetTemplate {
    /// Sheet width in mm.
    pub width_mm: f64,
    /// Sheet height in mm.
    pub height_mm: f64,
}

impl SheetTemplate {
    /// Standard title-block width in millimeters.
    pub const TITLE_BLOCK_WIDTH_MM: f64 = 180.0;
    /// Standard title-block height in millimeters.
    pub const TITLE_BLOCK_HEIGHT_MM: f64 = 60.0;

    /// Build a template for the given sheet.
    pub fn for_sheet(sheet: &Sheet) -> Self {
        let (w, h) = sheet.dimensions_mm();
        Self {
            width_mm: w,
            height_mm: h,
        }
    }

    /// All line segments needed to draw the title-block rectangle +
    /// the three internal divider lines (one per row).
    ///
    /// The rectangle is anchored bottom-right. The divider lines split
    /// the box into four equal-height rows used by
    /// [`Self::title_block_text_positions`].
    pub fn title_block_edges(&self) -> Vec<[(f64, f64); 2]> {
        let w = Self::TITLE_BLOCK_WIDTH_MM;
        let h = Self::TITLE_BLOCK_HEIGHT_MM;
        // Bottom-right corner: at (width_mm, 0).
        let x0 = self.width_mm - w;
        let y0 = 0.0;
        let x1 = self.width_mm;
        let y1 = h;
        let mut out: Vec<[(f64, f64); 2]> = vec![
            [(x0, y0), (x1, y0)], // bottom
            [(x1, y0), (x1, y1)], // right
            [(x1, y1), (x0, y1)], // top
            [(x0, y1), (x0, y0)], // left
        ];
        // Three horizontal dividers (four rows).
        for i in 1..4 {
            let y = y0 + h * (i as f64) / 4.0;
            out.push([(x0, y), (x1, y)]);
        }
        out
    }

    /// Suggested positions + content for the four title-block text
    /// fields. Returns `(x_mm, y_mm, text)` tuples — the export
    /// pipeline picks a font size appropriate to the renderer.
    pub fn title_block_text_positions(&self, sheet: &Sheet) -> Vec<(f64, f64, String)> {
        let w = Self::TITLE_BLOCK_WIDTH_MM;
        let h = Self::TITLE_BLOCK_HEIGHT_MM;
        let x0 = self.width_mm - w;
        let row_h = h / 4.0;
        let pad_x = 4.0;
        let pad_y = row_h / 2.0;
        // Rows from top to bottom: title / author / revision / date.
        let date = today_iso();
        vec![
            (
                x0 + pad_x,
                h - row_h * 0.5 + pad_y - row_h,
                format!("Title: {}", sheet.title),
            ),
            (
                x0 + pad_x,
                h - row_h * 1.5 + pad_y - row_h,
                format!("Drawn by: {}", sheet.author),
            ),
            (
                x0 + pad_x,
                h - row_h * 2.5 + pad_y - row_h,
                format!("Rev: {}", sheet.revision),
            ),
            (
                x0 + pad_x,
                h - row_h * 3.5 + pad_y - row_h,
                format!("Date: {date}"),
            ),
        ]
    }
}

/// Cheap fallback "today" formatter — we don't pull `chrono` for one
/// stamp. Returns a placeholder when the system clock is unavailable.
fn today_iso() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = secs / 86_400;
    // Civil-from-days (Hatcher's algorithm). Good enough for a stamp;
    // no leap-second / timezone fanciness needed.
    let z = days as i64 + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = y + if m <= 2 { 1 } else { 0 };
    format!("{year:04}-{m:02}-{d:02}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_series_dimensions_match_iso_216_landscape() {
        assert_eq!(SheetSize::A4.dimensions_mm(), (297.0, 210.0));
        assert_eq!(SheetSize::A3.dimensions_mm(), (420.0, 297.0));
        assert_eq!(SheetSize::A2.dimensions_mm(), (594.0, 420.0));
        assert_eq!(SheetSize::A1.dimensions_mm(), (841.0, 594.0));
        assert_eq!(SheetSize::A0.dimensions_mm(), (1189.0, 841.0));
    }

    #[test]
    fn custom_size_round_trips_dimensions() {
        let s = SheetSize::Custom {
            width: 500.0,
            height: 300.0,
        };
        assert_eq!(s.dimensions_mm(), (500.0, 300.0));
    }

    #[test]
    fn constructors_set_metadata() {
        let s = Sheet::a4_landscape("Bracket", "A. Engineer", "A");
        assert_eq!(s.title, "Bracket");
        assert_eq!(s.author, "A. Engineer");
        assert_eq!(s.revision, "A");
        assert_eq!(s.dimensions_mm(), (297.0, 210.0));
    }

    #[test]
    fn title_block_edges_returns_outer_rect_and_three_dividers() {
        let s = Sheet::a3_landscape("X", "Y", "Z");
        let t = SheetTemplate::for_sheet(&s);
        let edges = t.title_block_edges();
        // 4 outer + 3 internal dividers = 7 segments.
        assert_eq!(edges.len(), 7);
    }

    #[test]
    fn title_block_text_positions_emits_four_fields() {
        let s = Sheet::a4_landscape("Foo", "Bar", "v1");
        let t = SheetTemplate::for_sheet(&s);
        let texts = t.title_block_text_positions(&s);
        assert_eq!(texts.len(), 4);
        assert!(texts[0].2.contains("Foo"));
        assert!(texts[1].2.contains("Bar"));
        assert!(texts[2].2.contains("v1"));
        assert!(texts[3].2.starts_with("Date:"));
    }

    #[test]
    fn date_helper_is_iso_format() {
        let d = today_iso();
        // YYYY-MM-DD = 10 chars.
        assert_eq!(d.len(), 10);
        assert_eq!(&d[4..5], "-");
        assert_eq!(&d[7..8], "-");
    }
}
