//! Revision history block.
//!
//! Engineering drawings carry a revision history table — one row per
//! drawing revision with columns for the rev letter / date /
//! description / who made the change / who approved it. Convention
//! places the table near the title block (typically above it for a
//! standard A-series sheet).
//!
//! [`RevisionEntry`] is one row; [`RevisionBlock`] is the table. The
//! table renders to line + label segments via [`RevisionBlock::render`]
//! suitable for inclusion in any export pipeline (SVG / DXF / PDF).
//!
//! Auto-numbering: new entries appended via [`RevisionBlock::add_entry`]
//! get the next alphabetic letter (`A`, `B`, …) for their `rev` field
//! when constructed with [`RevisionEntry::next`]. Callers that supply
//! their own rev codes can use the basic [`RevisionEntry::new`].

use serde::{Deserialize, Serialize};

/// Output of [`RevisionBlock::render`] (and the symmetric BOM
/// renderer): a pair of `(line_segments, labels)`, where each label is
/// `(x_mm, y_mm, text)`.
pub type RenderedTable = (Vec<[(f64, f64); 2]>, Vec<(f64, f64, String)>);

/// One row in the revision history.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RevisionEntry {
    /// Revision letter / code (e.g. `"A"`, `"B"`, `"01"`).
    pub rev: String,
    /// ISO-8601 date string when the revision was made.
    pub date: String,
    /// Human-readable description of the change.
    pub description: String,
    /// Who made the change (initials or name).
    pub by: String,
    /// Who approved the change (initials or name). Blank until
    /// approval.
    pub approved: String,
}

impl RevisionEntry {
    /// Manually-numbered entry.
    pub fn new(rev: &str, date: &str, description: &str, by: &str, approved: &str) -> Self {
        Self {
            rev: rev.into(),
            date: date.into(),
            description: description.into(),
            by: by.into(),
            approved: approved.into(),
        }
    }

    /// Auto-numbered entry: `rev` is set to the next letter after the
    /// last entry already in `existing` (A → B → … → Z → AA → AB → …).
    /// `existing` typically comes from a [`RevisionBlock::entries`]
    /// snapshot.
    pub fn next(
        existing: &[RevisionEntry],
        date: &str,
        description: &str,
        by: &str,
        approved: &str,
    ) -> Self {
        Self::new(&next_letter(existing), date, description, by, approved)
    }
}

fn next_letter(existing: &[RevisionEntry]) -> String {
    let n = existing.len();
    if n == 0 {
        return "A".into();
    }
    // Counting in base-26: 0 → A, 25 → Z, 26 → AA, …
    let mut k = n; // n = number of *existing* entries; the next one is index n.
    let mut out: Vec<char> = Vec::new();
    loop {
        let digit = (k % 26) as u8;
        out.push((b'A' + digit) as char);
        k /= 26;
        if k == 0 {
            break;
        }
        k -= 1;
    }
    out.reverse();
    out.into_iter().collect()
}

/// A complete revision-history table.
#[derive(Clone, Debug, PartialEq, Default, Serialize, Deserialize)]
pub struct RevisionBlock {
    /// Rows in newest-on-top order (the rendered table draws row 0 at
    /// the top, with the header above it).
    pub entries: Vec<RevisionEntry>,
    /// Position of the table's **lower-left corner** on the sheet
    /// (mm). The drawing standard places this above the title block
    /// — callers can compute that with
    /// [`RevisionBlock::standard_position`].
    pub position: [f64; 2],
}

impl RevisionBlock {
    /// Empty block at the given position.
    pub fn new(position: [f64; 2]) -> Self {
        Self {
            entries: Vec::new(),
            position,
        }
    }

    /// Append an entry to the table. Returns the new entry's row
    /// index (0-based).
    pub fn add_entry(&mut self, entry: RevisionEntry) -> usize {
        let id = self.entries.len();
        self.entries.push(entry);
        id
    }

    /// Standard drawing position: just above the title block, anchored
    /// to the bottom-right corner of the sheet. The title block is
    /// 60 mm tall + a 2 mm gap = 62 mm above the sheet bottom.
    pub fn standard_position(sheet_width_mm: f64) -> [f64; 2] {
        let table_w = Self::standard_total_width();
        [
            sheet_width_mm - table_w,
            crate::sheet::SheetTemplate::TITLE_BLOCK_HEIGHT_MM + 2.0,
        ]
    }

    /// Standard column widths in millimeters: Rev / Date / Description
    /// / By / Approved.
    pub fn standard_column_widths() -> [f64; 5] {
        [12.0, 24.0, 80.0, 16.0, 16.0]
    }

    /// Total table width (sum of [`Self::standard_column_widths`]).
    pub fn standard_total_width() -> f64 {
        Self::standard_column_widths().iter().sum()
    }

    /// Render the table as `(grid_segments, labels)` in sheet mm. The
    /// header row sits at the *top* of the table; entries below.
    pub fn render(&self) -> RenderedTable {
        let col_widths = Self::standard_column_widths();
        let row_h = 6.0;
        let header_h = 6.0;
        let total_width: f64 = col_widths.iter().sum();
        let n_rows = self.entries.len();
        let total_height = header_h + row_h * n_rows as f64;
        let [ox, oy] = self.position;
        let mut grid: Vec<[(f64, f64); 2]> = vec![
            // Outer rectangle.
            [(ox, oy), (ox + total_width, oy)],
            [
                (ox + total_width, oy),
                (ox + total_width, oy + total_height),
            ],
            [
                (ox + total_width, oy + total_height),
                (ox, oy + total_height),
            ],
            [(ox, oy + total_height), (ox, oy)],
        ];
        // Column dividers.
        let mut acc = 0.0;
        for w in &col_widths[..col_widths.len() - 1] {
            acc += w;
            grid.push([(ox + acc, oy), (ox + acc, oy + total_height)]);
        }
        // Header divider (always present, separates header from body).
        let header_y = oy + total_height - header_h;
        grid.push([(ox, header_y), (ox + total_width, header_y)]);
        // Between-data-row dividers.
        for i in 1..n_rows {
            let y = header_y - row_h * i as f64;
            grid.push([(ox, y), (ox + total_width, y)]);
        }

        // Labels.
        let mut labels: Vec<(f64, f64, String)> = Vec::new();
        let pad_x = 1.0;
        let baseline_y_offset = row_h * 0.3;
        let headers = ["Rev", "Date", "Description", "By", "Approved"];
        let mut cx = ox;
        for (i, h) in headers.iter().enumerate() {
            labels.push((cx + pad_x, header_y + baseline_y_offset, (*h).into()));
            cx += col_widths[i];
        }
        for (row_idx, e) in self.entries.iter().enumerate() {
            let y = header_y - row_h * (row_idx + 1) as f64 + baseline_y_offset;
            let cells = [
                e.rev.clone(),
                e.date.clone(),
                e.description.clone(),
                e.by.clone(),
                e.approved.clone(),
            ];
            let mut cx = ox;
            for (i, txt) in cells.iter().enumerate() {
                labels.push((cx + pad_x, y, txt.clone()));
                cx += col_widths[i];
            }
        }
        (grid, labels)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Auto-letter cycles A→B→C…
    #[test]
    fn next_letter_basic_sequence() {
        let mut blk = RevisionBlock::new([100.0, 100.0]);
        blk.add_entry(RevisionEntry::next(
            &blk.entries,
            "2026-05-23",
            "initial",
            "GH",
            "",
        ));
        assert_eq!(blk.entries[0].rev, "A");
        blk.add_entry(RevisionEntry::next(
            &blk.entries,
            "2026-05-24",
            "fix dim",
            "GH",
            "BS",
        ));
        assert_eq!(blk.entries[1].rev, "B");
        blk.add_entry(RevisionEntry::next(
            &blk.entries,
            "2026-05-25",
            "tweak",
            "GH",
            "BS",
        ));
        assert_eq!(blk.entries[2].rev, "C");
    }

    /// Auto-letter wraps Z→AA at 26 entries.
    #[test]
    fn next_letter_wraps_to_aa_after_z() {
        let mut existing: Vec<RevisionEntry> = (0..26)
            .map(|i| RevisionEntry::new(&(b'A' + i as u8).to_string(), "2026-01-01", "x", "x", "x"))
            .collect();
        // The 27th entry (index 26) should be "AA".
        let n = RevisionEntry::next(&existing, "2026-05-23", "next", "G", "");
        assert_eq!(n.rev, "AA");
        existing.push(n);
        // The 28th (index 27) should be "AB".
        let m = RevisionEntry::next(&existing, "2026-05-23", "next", "G", "");
        assert_eq!(m.rev, "AB");
    }

    /// Manual `new` does not consult prior entries.
    #[test]
    fn manual_new_keeps_caller_rev() {
        let e = RevisionEntry::new("01", "2026-05-23", "init", "GH", "");
        assert_eq!(e.rev, "01");
    }

    /// add_entry appends and returns the new index.
    #[test]
    fn add_entry_returns_sequential_indices() {
        let mut blk = RevisionBlock::new([100.0, 100.0]);
        let i = blk.add_entry(RevisionEntry::new("A", "2026-05-23", "a", "G", ""));
        let j = blk.add_entry(RevisionEntry::new("B", "2026-05-24", "b", "G", ""));
        assert_eq!(i, 0);
        assert_eq!(j, 1);
        assert_eq!(blk.entries.len(), 2);
    }

    /// standard_position places the block above the title block.
    #[test]
    fn standard_position_above_title_block() {
        let pos = RevisionBlock::standard_position(297.0); // A4 landscape
                                                           // Y should be title-block height + 2 mm gap = 62.
        assert!((pos[1] - 62.0).abs() < 1e-6);
        // X should be sheet_width - table_width.
        let expected_x = 297.0 - RevisionBlock::standard_total_width();
        assert!((pos[0] - expected_x).abs() < 1e-6);
    }

    /// render emits the expected grid + label counts.
    #[test]
    fn render_emits_grid_and_labels() {
        let mut blk = RevisionBlock::new([10.0, 10.0]);
        blk.add_entry(RevisionEntry::new(
            "A",
            "2026-05-23",
            "initial release",
            "GH",
            "BS",
        ));
        blk.add_entry(RevisionEntry::new(
            "B",
            "2026-05-24",
            "fix dim chain",
            "GH",
            "BS",
        ));
        let (grid, labels) = blk.render();
        // 4 outer + 4 column dividers + 1 header + 1 between-row = 10.
        assert_eq!(grid.len(), 10);
        // 5 headers + 5 cells × 2 rows = 15.
        assert_eq!(labels.len(), 15);
        // Find the "Rev" header.
        assert!(labels.iter().any(|(_, _, t)| t == "Rev"));
        // Find the "A" entry rev label.
        assert!(labels.iter().any(|(_, _, t)| t == "A"));
        assert!(labels.iter().any(|(_, _, t)| t == "initial release"));
    }

    /// Empty block renders just the header row + outer rectangle.
    #[test]
    fn empty_block_renders_only_header() {
        let blk = RevisionBlock::new([0.0, 0.0]);
        let (grid, labels) = blk.render();
        // 4 outer + 4 column dividers + 1 header divider = 9.
        assert_eq!(grid.len(), 9);
        assert_eq!(labels.len(), 5); // 5 header labels only
    }
}
