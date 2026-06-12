//! Minimal PDF 1.4 export.
//!
//! Writes a single-page PDF with all the drawing's visible/hidden
//! edges, dimensions, and title block as line + text operators in a
//! single content stream. No external crate dependency — the PDF
//! format is reasonable to drive by hand for line-art output.
//!
//! Structure:
//!
//! ```text
//! %PDF-1.4
//! 1 0 obj  ← Catalog
//! 2 0 obj  ← Pages
//! 3 0 obj  ← Page (refs Resources + Contents)
//! 4 0 obj  ← Resources (Font dict)
//! 5 0 obj  ← Contents stream (the actual drawing ops)
//! 6 0 obj  ← Font (Helvetica)
//! xref
//! trailer
//! ```
//!
//! Coordinates are converted from mm to PDF points
//! (1 mm = 2.834645669… pt). PDF origin is bottom-left, same as our
//! sheet space — no y-flip.

use std::path::Path;

use crate::error::TechDrawError;
use crate::sheet::SheetTemplate;
use crate::Drawing;

/// Millimeters per PDF point. PDF uses 1/72 inch = 25.4/72 mm
/// ≈ 0.3527777… mm per point. The inverse is what we apply to mm to
/// get points.
const MM_TO_PT: f64 = 72.0 / 25.4;

/// Render `drawing` to a PDF file at `path`.
pub fn write(drawing: &Drawing, path: &Path) -> Result<(), TechDrawError> {
    let bytes = render_bytes(drawing);
    // R30: PDF bytes are already fully materialised in `bytes`; publish
    // them atomically (sidecar → fsync → rename) so a torn write can't
    // leave a corrupt PDF a viewer would reject.
    valenx_core::io_caps::atomic_write_bytes(path, &bytes)?;
    Ok(())
}

/// Render `drawing` to in-memory PDF bytes.
pub fn render_bytes(drawing: &Drawing) -> Vec<u8> {
    let (w_mm, h_mm) = drawing.sheet.dimensions_mm();
    let w_pt = w_mm * MM_TO_PT;
    let h_pt = h_mm * MM_TO_PT;

    let content = build_content_stream(drawing, h_mm);
    let content_bytes = content.as_bytes();

    let mut out = Vec::new();
    let mut offsets: Vec<usize> = Vec::new();
    let pdf = |s: &str, out: &mut Vec<u8>| out.extend_from_slice(s.as_bytes());

    pdf("%PDF-1.4\n%\u{00E2}\u{00E3}\u{00CF}\u{00D3}\n", &mut out);

    // 1: Catalog
    offsets.push(out.len());
    pdf(
        "1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n",
        &mut out,
    );

    // 2: Pages
    offsets.push(out.len());
    pdf(
        "2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n",
        &mut out,
    );

    // 3: Page
    offsets.push(out.len());
    pdf(
        &format!(
            "3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {w_pt:.3} {h_pt:.3}] \
             /Resources 4 0 R /Contents 5 0 R >>\nendobj\n"
        ),
        &mut out,
    );

    // 4: Resources (Font)
    offsets.push(out.len());
    pdf("4 0 obj\n<< /Font << /F1 6 0 R >> >>\nendobj\n", &mut out);

    // 5: Contents stream
    offsets.push(out.len());
    pdf(
        &format!("5 0 obj\n<< /Length {} >>\nstream\n", content_bytes.len()),
        &mut out,
    );
    out.extend_from_slice(content_bytes);
    pdf("\nendstream\nendobj\n", &mut out);

    // 6: Font (Helvetica)
    offsets.push(out.len());
    pdf(
        "6 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>\nendobj\n",
        &mut out,
    );

    // xref
    let xref_pos = out.len();
    pdf(&format!("xref\n0 {n}\n", n = offsets.len() + 1), &mut out);
    pdf("0000000000 65535 f \n", &mut out);
    for off in &offsets {
        pdf(&format!("{off:010} 00000 n \n"), &mut out);
    }
    pdf(
        &format!(
            "trailer\n<< /Size {n} /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF",
            n = offsets.len() + 1
        ),
        &mut out,
    );
    out
}

/// Build the contents stream — concatenation of "move-to + line-to +
/// stroke" operators plus "Tj" text operators.
fn build_content_stream(drawing: &Drawing, _h_mm: f64) -> String {
    let mut s = String::new();
    // 0.4 pt stroke for sheet border.
    s.push_str(&format!(
        "0.4 w\n0 0 {} {} re\nS\n",
        drawing.sheet.dimensions_mm().0 * MM_TO_PT,
        drawing.sheet.dimensions_mm().1 * MM_TO_PT
    ));

    // Views.
    for view in &drawing.views {
        // Place visible edges with translate + scale.
        s.push_str("q\n");
        s.push_str(&format!(
            "{sc} 0 0 {sc} {tx:.3} {ty:.3} cm\n",
            sc = view.scale * MM_TO_PT,
            tx = view.position[0] * MM_TO_PT,
            ty = view.position[1] * MM_TO_PT,
        ));
        s.push_str("0.5 w\n");
        for seg in &view.visible_edges {
            s.push_str(&format!(
                "{:.4} {:.4} m {:.4} {:.4} l S\n",
                seg[0].0, seg[0].1, seg[1].0, seg[1].1
            ));
        }
        // Hidden edges — dashed (PDF: `[2 1] 0 d`, set dash pattern).
        s.push_str("0.3 w\n[2 1] 0 d\n");
        for seg in &view.hidden_edges {
            s.push_str(&format!(
                "{:.4} {:.4} m {:.4} {:.4} l S\n",
                seg[0].0, seg[0].1, seg[1].0, seg[1].1
            ));
        }
        s.push_str("[] 0 d\nQ\n"); // clear dash + restore graphics state
    }

    // Dimensions (in sheet space, point units). Phase 18B: chains
    // expand inline so the renderer treats them the same way.
    s.push_str("0.3 w\n");
    let mut all_dims: Vec<crate::Dimension> = drawing.dimensions.clone();
    for chain in &drawing.dim_chains {
        all_dims.extend(chain.expand());
    }
    for dim in &all_dims {
        let (segs, (lx, ly, label)) = dim.render_segments();
        for seg in &segs {
            s.push_str(&format!(
                "{:.4} {:.4} m {:.4} {:.4} l S\n",
                seg[0].0 * MM_TO_PT,
                seg[0].1 * MM_TO_PT,
                seg[1].0 * MM_TO_PT,
                seg[1].1 * MM_TO_PT,
            ));
        }
        s.push_str(&text_op(lx * MM_TO_PT, ly * MM_TO_PT, 8.0, &label));
    }

    // Phase 18C-F annotations — minimal line-and-text rendering.
    // Balloons: circle approximated as 8-segment polygon, leader to target.
    for b in &drawing.balloons {
        let cx = b.position[0] * MM_TO_PT;
        let cy = b.position[1] * MM_TO_PT;
        let r = b.radius * MM_TO_PT;
        emit_polygon_approx(&mut s, cx, cy, r, 16);
        s.push_str(&format!(
            "{:.4} {:.4} m {:.4} {:.4} l S\n",
            cx,
            cy,
            b.target_point[0] * MM_TO_PT,
            b.target_point[1] * MM_TO_PT,
        ));
        s.push_str(&text_op(cx - r * 0.4, cy - r * 0.4, 8.0, &b.number));
    }
    // Leaders: line + small inline label.
    for l in &drawing.leaders {
        s.push_str(&format!(
            "{:.4} {:.4} m {:.4} {:.4} l S\n",
            l.start[0] * MM_TO_PT,
            l.start[1] * MM_TO_PT,
            l.end[0] * MM_TO_PT,
            l.end[1] * MM_TO_PT,
        ));
        if !l.text.is_empty() {
            s.push_str(&text_op(
                l.start[0] * MM_TO_PT,
                l.start[1] * MM_TO_PT,
                8.0,
                &l.text,
            ));
        }
    }
    // Welds: reference line + leader + size label.
    for w in &drawing.welds {
        let px = w.position[0] * MM_TO_PT;
        let py = w.position[1] * MM_TO_PT;
        s.push_str(&format!(
            "{:.4} {:.4} m {:.4} {:.4} l S\n",
            px,
            py,
            px + 14.0 * MM_TO_PT,
            py
        ));
        s.push_str(&format!(
            "{:.4} {:.4} m {:.4} {:.4} l S\n",
            px,
            py,
            w.arrow_target[0] * MM_TO_PT,
            w.arrow_target[1] * MM_TO_PT,
        ));
        let label = if w.size.is_empty() {
            w.weld_type.label().to_string()
        } else {
            format!("{} {}", w.weld_type.label(), w.size)
        };
        s.push_str(&text_op(px, py + 2.0, 8.0, &label));
    }
    // Surface finishes: V triangle + Ra label.
    for sf in &drawing.surface_finishes {
        let px = sf.position[0] * MM_TO_PT;
        let py = sf.position[1] * MM_TO_PT;
        let w = 2.0 * MM_TO_PT;
        let hgt = 5.0 * MM_TO_PT;
        s.push_str(&format!(
            "{:.4} {:.4} m {:.4} {:.4} l S\n",
            px,
            py,
            px - w,
            py + hgt
        ));
        s.push_str(&format!(
            "{:.4} {:.4} m {:.4} {:.4} l S\n",
            px,
            py,
            px + w,
            py + hgt
        ));
        if sf.roughness_value > 0.0 {
            s.push_str(&text_op(
                px - w,
                py + hgt + 1.0,
                8.0,
                &format!("Ra {:.2}", sf.roughness_value),
            ));
        }
    }
    // GD&T frames — rectangle outline + glyph + tolerance text.
    for g in &drawing.gdt {
        let px = g.position[0] * MM_TO_PT;
        let py = g.position[1] * MM_TO_PT;
        let cell_h = 5.0 * MM_TO_PT;
        let frame_w = (10.0 + 3.5 * g.datums.len() as f64) * MM_TO_PT;
        s.push_str(&format!("{px:.4} {py:.4} {frame_w:.4} {cell_h:.4} re S\n"));
        let mut datum_chunks = String::new();
        for d in &g.datums {
            datum_chunks.push(' ');
            datum_chunks.push_str(&d.letter);
            datum_chunks.push_str(d.modifier.glyph());
        }
        let label = format!(
            "{} {}{}{}",
            g.geometric_characteristic.label(),
            g.tolerance_value,
            g.material_condition.glyph(),
            datum_chunks
        );
        s.push_str(&text_op(px + 1.0, py + 1.5, 8.0, &label));
    }
    for d in &drawing.datums {
        let px = d.position[0] * MM_TO_PT;
        let py = d.position[1] * MM_TO_PT;
        let bw = 5.0 * MM_TO_PT;
        s.push_str(&format!(
            "{:.4} {:.4} {bw:.4} {bw:.4} re S\n",
            px - bw * 0.5,
            py - bw * 0.5
        ));
        s.push_str(&text_op(px - 1.0, py - 1.5, 8.0, &d.letter));
        s.push_str(&format!(
            "{:.4} {:.4} m {:.4} {:.4} l S\n",
            px,
            py,
            d.leader_target[0] * MM_TO_PT,
            d.leader_target[1] * MM_TO_PT,
        ));
    }

    // Phase 19 — BOM tables.
    s.push_str("0.3 w\n");
    for bp in &drawing.bom_placements {
        let (grid, labels) = bp.bom.render_table(bp.origin);
        for seg in &grid {
            s.push_str(&format!(
                "{:.4} {:.4} m {:.4} {:.4} l S\n",
                seg[0].0 * MM_TO_PT,
                seg[0].1 * MM_TO_PT,
                seg[1].0 * MM_TO_PT,
                seg[1].1 * MM_TO_PT,
            ));
        }
        for (x, y, txt) in labels {
            s.push_str(&text_op(x * MM_TO_PT, y * MM_TO_PT, 7.0, &txt));
        }
    }

    // Phase 19 — Revision blocks.
    for blk in &drawing.revision_blocks {
        let (grid, labels) = blk.render();
        for seg in &grid {
            s.push_str(&format!(
                "{:.4} {:.4} m {:.4} {:.4} l S\n",
                seg[0].0 * MM_TO_PT,
                seg[0].1 * MM_TO_PT,
                seg[1].0 * MM_TO_PT,
                seg[1].1 * MM_TO_PT,
            ));
        }
        for (x, y, txt) in labels {
            s.push_str(&text_op(x * MM_TO_PT, y * MM_TO_PT, 7.0, &txt));
        }
    }

    // Phase 19 — Detail views (bubble + magnified content + caption).
    for dv in &drawing.detail_views {
        if let Some(parent) = drawing.views.get(dv.parent_view_idx) {
            // Bubble in parent's local frame: apply translate + scale
            // through a PDF graphics-state push.
            s.push_str("q\n");
            s.push_str(&format!(
                "{sc} 0 0 {sc} {tx:.3} {ty:.3} cm\n",
                sc = parent.scale * MM_TO_PT,
                tx = parent.position[0] * MM_TO_PT,
                ty = parent.position[1] * MM_TO_PT,
            ));
            s.push_str("0.4 w\n");
            for seg in dv.bubble_segments() {
                s.push_str(&format!(
                    "{:.4} {:.4} m {:.4} {:.4} l S\n",
                    seg[0].0, seg[0].1, seg[1].0, seg[1].1,
                ));
            }
            s.push_str("Q\n");
            // Label near the bubble.
            let tick_x = parent.position[0]
                + (dv.center[0] + dv.radius * std::f64::consts::FRAC_1_SQRT_2 + 1.5) * parent.scale;
            let tick_y = parent.position[1]
                + (dv.center[1] + dv.radius * std::f64::consts::FRAC_1_SQRT_2 + 1.5) * parent.scale;
            s.push_str(&text_op(
                tick_x * MM_TO_PT,
                tick_y * MM_TO_PT,
                10.0,
                &dv.label,
            ));
            // Magnified content drawn at the detail's sheet position.
            let magnified = dv.clip_and_magnify(&parent.visible_edges);
            s.push_str("q\n");
            s.push_str(&format!(
                "{sc} 0 0 {sc} {tx:.3} {ty:.3} cm\n",
                sc = MM_TO_PT,
                tx = dv.position[0] * MM_TO_PT,
                ty = dv.position[1] * MM_TO_PT,
            ));
            s.push_str("0.5 w\n");
            for seg in &magnified {
                s.push_str(&format!(
                    "{:.4} {:.4} m {:.4} {:.4} l S\n",
                    seg[0].0, seg[0].1, seg[1].0, seg[1].1,
                ));
            }
            s.push_str("Q\n");
            // Caption.
            s.push_str(&text_op(
                dv.position[0] * MM_TO_PT,
                (dv.position[1] - 8.0) * MM_TO_PT,
                10.0,
                &dv.detail_caption(),
            ));
        }
    }

    // Title block.
    let tpl = SheetTemplate::for_sheet(&drawing.sheet);
    s.push_str("0.4 w\n");
    for seg in &tpl.title_block_edges() {
        s.push_str(&format!(
            "{:.4} {:.4} m {:.4} {:.4} l S\n",
            seg[0].0 * MM_TO_PT,
            seg[0].1 * MM_TO_PT,
            seg[1].0 * MM_TO_PT,
            seg[1].1 * MM_TO_PT,
        ));
    }
    for (x, y, txt) in tpl.title_block_text_positions(&drawing.sheet) {
        s.push_str(&text_op(x * MM_TO_PT, y * MM_TO_PT, 10.0, &txt));
    }

    s
}

/// Emit a regular-polygon approximation of a circle as a series of
/// "move + line" PDF operators followed by stroke. Used by the
/// balloon renderer where the SVG variant draws an actual circle.
fn emit_polygon_approx(s: &mut String, cx: f64, cy: f64, r: f64, n_sides: usize) {
    if n_sides < 3 {
        return;
    }
    let two_pi = std::f64::consts::TAU;
    let mut first = (cx + r, cy);
    s.push_str(&format!("{:.4} {:.4} m\n", first.0, first.1));
    for i in 1..n_sides {
        let a = (i as f64) * two_pi / n_sides as f64;
        let x = cx + r * a.cos();
        let y = cy + r * a.sin();
        s.push_str(&format!("{x:.4} {y:.4} l\n"));
        if i == n_sides - 1 {
            first = (cx + r, cy);
        }
    }
    s.push_str(&format!("{:.4} {:.4} l S\n", first.0, first.1));
}

fn text_op(x_pt: f64, y_pt: f64, size_pt: f64, text: &str) -> String {
    // BT … ET wraps a text object. Tf selects font + size. Td sets
    // baseline. Tj shows the text literal (PDF strings: `(text)` —
    // we escape `\`, `(`, `)`).
    let esc: String = text
        .chars()
        .flat_map(|c| match c {
            '\\' => "\\\\".chars().collect::<Vec<_>>(),
            '(' => "\\(".chars().collect(),
            ')' => "\\)".chars().collect(),
            // PDF Type1 fonts handle ASCII reliably; map non-ASCII to '?'.
            c if c.is_ascii() => vec![c],
            _ => vec!['?'],
        })
        .collect();
    format!("BT /F1 {size_pt:.2} Tf {x_pt:.3} {y_pt:.3} Td ({esc}) Tj ET\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sheet::Sheet;
    use crate::view::{View, ViewKind};

    /// Task 29 — PDF round-trip test: writing a drawing produces a
    /// non-empty file that starts with the PDF magic bytes.
    #[test]
    fn write_pdf_round_trip() {
        let mut d = Drawing::new(Sheet::a4_landscape("Test", "A. Engineer", "A"));
        let mut v = View::new(ViewKind::Front, 1.0, [50.0, 50.0]);
        v.visible_edges = vec![[(0.0, 0.0), (10.0, 0.0)]];
        d.add_view(v);
        let tmp = std::env::temp_dir().join("valenx_techdraw_pdf.pdf");
        write(&d, &tmp).unwrap();
        let bytes = std::fs::read(&tmp).unwrap();
        assert!(!bytes.is_empty());
        assert!(
            bytes.starts_with(b"%PDF-1.4"),
            "should start with PDF magic"
        );
        // Should contain the EOF marker.
        let tail = &bytes[bytes.len().saturating_sub(8)..];
        assert!(
            String::from_utf8_lossy(tail).contains("%%EOF"),
            "should end with %%EOF, got: {:?}",
            String::from_utf8_lossy(tail),
        );
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn empty_drawing_pdf_is_valid() {
        let d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let bytes = render_bytes(&d);
        assert!(bytes.starts_with(b"%PDF-1.4"));
        let s = String::from_utf8_lossy(&bytes);
        assert!(s.contains("xref"));
        assert!(s.contains("trailer"));
        assert!(s.contains("startxref"));
    }

    #[test]
    fn pdf_includes_title_text_in_content_stream() {
        let d = Drawing::new(Sheet::a4_landscape("UniqueTitle", "Y", "Z"));
        let bytes = render_bytes(&d);
        let s = String::from_utf8_lossy(&bytes);
        assert!(s.contains("UniqueTitle"));
    }

    /// Phase 19 — PDF content stream includes BOM headers, revision
    /// headers, and the detail caption.
    #[test]
    fn pdf_includes_phase19_table_text() {
        use crate::bom::{Bom, BomItem};
        use crate::detail_view::DetailView;
        use crate::revision_block::{RevisionBlock, RevisionEntry};

        let mut d = Drawing::new(Sheet::a3_landscape("X", "Y", "Z"));
        let mut v = View::new(ViewKind::Front, 1.0, [50.0, 100.0]);
        v.visible_edges = vec![[(0.0, 0.0), (20.0, 0.0)]];
        d.add_view(v);
        let mut bom = Bom::new();
        bom.add(BomItem::full("UniqBracket", 2, "P-001", "desc", "Al"));
        bom.renumber_items();
        d.add_bom_placement(bom, [200.0, 100.0]);
        let mut blk = RevisionBlock::new([10.0, 70.0]);
        blk.add_entry(RevisionEntry::new("A", "2026-05-23", "init", "GH", ""));
        d.add_revision_block(blk);
        d.add_detail_view(DetailView::new(
            0,
            [10.0, 5.0],
            4.0,
            [250.0, 150.0],
            2.0,
            "Q",
        ));
        let bytes = render_bytes(&d);
        let s = String::from_utf8_lossy(&bytes);
        assert!(s.contains("Part No.") || s.contains("UniqBracket"));
        assert!(s.contains("Rev"));
        assert!(s.contains("Detail Q"));
    }
}
