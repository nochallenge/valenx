//! SVG 1.1 export.
//!
//! Writes a plain-text SVG file with:
//! - One `<g>` per view (translated to its on-sheet position, scaled).
//! - `<line>` for every visible edge (solid stroke).
//! - `<line stroke-dasharray>` for every hidden edge.
//! - `<line>` for the title-block frame.
//! - `<text>` for each title-block field.
//! - `<g class="dimension">` per dimension with its rendered segments
//!   and a `<text>` label.
//!
//! Coordinates flow as: sheet-mm directly into the SVG viewBox; the
//! `width` / `height` attributes use `mm` units so the file scales
//! correctly when imported into Inkscape / browsers that respect the
//! viewBox.

use std::path::Path;

use crate::error::TechDrawError;
use crate::sheet::SheetTemplate;
use crate::Drawing;

/// Render `drawing` to an SVG file at `path`.
pub fn write(drawing: &Drawing, path: &Path) -> Result<(), TechDrawError> {
    let svg = render(drawing);
    // R30: SVG text is already fully built in `svg`; publish atomically
    // (sidecar → fsync → rename) so importers never see a half-written
    // document.
    valenx_core::io_caps::atomic_write_str(path, &svg)?;
    Ok(())
}

/// Render `drawing` to an SVG string (handy for tests + UI previews).
pub fn render(drawing: &Drawing) -> String {
    let (w, h) = drawing.sheet.dimensions_mm();
    let mut s = String::new();
    s.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    s.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}mm\" height=\"{h}mm\" \
         viewBox=\"0 0 {w} {h}\" font-family=\"sans-serif\">\n"
    ));
    // Sheet border.
    s.push_str(&format!(
        "  <rect x=\"0\" y=\"0\" width=\"{w}\" height=\"{h}\" fill=\"white\" \
         stroke=\"black\" stroke-width=\"0.4\"/>\n"
    ));

    // Views.
    for view in &drawing.views {
        let (vx, vy) = (view.position[0], view.position[1]);
        let sc = view.scale;
        // SVG y grows downward — our drawing-mm uses bottom-left
        // origin like a paper sheet, so we flip y inside each view's
        // group by mapping sheet y → (h - y). Apply the same flip
        // inside the per-view transform so the edges land right-side-up.
        s.push_str(&format!(
            "  <g class=\"view view-{idx}\" transform=\"translate({vx} {flipped_y}) scale({sc} {neg_sc})\">\n",
            idx = view.id,
            flipped_y = h - vy,
            sc = sc,
            neg_sc = -sc,
        ));
        for seg in &view.visible_edges {
            s.push_str(&line(seg, "stroke=\"black\" stroke-width=\"0.5\""));
        }
        for seg in &view.hidden_edges {
            s.push_str(&line(
                seg,
                "stroke=\"black\" stroke-width=\"0.3\" stroke-dasharray=\"2 1\"",
            ));
        }
        s.push_str("  </g>\n");
    }

    // Dimensions (drawn in sheet space, no per-view transform). Also
    // expand any auto-dim chains (Phase 18B) and render them inline so
    // the rest of the pipeline stays untouched.
    s.push_str("  <g class=\"dimensions\">\n");
    let mut all_dims: Vec<crate::Dimension> = drawing.dimensions.clone();
    for chain in &drawing.dim_chains {
        all_dims.extend(chain.expand());
    }
    for dim in &all_dims {
        let (segs, (lx, ly, label)) = dim.render_segments();
        for seg in &segs {
            s.push_str(&line(seg, "stroke=\"black\" stroke-width=\"0.3\""));
        }
        s.push_str(&format!(
            "    <text x=\"{}\" y=\"{}\" font-size=\"3\" fill=\"black\">{}</text>\n",
            lx,
            // Flip y for SVG.
            h - ly,
            escape_xml(&label),
        ));
    }
    s.push_str("  </g>\n");

    // Phase 18C — balloons + leaders.
    if !drawing.balloons.is_empty() {
        s.push_str("  <g class=\"balloons\">\n");
        for b in &drawing.balloons {
            s.push_str(&crate::balloon::render_svg(b, h));
        }
        s.push_str("  </g>\n");
    }
    if !drawing.leaders.is_empty() {
        s.push_str("  <g class=\"leaders\">\n");
        for l in &drawing.leaders {
            s.push_str(&crate::leader::render_svg(l, h));
        }
        s.push_str("  </g>\n");
    }
    // Phase 18D — weld symbols.
    if !drawing.welds.is_empty() {
        s.push_str("  <g class=\"welds\">\n");
        for w in &drawing.welds {
            s.push_str(&crate::weld::render_svg(w, h));
        }
        s.push_str("  </g>\n");
    }
    // Phase 18E — surface-finish callouts.
    if !drawing.surface_finishes.is_empty() {
        s.push_str("  <g class=\"surface-finishes\">\n");
        for sf in &drawing.surface_finishes {
            s.push_str(&crate::surface_finish::render_svg(sf, h));
        }
        s.push_str("  </g>\n");
    }
    // Phase 18F — GD&T feature control frames + datum symbols.
    if !drawing.gdt.is_empty() {
        s.push_str("  <g class=\"gdt-frames\">\n");
        for g in &drawing.gdt {
            s.push_str(&crate::gdt::render_frame_svg(g, h));
        }
        s.push_str("  </g>\n");
    }
    if !drawing.datums.is_empty() {
        s.push_str("  <g class=\"datums\">\n");
        for d in &drawing.datums {
            s.push_str(&crate::gdt::render_datum_svg(d, h));
        }
        s.push_str("  </g>\n");
    }

    // Phase 19 — BOM tables.
    if !drawing.bom_placements.is_empty() {
        s.push_str("  <g class=\"bom-tables\">\n");
        for bp in &drawing.bom_placements {
            let (grid, labels) = bp.bom.render_table(bp.origin);
            for seg in &grid {
                let flipped = [(seg[0].0, h - seg[0].1), (seg[1].0, h - seg[1].1)];
                s.push_str(&line(&flipped, "stroke=\"black\" stroke-width=\"0.3\""));
            }
            for (x, y, txt) in labels {
                s.push_str(&format!(
                    "    <text x=\"{}\" y=\"{}\" font-size=\"2.5\" fill=\"black\">{}</text>\n",
                    x,
                    h - y,
                    escape_xml(&txt),
                ));
            }
        }
        s.push_str("  </g>\n");
    }

    // Phase 19 — Revision blocks.
    if !drawing.revision_blocks.is_empty() {
        s.push_str("  <g class=\"revision-blocks\">\n");
        for blk in &drawing.revision_blocks {
            let (grid, labels) = blk.render();
            for seg in &grid {
                let flipped = [(seg[0].0, h - seg[0].1), (seg[1].0, h - seg[1].1)];
                s.push_str(&line(&flipped, "stroke=\"black\" stroke-width=\"0.3\""));
            }
            for (x, y, txt) in labels {
                s.push_str(&format!(
                    "    <text x=\"{}\" y=\"{}\" font-size=\"2.5\" fill=\"black\">{}</text>\n",
                    x,
                    h - y,
                    escape_xml(&txt),
                ));
            }
        }
        s.push_str("  </g>\n");
    }

    // Phase 19 — Detail view bubbles drawn on the parent view + the
    // magnified detail itself drawn at the detail's `position` on the
    // sheet. The bubble is rendered in the parent's local frame
    // (translate + scale + y-flip applied via a <g> wrapper); the
    // detail's magnified content is rendered the same way at its
    // sheet placement.
    if !drawing.detail_views.is_empty() {
        s.push_str("  <g class=\"detail-views\">\n");
        for dv in &drawing.detail_views {
            // Bubble on the parent.
            if let Some(parent) = drawing.views.get(dv.parent_view_idx) {
                let (vx, vy) = (parent.position[0], parent.position[1]);
                let sc = parent.scale;
                s.push_str(&format!(
                    "    <g class=\"detail-bubble\" transform=\"translate({vx} {flipped_y}) scale({sc} {neg_sc})\">\n",
                    flipped_y = h - vy,
                    neg_sc = -sc,
                ));
                for seg in dv.bubble_segments() {
                    s.push_str(&line(&seg, "stroke=\"black\" stroke-width=\"0.4\""));
                }
                s.push_str("    </g>\n");
                // Bubble label in sheet coords near the leader tick.
                let tick_x = parent.position[0]
                    + (dv.center[0] + dv.radius * std::f64::consts::FRAC_1_SQRT_2 + 1.5) * sc;
                let tick_y = parent.position[1]
                    + (dv.center[1] + dv.radius * std::f64::consts::FRAC_1_SQRT_2 + 1.5) * sc;
                s.push_str(&format!(
                    "    <text x=\"{}\" y=\"{}\" font-size=\"4\" fill=\"black\">{}</text>\n",
                    tick_x,
                    h - tick_y,
                    escape_xml(&dv.label),
                ));
                // Magnified detail view.
                let magnified = dv.clip_and_magnify(&parent.visible_edges);
                let dx = dv.position[0];
                let dy = dv.position[1];
                s.push_str(&format!(
                    "    <g class=\"detail-magnified\" transform=\"translate({dx} {flipped_dy}) scale(1 -1)\">\n",
                    flipped_dy = h - dy,
                ));
                for seg in &magnified {
                    s.push_str(&line(seg, "stroke=\"black\" stroke-width=\"0.5\""));
                }
                s.push_str("    </g>\n");
                // Detail caption text.
                s.push_str(&format!(
                    "    <text x=\"{}\" y=\"{}\" font-size=\"4\" fill=\"black\">{}</text>\n",
                    dv.position[0],
                    h - (dv.position[1] - 8.0),
                    escape_xml(&dv.detail_caption()),
                ));
            }
        }
        s.push_str("  </g>\n");
    }

    // Title block.
    let tpl = SheetTemplate::for_sheet(&drawing.sheet);
    s.push_str("  <g class=\"title-block\">\n");
    for seg in &tpl.title_block_edges() {
        let flipped = [(seg[0].0, h - seg[0].1), (seg[1].0, h - seg[1].1)];
        s.push_str(&line(&flipped, "stroke=\"black\" stroke-width=\"0.4\""));
    }
    for (x, y, txt) in tpl.title_block_text_positions(&drawing.sheet) {
        s.push_str(&format!(
            "    <text x=\"{}\" y=\"{}\" font-size=\"3.5\" fill=\"black\">{}</text>\n",
            x,
            h - y,
            escape_xml(&txt),
        ));
    }
    s.push_str("  </g>\n");

    s.push_str("</svg>\n");
    s
}

fn line(seg: &[(f64, f64); 2], attrs: &str) -> String {
    format!(
        "    <line x1=\"{:.4}\" y1=\"{:.4}\" x2=\"{:.4}\" y2=\"{:.4}\" {attrs}/>\n",
        seg[0].0, seg[0].1, seg[1].0, seg[1].1
    )
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
    use crate::dimension::Dimension;
    use crate::sheet::Sheet;
    use crate::view::{View, ViewKind};

    /// Task 27 — SVG round-trip test: a drawing with one view + one
    /// dimension yields an SVG with the expected number of `<line>`
    /// elements.
    #[test]
    fn drawing_with_view_and_dimension_emits_expected_lines() {
        let mut d = Drawing::new(Sheet::a4_landscape("Test", "A. Engineer", "A"));
        let mut v = View::new(ViewKind::Front, 1.0, [50.0, 50.0]);
        // Manually populate edges so we don't depend on solid extraction.
        v.visible_edges = vec![[(0.0, 0.0), (10.0, 0.0)], [(10.0, 0.0), (10.0, 5.0)]];
        v.hidden_edges = vec![[(0.0, 0.0), (10.0, 5.0)]];
        d.add_view(v);
        d.add_dimension(Dimension::Linear {
            from: [50.0, 50.0],
            to: [60.0, 50.0],
            offset: 10.0,
            value: 10.0,
        });
        let svg = render(&d);
        // 2 visible edges + 1 hidden edge + 7 dimension lines + 4 outer
        // title-block + 3 title-block dividers + 1 sheet border (rect,
        // not line) = at least 16 `<line>` elements (border is `<rect>`).
        let line_count = svg.matches("<line").count();
        assert!(line_count >= 16, "got only {line_count} lines");
        assert!(svg.contains("<svg "));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("10.00")); // dimension value label
    }

    #[test]
    fn empty_drawing_renders_valid_svg() {
        let d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let svg = render(&d);
        assert!(svg.starts_with("<?xml"));
        assert!(svg.contains("<svg "));
        assert!(svg.ends_with("</svg>\n"));
    }

    #[test]
    fn writes_to_file() {
        let d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let tmp = std::env::temp_dir().join("valenx_techdraw_svg.svg");
        write(&d, &tmp).unwrap();
        assert!(tmp.exists());
        let bytes = std::fs::metadata(&tmp).unwrap().len();
        assert!(bytes > 100, "SVG should be non-trivial");
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn escapes_xml_special_chars_in_title() {
        let d = Drawing::new(Sheet::a4_landscape("X<>&Y", "A. Engineer", "A"));
        let svg = render(&d);
        assert!(svg.contains("X&lt;&gt;&amp;Y"));
    }

    /// Phase 18 — every new annotation kind shows up in the SVG.
    #[test]
    fn svg_renders_all_phase18_annotations() {
        use crate::balloon::Balloon;
        use crate::dim_chain::{DimChain, DimChainKind};
        use crate::gdt::{Datum, GdtSymbol, GeometricCharacteristic};
        use crate::leader::Leader;
        use crate::surface_finish::SurfaceFinish;
        use crate::weld::WeldSymbol;

        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        d.add_balloon(Balloon::new([20.0, 20.0], "1", [40.0, 40.0]));
        d.add_leader(Leader::new([10.0, 10.0], [30.0, 30.0], "L1"));
        d.add_weld(WeldSymbol::new_fillet([15.0, 15.0], [25.0, 25.0], "5"));
        d.add_surface_finish(SurfaceFinish::new([35.0, 35.0], 3.2));
        d.add_gdt(GdtSymbol::new(
            [45.0, 45.0],
            GeometricCharacteristic::Flatness,
            "0.1",
        ));
        d.add_datum(Datum::new([55.0, 55.0], "A", [60.0, 60.0]));
        let mut chain = DimChain::new(DimChainKind::Chain, 2.0);
        chain.entries = vec![[0.0, 0.0], [5.0, 0.0], [10.0, 0.0]];
        d.add_dim_chain(chain);
        let svg = render(&d);
        assert!(svg.contains("class=\"balloons\""));
        assert!(svg.contains("class=\"leaders\""));
        assert!(svg.contains("class=\"welds\""));
        assert!(svg.contains("class=\"surface-finishes\""));
        assert!(svg.contains("class=\"gdt-frames\""));
        assert!(svg.contains("class=\"datums\""));
        // Chain produced 2 dim labels.
        let labels = svg.matches("font-size=\"3\" fill=\"black\"").count();
        assert!(labels >= 2);
    }

    /// Phase 19 — SVG renders BOM tables + revision blocks + detail
    /// view bubble + magnified content.
    #[test]
    fn svg_renders_phase19_tables_and_detail_view() {
        use crate::bom::{Bom, BomItem};
        use crate::detail_view::DetailView;
        use crate::revision_block::{RevisionBlock, RevisionEntry};

        let mut d = Drawing::new(Sheet::a3_landscape("X", "Y", "Z"));
        let mut v = View::new(ViewKind::Front, 1.0, [50.0, 100.0]);
        v.visible_edges = vec![[(0.0, 0.0), (20.0, 0.0)], [(20.0, 0.0), (20.0, 10.0)]];
        d.add_view(v);

        let mut bom = Bom::new();
        bom.add(BomItem::full("Bracket", 2, "P-1", "Mounting bracket", "Al"));
        bom.renumber_items();
        d.add_bom_placement(bom, [200.0, 100.0]);

        let mut blk = RevisionBlock::new([10.0, 70.0]);
        blk.add_entry(RevisionEntry::new("A", "2026-05-23", "init", "GH", ""));
        d.add_revision_block(blk);

        d.add_detail_view(DetailView::new(0, [10.0, 5.0], 4.0, [250.0, 150.0], 2.0, "A"));

        let svg = render(&d);
        assert!(svg.contains("class=\"bom-tables\""));
        assert!(svg.contains("class=\"revision-blocks\""));
        assert!(svg.contains("class=\"detail-views\""));
        assert!(svg.contains("class=\"detail-bubble\""));
        assert!(svg.contains("class=\"detail-magnified\""));
        // Header text and a cell value should both appear.
        assert!(svg.contains("Part No."));
        assert!(svg.contains("Mounting bracket"));
        assert!(svg.contains("Rev"));
        // Detail caption uses the standard format.
        assert!(svg.contains("Detail A"));
    }
}
