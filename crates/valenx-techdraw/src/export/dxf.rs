//! AutoCAD R12 ASCII DXF export.
//!
//! Writes a minimal DXF file containing a HEADER section (with
//! `$ACADVER = AC1009`) and an ENTITIES section populated with LINE
//! and TEXT records. R12 is widely supported (FreeCAD, LibreCAD,
//! AutoCAD, QCAD, BricsCAD) and avoids the complexity of the post-R13
//! object/handle bookkeeping.
//!
//! No external crate dependency — DXF group-code text format is simple
//! enough to drive by hand for line + text output.
//!
//! Coordinates are emitted in millimeters (DXF units are dimensionless;
//! consumers configure interpretation).

use std::path::Path;

use crate::error::TechDrawError;
use crate::sheet::SheetTemplate;
use crate::Drawing;

/// Render `drawing` to a DXF file at `path`.
pub fn write(drawing: &Drawing, path: &Path) -> Result<(), TechDrawError> {
    let dxf = render(drawing);
    // R30: the whole DXF is already materialised in `dxf`, so route it
    // through the crash-safe atomic writer (sidecar → fsync → rename) —
    // no memory penalty, and a torn DXF a CAD tool would choke on is
    // impossible.
    valenx_core::io_caps::atomic_write_str(path, &dxf)?;
    Ok(())
}

/// Render `drawing` to a DXF string.
pub fn render(drawing: &Drawing) -> String {
    let mut out = String::new();
    // ---- HEADER section: minimum is $ACADVER (declares R12). ----
    out.push_str("0\nSECTION\n2\nHEADER\n");
    out.push_str("9\n$ACADVER\n1\nAC1009\n");
    out.push_str("0\nENDSEC\n");

    // ---- TABLES + BLOCKS sections (empty but expected by some readers). ----
    out.push_str("0\nSECTION\n2\nTABLES\n0\nENDSEC\n");
    out.push_str("0\nSECTION\n2\nBLOCKS\n0\nENDSEC\n");

    // ---- ENTITIES section. ----
    out.push_str("0\nSECTION\n2\nENTITIES\n");

    // Views: emit LINE entities translated to position + scaled.
    for view in &drawing.views {
        for seg in &view.visible_edges {
            push_line(
                &mut out,
                view.position[0] + seg[0].0 * view.scale,
                view.position[1] + seg[0].1 * view.scale,
                view.position[0] + seg[1].0 * view.scale,
                view.position[1] + seg[1].1 * view.scale,
                "0", // layer
            );
        }
        for seg in &view.hidden_edges {
            push_line(
                &mut out,
                view.position[0] + seg[0].0 * view.scale,
                view.position[1] + seg[0].1 * view.scale,
                view.position[0] + seg[1].0 * view.scale,
                view.position[1] + seg[1].1 * view.scale,
                "HIDDEN", // separate layer so callers can re-style as dashed
            );
        }
    }

    // Dimensions — use plain LINE + TEXT instead of full DIMENSION
    // entities (which require a STYLE table entry). The geometry is
    // identical and the layer makes them filterable.
    // Phase 18B — expand dim chains inline through the same pipeline.
    let mut all_dims: Vec<crate::Dimension> = drawing.dimensions.clone();
    for chain in &drawing.dim_chains {
        all_dims.extend(chain.expand());
    }
    for dim in &all_dims {
        let (segs, (lx, ly, label)) = dim.render_segments();
        for seg in &segs {
            push_line(&mut out, seg[0].0, seg[0].1, seg[1].0, seg[1].1, "DIM");
        }
        push_text(&mut out, lx, ly, 2.5, &label, "DIM");
    }

    // Phase 18 annotations — emit as LINE + TEXT on dedicated layers.
    // Balloons: 16-gon outline + leader + number text.
    for b in &drawing.balloons {
        push_circle_approx(
            &mut out,
            b.position[0],
            b.position[1],
            b.radius,
            16,
            "BALLOON",
        );
        push_line(
            &mut out,
            b.position[0],
            b.position[1],
            b.target_point[0],
            b.target_point[1],
            "BALLOON",
        );
        push_text(
            &mut out,
            b.position[0],
            b.position[1],
            2.5,
            &b.number,
            "BALLOON",
        );
    }
    // Leaders.
    for l in &drawing.leaders {
        push_line(
            &mut out, l.start[0], l.start[1], l.end[0], l.end[1], "LEADER",
        );
        if !l.text.is_empty() {
            push_text(
                &mut out,
                l.start[0],
                l.start[1] + 1.0,
                2.5,
                &l.text,
                "LEADER",
            );
        }
    }
    // Welds.
    for w in &drawing.welds {
        push_line(
            &mut out,
            w.position[0],
            w.position[1],
            w.position[0] + 14.0,
            w.position[1],
            "WELD",
        );
        push_line(
            &mut out,
            w.position[0],
            w.position[1],
            w.arrow_target[0],
            w.arrow_target[1],
            "WELD",
        );
        let label = format!("{} {}", w.weld_type.label(), w.size);
        push_text(
            &mut out,
            w.position[0],
            w.position[1] + 1.0,
            2.5,
            &label,
            "WELD",
        );
    }
    // Surface finish — V triangle (two lines) + Ra label.
    for sf in &drawing.surface_finishes {
        let px = sf.position[0];
        let py = sf.position[1];
        push_line(&mut out, px, py, px - 2.0, py + 5.0, "SURFACE");
        push_line(&mut out, px, py, px + 2.0, py + 5.0, "SURFACE");
        if sf.roughness_value > 0.0 {
            push_text(
                &mut out,
                px - 2.0,
                py + 6.0,
                2.5,
                &format!("Ra {:.2}", sf.roughness_value),
                "SURFACE",
            );
        }
    }
    // GD&T frames — rectangle (four edges) + label.
    for g in &drawing.gdt {
        let px = g.position[0];
        let py = g.position[1];
        let w = 10.0 + 3.5 * g.datums.len() as f64;
        let h = 5.0;
        push_line(&mut out, px, py, px + w, py, "GDT");
        push_line(&mut out, px + w, py, px + w, py + h, "GDT");
        push_line(&mut out, px + w, py + h, px, py + h, "GDT");
        push_line(&mut out, px, py + h, px, py, "GDT");
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
        push_text(&mut out, px + 1.0, py + 1.0, 2.5, &label, "GDT");
    }
    // Datums.
    for d in &drawing.datums {
        let bw = 5.0;
        let px = d.position[0];
        let py = d.position[1];
        push_line(
            &mut out,
            px - bw * 0.5,
            py - bw * 0.5,
            px + bw * 0.5,
            py - bw * 0.5,
            "DATUM",
        );
        push_line(
            &mut out,
            px + bw * 0.5,
            py - bw * 0.5,
            px + bw * 0.5,
            py + bw * 0.5,
            "DATUM",
        );
        push_line(
            &mut out,
            px + bw * 0.5,
            py + bw * 0.5,
            px - bw * 0.5,
            py + bw * 0.5,
            "DATUM",
        );
        push_line(
            &mut out,
            px - bw * 0.5,
            py + bw * 0.5,
            px - bw * 0.5,
            py - bw * 0.5,
            "DATUM",
        );
        push_text(&mut out, px - 1.0, py - 1.0, 2.5, &d.letter, "DATUM");
        push_line(
            &mut out,
            px,
            py,
            d.leader_target[0],
            d.leader_target[1],
            "DATUM",
        );
    }

    // Phase 19 — BOM tables.
    for bp in &drawing.bom_placements {
        let (grid, labels) = bp.bom.render_table(bp.origin);
        for seg in &grid {
            push_line(&mut out, seg[0].0, seg[0].1, seg[1].0, seg[1].1, "BOM");
        }
        for (x, y, txt) in labels {
            push_text(&mut out, x, y, 2.0, &txt, "BOM");
        }
    }

    // Phase 19 — Revision blocks.
    for blk in &drawing.revision_blocks {
        let (grid, labels) = blk.render();
        for seg in &grid {
            push_line(&mut out, seg[0].0, seg[0].1, seg[1].0, seg[1].1, "REVISION");
        }
        for (x, y, txt) in labels {
            push_text(&mut out, x, y, 2.0, &txt, "REVISION");
        }
    }

    // Phase 19 — Detail views (bubble on parent + magnified content).
    for dv in &drawing.detail_views {
        if let Some(parent) = drawing.views.get(dv.parent_view_idx) {
            // Bubble in sheet coords: apply parent's translate + scale.
            for seg in dv.bubble_segments() {
                let a = (
                    parent.position[0] + seg[0].0 * parent.scale,
                    parent.position[1] + seg[0].1 * parent.scale,
                );
                let b = (
                    parent.position[0] + seg[1].0 * parent.scale,
                    parent.position[1] + seg[1].1 * parent.scale,
                );
                push_line(&mut out, a.0, a.1, b.0, b.1, "DETAIL");
            }
            // Bubble label.
            let tick_x = parent.position[0]
                + (dv.center[0] + dv.radius * std::f64::consts::FRAC_1_SQRT_2 + 1.5) * parent.scale;
            let tick_y = parent.position[1]
                + (dv.center[1] + dv.radius * std::f64::consts::FRAC_1_SQRT_2 + 1.5) * parent.scale;
            push_text(&mut out, tick_x, tick_y, 3.0, &dv.label, "DETAIL");
            // Magnified content.
            let magnified = dv.clip_and_magnify(&parent.visible_edges);
            for seg in &magnified {
                let a = (dv.position[0] + seg[0].0, dv.position[1] + seg[0].1);
                let b = (dv.position[0] + seg[1].0, dv.position[1] + seg[1].1);
                push_line(&mut out, a.0, a.1, b.0, b.1, "DETAIL");
            }
            // Caption.
            push_text(
                &mut out,
                dv.position[0],
                dv.position[1] - 8.0,
                3.0,
                &dv.detail_caption(),
                "DETAIL",
            );
        }
    }

    // Title block.
    let tpl = SheetTemplate::for_sheet(&drawing.sheet);
    for seg in &tpl.title_block_edges() {
        push_line(
            &mut out,
            seg[0].0,
            seg[0].1,
            seg[1].0,
            seg[1].1,
            "TITLEBLOCK",
        );
    }
    for (x, y, txt) in tpl.title_block_text_positions(&drawing.sheet) {
        push_text(&mut out, x, y, 3.0, &txt, "TITLEBLOCK");
    }

    out.push_str("0\nENDSEC\n");
    out.push_str("0\nEOF\n");
    out
}

fn push_line(out: &mut String, x1: f64, y1: f64, x2: f64, y2: f64, layer: &str) {
    out.push_str("0\nLINE\n");
    out.push_str(&format!("8\n{layer}\n"));
    out.push_str(&format!("10\n{x1:.4}\n20\n{y1:.4}\n30\n0.0\n"));
    out.push_str(&format!("11\n{x2:.4}\n21\n{y2:.4}\n31\n0.0\n"));
}

/// Emit a closed n-sided polygon approximating a circle as a series
/// of LINE entities on `layer`. Used by the balloon renderer where
/// the SVG path uses an actual `<circle>`.
fn push_circle_approx(out: &mut String, cx: f64, cy: f64, r: f64, n_sides: usize, layer: &str) {
    if n_sides < 3 {
        return;
    }
    let two_pi = std::f64::consts::TAU;
    let pt = |i: usize| {
        let a = (i as f64) * two_pi / n_sides as f64;
        (cx + r * a.cos(), cy + r * a.sin())
    };
    for i in 0..n_sides {
        let (x1, y1) = pt(i);
        let (x2, y2) = pt((i + 1) % n_sides);
        push_line(out, x1, y1, x2, y2, layer);
    }
}

fn push_text(out: &mut String, x: f64, y: f64, height: f64, text: &str, layer: &str) {
    // DXF TEXT entity. Group 1 = the literal string; non-ASCII is
    // mapped to '?' to keep the line ASCII-clean (DXF accepts UTF-8
    // with `$DWGCODEPAGE` but we skip that table for the v1 writer).
    let ascii: String = text
        .chars()
        .map(|c| if c.is_ascii() { c } else { '?' })
        .collect();
    out.push_str("0\nTEXT\n");
    out.push_str(&format!("8\n{layer}\n"));
    out.push_str(&format!("10\n{x:.4}\n20\n{y:.4}\n30\n0.0\n"));
    out.push_str(&format!("40\n{height:.4}\n"));
    out.push_str(&format!("1\n{ascii}\n"));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dimension::Dimension;
    use crate::sheet::Sheet;
    use crate::view::{View, ViewKind};

    /// Task 31 — DXF round-trip test: writing a drawing produces a
    /// file whose header section contains `$ACADVER = AC1009`.
    #[test]
    fn dxf_round_trip_writes_acadver_header() {
        let mut d = Drawing::new(Sheet::a4_landscape("Test", "A. Engineer", "A"));
        let mut v = View::new(ViewKind::Front, 1.0, [50.0, 50.0]);
        v.visible_edges = vec![[(0.0, 0.0), (10.0, 0.0)]];
        d.add_view(v);
        d.add_dimension(Dimension::Linear {
            from: [0.0, 0.0],
            to: [10.0, 0.0],
            offset: 5.0,
            value: 10.0,
        });
        let tmp = std::env::temp_dir().join("valenx_techdraw_dxf.dxf");
        write(&d, &tmp).unwrap();
        let s = std::fs::read_to_string(&tmp).unwrap();
        // Header markers.
        assert!(s.contains("$ACADVER"));
        assert!(s.contains("AC1009"));
        // Section delimiters.
        assert!(s.contains("SECTION"));
        assert!(s.contains("ENDSEC"));
        assert!(s.ends_with("EOF\n") || s.contains("EOF"));
        // At least one LINE entity (the view edge).
        assert!(s.contains("\nLINE\n"));
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn empty_drawing_dxf_is_valid_structure() {
        let d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let dxf = render(&d);
        assert!(dxf.starts_with("0\nSECTION\n"));
        // Even an empty drawing has title-block LINE entities.
        assert!(dxf.contains("LINE"));
        // Always exactly one EOF.
        assert_eq!(dxf.matches("\nEOF\n").count(), 1);
    }

    #[test]
    fn hidden_edges_go_on_separate_layer() {
        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        let mut v = View::new(ViewKind::Front, 1.0, [10.0, 10.0]);
        v.hidden_edges = vec![[(0.0, 0.0), (1.0, 1.0)]];
        d.add_view(v);
        let s = render(&d);
        assert!(s.contains("HIDDEN"));
    }

    /// Phase 18 — annotations emit on dedicated layers.
    #[test]
    fn dxf_emits_phase18_annotation_layers() {
        use crate::balloon::Balloon;
        use crate::gdt::{Datum, GdtSymbol, GeometricCharacteristic};
        use crate::leader::Leader;
        use crate::surface_finish::SurfaceFinish;
        use crate::weld::WeldSymbol;

        let mut d = Drawing::new(Sheet::a4_landscape("X", "Y", "Z"));
        d.add_balloon(Balloon::new([10.0, 10.0], "1", [20.0, 20.0]));
        d.add_leader(Leader::new([10.0, 10.0], [20.0, 20.0], "L"));
        d.add_weld(WeldSymbol::new_fillet([10.0, 10.0], [20.0, 20.0], "5"));
        d.add_surface_finish(SurfaceFinish::new([30.0, 30.0], 1.6));
        d.add_gdt(GdtSymbol::new(
            [40.0, 40.0],
            GeometricCharacteristic::Position,
            "0.1",
        ));
        d.add_datum(Datum::new([50.0, 50.0], "A", [60.0, 60.0]));
        let s = render(&d);
        for layer in ["BALLOON", "LEADER", "WELD", "SURFACE", "GDT", "DATUM"] {
            assert!(s.contains(layer), "missing layer `{layer}` in DXF");
        }
    }

    /// Phase 19 — DXF emits BOM / REVISION / DETAIL layers.
    #[test]
    fn dxf_emits_phase19_layers() {
        use crate::bom::{Bom, BomItem};
        use crate::detail_view::DetailView;
        use crate::revision_block::{RevisionBlock, RevisionEntry};

        let mut d = Drawing::new(Sheet::a3_landscape("X", "Y", "Z"));
        let mut v = View::new(ViewKind::Front, 1.0, [50.0, 100.0]);
        v.visible_edges = vec![[(0.0, 0.0), (20.0, 0.0)]];
        d.add_view(v);
        let mut bom = Bom::new();
        bom.add(BomItem::full("Bracket", 2, "P-1", "desc", "Al"));
        bom.renumber_items();
        d.add_bom_placement(bom, [200.0, 100.0]);
        let mut blk = RevisionBlock::new([10.0, 70.0]);
        blk.add_entry(RevisionEntry::new("A", "2026-05-23", "init", "GH", ""));
        d.add_revision_block(blk);
        d.add_detail_view(DetailView::new(0, [10.0, 5.0], 4.0, [250.0, 150.0], 2.0, "A"));
        let s = render(&d);
        for layer in ["BOM", "REVISION", "DETAIL"] {
            assert!(s.contains(layer), "missing layer `{layer}` in DXF");
        }
    }
}
