//! 2D DNA / plasmid viewport — the central-panel rendering target for
//! DNA / RNA workbenches.
//!
//! # Research basis
//!
//! Examined tools: Benchling (web 2024), SnapGene 7, ApE 3.1, Geneious
//! Prime 2024, IGV 2.17, UGENE 48.
//!
//! Views implemented here are the ones present in every tool and most
//! useful for exploratory work:
//!
//! - **Linear annotated map** — horizontal backbone with coloured
//!   directional arrows for features (CDS, promoter, terminator, primer
//!   binding sites) and a bp ruler at the top. Scroll to zoom, drag to
//!   pan. Mirrors SnapGene's "Sequence map" + ApE's linear view.
//!
//! - **Circular plasmid map** — molecule drawn as a ring, with coloured
//!   arc overlays for features and radial labels. Total bp count in the
//!   centre. Mirrors SnapGene's circular view and Benchling's plasmid map.
//!
//! # Data source
//!
//! Driven by an optional [`valenx_bioseq::SeqRecord`]. When `None`, a
//! built-in demo plasmid (pDEMO, 1 500 bp) is shown so the viewport is
//! never blank on first open.
//!
//! # Wiring to the Genetics Workbench — next step after user review
//!
//! The plumbing point is the `seq` parameter of [`show`]. To wire the
//! live sequence: in `update.rs`, when building the 2D viewport call,
//! pass `self.genetics.sequence.parsed_record.as_deref()` (or however
//! the Sequence panel exposes its active SeqRecord).

pub mod layout;

use eframe::egui::{self, Align2, Color32, FontId, Painter, Pos2, Rect, Sense, Stroke, Vec2};
use valenx_bioseq::{
    record::{SeqFeature, SeqRecord},
    seq::{Seq, Topology},
    Location, SeqKind, Span, Strand,
};

use std::f32::consts::{FRAC_PI_2, TAU};

// ─────────────────────────────────────────────────────────────────────────────
// Public state
// ─────────────────────────────────────────────────────────────────────────────

/// Which sub-view the 2D DNA viewport renders.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DnaSubView {
    /// Horizontal annotated sequence track with pan + zoom.
    Linear,
    /// Circular plasmid ring diagram.
    Circular,
    /// Linear on top, circular below (the default).
    #[default]
    Both,
}

/// Persistent per-session state for the 2D DNA viewport.
///
/// One instance lives on [`crate::ValenxApp`] alongside the other
/// workbench states — it survives viewport-kind switches so pan/zoom are
/// remembered when the user comes back.
#[derive(Debug)]
pub struct Viewport2dState {
    /// Base-pair index at the left edge of the linear window (fractional
    /// for smooth scrolling). Initialised to 0.
    pub pan: f32,
    /// Bases per pixel. Lower = more zoomed in. Auto-set to fit the
    /// full sequence on the first frame; the user can scroll to override.
    pub bases_per_pixel: f32,
    /// Which sub-view to render.
    pub sub_view: DnaSubView,
    /// `true` → fit the whole sequence into the track on the next frame.
    /// Reset on first render; also set when the user clicks "⊞ Fit".
    pub needs_fit: bool,
    /// Draw individual base letters when very zoomed in (< 0.1 bp/px).
    pub show_bases: bool,
    /// Overlay vertical tick marks at restriction-enzyme cut sites.
    /// (Not yet connected to a live restriction map; reserved for the
    /// follow-up wiring pass.)
    pub show_restriction_sites: bool,
}

impl Default for Viewport2dState {
    fn default() -> Self {
        Self {
            pan: 0.0,
            bases_per_pixel: 1.0,
            sub_view: DnaSubView::Both,
            needs_fit: true,
            show_bases: true,
            show_restriction_sites: false,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Entry point
// ─────────────────────────────────────────────────────────────────────────────

/// Render the 2D DNA viewport into `ui`'s available space.
///
/// `seq` — the record to display, or `None` to use the built-in demo
/// plasmid ([`demo_record`]).
pub fn show(ui: &mut egui::Ui, state: &mut Viewport2dState, seq: Option<&SeqRecord>) {
    let demo = demo_record();
    let record = seq.unwrap_or(&demo);
    let seq_len = record.seq.len();

    let track_w = ui.available_width().max(1.0);

    if state.needs_fit || state.bases_per_pixel <= 0.0 {
        state.bases_per_pixel = layout::fit_zoom(seq_len, track_w);
        state.pan = 0.0;
        state.needs_fit = false;
    }

    draw_toolbar(ui, state, record);
    ui.separator();

    let available = ui.available_rect_before_wrap();
    let (resp, painter) = ui.allocate_painter(
        Vec2::new(available.width(), available.height()),
        Sense::click_and_drag(),
    );
    let canvas = resp.rect;

    painter.rect_filled(canvas, 2.0, Color32::from_gray(26));

    match state.sub_view {
        DnaSubView::Linear => {
            draw_linear_track(&painter, canvas, state, record);
            handle_linear_input(&resp, state, seq_len, canvas.width());
        }
        DnaSubView::Circular => {
            draw_circular_map(&painter, canvas, record);
        }
        DnaSubView::Both => {
            let split_y = canvas.top() + canvas.height() * 0.45;
            let linear_rect =
                Rect::from_min_max(canvas.min, Pos2::new(canvas.max.x, split_y - 1.0));
            let circ_rect =
                Rect::from_min_max(Pos2::new(canvas.min.x, split_y + 3.0), canvas.max);

            draw_linear_track(&painter, linear_rect, state, record);
            handle_linear_input(&resp, state, seq_len, canvas.width());
            draw_circular_map(&painter, circ_rect, record);

            // Divider between the two sub-views
            painter.line_segment(
                [
                    Pos2::new(canvas.left(), split_y + 1.0),
                    Pos2::new(canvas.right(), split_y + 1.0),
                ],
                Stroke::new(1.0, Color32::from_gray(55)),
            );
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Demo record (pDEMO — visible when no sequence is loaded)
// ─────────────────────────────────────────────────────────────────────────────

/// Build a compact demo plasmid (pDEMO, 1 500 bp) with representative
/// features so the viewport is not blank on first open.
pub fn demo_record() -> SeqRecord {
    let seq = Seq::with_topology(
        SeqKind::Dna,
        b"ATGCATGCATGC".repeat(125), // 1 500 bp, valid IUPAC DNA
        Topology::Circular,
    )
    .expect("demo sequence is valid");

    let mut rec = SeqRecord::new("pDEMO", seq);
    rec.description = "Demo plasmid shown when no sequence is loaded".to_string();

    let add = |rec: &mut SeqRecord,
               ft: &str,
               start: usize,
               end: usize,
               strand: Strand,
               label: &str| {
        rec.features.push(
            SeqFeature::new(ft, Location::Single(Span::with_strand(start, end, strand)))
                .with_qualifier("label", label),
        );
    };

    add(&mut rec, "promoter", 50, 100, Strand::Forward, "Ptac");
    add(&mut rec, "primer_bind", 88, 108, Strand::Forward, "M13 fwd");
    add(&mut rec, "CDS", 100, 700, Strand::Forward, "lacZ-α");
    add(&mut rec, "terminator", 700, 750, Strand::Forward, "T1");
    add(&mut rec, "rep_origin", 800, 1050, Strand::Forward, "pUC ori");
    add(&mut rec, "CDS", 1100, 1450, Strand::Reverse, "AmpR");

    rec
}

// ─────────────────────────────────────────────────────────────────────────────
// Toolbar
// ─────────────────────────────────────────────────────────────────────────────

fn draw_toolbar(ui: &mut egui::Ui, state: &mut Viewport2dState, record: &SeqRecord) {
    ui.horizontal_wrapped(|ui| {
        let topo = if record.seq.is_circular() {
            "circular"
        } else {
            "linear"
        };
        ui.label(
            egui::RichText::new(format!(
                "{} · {} bp · {}",
                record.id,
                record.seq.len(),
                topo
            ))
            .strong(),
        );

        ui.separator();
        ui.label("View:");
        ui.selectable_value(&mut state.sub_view, DnaSubView::Linear, "Linear");
        ui.selectable_value(&mut state.sub_view, DnaSubView::Circular, "Circular");
        ui.selectable_value(&mut state.sub_view, DnaSubView::Both, "Both");

        ui.separator();
        ui.checkbox(&mut state.show_bases, "Bases")
            .on_hover_text("Show individual base letters when zoomed in past 10 px/base.");
        ui.checkbox(&mut state.show_restriction_sites, "Cut sites")
            .on_hover_text("Overlay restriction-enzyme cut-site tick marks (stub — wiring to valenx-bioseq digest pending).");

        if ui
            .small_button("⊞ Fit")
            .on_hover_text("Fit the entire sequence into the linear track.")
            .clicked()
        {
            state.needs_fit = true;
        }
    });
}

// ─────────────────────────────────────────────────────────────────────────────
// Linear annotated track
// ─────────────────────────────────────────────────────────────────────────────

const RULER_H: f32 = 22.0;
const BACKBONE_Y_FRAC: f32 = 0.50;
const FEATURE_H: f32 = 13.0;
const FEATURE_GAP: f32 = 4.0;
const ARROW_TIP_W: f32 = 8.0;

fn draw_linear_track(
    painter: &Painter,
    rect: Rect,
    state: &Viewport2dState,
    record: &SeqRecord,
) {
    let seq_len = record.seq.len();

    let bp_x = |bp: f32| -> f32 { layout::bp_to_x(bp, state.pan, state.bases_per_pixel, rect.left()) };

    // --- Ruler background ------------------------------------------------
    let ruler_rect = Rect::from_min_max(rect.min, Pos2::new(rect.max.x, rect.min.y + RULER_H));
    painter.rect_filled(ruler_rect, 0.0, Color32::from_gray(38));

    // --- Ruler ticks and labels ------------------------------------------
    let tick_bp = layout::nice_tick_spacing_bp(state.bases_per_pixel, 80.0) as i64;
    let vis_start = state.pan.floor() as i64 - tick_bp;
    let vis_end = (state.pan + rect.width() * state.bases_per_pixel).ceil() as i64 + tick_bp;
    let first_tick = (vis_start / tick_bp) * tick_bp;
    let mut tick = first_tick;
    while tick <= vis_end {
        let bp_usize = tick as usize;
        if tick >= 0 && bp_usize <= seq_len {
            let x = bp_x(tick as f32);
            if x >= rect.left() && x <= rect.right() {
                painter.line_segment(
                    [
                        Pos2::new(x, ruler_rect.max.y - 8.0),
                        Pos2::new(x, ruler_rect.max.y),
                    ],
                    Stroke::new(1.0, Color32::from_gray(130)),
                );
                painter.text(
                    Pos2::new(x + 2.0, ruler_rect.min.y + 3.0),
                    Align2::LEFT_TOP,
                    format_bp(bp_usize),
                    FontId::monospace(9.0),
                    Color32::from_gray(165),
                );
            }
        }
        tick += tick_bp;
    }

    // --- Backbone --------------------------------------------------------
    let work_h = rect.height() - RULER_H;
    let backbone_y = rect.min.y + RULER_H + work_h * BACKBONE_Y_FRAC;
    let x0 = bp_x(0.0).clamp(rect.left(), rect.right());
    let x1 = bp_x(seq_len as f32).clamp(rect.left(), rect.right());

    painter.line_segment(
        [Pos2::new(x0, backbone_y), Pos2::new(x1, backbone_y)],
        Stroke::new(3.0, Color32::from_gray(88)),
    );
    // Sequence boundary ticks
    for &bx in &[x0, x1] {
        painter.line_segment(
            [Pos2::new(bx, backbone_y - 6.0), Pos2::new(bx, backbone_y + 6.0)],
            Stroke::new(1.5, Color32::from_gray(130)),
        );
    }

    // --- Features --------------------------------------------------------
    for feature in &record.features {
        for span in feature.location.spans() {
            let fx0 = bp_x(span.start as f32);
            let fx1 = bp_x(span.end as f32);
            // Skip features entirely off-screen
            if fx1 < rect.left() || fx0 > rect.right() {
                continue;
            }
            let fx0c = fx0.max(rect.left());
            let fx1c = fx1.min(rect.right());
            let fw = (fx1c - fx0c).max(1.0);

            let (r, g, b) = layout::feature_rgb(&feature.feature_type);
            let color = Color32::from_rgba_premultiplied(r, g, b, 210);

            // Forward strand: above backbone; reverse: below
            let y_mid = match span.strand {
                Strand::Reverse => backbone_y + FEATURE_H / 2.0 + FEATURE_GAP,
                _ => backbone_y - FEATURE_H / 2.0 - FEATURE_GAP,
            };
            let feat_top = y_mid - FEATURE_H / 2.0;
            let feat_bot = y_mid + FEATURE_H / 2.0;

            // Skip features whose row is outside the track rectangle
            if feat_bot < rect.top() + RULER_H || feat_top > rect.bottom() {
                continue;
            }

            draw_feature_arrow(
                painter,
                Pos2::new(fx0c, feat_top),
                Pos2::new(fx1c, feat_bot),
                span.strand,
                color,
            );

            // Label when the feature is wide enough to fit text
            if fw > 30.0 {
                painter.text(
                    Pos2::new((fx0c + fx1c) * 0.5, feat_top - 2.0),
                    Align2::CENTER_BOTTOM,
                    feature.label(),
                    FontId::proportional(9.5),
                    Color32::from_rgb(r, g, b),
                );
            }
        }
    }

    // --- Individual base letters when very zoomed in ---------------------
    if state.show_bases && state.bases_per_pixel < 0.1 {
        let first = vis_start.max(0) as usize;
        let last = (vis_end as usize).min(seq_len);
        for i in first..last {
            if let Some(base) = record.seq.get(i) {
                let bx = bp_x(i as f32 + 0.5);
                if bx >= rect.left() && bx <= rect.right() {
                    painter.text(
                        Pos2::new(bx, backbone_y + 1.0),
                        Align2::CENTER_CENTER,
                        std::str::from_utf8(&[base]).unwrap_or("?"),
                        FontId::monospace(9.0),
                        base_letter_color(base),
                    );
                }
            }
        }
    }
}

/// Draw a directional arrow block for one feature span.
///
/// The bounding box is `[top_left, bot_right]`. The arrowhead points in
/// the direction of `strand`.
fn draw_feature_arrow(
    painter: &Painter,
    top_left: Pos2,
    bot_right: Pos2,
    strand: Strand,
    color: Color32,
) {
    let w = bot_right.x - top_left.x;
    let mid_y = (top_left.y + bot_right.y) * 0.5;
    let tip_w = ARROW_TIP_W.min(w * 0.35);

    if w < 3.0 {
        // Hairline for very narrow features (e.g. restriction sites)
        painter.line_segment(
            [top_left, Pos2::new(top_left.x, bot_right.y)],
            Stroke::new(1.5, color),
        );
        return;
    }

    match strand {
        Strand::Forward => {
            let body_end_x = bot_right.x - tip_w;
            painter.rect_filled(
                Rect::from_min_max(top_left, Pos2::new(body_end_x, bot_right.y)),
                0.0,
                color,
            );
            // Arrowhead as three lines pointing right
            painter.line_segment(
                [Pos2::new(body_end_x, top_left.y), Pos2::new(bot_right.x, mid_y)],
                Stroke::new(2.0, color),
            );
            painter.line_segment(
                [Pos2::new(bot_right.x, mid_y), Pos2::new(body_end_x, bot_right.y)],
                Stroke::new(2.0, color),
            );
        }
        Strand::Reverse => {
            let body_start_x = top_left.x + tip_w;
            painter.rect_filled(
                Rect::from_min_max(Pos2::new(body_start_x, top_left.y), bot_right),
                0.0,
                color,
            );
            // Arrowhead as three lines pointing left
            painter.line_segment(
                [Pos2::new(body_start_x, top_left.y), Pos2::new(top_left.x, mid_y)],
                Stroke::new(2.0, color),
            );
            painter.line_segment(
                [Pos2::new(top_left.x, mid_y), Pos2::new(body_start_x, bot_right.y)],
                Stroke::new(2.0, color),
            );
        }
        Strand::Unknown => {
            painter.rect_filled(
                Rect::from_min_max(top_left, bot_right),
                2.0,
                color,
            );
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Linear track interaction (pan + zoom)
// ─────────────────────────────────────────────────────────────────────────────

fn handle_linear_input(
    resp: &egui::Response,
    state: &mut Viewport2dState,
    seq_len: usize,
    track_px: f32,
) {
    // Scroll wheel → zoom, anchored to mouse position
    if resp.hovered() {
        let scroll_y = resp.ctx.input(|i| i.smooth_scroll_delta.y);
        if scroll_y.abs() > 0.5 {
            let cursor_x = resp
                .ctx
                .input(|i| i.pointer.hover_pos())
                .map(|p| p.x)
                .unwrap_or(resp.rect.center().x);

            // Remember which bp is under the cursor
            let cursor_bp =
                layout::x_to_bp(cursor_x, state.pan, state.bases_per_pixel, resp.rect.left());

            // Zoom towards / away
            let factor = if scroll_y > 0.0 { 0.85 } else { 1.0 / 0.85 };
            state.bases_per_pixel =
                layout::clamp_zoom(state.bases_per_pixel * factor, seq_len);

            // Restore cursor_bp under cursor_x after zoom
            state.pan = cursor_bp - (cursor_x - resp.rect.left()) * state.bases_per_pixel;
        }
    }

    // Drag → pan
    if resp.dragged() {
        let dx = resp.drag_delta().x;
        state.pan -= dx * state.bases_per_pixel;
    }

    state.pan = layout::clamp_pan(state.pan, seq_len, state.bases_per_pixel, track_px);
}

// ─────────────────────────────────────────────────────────────────────────────
// Circular plasmid map
// ─────────────────────────────────────────────────────────────────────────────

fn draw_circular_map(painter: &Painter, rect: Rect, record: &SeqRecord) {
    let seq_len = record.seq.len();
    if seq_len == 0 {
        return;
    }

    let center = rect.center();
    let radius = (rect.width().min(rect.height()) * 0.34).max(30.0);

    // --- Ring backbone (approximated with line segments) -----------------
    let n_ring = 128usize;
    let ring_pts: Vec<Pos2> = (0..=n_ring)
        .map(|i| {
            let a = TAU * (i as f32 / n_ring as f32) - FRAC_PI_2;
            Pos2::new(center.x + radius * a.cos(), center.y + radius * a.sin())
        })
        .collect();
    let ring_stroke = Stroke::new(3.0, Color32::from_gray(90));
    for i in 0..n_ring {
        painter.line_segment([ring_pts[i], ring_pts[i + 1]], ring_stroke);
    }

    // --- Feature arcs ----------------------------------------------------
    let feat_r = radius + 8.0;
    let feat_thickness = 7.0;

    for feature in &record.features {
        for span in feature.location.spans() {
            if span.is_empty() {
                continue;
            }
            let (r, g, b) = layout::feature_rgb(&feature.feature_type);
            let feat_color = Color32::from_rgb(r, g, b);
            let arc_stroke = Stroke::new(feat_thickness, feat_color);

            let a_start = layout::bp_to_angle(span.start, seq_len);
            let a_end = layout::bp_to_angle(span.end, seq_len);

            // Segment count proportional to arc length, min 2
            let n_arc = ((span.len() as f32 / seq_len as f32 * 64.0) as usize).clamp(2, 64);

            for j in 0..n_arc {
                let t0 = j as f32 / n_arc as f32;
                let t1 = (j + 1) as f32 / n_arc as f32;
                let a0 = a_start + (a_end - a_start) * t0;
                let a1 = a_start + (a_end - a_start) * t1;
                let p0 =
                    Pos2::new(center.x + feat_r * a0.cos(), center.y + feat_r * a0.sin());
                let p1 =
                    Pos2::new(center.x + feat_r * a1.cos(), center.y + feat_r * a1.sin());
                painter.line_segment([p0, p1], arc_stroke);
            }

            // Label at arc midpoint, extending outward
            let a_mid = (a_start + a_end) * 0.5;
            let label_r = feat_r + feat_thickness + 14.0;
            let lx = center.x + label_r * a_mid.cos();
            let ly = center.y + label_r * a_mid.sin();
            let align = if a_mid.cos() >= 0.0 {
                Align2::LEFT_CENTER
            } else {
                Align2::RIGHT_CENTER
            };
            painter.text(
                Pos2::new(lx, ly),
                align,
                feature.label(),
                FontId::proportional(9.5),
                feat_color,
            );
        }
    }

    // --- Tick marks at 0%, 25%, 50%, 75% --------------------------------
    for frac in [0.0_f32, 0.25, 0.50, 0.75] {
        let bp_pos = (frac * seq_len as f32) as usize;
        let angle = layout::bp_to_angle(bp_pos, seq_len);
        let inner = Pos2::new(
            center.x + (radius - 6.0) * angle.cos(),
            center.y + (radius - 6.0) * angle.sin(),
        );
        let outer = Pos2::new(
            center.x + (radius + 6.0) * angle.cos(),
            center.y + (radius + 6.0) * angle.sin(),
        );
        painter.line_segment([inner, outer], Stroke::new(1.5, Color32::from_gray(155)));

        // Tick label inside the ring
        let label_r = radius - 16.0;
        painter.text(
            Pos2::new(
                center.x + label_r * angle.cos(),
                center.y + label_r * angle.sin(),
            ),
            Align2::CENTER_CENTER,
            format_bp(bp_pos),
            FontId::monospace(8.0),
            Color32::from_gray(140),
        );
    }

    // --- Centre label: record name + bp count ----------------------------
    painter.text(
        center,
        Align2::CENTER_CENTER,
        format!("{}\n{} bp", record.id, seq_len),
        FontId::proportional(11.0),
        Color32::from_gray(195),
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Small helpers
// ─────────────────────────────────────────────────────────────────────────────

fn format_bp(bp: usize) -> String {
    if bp >= 1_000_000 {
        format!("{:.1}M", bp as f64 / 1_000_000.0)
    } else if bp >= 1_000 {
        format!("{:.1}k", bp as f64 / 1_000.0)
    } else {
        format!("{bp}")
    }
}

fn base_letter_color(base: u8) -> Color32 {
    match base {
        b'A' => Color32::from_rgb(80, 200, 80),
        b'T' => Color32::from_rgb(200, 80, 80),
        b'G' => Color32::from_rgb(80, 80, 220),
        b'C' => Color32::from_rgb(220, 180, 40),
        _ => Color32::from_gray(140),
    }
}
