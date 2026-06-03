//! Viewport HUD overlay: the corner orientation gizmo + live world-cursor
//! coordinate readout. Painted on top of the 3D scene with egui's
//! `Painter` (screen-space, no depth). The axis screen-directions and the
//! cursor ray-pick come from `valenx_viz::scene` (pure, unit-tested); this
//! module only paints + hit-tests.

use eframe::egui;
use valenx_viz::{
    gizmo_axis_screen_dirs, intersect_ground_y0, project_point, ray_from_screen,
    snap_ground_point, GizmoFace, OrbitCamera,
};

const GIZMO_RADIUS: f32 = 44.0;
const GIZMO_MARGIN: f32 = 58.0;
const HIT_RADIUS: f32 = 12.0;

/// Per-axis colours: X red, Y green, Z blue.
const AXIS_COLORS: [egui::Color32; 3] = [
    egui::Color32::from_rgb(222, 82, 82),
    egui::Color32::from_rgb(108, 198, 108),
    egui::Color32::from_rgb(96, 142, 236),
];
const AXIS_LABELS: [&str; 3] = ["X", "Y", "Z"];

/// Draw the orientation gizmo + the cursor-coordinate readout. Returns the
/// gizmo face under `mouse` (if any) so the caller can snap the view on a
/// click.
pub fn draw_overlay(
    painter: &egui::Painter,
    rect: egui::Rect,
    camera: &OrbitCamera,
    mouse: Option<egui::Pos2>,
    grid_spacing: f32,
    snap_enabled: bool,
) -> Option<GizmoFace> {
    // --- corner orientation gizmo (top-right) ---
    let center = egui::pos2(rect.right() - GIZMO_MARGIN, rect.top() + GIZMO_MARGIN);
    let dirs = gizmo_axis_screen_dirs(camera);

    // Build the six axis tips (±X, ±Y, ±Z) with their view-space depth,
    // then paint back-to-front so the nearest tip lands on top.
    let mut tips: Vec<(usize, bool, egui::Pos2, f32)> = Vec::with_capacity(6);
    for (i, dir) in dirs.iter().enumerate() {
        let d = egui::vec2(dir.0[0], dir.0[1]);
        let depth = dir.1;
        tips.push((i, true, center + d * GIZMO_RADIUS, depth));
        tips.push((i, false, center - d * GIZMO_RADIUS, -depth));
    }
    tips.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));

    let mut hovered: Option<GizmoFace> = None;
    for &(i, positive, tip, depth) in &tips {
        let base = AXIS_COLORS[i];
        let toward = depth >= 0.0;
        let col = if toward {
            base
        } else {
            egui::Color32::from_rgba_unmultiplied(base.r(), base.g(), base.b(), 120)
        };
        if positive {
            painter.line_segment([center, tip], egui::Stroke::new(2.0, col));
            painter.circle_filled(tip, 8.0, col);
            painter.text(
                tip,
                egui::Align2::CENTER_CENTER,
                AXIS_LABELS[i],
                egui::FontId::proportional(11.0),
                egui::Color32::from_gray(20),
            );
        } else {
            painter.circle_filled(tip, 5.0, col);
        }
        if let Some(m) = mouse {
            if (m - tip).length() <= HIT_RADIUS {
                hovered = Some(face_for(i, positive));
            }
        }
    }

    // --- live world-cursor coordinate readout (bottom-left) ---
    // When snapping is on, the readout shows the snapped grid node and a
    // marker is drawn there, so the user sees exactly where a click would
    // land — the same lattice the GPU grid draws (Fusion-style snap).
    if let Some(m) = mouse {
        if rect.contains(m) {
            let r = ray_from_screen(
                camera,
                rect.width(),
                rect.height(),
                [m.x - rect.left(), m.y - rect.top()],
            );
            if let Some(p) = intersect_ground_y0(&r) {
                let snapping = snap_enabled && grid_spacing > 0.0;
                let shown = if snapping {
                    snap_ground_point(p, grid_spacing)
                } else {
                    p
                };
                if snapping {
                    // Project the snapped world node back to the screen and
                    // mark it. project_point's origin is the rect top-left.
                    if let Some(sp) = project_point(
                        camera,
                        rect.width(),
                        rect.height(),
                        [shown.x, shown.y, shown.z],
                    ) {
                        draw_snap_marker(painter, rect.min + egui::vec2(sp.x, sp.y));
                    }
                }
                let label = if snapping { "snap  " } else { "cursor" };
                painter.text(
                    egui::pos2(rect.left() + 8.0, rect.bottom() - 42.0),
                    egui::Align2::LEFT_BOTTOM,
                    format!(
                        "{label} · x {:.3}  y {:.3}  z {:.3}",
                        shown.x, shown.y, shown.z
                    ),
                    egui::FontId::monospace(11.0),
                    egui::Color32::from_gray(if snapping { 185 } else { 150 }),
                );
            }
        }
    }

    hovered
}

/// A small diamond + centre dot marking the grid node the cursor snaps
/// to. Soft cyan so it reads over both the dark grid and lit geometry
/// without shouting. Built from line segments (no `PathStroke`
/// dependency) to match the gizmo's drawing style.
fn draw_snap_marker(painter: &egui::Painter, pos: egui::Pos2) {
    let c = egui::Color32::from_rgb(120, 200, 255);
    let s = 5.5;
    let top = pos + egui::vec2(0.0, -s);
    let right = pos + egui::vec2(s, 0.0);
    let bot = pos + egui::vec2(0.0, s);
    let left = pos + egui::vec2(-s, 0.0);
    let stroke = egui::Stroke::new(1.4, c);
    painter.line_segment([top, right], stroke);
    painter.line_segment([right, bot], stroke);
    painter.line_segment([bot, left], stroke);
    painter.line_segment([left, top], stroke);
    painter.circle_filled(pos, 1.5, c);
}

/// Map a clicked axis tip (axis index 0/1/2 = X/Y/Z, `positive` sign) to
/// the gizmo face / canonical view it snaps to.
fn face_for(axis: usize, positive: bool) -> GizmoFace {
    match (axis, positive) {
        (0, true) => GizmoFace::Right,   // +X
        (0, false) => GizmoFace::Left,   // -X
        (1, true) => GizmoFace::Top,     // +Y
        (1, false) => GizmoFace::Bottom, // -Y
        (2, true) => GizmoFace::Front,   // +Z
        _ => GizmoFace::Back,            // -Z
    }
}
