//! 2D drafting workbench — LibreCAD-style 2D CAD on `valenx-librecad-2d`.
//!
//! A right-side panel that holds a [`Drawing2D`] and renders its entities on a
//! pan/zoom 2D canvas (an egui painter): lines, circles, arcs, polylines, and
//! text. Add primitives from the form, and export the drawing to DXF. The
//! drawing model + DXF round-trip are headless-testable; the canvas is the
//! interactive view.

use eframe::egui;

use valenx_librecad_2d::{dxf, Drawing2D, Entity2D};

use crate::ValenxApp;

/// Persistent state for the 2D drafting workbench.
pub struct Draft2dWorkbenchState {
    drawing: Drawing2D,
    /// Drawing-space point shown at the canvas centre.
    pan: egui::Vec2,
    /// Pixels per drawing unit.
    scale: f32,
    /// Line input: [x1, y1, x2, y2].
    line: [f64; 4],
    /// Circle input: [cx, cy, r].
    circle: [f64; 3],
    status: String,
}

impl Default for Draft2dWorkbenchState {
    fn default() -> Self {
        Self {
            drawing: demo_drawing(),
            pan: egui::Vec2::new(30.0, 20.0),
            scale: 6.0,
            line: [0.0, 0.0, 50.0, 30.0],
            circle: [30.0, 20.0, 10.0],
            status: String::new(),
        }
    }
}

/// A small sample drawing: a closed rectangle, a circle, a diagonal, a label.
fn demo_drawing() -> Drawing2D {
    let mut d = Drawing2D::new();
    d.add(Entity2D::Polyline {
        layer: "0".to_string(),
        vertices: vec![[0.0, 0.0], [60.0, 0.0], [60.0, 40.0], [0.0, 40.0]],
        closed: true,
    });
    d.add(Entity2D::Circle {
        layer: "0".to_string(),
        centre: [30.0, 20.0],
        radius: 12.0,
    });
    d.add(Entity2D::Line {
        layer: "0".to_string(),
        a: [0.0, 0.0],
        b: [60.0, 40.0],
    });
    d.add(Entity2D::Text {
        layer: "0".to_string(),
        position: [1.0, 42.0],
        height: 4.0,
        text: "valenx draft".to_string(),
    });
    d
}

/// Convert a [`Drawing2D`] into the plain-data [`crate::Draft2dView`] the
/// workspace-tile painter consumes — drop the CAD layer / DXF metadata, keep
/// only the geometry, and compute the overall extent so the painter can fit it
/// to the tile. `Text` entities are skipped (the tile draws geometry only).
fn drawing_to_view(drawing: &Drawing2D) -> crate::Draft2dView {
    let mut entities = Vec::new();
    let mut min = [f64::INFINITY; 2];
    let mut max = [f64::NEG_INFINITY; 2];
    let mut grow = |p: &[f64; 2]| {
        min[0] = min[0].min(p[0]);
        min[1] = min[1].min(p[1]);
        max[0] = max[0].max(p[0]);
        max[1] = max[1].max(p[1]);
    };
    for e in &drawing.entities {
        match e {
            Entity2D::Line { a, b, .. } => {
                grow(a);
                grow(b);
                entities.push(crate::Draft2dEntity::Line { a: *a, b: *b });
            }
            Entity2D::Circle { centre, radius, .. } => {
                grow(&[centre[0] - radius, centre[1] - radius]);
                grow(&[centre[0] + radius, centre[1] + radius]);
                entities.push(crate::Draft2dEntity::Circle {
                    centre: *centre,
                    radius: *radius,
                });
            }
            Entity2D::Arc {
                centre,
                radius,
                start_angle_deg,
                end_angle_deg,
                ..
            } => {
                // Bound by the full circle (a safe superset of the arc).
                grow(&[centre[0] - radius, centre[1] - radius]);
                grow(&[centre[0] + radius, centre[1] + radius]);
                entities.push(crate::Draft2dEntity::Arc {
                    centre: *centre,
                    radius: *radius,
                    start_angle_deg: *start_angle_deg,
                    end_angle_deg: *end_angle_deg,
                });
            }
            Entity2D::Polyline {
                vertices, closed, ..
            } => {
                for v in vertices {
                    grow(v);
                }
                entities.push(crate::Draft2dEntity::Polyline {
                    vertices: vertices.clone(),
                    closed: *closed,
                });
            }
            // Text + any future variant: geometry-only tile, skip.
            _ => {}
        }
    }
    // A drawing with no boundable geometry leaves the default unit box.
    if !min[0].is_finite() {
        min = [0.0, 0.0];
        max = [1.0, 1.0];
    }
    crate::Draft2dView {
        entities,
        bounds: (min, max),
        lines: vec![
            format!("{} entities", drawing.entities.len()),
            format!(
                "extent {:.0} × {:.0} units",
                (max[0] - min[0]).abs(),
                (max[1] - min[1]).abs()
            ),
            "valenx-librecad-2d · DXF".into(),
        ],
    }
}

/// Build the agent-bridge **`draft2d` product** — the canonical demo 2-D CAD
/// drawing (a closed rectangle, a circle, a diagonal) as a
/// [`crate::Workspace2dKind::Draft2d`] painted by the tile's 2-D branch
/// (mirroring `rcbeam` / `dna`). A 2-D drawing, NOT a mesh: `mesh: None`,
/// `kind2d: Some(Draft2d(..))`.
pub(crate) fn draft2d_product() -> crate::WorkspaceProduct {
    let drawing = demo_drawing();
    let view = drawing_to_view(&drawing);
    let lines = view.lines.clone();
    crate::WorkspaceProduct {
        title: "2-D Drawing".into(),
        lines,
        mesh: None,
        vertex_colors: None,
        camera: valenx_viz::OrbitCamera::default(),
        kind2d: Some(crate::Workspace2dKind::Draft2d(view)),
        last_export: None,
        image: None,
        image_texture: None,
    }
}

/// Draw the 2D drafting workbench (a no-op unless toggled on via
/// View → 2D Drafting).
pub fn draw_draft2d_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_draft2d_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_draft2d_workbench",
        "2D Drafting",
        |app, ui| {
            ui.label(
                egui::RichText::new("LibreCAD-style · valenx-librecad-2d")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.draft2d;

            ui.horizontal(|ui| {
                if ui.button("Demo").clicked() {
                    s.drawing = demo_drawing();
                }
                if ui.button("Clear").clicked() {
                    s.drawing = Drawing2D::new();
                }
                if ui.button("Export DXF").clicked() {
                    let text = dxf::serialise(&s.drawing);
                    s.status = format!(
                        "DXF: {} entities, {} chars",
                        s.drawing.entities.len(),
                        text.len()
                    );
                }
            });
            ui.horizontal(|ui| {
                ui.label("line");
                for v in &mut s.line {
                    ui.add(egui::DragValue::new(v).speed(1.0));
                }
                if ui.small_button("+").clicked() {
                    s.drawing.add(Entity2D::Line {
                        layer: "0".to_string(),
                        a: [s.line[0], s.line[1]],
                        b: [s.line[2], s.line[3]],
                    });
                }
            });
            ui.horizontal(|ui| {
                ui.label("circle");
                for v in &mut s.circle {
                    ui.add(egui::DragValue::new(v).speed(1.0));
                }
                if ui.small_button("+").clicked() {
                    s.drawing.add(Entity2D::Circle {
                        layer: "0".to_string(),
                        centre: [s.circle[0], s.circle[1]],
                        radius: s.circle[2].max(0.1),
                    });
                }
            });
            ui.horizontal(|ui| {
                ui.label(format!("{} entities", s.drawing.entities.len()));
                ui.separator();
                ui.label("zoom");
                ui.add(
                    egui::DragValue::new(&mut s.scale)
                        .speed(0.2)
                        .range(0.5..=60.0),
                );
            });
            if !s.status.is_empty() {
                ui.label(egui::RichText::new(&s.status).small().weak());
            }
            ui.separator();

            // Canvas — pan with drag, draw the entities.
            let (resp, painter) = ui.allocate_painter(
                egui::vec2(ui.available_width(), ui.available_height().max(200.0)),
                egui::Sense::drag(),
            );
            let rect = resp.rect;
            painter.rect_filled(rect, 2.0, egui::Color32::from_gray(18));
            if resp.dragged() {
                let d = resp.drag_delta();
                s.pan -= egui::vec2(d.x / s.scale, -d.y / s.scale);
            }
            draw_entities(&painter, rect, s.pan, s.scale, &s.drawing);
        },
    );
    if close {
        app.show_draft2d_workbench = false;
    }
}

/// Render every entity in `drawing` onto `painter`, mapping drawing units to
/// screen pixels (CAD y-up → screen y-down).
fn draw_entities(
    painter: &egui::Painter,
    rect: egui::Rect,
    pan: egui::Vec2,
    scale: f32,
    drawing: &Drawing2D,
) {
    let c = rect.center();
    let to_screen = |p: &[f64; 2]| -> egui::Pos2 {
        egui::pos2(
            c.x + (p[0] as f32 - pan.x) * scale,
            c.y - (p[1] as f32 - pan.y) * scale,
        )
    };
    let stroke = egui::Stroke::new(1.5, egui::Color32::from_rgb(120, 200, 255));
    for e in &drawing.entities {
        match e {
            Entity2D::Line { a, b, .. } => {
                painter.line_segment([to_screen(a), to_screen(b)], stroke);
            }
            Entity2D::Circle { centre, radius, .. } => {
                painter.circle_stroke(to_screen(centre), *radius as f32 * scale, stroke);
            }
            Entity2D::Arc {
                centre,
                radius,
                start_angle_deg,
                end_angle_deg,
                ..
            } => {
                let cc = to_screen(centre);
                let r = *radius as f32 * scale;
                let a0 = start_angle_deg.to_radians() as f32;
                let a1 = end_angle_deg.to_radians() as f32;
                let n = 48;
                let mut prev: Option<egui::Pos2> = None;
                for i in 0..=n {
                    let t = a0 + (a1 - a0) * (i as f32 / n as f32);
                    let p = egui::pos2(cc.x + r * t.cos(), cc.y - r * t.sin());
                    if let Some(pp) = prev {
                        painter.line_segment([pp, p], stroke);
                    }
                    prev = Some(p);
                }
            }
            Entity2D::Polyline {
                vertices, closed, ..
            } => {
                for w in vertices.windows(2) {
                    painter.line_segment([to_screen(&w[0]), to_screen(&w[1])], stroke);
                }
                if *closed && vertices.len() > 2 {
                    painter.line_segment(
                        [to_screen(vertices.last().unwrap()), to_screen(&vertices[0])],
                        stroke,
                    );
                }
            }
            Entity2D::Text {
                position,
                height,
                text,
                ..
            } => {
                painter.text(
                    to_screen(position),
                    egui::Align2::LEFT_BOTTOM,
                    text,
                    egui::FontId::proportional((*height as f32 * scale).clamp(8.0, 32.0)),
                    egui::Color32::from_gray(210),
                );
            }
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn demo_drawing_has_entities() {
        let d = demo_drawing();
        assert!(
            d.entities.len() >= 4,
            "demo has geometry: {}",
            d.entities.len()
        );
    }

    #[test]
    fn drawing_round_trips_through_dxf() {
        let d = demo_drawing();
        let n = d.entities.len();
        let text = dxf::serialise(&d);
        let parsed = dxf::parse(&text).expect("DXF parses");
        assert_eq!(
            parsed.entities.len(),
            n,
            "entity count preserved through DXF round-trip"
        );
    }
}
