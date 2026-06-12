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

/// Draw the 2D drafting workbench (a no-op unless toggled on via
/// View → 2D Drafting).
pub fn draw_draft2d_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_draft2d_workbench {
        return;
    }
    egui::SidePanel::right("valenx_draft2d_workbench")
        .resizable(true)
        .default_width(440.0)
        .width_range(360.0..=760.0)
        .show(ctx, |ui| {
            ui.heading("2D Drafting");
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
        });
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
