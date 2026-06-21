//! Interior-design workbench — rooms + furniture on `valenx-interior`.
//!
//! A right-side 2-D **floor-plan** canvas: a room polygon drawn as walls, with
//! furniture placements rendered as labelled rectangles. A palette selects the
//! furniture kind; "Place" drops it at the room centre. The scene model
//! (`InteriorPanelState`) is headless-testable; the canvas is the view.

use eframe::egui;

use nalgebra::{Vector2, Vector3};
use valenx_interior::{Furniture, InteriorPanelState, Room};

use crate::ValenxApp;

/// Persistent state for the interior workbench.
pub struct InteriorWorkbenchState {
    scene: InteriorPanelState,
    pan: egui::Vec2,
    scale: f32,
    selected: Furniture,
}

impl Default for InteriorWorkbenchState {
    fn default() -> Self {
        Self {
            scene: demo_scene(),
            pan: egui::Vec2::new(3.0, 2.0),
            scale: 48.0,
            selected: Furniture::Table,
        }
    }
}

/// A 6 m × 4 m demo room with a couple of pieces placed.
fn demo_scene() -> InteriorPanelState {
    let mut scene = InteriorPanelState::new();
    let mut room = Room::new("room", "Living Room", 2.5);
    room.floor_polygon = vec![
        Vector2::new(0.0, 0.0),
        Vector2::new(6.0, 0.0),
        Vector2::new(6.0, 4.0),
        Vector2::new(0.0, 4.0),
    ];
    let _ = scene.add_room(room);
    scene.select(Furniture::Sofa);
    let _ = scene.click_to_place(Vector3::new(1.5, 1.0, 0.0), "room");
    scene.select(Furniture::Table);
    let _ = scene.click_to_place(Vector3::new(3.0, 2.0, 0.0), "room");
    scene
}

/// Centroid of the first room's floor polygon (drawing units).
fn room_centre(scene: &InteriorPanelState) -> Vector3<f64> {
    if let Some(room) = scene.rooms.first() {
        let poly = &room.floor_polygon;
        if !poly.is_empty() {
            let n = poly.len() as f64;
            let cx = poly.iter().map(|p| p.x).sum::<f64>() / n;
            let cy = poly.iter().map(|p| p.y).sum::<f64>() / n;
            return Vector3::new(cx, cy, 0.0);
        }
    }
    Vector3::new(0.0, 0.0, 0.0)
}

/// Convert an [`InteriorPanelState`] into the plain-data
/// [`crate::FloorPlanView`] the workspace-tile painter consumes — room wall
/// polygons + furniture footprints (centre + size + label) in metres, plus the
/// overall room extent so the painter can fit the plan to the tile.
fn scene_to_view(scene: &InteriorPanelState) -> crate::FloorPlanView {
    let mut rooms: Vec<Vec<[f64; 2]>> = Vec::new();
    let mut min = [f64::INFINITY; 2];
    let mut max = [f64::NEG_INFINITY; 2];
    for room in &scene.rooms {
        let poly: Vec<[f64; 2]> = room
            .floor_polygon
            .iter()
            .map(|p| {
                min[0] = min[0].min(p.x);
                min[1] = min[1].min(p.y);
                max[0] = max[0].max(p.x);
                max[1] = max[1].max(p.y);
                [p.x, p.y]
            })
            .collect();
        if !poly.is_empty() {
            rooms.push(poly);
        }
    }
    let furniture: Vec<crate::FloorPlanItem> = scene
        .placements
        .iter()
        .map(|p| {
            let size = p.kind.default_size();
            crate::FloorPlanItem {
                centre: [p.position.x, p.position.y],
                size: [size.x, size.y],
                label: p.kind.name().to_string(),
            }
        })
        .collect();
    if !min[0].is_finite() {
        min = [0.0, 0.0];
        max = [1.0, 1.0];
    }
    crate::FloorPlanView {
        rooms,
        furniture,
        bounds: (min, max),
        lines: vec![
            format!(
                "{} room · {} pieces",
                scene.rooms.len(),
                scene.placements.len()
            ),
            format!(
                "{:.1} × {:.1} m",
                (max[0] - min[0]).abs(),
                (max[1] - min[1]).abs()
            ),
            "valenx-interior · floor plan".into(),
        ],
    }
}

/// Build the agent-bridge **`interior` product** — the canonical demo floor
/// plan (a 6 m × 4 m living room with a sofa + table) as a
/// [`crate::Workspace2dKind::FloorPlan`] painted by the tile's 2-D branch
/// (mirroring `rcbeam` / `dna`). A 2-D drawing, NOT a mesh: `mesh: None`,
/// `kind2d: Some(FloorPlan(..))`.
pub(crate) fn interior_product() -> crate::WorkspaceProduct {
    let scene = demo_scene();
    let view = scene_to_view(&scene);
    let lines = view.lines.clone();
    crate::WorkspaceProduct {
        title: "Floor Plan".into(),
        lines,
        mesh: None,
        vertex_colors: None,
        camera: valenx_viz::OrbitCamera::default(),
        kind2d: Some(crate::Workspace2dKind::FloorPlan(view)),
        last_export: None,
        image: None,
        image_texture: None,
    }
}

/// Draw the interior workbench (a no-op unless toggled on via View → Interior).
pub fn draw_interior_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_interior_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_interior_workbench",
        "Interior Design",
        |app, ui| {
            ui.label(
                egui::RichText::new("floor plan + furniture · valenx-interior")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.interior;
            ui.label(egui::RichText::new("Furniture").strong());
            ui.horizontal_wrapped(|ui| {
                for &f in Furniture::all() {
                    ui.selectable_value(&mut s.selected, f, f.name());
                }
            });
            ui.horizontal(|ui| {
                if ui.button("Place at centre").clicked() {
                    let rid = s.scene.rooms.first().map(|r| r.id.clone());
                    if let Some(rid) = rid {
                        let centre = room_centre(&s.scene);
                        let f = s.selected;
                        s.scene.select(f);
                        let _ = s.scene.click_to_place(centre, &rid);
                    }
                }
                if ui.button("Reset").clicked() {
                    s.scene = demo_scene();
                }
                ui.label(format!("{} pieces", s.scene.placements.len()));
            });
            ui.separator();

            let (resp, painter) = ui.allocate_painter(
                egui::vec2(ui.available_width(), ui.available_height().max(220.0)),
                egui::Sense::drag(),
            );
            let rect = resp.rect;
            painter.rect_filled(rect, 2.0, egui::Color32::from_gray(18));
            if resp.dragged() {
                let d = resp.drag_delta();
                s.pan -= egui::vec2(d.x / s.scale, -d.y / s.scale);
            }
            draw_plan(&painter, rect, s.pan, s.scale, &s.scene);
        },
    );
    if close {
        app.show_interior_workbench = false;
    }
}

fn draw_plan(
    painter: &egui::Painter,
    rect: egui::Rect,
    pan: egui::Vec2,
    scale: f32,
    scene: &InteriorPanelState,
) {
    let c = rect.center();
    let to_screen = |x: f64, y: f64| -> egui::Pos2 {
        egui::pos2(
            c.x + (x as f32 - pan.x) * scale,
            c.y - (y as f32 - pan.y) * scale,
        )
    };
    let wall = egui::Stroke::new(2.0, egui::Color32::from_rgb(180, 190, 210));
    for room in &scene.rooms {
        let poly = &room.floor_polygon;
        for w in poly.windows(2) {
            painter.line_segment([to_screen(w[0].x, w[0].y), to_screen(w[1].x, w[1].y)], wall);
        }
        if poly.len() > 2 {
            let last = poly[poly.len() - 1];
            painter.line_segment(
                [to_screen(last.x, last.y), to_screen(poly[0].x, poly[0].y)],
                wall,
            );
        }
    }
    for p in &scene.placements {
        let size = p.kind.default_size();
        let centre = to_screen(p.position.x, p.position.y);
        let r = egui::Rect::from_center_size(
            centre,
            egui::vec2(size.x as f32 * scale, size.y as f32 * scale),
        );
        painter.rect_filled(r, 2.0, egui::Color32::from_rgb(70, 110, 90));
        painter.rect_stroke(
            r,
            2.0,
            egui::Stroke::new(1.0, egui::Color32::from_gray(200)),
        );
        painter.text(
            centre,
            egui::Align2::CENTER_CENTER,
            p.kind.name(),
            egui::FontId::proportional(11.0),
            egui::Color32::from_gray(220),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn demo_scene_has_a_room_and_furniture() {
        let scene = demo_scene();
        assert_eq!(scene.rooms.len(), 1);
        assert!(scene.placements.len() >= 2, "demo places furniture");
    }

    #[test]
    fn placing_furniture_grows_the_scene() {
        let mut scene = demo_scene();
        let before = scene.placements.len();
        scene.select(Furniture::Chair);
        scene
            .click_to_place(Vector3::new(2.0, 2.0, 0.0), "room")
            .expect("place");
        assert_eq!(scene.placements.len(), before + 1);
    }
}
