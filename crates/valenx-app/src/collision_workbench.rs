//! The right-side **Collision Workbench** panel — native AABB geometry +
//! overlap/separation tests over `valenx-collision`.
//!
//! Mirrors the springs / gears / geomatics / piping workbenches: a
//! resizable [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_collision_workbench`, toggled from the View
//! menu. The form takes two axis-aligned bounding boxes (each a min and a
//! max corner); the "Compute" button reports each box's volume, surface
//! area, space diagonal, and inradius, plus whether the pair overlaps and
//! their L2 separation, as a monospace readout.

use eframe::egui;
use nalgebra::Vector3;

use valenx_collision::{distance, intersect, Aabb};
use valenx_viz::{project_point, OrbitCamera, ViewDirection};

use crate::ValenxApp;

/// Persistent form + result state for the Collision Workbench.
pub struct CollisionWorkbenchState {
    /// Box A minimum corner `[x, y, z]`.
    a_min: [f64; 3],
    /// Box A maximum corner `[x, y, z]`.
    a_max: [f64; 3],
    /// Box B minimum corner `[x, y, z]`.
    b_min: [f64; 3],
    /// Box B maximum corner `[x, y, z]`.
    b_max: [f64; 3],
    /// Formatted readout (empty until the first compute).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
}

impl Default for CollisionWorkbenchState {
    fn default() -> Self {
        // A 10×20×30 box at the origin and a disjoint 10×10×10 box offset
        // along +x (a 10-unit gap), so the default shows both the per-box
        // scalars and a non-trivial separation.
        Self {
            a_min: [0.0, 0.0, 0.0],
            a_max: [10.0, 20.0, 30.0],
            b_min: [20.0, 0.0, 0.0],
            b_max: [30.0, 10.0, 10.0],
            result: String::new(),
            error: None,
        }
    }
}

/// Draw the Collision Workbench right-side panel. A no-op when the
/// `show_collision_workbench` toggle is off.
pub fn draw_collision_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_collision_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_collision_workbench",
        "Collision",
        |app, ui| {
            ui.label(
                egui::RichText::new("native AABB geometry + overlap test · valenx-collision")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.collision;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    corner_rows(ui, "Box A", &mut s.a_min, &mut s.a_max);
                    ui.add_space(4.0);
                    corner_rows(ui, "Box B", &mut s.b_min, &mut s.b_max);

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("\u{25B6} Compute").strong())
                        .clicked()
                    {
                        run_collision(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Geometry").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }

                    // Live AABB wireframes (box A cyan · box B orange, iso).
                    if let Some(edges) = preview_boxes(s) {
                        ui.add_space(6.0);
                        ui.label(
                            egui::RichText::new("Box preview (A cyan \u{00B7} B orange)").strong(),
                        );
                        draw_boxes_preview(ui, &edges);
                    }
                });
        },
    );
    if close {
        app.show_collision_workbench = false;
    }
}

/// Render a labelled `min` row and `max` row of three `DragValue`s for one box.
fn corner_rows(ui: &mut egui::Ui, label: &str, min: &mut [f64; 3], max: &mut [f64; 3]) {
    ui.label(egui::RichText::new(format!("{label} (min / max corners)")).strong());
    ui.horizontal(|ui| {
        ui.label("min");
        for v in min.iter_mut() {
            ui.add(egui::DragValue::new(v).speed(0.5));
        }
    });
    ui.horizontal(|ui| {
        ui.label("max");
        for v in max.iter_mut() {
            ui.add(egui::DragValue::new(v).speed(0.5));
        }
    });
}

/// The 12 edges of an axis-aligned box, each as a pair of corner points.
/// Pure geometry over the eight corners — no validation, never panics
/// (a degenerate `min == max` box yields 12 zero-length edges).
fn box_edges(min: [f64; 3], max: [f64; 3]) -> Vec<[Vector3<f64>; 2]> {
    // Corner `i` takes `max` on axis `k` when bit `k` of `i` is set.
    let corner = |i: usize| {
        Vector3::new(
            if (i & 1) == 0 { min[0] } else { max[0] },
            if (i & 2) == 0 { min[1] } else { max[1] },
            if (i & 4) == 0 { min[2] } else { max[2] },
        )
    };
    // Each edge joins two corners differing in exactly one axis bit:
    // four x-aligned, then four y-aligned, then four z-aligned.
    const EDGES: [(usize, usize); 12] = [
        (0, 1),
        (2, 3),
        (4, 5),
        (6, 7),
        (0, 2),
        (1, 3),
        (4, 6),
        (5, 7),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];
    EDGES.iter().map(|&(a, b)| [corner(a), corner(b)]).collect()
}

/// Build the coloured edge list for both AABBs (box A cyan, box B orange),
/// best-effort `None` for non-finite coordinates (the same precondition as
/// [`run_collision`]). Drives the live iso wireframe preview.
fn preview_boxes(s: &CollisionWorkbenchState) -> Option<Vec<([Vector3<f64>; 2], egui::Color32)>> {
    let all_finite = s
        .a_min
        .iter()
        .chain(s.a_max.iter())
        .chain(s.b_min.iter())
        .chain(s.b_max.iter())
        .all(|v| v.is_finite());
    if !all_finite {
        return None;
    }
    let cyan = egui::Color32::from_rgb(120, 200, 255);
    let orange = egui::Color32::from_rgb(255, 180, 90);
    let mut edges: Vec<([Vector3<f64>; 2], egui::Color32)> = Vec::with_capacity(24);
    edges.extend(box_edges(s.a_min, s.a_max).into_iter().map(|e| (e, cyan)));
    edges.extend(box_edges(s.b_min, s.b_max).into_iter().map(|e| (e, orange)));
    Some(edges)
}

/// Draw a set of coloured 3-D edges as an iso-view wireframe in a
/// fixed-height canvas, the camera auto-framed to all endpoints. A segment
/// is painted only when both endpoints project in front of the camera, so
/// the render path never panics.
fn draw_boxes_preview(ui: &mut egui::Ui, edges: &[([Vector3<f64>; 2], egui::Color32)]) {
    let (response, painter) = ui.allocate_painter(
        egui::vec2(ui.available_width(), 200.0),
        egui::Sense::hover(),
    );
    let rect = response.rect;

    let mut min = [f32::INFINITY; 3];
    let mut max = [f32::NEG_INFINITY; 3];
    for (edge, _) in edges {
        for p in edge {
            for k in 0..3 {
                let v = p[k] as f32;
                min[k] = min[k].min(v);
                max[k] = max[k].max(v);
            }
        }
    }

    let mut cam = OrbitCamera::default();
    cam.set_view(ViewDirection::Iso);
    cam.frame_bounds(min, max);

    let (w, h) = (rect.width(), rect.height());
    for (edge, color) in edges {
        let a = project_point(
            &cam,
            w,
            h,
            [edge[0].x as f32, edge[0].y as f32, edge[0].z as f32],
        );
        let b = project_point(
            &cam,
            w,
            h,
            [edge[1].x as f32, edge[1].y as f32, edge[1].z as f32],
        );
        if let (Some(a), Some(b)) = (a, b) {
            painter.line_segment(
                [
                    egui::pos2(rect.min.x + a.x, rect.min.y + a.y),
                    egui::pos2(rect.min.x + b.x, rect.min.y + b.y),
                ],
                egui::Stroke::new(1.5, *color),
            );
        }
    }
}

/// Build the two [`Aabb`]s from the form, compute their geometry + the
/// pairwise overlap/separation, and format the readout. Extracted from the
/// draw closure so it is unit-testable.
fn run_collision(s: &mut CollisionWorkbenchState) {
    s.error = None;

    for v in s
        .a_min
        .iter()
        .chain(s.a_max.iter())
        .chain(s.b_min.iter())
        .chain(s.b_max.iter())
    {
        if !v.is_finite() {
            s.error = Some("all box coordinates must be finite".into());
            return;
        }
    }

    let a = Aabb {
        min: Vector3::new(s.a_min[0], s.a_min[1], s.a_min[2]),
        max: Vector3::new(s.a_max[0], s.a_max[1], s.a_max[2]),
    };
    let b = Aabb {
        min: Vector3::new(s.b_min[0], s.b_min[1], s.b_min[2]),
        max: Vector3::new(s.b_max[0], s.b_max[1], s.b_max[2]),
    };

    let hit = intersect(&a, &b);
    let sep = distance(&a, &b);

    s.result = format!(
        "box A : vol {:.3}  surf {:.3}  diag {:.4}  inr {:.3}\n\
         box B : vol {:.3}  surf {:.3}  diag {:.4}  inr {:.3}\n\n\
         intersect  : {}\n\
         separation : {:.4}   (L2 gap; 0 when overlapping)",
        a.volume(),
        a.surface_area(),
        a.diagonal(),
        a.inradius(),
        b.volume(),
        b.surface_area(),
        b.diagonal(),
        b.inradius(),
        if hit {
            "yes \u{2014} boxes overlap"
        } else {
            "no \u{2014} disjoint"
        },
        sep,
    );
}

/// Build the **Collision (AABB)** result card for the Workbench+Agent bridge — a
/// DATA-ONLY [`crate::WorkspaceProduct`] (`mesh: None`) whose `lines` are the
/// genuine AABB overlap test + separation ([`run_collision`]) for the canonical
/// default pair of boxes (the default is a disjoint pair with a 10-unit gap).
/// Registered as the `"collision"` producer in
/// [`crate::products_registry::lookup`]; the tile renders it as a text card, not
/// a 3-D view.
pub(crate) fn collision_product() -> crate::WorkspaceProduct {
    let mut s = CollisionWorkbenchState::default();
    run_collision(&mut s);
    crate::WorkspaceProduct {
        title: "Collision (AABB)".into(),
        lines: crate::products_registry::lines_from_readout(&s.result),
        mesh: None,
        vertex_colors: None,
        camera: Default::default(),
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = CollisionWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn compute_default_disjoint_boxes() {
        let mut s = CollisionWorkbenchState::default();
        run_collision(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        assert!(s.result.contains("box A"));
        assert!(s.result.contains("intersect"));
        assert!(s.result.contains("separation"));
        // Recompute via the backend to confirm the defaults: A = 10×20×30
        // → volume 6000; the boxes are disjoint with a 10-unit x-gap.
        let a = Aabb {
            min: Vector3::new(0.0, 0.0, 0.0),
            max: Vector3::new(10.0, 20.0, 30.0),
        };
        let b = Aabb {
            min: Vector3::new(20.0, 0.0, 0.0),
            max: Vector3::new(30.0, 10.0, 10.0),
        };
        assert!((a.volume() - 6000.0).abs() < 1e-9);
        assert!(!intersect(&a, &b));
        assert!((distance(&a, &b) - 10.0).abs() < 1e-9);
        // The readout reflects them.
        assert!(s.result.contains("6000.000"));
        assert!(s.result.contains("disjoint"));
    }

    #[test]
    fn compute_overlapping_boxes_report_zero_separation() {
        // Box B now straddles A → they share volume.
        let mut s = CollisionWorkbenchState {
            b_min: [5.0, 5.0, 5.0],
            b_max: [15.0, 15.0, 15.0],
            ..Default::default()
        };
        run_collision(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("overlap"));
        assert!(s.result.contains("separation : 0.0000"));
    }

    #[test]
    fn compute_rejects_non_finite_coords() {
        let mut s = CollisionWorkbenchState {
            a_max: [f64::NAN, 20.0, 30.0],
            ..Default::default()
        };
        run_collision(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());
    }

    #[test]
    fn box_edges_unit_cube_has_twelve_unit_edges() {
        let edges = box_edges([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(edges.len(), 12);
        // Every edge of the unit cube spans exactly one unit.
        for [p, q] in &edges {
            assert!(((q - p).norm() - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn box_edges_degenerate_is_twelve_zero_length_edges() {
        // min == max → all eight corners coincide; no panic, 12 zero edges.
        let edges = box_edges([2.0, 2.0, 2.0], [2.0, 2.0, 2.0]);
        assert_eq!(edges.len(), 12);
        assert!(edges.iter().all(|[p, q]| (q - p).norm() == 0.0));
    }

    #[test]
    fn preview_boxes_default_has_two_coloured_boxes() {
        let s = CollisionWorkbenchState::default();
        let edges = preview_boxes(&s).expect("finite default boxes preview");
        // 12 edges per box × 2 boxes.
        assert_eq!(edges.len(), 24);
        let cyan = egui::Color32::from_rgb(120, 200, 255);
        let orange = egui::Color32::from_rgb(255, 180, 90);
        assert_eq!(edges.iter().filter(|(_, c)| *c == cyan).count(), 12);
        assert_eq!(edges.iter().filter(|(_, c)| *c == orange).count(), 12);
        // The combined wireframe spans both boxes: A.x ∈ [0,10], B.x ∈ [20,30].
        let (mut xmin, mut xmax) = (f64::INFINITY, f64::NEG_INFINITY);
        for ([p, q], _) in &edges {
            for v in [p.x, q.x] {
                xmin = xmin.min(v);
                xmax = xmax.max(v);
            }
        }
        assert!(xmin.abs() < 1e-9);
        assert!((xmax - 30.0).abs() < 1e-9);
    }

    #[test]
    fn preview_boxes_none_for_non_finite() {
        let s = CollisionWorkbenchState {
            a_max: [f64::NAN, 20.0, 30.0],
            ..Default::default()
        };
        assert!(preview_boxes(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    /// Render the whole workbench panel once in a headless egui context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_collision_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_collision_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_collision_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_collision_workbench = true;
        run_collision(&mut app.collision);
        app.collision.error = Some("invalid box".to_string());
        draw_workbench(&mut app);
    }
}
