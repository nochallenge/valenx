//! The right-side **Collision Workbench** panel — native AABB geometry +
//! overlap/separation tests over `valenx-collision`.
//!
//! Mirrors the springs / gears / geomatics / piping workbenches: a
//! resizable [`egui::SidePanel`] gated on
//! [`crate::ValenxApp::show_collision_workbench`], toggled from the View
//! menu. The form takes two axis-aligned bounding boxes (each a min and a
//! max corner); the "Compute" button reports each box's volume, surface
//! area, space diagonal, and inradius, plus whether the pair overlaps and
//! their L2 separation, as a monospace readout.

use eframe::egui;
use nalgebra::Vector3;

use valenx_collision::{distance, intersect, Aabb};

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

    egui::SidePanel::right("valenx_collision_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("Collision");
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
                });
        });
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
