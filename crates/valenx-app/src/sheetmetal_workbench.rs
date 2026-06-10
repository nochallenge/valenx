//! The right-side **Sheet Metal Workbench** panel — bend allowance and
//! bend deduction over `valenx-sheet-metal`.
//!
//! Mirrors the springs / gears / geomatics / piping / collision
//! workbenches: a resizable [`egui::SidePanel`] gated on
//! [`crate::ValenxApp::show_sheetmetal_workbench`], toggled from the View
//! menu. The form takes a bend's thickness, inside radius, angle, and
//! k-factor; the "Compute" button reports the neutral-axis radius, the
//! bend allowance (the neutral-axis arc length), and the bend deduction
//! (the flat-blank correction), as a monospace readout.

use eframe::egui;

use valenx_sheet_metal::Bend;

use crate::ValenxApp;

/// Persistent form + result state for the Sheet Metal Workbench.
pub struct SheetmetalWorkbenchState {
    /// Sheet thickness (mm).
    thickness_mm: f64,
    /// Inside (concave-side) bend radius (mm).
    inside_radius_mm: f64,
    /// Bend angle (degrees) — the swept angle of the bend.
    bend_angle_deg: f64,
    /// K-factor — the neutral-axis fraction (0.33–0.5 typical).
    k_factor: f64,
    /// Formatted bend readout (empty until the first compute).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
}

impl Default for SheetmetalWorkbenchState {
    fn default() -> Self {
        // A 90° bend in 1 mm stock, 1 mm inside radius, k = 0.44 — the
        // canonical worked example (bend allowance 1.44·π/2 ≈ 2.2619 mm).
        Self {
            thickness_mm: 1.0,
            inside_radius_mm: 1.0,
            bend_angle_deg: 90.0,
            k_factor: 0.44,
            result: String::new(),
            error: None,
        }
    }
}

/// Draw the Sheet Metal Workbench right-side panel. A no-op when the
/// `show_sheetmetal_workbench` toggle is off.
pub fn draw_sheetmetal_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_sheetmetal_workbench {
        return;
    }

    egui::SidePanel::right("valenx_sheetmetal_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("Sheet Metal");
            ui.label(
                egui::RichText::new("native bend allowance / deduction · valenx-sheet-metal")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.sheetmetal;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Bend parameters").strong());
                    ui.horizontal(|ui| {
                        ui.label("thickness (mm)");
                        ui.add(egui::DragValue::new(&mut s.thickness_mm).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("inside radius (mm)");
                        ui.add(egui::DragValue::new(&mut s.inside_radius_mm).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("bend angle (\u{00B0})");
                        ui.add(egui::DragValue::new(&mut s.bend_angle_deg).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("k-factor");
                        ui.add(egui::DragValue::new(&mut s.k_factor).speed(0.01));
                    });
                    ui.label(
                        egui::RichText::new("k-factor = neutral-axis fraction (0.33–0.5 typical)")
                            .weak()
                            .small(),
                    );

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("\u{25B6} Compute").strong())
                        .clicked()
                    {
                        run_sheetmetal(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Bend").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });
}

/// Build a [`Bend`] from the form, compute the bend allowance / deduction,
/// and format the readout. Extracted from the draw closure so it is
/// unit-testable. The bend-line endpoints do not affect either scalar, so
/// a unit line is used.
fn run_sheetmetal(s: &mut SheetmetalWorkbenchState) {
    s.error = None;

    if !(s.thickness_mm.is_finite() && s.thickness_mm > 0.0) {
        s.error = Some("thickness must be positive".into());
        return;
    }
    if !(s.inside_radius_mm.is_finite() && s.inside_radius_mm >= 0.0) {
        s.error = Some("inside radius must be \u{2265} 0".into());
        return;
    }
    if !(s.bend_angle_deg.is_finite() && s.bend_angle_deg > 0.0 && s.bend_angle_deg <= 180.0) {
        s.error = Some("bend angle must be in (0, 180]".into());
        return;
    }
    if !(s.k_factor.is_finite() && (0.0..=1.0).contains(&s.k_factor)) {
        s.error = Some("k-factor must be in [0, 1]".into());
        return;
    }

    let bend = Bend::new(
        [0.0, 0.0],
        [1.0, 0.0],
        s.bend_angle_deg.to_radians(),
        s.inside_radius_mm,
    );
    let r_neutral = s.inside_radius_mm + s.k_factor * s.thickness_mm;
    let ba = bend.bend_allowance(s.thickness_mm, s.k_factor);
    let bd = bend.bend_deduction(s.thickness_mm, s.k_factor);

    s.result = format!(
        "thickness      : {:.3} mm\n\
         inside radius  : {:.3} mm\n\
         bend angle     : {:.1}\u{00B0}\n\
         k-factor       : {:.3}\n\n\
         neutral radius : {:.4} mm   (r_i + k\u{00B7}t)\n\
         bend allowance : {:.4} mm   (r_n\u{00B7}\u{03B8}; neutral-axis arc)\n\
         bend deduction : {:.4} mm   (2\u{00B7}OSSB \u{2212} BA; flat-blank)",
        s.thickness_mm,
        s.inside_radius_mm,
        s.bend_angle_deg,
        s.k_factor,
        r_neutral,
        ba,
        bd,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = SheetmetalWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn compute_default_90deg_bend() {
        let mut s = SheetmetalWorkbenchState::default();
        run_sheetmetal(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        assert!(s.result.contains("neutral radius"));
        assert!(s.result.contains("bend allowance"));
        assert!(s.result.contains("bend deduction"));
        // 90° bend, r_i = 1, t = 1, k = 0.44: r_n = 1.44, BA = 1.44·π/2 ≈
        // 2.2619 mm. Recompute via the backend to confirm the default.
        let bend = Bend::new([0.0, 0.0], [1.0, 0.0], 90_f64.to_radians(), 1.0);
        assert!((bend.bend_allowance(1.0, 0.44) - 2.261_95).abs() < 1e-3);
        assert!(s.result.contains("2.2619"));
    }

    #[test]
    fn compute_rejects_bad_inputs() {
        for bad in [
            SheetmetalWorkbenchState {
                thickness_mm: 0.0,
                ..Default::default()
            },
            SheetmetalWorkbenchState {
                bend_angle_deg: 0.0,
                ..Default::default()
            },
            SheetmetalWorkbenchState {
                bend_angle_deg: 200.0,
                ..Default::default()
            },
            SheetmetalWorkbenchState {
                k_factor: 1.5,
                ..Default::default()
            },
        ] {
            let mut s = bad;
            run_sheetmetal(&mut s);
            assert!(s.error.is_some());
            assert!(s.result.is_empty());
        }
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
            draw_sheetmetal_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_sheetmetal_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_sheetmetal_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_sheetmetal_workbench = true;
        run_sheetmetal(&mut app.sheetmetal);
        app.sheetmetal.error = Some("invalid bend parameters".to_string());
        draw_workbench(&mut app);
    }
}
