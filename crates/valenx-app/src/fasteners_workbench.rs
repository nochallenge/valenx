//! The right-side **Fasteners Workbench** panel — ISO 4017 hex-bolt
//! dimensions over `valenx-fasteners`.
//!
//! Mirrors the springs / gears / … / fields workbenches: a resizable
//! [`egui::SidePanel`] gated on [`crate::ValenxApp::show_fasteners_workbench`],
//! toggled from the View menu. The form picks a standard metric bolt size
//! from the ISO 4017 hex table; the "Compute" button reports the width
//! across flats, head height, pitch diameter, and tensile stress area, as
//! a monospace readout.

use eframe::egui;

use valenx_fasteners::bolt::iso4017_hex_table;

use crate::ValenxApp;

/// Persistent form + result state for the Fasteners Workbench.
pub struct FastenersWorkbenchState {
    /// Selected ISO metric nominal designation (e.g. `"M6"`).
    nominal: String,
    /// Formatted dimension readout (empty until the first compute).
    result: String,
    /// Validation / lookup error, if any.
    error: Option<String>,
}

impl Default for FastenersWorkbenchState {
    fn default() -> Self {
        Self {
            nominal: "M6".to_string(),
            result: String::new(),
            error: None,
        }
    }
}

/// Draw the Fasteners Workbench right-side panel. A no-op when the
/// `show_fasteners_workbench` toggle is off.
pub fn draw_fasteners_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_fasteners_workbench {
        return;
    }

    egui::SidePanel::right("valenx_fasteners_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("Fasteners");
            ui.label(
                egui::RichText::new("ISO 4017 hex-bolt dimensions · valenx-fasteners")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.fasteners;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Bolt size (ISO 4017 hex)").strong());
                    let table = iso4017_hex_table();
                    egui::ComboBox::from_id_source("valenx_fasteners_nominal")
                        .selected_text(s.nominal.clone())
                        .show_ui(ui, |ui| {
                            for spec in &table {
                                ui.selectable_value(
                                    &mut s.nominal,
                                    spec.nominal.clone(),
                                    spec.nominal.as_str(),
                                );
                            }
                        });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("\u{25B6} Compute").strong())
                        .clicked()
                    {
                        run_fasteners(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Dimensions").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });
}

/// Look up the selected bolt in the ISO 4017 table and format its
/// dimensions. Extracted from the draw closure so it is unit-testable.
fn run_fasteners(s: &mut FastenersWorkbenchState) {
    s.error = None;

    let table = iso4017_hex_table();
    let spec = match table.iter().find(|b| b.nominal == s.nominal) {
        Some(b) => b,
        None => {
            s.error = Some(format!("unknown bolt size '{}'", s.nominal));
            return;
        }
    };

    s.result = format!(
        "designation         : {}  (ISO 4017 hex)\n\
         width across flats  : {:.3} mm\n\
         head height         : {:.3} mm\n\
         pitch diameter      : {:.4} mm\n\
         tensile stress area : {:.3} mm\u{00B2}",
        s.nominal,
        spec.width_across_flats_mm(),
        spec.head_height_mm(),
        spec.pitch_diameter_mm(),
        spec.tensile_stress_area_mm2(),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = FastenersWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn compute_default_m6() {
        let mut s = FastenersWorkbenchState::default();
        run_fasteners(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        assert!(s.result.contains("width across flats"));
        assert!(s.result.contains("pitch diameter"));
        assert!(s.result.contains("tensile stress area"));
        // M6: pitch diameter 6 − 0.6495·1 = 5.3505 mm; tensile stress area ≈
        // 20.12 mm². Recompute via the backend table to confirm the default.
        let m6 = iso4017_hex_table()
            .into_iter()
            .find(|b| b.nominal == "M6")
            .expect("M6 is in the ISO 4017 table");
        assert!((m6.pitch_diameter_mm() - 5.3505).abs() < 1e-4);
        assert!((m6.tensile_stress_area_mm2() - 20.12).abs() < 0.05);
        assert!(s.result.contains("5.3505"));
    }

    #[test]
    fn compute_rejects_unknown_size() {
        let mut s = FastenersWorkbenchState {
            nominal: "M999".to_string(),
            ..Default::default()
        };
        run_fasteners(&mut s);
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
            draw_fasteners_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_fasteners_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fasteners_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fasteners_workbench = true;
        run_fasteners(&mut app.fasteners);
        app.fasteners.error = Some("invalid bolt".to_string());
        draw_workbench(&mut app);
    }
}
