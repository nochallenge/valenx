//! The right-side **Gears Workbench** panel — native involute-gear design
//! over `valenx-gears`.
//!
//! Mirrors the springs / CFD / FEM workbenches: a resizable
//! [`egui::SidePanel`] gated on [`crate::ValenxApp::show_gears_workbench`],
//! toggled from the View menu. The form drives a [`valenx_gears::GearSpec`];
//! the "Analyze" button computes the design scalars — circular pitch and the
//! pitch / base / addendum / dedendum diameters — plus the meshing gear ratio
//! against a mating tooth count, and renders them as a monospace readout.

use eframe::egui;

use valenx_gears::{circular_pitch_mm, gear_ratio, GearKind, GearSpec};

use crate::ValenxApp;

/// Persistent form + result state for the Gears Workbench.
pub struct GearsWorkbenchState {
    /// Gear family (spur / helical / bevel / worm).
    kind: GearKind,
    /// Module `m` (mm) — pitch diameter = `m × teeth`.
    module_mm: f64,
    /// Tooth count `z`.
    teeth: u32,
    /// Pressure angle (degrees) — standard 20°.
    pressure_angle_deg: f64,
    /// Helix angle (degrees) — 0 for spur, ~20–30 for helical.
    helix_angle_deg: f64,
    /// Face width (mm).
    face_width_mm: f64,
    /// Mating gear's tooth count, for the meshing ratio.
    mate_teeth: u32,
    /// Formatted design readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
}

impl Default for GearsWorkbenchState {
    fn default() -> Self {
        // Mirror valenx-gears' standard 20-tooth, 1-module, 20° spur gear.
        Self {
            kind: GearKind::Spur,
            module_mm: 1.0,
            teeth: 20,
            pressure_angle_deg: 20.0,
            helix_angle_deg: 0.0,
            face_width_mm: 10.0,
            mate_teeth: 40,
            result: String::new(),
            error: None,
        }
    }
}

/// Draw the Gears Workbench right-side panel. A no-op when the
/// `show_gears_workbench` toggle is off.
pub fn draw_gears_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_gears_workbench {
        return;
    }

    egui::SidePanel::right("valenx_gears_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("Gears");
            ui.label(
                egui::RichText::new("native involute-gear design · valenx-gears")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.gears;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Family").strong());
                    ui.horizontal_wrapped(|ui| {
                        ui.radio_value(&mut s.kind, GearKind::Spur, "Spur");
                        ui.radio_value(&mut s.kind, GearKind::Helical, "Helical");
                        ui.radio_value(&mut s.kind, GearKind::Bevel, "Bevel");
                        ui.radio_value(&mut s.kind, GearKind::Worm, "Worm");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Geometry").strong());
                    ui.horizontal(|ui| {
                        ui.label("module m (mm)");
                        ui.add(egui::DragValue::new(&mut s.module_mm).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("teeth z");
                        ui.add(egui::DragValue::new(&mut s.teeth).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("pressure angle (°)");
                        ui.add(egui::DragValue::new(&mut s.pressure_angle_deg).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("helix angle (°)");
                        ui.add(egui::DragValue::new(&mut s.helix_angle_deg).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("face width (mm)");
                        ui.add(egui::DragValue::new(&mut s.face_width_mm).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Mesh").strong());
                    ui.horizontal(|ui| {
                        ui.label("mating teeth");
                        ui.add(egui::DragValue::new(&mut s.mate_teeth).speed(1.0));
                    });

                    // Live hint: the pitch (reference) diameter d = m·z.
                    if s.module_mm > 0.0 {
                        let pd = s.module_mm * s.teeth as f64;
                        ui.label(
                            egui::RichText::new(format!("pitch diameter d ≈ {pd:.2} mm  (m·z)"))
                                .weak()
                                .small(),
                        );
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_gears(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Design").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });
}

/// Build a [`GearSpec`] from the form, validate it, and format the
/// design-scalar readout. Extracted from the draw closure so it is
/// unit-testable.
fn run_gears(s: &mut GearsWorkbenchState) {
    s.error = None;

    if !(s.module_mm.is_finite() && s.module_mm > 0.0) {
        s.error = Some("module must be positive".into());
        return;
    }
    if s.teeth == 0 {
        s.error = Some("tooth count must be at least 1".into());
        return;
    }
    // The base circle d·cos(α) needs 0 < α < 90 (cos α > 0).
    if !(s.pressure_angle_deg.is_finite() && s.pressure_angle_deg > 0.0 && s.pressure_angle_deg < 90.0)
    {
        s.error = Some("pressure angle must be between 0° and 90°".into());
        return;
    }
    if !(s.helix_angle_deg.is_finite() && s.helix_angle_deg >= 0.0 && s.helix_angle_deg < 90.0) {
        s.error = Some("helix angle must be between 0° and 90°".into());
        return;
    }
    if !(s.face_width_mm.is_finite() && s.face_width_mm > 0.0) {
        s.error = Some("face width must be positive".into());
        return;
    }
    if s.mate_teeth == 0 {
        s.error = Some("mating tooth count must be at least 1".into());
        return;
    }

    let spec = GearSpec {
        kind: s.kind,
        module_mm: s.module_mm,
        teeth: s.teeth,
        pressure_angle_deg: s.pressure_angle_deg,
        helix_angle_deg: s.helix_angle_deg,
        face_width_mm: s.face_width_mm,
    };

    let p = circular_pitch_mm(s.module_mm);
    // Convention: this gear is the driver, the mating gear the driven, so the
    // ratio is mate ÷ this (> 1 reduces speed / multiplies torque).
    let ratio = gear_ratio(s.mate_teeth, s.teeth);

    s.result = format!(
        "family        : {}\n\
         module m      : {:.3} mm\n\
         teeth z       : {}\n\
         pressure angle: {:.2}\u{00B0}\n\
         helix angle   : {:.2}\u{00B0}\n\
         face width    : {:.3} mm\n\
         mating teeth  : {}\n\n\
         circular pitch: {:.4} mm  (\u{03C0}\u{00B7}m)\n\
         pitch diameter: {:.3} mm  (m\u{00B7}z)\n\
         base diameter : {:.3} mm  (d\u{00B7}cos \u{03B1})\n\
         addendum dia. : {:.3} mm  (d + 2m)\n\
         dedendum dia. : {:.3} mm  (d \u{2212} 2.5m)\n\
         gear ratio    : {:.3}  (mate \u{00F7} this; >1 reduces speed)",
        s.kind.label(),
        s.module_mm,
        s.teeth,
        s.pressure_angle_deg,
        s.helix_angle_deg,
        s.face_width_mm,
        s.mate_teeth,
        p,
        spec.pitch_diameter_mm(),
        spec.base_diameter_mm(),
        spec.addendum_diameter_mm(),
        spec.dedendum_diameter_mm(),
        ratio,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = GearsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_spur_gear() {
        let mut s = GearsWorkbenchState::default();
        run_gears(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        // The readout names the core design scalars.
        assert!(s.result.contains("circular pitch"));
        assert!(s.result.contains("pitch diameter"));
        assert!(s.result.contains("base diameter"));
        assert!(s.result.contains("gear ratio"));
        // Standard m=1, z=20: pitch diameter d = m·z = 20.
        assert!(s.result.contains("20.000 mm"));
        // base = d·cos20° = 20 × 0.9396926 ≈ 18.794.
        assert!(s.result.contains("18.794"));
        // ratio = mate ÷ this = 40 ÷ 20 = 2.
        assert!(s.result.contains("2.000"));
    }

    #[test]
    fn analyze_rejects_zero_teeth() {
        let mut s = GearsWorkbenchState {
            teeth: 0,
            ..Default::default()
        };
        run_gears(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());
    }

    #[test]
    fn analyze_rejects_bad_module_angle_and_mate() {
        for bad in [
            GearsWorkbenchState {
                module_mm: 0.0,
                ..Default::default()
            },
            GearsWorkbenchState {
                pressure_angle_deg: 90.0,
                ..Default::default()
            },
            GearsWorkbenchState {
                face_width_mm: 0.0,
                ..Default::default()
            },
            GearsWorkbenchState {
                mate_teeth: 0,
                ..Default::default()
            },
        ] {
            let mut s = bad;
            run_gears(&mut s);
            assert!(s.error.is_some());
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
            draw_gears_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_gears_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_gears_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_gears_workbench = true;
        run_gears(&mut app.gears);
        app.gears.error = Some("invalid gear parameters".to_string());
        draw_workbench(&mut app);
    }
}
