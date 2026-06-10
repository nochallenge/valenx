//! The right-side **Piping Workbench** panel — native pipe-section sizing
//! over `valenx-piping`.
//!
//! Mirrors the springs / gears / geomatics workbenches: a resizable
//! [`egui::SidePanel`] gated on [`crate::ValenxApp::show_piping_workbench`],
//! toggled from the View menu. The form picks an NPS designation, a wall
//! schedule, a material, and a length; the "Compute" button reports the
//! outer / inner diameters and the flow / metal cross-section areas, the
//! wetted perimeter, and the external surface area, as a monospace readout.

use eframe::egui;
use nalgebra::Vector3;

use valenx_piping::{Material, PipeSection, PipingError, Schedule};

use crate::ValenxApp;

/// The NPS designations the OD lookup table knows (`valenx_piping::dims`).
const NPS_OPTIONS: &[&str] = &[
    "1/8", "1/4", "3/8", "1/2", "3/4", "1", "1-1/4", "1-1/2", "2", "2-1/2", "3", "3-1/2", "4", "5",
    "6", "8", "10", "12", "14", "16", "18", "20", "24", "30", "36",
];

/// Persistent form + result state for the Piping Workbench.
pub struct PipingWorkbenchState {
    /// NPS designation (e.g. `"2"`), keyed into the OD table.
    nominal_size: String,
    /// Wall schedule (Sch 40 / 80 / 160).
    schedule: Schedule,
    /// Pipe material (BOM metadata; does not change the geometry).
    material: Material,
    /// Section length (mm).
    length_mm: f64,
    /// Formatted section readout (empty until the first compute).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
}

impl Default for PipingWorkbenchState {
    fn default() -> Self {
        Self {
            nominal_size: "2".to_string(),
            schedule: Schedule::Sch40,
            material: Material::CarbonSteel,
            length_mm: 1000.0,
            result: String::new(),
            error: None,
        }
    }
}

/// Draw the Piping Workbench right-side panel. A no-op when the
/// `show_piping_workbench` toggle is off.
pub fn draw_piping_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_piping_workbench {
        return;
    }

    egui::SidePanel::right("valenx_piping_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("Piping");
            ui.label(
                egui::RichText::new("native pipe-section sizing · valenx-piping")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.piping;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Nominal pipe size").strong());
                    egui::ComboBox::from_id_source("valenx_piping_nps")
                        .selected_text(format!("NPS {}", s.nominal_size))
                        .show_ui(ui, |ui| {
                            for nps in NPS_OPTIONS {
                                ui.selectable_value(&mut s.nominal_size, (*nps).to_string(), *nps);
                            }
                        });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Schedule").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.schedule, Schedule::Sch40, "Sch 40");
                        ui.radio_value(&mut s.schedule, Schedule::Sch80, "Sch 80");
                        ui.radio_value(&mut s.schedule, Schedule::Sch160, "Sch 160");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Material").strong());
                    ui.horizontal_wrapped(|ui| {
                        ui.radio_value(&mut s.material, Material::CarbonSteel, "Carbon steel");
                        ui.radio_value(&mut s.material, Material::StainlessSteel, "Stainless");
                        ui.radio_value(&mut s.material, Material::Copper, "Copper");
                        ui.radio_value(&mut s.material, Material::Pvc, "PVC");
                        ui.radio_value(&mut s.material, Material::Pex, "PEX");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Geometry").strong());
                    ui.horizontal(|ui| {
                        ui.label("length (mm)");
                        ui.add(egui::DragValue::new(&mut s.length_mm).speed(10.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("\u{25B6} Compute").strong())
                        .clicked()
                    {
                        run_piping(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Section").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });
}

/// Build a [`PipeSection`] from the form, run the section calculations,
/// and format the readout. Extracted from the draw closure so it is
/// unit-testable. Surfaces any backend [`PipingError`] (unknown NPS, or a
/// schedule with no tabulated wall thickness) as the error string.
fn run_piping(s: &mut PipingWorkbenchState) {
    s.error = None;

    if !(s.length_mm.is_finite() && s.length_mm > 0.0) {
        s.error = Some("length must be positive".into());
        return;
    }

    // A horizontal section from the origin; only its length feeds the
    // external-surface calculation, so the orientation is irrelevant.
    let section = PipeSection::new(
        Vector3::zeros(),
        Vector3::new(s.length_mm, 0.0, 0.0),
        s.nominal_size.clone(),
        s.schedule,
        s.material,
    );

    let compute = || -> Result<(f64, f64, f64, f64, f64, f64), PipingError> {
        Ok((
            section.outer_diameter_mm()?,
            section.inner_diameter_mm()?,
            section.flow_area_mm2()?,
            section.metal_cross_section_mm2()?,
            section.wetted_perimeter_mm()?,
            section.external_surface_area_mm2()?,
        ))
    };
    let (od, id, flow, metal, wetted, external) = match compute() {
        Ok(t) => t,
        Err(e) => {
            s.error = Some(e.to_string());
            return;
        }
    };

    s.result = format!(
        "NPS            : {}  ({:?})\n\
         material       : {:?}\n\
         length         : {:.1} mm\n\n\
         outer diameter : {:.3} mm\n\
         inner diameter : {:.3} mm   (OD \u{2212} 2\u{00B7}wall)\n\
         flow area      : {:.2} mm\u{00B2}   (\u{03C0}\u{00B7}ID\u{00B2}/4)\n\
         metal section  : {:.2} mm\u{00B2}   (\u{03C0}\u{00B7}(OD\u{00B2}\u{2212}ID\u{00B2})/4)\n\
         wetted perim.  : {:.3} mm   (\u{03C0}\u{00B7}ID)\n\
         external surf. : {:.0} mm\u{00B2}   (\u{03C0}\u{00B7}OD\u{00B7}L)",
        s.nominal_size,
        s.schedule,
        s.material,
        s.length_mm,
        od,
        id,
        flow,
        metal,
        wetted,
        external,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = PipingWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn compute_default_nps2_sch40() {
        let mut s = PipingWorkbenchState::default();
        run_piping(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        assert!(s.result.contains("outer diameter"));
        assert!(s.result.contains("inner diameter"));
        assert!(s.result.contains("flow area"));
        assert!(s.result.contains("wetted perim"));
        // NPS 2 outer diameter = 2.375 in × 25.4 = 60.325 mm; the bore is
        // strictly smaller. Recompute via the backend to confirm the
        // default is the intended worked example.
        let section = PipeSection::new(
            Vector3::zeros(),
            Vector3::new(1000.0, 0.0, 0.0),
            "2",
            Schedule::Sch40,
            Material::CarbonSteel,
        );
        assert!((section.outer_diameter_mm().unwrap() - 60.325).abs() < 1e-6);
        assert!(section.inner_diameter_mm().unwrap() < section.outer_diameter_mm().unwrap());
    }

    #[test]
    fn compute_rejects_bad_length_and_unknown_nps() {
        // Non-positive length → validation error before any backend call.
        let mut s = PipingWorkbenchState {
            length_mm: 0.0,
            ..Default::default()
        };
        run_piping(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());
        // Unknown NPS → the backend's UnknownNps error is surfaced.
        let mut s2 = PipingWorkbenchState {
            nominal_size: "999".to_string(),
            ..Default::default()
        };
        run_piping(&mut s2);
        assert!(s2.error.is_some());
        assert!(s2.result.is_empty());
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
            draw_piping_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_piping_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_piping_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_piping_workbench = true;
        run_piping(&mut app.piping);
        app.piping.error = Some("invalid pipe parameters".to_string());
        draw_workbench(&mut app);
    }
}
