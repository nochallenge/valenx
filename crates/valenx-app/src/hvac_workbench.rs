//! HVAC workbench — duct sizing + pressure drop on `valenx-hvac`.
//!
//! A right-side calculation panel: pick a round or rectangular duct section,
//! set the run length and air velocity, and compute the **Darcy–Weisbach
//! pressure drop** plus the hydraulic diameter, cross-sectional area, and a
//! CFM-based recommended duct size. Pure calculation — fully headless-testable.

use eframe::egui;

use valenx_hvac::{duct::CrossSection, flow, pressure_drop};

use crate::ValenxApp;

#[derive(Clone, Copy, PartialEq, Eq)]
enum DuctShape {
    Round,
    Rect,
}

/// Persistent state for the HVAC workbench.
pub struct HvacWorkbenchState {
    shape: DuctShape,
    diameter_mm: f64,
    width_mm: f64,
    height_mm: f64,
    length_m: f64,
    velocity_ms: f64,
    friction: f64,
    cfm: f64,
    max_velocity_fpm: f64,
    results: Option<HvacResults>,
}

impl Default for HvacWorkbenchState {
    fn default() -> Self {
        Self {
            shape: DuctShape::Round,
            diameter_mm: 200.0,
            width_mm: 300.0,
            height_mm: 200.0,
            length_m: 10.0,
            velocity_ms: 5.0,
            friction: 0.02,
            cfm: 500.0,
            max_velocity_fpm: 700.0,
            results: None,
        }
    }
}

struct HvacResults {
    hydraulic_diameter_mm: f64,
    area_m2: f64,
    pressure_drop_pa: f64,
    duct_w_in: f64,
    duct_h_in: f64,
}

/// Compute the duct hydraulics for the current settings.
fn run_hvac(s: &HvacWorkbenchState) -> HvacResults {
    let cs = match s.shape {
        DuctShape::Round => CrossSection::Round {
            d: s.diameter_mm.max(10.0),
        },
        DuctShape::Rect => CrossSection::Rect {
            w: s.width_mm.max(10.0),
            h: s.height_mm.max(10.0),
        },
    };
    let dh_mm = cs.hydraulic_diameter_mm();
    let area_m2 = cs.area_mm2() * 1.0e-6;
    let pressure_drop_pa = pressure_drop::darcy_weisbach(
        dh_mm / 1000.0,
        s.length_m.max(0.1),
        s.velocity_ms.max(0.1),
        s.friction.max(0.001),
    );
    let (duct_w_in, duct_h_in) =
        flow::cfm_to_duct_size(s.cfm.max(1.0), s.max_velocity_fpm.max(50.0));
    HvacResults {
        hydraulic_diameter_mm: dh_mm,
        area_m2,
        pressure_drop_pa,
        duct_w_in,
        duct_h_in,
    }
}

/// Draw the HVAC workbench (a no-op unless toggled on via View → HVAC).
pub fn draw_hvac_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_hvac_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_hvac_workbench",
        "HVAC",
        |app, ui| {
            ui.label(
                egui::RichText::new("duct sizing + pressure drop · valenx-hvac")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.hvac;
            ui.horizontal(|ui| {
                ui.selectable_value(&mut s.shape, DuctShape::Round, "Round");
                ui.selectable_value(&mut s.shape, DuctShape::Rect, "Rectangular");
            });
            egui::Grid::new("hvac_params")
                .num_columns(2)
                .show(ui, |ui| {
                    match s.shape {
                        DuctShape::Round => {
                            ui.label("diameter (mm)");
                            ui.add(
                                egui::DragValue::new(&mut s.diameter_mm)
                                    .speed(2.0)
                                    .range(10.0..=2000.0),
                            );
                            ui.end_row();
                        }
                        DuctShape::Rect => {
                            ui.label("width (mm)");
                            ui.add(
                                egui::DragValue::new(&mut s.width_mm)
                                    .speed(2.0)
                                    .range(10.0..=2000.0),
                            );
                            ui.end_row();
                            ui.label("height (mm)");
                            ui.add(
                                egui::DragValue::new(&mut s.height_mm)
                                    .speed(2.0)
                                    .range(10.0..=2000.0),
                            );
                            ui.end_row();
                        }
                    }
                    ui.label("length (m)");
                    ui.add(
                        egui::DragValue::new(&mut s.length_m)
                            .speed(0.5)
                            .range(0.1..=500.0),
                    );
                    ui.end_row();
                    ui.label("velocity (m/s)");
                    ui.add(
                        egui::DragValue::new(&mut s.velocity_ms)
                            .speed(0.1)
                            .range(0.1..=30.0),
                    );
                    ui.end_row();
                    ui.label("friction factor");
                    ui.add(
                        egui::DragValue::new(&mut s.friction)
                            .speed(0.001)
                            .range(0.001..=0.1),
                    );
                    ui.end_row();
                    ui.label("flow (CFM)");
                    ui.add(
                        egui::DragValue::new(&mut s.cfm)
                            .speed(10.0)
                            .range(1.0..=100000.0),
                    );
                    ui.end_row();
                    ui.label("max vel (FPM)");
                    ui.add(
                        egui::DragValue::new(&mut s.max_velocity_fpm)
                            .speed(10.0)
                            .range(50.0..=5000.0),
                    );
                    ui.end_row();
                });
            ui.separator();
            if ui.button("▶ Compute").clicked() {
                let r = run_hvac(s);
                s.results = Some(r);
            }
            if let Some(r) = &s.results {
                ui.separator();
                ui.label(
                    egui::RichText::new(format!(
                        "hydraulic Ø {:.1} mm\narea {:.4} m²\npressure drop {:.2} Pa\nCFM duct size {:.1} × {:.1} in",
                        r.hydraulic_diameter_mm, r.area_m2, r.pressure_drop_pa, r.duct_w_in, r.duct_h_in,
                    ))
                    .monospace()
                    .small(),
                );
            }
        },
    );
    if close {
        app.show_hvac_workbench = false;
    }
}

/// Build the **HVAC** result card for the Workbench+Agent bridge — a DATA-ONLY
/// [`crate::WorkspaceProduct`] (`mesh: None`) whose `lines` are the genuine duct
/// hydraulics ([`run_hvac`]) for the canonical default duct (a 200 mm round duct,
/// 10 m long at 5 m/s): hydraulic diameter, area, Darcy–Weisbach pressure drop,
/// and the CFM-recommended duct size. Registered as the `"hvac"` producer in
/// [`crate::products_registry::lookup`]; the tile renders it as a text card, not
/// a 3-D view. The rows mirror the panel's own readout format.
pub(crate) fn hvac_product() -> crate::WorkspaceProduct {
    let s = HvacWorkbenchState::default();
    let r = run_hvac(&s);
    let readout = format!(
        "hydraulic Ø {:.1} mm\narea {:.4} m²\npressure drop {:.2} Pa\nCFM duct size {:.1} × {:.1} in",
        r.hydraulic_diameter_mm, r.area_m2, r.pressure_drop_pa, r.duct_w_in, r.duct_h_in,
    );
    crate::WorkspaceProduct {
        title: "HVAC".into(),
        lines: crate::products_registry::lines_from_readout(&readout),
        mesh: None,
        vertex_colors: None,
        camera: Default::default(),
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        animation: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pressure_drop_is_positive_and_grows_with_length() {
        let r1 = run_hvac(&HvacWorkbenchState::default());
        assert!(
            r1.pressure_drop_pa > 0.0,
            "ΔP positive: {}",
            r1.pressure_drop_pa
        );
        assert!(r1.hydraulic_diameter_mm > 0.0 && r1.area_m2 > 0.0);
        let r2 = run_hvac(&HvacWorkbenchState {
            length_m: 20.0,
            ..Default::default()
        });
        assert!(
            r2.pressure_drop_pa > r1.pressure_drop_pa,
            "ΔP grows with duct length"
        );
    }

    #[test]
    fn rectangular_hydraulic_diameter_matches_4a_over_p() {
        let r = run_hvac(&HvacWorkbenchState {
            shape: DuctShape::Rect,
            width_mm: 200.0,
            height_mm: 100.0,
            ..Default::default()
        });
        // 4·A/P = 4·200·100 / (2·300) = 133.33 mm.
        assert!(
            (r.hydraulic_diameter_mm - 133.333).abs() < 0.1,
            "got {}",
            r.hydraulic_diameter_mm
        );
    }
}
