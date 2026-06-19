//! The right-side **Battery ECM Workbench** panel — native first-order
//! Thevenin equivalent-circuit cell analysis over `valenx-battery-ecm`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_batteryecm_workbench`,
//! toggled from the View menu. The form sets a cell's static parameters
//! (series resistance `R0`, polarisation pair `R1`/`C1`, capacity), a state
//! of charge and a load current; "Analyze" reads the open-circuit voltage
//! from the OCV-SoC table and reports the terminal voltage under load
//! `V = OCV(SoC) − I·R0 − V_rc`, and "Show 3-D" loads a representative
//! cylindrical cell solid into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_battery_ecm::ocv::OcvSocTable;
use valenx_battery_ecm::thevenin::{CellParams, CellState};

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// How the cell is loaded for the analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LoadMode {
    /// Conventional current leaves the positive terminal (`I > 0`); the
    /// terminal voltage drops below OCV.
    Discharge,
    /// No external load (`I = 0`); the terminal voltage equals OCV when
    /// the cell is relaxed.
    Rest,
    /// Current is forced into the cell (`I < 0`); the terminal voltage
    /// rises above OCV.
    Charge,
}

/// Persistent form + result state for the Battery ECM Workbench.
pub struct BatteryEcmWorkbenchState {
    /// Series ohmic resistance `R0` (ohms).
    r0_ohm: f64,
    /// Polarisation resistance `R1` (ohms) of the RC pair.
    r1_ohm: f64,
    /// Polarisation capacitance `C1` (farads) of the RC pair.
    c1_farad: f64,
    /// Usable cell capacity (ampere-hours).
    capacity_ah: f64,
    /// State of charge to evaluate at (fraction in `[0, 1]`).
    soc: f64,
    /// Load-current magnitude (amperes); the sign is set by the mode.
    current_a: f64,
    /// Whether the cell is discharging, resting or charging.
    mode: LoadMode,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D cell solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for BatteryEcmWorkbenchState {
    fn default() -> Self {
        // A representative ~2 Ah cylindrical Li-ion-like cell: 10 mOhm
        // ohmic, a 15 mOhm / 2000 F polarisation pair, on a smooth
        // 3.0-4.2 V rest-voltage curve, at half charge under a 2 A
        // discharge. OCV(0.5) = 3.70 V, so the steady terminal voltage
        // is 3.70 - 2*(0.010 + 0.015) = 3.65 V.
        Self {
            r0_ohm: 0.010,
            r1_ohm: 0.015,
            c1_farad: 2000.0,
            capacity_ah: 2.0,
            soc: 0.5,
            current_a: 2.0,
            mode: LoadMode::Discharge,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Battery ECM Workbench right-side panel. A no-op when the
/// `show_batteryecm_workbench` toggle is off.
pub fn draw_batteryecm_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_batteryecm_workbench {
        return;
    }

    egui::SidePanel::right("valenx_batteryecm_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Battery ECM",
                "native first-order Thevenin cell voltage · valenx-battery-ecm",
            ) {
                app.show_batteryecm_workbench = false;
            }

            let s = &mut app.batteryecm;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Cell parameters").strong());
                    ui.horizontal(|ui| {
                        ui.label("R0 ohmic (Ω)");
                        ui.add(egui::DragValue::new(&mut s.r0_ohm).speed(0.001));
                    });
                    ui.horizontal(|ui| {
                        ui.label("R1 polar. (Ω)");
                        ui.add(egui::DragValue::new(&mut s.r1_ohm).speed(0.001));
                    });
                    ui.horizontal(|ui| {
                        ui.label("C1 polar. (F)");
                        ui.add(egui::DragValue::new(&mut s.c1_farad).speed(50.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("capacity (Ah)");
                        ui.add(egui::DragValue::new(&mut s.capacity_ah).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Operating point").strong());
                    ui.horizontal(|ui| {
                        ui.label("state of charge");
                        ui.add(egui::DragValue::new(&mut s.soc).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("current |I| (A)");
                        ui.add(egui::DragValue::new(&mut s.current_a).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Load").strong());
                    ui.radio_value(&mut s.mode, LoadMode::Discharge, "discharge (I > 0)");
                    ui.radio_value(&mut s.mode, LoadMode::Rest, "rest (I = 0)");
                    ui.radio_value(&mut s.mode, LoadMode::Charge, "charge (I < 0)");

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_batteryecm(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative cylindrical cell (body + positive terminal nub) as a 3-D solid and load it into the central viewport to orbit",
                        )
                        .clicked()
                    {
                        s.show_3d_request = true;
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Terminal voltage").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });

    // Serviced after the panel draws (the `&mut app.batteryecm` borrow is
    // released here): build the cell's 3-D solid and load it.
    if app.batteryecm.show_3d_request {
        app.batteryecm.show_3d_request = false;
        load_cell_3d(app);
    }
}

/// Validate the form, evaluate the cell and format the readout.
fn run_batteryecm(s: &mut BatteryEcmWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The signed load current (amperes) implied by the magnitude and mode:
/// positive on discharge, zero at rest, negative on charge.
fn signed_current(s: &BatteryEcmWorkbenchState) -> f64 {
    match s.mode {
        LoadMode::Discharge => s.current_a,
        LoadMode::Rest => 0.0,
        LoadMode::Charge => -s.current_a,
    }
}

/// A small but realistic monotone OCV-SoC curve (3.0 V empty to 4.2 V
/// full) for a generic lithium-ion-like cell. Extracted so it is shared
/// by the readout and unit-testable.
fn ocv_table() -> Result<OcvSocTable, String> {
    OcvSocTable::new(
        vec![0.00, 0.10, 0.30, 0.50, 0.70, 0.90, 1.00],
        vec![3.00, 3.45, 3.60, 3.70, 3.85, 4.05, 4.20],
    )
    .map_err(|e| e.to_string())
}

/// Build the validated [`CellState`] for the current form (fully relaxed,
/// `V_rc = 0`). The single source of truth both the readout and the 3-D
/// gate need. Extracted so it is unit-testable and shared.
fn cell_state(s: &BatteryEcmWorkbenchState) -> Result<CellState, String> {
    let params = CellParams::new(s.r0_ohm, s.r1_ohm, s.c1_farad, s.capacity_ah)
        .map_err(|e| e.to_string())?;
    CellState::new(params, ocv_table()?, s.soc).map_err(|e| e.to_string())
}

/// Evaluate the cell and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &BatteryEcmWorkbenchState) -> Result<String, String> {
    let cell = cell_state(s)?;
    let i = signed_current(s);

    let ocv = cell.ocv();
    // Fully relaxed (V_rc = 0), so the loaded terminal voltage is the
    // instantaneous OCV - I*R0 drop.
    let v_terminal = cell.terminal_voltage(i);
    let ir0_drop = i * s.r0_ohm;
    // Diffusion over-potential's analytic steady state under a sustained
    // current I is I*R1 (the RC pair charges to I*R1 as t -> infinity).
    let v_rc_steady = i * s.r1_ohm;
    let v_steady = ocv - ir0_drop - v_rc_steady;
    let tau = s.r1_ohm * s.c1_farad;

    Ok(format!(
        "R0 / R1 / C1    : {:.4} Ω / {:.4} Ω / {:.0} F\n\
         capacity        : {:.2} Ah\n\
         RC time const τ : {:.1} s\n\
         state of charge : {:.3}\n\
         load current I  : {:.3} A\n\n\
         OCV(SoC)        : {:.4} V\n\
         I·R0 drop       : {:.4} V\n\
         I·R1 (RC final) : {:.4} V\n\
         V terminal (t=0): {:.4} V\n\
         V terminal (∞)  : {:.4} V",
        s.r0_ohm,
        s.r1_ohm,
        s.c1_farad,
        s.capacity_ah,
        tau,
        s.soc,
        i,
        ocv,
        ir0_drop,
        v_rc_steady,
        v_terminal,
        v_steady,
    ))
}

/// Append a closed cylinder (a `segments`-sided prism) about the `z` axis
/// to the buffers: centre `c`, radius `r`, half-height `hz`. Both end caps
/// are triangle-fanned.
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    r: f64,
    hz: f64,
    segments: usize,
) {
    use std::f64::consts::TAU;

    let base = nodes.len();
    // Bottom-rim ring, then top-rim ring, then the two cap centres.
    for k in 0..segments {
        let a = TAU * (k as f64) / (segments as f64);
        let (x, y) = (r * a.cos(), r * a.sin());
        nodes.push(c + Vector3::new(x, y, -hz));
    }
    for k in 0..segments {
        let a = TAU * (k as f64) / (segments as f64);
        let (x, y) = (r * a.cos(), r * a.sin());
        nodes.push(c + Vector3::new(x, y, hz));
    }
    let bottom_centre = nodes.len();
    nodes.push(c + Vector3::new(0.0, 0.0, -hz));
    let top_centre = nodes.len();
    nodes.push(c + Vector3::new(0.0, 0.0, hz));

    for k in 0..segments {
        let k1 = (k + 1) % segments;
        let b0 = base + k;
        let b1 = base + k1;
        let t0 = base + segments + k;
        let t1 = base + segments + k1;
        // Side wall: two triangles per segment.
        tris.extend_from_slice(&[b0, b1, t1, b0, t1, t0]);
        // Bottom cap fan (outward) and top cap fan.
        tris.extend_from_slice(&[bottom_centre, b1, b0]);
        tris.extend_from_slice(&[top_centre, t0, t1]);
    }
}

/// Build a representative cylindrical cell as a triangle [`Mesh`] — the
/// can body plus a small positive-terminal nub, both swept about the `z`
/// axis. Representative geometry (not to scale; the voltage numbers are
/// the `valenx-battery-ecm` result). `None` for an invalid configuration.
fn cell_solid_mesh(s: &BatteryEcmWorkbenchState) -> Option<Mesh> {
    cell_state(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Main can body (a tall cylinder centred on the origin).
    push_cylinder(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.0),
        0.18,
        0.35,
        32,
    );
    // Positive-terminal nub on top (a short, narrow cylinder).
    push_cylinder(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.39),
        0.06,
        0.04,
        24,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-battery-ecm");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D cell solid and load it into the central viewport.
fn load_cell_3d(app: &mut ValenxApp) {
    let Some(mesh) = cell_solid_mesh(&app.batteryecm) else {
        app.batteryecm.error =
            Some("cell parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<cell>/valenx-battery-ecm"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = BatteryEcmWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_ocv_and_terminal_voltage() {
        let mut s = BatteryEcmWorkbenchState::default();
        run_batteryecm(&mut s);
        assert!(
            s.error.is_none(),
            "default cell should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("OCV(SoC)"));
        assert!(s.result.contains("V terminal (t=0)"));
        // OCV at SoC=0.5 is the tabulated 3.70 V.
        assert!(s.result.contains("3.7000"));
        // V terminal at t=0 = 3.70 - 2*0.010 = 3.68 V.
        assert!(s.result.contains("3.6800"));
        // V terminal at steady state = 3.70 - 2*(0.010+0.015) = 3.65 V.
        assert!(s.result.contains("3.6500"));
    }

    #[test]
    fn analyze_rejects_non_positive_resistance() {
        let mut s = BatteryEcmWorkbenchState {
            r0_ohm: 0.0,
            ..Default::default()
        };
        run_batteryecm(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn terminal_voltage_equals_ocv_at_rest_and_drops_by_ir0_on_discharge() {
        // Ground truth: at I = 0 (relaxed) terminal V == OCV(SoC), and a
        // discharge current drops it by exactly I*R0 (V_rc = 0 at t=0).
        let s_rest = BatteryEcmWorkbenchState {
            mode: LoadMode::Rest,
            ..Default::default()
        };
        let cell = cell_state(&s_rest).unwrap();
        let ocv = cell.ocv();
        // OCV(0.5) is the tabulated breakpoint, exactly 3.70 V.
        assert!((ocv - 3.70).abs() < 1e-12, "ocv = {ocv}");
        assert!((cell.terminal_voltage(signed_current(&s_rest)) - ocv).abs() < 1e-12);

        let s_dis = BatteryEcmWorkbenchState::default(); // discharge, 2 A
        let i = signed_current(&s_dis);
        assert!((i - 2.0).abs() < 1e-12);
        let v = cell_state(&s_dis).unwrap().terminal_voltage(i);
        let expected = ocv - i * s_dis.r0_ohm; // 3.70 - 2*0.010 = 3.68
        assert!((v - expected).abs() < 1e-12, "v = {v}");
        assert!((v - 3.68).abs() < 1e-12, "v = {v}");
        assert!(v < ocv);
    }

    #[test]
    fn charge_mode_raises_terminal_voltage_above_ocv() {
        let s = BatteryEcmWorkbenchState {
            mode: LoadMode::Charge,
            ..Default::default()
        };
        let cell = cell_state(&s).unwrap();
        let i = signed_current(&s); // -2 A
        assert!((i + 2.0).abs() < 1e-12);
        let v = cell.terminal_voltage(i);
        assert!(v > cell.ocv());
        // 3.70 - (-2)*0.010 = 3.72 V.
        assert!((v - 3.72).abs() < 1e-12, "v = {v}");
    }

    #[test]
    fn cell_mesh_for_default_is_nonempty_and_in_range() {
        let s = BatteryEcmWorkbenchState::default();
        let mesh = cell_solid_mesh(&s).expect("default cell yields a solid");
        assert!(mesh.nodes.len() > 8, "expected can body + terminal nub");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn cell_mesh_none_for_invalid() {
        let s = BatteryEcmWorkbenchState {
            capacity_ah: 0.0,
            ..Default::default()
        };
        assert!(cell_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_batteryecm_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_batteryecm_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_batteryecm_workbench = true;
        run_batteryecm(&mut app.batteryecm);
        draw_workbench(&mut app);
    }
}
