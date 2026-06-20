//! The right-side **Refrigeration Workbench** panel — native
//! vapor-compression / Carnot refrigeration analysis over
//! `valenx-refrigeration`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_refrigeration_workbench`,
//! toggled from the View menu. The form sets the two reservoir temperatures
//! (in kelvin, for the Carnot limit) and the three cycle-corner specific
//! enthalpies of a single-stage R-134a-style cycle, plus a cooling duty.
//! "Analyze" builds a [`Cycle`], reports the specific refrigerating effect,
//! compressor work and condenser rejection, the actual cooling/heating COP,
//! the reversible Carnot COP and the second-law efficiency, and scales the
//! cycle to the requested duty; "Show 3-D" loads a representative compressor
//! solid into the central viewport. A radio selects whether the headline COP
//! is read as a refrigerator (cooling) or a heat pump (heating).

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_refrigeration::carnot::{
    carnot_cop_cool, carnot_cop_heat, second_law_efficiency_cool, second_law_efficiency_heat,
};
use valenx_refrigeration::cycle::Cycle;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Whether the headline coefficient of performance is read as a
/// refrigerator (the evaporator cooling effect is useful) or a heat pump
/// (the condenser heat rejection is useful).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Duty {
    /// Refrigerator / air conditioner: the useful effect is the cooling.
    Cool,
    /// Heat pump: the useful effect is the heat delivered at the condenser.
    Heat,
}

/// Persistent form + result state for the Refrigeration Workbench.
pub struct RefrigerationWorkbenchState {
    /// Cold-reservoir (refrigerated space) absolute temperature `Tc` (K).
    t_cold_k: f64,
    /// Hot-reservoir (heat-rejection) absolute temperature `Th` (K).
    t_hot_k: f64,
    /// Specific enthalpy at the evaporator outlet / compressor inlet
    /// (state 1), kJ/kg.
    h1: f64,
    /// Specific enthalpy at the compressor outlet / condenser inlet
    /// (state 2), kJ/kg.
    h2: f64,
    /// Specific enthalpy at the condenser outlet (state 3), kJ/kg. The
    /// throttle is isenthalpic so the evaporator-inlet enthalpy `h4 = h3`.
    h3: f64,
    /// Target cooling capacity (refrigeration duty) `Q_evap` (kW).
    duty_kw: f64,
    /// Which effect the headline COP reports.
    duty: Duty,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D compressor solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for RefrigerationWorkbenchState {
    fn default() -> Self {
        // A domestic R-134a refrigerator: evaporating around -20 C
        // (253.15 K), condensing around 31 C (304.15 K), with the classic
        // textbook cycle enthalpies (Cengel & Boles, ideal -20 C / 0.8 MPa
        // R-134a). Cooling COP ~3.97 against a Carnot limit ~4.96, i.e. a
        // second-law efficiency near 80 %; a 2 kW cooling duty needs about
        // 0.5 kW of compressor power.
        Self {
            t_cold_k: 253.15,
            t_hot_k: 304.15,
            h1: 239.16,
            h2: 275.39,
            h3: 95.47,
            duty_kw: 2.0,
            duty: Duty::Cool,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Refrigeration Workbench right-side panel. A no-op when the
/// `show_refrigeration_workbench` toggle is off.
pub fn draw_refrigeration_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_refrigeration_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(app, ctx, "valenx_refrigeration_workbench", "Refrigeration", |app, ui| {
            ui.label(egui::RichText::new("native vapor-compression / Carnot COP · valenx-refrigeration").weak().small());
            ui.separator();

            let s = &mut app.refrigeration;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Reservoirs (Carnot limit)").strong());
                    ui.horizontal(|ui| {
                        ui.label("cold Tc (K)");
                        ui.add(egui::DragValue::new(&mut s.t_cold_k).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("hot Th (K)");
                        ui.add(egui::DragValue::new(&mut s.t_hot_k).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Cycle enthalpies (kJ/kg)").strong());
                    ui.horizontal(|ui| {
                        ui.label("h1 evap. outlet");
                        ui.add(egui::DragValue::new(&mut s.h1).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("h2 comp. outlet");
                        ui.add(egui::DragValue::new(&mut s.h2).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("h3 cond. outlet");
                        ui.add(egui::DragValue::new(&mut s.h3).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Duty").strong());
                    ui.horizontal(|ui| {
                        ui.label("cooling load (kW)");
                        ui.add(egui::DragValue::new(&mut s.duty_kw).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.duty, Duty::Cool, "refrigerator");
                        ui.radio_value(&mut s.duty, Duty::Heat, "heat pump");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_refrigeration(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative compressor (a capped cylinder with a drive shaft) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Performance").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        }, );
    if close { app.show_refrigeration_workbench = false; }

    // Serviced after the panel draws (the `&mut app.refrigeration` borrow is
    // released here): build the compressor's 3-D solid and load it.
    if app.refrigeration.show_3d_request {
        app.refrigeration.show_3d_request = false;
        load_unit_3d(app);
    }
}

/// Validate the form, evaluate the cycle and format the readout.
fn run_refrigeration(s: &mut RefrigerationWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the refrigeration cycle and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &RefrigerationWorkbenchState) -> Result<String, String> {
    let cycle = Cycle::new(s.h1, s.h2, s.h3).map_err(|e| e.to_string())?;
    let report = cycle
        .report_for_duty(s.duty_kw)
        .map_err(|e| e.to_string())?;

    let cop_cool = report.cop_cool;
    let cop_heat = report.cop_heat;
    let carnot_cool = carnot_cop_cool(s.t_cold_k, s.t_hot_k).map_err(|e| e.to_string())?;
    let carnot_heat = carnot_cop_heat(s.t_cold_k, s.t_hot_k).map_err(|e| e.to_string())?;
    let eta_cool =
        second_law_efficiency_cool(cop_cool, s.t_cold_k, s.t_hot_k).map_err(|e| e.to_string())?;
    let eta_heat =
        second_law_efficiency_heat(cop_heat, s.t_cold_k, s.t_hot_k).map_err(|e| e.to_string())?;

    // The application-relevant figures depend on which effect is "useful".
    let (mode, cop, carnot, eta) = match s.duty {
        Duty::Cool => ("refrigerator", cop_cool, carnot_cool, eta_cool),
        Duty::Heat => ("heat pump", cop_heat, carnot_heat, eta_heat),
    };

    let lift = s.t_hot_k - s.t_cold_k;
    // `report_for_duty` always populates the flow fields; default to keep
    // the formatter total.
    let m_dot = report.mass_flow.unwrap_or(0.0);
    let w_dot = report.compressor_power.unwrap_or(0.0);
    let q_rej = report.heat_rejection_rate.unwrap_or(0.0);

    Ok(format!(
        "mode            : {mode}\n\
         reservoirs Tc/Th: {:.2} / {:.2} K  (lift {:.2} K)\n\
         enthalpies 1/2/3: {:.2} / {:.2} / {:.2} kJ/kg\n\n\
         refrig. effect  : {:.2} kJ/kg\n\
         compressor work : {:.2} kJ/kg\n\
         heat rejected   : {:.2} kJ/kg\n\n\
         COP (actual)    : {:.3}\n\
         COP (Carnot)    : {:.3}\n\
         2nd-law eff.    : {:.1} %\n\n\
         cooling duty    : {:.2} kW\n\
         mass flow       : {:.4} kg/s\n\
         compressor power: {:.3} kW\n\
         heat reject rate: {:.3} kW",
        s.t_cold_k,
        s.t_hot_k,
        lift,
        s.h1,
        s.h2,
        s.h3,
        report.refrigerating_effect,
        report.compressor_work,
        report.heat_rejected,
        cop,
        carnot,
        eta * 100.0,
        s.duty_kw,
        m_dot,
        w_dot,
        q_rej,
    ))
}

/// Append a capped cylinder along the x-axis, double-sided, centred at
/// `center` with half-length `half_len` and `radius`, `seg` segments.
fn push_cyl_x(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    center: Vector3<f64>,
    half_len: f64,
    radius: f64,
    seg: usize,
) {
    let (x0, x1) = (center.x - half_len, center.x + half_len);
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(
            x0,
            center.y + radius * a.cos(),
            center.z + radius * a.sin(),
        ));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(
            x1,
            center.y + radius * a.cos(),
            center.z + radius * a.sin(),
        ));
    }
    let cap0 = nodes.len();
    nodes.push(Vector3::new(x0, center.y, center.z));
    let cap1 = nodes.len();
    nodes.push(Vector3::new(x1, center.y, center.z));
    for j in 0..seg {
        let jn = (j + 1) % seg;
        // Side wall (double-sided).
        tris.extend_from_slice(&[
            lo + j,
            hi + j,
            hi + jn,
            lo + j,
            hi + jn,
            lo + jn,
            lo + j,
            hi + jn,
            hi + j,
            lo + j,
            lo + jn,
            hi + jn,
        ]);
        // Caps (double-sided fans).
        tris.extend_from_slice(&[cap0, lo + jn, lo + j, cap0, lo + j, lo + jn]);
        tris.extend_from_slice(&[cap1, hi + j, hi + jn, cap1, hi + jn, hi + j]);
    }
}

/// Build the compressor as a triangle [`Mesh`] — a cylindrical hermetic
/// shell with a stub drive shaft protruding from the front face.
/// Representative geometry (not to scale; the performance numbers are the
/// `valenx-refrigeration` result). `None` for an invalid cycle.
fn compressor_solid_mesh(s: &RefrigerationWorkbenchState) -> Option<Mesh> {
    Cycle::new(s.h1, s.h2, s.h3).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    // Hermetic shell / body.
    push_cyl_x(&mut nodes, &mut tris, Vector3::zeros(), 1.0, 0.9, 24);
    // Drive shaft, protruding from the +x face.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(1.4, 0.0, 0.0),
        0.5,
        0.15,
        16,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-refrigeration");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D compressor solid and load it into the central viewport.
fn load_unit_3d(app: &mut ValenxApp) {
    let Some(mesh) = compressor_solid_mesh(&app.refrigeration) else {
        app.refrigeration.error =
            Some("cycle parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<unit>/valenx-refrigeration"),
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
        let s = RefrigerationWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_cop_and_efficiency() {
        let mut s = RefrigerationWorkbenchState::default();
        run_refrigeration(&mut s);
        assert!(
            s.error.is_none(),
            "default cycle should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("COP (actual)"));
        assert!(s.result.contains("COP (Carnot)"));
        assert!(s.result.contains("2nd-law eff."));
        // The R-134a textbook cooling COP is ~3.966 against a Carnot limit
        // of ~4.964 for the 253.15 / 304.15 K reservoirs.
        assert!(s.result.contains("3.966"));
        assert!(s.result.contains("4.964"));
        // 2 kW of cooling on this cycle needs ~0.504 kW of compressor power.
        assert!(s.result.contains("0.504"));
    }

    #[test]
    fn heat_pump_mode_reports_plus_one_cop() {
        // Switching the duty to a heat pump reports the heating COP, which
        // is the cooling COP plus one (~4.966).
        let mut s = RefrigerationWorkbenchState {
            duty: Duty::Heat,
            ..Default::default()
        };
        run_refrigeration(&mut s);
        assert!(s.error.is_none(), "{:?}", s.error);
        assert!(s.result.contains("heat pump"));
        assert!(s.result.contains("4.966"));
    }

    #[test]
    fn analyze_rejects_unordered_enthalpies() {
        // h2 <= h1 leaves no compressor work; the cycle constructor rejects
        // it and the workbench surfaces the error.
        let mut s = RefrigerationWorkbenchState {
            h2: 239.16,
            ..Default::default()
        };
        run_refrigeration(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn carnot_cool_ground_truth() {
        // Ground truth: the reversible cooling COP between Tc and Th is
        // Tc / (Th - Tc). For 253.15 K and 304.15 K the lift is 51 K, so
        // COP = 253.15 / 51 = 4.96372... — verify the crate agrees with the
        // hand-computed value.
        let tc = 253.15;
        let th = 304.15;
        let carnot = carnot_cop_cool(tc, th).unwrap();
        let hand = tc / (th - tc);
        assert!((carnot - hand).abs() < 1e-12, "carnot={carnot}");
        assert!(
            (carnot - 4.963_725_490_196_081_f64).abs() < 1e-9,
            "carnot={carnot}"
        );
    }

    #[test]
    fn unit_mesh_for_default_is_nonempty_and_in_range() {
        let s = RefrigerationWorkbenchState::default();
        let mesh = compressor_solid_mesh(&s).expect("default cycle yields a solid");
        assert!(mesh.nodes.len() > 8, "expected shell + shaft");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn unit_mesh_none_for_invalid() {
        let s = RefrigerationWorkbenchState {
            h2: 239.16,
            ..Default::default()
        };
        assert!(compressor_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_refrigeration_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_refrigeration_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_refrigeration_workbench = true;
        run_refrigeration(&mut app.refrigeration);
        draw_workbench(&mut app);
    }
}
