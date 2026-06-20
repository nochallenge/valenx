//! The right-side **Resistor Network Workbench** panel — native
//! closed-form DC resistor-network analysis over `valenx-resistor-network`.
//!
//! Mirrors the Fatigue / Heat Transfer / Antenna workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_resistornetwork_workbench`,
//! toggled from the View menu. The form sets three resistors and a source,
//! and a network mode (series equivalent, parallel equivalent, or a
//! two-resistor voltage divider). "Analyze" reports the equivalent
//! resistance (series / parallel) or the divider output voltage and the
//! current split, and "Show 3-D" loads a representative row of axial
//! resistor solids into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_resistor_network::{
    current_divider_i1, current_divider_i2, parallel, series, voltage_divider,
};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which closed-form network relation the workbench evaluates.
#[derive(Clone, Copy, PartialEq, Eq)]
enum NetworkMode {
    /// Series equivalent resistance `R_eq = R1 + R2 + R3`.
    Series,
    /// Parallel equivalent resistance `1/R_eq = 1/R1 + 1/R2 + 1/R3`.
    Parallel,
    /// Two-resistor voltage divider `Vout = Vin * R2 / (R1 + R2)`.
    Divider,
}

/// Persistent form + result state for the Resistor Network Workbench.
pub struct ResistorNetworkWorkbenchState {
    /// Which network relation to evaluate.
    mode: NetworkMode,
    /// First resistor `R1` (ohm).
    r1: f64,
    /// Second resistor `R2` (ohm).
    r2: f64,
    /// Third resistor `R3` (ohm) — used by the series / parallel modes.
    r3: f64,
    /// Source voltage `Vin` (V) for the voltage-divider mode.
    vin: f64,
    /// Source current `I_in` (A) for the divider current-split readout.
    i_in: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D resistor row (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for ResistorNetworkWorkbenchState {
    fn default() -> Self {
        // Three textbook resistors 100 / 220 / 330 ohm. In series that is
        // 650 ohm; in parallel ~56.9 ohm. With Vin = 12 V the R1/R2 divider
        // gives 12 * 220 / 320 = 8.25 V.
        Self {
            mode: NetworkMode::Series,
            r1: 100.0,
            r2: 220.0,
            r3: 330.0,
            vin: 12.0,
            i_in: 1.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Resistor Network Workbench right-side panel. A no-op when the
/// `show_resistornetwork_workbench` toggle is off.
pub fn draw_resistornetwork_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_resistornetwork_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_resistornetwork_workbench",
        "Resistor Network",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native closed-form DC resistor-network analysis · valenx-resistor-network",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.resistornetwork;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Network mode").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.mode, NetworkMode::Series, "Series");
                        ui.radio_value(&mut s.mode, NetworkMode::Parallel, "Parallel");
                        ui.radio_value(&mut s.mode, NetworkMode::Divider, "Divider");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Resistors").strong());
                    ui.horizontal(|ui| {
                        ui.label("R1 (Ω)");
                        ui.add(egui::DragValue::new(&mut s.r1).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("R2 (Ω)");
                        ui.add(egui::DragValue::new(&mut s.r2).speed(5.0));
                    });
                    if s.mode != NetworkMode::Divider {
                        ui.horizontal(|ui| {
                            ui.label("R3 (Ω)");
                            ui.add(egui::DragValue::new(&mut s.r3).speed(5.0));
                        });
                    }

                    if s.mode == NetworkMode::Divider {
                        ui.add_space(4.0);
                        ui.label(egui::RichText::new("Source").strong());
                        ui.horizontal(|ui| {
                            ui.label("Vin (V)");
                            ui.add(egui::DragValue::new(&mut s.vin).speed(0.5));
                        });
                        ui.horizontal(|ui| {
                            ui.label("I_in (A)");
                            ui.add(egui::DragValue::new(&mut s.i_in).speed(0.1));
                        });
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_resistornetwork(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative row of axial resistor solids (a body cylinder with two wire leads each) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Network").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_resistornetwork_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.resistornetwork` borrow
    // is released here): build the resistor row's 3-D solid and load it.
    if app.resistornetwork.show_3d_request {
        app.resistornetwork.show_3d_request = false;
        load_network_3d(app);
    }
}

/// Validate the form, evaluate the network and format the readout.
fn run_resistornetwork(s: &mut ResistorNetworkWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the selected network relation and format the full readout,
/// mapping any domain error to a display string. Extracted so it is
/// unit-testable.
fn compute(s: &ResistorNetworkWorkbenchState) -> Result<String, String> {
    match s.mode {
        NetworkMode::Series => {
            let r_eq = series(&[s.r1, s.r2, s.r3]).map_err(|e| e.to_string())?;
            Ok(format!(
                "mode            : series\n\
                 R1 / R2 / R3    : {r1:.1} / {r2:.1} / {r3:.1} Ω\n\n\
                 R equivalent    : {r_eq:.2} Ω\n\
                 (R_eq = R1 + R2 + R3)",
                r1 = s.r1,
                r2 = s.r2,
                r3 = s.r3,
            ))
        }
        NetworkMode::Parallel => {
            let r_eq = parallel(&[s.r1, s.r2, s.r3]).map_err(|e| e.to_string())?;
            Ok(format!(
                "mode            : parallel\n\
                 R1 / R2 / R3    : {r1:.1} / {r2:.1} / {r3:.1} Ω\n\n\
                 R equivalent    : {r_eq:.2} Ω\n\
                 (1/R_eq = 1/R1 + 1/R2 + 1/R3)",
                r1 = s.r1,
                r2 = s.r2,
                r3 = s.r3,
            ))
        }
        NetworkMode::Divider => {
            let vout = voltage_divider(s.vin, s.r1, s.r2).map_err(|e| e.to_string())?;
            let i1 = current_divider_i1(s.i_in, s.r1, s.r2).map_err(|e| e.to_string())?;
            let i2 = current_divider_i2(s.i_in, s.r1, s.r2).map_err(|e| e.to_string())?;
            Ok(format!(
                "mode            : voltage divider\n\
                 R1 / R2         : {r1:.1} / {r2:.1} Ω\n\
                 Vin / I_in      : {vin:.2} V / {i_in:.2} A\n\n\
                 Vout across R2  : {vout:.3} V\n\
                 (Vout = Vin·R2/(R1+R2))\n\
                 I through R1    : {i1:.4} A\n\
                 I through R2    : {i2:.4} A",
                r1 = s.r1,
                r2 = s.r2,
                vin = s.vin,
                i_in = s.i_in,
            ))
        }
    }
}

/// Append an axis-aligned cylinder along the `x` axis (centre `c`, radius
/// `r`, half-length `hx`, `seg` facets) to the buffers as a closed
/// triangle tube with end caps.
fn push_x_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    r: f64,
    hx: f64,
    seg: usize,
) {
    let base = nodes.len();
    // Two rings of `seg` vertices at -hx and +hx.
    for ring in 0..2 {
        let x = c.x + if ring == 0 { -hx } else { hx };
        for k in 0..seg {
            let theta = TAU * (k as f64) / (seg as f64);
            nodes.push(Vector3::new(
                x,
                c.y + r * theta.cos(),
                c.z + r * theta.sin(),
            ));
        }
    }
    // Side wall: quad between ring0[k],ring0[k+1],ring1[k+1],ring1[k].
    for k in 0..seg {
        let kn = (k + 1) % seg;
        let a = base + k;
        let b = base + kn;
        let cc = base + seg + kn;
        let d = base + seg + k;
        tris.extend_from_slice(&[a, b, cc, a, cc, d]);
    }
    // End-cap centres.
    let c0 = nodes.len();
    nodes.push(Vector3::new(c.x - hx, c.y, c.z));
    let c1 = nodes.len();
    nodes.push(Vector3::new(c.x + hx, c.y, c.z));
    for k in 0..seg {
        let kn = (k + 1) % seg;
        // -x cap (faces outward in -x).
        tris.extend_from_slice(&[c0, base + kn, base + k]);
        // +x cap (faces outward in +x).
        tris.extend_from_slice(&[c1, base + seg + k, base + seg + kn]);
    }
}

/// How many resistor bodies the current mode draws (three for the
/// series / parallel reductions, two for the divider).
fn body_count(mode: NetworkMode) -> usize {
    match mode {
        NetworkMode::Series | NetworkMode::Parallel => 3,
        NetworkMode::Divider => 2,
    }
}

/// Build a representative row of axial resistor solids as a triangle
/// [`Mesh`] — one per network element, each a fat body cylinder along the
/// `x` axis with two thin wire leads, spaced along `y`. Representative
/// geometry (not to scale; the analysis numbers are the
/// `valenx-resistor-network` result). `None` for an invalid configuration.
fn network_solid_mesh(s: &ResistorNetworkWorkbenchState) -> Option<Mesh> {
    // Reject the same out-of-domain inputs the active mode's analyze would.
    compute(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let seg = 20;
    let n = body_count(s.mode);
    let pitch = 0.6;
    let y0 = -pitch * (n as f64 - 1.0) / 2.0;
    for i in 0..n {
        let y = y0 + pitch * (i as f64);
        let c = Vector3::new(0.0, y, 0.5);
        // Resistor body.
        push_x_cylinder(&mut nodes, &mut tris, c, 0.12, 0.3, seg);
        // Wire lead on the -x side.
        push_x_cylinder(
            &mut nodes,
            &mut tris,
            Vector3::new(-0.5, y, 0.5),
            0.03,
            0.2,
            seg,
        );
        // Wire lead on the +x side.
        push_x_cylinder(
            &mut nodes,
            &mut tris,
            Vector3::new(0.5, y, 0.5),
            0.03,
            0.2,
            seg,
        );
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-resistor-network");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D resistor row solid and load it into the central viewport.
fn load_network_3d(app: &mut ValenxApp) {
    let Some(mesh) = network_solid_mesh(&app.resistornetwork) else {
        app.resistornetwork.error =
            Some("network parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<network>/valenx-resistor-network"),
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
        let s = ResistorNetworkWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_series_equivalent() {
        let mut s = ResistorNetworkWorkbenchState::default();
        run_resistornetwork(&mut s);
        assert!(
            s.error.is_none(),
            "default series network should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("R equivalent"));
        // 100 + 220 + 330 = 650 ohm.
        assert!(s.result.contains("650.00"));
    }

    #[test]
    fn analyze_parallel_reports_equivalent_below_smallest() {
        let mut s = ResistorNetworkWorkbenchState {
            mode: NetworkMode::Parallel,
            ..Default::default()
        };
        run_resistornetwork(&mut s);
        assert!(s.error.is_none(), "parallel network should analyze");
        // 100 || 220 || 330 = 1 / (1/100 + 1/220 + 1/330) = 56.90 ohm.
        assert!(s.result.contains("56.90"));
    }

    #[test]
    fn analyze_divider_reports_output_voltage() {
        let mut s = ResistorNetworkWorkbenchState {
            mode: NetworkMode::Divider,
            ..Default::default()
        };
        run_resistornetwork(&mut s);
        assert!(s.error.is_none(), "divider should analyze");
        assert!(s.result.contains("Vout across R2"));
        // 12 V * 220 / (100 + 220) = 8.25 V.
        assert!(s.result.contains("8.250"));
    }

    #[test]
    fn analyze_rejects_zero_resistor() {
        let mut s = ResistorNetworkWorkbenchState {
            r1: 0.0,
            ..Default::default()
        };
        run_resistornetwork(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn parallel_of_two_equal_is_half_each() {
        // Ground truth: two equal resistors R in parallel give exactly R/2.
        let r: f64 = 1000.0;
        let r_eq = parallel(&[r, r]).expect("valid");
        assert!((r_eq - r / 2.0).abs() < 1e-9, "got {r_eq}");
        // And the hand-computed voltage divider: Vout = Vin*R2/(R1+R2).
        let vout = voltage_divider(12.0, 100.0, 220.0).expect("valid");
        let expected = 12.0 * 220.0 / (100.0 + 220.0);
        assert!((vout - expected).abs() < 1e-9, "got {vout}");
        assert!((vout - 8.25).abs() < 1e-9, "got {vout}");
    }

    #[test]
    fn network_mesh_for_default_is_nonempty_and_in_range() {
        let s = ResistorNetworkWorkbenchState::default();
        let mesh = network_solid_mesh(&s).expect("default network yields a solid");
        assert!(mesh.nodes.len() > 8, "expected a row of resistor bodies");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn network_mesh_none_for_invalid() {
        let s = ResistorNetworkWorkbenchState {
            r1: 0.0,
            ..Default::default()
        };
        assert!(network_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_resistornetwork_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_resistornetwork_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_resistornetwork_workbench = true;
        run_resistornetwork(&mut app.resistornetwork);
        draw_workbench(&mut app);
    }
}
