//! The right-side **Op-Amp Workbench** panel — native ideal closed-loop
//! op-amp analysis over `valenx-opamp`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_opamp_workbench`,
//! toggled from the View menu. The form picks an inverting or
//! non-inverting topology, sets the input / feedback resistors, an input
//! voltage and a single-pole gain-bandwidth product; "Analyze" reports the
//! closed-loop gain (`-Rf/Rin` or `1 + Rf/Rin`), the gain magnitude, the
//! output voltage and the GBW-limited closed-loop bandwidth, and
//! "Show 3-D" loads a representative 8-pin DIP op-amp IC solid into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use valenx_opamp::{Gbw, Inverting, NonInverting};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which closed-loop topology the workbench analyses.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Topology {
    /// Inverting amplifier, gain `-Rf/Rin`.
    Inverting,
    /// Non-inverting amplifier, gain `1 + Rf/Rin`.
    NonInverting,
}

/// Persistent form + result state for the Op-Amp Workbench.
pub struct OpAmpWorkbenchState {
    /// The closed-loop topology to analyse.
    topology: Topology,
    /// Input resistance `Rin` (ohms).
    r_in_ohm: f64,
    /// Feedback resistance `Rf` (ohms).
    r_f_ohm: f64,
    /// Input signal voltage `Vin` (volts).
    v_in_v: f64,
    /// Single-pole gain-bandwidth product `GBW` (hertz).
    gbw_hz: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D op-amp solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for OpAmpWorkbenchState {
    fn default() -> Self {
        // A textbook inverting stage: Rin = 1 kΩ, Rf = 10 kΩ → gain = -10,
        // with a 0.1 V input (→ -1.0 V out) and a 1 MHz GBW part, so the
        // closed-loop bandwidth is 1 MHz / 10 = 100 kHz.
        Self {
            topology: Topology::Inverting,
            r_in_ohm: 1000.0,
            r_f_ohm: 10000.0,
            v_in_v: 0.1,
            gbw_hz: 1.0e6,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Op-Amp Workbench right-side panel. A no-op when the
/// `show_opamp_workbench` toggle is off.
pub fn draw_opamp_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_opamp_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(app, ctx, "valenx_opamp_workbench", "Op-Amp", |app, ui| {
            ui.label(egui::RichText::new("native ideal closed-loop op-amp gain & bandwidth · valenx-opamp").weak().small());
            ui.separator();

            let s = &mut app.opamp;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Topology").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.topology, Topology::Inverting, "inverting");
                        ui.radio_value(
                            &mut s.topology,
                            Topology::NonInverting,
                            "non-inverting",
                        );
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Resistors").strong());
                    ui.horizontal(|ui| {
                        ui.label("input Rin (Ω)");
                        ui.add(egui::DragValue::new(&mut s.r_in_ohm).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("feedback Rf (Ω)");
                        ui.add(egui::DragValue::new(&mut s.r_f_ohm).speed(10.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Signal & bandwidth").strong());
                    ui.horizontal(|ui| {
                        ui.label("input Vin (V)");
                        ui.add(egui::DragValue::new(&mut s.v_in_v).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("GBW (Hz)");
                        ui.add(egui::DragValue::new(&mut s.gbw_hz).speed(1000.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_opamp(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative 8-pin DIP op-amp IC (body, notch tab and pin rows) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Closed-loop response").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        }, );
    if close {
        app.show_opamp_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.opamp` borrow is
    // released here): build the op-amp's 3-D solid and load it.
    if app.opamp.show_3d_request {
        app.opamp.show_3d_request = false;
        load_opamp_3d(app);
    }
}

/// Validate the form, evaluate the stage and format the readout.
fn run_opamp(s: &mut OpAmpWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The closed-loop gain (signed) and gain magnitude for the configured
/// topology — the quantities both the readout and the 3-D gate need.
/// Extracted so it is unit-testable and shared.
fn gains(s: &OpAmpWorkbenchState) -> Result<(f64, f64), String> {
    match s.topology {
        Topology::Inverting => {
            let amp = Inverting::new(s.r_in_ohm, s.r_f_ohm).map_err(|e| e.to_string())?;
            Ok((amp.gain(), amp.gain_magnitude()))
        }
        Topology::NonInverting => {
            let amp = NonInverting::new(s.r_in_ohm, s.r_f_ohm).map_err(|e| e.to_string())?;
            Ok((amp.gain(), amp.gain_magnitude()))
        }
    }
}

/// Evaluate the stage and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &OpAmpWorkbenchState) -> Result<String, String> {
    let (gain, gain_mag) = gains(s)?;
    let v_out = gain * s.v_in_v;
    let gbw = Gbw::new(s.gbw_hz).map_err(|e| e.to_string())?;
    let bw = gbw
        .closed_loop_bandwidth(gain_mag)
        .map_err(|e| e.to_string())?;

    let topo = match s.topology {
        Topology::Inverting => "inverting (G = -Rf/Rin)",
        Topology::NonInverting => "non-inverting (G = 1 + Rf/Rin)",
    };

    Ok(format!(
        "topology        : {topo}\n\
         Rin / Rf        : {:.1} / {:.1} Ω\n\
         input Vin       : {:.4} V\n\n\
         closed-loop gain: {:.4}\n\
         gain magnitude  : {:.4}\n\
         output Vout     : {:.4} V\n\
         GBW             : {:.1} Hz\n\
         closed-loop BW  : {:.1} Hz",
        s.r_in_ohm, s.r_f_ohm, s.v_in_v, gain, gain_mag, v_out, s.gbw_hz, bw,
    ))
}

/// Append an outward-facing box (centre `c`, half-extents `h`) to the
/// buffers.
fn push_box(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    h: Vector3<f64>,
) {
    let base = nodes.len();
    let signs = [
        (-1.0, -1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    ];
    for (sx, sy, sz) in signs {
        nodes.push(c + Vector3::new(sx * h.x, sy * h.y, sz * h.z));
    }
    let faces = [
        [1, 2, 6, 5],
        [0, 4, 7, 3],
        [3, 7, 6, 2],
        [0, 1, 5, 4],
        [4, 5, 6, 7],
        [0, 3, 2, 1],
    ];
    for f in faces {
        tris.extend_from_slice(&[
            base + f[0],
            base + f[1],
            base + f[2],
            base + f[0],
            base + f[2],
            base + f[3],
        ]);
    }
}

/// Build a representative 8-pin DIP op-amp IC as a triangle [`Mesh`] — a
/// rectangular body with a small pin-1 notch tab and two rows of four
/// pins. Representative geometry (not to scale; the gain / bandwidth
/// numbers are the `valenx-opamp` result). `None` for an invalid
/// configuration.
fn opamp_solid_mesh(s: &OpAmpWorkbenchState) -> Option<Mesh> {
    gains(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Main IC body (DIP package), long in y, raised in z above the pins.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.25),
        Vector3::new(0.30, 0.60, 0.12),
    );
    // Pin-1 notch tab on one end (a small half-cylinder is approximated by
    // a thin box centred on the -y end face).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, -0.60, 0.30),
        Vector3::new(0.07, 0.03, 0.04),
    );

    // Two rows of four pins along the x edges, descending below the body.
    let pin_half = Vector3::new(0.04, 0.05, 0.14);
    let ys = [-0.39, -0.13, 0.13, 0.39];
    for &y in &ys {
        for &x in &[-0.34, 0.34] {
            push_box(&mut nodes, &mut tris, Vector3::new(x, y, 0.05), pin_half);
        }
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-opamp");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D op-amp solid and load it into the central viewport.
fn load_opamp_3d(app: &mut ValenxApp) {
    let Some(mesh) = opamp_solid_mesh(&app.opamp) else {
        app.opamp.error = Some("op-amp parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<opamp>/valenx-opamp"),
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
        let s = OpAmpWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_gain_and_bandwidth() {
        let mut s = OpAmpWorkbenchState::default();
        run_opamp(&mut s);
        assert!(
            s.error.is_none(),
            "default stage should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("closed-loop gain"));
        assert!(s.result.contains("output Vout"));
        assert!(s.result.contains("closed-loop BW"));
        // Rin = 1k, Rf = 10k → inverting gain -10, |gain| 10.
        assert!(s.result.contains("-10.0000"));
        // Vout = -10 * 0.1 = -1.0 V.
        assert!(s.result.contains("-1.0000 V"));
        // GBW 1 MHz / |gain| 10 = 100000 Hz.
        assert!(s.result.contains("100000.0 Hz"));
    }

    #[test]
    fn analyze_non_inverting_gain_is_one_plus_ratio() {
        // Ground truth: non-inverting gain = 1 + Rf/Rin = 1 + 10k/1k = 11.
        let mut s = OpAmpWorkbenchState {
            topology: Topology::NonInverting,
            ..Default::default()
        };
        run_opamp(&mut s);
        assert!(s.error.is_none(), "{:?}", s.error);
        let (gain, gain_mag) = gains(&s).unwrap();
        let expected: f64 = 1.0 + 10000.0 / 1000.0;
        assert!((gain - expected).abs() < 1e-12, "gain = {gain}");
        assert!((gain_mag - 11.0).abs() < 1e-12, "|gain| = {gain_mag}");
        assert!(s.result.contains("11.0000"));
    }

    #[test]
    fn inverting_gain_is_negative_ratio_ground_truth() {
        // Ground truth: inverting gain = -Rf/Rin = -10k/1k = -10, hand-computed.
        let s = OpAmpWorkbenchState::default();
        let (gain, gain_mag) = gains(&s).unwrap();
        let expected: f64 = -10000.0 / 1000.0;
        assert!((gain - expected).abs() < 1e-12, "gain = {gain}");
        assert!((gain - (-10.0)).abs() < 1e-12, "gain = {gain}");
        assert!((gain_mag - 10.0).abs() < 1e-12, "|gain| = {gain_mag}");
    }

    #[test]
    fn analyze_rejects_zero_resistor() {
        let mut s = OpAmpWorkbenchState {
            r_in_ohm: 0.0,
            ..Default::default()
        };
        run_opamp(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn opamp_mesh_for_default_is_nonempty_and_in_range() {
        let s = OpAmpWorkbenchState::default();
        let mesh = opamp_solid_mesh(&s).expect("default op-amp yields a solid");
        // Body + notch + 8 pins = 10 boxes × 8 verts = 80 nodes.
        assert!(mesh.nodes.len() > 8, "expected body + notch + pins");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_opamp_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_opamp_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_opamp_workbench = true;
        run_opamp(&mut app.opamp);
        draw_workbench(&mut app);
    }
}
