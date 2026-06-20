//! The right-side **Rectifier Workbench** panel — native ideal-diode
//! rectifier figures of merit over `valenx-rectifier`.
//!
//! Mirrors the Capacitor / Heat Transfer / Antenna workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_rectifier_workbench`,
//! toggled from the View menu. The form picks a topology (half-wave or
//! full-wave), a sinusoidal peak input voltage, and a capacitor-input
//! filter (load current, mains frequency, reservoir capacitance);
//! "Analyze" reports the average (DC) and RMS output, the dimensionless
//! ripple factor / form factor / rectification efficiency, and the
//! peak-to-peak filter ripple, and "Show 3-D" loads a representative
//! diode-bridge-plus-smoothing-capacitor solid into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_rectifier::{
    capacitor_ripple_pp_for, form_factor, rectification_efficiency, ripple_factor, vdc, vrms,
    Topology,
};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Rectifier Workbench.
pub struct RectifierWorkbenchState {
    /// Which diode topology is analysed.
    topology: Topology,
    /// Sinusoidal peak input voltage `Vpeak` (volts).
    v_peak: f64,
    /// Steady load current `I` drawn from the filter (amperes).
    load_current_a: f64,
    /// Mains (line) frequency `f_mains` (hertz).
    mains_freq_hz: f64,
    /// Reservoir capacitance `C` of the input filter (farads).
    cap_farads: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D rectifier solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for RectifierWorkbenchState {
    fn default() -> Self {
        // A full-wave bridge off a 12 V peak winding feeding a 0.5 A load
        // through a 2200 uF reservoir on a 60 Hz mains (so the ripple is at
        // 120 Hz): Vdc ~ 7.64 V, Vrms ~ 8.49 V, eta ~ 81.1%, ripple ~ 1.89 V.
        Self {
            topology: Topology::FullWave,
            v_peak: 12.0,
            load_current_a: 0.5,
            mains_freq_hz: 60.0,
            cap_farads: 2200.0e-6,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Rectifier Workbench right-side panel. A no-op when the
/// `show_rectifier_workbench` toggle is off.
pub fn draw_rectifier_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_rectifier_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_rectifier_workbench",
        "Rectifier",
        |app, ui| {
            ui.label(egui::RichText::new("native ideal-diode rectifier figures + capacitor-filter ripple · valenx-rectifier").weak().small());
            ui.separator();

            let s = &mut app.rectifier;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Topology").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.topology, Topology::HalfWave, "half-wave");
                        ui.radio_value(&mut s.topology, Topology::FullWave, "full-wave");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Input").strong());
                    ui.horizontal(|ui| {
                        ui.label("peak voltage Vpeak (V)");
                        ui.add(egui::DragValue::new(&mut s.v_peak).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Capacitor filter").strong());
                    ui.horizontal(|ui| {
                        ui.label("load current I (A)");
                        ui.add(egui::DragValue::new(&mut s.load_current_a).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("mains frequency f (Hz)");
                        ui.add(egui::DragValue::new(&mut s.mains_freq_hz).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("capacitance C (F)");
                        ui.add(egui::DragValue::new(&mut s.cap_farads).speed(0.0001));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_rectifier(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative diode-bridge package with its smoothing reservoir capacitor as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Rectifier figures").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_rectifier_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.rectifier` borrow is
    // released here): build the rectifier's 3-D solid and load it.
    if app.rectifier.show_3d_request {
        app.rectifier.show_3d_request = false;
        load_rectifier_3d(app);
    }
}

/// Validate the form, evaluate the rectifier and format the readout.
fn run_rectifier(s: &mut RectifierWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the rectifier and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &RectifierWorkbenchState) -> Result<String, String> {
    let v_dc = vdc(s.topology, s.v_peak).map_err(|e| e.to_string())?;
    let v_rms = vrms(s.topology, s.v_peak).map_err(|e| e.to_string())?;
    let r = ripple_factor(s.topology).map_err(|e| e.to_string())?;
    let ff = form_factor(s.topology).map_err(|e| e.to_string())?;
    let eta = rectification_efficiency(s.topology).map_err(|e| e.to_string())?;
    let ripple_hz = s.mains_freq_hz * s.topology.ripple_frequency_multiplier();
    let v_ripple =
        capacitor_ripple_pp_for(s.topology, s.load_current_a, s.mains_freq_hz, s.cap_farads)
            .map_err(|e| e.to_string())?;

    let topology = match s.topology {
        Topology::HalfWave => "half-wave",
        Topology::FullWave => "full-wave",
    };

    Ok(format!(
        "topology        : {topology}\n\
         peak voltage    : {:.2} V\n\
         load / mains    : {:.3} A / {:.1} Hz\n\
         capacitance     : {:.1} µF\n\n\
         Vdc (average)   : {v_dc:.4} V\n\
         Vrms            : {v_rms:.4} V\n\
         ripple factor r : {r:.4}\n\
         form factor FF  : {ff:.4}\n\
         efficiency η    : {:.2} %\n\
         ripple freq     : {ripple_hz:.1} Hz\n\
         ripple Vr (p-p) : {v_ripple:.4} V",
        s.v_peak,
        s.load_current_a,
        s.mains_freq_hz,
        s.cap_farads * 1.0e6,
        eta * 100.0,
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

/// Append a closed solid cylinder of radius `r` standing along the `z`
/// axis from `z = z0` to `z = z0 + len`, centred on `(cx, cy)`, to the
/// buffers (both end caps and the side wall). `seg` is the number of
/// angular segments.
#[allow(clippy::too_many_arguments)]
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    cx: f64,
    cy: f64,
    z0: f64,
    r: f64,
    len: f64,
    seg: usize,
) {
    let base = nodes.len();
    let top = z0 + len;
    // Per angular step i: side-top node then side-bottom node.
    for i in 0..seg {
        let a = std::f64::consts::TAU * (i as f64) / (seg as f64);
        let (sin, cos) = a.sin_cos();
        nodes.push(Vector3::new(cx + r * cos, cy + r * sin, top));
        nodes.push(Vector3::new(cx + r * cos, cy + r * sin, z0));
    }
    // Two cap-centre nodes appended after the rings.
    let top_centre = nodes.len();
    nodes.push(Vector3::new(cx, cy, top));
    let bot_centre = nodes.len();
    nodes.push(Vector3::new(cx, cy, z0));
    for i in 0..seg {
        let j = (i + 1) % seg;
        let (t, b) = (2 * i, 2 * i + 1);
        let (nt, nb) = (2 * j, 2 * j + 1);
        // Side wall.
        tris.extend_from_slice(&[
            base + t,
            base + b,
            base + nb,
            base + t,
            base + nb,
            base + nt,
        ]);
        // Top cap fan.
        tris.extend_from_slice(&[base + t, base + nt, top_centre]);
        // Bottom cap fan.
        tris.extend_from_slice(&[base + b, bot_centre, base + nb]);
    }
}

/// Build the rectifier as a triangle [`Mesh`] — a diode-bridge package
/// (a box) sitting on a base next to its smoothing reservoir capacitor
/// (a standing cylinder). Representative geometry (not to scale; the
/// figures of merit are the `valenx-rectifier` result). `None` for an
/// invalid configuration (a non-positive peak voltage).
fn rectifier_solid_mesh(s: &RectifierWorkbenchState) -> Option<Mesh> {
    // Gate the geometry on a valid, analysable rectifier.
    vdc(s.topology, s.v_peak).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Diode-bridge package body.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.45, 0.0, 0.25),
        Vector3::new(0.35, 0.45, 0.2),
    );
    // Smoothing reservoir capacitor (standing cylinder).
    push_cylinder(&mut nodes, &mut tris, 0.55, 0.0, 0.1, 0.35, 0.9, 48);
    // Base board.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(1.1, 0.6, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-rectifier");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D rectifier solid and load it into the central viewport.
fn load_rectifier_3d(app: &mut ValenxApp) {
    let Some(mesh) = rectifier_solid_mesh(&app.rectifier) else {
        app.rectifier.error =
            Some("rectifier parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<rectifier>/valenx-rectifier"),
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
    use std::f64::consts::PI;

    #[test]
    fn default_state_is_idle() {
        let s = RectifierWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_outputs_and_efficiency() {
        let mut s = RectifierWorkbenchState::default();
        run_rectifier(&mut s);
        assert!(
            s.error.is_none(),
            "default rectifier should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("Vdc (average)"));
        assert!(s.result.contains("ripple factor r"));
        assert!(s.result.contains("efficiency η"));
        // Full-wave, 12 V peak: Vdc = 2*12/pi ~ 7.6394 V, Vrms = 12/sqrt2
        // ~ 8.4853 V, eta ~ 81.06%; 0.5 A / 2200 uF at 120 Hz -> ~1.8939 V.
        assert!(s.result.contains("7.6394"));
        assert!(s.result.contains("8.4853"));
        assert!(s.result.contains("81.06"));
        assert!(s.result.contains("120.0"));
        assert!(s.result.contains("1.8939"));
    }

    #[test]
    fn half_wave_reports_its_distinct_figures() {
        let mut s = RectifierWorkbenchState {
            topology: Topology::HalfWave,
            ..Default::default()
        };
        run_rectifier(&mut s);
        assert!(s.error.is_none(), "half-wave should analyze: {:?}", s.error);
        assert!(s.result.contains("half-wave"));
        // Half-wave, 12 V peak: Vdc = 12/pi ~ 3.8197 V; efficiency 40.53%;
        // ripple at the mains frequency (60 Hz, not doubled).
        assert!(s.result.contains("3.8197"));
        assert!(s.result.contains("40.53"));
        assert!(s.result.contains("60.0"));
    }

    #[test]
    fn analyze_rejects_non_positive_peak() {
        let mut s = RectifierWorkbenchState {
            v_peak: 0.0,
            ..Default::default()
        };
        run_rectifier(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn compute_matches_full_wave_vdc_ground_truth() {
        // Ground truth: the average output of an ideal full-wave rectifier
        // is Vdc = 2 * Vpeak / pi — hand-computed here and checked against
        // the crate, then confirmed to appear in the readout.
        let s = RectifierWorkbenchState::default();
        let expected: f64 = 2.0 * s.v_peak / PI;
        let api = vdc(Topology::FullWave, s.v_peak).unwrap();
        assert!(
            (api - expected).abs() < 1e-12,
            "api={api} expected={expected}"
        );
        let r = compute(&s).expect("default analyzes");
        assert!(r.contains(&format!("{expected:.4}")));
    }

    #[test]
    fn rectifier_mesh_for_default_is_nonempty_and_in_range() {
        let s = RectifierWorkbenchState::default();
        let mesh = rectifier_solid_mesh(&s).expect("default rectifier yields a solid");
        assert!(mesh.nodes.len() > 8, "expected bridge + capacitor + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn rectifier_mesh_none_for_invalid() {
        let s = RectifierWorkbenchState {
            v_peak: -1.0,
            ..Default::default()
        };
        assert!(rectifier_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_rectifier_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_rectifier_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_rectifier_workbench = true;
        run_rectifier(&mut app.rectifier);
        draw_workbench(&mut app);
    }
}
