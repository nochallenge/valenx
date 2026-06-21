//! The right-side **RC / RLC Filter Workbench** panel — native analog
//! first-order RC and series-RLC resonant filter analysis over
//! `valenx-filter`.
//!
//! Mirrors the Heat Transfer / Combustion / Antenna workbenches: a
//! resizable [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_filter_workbench`, toggled from the View menu.
//! The form picks a filter type — an RC low-pass, an RC high-pass, or a
//! series RLC band-pass — and its component values. "Analyze" reports, for
//! the RC kinds, the `-3 dB` cutoff `fc = 1/(2 pi R C)`, the time constant
//! `tau = R C`, and the magnitude (linear and in dB) and phase at a chosen
//! probe frequency; for the RLC kind, the resonant frequency
//! `f0 = 1/(2 pi sqrt(L C))`, the quality factor `Q`, the `-3 dB`
//! bandwidth, and both the narrow-band and the exact half-power band
//! edges. "Show 3-D" loads a representative R / L / C component trio (a
//! row of cylindrical component bodies) into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_filter::rc::{RcFilter, RcKind};
use valenx_filter::rlc::RlcCircuit;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which filter topology the form is configured for, selecting both the
/// analysis path and the set of component fields that matter.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum FilterMode {
    /// First-order RC low-pass (output across the capacitor).
    RcLowPass,
    /// First-order RC high-pass (output across the resistor).
    RcHighPass,
    /// Series RLC resonant band-pass.
    Rlc,
}

/// Persistent form + result state for the RC / RLC Filter Workbench.
pub struct FilterWorkbenchState {
    /// Filter topology being analysed.
    mode: FilterMode,
    /// Resistance `R`, in ohms (used by every mode).
    resistance_ohm: f64,
    /// Capacitance `C`, in farads (used by every mode).
    capacitance_f: f64,
    /// Inductance `L`, in henries (used by the RLC mode only).
    inductance_h: f64,
    /// Probe frequency for the RC magnitude / phase readout, in hertz
    /// (used by the RC modes only).
    probe_hz: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D component solids (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for FilterWorkbenchState {
    fn default() -> Self {
        // The canonical textbook RC low-pass: R = 1 kOhm, C = 1 uF, giving
        // fc = 1/(2*pi*1e3*1e-6) ~ 159.1549 Hz, tau = 1 ms. Probed at the
        // cutoff, where the gain is exactly 1/sqrt(2) = -3.0103 dB and the
        // phase is -45 deg. (In RLC mode the same L / C resonate at
        // ~1591.5 Hz; with the default 1 kΩ R that is a broad, low-Q
        // (~0.1) band-pass — lower R for a sharper resonance.)
        Self {
            mode: FilterMode::RcLowPass,
            resistance_ohm: 1_000.0,
            capacitance_f: 1e-6,
            inductance_h: 10e-3,
            probe_hz: 159.154_943_091_895_34,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the RC / RLC Filter Workbench right-side panel. A no-op when the
/// `show_filter_workbench` toggle is off.
pub fn draw_filter_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_filter_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_filter_workbench",
        "RC / RLC Filter",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native RC / series-RLC analog filter response · valenx-filter",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.filter;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Filter type").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.mode, FilterMode::RcLowPass, "RC low-pass");
                        ui.radio_value(&mut s.mode, FilterMode::RcHighPass, "RC high-pass");
                    });
                    ui.radio_value(&mut s.mode, FilterMode::Rlc, "Series RLC band-pass");

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Components").strong());
                    ui.horizontal(|ui| {
                        ui.label("resistance R (Ω)");
                        ui.add(egui::DragValue::new(&mut s.resistance_ohm).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("capacitance C (F)");
                        ui.add(
                            egui::DragValue::new(&mut s.capacitance_f)
                                .speed(1e-8)
                                .max_decimals(12),
                        );
                    });
                    if s.mode == FilterMode::Rlc {
                        ui.horizontal(|ui| {
                            ui.label("inductance L (H)");
                            ui.add(
                                egui::DragValue::new(&mut s.inductance_h)
                                    .speed(1e-4)
                                    .max_decimals(9),
                            );
                        });
                    } else {
                        ui.add_space(4.0);
                        ui.label(egui::RichText::new("Probe").strong());
                        ui.horizontal(|ui| {
                            ui.label("frequency (Hz)");
                            ui.add(egui::DragValue::new(&mut s.probe_hz).speed(1.0));
                        });
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_filter(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative R / L / C component trio (a row of cylindrical component bodies) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Response").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_filter_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.filter` borrow is
    // released here): build the component trio's 3-D solid and load it.
    if app.filter.show_3d_request {
        app.filter.show_3d_request = false;
        load_filter_3d(app);
    }
}

/// Validate the form, evaluate the filter and format the readout.
fn run_filter(s: &mut FilterWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the configured filter and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &FilterWorkbenchState) -> Result<String, String> {
    match s.mode {
        FilterMode::RcLowPass => compute_rc(s, RcKind::LowPass, "RC low-pass"),
        FilterMode::RcHighPass => compute_rc(s, RcKind::HighPass, "RC high-pass"),
        FilterMode::Rlc => compute_rlc(s),
    }
}

/// Build and format the RC (low- or high-pass) readout: cutoff, time
/// constant, and the magnitude / phase at the probe frequency.
fn compute_rc(s: &FilterWorkbenchState, kind: RcKind, label: &str) -> Result<String, String> {
    let filter =
        RcFilter::new(s.resistance_ohm, s.capacitance_f, kind).map_err(|e| e.to_string())?;
    let fc = filter.cutoff_hz();
    let tau = filter.time_constant_s();
    let resp = filter.response(s.probe_hz).map_err(|e| e.to_string())?;

    Ok(format!(
        "topology        : {label}\n\
         R / C           : {r:.1} Ω / {c:.3} µF\n\
         cutoff fc       : {fc:.4} Hz\n\
         time const τ    : {tau_ms:.4} ms\n\n\
         probe frequency : {probe:.4} Hz\n\
         magnitude |H|   : {mag:.4}\n\
         magnitude       : {db:.4} dB\n\
         phase           : {phase:.2} °",
        r = s.resistance_ohm,
        c = s.capacitance_f * 1e6,
        tau_ms = tau * 1e3,
        probe = s.probe_hz,
        mag = resp.magnitude,
        db = resp.magnitude_db(),
        phase = resp.phase_deg(),
    ))
}

/// Build and format the series-RLC readout: resonant frequency, quality
/// factor, `-3 dB` bandwidth, and both the narrow-band and the exact
/// half-power band edges.
fn compute_rlc(s: &FilterWorkbenchState) -> Result<String, String> {
    let circuit = RlcCircuit::new(s.resistance_ohm, s.inductance_h, s.capacitance_f)
        .map_err(|e| e.to_string())?;
    let f0 = circuit.resonant_hz();
    let q = circuit.quality_factor();
    let bw = circuit.bandwidth_hz();

    Ok(format!(
        "topology        : series RLC band-pass\n\
         R / L / C       : {r:.1} Ω / {l:.3} mH / {c:.3} µF\n\
         resonant f0     : {f0:.4} Hz\n\
         quality factor Q: {q:.4}\n\
         bandwidth BW    : {bw:.4} Hz\n\n\
         edges (approx, f0 ± BW/2)\n\
         lower / upper   : {lo:.4} / {hi:.4} Hz\n\n\
         edges (exact, half-power)\n\
         lower / upper   : {lo_x:.4} / {hi_x:.4} Hz",
        r = s.resistance_ohm,
        l = s.inductance_h * 1e3,
        c = s.capacitance_f * 1e6,
        lo = circuit.lower_cutoff_hz(),
        hi = circuit.upper_cutoff_hz(),
        lo_x = circuit.lower_cutoff_exact_hz(),
        hi_x = circuit.upper_cutoff_exact_hz(),
    ))
}

/// Append a (double-sided) cylinder whose axis runs along `+x`, spanning
/// `base.x ..= base.x + length` with circle centre `(base.y, base.z)`.
/// Used for each R / L / C component body in the 3-D trio.
fn push_cyl_x(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    length: f64,
    r: f64,
    seg: usize,
) {
    let (x0, x1) = (base.x, base.x + length);
    let left = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x0, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    let right = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x1, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[
            left + j,
            right + j,
            right + jn,
            left + j,
            right + jn,
            left + jn,
        ]);
    }
}

/// Build the filter as a triangle [`Mesh`] — a representative row of
/// cylindrical component bodies: a resistor, an inductor, and a capacitor
/// along the `x` axis. The inductor body is shown only for the RLC mode
/// (the RC kinds have no inductor). Representative geometry (not to scale;
/// the response numbers are the `valenx-filter` result). `None` for an
/// invalid configuration.
fn filter_solid_mesh(s: &FilterWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a configuration the analysis itself accepts.
    let valid = match s.mode {
        FilterMode::RcLowPass | FilterMode::RcHighPass => {
            RcFilter::low_pass(s.resistance_ohm, s.capacitance_f).is_ok()
        }
        FilterMode::Rlc => {
            RlcCircuit::new(s.resistance_ohm, s.inductance_h, s.capacitance_f).is_ok()
        }
    };
    if !valid {
        return None;
    }

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Resistor body (a stout cylinder near the left).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-1.2, 0.0, 0.0),
        0.7,
        0.18,
        48,
    );
    // Capacitor body (a wider, shorter cylinder on the right).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(0.5, 0.0, 0.0),
        0.5,
        0.26,
        48,
    );
    // Inductor body (a long, slender coil-former) — RLC mode only.
    if s.mode == FilterMode::Rlc {
        push_cyl_x(
            &mut nodes,
            &mut tris,
            Vector3::new(-0.4, 0.0, 0.0),
            0.8,
            0.12,
            48,
        );
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-filter");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D component trio and load it into the central viewport.
fn load_filter_3d(app: &mut ValenxApp) {
    let Some(mesh) = filter_solid_mesh(&app.filter) else {
        app.filter.error =
            Some("filter parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<filter>/valenx-filter"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"filter"}`** product: the canonical
/// passive filter built as a 3-D solid, paired with the workbench's own
/// `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`FilterWorkbenchState::default`].
pub(crate) fn filter_product() -> crate::WorkspaceProduct {
    let s = FilterWorkbenchState::default();
    let mesh = filter_solid_mesh(&s).expect("canonical filter ⇒ solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<filter>/valenx-filter");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical filter ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Filter (passive)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = FilterWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_rc_reports_cutoff_and_minus_3db() {
        let mut s = FilterWorkbenchState::default();
        run_filter(&mut s);
        assert!(
            s.error.is_none(),
            "default RC low-pass should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("cutoff fc"));
        assert!(s.result.contains("time const"));
        // R = 1 kOhm, C = 1 uF -> fc ~ 159.1549 Hz, tau = 1 ms.
        assert!(s.result.contains("159.1549"));
        assert!(s.result.contains("1.0000 ms"));
        // Probed at the cutoff: gain is exactly 1/sqrt(2) = -3.0103 dB at
        // a -45 deg phase.
        assert!(s.result.contains("-3.0103 dB"));
        assert!(s.result.contains("-45.00 °"));
    }

    #[test]
    fn analyze_rlc_reports_resonance_q_and_bandwidth() {
        let mut s = FilterWorkbenchState {
            mode: FilterMode::Rlc,
            // R = 10 Ω (not the 1 kΩ RC default) gives a Q = 10 resonator.
            resistance_ohm: 10.0,
            ..Default::default()
        };
        run_filter(&mut s);
        assert!(
            s.error.is_none(),
            "default RLC should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("resonant f0"));
        assert!(s.result.contains("quality factor Q"));
        // R = 10 Ohm, L = 10 mH, C = 1 uF -> f0 ~ 1591.5494 Hz, Q = 10,
        // BW ~ 159.1549 Hz.
        assert!(s.result.contains("1591.5494"));
        assert!(s.result.contains("10.0000"));
        assert!(s.result.contains("159.1549"));
        // The exact half-power edges differ from the narrow-band ones.
        assert!(s.result.contains("1513.9602"));
        assert!(s.result.contains("1673.1151"));
    }

    #[test]
    fn analyze_rejects_zero_capacitance() {
        let mut s = FilterWorkbenchState {
            capacitance_f: 0.0,
            ..Default::default()
        };
        run_filter(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());
    }

    #[test]
    fn ground_truth_rc_cutoff_and_half_power_db() {
        // Hand-computed RC corner frequency and half-power level:
        //   fc = 1/(2*pi*R*C) with R = 1 kOhm, C = 1 uF
        //   |H(fc)| = 1/sqrt(2), so 20*log10|H| = -3.0103 dB.
        let r = 1_000.0_f64;
        let c = 1e-6_f64;
        let filter = RcFilter::low_pass(r, c).unwrap();
        let fc_expected = 1.0 / (2.0 * std::f64::consts::PI * r * c);
        assert!((filter.cutoff_hz() - fc_expected).abs() < 1e-9);
        assert!((filter.cutoff_hz() - 159.154_943_091_895_34).abs() < 1e-6);
        let db = filter.response(filter.cutoff_hz()).unwrap().magnitude_db();
        assert!((db - (-3.010_299_956_639_812_f64)).abs() < 1e-9);
    }

    #[test]
    fn filter_mesh_for_default_is_nonempty_and_in_range() {
        let s = FilterWorkbenchState::default();
        let mesh = filter_solid_mesh(&s).expect("default RC yields a solid");
        assert!(mesh.nodes.len() > 8, "expected component bodies");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
        // The RLC mode adds the inductor body, so it has strictly more nodes.
        let rlc = FilterWorkbenchState {
            mode: FilterMode::Rlc,
            ..Default::default()
        };
        let rlc_mesh = filter_solid_mesh(&rlc).expect("default RLC yields a solid");
        assert!(rlc_mesh.nodes.len() > mesh.nodes.len());
    }

    #[test]
    fn filter_mesh_none_for_invalid() {
        let s = FilterWorkbenchState {
            resistance_ohm: 0.0,
            ..Default::default()
        };
        assert!(filter_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_filter_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_filter_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_filter_workbench = true;
        run_filter(&mut app.filter);
        draw_workbench(&mut app);
    }
}
