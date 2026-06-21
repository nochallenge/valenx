//! The right-side **Strain Gauge Workbench** panel — native Wheatstone-
//! bridge strain-gauge analysis over `valenx-straingauge`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_straingauge_workbench`,
//! toggled from the View menu. The form sets a gauge factor, an applied
//! strain (in microstrain), a bridge configuration (quarter / half / full),
//! an excitation voltage, a nominal gauge resistance and the specimen's
//! Young's modulus; "Analyze" reports the fractional resistance change
//! `ΔR/R = GF·ε`, the absolute `ΔR`, the normalised bridge output
//! `Vout/Vin = (N/4)·GF·ε` and the absolute output voltage, and the
//! uniaxial Hooke's-law stress `σ = E·ε`, and "Show 3-D" loads a
//! representative strained-specimen beam (with a gauge patch) into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_straingauge::{microstrain, stress, Bridge, BridgeConfig, Gauge};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Strain Gauge Workbench.
pub struct StrainGaugeWorkbenchState {
    /// Dimensionless gauge factor `GF = (ΔR/R)/ε`. Must be `> 0`.
    gauge_factor: f64,
    /// Applied strain in *microstrain* (`µε`, units of `1e-6`); signed
    /// positive in tension, negative in compression.
    strain_microstrain: f64,
    /// Wheatstone-bridge wiring (quarter / half / full).
    config: BridgeConfig,
    /// Bridge excitation voltage `Vin` (V). Must be `> 0`.
    excitation_voltage_v: f64,
    /// Nominal gauge resistance `R` (Ω), e.g. 120 or 350. Must be `> 0`.
    nominal_resistance_ohm: f64,
    /// Specimen Young's modulus `E` (GPa) for the Hooke's-law stress.
    youngs_modulus_gpa: f64,
    /// Formatted readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D specimen solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for StrainGaugeWorkbenchState {
    fn default() -> Self {
        // A standard constantan foil gauge (GF = 2.0) on a steel specimen
        // (E = 200 GPa) at 1000 µε tension, 350 Ω, 5 V quarter bridge:
        // ΔR/R = 2e-3, Vout/Vin = 5e-4, Vout = 2.5 mV, σ = 200 MPa.
        Self {
            gauge_factor: 2.0,
            strain_microstrain: 1000.0,
            config: BridgeConfig::Quarter,
            excitation_voltage_v: 5.0,
            nominal_resistance_ohm: 350.0,
            youngs_modulus_gpa: 200.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Strain Gauge Workbench right-side panel. A no-op when the
/// `show_straingauge_workbench` toggle is off.
pub fn draw_straingauge_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_straingauge_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_straingauge_workbench",
        "Strain Gauge",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native Wheatstone-bridge strain-gauge analysis · valenx-straingauge",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.straingauge;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Gauge").strong());
                    ui.horizontal(|ui| {
                        ui.label("gauge factor GF");
                        ui.add(egui::DragValue::new(&mut s.gauge_factor).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("nominal R (Ω)");
                        ui.add(egui::DragValue::new(&mut s.nominal_resistance_ohm).speed(5.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Loading").strong());
                    ui.horizontal(|ui| {
                        ui.label("strain (µε)");
                        ui.add(egui::DragValue::new(&mut s.strain_microstrain).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Young's modulus E (GPa)");
                        ui.add(egui::DragValue::new(&mut s.youngs_modulus_gpa).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Bridge").strong());
                    ui.horizontal(|ui| {
                        ui.label("excitation Vin (V)");
                        ui.add(egui::DragValue::new(&mut s.excitation_voltage_v).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.config, BridgeConfig::Quarter, "Quarter");
                        ui.radio_value(&mut s.config, BridgeConfig::Half, "Half");
                        ui.radio_value(&mut s.config, BridgeConfig::Full, "Full");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_straingauge(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative strained-specimen beam (with a bonded gauge patch) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Bridge output").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_straingauge_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.straingauge` borrow is
    // released here): build the specimen's 3-D solid and load it.
    if app.straingauge.show_3d_request {
        app.straingauge.show_3d_request = false;
        load_specimen_3d(app);
    }
}

/// Validate the form, evaluate the gauge / bridge and format the readout.
fn run_straingauge(s: &mut StrainGaugeWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Bridge`] for this configuration — the gauge plus
/// wiring both the readout and the 3-D gate need. Extracted so it is
/// unit-testable and shared.
fn bridge(s: &StrainGaugeWorkbenchState) -> Result<Bridge, String> {
    let gauge = Gauge::new(s.gauge_factor).map_err(|e| e.to_string())?;
    Ok(Bridge::new(gauge, s.config))
}

/// Evaluate the gauge / bridge / stress chain and format the full readout,
/// mapping any domain error to a display string. Extracted so it is
/// unit-testable.
fn compute(s: &StrainGaugeWorkbenchState) -> Result<String, String> {
    let b = bridge(s)?;
    // Microstrain (µε) -> pure ratio (m/m); GPa -> Pa for E.
    let strain = microstrain(s.strain_microstrain);
    let e_pa = s.youngs_modulus_gpa * 1.0e9;

    let frac = b
        .gauge
        .fractional_resistance_change(strain)
        .map_err(|e| e.to_string())?;
    let delta_r = b
        .gauge
        .resistance_change_ohm(s.nominal_resistance_ohm, strain)
        .map_err(|e| e.to_string())?;
    let ratio = b.output_ratio(strain).map_err(|e| e.to_string())?;
    let vout = b
        .output_voltage(s.excitation_voltage_v, strain)
        .map_err(|e| e.to_string())?;
    let sigma_pa = stress(e_pa, strain).map_err(|e| e.to_string())?;
    let sigma_mpa = sigma_pa / 1.0e6;
    let vout_mv = vout * 1.0e3;

    Ok(format!(
        "gauge factor GF : {:.3}\n\
         applied strain  : {:.1} µε ({:.3e} m/m)\n\
         nominal R       : {:.1} Ω\n\
         config          : {} (N = {})\n\
         excitation Vin  : {:.2} V\n\n\
         ΔR/R = GF·ε     : {:.3e}\n\
         ΔR              : {:.4} Ω\n\
         Vout/Vin        : {:.3e}\n\
         Vout            : {:.3} mV\n\
         stress σ = E·ε  : {:.2} MPa",
        s.gauge_factor,
        s.strain_microstrain,
        strain,
        s.nominal_resistance_ohm,
        s.config.label(),
        s.config.active_arms(),
        s.excitation_voltage_v,
        frac,
        delta_r,
        ratio,
        vout_mv,
        sigma_mpa,
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

/// Build the strained specimen as a triangle [`Mesh`] — a slender beam
/// (long in `x`) with a small bonded gauge patch on its top face.
/// Representative geometry (not to scale; the bridge / stress numbers are
/// the `valenx-straingauge` result). `None` for an invalid configuration.
fn specimen_solid_mesh(s: &StrainGaugeWorkbenchState) -> Option<Mesh> {
    bridge(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Specimen beam (long in x, thin in z).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.18, 0.06),
    );
    // Bonded gauge patch on the top (+z) face, near mid-span.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.075),
        Vector3::new(0.12, 0.08, 0.015),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-straingauge");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D specimen solid and load it into the central viewport.
fn load_specimen_3d(app: &mut ValenxApp) {
    let Some(mesh) = specimen_solid_mesh(&app.straingauge) else {
        app.straingauge.error =
            Some("gauge parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<specimen>/valenx-straingauge"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical straingauge workbench as a 3-D solid plus its
/// `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn straingauge_product() -> crate::WorkspaceProduct {
    let s = StrainGaugeWorkbenchState::default();
    let mesh = specimen_solid_mesh(&s).expect("canonical straingauge ⇒ specimen solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<straingauge>/valenx-specimen");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical straingauge ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Strain gauge (Wheatstone bridge)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = StrainGaugeWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_bridge_output_and_stress() {
        let mut s = StrainGaugeWorkbenchState::default();
        run_straingauge(&mut s);
        assert!(
            s.error.is_none(),
            "default gauge should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("ΔR/R"));
        assert!(s.result.contains("Vout/Vin"));
        assert!(s.result.contains("stress σ"));
        // GF = 2, ε = 1000 µε = 1e-3, quarter bridge:
        // ΔR/R = 2e-3, Vout/Vin = GF·ε/4 = 5e-4, σ = 200 GPa · 1e-3 = 200 MPa.
        assert!(s.result.contains("2.000e-3"));
        assert!(s.result.contains("5.000e-4"));
        assert!(s.result.contains("200.00 MPa"));
    }

    #[test]
    fn analyze_rejects_non_positive_gauge_factor() {
        let mut s = StrainGaugeWorkbenchState {
            gauge_factor: 0.0,
            ..Default::default()
        };
        run_straingauge(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn quarter_bridge_output_is_gf_strain_over_four() {
        // Ground truth: a quarter bridge puts out Vout/Vin = GF · ε / 4.
        // GF = 2, ε = 1e-3  ⇒  5e-4, hand-computed.
        let gauge = Gauge::new(2.0).unwrap();
        let q = Bridge::new(gauge, BridgeConfig::Quarter);
        let strain = microstrain(1000.0);
        let ratio = q.output_ratio(strain).unwrap();
        let want = 2.0 * strain / 4.0;
        assert!((ratio - 5.0e-4).abs() < 1e-15, "got {ratio}");
        assert!((ratio - want).abs() < 1e-15);
    }

    #[test]
    fn full_bridge_quadruples_quarter_output() {
        // N/4 gain: a full bridge (N = 4) is four times the quarter
        // bridge (N = 1) at the same strain.
        let quarter = StrainGaugeWorkbenchState {
            config: BridgeConfig::Quarter,
            ..Default::default()
        };
        let full = StrainGaugeWorkbenchState {
            config: BridgeConfig::Full,
            ..Default::default()
        };
        let bq = bridge(&quarter).unwrap();
        let bf = bridge(&full).unwrap();
        let strain = microstrain(quarter.strain_microstrain);
        let rq = bq.output_ratio(strain).unwrap();
        let rf = bf.output_ratio(strain).unwrap();
        assert!((rf - 4.0 * rq).abs() < 1e-15, "rq={rq} rf={rf}");
    }

    #[test]
    fn specimen_mesh_for_default_is_nonempty_and_in_range() {
        let s = StrainGaugeWorkbenchState::default();
        let mesh = specimen_solid_mesh(&s).expect("default specimen yields a solid");
        assert!(mesh.nodes.len() > 8, "expected beam + gauge patch");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn specimen_mesh_none_for_invalid() {
        let s = StrainGaugeWorkbenchState {
            gauge_factor: -1.0,
            ..Default::default()
        };
        assert!(specimen_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_straingauge_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_straingauge_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_straingauge_workbench = true;
        run_straingauge(&mut app.straingauge);
        draw_workbench(&mut app);
    }
}
