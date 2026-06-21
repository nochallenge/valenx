//! The right-side **Thermocouple Workbench** panel — native linear-Seebeck
//! thermoelectric analysis over `valenx-thermocouple`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_thermocouple_workbench`,
//! toggled from the View menu. The form picks a named thermocouple type
//! (K / J / T / E) and sets the hot (measurement) and cold (reference)
//! junction temperatures plus a measured voltage; "Analyze" reports the raw
//! and cold-junction-compensated EMF for the forward direction and the
//! hot-junction temperature recovered from the measured voltage, and
//! "Show 3-D probe" loads a representative two-wire probe with its junction
//! bead into the central viewport.
//!
//! The model is the linear Seebeck approximation `EMF = S * (T_hot - T_cold)`
//! with a single representative sensitivity per type — research/educational
//! grade, not the full NIST ITS-90 reference functions.

use std::f64::consts::PI;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_thermocouple::{TcType, Thermocouple};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Thermocouple Workbench.
pub struct ThermocoupleWorkbenchState {
    /// Selected standard thermocouple type (sets the Seebeck sensitivity).
    tc_type: TcType,
    /// Hot (measurement) junction temperature (deg C).
    t_hot_c: f64,
    /// Cold (reference) junction temperature (deg C).
    t_cold_c: f64,
    /// A measured voltage to invert back to a temperature (millivolts).
    emf_input_mv: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D probe solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for ThermocoupleWorkbenchState {
    fn default() -> Self {
        // A type-K probe reading a 100 C process against a 25 C reference:
        // raw EMF = 41 uV/C * 75 C = 3.075 mV, compensated = 41 uV/C * 100 C
        // = 4.100 mV. The 4.100 mV inverse input recovers 125 C from a 25 C
        // reference (25 + 4.100 mV / 41 uV/C).
        Self {
            tc_type: TcType::K,
            t_hot_c: 100.0,
            t_cold_c: 25.0,
            emf_input_mv: 4.100,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Thermocouple Workbench right-side panel. A no-op when the
/// `show_thermocouple_workbench` toggle is off.
pub fn draw_thermocouple_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_thermocouple_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_thermocouple_workbench",
        "Thermocouple",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native linear-Seebeck EMF + cold-junction compensation · valenx-thermocouple",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.thermocouple;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Thermocouple type").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.tc_type, TcType::K, "K");
                        ui.radio_value(&mut s.tc_type, TcType::J, "J");
                        ui.radio_value(&mut s.tc_type, TcType::T, "T");
                        ui.radio_value(&mut s.tc_type, TcType::E, "E");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Junctions").strong());
                    ui.horizontal(|ui| {
                        ui.label("hot / measurement (°C)");
                        ui.add(egui::DragValue::new(&mut s.t_hot_c).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("cold / reference (°C)");
                        ui.add(egui::DragValue::new(&mut s.t_cold_c).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Inverse (V → T)").strong());
                    ui.horizontal(|ui| {
                        ui.label("measured EMF (mV)");
                        ui.add(egui::DragValue::new(&mut s.emf_input_mv).speed(0.05));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_thermocouple(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D probe").strong())
                        .on_hover_text(
                            "Build a representative thermocouple probe (two dissimilar-metal wires meeting at a junction bead) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Thermoelectric response").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_thermocouple_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.thermocouple` borrow is
    // released here): build the probe's 3-D solid and load it.
    if app.thermocouple.show_3d_request {
        app.thermocouple.show_3d_request = false;
        load_probe_3d(app);
    }
}

/// Validate the form, evaluate the thermocouple and format the readout.
fn run_thermocouple(s: &mut ThermocoupleWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the thermocouple and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
///
/// Reports the forward raw and cold-junction-compensated EMF for the
/// hot/cold junction pair and the hot-junction temperature recovered from
/// the separately entered measured voltage. Voltages are shown in
/// millivolts; the crate works in volts.
fn compute(s: &ThermocoupleWorkbenchState) -> Result<String, String> {
    let tc = Thermocouple::of_type(s.tc_type);
    let sens_uv_per_c = tc.sensitivity_v_per_c() * 1.0e6;

    let emf_v = tc.emf(s.t_hot_c, s.t_cold_c).map_err(|e| e.to_string())?;
    let emf_comp_v = tc
        .emf_compensated(s.t_hot_c, s.t_cold_c)
        .map_err(|e| e.to_string())?;
    let dt = s.t_hot_c - s.t_cold_c;

    // Invert the measured voltage (entered in mV) against the cold junction.
    let emf_input_v = s.emf_input_mv * 1.0e-3;
    let t_recovered = tc
        .temperature_from_emf(emf_input_v, s.t_cold_c)
        .map_err(|e| e.to_string())?;

    let type_label = tc_type_label(s.tc_type);

    Ok(format!(
        "type            : {type_label}\n\
         sensitivity S   : {sens_uv_per_c:.1} µV/°C\n\
         hot / cold      : {:.1} / {:.1} °C\n\
         delta T         : {dt:.1} °C\n\n\
         raw EMF         : {:.3} mV\n\
         compensated EMF : {:.3} mV\n\n\
         measured EMF    : {:.3} mV\n\
         recovered T_hot : {t_recovered:.2} °C",
        s.t_hot_c,
        s.t_cold_c,
        emf_v * 1.0e3,
        emf_comp_v * 1.0e3,
        s.emf_input_mv,
    ))
}

/// Short display name for a thermocouple type.
fn tc_type_label(kind: TcType) -> &'static str {
    match kind {
        TcType::K => "K (chromel/alumel)",
        TcType::J => "J (iron/constantan)",
        TcType::T => "T (copper/constantan)",
        TcType::E => "E (chromel/constantan)",
    }
}

/// Append an axial cylinder (a wire) running along `+z` from `base` for
/// `length`, with the given `radius`, approximated by a `seg`-sided prism
/// (side quads split into triangles, plus the two end caps as fans).
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    length: f64,
    radius: f64,
    seg: usize,
) {
    let ring0 = nodes.len();
    // Bottom then top ring of `seg` vertices.
    for k in 0..seg {
        let a = 2.0 * PI * (k as f64) / (seg as f64);
        let (x, y) = (radius * a.cos(), radius * a.sin());
        nodes.push(base + Vector3::new(x, y, 0.0));
    }
    let ring1 = nodes.len();
    for k in 0..seg {
        let a = 2.0 * PI * (k as f64) / (seg as f64);
        let (x, y) = (radius * a.cos(), radius * a.sin());
        nodes.push(base + Vector3::new(x, y, length));
    }
    // Side wall.
    for k in 0..seg {
        let kn = (k + 1) % seg;
        let b0 = ring0 + k;
        let b1 = ring0 + kn;
        let t0 = ring1 + k;
        let t1 = ring1 + kn;
        tris.extend_from_slice(&[b0, b1, t1, b0, t1, t0]);
    }
    // End caps via a central vertex fan.
    let cap_bot = nodes.len();
    nodes.push(base);
    let cap_top = nodes.len();
    nodes.push(base + Vector3::new(0.0, 0.0, length));
    for k in 0..seg {
        let kn = (k + 1) % seg;
        tris.extend_from_slice(&[cap_bot, ring0 + kn, ring0 + k]);
        tris.extend_from_slice(&[cap_top, ring1 + k, ring1 + kn]);
    }
}

/// Append a UV sphere of the given `radius` centred at `c` (the junction
/// bead), with `seg` longitudes and `seg` latitude bands.
fn push_sphere(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    radius: f64,
    seg: usize,
) {
    let base = nodes.len();
    // Rows i = 0..=seg are latitudes from south to north pole.
    for i in 0..=seg {
        let theta = PI * (i as f64) / (seg as f64);
        let (st, ct) = (theta.sin(), theta.cos());
        for j in 0..=seg {
            let phi = 2.0 * PI * (j as f64) / (seg as f64);
            let (sp, cp) = (phi.sin(), phi.cos());
            nodes.push(c + Vector3::new(radius * st * cp, radius * st * sp, radius * ct));
        }
    }
    let stride = seg + 1;
    for i in 0..seg {
        for j in 0..seg {
            let a = base + i * stride + j;
            let b = base + i * stride + j + 1;
            let cc = base + (i + 1) * stride + j;
            let d = base + (i + 1) * stride + j + 1;
            tris.extend_from_slice(&[a, b, d, a, d, cc]);
        }
    }
}

/// Build the thermocouple probe as a triangle [`Mesh`] — two dissimilar
/// wire cylinders running up to a junction bead where they meet, plus the
/// bead itself. Representative geometry (not to scale; the thermoelectric
/// numbers are the `valenx-thermocouple` result). `None` for an invalid
/// configuration.
fn probe_solid_mesh(s: &ThermocoupleWorkbenchState) -> Option<Mesh> {
    // Gate the geometry on the same validity the readout needs.
    let tc = Thermocouple::of_type(s.tc_type);
    tc.emf(s.t_hot_c, s.t_cold_c).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let seg = 16usize;
    let wire_r = 0.05;
    let wire_len = 1.0;
    let bead_r = 0.12;
    // The bead sits at the top of the wires; the two wires are laterally
    // offset so they visibly converge into it.
    let bead = Vector3::new(0.0, 0.0, wire_len);
    push_cylinder(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.08, 0.0, 0.0),
        wire_len,
        wire_r,
        seg,
    );
    push_cylinder(
        &mut nodes,
        &mut tris,
        Vector3::new(0.08, 0.0, 0.0),
        wire_len,
        wire_r,
        seg,
    );
    push_sphere(&mut nodes, &mut tris, bead, bead_r, seg);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-thermocouple");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D probe solid and load it into the central viewport.
fn load_probe_3d(app: &mut ValenxApp) {
    let Some(mesh) = probe_solid_mesh(&app.thermocouple) else {
        app.thermocouple.error =
            Some("thermocouple parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<probe>/valenx-thermocouple"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical thermocouple workbench as a 3-D solid
/// plus its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn thermocouple_product() -> crate::WorkspaceProduct {
    let s = ThermocoupleWorkbenchState::default();
    let mesh = probe_solid_mesh(&s).expect("canonical thermocouple ⇒ probe solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<thermocouple>/valenx-probe");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical thermocouple ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Thermocouple (Seebeck EMF)".into(),
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
        let s = ThermocoupleWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_emf_and_recovered_temperature() {
        let mut s = ThermocoupleWorkbenchState::default();
        run_thermocouple(&mut s);
        assert!(
            s.error.is_none(),
            "default probe should analyze: {:?}",
            s.error
        );
        // Type K, 100 -> 25 C: raw EMF = 41 uV/C * 75 C = 3.075 mV;
        // compensated = 41 uV/C * 100 C = 4.100 mV; the 4.100 mV inverse
        // input recovers 125.00 C against the 25 C reference.
        assert!(s.result.contains("raw EMF         : 3.075 mV"));
        assert!(s.result.contains("compensated EMF : 4.100 mV"));
        assert!(s.result.contains("recovered T_hot : 125.00 °C"));
    }

    #[test]
    fn ground_truth_equal_junctions_give_zero_emf() {
        // Ground truth: the linear Seebeck EMF vanishes when the hot and
        // cold junctions sit at the same temperature, S * (T - T) = 0.
        let mut s = ThermocoupleWorkbenchState {
            t_hot_c: 50.0,
            t_cold_c: 50.0,
            ..ThermocoupleWorkbenchState::default()
        };
        run_thermocouple(&mut s);
        assert!(s.error.is_none(), "equal junctions analyze: {:?}", s.error);
        assert!(s.result.contains("raw EMF         : 0.000 mV"));
        assert!(s.result.contains("delta T         : 0.0 °C"));
    }

    #[test]
    fn type_selection_changes_sensitivity_and_emf() {
        // Switching K -> E (61 uV/C) raises both the reported sensitivity
        // and the EMF for the same junction pair.
        let mut k = ThermocoupleWorkbenchState::default();
        run_thermocouple(&mut k);
        assert!(k.result.contains("sensitivity S   : 41.0 µV/°C"));

        let mut e = ThermocoupleWorkbenchState {
            tc_type: TcType::E,
            ..ThermocoupleWorkbenchState::default()
        };
        run_thermocouple(&mut e);
        assert!(e.result.contains("sensitivity S   : 61.0 µV/°C"));
        // E: raw EMF = 61 uV/C * 75 C = 4.575 mV.
        assert!(e.result.contains("raw EMF         : 4.575 mV"));
    }

    #[test]
    fn probe_mesh_for_default_is_nonempty_and_in_range() {
        let s = ThermocoupleWorkbenchState::default();
        let mesh = probe_solid_mesh(&s).expect("default probe yields a solid");
        assert!(mesh.nodes.len() > 8, "expected two wires + a bead");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn analyze_rejects_non_finite_temperature() {
        let mut s = ThermocoupleWorkbenchState {
            t_hot_c: f64::NAN,
            ..ThermocoupleWorkbenchState::default()
        };
        run_thermocouple(&mut s);
        assert!(s.error.is_some());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_thermocouple_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_thermocouple_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_thermocouple_workbench = true;
        run_thermocouple(&mut app.thermocouple);
        draw_workbench(&mut app);
    }
}
