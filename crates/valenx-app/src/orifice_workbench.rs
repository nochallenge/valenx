//! The right-side **Orifice Meter Workbench** panel — native
//! differential-pressure flow-meter sizing over `valenx-orifice`.
//!
//! Mirrors the Pipe Flow / Heat Transfer workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_orifice_workbench`,
//! toggled from the View menu. The form sets a constriction geometry (a
//! throat bore `d` inside a pipe `D`), a meter family (orifice plate /
//! flow nozzle / Venturi tube) that fixes the discharge coefficient, and a
//! fluid density with a measured pressure drop; "Analyze" evaluates the
//! incompressible Bernoulli working equation — diameter ratio `beta`,
//! throat area, velocity-of-approach factor, flow coefficient, volumetric
//! and mass flow, and the permanent (unrecovered) pressure loss — and
//! "Show 3-D" loads a representative pipe with an orifice-plate restriction
//! into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_orifice::{Meter, MeterGeometry, MeterKind};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Orifice Meter Workbench.
pub struct OrificeWorkbenchState {
    /// Throat / bore diameter `d` (m).
    throat_diameter_m: f64,
    /// Upstream pipe internal diameter `D` (m).
    pipe_diameter_m: f64,
    /// Meter family, which fixes the representative discharge coefficient.
    kind: MeterKind,
    /// Fluid density `rho` (kg/m^3).
    density_kg_m3: f64,
    /// Measured pressure drop `dP` across the constriction (Pa).
    pressure_drop_pa: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D solid (serviced after the panel
    /// draws).
    show_3d_request: bool,
}

impl Default for OrificeWorkbenchState {
    fn default() -> Self {
        // The crate's canonical case: a 50 mm orifice bore in a 100 mm pipe
        // (beta = 0.5), water (rho = 1000), 50 kPa across the plate ->
        // Q ~ 0.01237 m^3/s (12.37 L/s), mdot ~ 12.37 kg/s.
        Self {
            throat_diameter_m: 0.05,
            pipe_diameter_m: 0.10,
            kind: MeterKind::OrificePlate,
            density_kg_m3: 1000.0,
            pressure_drop_pa: 50_000.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Orifice Meter Workbench right-side panel. A no-op when the
/// `show_orifice_workbench` toggle is off.
pub fn draw_orifice_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_orifice_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_orifice_workbench",
        "Orifice Meter",
        |app, ui| {
            ui.label(
                egui::RichText::new("native incompressible dP flow-meter sizing · valenx-orifice")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.orifice;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Geometry").strong());
                    ui.horizontal(|ui| {
                        ui.label("throat bore d (m)");
                        ui.add(egui::DragValue::new(&mut s.throat_diameter_m).speed(0.001));
                    });
                    ui.horizontal(|ui| {
                        ui.label("pipe bore D (m)");
                        ui.add(egui::DragValue::new(&mut s.pipe_diameter_m).speed(0.001));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Meter type").strong());
                    ui.radio_value(&mut s.kind, MeterKind::OrificePlate, "orifice plate (Cd≈0.61)");
                    ui.radio_value(&mut s.kind, MeterKind::FlowNozzle, "flow nozzle (Cd≈0.97)");
                    ui.radio_value(&mut s.kind, MeterKind::VenturiTube, "Venturi tube (Cd≈0.98)");

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Fluid & differential").strong());
                    ui.horizontal(|ui| {
                        ui.label("density ρ (kg/m³)");
                        ui.add(egui::DragValue::new(&mut s.density_kg_m3).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("pressure drop dP (Pa)");
                        ui.add(egui::DragValue::new(&mut s.pressure_drop_pa).speed(100.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_orifice(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative pipe run with an orifice-plate restriction (an annular bore across the bore) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Flow").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_orifice_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.orifice` borrow is
    // released here): build the meter's 3-D solid and load it.
    if app.orifice.show_3d_request {
        app.orifice.show_3d_request = false;
        load_orifice_3d(app);
    }
}

/// Validate the form, evaluate the meter and format the readout.
fn run_orifice(s: &mut OrificeWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Meter`] for the current form, the value both the
/// readout and the 3-D gate need. Extracted so it is unit-testable and
/// shared.
fn build_meter(s: &OrificeWorkbenchState) -> Result<Meter, String> {
    let geometry =
        MeterGeometry::new(s.throat_diameter_m, s.pipe_diameter_m).map_err(|e| e.to_string())?;
    Meter::with_typical_cd(geometry, s.kind).map_err(|e| e.to_string())
}

/// Evaluate the meter and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &OrificeWorkbenchState) -> Result<String, String> {
    let meter = build_meter(s)?;
    let g = meter.geometry();
    let beta = g.beta();
    let area = g.throat_area();
    let approach = g.velocity_of_approach_factor();
    let k = meter.flow_coefficient();
    let cd = meter.discharge_coefficient();

    let q = meter
        .flow_rate(s.density_kg_m3, s.pressure_drop_pa)
        .map_err(|e| e.to_string())?;
    let mdot = meter
        .mass_flow_rate(s.density_kg_m3, s.pressure_drop_pa)
        .map_err(|e| e.to_string())?;
    let loss_ratio = meter.permanent_pressure_loss_ratio();
    // The permanent loss as a fraction of the measured differential.
    let loss_pa = loss_ratio * s.pressure_drop_pa;

    let q_lps = q * 1000.0;
    let loss_frac = loss_ratio * 100.0;
    let loss_kpa = loss_pa / 1000.0;
    Ok(format!(
        "meter type      : {kind}\n\
         throat / pipe d : {d:.4} / {pipe:.4} m\n\
         density ρ       : {rho:.1} kg/m³\n\
         pressure drop   : {dp_kpa:.2} kPa\n\n\
         beta ratio d/D  : {beta:.4}\n\
         throat area A   : {area:.6} m²\n\
         approach E      : {approach:.4}\n\
         discharge Cd    : {cd:.3}\n\
         flow coeff K    : {k:.4}\n\
         flow rate Q     : {q:.6} m³/s  ({q_lps:.2} L/s)\n\
         mass flow ṁ     : {mdot:.4} kg/s\n\
         perm loss ratio : {loss_frac:.1} %\n\
         perm loss dΩ    : {loss_kpa:.2} kPa",
        kind = s.kind.code(),
        d = s.throat_diameter_m,
        pipe = s.pipe_diameter_m,
        rho = s.density_kg_m3,
        dp_kpa = s.pressure_drop_pa / 1000.0,
    ))
}

/// Append a (double-sided) cylinder whose axis runs along `+x`, spanning
/// `base.x ..= base.x + length` with circle centre `(base.y, base.z)`.
fn push_cyl_x(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    length: f64,
    r: f64,
    seg: usize,
) {
    let (x0, x1) = (base.x, base.x + length);
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x0, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x1, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[
            lo + j,
            hi + j,
            hi + jn,
            lo + j,
            hi + jn,
            lo + jn,
            lo + j,
            lo + jn,
            hi + jn,
        ]);
    }
}

/// Append the orifice plate as a flat annular disc across the pipe — an
/// outer rim of radius `r_out` (the pipe bore) with a central bore of
/// radius `r_in` (the throat), centred on the pipe axis at `x = cx` and of
/// axial half-thickness `hx`. Triangulated as a closed solid (the two
/// faces, the outer rim wall, and the inner bore wall) so the annular
/// restriction reads clearly. `seg` is the number of angular segments.
fn push_orifice_plate(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    cx: f64,
    r_out: f64,
    r_in: f64,
    hx: f64,
    seg: usize,
) {
    let base = nodes.len();
    // Order per angular step i: outer-front, outer-back, inner-front,
    // inner-back, where front is `-x` and back is `+x`.
    for i in 0..seg {
        let a = TAU * (i as f64) / (seg as f64);
        let (sin, cos) = a.sin_cos();
        nodes.push(Vector3::new(cx - hx, r_out * cos, r_out * sin));
        nodes.push(Vector3::new(cx + hx, r_out * cos, r_out * sin));
        nodes.push(Vector3::new(cx - hx, r_in * cos, r_in * sin));
        nodes.push(Vector3::new(cx + hx, r_in * cos, r_in * sin));
    }
    let mut quad = |a: usize, b: usize, c: usize, d: usize| {
        tris.extend_from_slice(&[base + a, base + b, base + c, base + a, base + c, base + d]);
    };
    for i in 0..seg {
        let j = (i + 1) % seg;
        let (of, ob, inf, ib) = (4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3);
        let (nof, nob, ninf, nib) = (4 * j, 4 * j + 1, 4 * j + 2, 4 * j + 3);
        // Front face annulus, back face annulus, outer rim, inner bore.
        quad(inf, of, nof, ninf);
        quad(ob, ib, nib, nob);
        quad(of, ob, nob, nof);
        quad(ib, inf, ninf, nib);
    }
}

/// Build the meter as a triangle [`Mesh`] — a horizontal pipe barrel with a
/// flat orifice plate (an annular restriction whose central bore is the
/// throat) set across the bore at mid-length. Representative geometry (the
/// pipe is drawn at a fixed visual size; only the plate's beta ratio is to
/// scale, and the flow numbers are the `valenx-orifice` result). `None` for
/// an invalid configuration.
fn meter_solid_mesh(s: &OrificeWorkbenchState) -> Option<Mesh> {
    let meter = build_meter(s).ok()?;
    let beta = meter.geometry().beta();

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Fixed visual pipe radius; the throat bore tracks the real beta ratio.
    let r_pipe = 0.25;
    let r_throat = beta * r_pipe;

    // Pipe barrel along +x.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-1.0, 0.0, 0.0),
        2.0,
        r_pipe,
        40,
    );
    // Orifice plate across the bore at mid-length.
    push_orifice_plate(&mut nodes, &mut tris, 0.0, r_pipe, r_throat, 0.02, 40);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-orifice");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D meter solid and load it into the central viewport.
fn load_orifice_3d(app: &mut ValenxApp) {
    let Some(mesh) = meter_solid_mesh(&app.orifice) else {
        app.orifice.error =
            Some("meter parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<meter>/valenx-orifice"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"orifice"}`** product: the canonical
/// orifice-plate flow meter built as a 3-D solid, paired with the workbench's
/// own `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`OrificeWorkbenchState::default`].
pub(crate) fn orifice_product() -> crate::WorkspaceProduct {
    let s = OrificeWorkbenchState::default();
    let mesh = meter_solid_mesh(&s).expect("canonical orifice meter ⇒ solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<orifice>/valenx-orifice");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical orifice meter ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Orifice meter (DP sizing)".into(),
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
        let s = OrificeWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_flow_and_mass_flow() {
        let mut s = OrificeWorkbenchState::default();
        run_orifice(&mut s);
        assert!(
            s.error.is_none(),
            "default meter should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("flow rate Q"));
        assert!(s.result.contains("mass flow"));
        assert!(s.result.contains("beta ratio"));
        // beta = d/D = 0.5; Cd = 0.61 for the orifice plate.
        assert!(s.result.contains("0.5000"));
        assert!(s.result.contains("0.610"));
        // The crate's canonical Q for this case: 0.012370... m^3/s
        // (12.37 L/s), mdot 12.37 kg/s.
        assert!(s.result.contains("0.012370"));
        assert!(s.result.contains("12.37"));
    }

    #[test]
    fn analyze_rejects_throat_not_smaller_than_pipe() {
        let mut s = OrificeWorkbenchState {
            throat_diameter_m: 0.10,
            pipe_diameter_m: 0.10,
            ..Default::default()
        };
        run_orifice(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn flow_matches_bernoulli_ground_truth() {
        // Ground truth: the incompressible dP working equation
        // Q = Cd * A * sqrt(2 * dP / (rho * (1 - beta^4))), hand-computed
        // for d=0.05, D=0.10, Cd=0.61, rho=1000, dP=50000.
        let cd = 0.61_f64;
        let d = 0.05_f64;
        let big_d = 0.10_f64;
        let rho = 1000.0_f64;
        let dp = 50_000.0_f64;
        let area = std::f64::consts::PI * d * d / 4.0;
        let beta = d / big_d;
        let beta4 = beta * beta * beta * beta;
        let q_expected = cd * area * (2.0 * dp / (rho * (1.0 - beta4))).sqrt();

        let meter = build_meter(&OrificeWorkbenchState::default()).expect("valid meter");
        let q = meter.flow_rate(rho, dp).expect("flow");
        assert!(
            (q - q_expected).abs() < 1e-12,
            "Q = {q}, hand calc {q_expected}"
        );
        // And the independently known magnitude.
        assert!((q - 0.012_370_124_961_719_516).abs() < 1e-12, "Q = {q}");
    }

    #[test]
    fn venturi_passes_more_flow_than_orifice_in_readout() {
        // At the same geometry / fluid / dP, the Venturi (Cd ~ 0.98) must
        // report a larger flow than the orifice plate (Cd ~ 0.61).
        let orifice = OrificeWorkbenchState {
            kind: MeterKind::OrificePlate,
            ..Default::default()
        };
        let venturi = OrificeWorkbenchState {
            kind: MeterKind::VenturiTube,
            ..Default::default()
        };
        let q_o = build_meter(&orifice)
            .unwrap()
            .flow_rate(orifice.density_kg_m3, orifice.pressure_drop_pa)
            .unwrap();
        let q_v = build_meter(&venturi)
            .unwrap()
            .flow_rate(venturi.density_kg_m3, venturi.pressure_drop_pa)
            .unwrap();
        assert!(q_v > q_o, "venturi flow {q_v} > orifice {q_o}");
    }

    #[test]
    fn meter_mesh_for_default_is_nonempty_and_in_range() {
        let s = OrificeWorkbenchState::default();
        let mesh = meter_solid_mesh(&s).expect("default meter yields a solid");
        assert!(mesh.nodes.len() > 8, "expected pipe barrel + orifice plate");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn meter_mesh_none_for_invalid() {
        let s = OrificeWorkbenchState {
            throat_diameter_m: 0.20,
            pipe_diameter_m: 0.10,
            ..Default::default()
        };
        assert!(meter_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_orifice_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_orifice_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_orifice_workbench = true;
        run_orifice(&mut app.orifice);
        draw_workbench(&mut app);
    }
}
