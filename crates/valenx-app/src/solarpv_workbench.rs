//! The right-side **Solar PV Workbench** panel — native single-diode
//! photovoltaic cell performance over `valenx-solarpv`.
//!
//! Mirrors the Wind Turbine / Rail / Drone workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_solarpv_workbench`,
//! toggled from the View menu. The form drives a
//! [`valenx_solarpv::SingleDiode`] cell; "Analyze" reports the
//! short-circuit current, open-circuit voltage, the maximum-power point
//! (Vmp / Imp / Pmax), the fill factor and the conversion efficiency, and
//! "Show 3-D panel" loads a tilted-module-on-a-post solid into the central
//! viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_solarpv::{
    celsius_to_kelvin, efficiency, fill_factor, max_power_point, SingleDiode,
    STC_IRRADIANCE_W_PER_M2,
};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Number of I-V samples used by the maximum-power-point search.
const MPP_SAMPLES: usize = 256;

/// Persistent form + result state for the Solar PV Workbench.
pub struct SolarPvWorkbenchState {
    /// Photocurrent `Iph` (A).
    photocurrent_a: f64,
    /// Diode reverse-saturation current `I0` (A).
    saturation_current_a: f64,
    /// Diode ideality factor `n`.
    ideality: f64,
    /// Cell temperature (deg C).
    temperature_c: f64,
    /// Series resistance `Rs` (ohm).
    series_ohms: f64,
    /// Shunt resistance `Rsh` (ohm).
    shunt_ohms: f64,
    /// Cell area (m^2) for the efficiency.
    cell_area_m2: f64,
    /// Plane-of-array irradiance (W/m^2) for the efficiency.
    irradiance_w_per_m2: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D panel solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for SolarPvWorkbenchState {
    fn default() -> Self {
        // A representative crystalline-silicon cell at Standard Test
        // Conditions (1000 W/m^2, 25 deg C), ~156 mm pseudo-square.
        Self {
            photocurrent_a: 9.0,
            saturation_current_a: 1.0e-9,
            ideality: 1.3,
            temperature_c: 25.0,
            series_ohms: 0.005,
            shunt_ohms: 15.0,
            cell_area_m2: 0.0244,
            irradiance_w_per_m2: STC_IRRADIANCE_W_PER_M2,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Solar PV Workbench right-side panel. A no-op when the
/// `show_solarpv_workbench` toggle is off.
pub fn draw_solarpv_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_solarpv_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_solarpv_workbench",
        "Solar PV",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native single-diode photovoltaic cell performance · valenx-solarpv",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.solarpv;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Single-diode cell").strong());
                    ui.horizontal(|ui| {
                        ui.label("photocurrent Iph (A)");
                        ui.add(egui::DragValue::new(&mut s.photocurrent_a).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("sat. current I0 (A)");
                        ui.add(egui::DragValue::new(&mut s.saturation_current_a).speed(1.0e-10));
                    });
                    ui.horizontal(|ui| {
                        ui.label("ideality n");
                        ui.add(egui::DragValue::new(&mut s.ideality).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("temperature (°C)");
                        ui.add(egui::DragValue::new(&mut s.temperature_c).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("series Rs (Ω)");
                        ui.add(egui::DragValue::new(&mut s.series_ohms).speed(0.001));
                    });
                    ui.horizontal(|ui| {
                        ui.label("shunt Rsh (Ω)");
                        ui.add(egui::DragValue::new(&mut s.shunt_ohms).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Efficiency basis").strong());
                    ui.horizontal(|ui| {
                        ui.label("cell area (m²)");
                        ui.add(egui::DragValue::new(&mut s.cell_area_m2).speed(0.001));
                    });
                    ui.horizontal(|ui| {
                        ui.label("irradiance (W/m²)");
                        ui.add(egui::DragValue::new(&mut s.irradiance_w_per_m2).speed(5.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_cell(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D panel").strong())
                        .on_hover_text(
                            "Build a tilted PV module on a support post as a 3-D solid and load it into the central viewport to orbit",
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
        },
    );
    if close {
        app.show_solarpv_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.solarpv` borrow is
    // released here): build the panel's 3-D solid and load it.
    if app.solarpv.show_3d_request {
        app.solarpv.show_3d_request = false;
        load_panel_3d(app);
    }
}

/// Build a validated [`SingleDiode`] cell from the form, mapping the domain
/// error to a display string.
fn build_cell(s: &SolarPvWorkbenchState) -> Result<SingleDiode, String> {
    SingleDiode::new(
        s.photocurrent_a,
        s.saturation_current_a,
        s.ideality,
        celsius_to_kelvin(s.temperature_c),
        s.series_ohms,
        s.shunt_ohms,
    )
    .map_err(|e| e.to_string())
}

/// Validate the form, compute the cell performance and format the readout.
fn run_cell(s: &mut SolarPvWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Compute the full performance readout, mapping any domain error to a
/// display string. Extracted so it is unit-testable.
fn compute(s: &SolarPvWorkbenchState) -> Result<String, String> {
    let cell = build_cell(s)?;
    let isc = cell.isc().map_err(|e| e.to_string())?;
    let voc = cell.voc().map_err(|e| e.to_string())?;
    let mpp = max_power_point(&cell, MPP_SAMPLES).map_err(|e| e.to_string())?;
    let ff = fill_factor(&cell, MPP_SAMPLES).map_err(|e| e.to_string())?;
    let eff = efficiency(&cell, s.irradiance_w_per_m2, s.cell_area_m2, MPP_SAMPLES)
        .map_err(|e| e.to_string())?;
    Ok(format!(
        "photocurrent Iph: {:.2} A\n\
         sat current I0  : {:.2e} A\n\
         ideality n      : {:.2}\n\
         cell temp       : {:.1} °C\n\
         Rs / Rsh        : {:.4} / {:.1} Ω\n\n\
         short-circ Isc  : {:.2} A\n\
         open-circ Voc   : {:.3} V\n\
         max-power Vmp   : {:.3} V\n\
         max-power Imp   : {:.2} A\n\
         max power Pmax  : {:.2} W\n\
         fill factor FF  : {:.3}\n\
         efficiency      : {:.1} %  (at {:.0} W/m², {:.4} m²)",
        s.photocurrent_a,
        s.saturation_current_a,
        s.ideality,
        s.temperature_c,
        s.series_ohms,
        s.shunt_ohms,
        isc,
        voc,
        mpp.v_mp,
        mpp.i_mp,
        mpp.p_max,
        ff,
        eff * 100.0,
        s.irradiance_w_per_m2,
        s.cell_area_m2,
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

/// Append a flat (double-sided) quad `a-b-c-d` (used for the panel face).
fn push_quad(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    a: Vector3<f64>,
    b: Vector3<f64>,
    c: Vector3<f64>,
    d: Vector3<f64>,
) {
    let base = nodes.len();
    nodes.push(a);
    nodes.push(b);
    nodes.push(c);
    nodes.push(d);
    tris.extend_from_slice(&[
        base,
        base + 1,
        base + 2,
        base,
        base + 2,
        base + 3,
        base,
        base + 2,
        base + 1,
        base,
        base + 3,
        base + 2,
    ]);
}

/// Build a tilted PV module on a support post as a triangle [`Mesh`] — the
/// module face is a flat panel tilted 30 deg from horizontal toward +x, on a
/// central post and a ground base. Representative geometry (the electrical
/// performance is the single-diode cell). `None` for an invalid cell / area.
fn panel_solid_mesh(s: &SolarPvWorkbenchState) -> Option<Mesh> {
    build_cell(s).ok()?;
    if !(s.cell_area_m2.is_finite() && s.cell_area_m2 > 0.0) {
        return None;
    }
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Module face: a flat rectangle tilted 30 deg from horizontal, facing +x.
    let beta = 30.0_f64.to_radians();
    let up_slope = Vector3::new(beta.cos(), 0.0, beta.sin());
    let cross = Vector3::new(0.0, 1.0, 0.0);
    let centre = Vector3::new(0.0, 0.0, 1.3);
    let (half_len, half_wid) = (1.1, 0.8);
    let c0 = centre - up_slope * half_len - cross * half_wid;
    let c1 = centre + up_slope * half_len - cross * half_wid;
    let c2 = centre + up_slope * half_len + cross * half_wid;
    let c3 = centre - up_slope * half_len + cross * half_wid;
    push_quad(&mut nodes, &mut tris, c0, c1, c2, c3);

    // Support post.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.2, 0.0, 0.6),
        Vector3::new(0.1, 0.1, 0.6),
    );
    // Ground base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.2, 0.0, 0.06),
        Vector3::new(0.45, 0.45, 0.06),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-solarpv");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D panel solid and load it into the central viewport.
fn load_panel_3d(app: &mut ValenxApp) {
    let Some(mesh) = panel_solid_mesh(&app.solarpv) else {
        app.solarpv.error =
            Some("cell / area parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<panel>/valenx-solarpv"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"solarpv"}`** product: the canonical PV
/// panel built as a 3-D solid, paired with the workbench's own `compute()`
/// readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`SolarPvWorkbenchState::default`].
pub(crate) fn solarpv_product() -> crate::WorkspaceProduct {
    let s = SolarPvWorkbenchState::default();
    let mesh = panel_solid_mesh(&s).expect("canonical PV panel ⇒ solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<panel>/valenx-solarpv");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical PV panel ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Solar PV (panel)".into(),
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
        let s = SolarPvWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_mpp_and_efficiency() {
        let mut s = SolarPvWorkbenchState::default();
        run_cell(&mut s);
        assert!(
            s.error.is_none(),
            "default cell should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("Pmax"));
        assert!(s.result.contains("efficiency"));
        assert!(s.result.contains("Voc"));
    }

    #[test]
    fn analyze_rejects_nonpositive_area() {
        let mut s = SolarPvWorkbenchState {
            cell_area_m2: 0.0,
            ..Default::default()
        };
        run_cell(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_nonpositive_irradiance() {
        let mut s = SolarPvWorkbenchState {
            irradiance_w_per_m2: 0.0,
            ..Default::default()
        };
        run_cell(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn panel_mesh_for_default_is_nonempty_and_in_range() {
        let s = SolarPvWorkbenchState::default();
        let mesh = panel_solid_mesh(&s).expect("default panel yields a solid");
        assert!(mesh.nodes.len() > 8, "expected module + post + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn panel_mesh_none_for_invalid() {
        let s = SolarPvWorkbenchState {
            cell_area_m2: 0.0,
            ..Default::default()
        };
        assert!(panel_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_solarpv_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_solarpv_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_solarpv_workbench = true;
        run_cell(&mut app.solarpv);
        draw_workbench(&mut app);
    }
}
