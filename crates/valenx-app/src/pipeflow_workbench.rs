//! The right-side **Pipe Flow Workbench** panel — native Darcy-Weisbach
//! pipe-flow analysis over `valenx-pipeflow`.
//!
//! Mirrors the Heat Exchanger / Pump workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_pipeflow_workbench`,
//! toggled from the View menu. The form sets the fluid (density,
//! viscosity), the pipe (diameter, length, absolute roughness) and the
//! bulk velocity; "Analyze" runs [`valenx_pipeflow::headloss::solve_pipe`]
//! and reports the Reynolds number, flow regime, Darcy friction factor,
//! volumetric flow, head loss and pressure drop, and "Show 3-D pipe"
//! loads a flanged pipe-run solid into the central viewport.

use std::f64::consts::{PI, TAU};
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_pipeflow::headloss::solve_pipe;
use valenx_pipeflow::reynolds::FlowRegime;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Pipe Flow Workbench.
pub struct PipeFlowWorkbenchState {
    /// Fluid density `rho` (kg/m³).
    density_kg_m3: f64,
    /// Fluid dynamic viscosity `mu` (Pa·s).
    viscosity_pa_s: f64,
    /// Internal pipe diameter `D` (m).
    diameter_m: f64,
    /// Pipe run length `L` (m).
    length_m: f64,
    /// Absolute wall roughness `epsilon` (m).
    roughness_m: f64,
    /// Bulk (area-averaged) velocity `V` (m/s).
    velocity_m_s: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D pipe solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for PipeFlowWorkbenchState {
    fn default() -> Self {
        // Water at 20 °C through 100 m of 100 mm commercial-steel pipe
        // (ε ≈ 0.046 mm) at 2 m/s — Re ≈ 2.0e5, turbulent.
        Self {
            density_kg_m3: 998.0,
            viscosity_pa_s: 1.002e-3,
            diameter_m: 0.1,
            length_m: 100.0,
            roughness_m: 4.6e-5,
            velocity_m_s: 2.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Pipe Flow Workbench right-side panel. A no-op when the
/// `show_pipeflow_workbench` toggle is off.
pub fn draw_pipeflow_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_pipeflow_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_pipeflow_workbench",
        "Pipe Flow",
        |app, ui| {
            ui.label(
                egui::RichText::new("native Darcy-Weisbach pipe-flow analysis · valenx-pipeflow")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.pipeflow;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Fluid").strong());
                    ui.horizontal(|ui| {
                        ui.label("density ρ (kg/m³)");
                        ui.add(egui::DragValue::new(&mut s.density_kg_m3).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("viscosity μ (Pa·s)");
                        ui.add(egui::DragValue::new(&mut s.viscosity_pa_s).speed(1.0e-4));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Pipe").strong());
                    ui.horizontal(|ui| {
                        ui.label("diameter D (m)");
                        ui.add(egui::DragValue::new(&mut s.diameter_m).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("length L (m)");
                        ui.add(egui::DragValue::new(&mut s.length_m).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("roughness ε (m)");
                        ui.add(egui::DragValue::new(&mut s.roughness_m).speed(1.0e-5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Flow").strong());
                    ui.horizontal(|ui| {
                        ui.label("velocity V (m/s)");
                        ui.add(egui::DragValue::new(&mut s.velocity_m_s).speed(0.05));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_pipe(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D pipe").strong())
                        .on_hover_text(
                            "Build a flanged pipe run on saddle supports as a 3-D solid (representative geometry, not to scale) and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Result").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_pipeflow_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.pipeflow` borrow is
    // released here): build the pipe's 3-D solid and load it.
    if app.pipeflow.show_3d_request {
        app.pipeflow.show_3d_request = false;
        load_pipe_3d(app);
    }
}

/// Validate the form, solve the pipe flow and format the readout.
fn run_pipe(s: &mut PipeFlowWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Solve the pipe flow and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &PipeFlowWorkbenchState) -> Result<String, String> {
    let rel_roughness = s.roughness_m / s.diameter_m;
    let r = solve_pipe(
        s.density_kg_m3,
        s.viscosity_pa_s,
        s.diameter_m,
        s.length_m,
        rel_roughness,
        s.velocity_m_s,
    )
    .map_err(|e| e.to_string())?;

    let regime = match r.friction.regime {
        FlowRegime::Laminar => "Laminar",
        FlowRegime::Transitional => "Transitional",
        FlowRegime::Turbulent => "Turbulent",
    };
    let area_m2 = 0.25 * PI * s.diameter_m * s.diameter_m;
    let flow_m3s = s.velocity_m_s * area_m2;

    Ok(format!(
        "pipe D / L      : {:.3} m / {:.1} m\n\
         roughness ε/D   : {:.2e}\n\
         velocity V      : {:.2} m/s\n\n\
         Reynolds Re     : {:.3e}\n\
         regime          : {}\n\
         friction factor : {:.4}\n\
         flow rate Q     : {:.4} m³/s  ({:.1} L/s)\n\
         head loss hf    : {:.3} m\n\
         pressure drop   : {:.2} kPa",
        s.diameter_m,
        s.length_m,
        rel_roughness,
        s.velocity_m_s,
        r.friction.reynolds,
        regime,
        r.friction.friction_factor,
        flow_m3s,
        flow_m3s * 1000.0,
        r.head_loss_m,
        r.pressure_drop_pa / 1000.0,
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
            hi + jn,
            hi + j,
            lo + j,
            lo + jn,
            hi + jn,
        ]);
    }
}

/// Build the pipe run as a triangle [`Mesh`] — a horizontal pipe with a
/// raised flange at each end on two saddle supports and a base plate.
/// Representative geometry (not to scale; the physics is the
/// `valenx-pipeflow` Darcy-Weisbach solve). `None` for an inconsistent
/// configuration.
fn pipe_solid_mesh(s: &PipeFlowWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a flow that actually solves.
    let rel_roughness = s.roughness_m / s.diameter_m;
    solve_pipe(
        s.density_kg_m3,
        s.viscosity_pa_s,
        s.diameter_m,
        s.length_m,
        rel_roughness,
        s.velocity_m_s,
    )
    .ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    let axis_z = 0.55;

    // Pipe barrel.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-1.2, 0.0, axis_z),
        2.4,
        0.22,
        28,
    );
    // End flanges (short fat rings).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-1.25, 0.0, axis_z),
        0.12,
        0.34,
        28,
    );
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(1.13, 0.0, axis_z),
        0.12,
        0.34,
        28,
    );
    // Saddle supports under the pipe.
    for &x in &[-0.6, 0.6] {
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(x, 0.0, 0.24),
            Vector3::new(0.1, 0.28, 0.24),
        );
    }
    // Base plate.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.04),
        Vector3::new(1.0, 0.4, 0.04),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-pipeflow");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D pipe solid and load it into the central viewport.
fn load_pipe_3d(app: &mut ValenxApp) {
    let Some(mesh) = pipe_solid_mesh(&app.pipeflow) else {
        app.pipeflow.error =
            Some("flow parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<pipe>/valenx-pipeflow"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"pipeflow"}`** product: the canonical
/// Darcy-Weisbach pipe run built as a 3-D solid, paired with the workbench's
/// own `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`PipeFlowWorkbenchState::default`].
pub(crate) fn pipeflow_product() -> crate::WorkspaceProduct {
    let s = PipeFlowWorkbenchState::default();
    let mesh = pipe_solid_mesh(&s).expect("canonical pipe flow ⇒ solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<pipeflow>/valenx-pipe");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical pipe flow ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Pipe flow (Darcy-Weisbach)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        animation: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = PipeFlowWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_reynolds_and_headloss() {
        let mut s = PipeFlowWorkbenchState::default();
        run_pipe(&mut s);
        assert!(
            s.error.is_none(),
            "default pipe should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("Reynolds"));
        assert!(s.result.contains("friction factor"));
        assert!(s.result.contains("head loss"));
        assert!(s.result.contains("pressure drop"));
        // Default 2 m/s through 100 mm water is firmly turbulent.
        assert!(s.result.contains("Turbulent"));
    }

    #[test]
    fn analyze_rejects_zero_diameter() {
        let mut s = PipeFlowWorkbenchState {
            diameter_m: 0.0,
            ..Default::default()
        };
        run_pipe(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn laminar_friction_is_64_over_re_at_low_velocity() {
        // Ground truth: below the transition Reynolds number the Darcy
        // friction factor is exactly 64/Re. At 0.02 m/s, Re ≈ 1990.
        let r = solve_pipe(998.0, 1.002e-3, 0.1, 100.0, 4.6e-5 / 0.1, 0.02).unwrap();
        assert_eq!(r.friction.regime, FlowRegime::Laminar);
        let expected_f = 64.0 / r.friction.reynolds;
        assert!(
            (r.friction.friction_factor - expected_f).abs() < 1e-9,
            "f={} expected {}",
            r.friction.friction_factor,
            expected_f
        );
    }

    #[test]
    fn pipe_mesh_for_default_is_nonempty_and_in_range() {
        let s = PipeFlowWorkbenchState::default();
        let mesh = pipe_solid_mesh(&s).expect("default pipe yields a solid");
        assert!(mesh.nodes.len() > 8, "expected barrel + flanges + supports");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn pipe_mesh_none_for_invalid() {
        let s = PipeFlowWorkbenchState {
            diameter_m: 0.0,
            ..Default::default()
        };
        assert!(pipe_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_pipeflow_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_pipeflow_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_pipeflow_workbench = true;
        run_pipe(&mut app.pipeflow);
        draw_workbench(&mut app);
    }
}
