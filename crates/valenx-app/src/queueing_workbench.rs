//! The right-side **Queueing Workbench** panel — native single-server
//! M/M/1 steady-state analysis over `valenx-queueing`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_queueing_workbench`,
//! toggled from the View menu. The form sets the arrival rate `λ` and (in
//! "given μ" mode) the service rate `μ`, or (in "size μ for a target
//! response time" mode) a target mean time-in-system `W`. "Analyze"
//! validates the load, solves the closed-form M/M/1 chain and reports the
//! utilization `ρ`, the mean numbers in system / in queue (`L`, `Lq`), the
//! mean times in system / waiting (`W`, `Wq`) and the idle probability
//! `P(0)`; "Show 3-D" loads a representative queue — a row of waiting
//! "customer" boxes feeding one server box — into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_queueing::Mm1;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which M/M/1 quantity the form supplies — either both rates directly, or
/// the arrival rate plus a target response time from which the required
/// service rate is sized via
/// [`Mm1::service_rate_for_mean_response_time`].
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Analysis {
    /// Supply both the arrival rate `λ` and the service rate `μ` and solve
    /// the queue directly.
    #[default]
    GivenRates,
    /// Supply the arrival rate `λ` and a target mean time-in-system `W`,
    /// size the single-server rate `μ = λ + 1/W` that meets it, then solve.
    SizeServiceRate,
}

/// Persistent form + result state for the Queueing Workbench.
pub struct QueueingWorkbenchState {
    /// Which quantity to supply (rates directly, or size `μ` for a target
    /// response time).
    analysis: Analysis,
    /// Mean arrival rate `λ` (customers per unit time).
    lambda: f64,
    /// Mean service rate `μ` (customers per unit time) — used in the
    /// `GivenRates` mode.
    mu: f64,
    /// Target mean time in system `W` (per-unit-time) — used in the
    /// `SizeServiceRate` mode.
    target_w: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D queue solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for QueueingWorkbenchState {
    fn default() -> Self {
        // λ = 8, μ = 10  =>  ρ = 0.8, L = ρ/(1-ρ) = 4 customers in system.
        Self {
            analysis: Analysis::GivenRates,
            lambda: 8.0,
            mu: 10.0,
            target_w: 0.5,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Queueing Workbench right-side panel. A no-op when the
/// `show_queueing_workbench` toggle is off.
pub fn draw_queueing_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_queueing_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_queueing_workbench",
        "Queueing (M/M/1)",
        |app, ui| {
            ui.label(
                egui::RichText::new("native single-server steady-state queue · valenx-queueing")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.queueing;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Mode").strong());
                    ui.radio_value(&mut s.analysis, Analysis::GivenRates, "given μ");
                    ui.radio_value(
                        &mut s.analysis,
                        Analysis::SizeServiceRate,
                        "size μ for target W",
                    );

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Rates").strong());
                    ui.horizontal(|ui| {
                        ui.label("arrival rate λ");
                        ui.add(egui::DragValue::new(&mut s.lambda).speed(0.25));
                    });
                    match s.analysis {
                        Analysis::GivenRates => {
                            ui.horizontal(|ui| {
                                ui.label("service rate μ");
                                ui.add(egui::DragValue::new(&mut s.mu).speed(0.25));
                            });
                        }
                        Analysis::SizeServiceRate => {
                            ui.horizontal(|ui| {
                                ui.label("target time in system W");
                                ui.add(egui::DragValue::new(&mut s.target_w).speed(0.05));
                            });
                        }
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_queueing(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative single-server queue (a row of waiting \"customer\" boxes feeding one server box) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Steady state").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_queueing_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.queueing` borrow is
    // released here): build the queue's 3-D solid and load it.
    if app.queueing.show_3d_request {
        app.queueing.show_3d_request = false;
        load_queue_3d(app);
    }
}

/// Validate the form, solve the queue and format the readout.
fn run_queueing(s: &mut QueueingWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Mm1`] queue the form describes — either from both
/// rates directly, or from the arrival rate plus the service rate sized to
/// hit the target response time. Extracted so it is unit-testable and
/// shared with the 3-D gate.
fn build_queue(s: &QueueingWorkbenchState) -> Result<Mm1, String> {
    let mu = match s.analysis {
        Analysis::GivenRates => s.mu,
        Analysis::SizeServiceRate => Mm1::service_rate_for_mean_response_time(s.lambda, s.target_w)
            .map_err(|e| e.to_string())?,
    };
    Mm1::new(s.lambda, mu).map_err(|e| e.to_string())
}

/// Solve the queue and format the full steady-state readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &QueueingWorkbenchState) -> Result<String, String> {
    let q = build_queue(s)?;
    Ok(format!(
        "arrival rate λ  : {lambda:.4}\n\
         service rate μ  : {mu:.4}\n\
         utilization ρ   : {rho:.4}\n\n\
         L  (in system)  : {l:.4}\n\
         Lq (in queue)   : {lq:.4}\n\
         W  (time sys)   : {w:.4}\n\
         Wq (wait queue) : {wq:.4}\n\
         P(0) idle       : {p0:.4}\n\
         Little residual : {residual:.2e}",
        lambda = q.lambda(),
        mu = q.mu(),
        rho = q.rho(),
        l = q.l(),
        lq = q.lq(),
        w = q.w(),
        wq = q.wq(),
        p0 = q.p0(),
        residual = q.little_residual(),
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

/// Build the queue as a triangle [`Mesh`] — a row of waiting "customer"
/// boxes (the line) feeding a larger server box at the head, on a base.
/// Representative geometry (not to scale; the steady-state numbers are the
/// `valenx-queueing` result). `None` for an invalid / unstable
/// configuration.
fn queue_solid_mesh(s: &QueueingWorkbenchState) -> Option<Mesh> {
    build_queue(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // A fixed, representative line of waiting customers along -x, with the
    // server at the head (+x). Independent of scale.
    let customers = 5;
    for i in 0..customers {
        let x = -0.3 - (i as f64) * 0.32;
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(x, 0.0, 0.18),
            Vector3::new(0.11, 0.11, 0.11),
        );
    }
    // Server box at the head of the line.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.18, 0.0, 0.24),
        Vector3::new(0.18, 0.2, 0.18),
    );
    // Base / floor under the whole queue.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.7, 0.0, 0.02),
        Vector3::new(1.1, 0.3, 0.02),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-queueing");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D queue solid and load it into the central viewport.
fn load_queue_3d(app: &mut ValenxApp) {
    let Some(mesh) = queue_solid_mesh(&app.queueing) else {
        app.queueing.error =
            Some("queue parameters are invalid or unstable — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<queue>/valenx-queueing"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"queueing"}`** product: the canonical
/// queue built as a 3-D solid, paired with the workbench's own `compute()`
/// readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`QueueingWorkbenchState::default`].
pub(crate) fn queueing_product() -> crate::WorkspaceProduct {
    let s = QueueingWorkbenchState::default();
    let mesh = queue_solid_mesh(&s).expect("canonical queue ⇒ solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<queue>/valenx-queueing");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical queue ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Queueing (M/M/c)".into(),
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
        let s = QueueingWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_rho_and_l() {
        let mut s = QueueingWorkbenchState::default();
        run_queueing(&mut s);
        assert!(
            s.error.is_none(),
            "default queue should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("utilization ρ"));
        assert!(s.result.contains("L  (in system)"));
        // Ground truth at λ = 8, μ = 10: ρ = λ/μ = 0.8, and
        // L = ρ/(1-ρ) = 0.8/0.2 = 4 customers in system. Both print at
        // the readout's :.4 precision.
        assert!(s.result.contains("0.8000"));
        assert!(s.result.contains("4.0000"));
    }

    #[test]
    fn analyze_rejects_unstable_load() {
        // λ >= μ has no finite steady state and must be reported.
        let mut s = QueueingWorkbenchState {
            lambda: 10.0,
            mu: 8.0,
            ..Default::default()
        };
        run_queueing(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn ground_truth_rho_and_l_at_lambda8_mu10() {
        // Hand-computed M/M/1 ground truth, independent of the formatter:
        // ρ = λ/μ = 8/10 = 0.8, L = ρ/(1-ρ) = 4, P(0) = 1-ρ = 0.2.
        let q = Mm1::new(8.0, 10.0).unwrap();
        assert!((q.rho() - 0.8).abs() < 1e-12);
        assert!((q.l() - 4.0).abs() < 1e-12);
        assert!((q.p0() - 0.2).abs() < 1e-12);
    }

    #[test]
    fn size_service_rate_mode_meets_target_w() {
        // In sizing mode μ = λ + 1/W, so the solved queue reproduces W.
        let s = QueueingWorkbenchState {
            analysis: Analysis::SizeServiceRate,
            lambda: 8.0,
            target_w: 0.5,
            ..Default::default()
        };
        let q = build_queue(&s).expect("sizing mode yields a stable queue");
        // μ = 8 + 1/0.5 = 10, and W = 1/(μ-λ) = 1/2 = 0.5.
        assert!((q.mu() - 10.0).abs() < 1e-12);
        assert!((q.w() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn queue_mesh_for_default_is_nonempty_and_in_range() {
        let s = QueueingWorkbenchState::default();
        let mesh = queue_solid_mesh(&s).expect("default queue yields a solid");
        assert!(mesh.nodes.len() > 8, "expected customers + server + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn queue_mesh_none_for_unstable() {
        let s = QueueingWorkbenchState {
            lambda: 10.0,
            mu: 8.0,
            ..Default::default()
        };
        assert!(queue_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_queueing_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_queueing_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_queueing_workbench = true;
        run_queueing(&mut app.queueing);
        draw_workbench(&mut app);
    }
}
