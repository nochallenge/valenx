//! The right-side **Pipe Network Workbench** panel — native looped
//! pipe-network flow balancing over `valenx-pipenetwork`.
//!
//! Mirrors the Pipe Flow / Heat Transfer workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_pipenetwork_workbench`,
//! toggled from the View menu. The form sets a single closed loop of two
//! parallel pipes (A → B) by their resistance coefficients `k0`, `k1` and a
//! total inflow `Q`; "Analyze" splits `Q` as a continuity-respecting initial
//! guess and runs the **Hardy-Cross** iteration ([`valenx_pipenetwork::Network::solve`])
//! to balance the loop, reporting the converged flow in each pipe, the
//! per-pipe head losses, the residual loop head loss and the solver's
//! iteration count, and "Show 3-D network" loads a representative
//! parallel-loop solid (two pipes between two junction blocks) into the
//! central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_pipenetwork::{Loop, LoopMember, Network, Pipe, SolveConfig};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Pipe Network Workbench.
pub struct PipeNetworkWorkbenchState {
    /// Resistance coefficient `k0` of pipe 0 in the loss law `h = k q |q|`.
    k0: f64,
    /// Resistance coefficient `k1` of pipe 1.
    k1: f64,
    /// Total volumetric inflow `Q` entering node A and drawn at node B
    /// (m³/s); split between the two parallel pipes.
    total_inflow: f64,
    /// Convergence tolerance on the largest loop correction `|dQ|`.
    tolerance: f64,
    /// Maximum Hardy-Cross sweeps before giving up.
    max_iterations: usize,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D network solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for PipeNetworkWorkbenchState {
    fn default() -> Self {
        // The crate's canonical parallel-loop example: two pipes A → B with
        // k0 = 1, k1 = 4 and a total inflow Q = 3. Balancing k0 q0² = k1 q1²
        // under q0 + q1 = 3 gives the exact split q0 = 2, q1 = 1
        // (1·2² = 4·1² = 4). Hardy-Cross from the even guess (1.5, 1.5)
        // converges to it in a handful of sweeps.
        Self {
            k0: 1.0,
            k1: 4.0,
            total_inflow: 3.0,
            tolerance: 1e-9,
            max_iterations: 200,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Pipe Network Workbench right-side panel. A no-op when the
/// `show_pipenetwork_workbench` toggle is off.
pub fn draw_pipenetwork_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_pipenetwork_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_pipenetwork_workbench",
        "Pipe Network",
        |app, ui| {
            ui.label(
                egui::RichText::new("native pipe-network flow balancing · valenx-pipenetwork")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.pipenetwork;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Parallel loop A → B").strong());
                    ui.horizontal(|ui| {
                        ui.label("pipe 0 resistance k₀");
                        ui.add(egui::DragValue::new(&mut s.k0).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("pipe 1 resistance k₁");
                        ui.add(egui::DragValue::new(&mut s.k1).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("total inflow Q (m³/s)");
                        ui.add(egui::DragValue::new(&mut s.total_inflow).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Solver (Hardy-Cross)").strong());
                    ui.horizontal(|ui| {
                        ui.label("tolerance |dQ|");
                        ui.add(egui::DragValue::new(&mut s.tolerance).speed(1.0e-9));
                    });
                    ui.horizontal(|ui| {
                        ui.label("max iterations");
                        ui.add(egui::DragValue::new(&mut s.max_iterations).speed(1.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_network(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D network").strong())
                        .on_hover_text(
                            "Build the parallel loop as two pipes between two junction blocks as a 3-D solid (representative geometry, not to scale) and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Balanced flows").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_pipenetwork_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.pipenetwork` borrow is
    // released here): build the network's 3-D solid and load it.
    if app.pipenetwork.show_3d_request {
        app.pipenetwork.show_3d_request = false;
        load_network_3d(app);
    }
}

/// Validate the form, balance the network and format the readout.
fn run_network(s: &mut PipeNetworkWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Assemble the two-pipe parallel loop from the form, returning the
/// constructed (unsolved) [`Network`]. The two pipes share node A (inflow)
/// and node B (outflow); pipe 0 is traversed forward and pipe 1 reverse
/// around the single loop. The initial flows split `Q` evenly, which
/// satisfies continuity at both nodes. Extracted so it is shared by the
/// solver and the 3-D gate.
fn build_network(s: &PipeNetworkWorkbenchState) -> Result<Network, String> {
    let guess = s.total_inflow / 2.0;
    let pipes = vec![
        Pipe::new(s.k0, guess).map_err(|e| e.to_string())?,
        Pipe::new(s.k1, guess).map_err(|e| e.to_string())?,
    ];
    let lp = Loop::new(
        "loop-AB",
        vec![LoopMember::forward(0), LoopMember::reverse(1)],
    )
    .map_err(|e| e.to_string())?;
    Network::new(pipes, vec![lp]).map_err(|e| e.to_string())
}

/// Build, balance and format the full readout, mapping any domain error to a
/// display string. Extracted so it is unit-testable.
fn compute(s: &PipeNetworkWorkbenchState) -> Result<String, String> {
    let mut net = build_network(s)?;
    let config = SolveConfig {
        tolerance: s.tolerance,
        max_iterations: s.max_iterations,
    };
    let report = net.solve(&config).map_err(|e| e.to_string())?;

    let q0 = net.pipes()[0].q;
    let q1 = net.pipes()[1].q;
    let h0 = net.pipes()[0].head_loss();
    let h1 = net.pipes()[1].head_loss();
    // Signed head loss summed around the loop — driven to ~0 at balance.
    let loop_residual = net.loop_head_loss(0).unwrap_or(f64::NAN);

    Ok(format!(
        "resistances k    : {:.3} / {:.3}\n\
         total inflow Q   : {:.3} m³/s\n\n\
         pipe 0 flow q₀   : {:.4} m³/s\n\
         pipe 1 flow q₁   : {:.4} m³/s\n\
         sum q₀ + q₁      : {:.4} m³/s\n\n\
         pipe 0 head loss : {:.4} m\n\
         pipe 1 head loss : {:.4} m\n\
         loop residual    : {:.2e} m\n\n\
         iterations       : {}\n\
         final |dQ|       : {:.2e}",
        s.k0,
        s.k1,
        s.total_inflow,
        q0,
        q1,
        q0 + q1,
        h0,
        h1,
        loop_residual,
        report.iterations,
        report.final_residual,
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

/// Build the parallel loop as a triangle [`Mesh`] — two horizontal pipes
/// (the loop's two branches) running between two junction blocks (node A on
/// the left, node B on the right), with a base plate. Representative
/// geometry (not to scale; the flows are the `valenx-pipenetwork`
/// Hardy-Cross result). `None` for a configuration that does not build /
/// solve.
fn network_solid_mesh(s: &PipeNetworkWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a network that actually assembles and balances.
    let mut net = build_network(s).ok()?;
    let config = SolveConfig {
        tolerance: s.tolerance,
        max_iterations: s.max_iterations,
    };
    net.solve(&config).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    let z = 0.5;

    // Junction A (left) and junction B (right) blocks.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-1.1, 0.0, z),
        Vector3::new(0.12, 0.5, 0.16),
    );
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(1.1, 0.0, z),
        Vector3::new(0.12, 0.5, 0.16),
    );
    // Two parallel pipe branches between the junctions (the loop's pipe 0
    // on the +y side, pipe 1 on the -y side).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-1.0, 0.32, z),
        2.0,
        0.1,
        24,
    );
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-1.0, -0.32, z),
        2.0,
        0.1,
        24,
    );
    // Base plate.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.04),
        Vector3::new(1.3, 0.6, 0.04),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-pipenetwork");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D network solid and load it into the central viewport.
fn load_network_3d(app: &mut ValenxApp) {
    let Some(mesh) = network_solid_mesh(&app.pipenetwork) else {
        app.pipenetwork.error =
            Some("network parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<network>/valenx-pipenetwork"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"pipenetwork"}`** product: the canonical
/// looped pipe network built as a 3-D solid, paired with the workbench's own
/// `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`PipeNetworkWorkbenchState::default`].
pub(crate) fn pipenetwork_product() -> crate::WorkspaceProduct {
    let s = PipeNetworkWorkbenchState::default();
    let mesh = network_solid_mesh(&s).expect("canonical pipe network ⇒ solid builds");
    let loaded =
        crate::products_registry::loaded_mesh_from(mesh, "<pipenetwork>/valenx-pipenetwork");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical pipe network ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Pipe network (loop balancing)".into(),
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
        let s = PipeNetworkWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_flows_and_iterations() {
        let mut s = PipeNetworkWorkbenchState::default();
        run_network(&mut s);
        assert!(
            s.error.is_none(),
            "default network should balance: {:?}",
            s.error
        );
        assert!(s.result.contains("pipe 0 flow"));
        assert!(s.result.contains("pipe 1 flow"));
        assert!(s.result.contains("loop residual"));
        assert!(s.result.contains("iterations"));
        // Canonical k=(1,4), Q=3 split balances to q0=2, q1=1.
        assert!(s.result.contains("2.0000"));
        assert!(s.result.contains("1.0000"));
    }

    #[test]
    fn analyze_rejects_negative_resistance() {
        let mut s = PipeNetworkWorkbenchState {
            k0: -1.0,
            ..Default::default()
        };
        run_network(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn balanced_loop_conserves_mass_and_has_zero_head_loss() {
        // Ground truth: at the Hardy-Cross balance the two branch flows must
        // (a) sum to the total inflow Q (node mass conservation, preserved
        // exactly from the even initial guess), and (b) leave zero net head
        // loss around the loop (k0 q0² = k1 q1²).
        let s = PipeNetworkWorkbenchState::default();
        let mut net = build_network(&s).unwrap();
        net.solve(&SolveConfig::default()).unwrap();
        let q0 = net.pipes()[0].q;
        let q1 = net.pipes()[1].q;
        // Continuity: q0 + q1 == Q.
        assert!((q0 + q1 - s.total_inflow).abs() < 1e-7, "q0+q1={}", q0 + q1);
        // Loop balance: net signed head loss ~ 0, i.e. k0 q0² == k1 q1².
        assert!(net.loop_head_loss(0).unwrap().abs() < 1e-7);
        assert!((s.k0 * q0 * q0 - s.k1 * q1 * q1).abs() < 1e-6);
        // And the analytic split for k=(1,4), Q=3 is exactly (2, 1).
        assert!((q0 - 2.0).abs() < 1e-7);
        assert!((q1 - 1.0).abs() < 1e-7);
    }

    #[test]
    fn network_mesh_for_default_is_nonempty_and_in_range() {
        let s = PipeNetworkWorkbenchState::default();
        let mesh = network_solid_mesh(&s).expect("default network yields a solid");
        assert!(
            mesh.nodes.len() > 8,
            "expected two junctions + two pipes + base"
        );
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn network_mesh_none_for_invalid() {
        let s = PipeNetworkWorkbenchState {
            k0: -1.0,
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
            draw_pipenetwork_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_pipenetwork_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_pipenetwork_workbench = true;
        run_network(&mut app.pipenetwork);
        draw_workbench(&mut app);
    }
}
