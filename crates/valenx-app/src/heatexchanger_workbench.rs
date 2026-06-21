//! The right-side **Heat Exchanger Workbench** panel — native
//! effectiveness-NTU analysis over `valenx-heatexchanger`.
//!
//! Mirrors the Pump / Solar PV workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_heatexchanger_workbench`,
//! toggled from the View menu. The form drives an
//! [`valenx_heatexchanger::NtuProblem`] (stream heat-capacity rates, the
//! overall conductance UA, inlet temperatures and the flow arrangement);
//! "Analyze" reports the capacity ratio, NTU, effectiveness, the maximum
//! and actual duty, both outlet temperatures, the log-mean temperature
//! difference (LMTD) driving the duty and the conductance `UA = Q / LMTD`
//! it implies (a self-consistency cross-check against the input UA), and
//! "Show 3-D exchanger" loads a shell-and-tube solid into the central
//! viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_heatexchanger::lmtd::lmtd;
use valenx_heatexchanger::ntu::{solve, NtuProblem};
use valenx_heatexchanger::{FlowArrangement, TerminalTemperatures};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Number of tubes drawn in the 3-D bundle (one central + a hex ring).
const TUBE_COUNT: usize = 7;

/// Persistent form + result state for the Heat Exchanger Workbench.
pub struct HeatExchangerWorkbenchState {
    /// Hot-stream heat-capacity rate `Ch = m_dot · cp` (W/K).
    c_hot_w_per_k: f64,
    /// Cold-stream heat-capacity rate `Cc = m_dot · cp` (W/K).
    c_cold_w_per_k: f64,
    /// Overall conductance `UA = U · A` (W/K).
    ua_w_per_k: f64,
    /// Hot-stream inlet temperature (°C).
    hot_in_c: f64,
    /// Cold-stream inlet temperature (°C).
    cold_in_c: f64,
    /// Flow arrangement (counter- or parallel-flow).
    arrangement: FlowArrangement,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D exchanger solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for HeatExchangerWorkbenchState {
    fn default() -> Self {
        // A representative liquid-liquid exchanger: hot 90 °C / 2 kW·K⁻¹
        // against cold 20 °C / 3 kW·K⁻¹ through UA = 5 kW/K, which gives
        // NTU = 2.5, ε ≈ 0.80 counterflow.
        Self {
            c_hot_w_per_k: 2000.0,
            c_cold_w_per_k: 3000.0,
            ua_w_per_k: 5000.0,
            hot_in_c: 90.0,
            cold_in_c: 20.0,
            arrangement: FlowArrangement::Counterflow,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Heat Exchanger Workbench right-side panel. A no-op when the
/// `show_heatexchanger_workbench` toggle is off.
pub fn draw_heatexchanger_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_heatexchanger_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_heatexchanger_workbench",
        "Heat Exchanger",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native effectiveness-NTU exchanger analysis · valenx-heatexchanger",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.heatexchanger;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Stream capacity rates").strong());
                    ui.horizontal(|ui| {
                        ui.label("hot  Ch (W/K)");
                        ui.add(egui::DragValue::new(&mut s.c_hot_w_per_k).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("cold Cc (W/K)");
                        ui.add(egui::DragValue::new(&mut s.c_cold_w_per_k).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("conductance UA (W/K)");
                        ui.add(egui::DragValue::new(&mut s.ua_w_per_k).speed(50.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Inlet temperatures").strong());
                    ui.horizontal(|ui| {
                        ui.label("hot  in (°C)");
                        ui.add(egui::DragValue::new(&mut s.hot_in_c).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("cold in (°C)");
                        ui.add(egui::DragValue::new(&mut s.cold_in_c).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Arrangement").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(
                            &mut s.arrangement,
                            FlowArrangement::Counterflow,
                            "Counterflow",
                        );
                        ui.radio_value(
                            &mut s.arrangement,
                            FlowArrangement::ParallelFlow,
                            "Parallel",
                        );
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_exchanger(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D exchanger").strong())
                        .on_hover_text(
                            "Build a shell-and-tube exchanger (shell, tube bundle, channel heads and nozzles) as a 3-D solid and load it into the central viewport to orbit",
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
        app.show_heatexchanger_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.heatexchanger` borrow
    // is released here): build the exchanger's 3-D solid and load it.
    if app.heatexchanger.show_3d_request {
        app.heatexchanger.show_3d_request = false;
        load_exchanger_3d(app);
    }
}

/// Validate the form, solve the exchanger and format the readout.
fn run_exchanger(s: &mut HeatExchangerWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build a validated [`NtuProblem`] from the form, mapping the domain
/// error to a display string.
fn build_problem(s: &HeatExchangerWorkbenchState) -> Result<NtuProblem, String> {
    NtuProblem::new(
        s.c_hot_w_per_k,
        s.c_cold_w_per_k,
        s.ua_w_per_k,
        s.hot_in_c,
        s.cold_in_c,
    )
    .map_err(|e| e.to_string())
}

/// Solve the exchanger and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &HeatExchangerWorkbenchState) -> Result<String, String> {
    let problem = build_problem(s)?;
    let r = solve(&problem, s.arrangement).map_err(|e| e.to_string())?;
    // The four solved terminal temperatures close the loop to the LMTD
    // method: build the validated terminal set and take the log-mean
    // driving temperature difference for this arrangement.
    let temps = TerminalTemperatures::new(s.hot_in_c, r.hot_out, s.cold_in_c, r.cold_out)
        .map_err(|e| e.to_string())?;
    let lmtd_c = lmtd(&temps, s.arrangement).map_err(|e| e.to_string())?;
    // Conductance implied by the duty over that LMTD (Q = UA · LMTD).
    // For a self-consistent solve this round-trips back to the input UA,
    // so it doubles as a cross-check on the effectiveness-NTU result.
    let ua_from_lmtd = r.q_w / lmtd_c;
    Ok(format!(
        "arrangement     : {}\n\
         Cmin / Cmax     : {:.0} / {:.0} W/K\n\
         capacity ratio  : {:.3}\n\
         NTU             : {:.3}\n\
         effectiveness ε : {:.1} %\n\n\
         max duty Q_max  : {:.2} kW\n\
         duty Q          : {:.2} kW\n\
         hot  in → out   : {:.1} → {:.1} °C\n\
         cold in → out   : {:.1} → {:.1} °C\n\
         LMTD ΔT_lm      : {:.2} °C\n\
         UA = Q / LMTD   : {:.0} W/K",
        s.arrangement.label(),
        problem.c_min(),
        problem.c_max(),
        r.capacity_ratio,
        r.ntu,
        r.effectiveness * 100.0,
        r.q_max_w / 1000.0,
        r.q_w / 1000.0,
        s.hot_in_c,
        r.hot_out,
        s.cold_in_c,
        r.cold_out,
        lmtd_c,
        ua_from_lmtd,
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

/// Append a (double-sided) cylinder whose axis runs along `+z`, spanning
/// `base.z ..= base.z + height` with circle centre `(base.x, base.y)`.
fn push_cyl_z(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    height: f64,
    r: f64,
    seg: usize,
) {
    let (z0, z1) = (base.z, base.z + height);
    let bot = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(base.x + r * a.cos(), base.y + r * a.sin(), z0));
    }
    let top = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(base.x + r * a.cos(), base.y + r * a.sin(), z1));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[
            bot + j,
            top + j,
            top + jn,
            bot + j,
            top + jn,
            bot + jn,
            bot + j,
            top + jn,
            top + j,
            bot + j,
            bot + jn,
            top + jn,
        ]);
    }
}

/// Build the exchanger as a triangle [`Mesh`] — a horizontal shell with a
/// tube bundle (one central tube + a hex ring), a channel head at each
/// end, two shell-side nozzles and a saddle base. Representative geometry
/// (the performance is the `valenx-heatexchanger` effectiveness-NTU
/// result). `None` for an inconsistent configuration.
fn exchanger_solid_mesh(s: &HeatExchangerWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a consistent exchanger configuration.
    build_problem(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    let axis_z = 0.6;

    // Shell — the large horizontal cylinder.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-1.0, 0.0, axis_z),
        2.0,
        0.4,
        28,
    );
    // Tube bundle — one central tube plus a hex ring, running the full
    // length and protruding slightly past the heads.
    for k in 0..TUBE_COUNT {
        let (cy, cz) = if k == 0 {
            (0.0, axis_z)
        } else {
            let a = (k - 1) as f64 / 6.0 * TAU;
            (0.22 * a.cos(), axis_z + 0.22 * a.sin())
        };
        push_cyl_x(
            &mut nodes,
            &mut tris,
            Vector3::new(-1.1, cy, cz),
            2.2,
            0.06,
            12,
        );
    }
    // Channel heads — a slightly larger short cylinder at each end.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-1.18, 0.0, axis_z),
        0.18,
        0.42,
        28,
    );
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(1.0, 0.0, axis_z),
        0.18,
        0.42,
        28,
    );
    // Shell-side nozzles — one up near each end.
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.7, 0.0, axis_z + 0.35),
        0.4,
        0.1,
        16,
    );
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.7, 0.0, axis_z + 0.35),
        0.4,
        0.1,
        16,
    );
    // Saddle base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.06),
        Vector3::new(0.9, 0.45, 0.06),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-heatexchanger");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D exchanger solid and load it into the central viewport.
fn load_exchanger_3d(app: &mut ValenxApp) {
    let Some(mesh) = exchanger_solid_mesh(&app.heatexchanger) else {
        app.heatexchanger.error =
            Some("exchanger configuration is inconsistent — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<exchanger>/valenx-heatexchanger"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical heat-exchanger workbench as a 3-D solid
/// plus its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn heatexchanger_product() -> crate::WorkspaceProduct {
    let s = HeatExchangerWorkbenchState::default();
    let mesh =
        exchanger_solid_mesh(&s).expect("canonical heat exchanger ⇒ shell-tube solid builds");
    let loaded =
        crate::products_registry::loaded_mesh_from(mesh, "<heatexchanger>/valenx-exchanger");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical heat exchanger ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Heat exchanger (LMTD/duty)".into(),
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
        let s = HeatExchangerWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_effectiveness_and_outlets() {
        let mut s = HeatExchangerWorkbenchState::default();
        run_exchanger(&mut s);
        assert!(
            s.error.is_none(),
            "default exchanger should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("effectiveness"));
        assert!(s.result.contains("NTU"));
        assert!(s.result.contains("hot  in → out"));
        assert!(s.result.contains("cold in → out"));
    }

    #[test]
    fn analyze_rejects_hot_inlet_below_cold() {
        let mut s = HeatExchangerWorkbenchState {
            hot_in_c: 20.0,
            cold_in_c: 90.0,
            ..Default::default()
        };
        run_exchanger(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_default_reports_lmtd_and_implied_ua() {
        // Ground truth for the default counterflow exchanger
        // (Ch = 2000, Cc = 3000, UA = 5000 W/K, hot 90 °C, cold 20 °C):
        //   Cr = 2/3, NTU = 2.5, ε = 0.796040…, Q = 111.4456… kW,
        //   hot_out = 34.2772 °C, cold_out = 57.1485 °C, so the two
        //   counterflow approaches are
        //     dT1 = 90 − 57.1485 = 32.8515,  dT2 = 34.2772 − 20 = 14.2772,
        //   giving LMTD = (dT1 − dT2)/ln(dT1/dT2) = 22.2891… °C and a
        //   back-computed UA = Q/LMTD = 5000 W/K (it must round-trip to
        //   the input UA, the Q = UA·LMTD self-consistency identity).
        let s = HeatExchangerWorkbenchState::default();
        let out = compute(&s).expect("default exchanger computes");
        assert!(
            out.contains("LMTD ΔT_lm      : 22.29 °C"),
            "LMTD readout wrong:\n{out}"
        );
        assert!(
            out.contains("UA = Q / LMTD   : 5000 W/K"),
            "implied-UA readout wrong:\n{out}"
        );

        // Independent hand check of the LMTD against the closed form,
        // recomputing the duty from ε·q_max rather than reusing the
        // formatted string.
        let problem = build_problem(&s).unwrap();
        let r = solve(&problem, s.arrangement).unwrap();
        let dt1 = s.hot_in_c - r.cold_out;
        let dt2 = r.hot_out - s.cold_in_c;
        let expected_lmtd = (dt1 - dt2) / (dt1 / dt2).ln();
        assert!(
            (expected_lmtd - 22.289_126_442_5).abs() < 1e-6,
            "expected LMTD ≈ 22.2891 °C, got {expected_lmtd}"
        );
        // Q = UA · LMTD must hold to round-off.
        assert!(
            (r.q_w - problem.ua_w_per_k * expected_lmtd).abs() < 1e-6,
            "Q = UA·LMTD identity broken: Q = {}, UA·LMTD = {}",
            r.q_w,
            problem.ua_w_per_k * expected_lmtd
        );
    }

    #[test]
    fn counterflow_is_at_least_as_effective_as_parallel() {
        // Ground-truth physics: for 0 < Cr <= 1 a counterflow exchanger
        // is never less effective than the parallel-flow one at equal NTU.
        let problem = NtuProblem::new(2000.0, 3000.0, 5000.0, 90.0, 20.0).unwrap();
        let cf = solve(&problem, FlowArrangement::Counterflow).unwrap();
        let pf = solve(&problem, FlowArrangement::ParallelFlow).unwrap();
        assert!(
            cf.effectiveness >= pf.effectiveness - 1e-12,
            "counterflow ε {} should be >= parallel ε {}",
            cf.effectiveness,
            pf.effectiveness
        );
    }

    #[test]
    fn exchanger_mesh_for_default_is_nonempty_and_in_range() {
        let s = HeatExchangerWorkbenchState::default();
        let mesh = exchanger_solid_mesh(&s).expect("default exchanger yields a solid");
        assert!(
            mesh.nodes.len() > 8,
            "expected shell + tubes + heads + base"
        );
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn exchanger_mesh_none_for_invalid() {
        let s = HeatExchangerWorkbenchState {
            hot_in_c: 20.0,
            cold_in_c: 90.0,
            ..Default::default()
        };
        assert!(exchanger_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_heatexchanger_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_heatexchanger_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_heatexchanger_workbench = true;
        run_exchanger(&mut app.heatexchanger);
        draw_workbench(&mut app);
    }
}
