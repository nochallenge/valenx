//! The right-side **Population Dynamics Workbench** panel — native
//! closed-form epidemiology / ecology over `valenx-popdynamics`.
//!
//! Mirrors the Heat Transfer / Pharmacokinetics workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_popdynamics_workbench`,
//! toggled from the View menu. A model selector picks one of three textbook
//! models — SIR epidemic, logistic growth, or Lotka-Volterra predator-prey —
//! and only that model's inputs are shown. "Analyze" reports the model's
//! closed-form diagnostics (SIR's `R0` and herd-immunity threshold; the
//! logistic population at a chosen time and its carrying capacity; the
//! Lotka-Volterra coexistence equilibrium and conserved quantity), and
//! "Show 3-D" loads a representative trajectory curve as a swept ribbon into
//! the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_popdynamics::{Logistic, LotkaVolterra, Sir, SirState};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which population-dynamics model the workbench evaluates.
#[derive(Debug, Clone, Copy, PartialEq)]
enum PopModel {
    /// Kermack-McKendrick SIR epidemic (`R0`, herd-immunity threshold).
    Sir,
    /// Single-species logistic (Verhulst) growth (closed-form `N(t)`).
    Logistic,
    /// Two-species Lotka-Volterra predator-prey (equilibrium, conserved `H`).
    LotkaVolterra,
}

/// Persistent form + result state for the Population Dynamics Workbench.
pub struct PopDynamicsWorkbenchState {
    /// Which model is active.
    model: PopModel,

    // --- SIR -----------------------------------------------------------
    /// SIR transmission rate `beta` (per day).
    sir_beta: f64,
    /// SIR recovery rate `gamma` (per day); `1/gamma` is the infectious period.
    sir_gamma: f64,
    /// Initial infectious seed `I0` (people), for the representative curve.
    sir_i0: f64,
    /// Total population `N` (people), for the representative curve.
    sir_n: f64,

    // --- Logistic ------------------------------------------------------
    /// Logistic intrinsic per-capita growth rate `r` (per unit time).
    log_r: f64,
    /// Logistic carrying capacity `K`.
    log_k: f64,
    /// Logistic initial population `N0`.
    log_n0: f64,
    /// Time `t` at which to report the closed-form `N(t)`.
    log_t: f64,

    // --- Lotka-Volterra ------------------------------------------------
    /// Prey birth rate `alpha`.
    lv_alpha: f64,
    /// Predation rate `beta`.
    lv_beta: f64,
    /// Predator death rate `gamma`.
    lv_gamma: f64,
    /// Predator reproduction-per-prey rate `delta`.
    lv_delta: f64,
    /// Initial prey population.
    lv_prey0: f64,
    /// Initial predator population.
    lv_pred0: f64,

    /// Formatted diagnostics readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D curve ribbon (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for PopDynamicsWorkbenchState {
    fn default() -> Self {
        // Realistic textbook defaults.
        //   SIR:  beta=0.5, gamma=0.2/day -> R0 = 2.5, HIT = 0.60.
        //   Logistic: r=0.7, K=1000, N0=10; N(10) ~ 911.
        //   Lotka-Volterra: classic Hare/Lynx-style rates.
        Self {
            model: PopModel::Sir,

            sir_beta: 0.5,
            sir_gamma: 0.2,
            sir_i0: 1.0,
            sir_n: 1000.0,

            log_r: 0.7,
            log_k: 1000.0,
            log_n0: 10.0,
            log_t: 10.0,

            lv_alpha: 1.0,
            lv_beta: 0.1,
            lv_gamma: 1.5,
            lv_delta: 0.075,
            lv_prey0: 10.0,
            lv_pred0: 5.0,

            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Population Dynamics Workbench right-side panel. A no-op when the
/// `show_popdynamics_workbench` toggle is off.
pub fn draw_popdynamics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_popdynamics_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_popdynamics_workbench",
        "Population Dynamics",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native closed-form epidemiology / ecology · valenx-popdynamics",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.popdynamics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Model").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.model, PopModel::Sir, "SIR");
                        ui.radio_value(&mut s.model, PopModel::Logistic, "Logistic");
                        ui.radio_value(
                            &mut s.model,
                            PopModel::LotkaVolterra,
                            "Lotka-Volterra",
                        );
                    });

                    ui.add_space(4.0);
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as
                    // its accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name, leaving it anonymous to a screen
                    // reader / AI driver otherwise); the hover text mirrors the
                    // caption for a mouse user. Applied across every model arm.
                    match s.model {
                        PopModel::Sir => {
                            ui.label(egui::RichText::new("SIR epidemic").strong());
                            ui.horizontal(|ui| {
                                let lbl = ui.label("transmission β (/day)");
                                ui.add(egui::DragValue::new(&mut s.sir_beta).speed(0.01))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("transmission β (/day)");
                            });
                            ui.horizontal(|ui| {
                                let lbl = ui.label("recovery γ (/day)");
                                ui.add(egui::DragValue::new(&mut s.sir_gamma).speed(0.01))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("recovery γ (/day)");
                            });
                            ui.horizontal(|ui| {
                                let lbl = ui.label("population N");
                                ui.add(egui::DragValue::new(&mut s.sir_n).speed(10.0))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("population N");
                            });
                            ui.horizontal(|ui| {
                                let lbl = ui.label("initial infectious I₀");
                                ui.add(egui::DragValue::new(&mut s.sir_i0).speed(1.0))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("initial infectious I₀");
                            });
                        }
                        PopModel::Logistic => {
                            ui.label(egui::RichText::new("Logistic growth").strong());
                            ui.horizontal(|ui| {
                                let lbl = ui.label("growth rate r");
                                ui.add(egui::DragValue::new(&mut s.log_r).speed(0.01))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("growth rate r");
                            });
                            ui.horizontal(|ui| {
                                let lbl = ui.label("carrying capacity K");
                                ui.add(egui::DragValue::new(&mut s.log_k).speed(10.0))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("carrying capacity K");
                            });
                            ui.horizontal(|ui| {
                                let lbl = ui.label("initial N₀");
                                ui.add(egui::DragValue::new(&mut s.log_n0).speed(1.0))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("initial N₀");
                            });
                            ui.horizontal(|ui| {
                                let lbl = ui.label("report at time t");
                                ui.add(egui::DragValue::new(&mut s.log_t).speed(0.5))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("report at time t");
                            });
                        }
                        PopModel::LotkaVolterra => {
                            ui.label(egui::RichText::new("Predator-prey").strong());
                            ui.horizontal(|ui| {
                                let lbl = ui.label("prey birth α");
                                ui.add(egui::DragValue::new(&mut s.lv_alpha).speed(0.01))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("prey birth α");
                            });
                            ui.horizontal(|ui| {
                                let lbl = ui.label("predation β");
                                ui.add(egui::DragValue::new(&mut s.lv_beta).speed(0.005))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("predation β");
                            });
                            ui.horizontal(|ui| {
                                let lbl = ui.label("predator death γ");
                                ui.add(egui::DragValue::new(&mut s.lv_gamma).speed(0.01))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("predator death γ");
                            });
                            ui.horizontal(|ui| {
                                let lbl = ui.label("predator growth δ");
                                ui.add(egui::DragValue::new(&mut s.lv_delta).speed(0.005))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("predator growth δ");
                            });
                            ui.horizontal(|ui| {
                                let lbl = ui.label("initial prey");
                                ui.add(egui::DragValue::new(&mut s.lv_prey0).speed(1.0))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("initial prey");
                            });
                            ui.horizontal(|ui| {
                                let lbl = ui.label("initial predator");
                                ui.add(egui::DragValue::new(&mut s.lv_pred0).speed(1.0))
                                    .labelled_by(lbl.id)
                                    .on_hover_text("initial predator");
                            });
                        }
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_popdynamics(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative trajectory curve (SIR infectious / logistic N(t) / Lotka-Volterra prey) as a swept 3-D ribbon and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Diagnostics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_popdynamics_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.popdynamics` borrow is
    // released here): build the trajectory ribbon and load it.
    if app.popdynamics.show_3d_request {
        app.popdynamics.show_3d_request = false;
        load_curve_3d(app);
    }
}

/// Validate the form, evaluate the active model and format the readout.
fn run_popdynamics(s: &mut PopDynamicsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the active model's closed-form diagnostics and format the
/// readout, mapping any domain error to a display string. Extracted so it
/// is unit-testable.
fn compute(s: &PopDynamicsWorkbenchState) -> Result<String, String> {
    match s.model {
        PopModel::Sir => {
            let m = Sir::new(s.sir_beta, s.sir_gamma).map_err(|e| e.to_string())?;
            let r0 = m.r0();
            let hit = m.herd_immunity_threshold();
            let grows = if r0 > 1.0 { "yes" } else { "no" };
            Ok(format!(
                "model           : SIR epidemic\n\
                 transmission β  : {beta:.3} /day\n\
                 recovery γ      : {gamma:.3} /day\n\
                 infectious period: {period:.2} day\n\n\
                 R0 = β/γ        : {r0:.3}\n\
                 herd immunity   : {hit_pct:.1} %  (1 − 1/R0)\n\
                 epidemic grows  : {grows}  (R0 > 1)",
                beta = s.sir_beta,
                gamma = s.sir_gamma,
                period = 1.0 / s.sir_gamma,
                hit_pct = hit * 100.0,
            ))
        }
        PopModel::Logistic => {
            let m = Logistic::new(s.log_r, s.log_k).map_err(|e| e.to_string())?;
            if !s.log_n0.is_finite() || s.log_n0 < 0.0 {
                return Err("initial population N₀ must be finite and non-negative".into());
            }
            let n_t = m.analytic(s.log_n0, s.log_t);
            Ok(format!(
                "model           : logistic growth\n\
                 growth rate r   : {r:.3}\n\
                 carrying cap. K : {k:.1}\n\
                 initial N₀      : {n0:.1}\n\
                 report time t   : {t:.2}\n\n\
                 N(t)            : {n_t:.2}\n\
                 fraction of K   : {frac_pct:.1} %\n\
                 N(∞) → K        : {k:.1}",
                r = s.log_r,
                k = s.log_k,
                n0 = s.log_n0,
                t = s.log_t,
                frac_pct = (n_t / s.log_k) * 100.0,
            ))
        }
        PopModel::LotkaVolterra => {
            let m = LotkaVolterra::new(s.lv_alpha, s.lv_beta, s.lv_gamma, s.lv_delta)
                .map_err(|e| e.to_string())?;
            let eq = m.equilibrium();
            let state = [s.lv_prey0, s.lv_pred0];
            let h = m.conserved_quantity(&state);
            Ok(format!(
                "model           : Lotka-Volterra\n\
                 prey birth α    : {alpha:.3}\n\
                 predation β     : {beta:.3}\n\
                 predator death γ: {gamma:.3}\n\
                 predator grow δ : {delta:.3}\n\n\
                 equilibrium x*  : {prey_eq:.3}  (γ/δ)\n\
                 equilibrium y*  : {pred_eq:.3}  (α/β)\n\
                 conserved H     : {h:.4}",
                alpha = s.lv_alpha,
                beta = s.lv_beta,
                gamma = s.lv_gamma,
                delta = s.lv_delta,
                prey_eq = eq[0],
                pred_eq = eq[1],
            ))
        }
    }
}

/// Sample the active model's representative trajectory as `(t, value)`
/// pairs over a fixed horizon: SIR infectious count, logistic `N(t)`, or
/// Lotka-Volterra prey count. `None` for an invalid configuration. The
/// quantities all come from `valenx-popdynamics` (closed-form for the
/// logistic, RK4 `simulate` for SIR / Lotka-Volterra).
fn curve_samples(s: &PopDynamicsWorkbenchState) -> Option<Vec<(f64, f64)>> {
    match s.model {
        PopModel::Sir => {
            let m = Sir::new(s.sir_beta, s.sir_gamma).ok()?;
            let i0 = s.sir_i0.clamp(0.0, s.sir_n.max(0.0));
            let init = SirState::new((s.sir_n - i0).max(0.0), i0, 0.0);
            let traj = m.simulate(init, 60.0, 0.2).ok()?;
            Some(traj.iter().map(|p| (p.t, p.y[1])).collect())
        }
        PopModel::Logistic => {
            let m = Logistic::new(s.log_r, s.log_k).ok()?;
            if !s.log_n0.is_finite() || s.log_n0 < 0.0 {
                return None;
            }
            let t_end = s.log_t.max(1.0) * 1.5;
            let steps = 240;
            Some(
                (0..=steps)
                    .map(|j| {
                        let t = t_end * j as f64 / steps as f64;
                        (t, m.analytic(s.log_n0, t))
                    })
                    .collect(),
            )
        }
        PopModel::LotkaVolterra => {
            let m = LotkaVolterra::new(s.lv_alpha, s.lv_beta, s.lv_gamma, s.lv_delta).ok()?;
            let traj = m.simulate([s.lv_prey0, s.lv_pred0], 30.0, 0.05).ok()?;
            Some(traj.iter().map(|p| (p.t, p.y[0])).collect())
        }
    }
}

/// Append one thin swept-ribbon quad strip following the sampled curve to
/// the buffers, normalising time to `x in [0, span_x]` and value to
/// `z in [0, span_z]`, extruded a half-width `hw` in `+/- y`. Two triangles
/// per segment.
fn push_curve_ribbon(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    samples: &[(f64, f64)],
    span_x: f64,
    span_z: f64,
    hw: f64,
) {
    if samples.len() < 2 {
        return;
    }
    let t_max = samples.last().map_or(1.0, |p| p.0).max(1e-9);
    let v_max = samples
        .iter()
        .map(|p| p.1)
        .fold(f64::MIN, f64::max)
        .max(1e-9);

    let base = nodes.len();
    for &(t, v) in samples {
        let x = (t / t_max) * span_x;
        let z = (v / v_max) * span_z;
        nodes.push(Vector3::new(x, -hw, z));
        nodes.push(Vector3::new(x, hw, z));
    }
    for j in 0..samples.len() - 1 {
        let a = base + 2 * j;
        let b = a + 1;
        let c = a + 2;
        let d = a + 3;
        // Two triangles (both windings so the thin strip is visible from
        // either side).
        tris.extend_from_slice(&[a, c, d, a, d, b, a, d, c, a, b, d]);
    }
}

/// Build the active model's representative trajectory as a triangle [`Mesh`]
/// — the curve swept into a thin ribbon, with a thin base strip along the
/// time axis. Representative geometry (not to scale; the diagnostics numbers
/// are the `valenx-popdynamics` result). `None` for an invalid configuration.
fn curve_ribbon_mesh(s: &PopDynamicsWorkbenchState) -> Option<Mesh> {
    let samples = curve_samples(s)?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let span_x = 2.0;
    let span_z = 1.0;
    push_curve_ribbon(&mut nodes, &mut tris, &samples, span_x, span_z, 0.04);
    // Base strip along the time axis (a flat ground line for the curve).
    push_curve_ribbon(
        &mut nodes,
        &mut tris,
        &[(0.0, 0.0), (1.0, 0.0)],
        span_x,
        span_z,
        0.06,
    );

    if tris.is_empty() {
        return None;
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-popdynamics");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D trajectory ribbon and load it into the central viewport.
fn load_curve_3d(app: &mut ValenxApp) {
    let Some(mesh) = curve_ribbon_mesh(&app.popdynamics) else {
        app.popdynamics.error =
            Some("model parameters are invalid — cannot build the 3-D curve".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<curve>/valenx-popdynamics"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The active model's representative trajectory as a 2-D **line-chart** series —
/// population vs time. The plotted quantity matches [`curve_samples`]: SIR
/// infectious count, logistic `N(t)`, or Lotka-Volterra prey count, all the
/// genuine `valenx-popdynamics` result (closed-form logistic / RK4 `simulate`
/// for the others). `None` for an invalid configuration.
fn trajectory_chart(s: &PopDynamicsWorkbenchState) -> Option<crate::ChartData> {
    let samples = curve_samples(s)?;
    let (label, y_label) = match s.model {
        PopModel::Sir => ("infectious I(t)", "population"),
        PopModel::Logistic => ("N(t)", "population N"),
        PopModel::LotkaVolterra => ("prey x(t)", "population"),
    };
    let points: Vec<[f64; 2]> = samples.iter().map(|&(t, v)| [t, v]).collect();
    Some(crate::ChartData {
        title: "Population trajectory".into(),
        x_label: "time".into(),
        y_label: y_label.into(),
        series: vec![crate::ChartSeries {
            label: label.into(),
            points,
            bars: false,
        }],
    })
}

/// The agent-bridge **`show_3d{kind:"popdynamics"}`** product: the canonical
/// model trajectory presented as a 2-D **line chart** (population vs time — see
/// [`trajectory_chart`]) paired with the workbench's own model-diagnostics
/// headline numbers. A population-vs-time trajectory reads far better as a
/// framed, auto-scaled curve than as a swept 3-D ribbon, so this product carries
/// `kind2d: Some(Chart(..))` and no mesh. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`PopDynamicsWorkbenchState::default`] (the default SIR
/// model). The readout rows mirror the panel's `compute()` readout.
pub(crate) fn popdynamics_product() -> crate::WorkspaceProduct {
    let s = PopDynamicsWorkbenchState::default();
    let chart = trajectory_chart(&s).expect("default SIR model ⇒ a trajectory chart");
    let readout = compute(&s).expect("default SIR model ⇒ a valid readout");
    let lines = crate::products_registry::lines_from_readout(&readout);
    crate::WorkspaceProduct {
        title: "Population Dynamics".into(),
        lines,
        mesh: None,
        vertex_colors: None,
        camera: valenx_viz::OrbitCamera::default(),
        kind2d: Some(crate::Workspace2dKind::Chart(chart)),
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
        let s = PopDynamicsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
        assert_eq!(s.model, PopModel::Sir);
    }

    #[test]
    fn analyze_sir_reports_r0_and_herd_immunity() {
        let mut s = PopDynamicsWorkbenchState::default();
        run_popdynamics(&mut s);
        assert!(
            s.error.is_none(),
            "default SIR should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("R0 = β/γ"));
        assert!(s.result.contains("herd immunity"));
        // beta=0.5, gamma=0.2 -> R0 = 2.5 exactly; HIT = 1 - 1/2.5 = 0.60.
        assert!(s.result.contains("2.500"));
        assert!(s.result.contains("60.0 %"));
        assert!(s.result.contains("epidemic grows  : yes"));
    }

    #[test]
    fn sir_ground_truth_r0_and_hit_hand_computed() {
        // GROUND TRUTH (hand-computed): R0 = beta/gamma and the
        // herd-immunity threshold HIT = 1 - 1/R0.
        let m = Sir::new(0.5, 0.2).unwrap();
        assert!((m.r0() - 2.5).abs() < 1e-12, "R0={}", m.r0());
        // 1 - 1/2.5 = 1 - 0.4 = 0.6.
        assert!((m.herd_immunity_threshold() - 0.6).abs() < 1e-12);
    }

    #[test]
    fn analyze_logistic_reports_endpoints() {
        let mut s = PopDynamicsWorkbenchState {
            model: PopModel::Logistic,
            ..Default::default()
        };
        run_popdynamics(&mut s);
        assert!(s.error.is_none(), "logistic should analyze: {:?}", s.error);
        assert!(s.result.contains("N(t)"));
        assert!(s.result.contains("N(∞) → K"));
        // GROUND TRUTH endpoints of the closed form N(t)=K/(1+((K-N0)/N0)e^-rt):
        //   N(0) = N0 = 10, and N(t) -> K = 1000 as t -> infinity.
        let m = Logistic::new(s.log_r, s.log_k).unwrap();
        assert!((m.analytic(s.log_n0, 0.0) - s.log_n0).abs() < 1e-9);
        // `analytic` uses the e^{r t} form, which overflows to NaN for very
        // large t; t = 50 is already well converged to K (within ~1e-10).
        assert!((m.analytic(s.log_n0, 50.0) - s.log_k).abs() < 1e-3);
    }

    #[test]
    fn analyze_lotka_volterra_reports_equilibrium_and_conserved() {
        let mut s = PopDynamicsWorkbenchState {
            model: PopModel::LotkaVolterra,
            ..Default::default()
        };
        run_popdynamics(&mut s);
        assert!(
            s.error.is_none(),
            "Lotka-Volterra should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("equilibrium x*"));
        assert!(s.result.contains("conserved H"));
        // x* = gamma/delta = 1.5/0.075 = 20, y* = alpha/beta = 1.0/0.1 = 10.
        let m = LotkaVolterra::new(s.lv_alpha, s.lv_beta, s.lv_gamma, s.lv_delta).unwrap();
        let eq = m.equilibrium();
        assert!((eq[0] - 20.0).abs() < 1e-9, "x*={}", eq[0]);
        assert!((eq[1] - 10.0).abs() < 1e-9, "y*={}", eq[1]);
    }

    #[test]
    fn analyze_rejects_bad_parameters_per_model() {
        // SIR: gamma must be strictly positive.
        let mut s = PopDynamicsWorkbenchState {
            sir_gamma: 0.0,
            ..Default::default()
        };
        run_popdynamics(&mut s);
        assert!(s.error.is_some());

        // Logistic: K must be strictly positive.
        let mut s = PopDynamicsWorkbenchState {
            model: PopModel::Logistic,
            log_k: 0.0,
            ..Default::default()
        };
        run_popdynamics(&mut s);
        assert!(s.error.is_some());

        // Lotka-Volterra: all four rates must be strictly positive.
        let mut s = PopDynamicsWorkbenchState {
            model: PopModel::LotkaVolterra,
            lv_alpha: 0.0,
            ..Default::default()
        };
        run_popdynamics(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn curve_mesh_for_each_model_is_nonempty_and_in_range() {
        for model in [PopModel::Sir, PopModel::Logistic, PopModel::LotkaVolterra] {
            let s = PopDynamicsWorkbenchState {
                model,
                ..Default::default()
            };
            let mesh = curve_ribbon_mesh(&s).expect("default model yields a ribbon");
            assert!(
                mesh.nodes.len() > 4,
                "expected a swept ribbon for {model:?}"
            );
            let n = mesh.nodes.len() as u32;
            for blk in &mesh.element_blocks {
                assert!(!blk.connectivity.is_empty());
                assert_eq!(blk.connectivity.len() % 3, 0);
                assert!(blk.connectivity.iter().all(|&i| i < n));
            }
        }
    }

    #[test]
    fn curve_mesh_none_for_invalid() {
        let s = PopDynamicsWorkbenchState {
            sir_gamma: 0.0,
            ..Default::default()
        };
        assert!(curve_ribbon_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_popdynamics_workbench(app, ctx);
        });
    }

    /// As `draw_workbench`, but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_popdynamics_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_popdynamics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_popdynamics_workbench = true;
        run_popdynamics(&mut app.popdynamics);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The default model is SIR, which shows four DragValues; each is a
        // SpinButton that must be `labelled_by` its caption (egui clears a
        // DragValue's own Name), so an AI / screen reader can find the control
        // by the caption text. Only the SIR fields are counted — the Logistic
        // and Lotka-Volterra fields are hidden behind a non-default selection.
        let mut app = ValenxApp::default();
        app.show_popdynamics_workbench = true;
        assert_eq!(app.popdynamics.model, PopModel::Sir, "default model is SIR");
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // SIR shows transmission β, recovery γ, population N, initial I₀.
        assert!(
            spin_buttons.len() >= 4,
            "expected the SIR numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every SIR DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["transmission β (/day)", "recovery γ (/day)", "population N"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        // The Analyze button stays a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::Button
                && n.name().is_some_and(|s| s.contains("Analyze"))),
            "the Analyze button is a named, invokable node"
        );
    }
}
