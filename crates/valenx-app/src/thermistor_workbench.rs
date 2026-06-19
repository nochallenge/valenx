//! The right-side **Thermistor Workbench** panel — native NTC
//! resistance-temperature analysis over `valenx-thermistor`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_thermistor_workbench`,
//! toggled from the View menu. The form picks one of two resistance-
//! temperature laws — the single-parameter beta model
//! `R = R0 * exp(beta (1/T - 1/T0))` or the three-coefficient
//! Steinhart-Hart model `1/T = A + B ln(R) + C ln(R)^3` — and a query
//! temperature and resistance. "Analyze" converts in both directions
//! (resistance at the query temperature, temperature at the query
//! resistance) and reports the local temperature coefficient of
//! resistance, and "Show 3-D" loads a representative thermistor bead
//! solid into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_thermistor::units::{celsius_to_kelvin, kelvin_to_celsius};
use valenx_thermistor::{BetaModel, SteinhartHart};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which resistance-temperature law the workbench evaluates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ThermistorModel {
    /// Single-parameter beta (`B`) law `R = R0 exp(beta (1/T - 1/T0))`.
    Beta,
    /// Three-coefficient Steinhart-Hart law `1/T = A + B ln R + C ln(R)^3`.
    Steinhart,
}

/// Persistent form + result state for the Thermistor Workbench.
pub struct ThermistorWorkbenchState {
    /// Which resistance-temperature law to evaluate.
    model: ThermistorModel,
    /// Beta model: reference resistance `R0` (ohms).
    r0_ohms: f64,
    /// Beta model: reference temperature `T0` (deg C).
    t0_c: f64,
    /// Beta model: material constant `beta` (`B`), in kelvin.
    beta_kelvin: f64,
    /// Steinhart-Hart `A` coefficient (1/K).
    sh_a: f64,
    /// Steinhart-Hart `B` coefficient (1/K).
    sh_b: f64,
    /// Steinhart-Hart `C` coefficient (1/K).
    sh_c: f64,
    /// Query temperature for resistance-at-T (deg C).
    query_t_c: f64,
    /// Query resistance for temperature-at-R (ohms).
    query_r_ohms: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D bead solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for ThermistorWorkbenchState {
    fn default() -> Self {
        // A 10 kohm-at-25C NTC with beta = 3950 K (the most common
        // hobbyist/industrial part). The Steinhart-Hart defaults are a
        // widely tabulated Vishay-style coefficient set for a ~10 kohm
        // part. Query at 50 C and 20 kohm.
        Self {
            model: ThermistorModel::Beta,
            r0_ohms: 10_000.0,
            t0_c: 25.0,
            beta_kelvin: 3950.0,
            sh_a: 1.009_249_522e-3,
            sh_b: 2.378_405_444e-4,
            sh_c: 2.019_202_697e-7,
            query_t_c: 50.0,
            query_r_ohms: 20_000.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Thermistor Workbench right-side panel. A no-op when the
/// `show_thermistor_workbench` toggle is off.
pub fn draw_thermistor_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_thermistor_workbench {
        return;
    }

    egui::SidePanel::right("valenx_thermistor_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Thermistor",
                "native NTC resistance-temperature models · valenx-thermistor",
            ) {
                app.show_thermistor_workbench = false;
            }

            let s = &mut app.thermistor;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Model").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.model, ThermistorModel::Beta, "Beta (B)");
                        ui.radio_value(
                            &mut s.model,
                            ThermistorModel::Steinhart,
                            "Steinhart-Hart",
                        );
                    });

                    ui.add_space(4.0);
                    match s.model {
                        ThermistorModel::Beta => {
                            ui.label(egui::RichText::new("Beta parameters").strong());
                            ui.horizontal(|ui| {
                                ui.label("R0 (Ω)");
                                ui.add(egui::DragValue::new(&mut s.r0_ohms).speed(100.0));
                            });
                            ui.horizontal(|ui| {
                                ui.label("T0 (°C)");
                                ui.add(egui::DragValue::new(&mut s.t0_c).speed(0.5));
                            });
                            ui.horizontal(|ui| {
                                ui.label("beta (K)");
                                ui.add(egui::DragValue::new(&mut s.beta_kelvin).speed(10.0));
                            });
                        }
                        ThermistorModel::Steinhart => {
                            ui.label(egui::RichText::new("Steinhart-Hart coefficients").strong());
                            ui.horizontal(|ui| {
                                ui.label("A (1/K)");
                                ui.add(
                                    egui::DragValue::new(&mut s.sh_a).speed(1e-5).max_decimals(12),
                                );
                            });
                            ui.horizontal(|ui| {
                                ui.label("B (1/K)");
                                ui.add(
                                    egui::DragValue::new(&mut s.sh_b).speed(1e-6).max_decimals(12),
                                );
                            });
                            ui.horizontal(|ui| {
                                ui.label("C (1/K)");
                                ui.add(
                                    egui::DragValue::new(&mut s.sh_c).speed(1e-9).max_decimals(15),
                                );
                            });
                        }
                    }

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Query").strong());
                    ui.horizontal(|ui| {
                        ui.label("temperature (°C)");
                        ui.add(egui::DragValue::new(&mut s.query_t_c).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("resistance (Ω)");
                        ui.add(egui::DragValue::new(&mut s.query_r_ohms).speed(100.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_thermistor(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative thermistor bead (a glass-bead body on two lead wires) as a 3-D solid and load it into the central viewport to orbit",
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
        });

    // Serviced after the panel draws (the `&mut app.thermistor` borrow is
    // released here): build the bead's 3-D solid and load it.
    if app.thermistor.show_3d_request {
        app.thermistor.show_3d_request = false;
        load_bead_3d(app);
    }
}

/// Validate the form, evaluate the model and format the readout.
fn run_thermistor(s: &mut ThermistorWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the selected model at the query temperature and resistance —
/// resistance-at-T, temperature-at-R, and the local temperature
/// coefficient of resistance — returning `(r_at_t, t_at_r_kelvin,
/// alpha_per_k)`. Extracted so it is unit-testable and shared with the
/// readout and the 3-D gate. Maps any domain error to a display string.
fn evaluate(s: &ThermistorWorkbenchState) -> Result<(f64, f64, f64), String> {
    let t_query_k = celsius_to_kelvin(s.query_t_c);
    match s.model {
        ThermistorModel::Beta => {
            let m = BetaModel::new(s.r0_ohms, celsius_to_kelvin(s.t0_c), s.beta_kelvin)
                .map_err(|e| e.to_string())?;
            let r_at_t = m.resistance_at(t_query_k).map_err(|e| e.to_string())?;
            let t_at_r = m
                .temperature_at(s.query_r_ohms)
                .map_err(|e| e.to_string())?;
            let alpha = m
                .temperature_coefficient_at(t_query_k)
                .map_err(|e| e.to_string())?;
            Ok((r_at_t, t_at_r, alpha))
        }
        ThermistorModel::Steinhart => {
            let m = SteinhartHart::new(s.sh_a, s.sh_b, s.sh_c).map_err(|e| e.to_string())?;
            let r_at_t = m.resistance_at(t_query_k).map_err(|e| e.to_string())?;
            let t_at_r = m
                .temperature_at(s.query_r_ohms)
                .map_err(|e| e.to_string())?;
            let alpha = m
                .temperature_coefficient_at(t_query_k)
                .map_err(|e| e.to_string())?;
            Ok((r_at_t, t_at_r, alpha))
        }
    }
}

/// Evaluate the model and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &ThermistorWorkbenchState) -> Result<String, String> {
    let (r_at_t, t_at_r_k, alpha) = evaluate(s)?;
    let t_at_r_c = kelvin_to_celsius(t_at_r_k);
    let model_name = match s.model {
        ThermistorModel::Beta => "beta (B)",
        ThermistorModel::Steinhart => "Steinhart-Hart",
    };

    Ok(format!(
        "model           : {model_name}\n\
         query T / R     : {:.2} °C / {:.1} Ω\n\n\
         R at query T    : {:.2} Ω\n\
         T at query R    : {:.2} °C ({:.2} K)\n\
         coeff alpha     : {:.6} /K\n\
         coeff alpha     : {:.4} %/K",
        s.query_t_c,
        s.query_r_ohms,
        r_at_t,
        t_at_r_c,
        t_at_r_k,
        alpha,
        alpha * 100.0,
    ))
}

/// Append an outward-facing capped cylinder (axis along `+z`, centred on
/// `(cx, cy)` and spanning `z0..z1`) of radius `r` with `seg` facets to
/// the buffers, as a fan-capped triangle tube.
#[allow(clippy::too_many_arguments)]
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    cx: f64,
    cy: f64,
    z0: f64,
    z1: f64,
    r: f64,
    seg: usize,
) {
    let base = nodes.len();
    // Two rings of `seg` rim vertices, then two cap centres.
    for &z in &[z0, z1] {
        for k in 0..seg {
            let theta = TAU * (k as f64) / (seg as f64);
            nodes.push(Vector3::new(cx + r * theta.cos(), cy + r * theta.sin(), z));
        }
    }
    let bottom_centre = nodes.len();
    nodes.push(Vector3::new(cx, cy, z0));
    let top_centre = nodes.len();
    nodes.push(Vector3::new(cx, cy, z1));

    for k in 0..seg {
        let next = (k + 1) % seg;
        let b0 = base + k;
        let b1 = base + next;
        let t0 = base + seg + k;
        let t1 = base + seg + next;
        // Side quad as two triangles.
        tris.extend_from_slice(&[b0, b1, t1, b0, t1, t0]);
        // Bottom cap fan.
        tris.extend_from_slice(&[bottom_centre, b1, b0]);
        // Top cap fan.
        tris.extend_from_slice(&[top_centre, t0, t1]);
    }
}

/// Build the thermistor as a triangle [`Mesh`] — an epoxy/glass bead body
/// (a faceted cylinder along `z`) on two thin lead wires below it.
/// Representative geometry (not to scale; the resistance-temperature
/// numbers are the `valenx-thermistor` result). `None` for a model
/// configuration the solver rejects.
fn bead_solid_mesh(s: &ThermistorWorkbenchState) -> Option<Mesh> {
    evaluate(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let seg = 24;
    // Bead body.
    push_cylinder(&mut nodes, &mut tris, 0.0, 0.0, 0.3, 0.9, 0.35, seg);
    // Two lead wires hanging below.
    push_cylinder(&mut nodes, &mut tris, -0.15, 0.0, -0.6, 0.32, 0.04, seg);
    push_cylinder(&mut nodes, &mut tris, 0.15, 0.0, -0.6, 0.32, 0.04, seg);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-thermistor");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D bead solid and load it into the central viewport.
fn load_bead_3d(app: &mut ValenxApp) {
    let Some(mesh) = bead_solid_mesh(&app.thermistor) else {
        app.thermistor.error =
            Some("model parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<bead>/valenx-thermistor"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = ThermistorWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_resistance_temperature_and_coefficient() {
        let mut s = ThermistorWorkbenchState::default();
        run_thermistor(&mut s);
        assert!(
            s.error.is_none(),
            "default beta model should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("R at query T"));
        assert!(s.result.contains("T at query R"));
        assert!(s.result.contains("alpha"));
        // 10k-at-25C, beta=3950 at 50 C: R = 3588.18 ohms.
        assert!(s.result.contains("3588.18"), "got:\n{}", s.result);
    }

    #[test]
    fn analyze_steinhart_default_reports_about_10k_room_resistance() {
        let mut s = ThermistorWorkbenchState {
            model: ThermistorModel::Steinhart,
            query_t_c: 25.0,
            ..Default::default()
        };
        run_thermistor(&mut s);
        assert!(
            s.error.is_none(),
            "default Steinhart-Hart model should analyze: {:?}",
            s.error
        );
        // The tabulated coefficient set is a ~10 kohm part near room temp.
        assert!(s.result.contains("R at query T"));
    }

    #[test]
    fn analyze_rejects_nonpositive_r0() {
        let mut s = ThermistorWorkbenchState {
            r0_ohms: 0.0,
            ..Default::default()
        };
        run_thermistor(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn beta_resistance_at_reference_equals_r0_exactly() {
        // Ground truth: the beta law R = R0 exp(beta (1/T - 1/T0)) has a
        // zero exponent at T == T0, so R(T0) == R0 exactly.
        let s = ThermistorWorkbenchState {
            query_t_c: 25.0, // == T0
            ..Default::default()
        };
        let (r_at_t, _t_at_r, _alpha) = evaluate(&s).expect("default beta evaluates");
        assert!((r_at_t - 10_000.0).abs() < 1e-6, "R(T0) = {r_at_t}");
    }

    #[test]
    fn beta_resistance_matches_hand_computation_at_50c() {
        // Independent hand value: R(323.15) = 10000 exp(3950 (1/323.15 -
        // 1/298.15)) = 3588.182582 ohms.
        let s = ThermistorWorkbenchState::default(); // query at 50 C
        let (r_at_t, _t_at_r, _alpha) = evaluate(&s).expect("default beta evaluates");
        let expected: f64 = 10_000.0 * (3950.0_f64 * (1.0 / 323.15 - 1.0 / 298.15)).exp();
        assert!((r_at_t - expected).abs() < 1e-9, "got {r_at_t}");
        assert!((r_at_t - 3588.182582).abs() < 1e-3, "got {r_at_t}");
    }

    #[test]
    fn bead_mesh_for_default_is_nonempty_and_in_range() {
        let s = ThermistorWorkbenchState::default();
        let mesh = bead_solid_mesh(&s).expect("default model yields a solid");
        assert!(mesh.nodes.len() > 8, "expected bead body + two leads");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn bead_mesh_none_for_invalid() {
        let s = ThermistorWorkbenchState {
            beta_kelvin: 0.0,
            ..Default::default()
        };
        assert!(bead_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_thermistor_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_thermistor_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_thermistor_workbench = true;
        run_thermistor(&mut app.thermistor);
        draw_workbench(&mut app);
    }
}
