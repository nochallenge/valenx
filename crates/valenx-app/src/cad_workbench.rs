//! Parametric-CAD workbench — named parameters drive sketch geometry.
//!
//! A right-side panel over `valenx-solvespace-3d`: define **named parameters**
//! with expressions (Fusion's "Change Parameters"), pick one to drive a
//! circle's radius, and Solve. The panel resolves every parameter (showing its
//! value or an error) and the constraint solver lands the circle on the
//! parameter-driven radius — the parameter → expression → constraint-target →
//! geometry loop, in the UI. Compute is synchronous (sub-millisecond).

use eframe::egui;

use valenx_solvespace_3d::{Constraint3D, ParameterTable, Sketch3D};

use crate::ValenxApp;

/// Persistent state for the parametric-CAD workbench.
pub struct CadWorkbenchState {
    /// Editable named parameters as (name, expression) rows.
    params: Vec<(String, String)>,
    /// Name of the parameter that drives the circle's radius.
    radius_param: String,
    results: Option<CadResults>,
}

impl Default for CadWorkbenchState {
    fn default() -> Self {
        Self {
            params: vec![
                ("base".to_string(), "4".to_string()),
                ("radius".to_string(), "base + 1".to_string()),
            ],
            radius_param: "radius".to_string(),
            results: None,
        }
    }
}

struct CadResults {
    /// Each parameter's resolved value or error message.
    resolved: Vec<(String, Result<f64, String>)>,
    /// The solved circle radius, if the sketch solved.
    solved_radius: Option<f64>,
    /// Solver / status message.
    status: String,
}

/// Resolve the parameters and solve a circle whose radius is driven by the
/// chosen parameter.
fn run_cad(s: &CadWorkbenchState) -> CadResults {
    let mut table = ParameterTable::new();
    for (n, e) in &s.params {
        let n = n.trim();
        if !n.is_empty() {
            table.set(n, e);
        }
    }
    let resolved: Vec<(String, Result<f64, String>)> = s
        .params
        .iter()
        .filter(|(n, _)| !n.trim().is_empty())
        .map(|(n, _)| (n.clone(), table.value(n.trim()).map_err(|e| e.to_string())))
        .collect();

    let (solved_radius, status) = match table.value(s.radius_param.trim()) {
        Ok(r) => {
            let mut sk = Sketch3D::new();
            let c = sk.add_point(0.0, 0.0, 0.0);
            let circle = sk.add_circle(c, 1.0, 0.0, 0.0, 1.0).expect("centre is a point");
            sk.add_constraint(Constraint3D::CircleRadius { circle, target: r });
            match sk.solve() {
                Ok(rep) => {
                    let solved = sk.circle_radius(circle);
                    (Some(solved), format!("{:?} — circle radius = {solved:.4}", rep.status))
                }
                Err(e) => (None, format!("solve error: {e}")),
            }
        }
        Err(e) => (None, format!("radius parameter '{}': {e}", s.radius_param.trim())),
    };

    CadResults { resolved, solved_radius, status }
}

/// Draw the parametric-CAD workbench (a no-op unless toggled on via
/// View → Parametric CAD).
pub fn draw_cad_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_cad_workbench {
        return;
    }
    egui::SidePanel::right("valenx_cad_workbench")
        .resizable(true)
        .default_width(340.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("Parametric CAD");
            ui.label(
                egui::RichText::new("named parameters · valenx-solvespace-3d")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.cad;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Parameters (name = expression)").strong());
                    let mut remove: Option<usize> = None;
                    for (i, (name, expr)) in s.params.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.add(
                                egui::TextEdit::singleline(name)
                                    .desired_width(80.0)
                                    .hint_text("name"),
                            );
                            ui.label("=");
                            ui.add(
                                egui::TextEdit::singleline(expr)
                                    .desired_width(130.0)
                                    .hint_text("expr"),
                            );
                            if ui.small_button("✕").clicked() {
                                remove = Some(i);
                            }
                        });
                    }
                    if let Some(i) = remove {
                        s.params.remove(i);
                    }
                    if ui.button("+ parameter").clicked() {
                        s.params.push((String::new(), String::new()));
                    }
                    ui.add_space(4.0);
                    ui.horizontal(|ui| {
                        ui.label("circle radius =");
                        ui.add(
                            egui::TextEdit::singleline(&mut s.radius_param)
                                .desired_width(100.0)
                                .hint_text("parameter"),
                        );
                    });
                    ui.separator();
                    if ui.button("▶ Solve").clicked() {
                        let res = run_cad(s);
                        s.results = Some(res);
                    }
                    if let Some(r) = &s.results {
                        ui.separator();
                        ui.label(egui::RichText::new("Resolved").strong());
                        for (name, val) in &r.resolved {
                            match val {
                                Ok(v) => ui.label(
                                    egui::RichText::new(format!("{name} = {v:.4}"))
                                        .monospace()
                                        .small(),
                                ),
                                Err(e) => ui.colored_label(
                                    egui::Color32::from_rgb(220, 120, 80),
                                    egui::RichText::new(format!("{name}: {e}")).small(),
                                ),
                            };
                        }
                        ui.add_space(4.0);
                        let color = if r.solved_radius.is_some() {
                            egui::Color32::from_rgb(80, 220, 120)
                        } else {
                            egui::Color32::from_rgb(220, 120, 80)
                        };
                        ui.colored_label(color, &r.status);
                    }
                });
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parameter_drives_circle_radius() {
        // Defaults: base = 4, radius = base + 1, circle radius driven by `radius`.
        let s = CadWorkbenchState::default();
        let r = run_cad(&s);
        let solved = r.solved_radius.expect("sketch solved");
        assert!((solved - 5.0).abs() < 1e-4, "radius = {solved}");
        let radius = r.resolved.iter().find(|(n, _)| n == "radius").unwrap();
        assert_eq!(radius.1.as_ref().ok().map(|v| (v * 1e4).round() / 1e4), Some(5.0));
    }

    #[test]
    fn cyclic_parameters_report_an_error_not_a_panic() {
        let s = CadWorkbenchState {
            params: vec![
                ("a".to_string(), "b + 1".to_string()),
                ("b".to_string(), "a + 1".to_string()),
            ],
            radius_param: "a".to_string(),
            results: None,
        };
        let r = run_cad(&s);
        assert!(r.solved_radius.is_none());
        assert!(r.status.to_lowercase().contains("cyclic"), "status: {}", r.status);
    }
}
