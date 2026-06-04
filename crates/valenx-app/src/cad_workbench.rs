//! Parametric-CAD workbench — named parameters drive sketch geometry *and* a
//! CSG feature tree.
//!
//! A right-side panel over `valenx-solvespace-3d`. One shared table of **named
//! parameters** (Fusion's "Change Parameters") feeds two consumers:
//!
//! 1. **Sketch** — pick a parameter to drive a circle's radius and Solve; the
//!    constraint solver lands the circle on the parameter-driven radius.
//! 2. **Feature tree (CSG)** — an ordered list of steps, each placing a
//!    primitive (box / cylinder) and combining it with the running body via
//!    New / Join / Cut / Intersect. Rebuild folds the tree into one solid,
//!    tessellates it, and pushes it into the central 3-D viewport. Edit a
//!    parameter, rebuild, and the whole model re-drives — a hole moves, a boss
//!    grows.
//!
//! Compute is synchronous: parameter resolution is sub-millisecond and a
//! handful of boolean ops on primitives is well under a frame.

use eframe::egui;

use valenx_solvespace_3d::{
    Constraint3D, Feature, FeatureTimeline, Op, ParameterTable, Sketch3D, Step,
};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which primitive a feature-tree step builds.
#[derive(Clone, Copy, PartialEq, Eq)]
enum FeatureKind {
    Box,
    Cylinder,
}

/// One UI-editable feature-tree step. Carries both the box and cylinder
/// dimension fields so toggling `kind` preserves whatever the user typed.
#[derive(Clone)]
struct UiStep {
    op: Op,
    kind: FeatureKind,
    dx: String,
    dy: String,
    dz: String,
    radius: String,
    height: String,
    x: String,
    y: String,
    z: String,
}

impl UiStep {
    fn new_box() -> Self {
        Self {
            op: Op::Join,
            kind: FeatureKind::Box,
            dx: "1".into(),
            dy: "1".into(),
            dz: "1".into(),
            radius: "0.5".into(),
            height: "1".into(),
            x: "0".into(),
            y: "0".into(),
            z: "0".into(),
        }
    }

    fn new_cylinder() -> Self {
        Self {
            op: Op::Cut,
            kind: FeatureKind::Cylinder,
            dx: "1".into(),
            dy: "1".into(),
            dz: "1".into(),
            radius: "0.25".into(),
            height: "2".into(),
            x: "0".into(),
            y: "0".into(),
            z: "0".into(),
        }
    }

    /// Translate into a solver-crate [`Step`].
    fn to_step(&self) -> Step {
        let feature = match self.kind {
            FeatureKind::Box => Feature::Box {
                dx: self.dx.clone(),
                dy: self.dy.clone(),
                dz: self.dz.clone(),
            },
            FeatureKind::Cylinder => Feature::Cylinder {
                radius: self.radius.clone(),
                height: self.height.clone(),
            },
        };
        Step::placed(self.op, feature, self.x.clone(), self.y.clone(), self.z.clone())
    }
}

/// The default feature tree: a unit box with a cylinder punched through it —
/// `valenx-cad`'s proven "punched cube" geometry, so the seed always rebuilds.
fn default_steps() -> Vec<UiStep> {
    vec![
        UiStep {
            op: Op::New,
            kind: FeatureKind::Box,
            dx: "size".into(),
            dy: "size".into(),
            dz: "size".into(),
            radius: "hole_r".into(),
            height: "hole_h".into(),
            x: "0".into(),
            y: "0".into(),
            z: "0".into(),
        },
        UiStep {
            op: Op::Cut,
            kind: FeatureKind::Cylinder,
            dx: "size".into(),
            dy: "size".into(),
            dz: "size".into(),
            radius: "hole_r".into(),
            height: "hole_h".into(),
            x: "size / 2".into(),
            y: "size / 2".into(),
            z: "-0.5".into(),
        },
    ]
}

/// Persistent state for the parametric-CAD workbench.
pub struct CadWorkbenchState {
    /// Editable named parameters as (name, expression) rows — shared by the
    /// sketch demo and the feature tree.
    params: Vec<(String, String)>,
    /// Name of the parameter that drives the circle's radius.
    radius_param: String,
    results: Option<CadResults>,
    /// Feature-tree steps, in build order.
    steps: Vec<UiStep>,
    /// Last rebuild outcome: `Ok(status)` or `Err(message)`.
    tree_status: Option<Result<String, String>>,
    /// Tessellated body waiting to be pushed into the viewport (deferred out
    /// of the panel borrow).
    rebuilt_mesh: Option<valenx_mesh::Mesh>,
    /// Set when a fresh rebuild needs pushing into the viewport.
    push_rebuild: bool,
    /// Snapshots from the last rebuild — `history[k]` is the running body
    /// after step k (the last entry is the final body). Drives the scrubber.
    history: Option<Vec<valenx_cad::Solid>>,
    /// 1-based step the history scrubber is showing (`1..=history.len()`).
    scrub: usize,
}

impl Default for CadWorkbenchState {
    fn default() -> Self {
        Self {
            params: vec![
                ("base".to_string(), "4".to_string()),
                ("radius".to_string(), "base + 1".to_string()),
                ("size".to_string(), "1".to_string()),
                ("hole_r".to_string(), "0.25".to_string()),
                ("hole_h".to_string(), "2".to_string()),
            ],
            radius_param: "radius".to_string(),
            results: None,
            steps: default_steps(),
            tree_status: None,
            rebuilt_mesh: None,
            push_rebuild: false,
            history: None,
            scrub: 1,
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

/// Build a [`ParameterTable`] from the editable rows, skipping blank names.
fn build_table(params: &[(String, String)]) -> ParameterTable {
    let mut table = ParameterTable::new();
    for (n, e) in params {
        let n = n.trim();
        if !n.is_empty() {
            table.set(n, e);
        }
    }
    table
}

/// Resolve the parameters and solve a circle whose radius is driven by the
/// chosen parameter.
fn run_cad(s: &CadWorkbenchState) -> CadResults {
    let table = build_table(&s.params);
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

/// Rebuild the feature tree against the parameters. Returns the per-step
/// snapshot solids (`snapshots[k]` = running body after step k; the last entry
/// is the final body) plus a one-line status, or an error message.
fn rebuild_tree(s: &CadWorkbenchState) -> Result<(Vec<valenx_cad::Solid>, String), String> {
    let table = build_table(&s.params);
    let mut tl = FeatureTimeline::new();
    for st in &s.steps {
        tl.push(st.to_step());
    }
    let model = tl.rebuild(&table).map_err(|e| e.to_string())?;
    let faces = model.body.faces();
    let edges = model.body.edges();
    let verts = model.body.vertices();
    let status = format!(
        "{faces} faces · {edges} edges · {verts} vertices · {} steps",
        s.steps.len()
    );
    Ok((model.snapshots, status))
}

/// Tessellate the running body at a given 1-based step for viewport display.
fn tessellate_step(
    history: &[valenx_cad::Solid],
    step_1based: usize,
) -> Result<valenx_mesh::Mesh, String> {
    let idx = step_1based
        .saturating_sub(1)
        .min(history.len().saturating_sub(1));
    let solid = history.get(idx).ok_or_else(|| "no such step".to_string())?;
    valenx_cad::solid_to_mesh(solid, valenx_cad::DEFAULT_TESS_TOLERANCE).map_err(|e| e.to_string())
}

fn op_label(op: Op) -> &'static str {
    match op {
        Op::New => "New",
        Op::Join => "Join",
        Op::Cut => "Cut",
        Op::Intersect => "Intersect",
    }
}

/// A narrow single-line editor for a dimension / placement expression.
fn dim_edit(ui: &mut egui::Ui, v: &mut String) {
    ui.add(egui::TextEdit::singleline(v).desired_width(52.0));
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

                    // ---- Sketch: parameter-driven circle radius ----
                    ui.add_space(6.0);
                    ui.label(egui::RichText::new("Sketch").strong());
                    ui.horizontal(|ui| {
                        ui.label("circle radius =");
                        ui.add(
                            egui::TextEdit::singleline(&mut s.radius_param)
                                .desired_width(100.0)
                                .hint_text("parameter"),
                        );
                    });
                    if ui.button("▶ Solve").clicked() {
                        let res = run_cad(s);
                        s.results = Some(res);
                    }
                    if let Some(r) = &s.results {
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
                        let color = if r.solved_radius.is_some() {
                            egui::Color32::from_rgb(80, 220, 120)
                        } else {
                            egui::Color32::from_rgb(220, 120, 80)
                        };
                        ui.colored_label(color, &r.status);
                    }

                    // ---- Feature tree (CSG) ----
                    ui.separator();
                    ui.label(egui::RichText::new("Feature tree (CSG)").strong());
                    ui.label(
                        egui::RichText::new(
                            "each step places a primitive and combines it with the running body",
                        )
                        .weak()
                        .small(),
                    );

                    let mut remove_step: Option<usize> = None;
                    for (i, st) in s.steps.iter_mut().enumerate() {
                        ui.group(|ui| {
                            ui.horizontal(|ui| {
                                ui.label(format!("{}.", i + 1));
                                egui::ComboBox::from_id_source(("cad_op", i))
                                    .selected_text(op_label(st.op))
                                    .width(92.0)
                                    .show_ui(ui, |ui| {
                                        for op in [Op::New, Op::Join, Op::Cut, Op::Intersect] {
                                            ui.selectable_value(&mut st.op, op, op_label(op));
                                        }
                                    });
                                egui::ComboBox::from_id_source(("cad_kind", i))
                                    .selected_text(match st.kind {
                                        FeatureKind::Box => "Box",
                                        FeatureKind::Cylinder => "Cylinder",
                                    })
                                    .width(92.0)
                                    .show_ui(ui, |ui| {
                                        ui.selectable_value(
                                            &mut st.kind,
                                            FeatureKind::Box,
                                            "Box",
                                        );
                                        ui.selectable_value(
                                            &mut st.kind,
                                            FeatureKind::Cylinder,
                                            "Cylinder",
                                        );
                                    });
                                if ui.small_button("✕").clicked() {
                                    remove_step = Some(i);
                                }
                            });
                            match st.kind {
                                FeatureKind::Box => {
                                    ui.horizontal(|ui| {
                                        ui.label("dx,dy,dz");
                                        dim_edit(ui, &mut st.dx);
                                        dim_edit(ui, &mut st.dy);
                                        dim_edit(ui, &mut st.dz);
                                    });
                                }
                                FeatureKind::Cylinder => {
                                    ui.horizontal(|ui| {
                                        ui.label("r, h");
                                        dim_edit(ui, &mut st.radius);
                                        dim_edit(ui, &mut st.height);
                                    });
                                }
                            }
                            ui.horizontal(|ui| {
                                ui.label("at x,y,z");
                                dim_edit(ui, &mut st.x);
                                dim_edit(ui, &mut st.y);
                                dim_edit(ui, &mut st.z);
                            });
                        });
                    }
                    if let Some(i) = remove_step {
                        s.steps.remove(i);
                    }
                    ui.horizontal(|ui| {
                        if ui.button("+ Box").clicked() {
                            s.steps.push(UiStep::new_box());
                        }
                        if ui.button("+ Cylinder").clicked() {
                            s.steps.push(UiStep::new_cylinder());
                        }
                    });
                    if ui.button("▶ Rebuild → viewport").clicked() {
                        match rebuild_tree(s) {
                            Ok((history, status)) => {
                                let k = history.len();
                                match tessellate_step(&history, k) {
                                    Ok(mesh) => {
                                        s.rebuilt_mesh = Some(mesh);
                                        s.push_rebuild = true;
                                        s.tree_status = Some(Ok(status));
                                        s.scrub = k;
                                        s.history = Some(history);
                                    }
                                    Err(e) => {
                                        s.history = None;
                                        s.rebuilt_mesh = None;
                                        s.push_rebuild = false;
                                        s.tree_status = Some(Err(e));
                                    }
                                }
                            }
                            Err(e) => {
                                s.history = None;
                                s.rebuilt_mesh = None;
                                s.push_rebuild = false;
                                s.tree_status = Some(Err(e));
                            }
                        }
                    }
                    if let Some(res) = &s.tree_status {
                        match res {
                            Ok(status) => ui.colored_label(
                                egui::Color32::from_rgb(80, 220, 120),
                                status,
                            ),
                            Err(e) => ui.colored_label(
                                egui::Color32::from_rgb(220, 120, 80),
                                format!("rebuild failed: {e}"),
                            ),
                        };
                    }

                    // History scrubber — roll the model back/forward through
                    // the per-step snapshots from the last rebuild, pushing the
                    // selected step into the viewport.
                    let n = s.history.as_ref().map_or(0, |h| h.len());
                    if n > 1 {
                        ui.add_space(2.0);
                        ui.label(egui::RichText::new("History").strong());
                        let resp = ui
                            .add(egui::Slider::new(&mut s.scrub, 1..=n).integer().text("step"));
                        if resp.changed() {
                            let scrub = s.scrub;
                            let mesh =
                                s.history.as_ref().and_then(|h| tessellate_step(h, scrub).ok());
                            if let Some(mesh) = mesh {
                                s.rebuilt_mesh = Some(mesh);
                                s.push_rebuild = true;
                            }
                        }
                        let label = s
                            .history
                            .as_ref()
                            .and_then(|h| h.get(s.scrub.saturating_sub(1)))
                            .map(|solid| {
                                format!("step {} / {n} — {} faces", s.scrub, solid.faces())
                            });
                        if let Some(label) = label {
                            ui.label(egui::RichText::new(label).small().monospace());
                        }
                    }
                });
        });

    // Deferred (outside the panel borrow): push the rebuilt solid's mesh into
    // the central 3-D viewport.
    if app.cad.push_rebuild {
        app.cad.push_rebuild = false;
        if let Some(mesh) = app.cad.rebuilt_mesh.take() {
            let quality = valenx_mesh::quality_report(&mesh);
            let aspect_hist =
                valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
            let skew_hist =
                valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
            app.stl = None;
            app.aero_field_overlay = None;
            app.mesh = Some(LoadedMesh {
                path: std::path::PathBuf::from("<cad>/feature-tree"),
                mesh,
                quality,
                aspect_hist,
                skew_hist,
            });
            app.frame_current_mesh();
        }
    }
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
            ..CadWorkbenchState::default()
        };
        let r = run_cad(&s);
        assert!(r.solved_radius.is_none());
        assert!(r.status.to_lowercase().contains("cyclic"), "status: {}", r.status);
    }

    #[test]
    fn feature_tree_rebuilds_punched_cube_with_history() {
        // The default tree is New box + Cut cylinder — the punched cube.
        let s = CadWorkbenchState::default();
        let (history, status) = rebuild_tree(&s).expect("default tree rebuilds");
        assert_eq!(history.len(), 2, "two steps → two snapshots");
        assert!(status.contains("faces"), "status: {status}");
        // The scrubber's intermediate history: step 1 is the bare box (6
        // faces), step 2 is the punched cube (more than 6).
        assert_eq!(history[0].faces(), 6, "first snapshot is the bare box");
        assert!(history[1].faces() > 6, "second snapshot is punched");
        // Every step tessellates to a non-empty viewport mesh.
        for k in 1..=history.len() {
            let mesh = tessellate_step(&history, k).expect("tessellate step");
            assert!(
                crate::mesh_loader::mesh_bounding_box(&mesh).is_some(),
                "step {k} should tessellate to a non-empty mesh"
            );
        }
    }

    #[test]
    fn feature_tree_reports_a_no_base_body_error() {
        // A lone Cut step has no body to cut from — surfaces as an error,
        // not a panic, and pushes nothing.
        let s = CadWorkbenchState {
            steps: vec![UiStep::new_cylinder()], // a single Cut step
            ..CadWorkbenchState::default()
        };
        let err = rebuild_tree(&s).expect_err("a lone Cut must fail");
        assert!(!err.is_empty());
    }
}
