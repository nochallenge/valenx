//! The right-side **Rail / Train Workbench** panel — native train running
//! resistance and tractive performance over `valenx-rail`.
//!
//! Mirrors the Drone / Springs / Marine workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_rail_workbench`,
//! toggled from the View menu. The form drives a [`valenx_rail::Train`];
//! "Analyze" reports the Davis running resistance, the grade resistance,
//! the net tractive force, the acceleration, the drawbar power and the
//! constant-tractive-effort balancing speed at the chosen speed and grade,
//! and "Show 3-D train" loads a boxcar solid (body + underframe + wheels)
//! into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_rail::Train;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Rail / Train Workbench.
pub struct RailWorkbenchState {
    /// Total mass `m` (kg).
    mass_kg: f64,
    /// Davis constant term `A` (N).
    davis_a_n: f64,
    /// Davis linear coefficient `B` (N*s/m).
    davis_b_n_s_per_m: f64,
    /// Davis quadratic (aerodynamic) coefficient `C` (N*s^2/m^2).
    davis_c_n_s2_per_m2: f64,
    /// Constant tractive effort `TE` (N).
    max_tractive_effort_n: f64,
    /// Operating speed `v` (m/s) the readout is evaluated at.
    speed_m_s: f64,
    /// Track grade (rise / run) as a percentage the readout is evaluated at.
    grade_percent: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D boxcar solid and load it into the
    /// central viewport (serviced after the panel draws).
    show_3d_request: bool,
}

impl Default for RailWorkbenchState {
    fn default() -> Self {
        // A heavy ~4000 t freight train on a 1.4 % grade at 20 m/s.
        Self {
            mass_kg: 4000000.0,
            davis_a_n: 30000.0,
            davis_b_n_s_per_m: 350.0,
            davis_c_n_s2_per_m2: 12.0,
            max_tractive_effort_n: 600000.0,
            speed_m_s: 20.0,
            grade_percent: 1.4,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Rail / Train Workbench right-side panel. A no-op when the
/// `show_rail_workbench` toggle is off.
pub fn draw_rail_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_rail_workbench {
        return;
    }

    egui::SidePanel::right("valenx_rail_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Rail / Train",
                "native train resistance + tractive effort · valenx-rail",
            ) {
                app.show_rail_workbench = false;
            }

            let s = &mut app.rail;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Train").strong());
                    ui.horizontal(|ui| {
                        ui.label("mass m (kg)");
                        ui.add(egui::DragValue::new(&mut s.mass_kg).speed(1000.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("tractive effort TE (N)");
                        ui.add(egui::DragValue::new(&mut s.max_tractive_effort_n).speed(1000.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Davis resistance R = A + B·v + C·v²").strong());
                    ui.horizontal(|ui| {
                        ui.label("A (N)");
                        ui.add(egui::DragValue::new(&mut s.davis_a_n).speed(100.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("B (N·s/m)");
                        ui.add(egui::DragValue::new(&mut s.davis_b_n_s_per_m).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("C (N·s²/m²)");
                        ui.add(egui::DragValue::new(&mut s.davis_c_n_s2_per_m2).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Operating point").strong());
                    ui.horizontal(|ui| {
                        ui.label("speed v (m/s)");
                        ui.add(egui::DragValue::new(&mut s.speed_m_s).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("grade (%)");
                        ui.add(egui::DragValue::new(&mut s.grade_percent).speed(0.05));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_train(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D train").strong())
                        .on_hover_text(
                            "Build a boxcar (body + underframe + wheels) as a 3-D solid and load it into the central viewport to orbit",
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
        });

    // Serviced after the panel draws (the `&mut app.rail` borrow is
    // released here): build the boxcar's 3-D solid and load it.
    if app.rail.show_3d_request {
        app.rail.show_3d_request = false;
        load_train_3d(app);
    }
}

/// Build a validated [`Train`] from the form, mapping the domain error to a
/// display string.
fn build_train(s: &RailWorkbenchState) -> Result<Train, String> {
    Train::new(
        s.mass_kg,
        s.davis_a_n,
        s.davis_b_n_s_per_m,
        s.davis_c_n_s2_per_m2,
        s.max_tractive_effort_n,
    )
    .map_err(|e| e.to_string())
}

/// Validate the form, compute the running / tractive performance at the
/// chosen speed and grade, and format the readout.
fn run_train(s: &mut RailWorkbenchState) {
    s.error = None;
    match build_train(s) {
        Ok(t) => {
            let v = s.speed_m_s.max(0.0);
            let g = s.grade_percent / 100.0;
            let bal = match t.balancing_speed(g) {
                Some(b) => format!("{:.1} m/s ({:.0} km/h)", b, b * 3.6),
                None => "— (cannot sustain a steady speed)".to_string(),
            };
            s.result = format!(
                "mass m         : {:.0} kg\n\
                 tractive effort: {:.0} N\n\
                 Davis A B C    : {:.0} / {:.0} / {:.1}\n\n\
                 at v = {:.1} m/s, grade = {:.2} %:\n\
                 running resist : {:.0} N\n\
                 grade resist   : {:.0} N\n\
                 total resist   : {:.0} N\n\
                 net tractive   : {:.0} N\n\
                 acceleration   : {:.4} m/s²\n\
                 drawbar power  : {:.0} kW\n\n\
                 max start grade: {:.2} %\n\
                 balancing speed: {}",
                s.mass_kg,
                s.max_tractive_effort_n,
                s.davis_a_n,
                s.davis_b_n_s_per_m,
                s.davis_c_n_s2_per_m2,
                v,
                s.grade_percent,
                t.running_resistance(v),
                t.grade_resistance(g),
                t.total_resistance(v, g),
                t.net_tractive_force(v, g),
                t.acceleration(v, g),
                t.drawbar_power(v, g) / 1000.0,
                t.max_startable_grade() * 100.0,
                bal,
            );
        }
        Err(e) => s.error = Some(e),
    }
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

/// Append a (double-sided) wheel — a short cylinder whose axis is the
/// cross-track `y` direction — of radius `r` and half-width `hw` at centre
/// `c`, with `seg` segments around the rim.
fn push_wheel(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    r: f64,
    hw: f64,
    seg: usize,
) {
    let left_center = nodes.len();
    nodes.push(c + Vector3::new(0.0, -hw, 0.0));
    let right_center = nodes.len();
    nodes.push(c + Vector3::new(0.0, hw, 0.0));
    let left_ring = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(c + Vector3::new(r * a.cos(), -hw, r * a.sin()));
    }
    let right_ring = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(c + Vector3::new(r * a.cos(), hw, r * a.sin()));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
        // Both caps, double-sided.
        tris.extend_from_slice(&[
            left_center,
            left_ring + j,
            left_ring + jn,
            left_center,
            left_ring + jn,
            left_ring + j,
        ]);
        tris.extend_from_slice(&[
            right_center,
            right_ring + j,
            right_ring + jn,
            right_center,
            right_ring + jn,
            right_ring + j,
        ]);
        // Tread (side wall) quad, double-sided.
        tris.extend_from_slice(&[
            left_ring + j,
            right_ring + j,
            right_ring + jn,
            left_ring + j,
            right_ring + jn,
            left_ring + jn,
            left_ring + j,
            right_ring + jn,
            right_ring + j,
            left_ring + j,
            left_ring + jn,
            right_ring + jn,
        ]);
    }
}

/// Build the train as a triangle [`Mesh`] — a boxcar body, an underframe and
/// eight wheels on four axles. `None` for an invalid (out-of-domain) train.
/// The geometry is representative (a generic freight car), not derived from
/// the resistance parameters.
fn train_solid_mesh(s: &RailWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D build on a physically valid train, like the other
    // workbenches.
    let _train = build_train(s).ok()?;
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    // Car body.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 2.6),
        Vector3::new(10.0, 1.5, 1.6),
    );
    // Underframe / chassis (slightly longer, narrower, thin).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 1.1),
        Vector3::new(10.5, 1.2, 0.3),
    );
    // Eight wheels: four axles, two wheels each, axes cross-track.
    let rw = 0.8;
    let hw = 0.25;
    for &x in &[-8.5, -6.0, 6.0, 8.5] {
        for &y in &[-1.25, 1.25] {
            push_wheel(&mut nodes, &mut tris, Vector3::new(x, y, rw), rw, hw, 16);
        }
    }
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-rail");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D boxcar solid and load it into the central viewport.
fn load_train_3d(app: &mut ValenxApp) {
    let Some(mesh) = train_solid_mesh(&app.rail) else {
        app.rail.error = Some("train parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<train>/valenx-rail"),
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
        let s = RailWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_train_reports_resistance_and_accel() {
        let mut s = RailWorkbenchState::default();
        run_train(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("running resist"));
        assert!(s.result.contains("acceleration"));
        assert!(s.result.contains("balancing speed"));
    }

    #[test]
    fn analyze_rejects_out_of_domain_inputs() {
        let mut zero_mass = RailWorkbenchState {
            mass_kg: 0.0,
            ..Default::default()
        };
        run_train(&mut zero_mass);
        assert!(zero_mass.error.is_some());
        let mut zero_c = RailWorkbenchState {
            davis_c_n_s2_per_m2: 0.0,
            ..Default::default()
        };
        run_train(&mut zero_c);
        assert!(zero_c.error.is_some());
    }

    #[test]
    fn train_mesh_for_default_is_nonempty_and_in_range() {
        let s = RailWorkbenchState::default();
        let mesh = train_solid_mesh(&s).expect("default train yields a solid");
        assert!(mesh.nodes.len() > 8, "expected body + underframe + wheels");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn train_mesh_none_for_invalid() {
        let s = RailWorkbenchState {
            mass_kg: 0.0,
            ..Default::default()
        };
        assert!(train_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_rail_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_rail_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_rail_workbench = true;
        run_train(&mut app.rail);
        draw_workbench(&mut app);
    }
}
