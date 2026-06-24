//! The right-side **Projectile Workbench** panel — native point-mass
//! ballistics over `valenx-projectile`.
//!
//! Mirrors the Heat Transfer / Antenna workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_projectile_workbench`,
//! toggled from the View menu. The form sets a launch condition (speed,
//! elevation angle, gravity) and a flight model — drag-free vacuum
//! (closed form) or quadratic aerodynamic drag (RK4-integrated). "Analyze"
//! reports the range, apex height, time of flight and the range-maximising
//! launch angle; "Show 3-D" loads a representative ground plane with the
//! launch-to-landing trajectory arc as a 3-D solid into the central
//! viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_projectile::drag::{optimal_drag_angle, DragShot, OptimizeConfig};
use valenx_projectile::vacuum::{optimal_vacuum_angle_rad, VacuumShot};

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// The flight model the workbench evaluates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FlightModel {
    /// Drag-free closed-form kinematics ([`valenx_projectile::vacuum`]).
    Vacuum,
    /// Quadratic aerodynamic drag, RK4-integrated
    /// ([`valenx_projectile::drag`]).
    Drag,
}

/// Persistent form + result state for the Projectile Workbench.
pub struct ProjectileWorkbenchState {
    /// Launch speed `v0` (m/s).
    speed_m_per_s: f64,
    /// Launch elevation angle `θ` (degrees), in `[0, 90]`.
    angle_deg: f64,
    /// Gravitational acceleration `g` (m/s²).
    gravity_m_per_s2: f64,
    /// Flight model: vacuum or quadratic drag.
    model: FlightModel,
    /// Lumped quadratic-drag coefficient `k = ρ·C_d·A/(2m)` (1/m), used
    /// only by the drag model.
    drag_coeff_per_m: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D trajectory solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for ProjectileWorkbenchState {
    fn default() -> Self {
        // A 30 m/s launch at 45° under standard gravity: the vacuum range
        // peaks here at v0^2/g ~ 91.8 m, apex ~ 22.9 m, flight ~ 4.3 s.
        Self {
            speed_m_per_s: 30.0,
            angle_deg: 45.0,
            gravity_m_per_s2: 9.80665,
            model: FlightModel::Vacuum,
            drag_coeff_per_m: 0.01,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Projectile Workbench right-side panel. A no-op when the
/// `show_projectile_workbench` toggle is off.
pub fn draw_projectile_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_projectile_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_projectile_workbench",
        "Projectile",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native point-mass ballistics (vacuum + quadratic drag) · valenx-projectile",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.projectile;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Launch").strong());
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as
                    // its accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name otherwise, leaving it anonymous to a
                    // screen reader / AI driver).
                    ui.horizontal(|ui| {
                        let l = ui.label("speed v0 (m/s)");
                        ui.add(egui::DragValue::new(&mut s.speed_m_per_s).speed(0.5))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("angle θ (°)");
                        ui.add(egui::DragValue::new(&mut s.angle_deg).speed(0.5))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("gravity g (m/s²)");
                        ui.add(egui::DragValue::new(&mut s.gravity_m_per_s2).speed(0.05))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Flight model").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.model, FlightModel::Vacuum, "vacuum");
                        ui.radio_value(&mut s.model, FlightModel::Drag, "drag");
                    });
                    if s.model == FlightModel::Drag {
                        ui.horizontal(|ui| {
                            let l = ui.label("drag k = ρ·Cd·A/2m (1/m)");
                            ui.add(egui::DragValue::new(&mut s.drag_coeff_per_m).speed(0.001))
                                .labelled_by(l.id);
                        });
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_projectile(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative ground plane with the launch-to-landing trajectory arc as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Trajectory").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_projectile_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.projectile` borrow is
    // released here): build the trajectory's 3-D solid and load it.
    if app.projectile.show_3d_request {
        app.projectile.show_3d_request = false;
        load_trajectory_3d(app);
    }
}

/// Validate the form, evaluate the trajectory and format the readout.
fn run_projectile(s: &mut ProjectileWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The landing range and apex height `(range_m, apex_m)` for the current
/// model, the quantities the 3-D arc geometry needs. Extracted so it is
/// unit-testable and shared with the readout.
fn range_and_apex(s: &ProjectileWorkbenchState) -> Result<(f64, f64), String> {
    let angle_rad = s.angle_deg.to_radians();
    match s.model {
        FlightModel::Vacuum => {
            let shot = VacuumShot::new(s.speed_m_per_s, angle_rad, s.gravity_m_per_s2)
                .map_err(|e| e.to_string())?;
            Ok((shot.range(), shot.apex_height()))
        }
        FlightModel::Drag => {
            let shot = DragShot::new(
                s.speed_m_per_s,
                angle_rad,
                s.gravity_m_per_s2,
                s.drag_coeff_per_m,
            )
            .map_err(|e| e.to_string())?;
            let traj = shot.integrate(1e-3, 5_000_000).map_err(|e| e.to_string())?;
            Ok((traj.range, traj.apex_height))
        }
    }
}

/// Evaluate the trajectory and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &ProjectileWorkbenchState) -> Result<String, String> {
    let angle_rad = s.angle_deg.to_radians();
    match s.model {
        FlightModel::Vacuum => {
            let shot = VacuumShot::new(s.speed_m_per_s, angle_rad, s.gravity_m_per_s2)
                .map_err(|e| e.to_string())?;
            // Closed-form drag-free optimum is exactly 45°, independent of
            // speed and gravity.
            let opt_deg = optimal_vacuum_angle_rad().to_degrees();
            Ok(format!(
                "model           : vacuum (drag-free)\n\
                 speed v0        : {:.2} m/s\n\
                 angle θ         : {:.2}°\n\
                 gravity g       : {:.4} m/s²\n\n\
                 range R         : {:.3} m\n\
                 apex height H   : {:.3} m\n\
                 time of flight  : {:.3} s\n\
                 time to apex    : {:.3} s\n\
                 optimal angle   : {:.2}°",
                s.speed_m_per_s,
                s.angle_deg,
                s.gravity_m_per_s2,
                shot.range(),
                shot.apex_height(),
                shot.time_of_flight(),
                shot.time_to_apex(),
                opt_deg,
            ))
        }
        FlightModel::Drag => {
            let shot = DragShot::new(
                s.speed_m_per_s,
                angle_rad,
                s.gravity_m_per_s2,
                s.drag_coeff_per_m,
            )
            .map_err(|e| e.to_string())?;
            let traj = shot.integrate(1e-3, 5_000_000).map_err(|e| e.to_string())?;
            // Numerically search for the drag-corrected optimum over a
            // sensible interior bracket (1°..89°).
            let cfg = OptimizeConfig::new(
                s.speed_m_per_s,
                s.gravity_m_per_s2,
                s.drag_coeff_per_m,
                1e-3,
                5_000_000,
            )
            .map_err(|e| e.to_string())?;
            let opt =
                optimal_drag_angle(&cfg, 1.0_f64.to_radians(), 89.0_f64.to_radians(), 1e-5, 200)
                    .map_err(|e| e.to_string())?;
            Ok(format!(
                "model           : quadratic drag (RK4)\n\
                 speed v0        : {:.2} m/s\n\
                 angle θ         : {:.2}°\n\
                 gravity g       : {:.4} m/s²\n\
                 drag k          : {:.4} 1/m\n\n\
                 range R         : {:.3} m\n\
                 apex height H   : {:.3} m\n\
                 time of flight  : {:.3} s\n\
                 time to apex    : {:.3} s\n\
                 RK4 steps       : {}\n\
                 optimal angle   : {:.2}°",
                s.speed_m_per_s,
                s.angle_deg,
                s.gravity_m_per_s2,
                s.drag_coeff_per_m,
                traj.range,
                traj.apex_height,
                traj.time_of_flight,
                traj.apex_time,
                traj.steps,
                opt.angle_deg(),
            ))
        }
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

/// Build the trajectory as a triangle [`Mesh`] — a flat ground plane plus
/// the launch-to-landing arc rendered as a swept tube of small box
/// segments. The arc samples the *parabolic* vacuum shape scaled to the
/// model's true `(range, apex)`, so the solid matches the analyzed numbers
/// for both models. `None` for an invalid configuration.
fn trajectory_solid_mesh(s: &ProjectileWorkbenchState) -> Option<Mesh> {
    let (range, apex) = range_and_apex(s).ok()?;
    if range <= 0.0 || apex <= 0.0 {
        return None;
    }

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Ground plane: a thin slab spanning the flight, centred under the arc.
    let half_span = 0.5 * range;
    let depth = (0.12 * range).clamp(0.5, 8.0);
    let slab_thickness = (0.01 * range).clamp(0.05, 0.5);
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(half_span, 0.0, -slab_thickness),
        Vector3::new(half_span * 1.05, depth, slab_thickness),
    );

    // Arc: sample the height profile y(x) = 4·H·(x/R)·(1 − x/R), the
    // parabola through (0,0), apex H at x = R/2, and (R,0). Lay a small
    // box at each sample to sweep a tube along the flight path.
    let segments = 48_usize;
    let tube = (0.012 * range).clamp(0.08, 0.6);
    let mut prev = Vector3::new(0.0, 0.0, 0.0);
    for i in 0..=segments {
        let frac = i as f64 / segments as f64;
        let x = frac * range;
        let z = 4.0 * apex * frac * (1.0 - frac);
        let p = Vector3::new(x, 0.0, z);
        if i > 0 {
            let mid = 0.5 * (prev + p);
            push_box(&mut nodes, &mut tris, mid, Vector3::new(tube, tube, tube));
        }
        prev = p;
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-projectile");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D trajectory solid and load it into the central viewport.
fn load_trajectory_3d(app: &mut ValenxApp) {
    let Some(mesh) = trajectory_solid_mesh(&app.projectile) else {
        app.projectile.error =
            Some("launch parameters give no flight — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<trajectory>/valenx-projectile"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical projectile-motion workbench as a 3-D
/// solid plus its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn projectile_product() -> crate::WorkspaceProduct {
    let s = ProjectileWorkbenchState::default();
    let mesh = trajectory_solid_mesh(&s).expect("canonical projectile ⇒ trajectory solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<projectile>/valenx-trajectory");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical projectile ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Projectile motion (range/apex)".into(),
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
        let s = ProjectileWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_range_apex_and_optimum() {
        let mut s = ProjectileWorkbenchState::default();
        run_projectile(&mut s);
        assert!(
            s.error.is_none(),
            "default vacuum shot should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("range R"));
        assert!(s.result.contains("apex height H"));
        assert!(s.result.contains("time of flight"));
        // 30 m/s at 45° under g0: R = v0^2/g ~ 91.774 m (3 d.p.).
        assert!(s.result.contains("91.774"));
        // The drag-free optimum is exactly 45°.
        assert!(s.result.contains("45.00°"));
    }

    #[test]
    fn analyze_rejects_zero_speed() {
        let mut s = ProjectileWorkbenchState {
            speed_m_per_s: 0.0,
            ..ProjectileWorkbenchState::default()
        };
        run_projectile(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn vacuum_range_matches_hand_computed_formula() {
        // Ground truth: vacuum range R = v0^2·sin(2θ)/g, maximised at 45°
        // where sin(2θ) = 1 so R_max = v0^2/g. Hand-computed for
        // v0 = 30 m/s, g = 9.80665 m/s²: R_max = 900 / 9.80665 = 91.7755… m.
        let v0 = 30.0_f64;
        let g = 9.80665_f64;
        let hand = v0 * v0 / g;
        let shot = VacuumShot::new(v0, optimal_vacuum_angle_rad(), g).unwrap();
        assert!((shot.range() - hand).abs() < 1e-9);
        // No other angle beats the 45° range across the quadrant.
        for deg in 0..=90 {
            let theta = (deg as f64).to_radians();
            let r = VacuumShot::new(v0, theta, g).unwrap().range();
            assert!(r <= shot.range() + 1e-9, "angle {deg}° beat the 45° range");
        }
    }

    #[test]
    fn drag_model_analyzes_and_shortens_range() {
        let mut s = ProjectileWorkbenchState {
            model: FlightModel::Drag,
            ..ProjectileWorkbenchState::default()
        };
        run_projectile(&mut s);
        assert!(s.error.is_none(), "drag shot should analyze: {:?}", s.error);
        assert!(s.result.contains("quadratic drag"));
        // Drag must shorten the flight below the 91.766 m vacuum range.
        let (range, _apex) = range_and_apex(&s).unwrap();
        let vac = VacuumShot::new(
            s.speed_m_per_s,
            s.angle_deg.to_radians(),
            s.gravity_m_per_s2,
        )
        .unwrap()
        .range();
        assert!(range < vac, "drag range {range} not below vacuum {vac}");
    }

    #[test]
    fn trajectory_mesh_for_default_is_nonempty_and_in_range() {
        let s = ProjectileWorkbenchState::default();
        let mesh = trajectory_solid_mesh(&s).expect("default shot yields a solid");
        assert!(mesh.nodes.len() > 8, "expected ground plane + arc segments");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
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
            draw_projectile_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_projectile_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_projectile_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_projectile_workbench = true;
        run_projectile(&mut app.projectile);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Each numeric DragValue is a SpinButton; each must be `labelled_by`
        // its caption (egui clears a DragValue's own Name), so an AI / screen
        // reader can find the control by the caption text.
        let mut app = ValenxApp::default();
        app.show_projectile_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 3,
            "expected the projectile numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every projectile DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["speed v0 (m/s)", "angle θ (°)", "gravity g (m/s²)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
