//! The right-side **Drone / Multirotor Workbench** panel — native
//! multirotor hover performance over `valenx-drone`.
//!
//! Mirrors the Springs / Marine workbenches: a resizable [`egui::SidePanel`]
//! gated on `crate::ValenxApp::show_drone_workbench`, toggled from the View
//! menu. The form drives a [`valenx_drone::Multirotor`]; "Analyze" reports
//! the disk loading, the hover induced velocity and the ideal / actual
//! hover power plus a battery-limited hover endurance, and "Show 3-D drone"
//! loads a hub + rotor-disk solid into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_drone::{Multirotor, SEA_LEVEL_AIR_DENSITY};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Drone / Multirotor Workbench.
pub struct DroneWorkbenchState {
    /// Number of lifting rotors `n`.
    rotor_count: u32,
    /// Rotor radius `R` (m).
    rotor_radius_m: f64,
    /// All-up mass `m` (kg).
    mass_kg: f64,
    /// Rotor figure of merit `FM` in `(0, 1]`.
    figure_of_merit: f64,
    /// Air density `rho` (kg/m^3).
    air_density: f64,
    /// Usable battery energy (Wh) for the hover-endurance estimate.
    battery_wh: f64,
    /// Formatted hover-performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D drone solid and load it into the
    /// central viewport (serviced after the panel draws).
    show_3d_request: bool,
}

impl Default for DroneWorkbenchState {
    fn default() -> Self {
        // A 1.5 kg quadcopter with 15 cm rotors at sea level.
        Self {
            rotor_count: 4,
            rotor_radius_m: 0.15,
            mass_kg: 1.5,
            figure_of_merit: 0.70,
            air_density: SEA_LEVEL_AIR_DENSITY,
            battery_wh: 50.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Drone / Multirotor Workbench right-side panel. A no-op when the
/// `show_drone_workbench` toggle is off.
pub fn draw_drone_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_drone_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_drone_workbench",
        "Drone / Multirotor",
        |app, ui| {
            ui.label(
                egui::RichText::new("native multirotor hover performance · valenx-drone")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.drone;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Rotors + airframe").strong());
                    ui.horizontal(|ui| {
                        ui.label("rotor count n");
                        ui.add(egui::DragValue::new(&mut s.rotor_count).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("rotor radius R (m)");
                        ui.add(egui::DragValue::new(&mut s.rotor_radius_m).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("all-up mass (kg)");
                        ui.add(egui::DragValue::new(&mut s.mass_kg).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("figure of merit");
                        ui.add(egui::DragValue::new(&mut s.figure_of_merit).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Air + battery").strong());
                    ui.horizontal(|ui| {
                        ui.label("air ρ (kg/m³)");
                        ui.add(egui::DragValue::new(&mut s.air_density).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("battery (Wh)");
                        ui.add(egui::DragValue::new(&mut s.battery_wh).speed(1.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_drone(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D drone").strong())
                        .on_hover_text(
                            "Build the hub + rotor disks as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Hover performance").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_drone_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.drone` borrow is
    // released here): build the drone's 3-D solid and load it.
    if app.drone.show_3d_request {
        app.drone.show_3d_request = false;
        load_drone_3d(app);
    }
}

/// Build a validated [`Multirotor`] from the form, mapping the domain error
/// to a display string.
fn build_drone(s: &DroneWorkbenchState) -> Result<Multirotor, String> {
    Multirotor::new(
        s.rotor_count,
        s.rotor_radius_m,
        s.mass_kg,
        s.figure_of_merit,
        s.air_density,
    )
    .map_err(|e| e.to_string())
}

/// Validate the form, compute the hover performance and format the readout.
fn run_drone(s: &mut DroneWorkbenchState) {
    s.error = None;
    match build_drone(s) {
        Ok(m) => {
            let p = m.hover_performance();
            let endurance = m.hover_endurance_minutes(s.battery_wh.max(0.0));
            s.result = format!(
                "rotors n       : {}\n\
                 rotor radius R : {:.3} m\n\
                 all-up mass m  : {:.2} kg\n\
                 figure of merit: {:.2}\n\n\
                 weight         : {:.1} N\n\
                 disk area      : {:.4} m²\n\
                 disk loading   : {:.1} N/m²\n\
                 induced vel vi : {:.2} m/s\n\
                 ideal power    : {:.0} W\n\
                 actual power   : {:.0} W  (P_ideal / FM)\n\
                 hover endurance: {:.1} min  ({:.0} Wh battery)",
                s.rotor_count,
                s.rotor_radius_m,
                s.mass_kg,
                s.figure_of_merit,
                p.weight_n,
                p.disk_area_m2,
                p.disk_loading_pa,
                p.induced_velocity_m_s,
                p.ideal_power_w,
                p.actual_power_w,
                endurance,
                s.battery_wh,
            );
        }
        Err(e) => s.error = Some(e),
    }
}

/// Append an outward-facing box (centre `c`, half-extents) to the buffers.
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

/// Append a flat (double-sided) disk of radius `r` in the horizontal plane
/// at centre `c`.
fn push_disk(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    r: f64,
    seg: usize,
) {
    let center = nodes.len();
    nodes.push(c);
    let ring = nodes.len();
    for j in 0..seg {
        let t = j as f64 / seg as f64 * TAU;
        nodes.push(c + Vector3::new(r * t.cos(), r * t.sin(), 0.0));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
        // Both windings so the disk shades from above and below.
        tris.extend_from_slice(&[center, ring + j, ring + jn, center, ring + jn, ring + j]);
    }
}

/// Append a flat (double-sided) arm strip of half-width `w` from `from` to
/// `to`, in the horizontal plane.
fn push_arm(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    from: Vector3<f64>,
    to: Vector3<f64>,
    w: f64,
) {
    let dir = to - from;
    if dir.norm() < 1e-9 {
        return;
    }
    let perp = Vector3::new(-dir.y, dir.x, 0.0).normalize() * w;
    let base = nodes.len();
    nodes.push(from + perp);
    nodes.push(from - perp);
    nodes.push(to - perp);
    nodes.push(to + perp);
    tris.extend_from_slice(&[
        base,
        base + 1,
        base + 2,
        base,
        base + 2,
        base + 3,
        base,
        base + 2,
        base + 1,
        base,
        base + 3,
        base + 2,
    ]);
}

/// Build the multirotor as a triangle [`Mesh`] — a central hub box with `n`
/// rotor disks on flat arms, evenly spaced around it. `None` for an invalid
/// configuration.
fn drone_solid_mesh(s: &DroneWorkbenchState) -> Option<Mesh> {
    let m = build_drone(s).ok()?;
    let n = m.rotor_count;
    let r = m.rotor_radius_m;
    let arm = r * 2.3;
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    // Central hub.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::zeros(),
        Vector3::new(r * 0.5, r * 0.5, r * 0.25),
    );
    for i in 0..n {
        let a = i as f64 / n as f64 * TAU;
        let c = Vector3::new(arm * a.cos(), arm * a.sin(), 0.0);
        push_arm(&mut nodes, &mut tris, Vector3::zeros(), c, r * 0.12);
        push_disk(&mut nodes, &mut tris, c, r, 16);
    }
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-drone");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D drone solid and load it into the central viewport.
fn load_drone_3d(app: &mut ValenxApp) {
    let Some(mesh) = drone_solid_mesh(&app.drone) else {
        app.drone.error = Some("drone parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<drone>/valenx-drone"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"drone"}`** product: the canonical hub +
/// rotor-disk multirotor solid (the panel's "Show 3-D drone" geometry) paired
/// with the workbench's own hover-performance headline numbers, at a fixed
/// 3/4 camera. Registered in [`crate::products_registry`]; the per-tool
/// builder the registry dispatches to. Pure — driven off
/// [`DroneWorkbenchState::default`].
///
/// The readout rows mirror the panel's `run_drone` "Hover performance"
/// readout (this workbench formats its readout into `result` rather than via a
/// shared `compute()` string fn, so the rows are taken from that here).
pub(crate) fn drone_product() -> crate::WorkspaceProduct {
    let mut s = DroneWorkbenchState::default();
    run_drone(&mut s);
    let mesh = drone_solid_mesh(&s).expect("default quadcopter ⇒ a 3-D solid");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<drone>/valenx-drone");
    let lines = crate::products_registry::lines_from_readout(&s.result);
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Drone / Multirotor".into(),
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
        let s = DroneWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_drone_reports_power_and_endurance() {
        let mut s = DroneWorkbenchState::default();
        run_drone(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("actual power"));
        assert!(s.result.contains("hover endurance"));
        assert!(s.result.contains("disk loading"));
    }

    #[test]
    fn analyze_rejects_zero_rotors_and_bad_fm() {
        let mut z = DroneWorkbenchState {
            rotor_count: 0,
            ..Default::default()
        };
        run_drone(&mut z);
        assert!(z.error.is_some());
        let mut fm = DroneWorkbenchState {
            figure_of_merit: 1.5,
            ..Default::default()
        };
        run_drone(&mut fm);
        assert!(fm.error.is_some());
    }

    #[test]
    fn drone_mesh_for_default_is_nonempty_and_in_range() {
        let s = DroneWorkbenchState::default();
        let mesh = drone_solid_mesh(&s).expect("default drone yields a solid");
        assert!(mesh.nodes.len() > 8, "expected hub + rotor disks");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn drone_mesh_none_for_invalid() {
        let s = DroneWorkbenchState {
            rotor_radius_m: 0.0,
            ..Default::default()
        };
        assert!(drone_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_drone_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_drone_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_drone_workbench = true;
        run_drone(&mut app.drone);
        draw_workbench(&mut app);
    }
}
