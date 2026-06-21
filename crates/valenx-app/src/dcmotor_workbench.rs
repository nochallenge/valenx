//! The right-side **DC Motor Workbench** panel — native brushed-DC-motor
//! performance over `valenx-dcmotor`.
//!
//! Mirrors the Truss / Wind Turbine workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_dcmotor_workbench`,
//! toggled from the View menu. The form drives a [`valenx_dcmotor::DcMotor`];
//! "Analyze" reports the stall torque / current, the no-load speed, the peak
//! output power and the operating point (current, torque, power, efficiency)
//! at a chosen shaft speed, and "Show 3-D motor" loads a can-plus-shaft solid
//! into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_dcmotor::DcMotor;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// rad/s -> rev/min.
fn rpm(omega: f64) -> f64 {
    omega * 60.0 / TAU
}

/// Persistent form + result state for the DC Motor Workbench.
pub struct DcMotorWorkbenchState {
    /// Armature resistance `R` (ohm).
    resistance_ohm: f64,
    /// Back-EMF constant `Ke` (V*s/rad).
    ke: f64,
    /// Torque constant `Kt` (N*m/A).
    kt: f64,
    /// Supply (terminal) voltage `V` (volts).
    supply_v: f64,
    /// Operating shaft speed `omega` (rad/s) for the operating-point readout.
    omega_rad_s: f64,
    /// Formatted readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D motor solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for DcMotorWorkbenchState {
    fn default() -> Self {
        // A small coherent-SI brushed motor (Kt == Ke), 12 V supply.
        Self {
            resistance_ohm: 1.0,
            ke: 0.05,
            kt: 0.05,
            supply_v: 12.0,
            omega_rad_s: 120.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the DC Motor Workbench right-side panel. A no-op when the
/// `show_dcmotor_workbench` toggle is off.
pub fn draw_dcmotor_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_dcmotor_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_dcmotor_workbench",
        "DC Motor",
        |app, ui| {
            ui.label(
                egui::RichText::new("native brushed-DC-motor performance · valenx-dcmotor")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.dcmotor;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Motor constants").strong());
                    ui.horizontal(|ui| {
                        ui.label("resistance R (Ω)");
                        ui.add(egui::DragValue::new(&mut s.resistance_ohm).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("back-EMF Ke (V·s/rad)");
                        ui.add(egui::DragValue::new(&mut s.ke).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("torque Kt (N·m/A)");
                        ui.add(egui::DragValue::new(&mut s.kt).speed(0.005));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Operating point").strong());
                    ui.horizontal(|ui| {
                        ui.label("supply V (V)");
                        ui.add(egui::DragValue::new(&mut s.supply_v).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("shaft ω (rad/s)");
                        ui.add(egui::DragValue::new(&mut s.omega_rad_s).speed(2.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_motor(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D motor").strong())
                        .on_hover_text(
                            "Build the motor can + output shaft as a 3-D solid and load it into the central viewport to orbit",
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
        app.show_dcmotor_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.dcmotor` borrow is
    // released here): build the motor's 3-D solid and load it.
    if app.dcmotor.show_3d_request {
        app.dcmotor.show_3d_request = false;
        load_motor_3d(app);
    }
}

/// Build a validated [`DcMotor`] from the form, mapping the domain error to a
/// display string.
fn build_motor(s: &DcMotorWorkbenchState) -> Result<DcMotor, String> {
    DcMotor::new(s.resistance_ohm, s.ke, s.kt).map_err(|e| e.to_string())
}

/// Validate the form, compute the motor performance and format the readout.
fn run_motor(s: &mut DcMotorWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Compute the full readout, mapping any domain error to a display string.
/// Extracted so it is unit-testable.
fn compute(s: &DcMotorWorkbenchState) -> Result<String, String> {
    let m = build_motor(s)?;
    let v = s.supply_v;
    let stall_t = m.stall_torque(v).map_err(|e| e.to_string())?;
    let stall_i = m.stall_current(v).map_err(|e| e.to_string())?;
    let no_load = m.no_load_speed(v).map_err(|e| e.to_string())?;
    let mpp = m.max_power_point(v).map_err(|e| e.to_string())?;
    let op = m
        .operating_point(v, s.omega_rad_s)
        .map_err(|e| e.to_string())?;
    Ok(format!(
        "resistance R : {:.3} Ω\n\
         Ke / Kt      : {:.4} / {:.4}\n\
         supply V     : {:.1} V\n\n\
         stall torque : {:.3} N·m\n\
         stall current: {:.1} A\n\
         no-load speed: {:.0} rad/s ({:.0} rpm)\n\n\
         max output P : {:.1} W  at {:.0} rpm, {:.3} N·m\n\n\
         at ω = {:.0} rad/s ({:.0} rpm):\n\
         current   : {:.2} A\n\
         torque    : {:.3} N·m\n\
         mech power: {:.1} W\n\
         efficiency: {:.1} %",
        s.resistance_ohm,
        s.ke,
        s.kt,
        v,
        stall_t,
        stall_i,
        no_load,
        rpm(no_load),
        mpp.mechanical_power_w,
        rpm(mpp.omega_rad_s),
        mpp.torque_nm,
        s.omega_rad_s,
        rpm(s.omega_rad_s),
        op.current_a,
        op.torque_nm,
        op.mechanical_power_w,
        op.efficiency * 100.0,
    ))
}

/// Append a capped cylinder along the x-axis, double-sided, centred at
/// `center` with half-length `half_len` and `radius`, `seg` segments.
fn push_cyl_x(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    center: Vector3<f64>,
    half_len: f64,
    radius: f64,
    seg: usize,
) {
    let (x0, x1) = (center.x - half_len, center.x + half_len);
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(
            x0,
            center.y + radius * a.cos(),
            center.z + radius * a.sin(),
        ));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(
            x1,
            center.y + radius * a.cos(),
            center.z + radius * a.sin(),
        ));
    }
    let cap0 = nodes.len();
    nodes.push(Vector3::new(x0, center.y, center.z));
    let cap1 = nodes.len();
    nodes.push(Vector3::new(x1, center.y, center.z));
    for j in 0..seg {
        let jn = (j + 1) % seg;
        // Side wall (double-sided).
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
        // Caps (double-sided fans).
        tris.extend_from_slice(&[cap0, lo + jn, lo + j, cap0, lo + j, lo + jn]);
        tris.extend_from_slice(&[cap1, hi + j, hi + jn, cap1, hi + jn, hi + j]);
    }
}

/// Build the motor as a triangle [`Mesh`] — a cylindrical can with an output
/// shaft protruding from the front face. Representative geometry (the
/// electrical performance is the single-machine model). `None` for an
/// invalid motor.
fn motor_solid_mesh(s: &DcMotorWorkbenchState) -> Option<Mesh> {
    build_motor(s).ok()?;
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    // Can / body.
    push_cyl_x(&mut nodes, &mut tris, Vector3::zeros(), 1.2, 1.0, 24);
    // Output shaft, protruding from the +x face.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(1.8, 0.0, 0.0),
        0.6,
        0.18,
        16,
    );
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-dcmotor");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D motor solid and load it into the central viewport.
fn load_motor_3d(app: &mut ValenxApp) {
    let Some(mesh) = motor_solid_mesh(&app.dcmotor) else {
        app.dcmotor.error =
            Some("motor parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<motor>/valenx-dcmotor"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical DC-motor workbench as a 3-D solid plus
/// its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn dcmotor_product() -> crate::WorkspaceProduct {
    let s = DcMotorWorkbenchState::default();
    let mesh = motor_solid_mesh(&s).expect("canonical DC motor ⇒ motor solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<dcmotor>/valenx-motor");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical DC motor ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "DC motor (torque/speed/efficiency)".into(),
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
        let s = DcMotorWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_stall_and_power() {
        let mut s = DcMotorWorkbenchState::default();
        run_motor(&mut s);
        assert!(s.error.is_none(), "default motor analyzes: {:?}", s.error);
        assert!(s.result.contains("stall torque"));
        assert!(s.result.contains("no-load speed"));
        assert!(s.result.contains("efficiency"));
    }

    #[test]
    fn analyze_rejects_nonpositive_resistance() {
        let mut s = DcMotorWorkbenchState {
            resistance_ohm: 0.0,
            ..Default::default()
        };
        run_motor(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn motor_mesh_for_default_is_nonempty_and_in_range() {
        let s = DcMotorWorkbenchState::default();
        let mesh = motor_solid_mesh(&s).expect("default motor yields a solid");
        assert!(mesh.nodes.len() > 8, "expected can + shaft");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn motor_mesh_none_for_invalid() {
        let s = DcMotorWorkbenchState {
            kt: 0.0,
            ..Default::default()
        };
        assert!(motor_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_dcmotor_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_dcmotor_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_dcmotor_workbench = true;
        run_motor(&mut app.dcmotor);
        draw_workbench(&mut app);
    }
}
