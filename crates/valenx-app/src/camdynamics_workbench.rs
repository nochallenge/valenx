//! The right-side **Cam Dynamics Workbench** panel — native closed-form
//! cam-follower rise kinematics over `valenx-camdynamics`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_camdynamics_workbench`,
//! toggled from the View menu. The form sets a follower rise (lift, rise
//! angle `beta`, motion law); "Analyze" builds a [`RiseProfile`] and reports
//! the endpoint displacements, the peak velocity and acceleration (and their
//! dimensionless coefficients), and "Show 3-D" loads a representative disc
//! cam (an extruded disc with an eccentric lift lobe) into the central
//! viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_camdynamics::{MotionLaw, RiseProfile};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Cam Dynamics Workbench.
pub struct CamDynamicsWorkbenchState {
    /// Follower lift `L` (total rise displacement, mm).
    lift_mm: f64,
    /// Rise angle `beta` over which the lift occurs (degrees).
    beta_deg: f64,
    /// The motion law connecting the rise.
    law: MotionLaw,
    /// Formatted kinematic readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D disc cam (serviced after the panel
    /// draws).
    show_3d_request: bool,
}

impl Default for CamDynamicsWorkbenchState {
    fn default() -> Self {
        // A 10 mm follower rise over a 90 deg cam interval using the smooth
        // cycloidal law — a representative automotive / machine-tool cam.
        Self {
            lift_mm: 10.0,
            beta_deg: 90.0,
            law: MotionLaw::Cycloidal,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Cam Dynamics Workbench right-side panel. A no-op when the
/// `show_camdynamics_workbench` toggle is off.
pub fn draw_camdynamics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_camdynamics_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_camdynamics_workbench",
        "Cam Dynamics",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native closed-form cam-follower rise kinematics · valenx-camdynamics",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.camdynamics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Rise").strong());
                    ui.horizontal(|ui| {
                        ui.label("lift L (mm)");
                        ui.add(egui::DragValue::new(&mut s.lift_mm).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("rise angle β (deg)");
                        ui.add(egui::DragValue::new(&mut s.beta_deg).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Motion law").strong());
                    ui.radio_value(
                        &mut s.law,
                        MotionLaw::SimpleHarmonic,
                        "simple-harmonic (SHM)",
                    );
                    ui.radio_value(&mut s.law, MotionLaw::Cycloidal, "cycloidal");

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_camdynamics(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative disc cam (an extruded disc with an eccentric lift lobe scaled by the configured lift) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Kinematics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_camdynamics_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.camdynamics` borrow is
    // released here): build the disc cam's 3-D solid and load it.
    if app.camdynamics.show_3d_request {
        app.camdynamics.show_3d_request = false;
        load_cam_3d(app);
    }
}

/// Validate the form, evaluate the rise profile and format the readout.
fn run_camdynamics(s: &mut CamDynamicsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`RiseProfile`] from the form (converting the rise
/// angle from degrees to radians), the quantity both the readout and the
/// 3-D gate need. Extracted so it is unit-testable and shared.
fn profile(s: &CamDynamicsWorkbenchState) -> Result<RiseProfile, String> {
    let beta_rad = s.beta_deg.to_radians();
    RiseProfile::new(s.lift_mm, beta_rad, s.law).map_err(|e| e.to_string())
}

/// Evaluate the rise profile and format the full kinematic readout, mapping
/// any domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &CamDynamicsWorkbenchState) -> Result<String, String> {
    let p = profile(s)?;
    let law = s.law.name();
    let beta_deg = s.beta_deg;
    let beta = p.beta();
    let lift = p.lift();
    let s_start = p.at(0.0).displacement;
    let s_end = p.at(beta).displacement;
    let s_mid = p.at(beta / 2.0).displacement;
    let v_peak = p.peak_velocity();
    let a_peak = p.peak_acceleration();
    // Dimensionless kinematic coefficients Cv = v_max·beta/L and
    // Ca = a_max·beta^2/L (only defined for a non-zero lift).
    let (cv, ca) = if lift > 0.0 {
        (v_peak * beta / lift, a_peak * beta * beta / lift)
    } else {
        (0.0, 0.0)
    };

    Ok(format!(
        "law             : {law}\n\
         lift L          : {lift:.3} mm\n\
         rise angle β    : {beta_deg:.2} deg ({beta:.4} rad)\n\n\
         displacement s(θ) [mm]\n\
         at start (θ=0)  : {s_start:.4}\n\
         at mid-rise     : {s_mid:.4}\n\
         at end  (θ=β)   : {s_end:.4}\n\n\
         peak velocity   : {v_peak:.4} mm/rad  (Cv {cv:.4})\n\
         peak accel      : {a_peak:.4} mm/rad²  (Ca {ca:.4})"
    ))
}

/// Append a closed disc-cam shell (an extruded disc whose radius is modulated
/// by a single eccentric lift lobe) to the buffers: a bottom ring and a top
/// ring connected by a side wall and capped with triangle fans. `seg` is the
/// angular segment count, `base_r` the base radius, `lobe` the lobe height,
/// `thick` the axial thickness.
fn push_disc_cam(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    seg: usize,
    base_r: f64,
    lobe: f64,
    thick: f64,
) {
    // Lobe profile: one smooth bump in [0, TAU) reaching `lobe` at angle 0,
    // using a raised-cosine so the rim closes continuously.
    let radius_at = |phi: f64| base_r + lobe * 0.5 * (1.0 + phi.cos());

    let base = nodes.len();
    // Bottom ring (z = 0), then top ring (z = thick).
    for ring in 0..2 {
        let z = ring as f64 * thick;
        for j in 0..seg {
            let phi = j as f64 / seg as f64 * TAU;
            let r = radius_at(phi);
            nodes.push(Vector3::new(r * phi.cos(), r * phi.sin(), z));
        }
    }

    // Side wall: connect the bottom ring to the top ring.
    for j in 0..seg {
        let jn = (j + 1) % seg;
        let b0 = base + j;
        let b1 = base + jn;
        let t0 = base + seg + j;
        let t1 = base + seg + jn;
        tris.extend_from_slice(&[b0, t0, t1, b0, t1, b1]);
    }

    // End caps: a triangle fan from each ring's centre point.
    let bottom_center = nodes.len();
    nodes.push(Vector3::new(0.0, 0.0, 0.0));
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[bottom_center, base + jn, base + j]);
    }
    let top_center = nodes.len();
    nodes.push(Vector3::new(0.0, 0.0, thick));
    for j in 0..seg {
        let jn = (j + 1) % seg;
        tris.extend_from_slice(&[top_center, base + seg + j, base + seg + jn]);
    }
}

/// Build the disc cam as a triangle [`Mesh`] — a representative extruded disc
/// whose rim carries a single eccentric lift lobe scaled by the configured
/// lift. Representative geometry (not a manufactured profile; the reported
/// numbers are the `valenx-camdynamics` rise kinematics). `None` for an
/// invalid configuration.
fn cam_solid_mesh(s: &CamDynamicsWorkbenchState) -> Option<Mesh> {
    let p = profile(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // A base disc with a single lobe whose height is the follower lift; the
    // base radius keeps the shape disc-like at a representative scale.
    let lift = p.lift();
    let base_r = (lift.max(1.0)) * 2.0;
    let thick = (lift.max(1.0)) * 0.6;
    push_disc_cam(&mut nodes, &mut tris, 64, base_r, lift, thick);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-camdynamics");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D disc cam solid and load it into the central viewport.
fn load_cam_3d(app: &mut ValenxApp) {
    let Some(mesh) = cam_solid_mesh(&app.camdynamics) else {
        app.camdynamics.error =
            Some("cam parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<cam>/valenx-camdynamics"),
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
        let s = CamDynamicsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_kinematics() {
        let mut s = CamDynamicsWorkbenchState::default();
        run_camdynamics(&mut s);
        assert!(
            s.error.is_none(),
            "default cam should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("peak velocity"));
        assert!(s.result.contains("peak accel"));
        assert!(s.result.contains("cycloidal"));
    }

    #[test]
    fn analyze_rejects_zero_rise_angle() {
        let mut s = CamDynamicsWorkbenchState {
            beta_deg: 0.0,
            ..Default::default()
        };
        run_camdynamics(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn ground_truth_endpoints_and_cycloidal_peak_velocity() {
        // Ground truth: the rise displacement is 0 at the start (θ=0) and
        // exactly the lift at the end (θ=β); and for the default cycloidal
        // law the peak velocity is the closed-form v_max = 2 L / beta.
        let s = CamDynamicsWorkbenchState::default();
        let p = profile(&s).expect("default profile builds");
        let start = p.at(0.0).displacement;
        let end = p.at(p.beta()).displacement;
        assert!(start.abs() < 1e-9, "s(0) should be 0, got {start}");
        assert!(
            (end - p.lift()).abs() < 1e-9,
            "s(beta) should equal lift {}, got {end}",
            p.lift()
        );
        let expected_v = 2.0 * p.lift() / p.beta();
        let got_v = p.peak_velocity();
        assert!(
            (got_v - expected_v).abs() < 1e-9,
            "cycloidal peak velocity should be {expected_v}, got {got_v}"
        );
    }

    #[test]
    fn cam_mesh_for_default_is_nonempty_and_in_range() {
        let s = CamDynamicsWorkbenchState::default();
        let mesh = cam_solid_mesh(&s).expect("default cam yields a solid");
        assert!(mesh.nodes.len() > 8, "expected two rings plus centres");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn cam_mesh_none_for_invalid() {
        let s = CamDynamicsWorkbenchState {
            beta_deg: 0.0,
            ..Default::default()
        };
        assert!(cam_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_camdynamics_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_camdynamics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_camdynamics_workbench = true;
        run_camdynamics(&mut app.camdynamics);
        draw_workbench(&mut app);
    }
}
