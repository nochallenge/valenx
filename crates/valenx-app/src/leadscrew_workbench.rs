//! The right-side **Lead Screw Workbench** panel — native closed-form
//! power-screw kinematics and statics over `valenx-leadscrew`.
//!
//! Mirrors the Heat Transfer / DC Motor workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_leadscrew_workbench`,
//! toggled from the View menu. The form sets the screw geometry (lead and
//! pitch diameter), a screw speed, a microstep count, a nut/thread friction
//! coefficient and a load — chosen as either an applied screw torque (and the
//! thrust it produces) or a target axial thrust (and the torque it needs).
//! "Analyze" reports the lead angle, the linear feed, the per-microstep
//! resolution, the static torque/thrust pairing, the friction-derived screw
//! efficiency and the self-locking / back-drive classification, and
//! "Show 3-D" loads a representative threaded screw shaft into the central
//! viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_leadscrew::LeadScrew;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which static quantity the form drives the screw from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LoadMode {
    /// Apply a screw torque; report the axial thrust it produces.
    ThrustFromTorque,
    /// Target an axial thrust; report the screw torque it needs.
    TorqueFromThrust,
}

/// Persistent form + result state for the Lead Screw Workbench.
pub struct LeadscrewWorkbenchState {
    /// Lead — axial advance per revolution (mm/rev).
    lead_mm: f64,
    /// Pitch (mean) diameter of the thread (mm).
    pitch_diameter_mm: f64,
    /// Screw speed (rev/min).
    rpm: f64,
    /// Microsteps per revolution (drive setting).
    microsteps_per_rev: u32,
    /// Nut/thread coefficient of friction `mu`.
    friction_mu: f64,
    /// Which static quantity is the input.
    mode: LoadMode,
    /// Applied screw torque (N·mm), used in [`LoadMode::ThrustFromTorque`].
    torque_n_mm: f64,
    /// Target axial thrust (N), used in [`LoadMode::TorqueFromThrust`].
    thrust_n: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D screw shaft (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for LeadscrewWorkbenchState {
    fn default() -> Self {
        // A common "T8x8" leadscrew: 8 mm lead on an 8 mm pitch diameter,
        // a 200-step motor at 16x microstepping (3200 microsteps/rev),
        // steel-on-bronze mu ~ 0.2. At 8 mm lead on an 8 mm screw the lead
        // angle is steep (~17.7 deg, tan ~ 0.318 > 0.2) so it back-drives,
        // and 200 N·mm of torque at the friction-implied efficiency drives
        // a useful axial thrust.
        Self {
            lead_mm: 8.0,
            pitch_diameter_mm: 8.0,
            rpm: 300.0,
            microsteps_per_rev: 3200,
            friction_mu: 0.2,
            mode: LoadMode::ThrustFromTorque,
            torque_n_mm: 200.0,
            thrust_n: 500.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Lead Screw Workbench right-side panel. A no-op when the
/// `show_leadscrew_workbench` toggle is off.
pub fn draw_leadscrew_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_leadscrew_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_leadscrew_workbench",
        "Lead Screw",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native closed-form power-screw kinematics & statics · valenx-leadscrew",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.leadscrew;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Screw geometry").strong());
                    ui.horizontal(|ui| {
                        ui.label("lead (mm/rev)");
                        ui.add(egui::DragValue::new(&mut s.lead_mm).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("pitch diameter (mm)");
                        ui.add(egui::DragValue::new(&mut s.pitch_diameter_mm).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Drive").strong());
                    ui.horizontal(|ui| {
                        ui.label("speed (rev/min)");
                        ui.add(egui::DragValue::new(&mut s.rpm).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("microsteps/rev");
                        ui.add(egui::DragValue::new(&mut s.microsteps_per_rev).speed(50.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Friction").strong());
                    ui.horizontal(|ui| {
                        ui.label("coefficient μ");
                        ui.add(egui::DragValue::new(&mut s.friction_mu).speed(0.005));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Load").strong());
                    ui.radio_value(&mut s.mode, LoadMode::ThrustFromTorque, "thrust from torque");
                    ui.radio_value(&mut s.mode, LoadMode::TorqueFromThrust, "torque from thrust");
                    match s.mode {
                        LoadMode::ThrustFromTorque => {
                            ui.horizontal(|ui| {
                                ui.label("torque (N·mm)");
                                ui.add(egui::DragValue::new(&mut s.torque_n_mm).speed(5.0));
                            });
                        }
                        LoadMode::TorqueFromThrust => {
                            ui.horizontal(|ui| {
                                ui.label("thrust (N)");
                                ui.add(egui::DragValue::new(&mut s.thrust_n).speed(10.0));
                            });
                        }
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_leadscrew(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative threaded screw shaft as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Screw performance").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_leadscrew_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.leadscrew` borrow is
    // released here): build the screw shaft's 3-D solid and load it.
    if app.leadscrew.show_3d_request {
        app.leadscrew.show_3d_request = false;
        load_screw_3d(app);
    }
}

/// Validate the form, evaluate the screw and format the readout.
fn run_leadscrew(s: &mut LeadscrewWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the screw and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &LeadscrewWorkbenchState) -> Result<String, String> {
    let screw = LeadScrew::new(s.lead_mm, s.pitch_diameter_mm).map_err(|e| e.to_string())?;

    let lead_angle_deg = screw.lead_angle_deg();
    let feed_mm_per_min = screw
        .linear_speed_mm_per_min(s.rpm)
        .map_err(|e| e.to_string())?;
    let feed_mm_per_s = screw
        .linear_speed_mm_per_s(s.rpm)
        .map_err(|e| e.to_string())?;
    let resolution_mm = screw
        .resolution_mm(s.microsteps_per_rev)
        .map_err(|e| e.to_string())?;

    // Friction-derived ideal raising efficiency and self-locking class.
    let eta = screw
        .screw_efficiency(s.friction_mu)
        .map_err(|e| e.to_string())?;
    let mu_crit = screw.critical_friction();
    let bd = screw.back_drive(s.friction_mu).map_err(|e| e.to_string())?;
    let lock_label = if bd.self_locking {
        "self-locking"
    } else {
        "back-drivable"
    };
    // Angular margin to the self-locking threshold: phi - lambda. Positive
    // when self-locking (with headroom), zero at the boundary, negative when
    // back-drivable — a single "how far from locking" number.
    let locking_margin_deg = bd.locking_margin_rad().to_degrees();

    // The static torque/thrust pairing, using the friction-derived eta as the
    // lumped efficiency so the two directions stay mutually consistent.
    let (torque_n_mm, thrust_n) = match s.mode {
        LoadMode::ThrustFromTorque => {
            let f = screw
                .thrust_n(s.torque_n_mm, eta)
                .map_err(|e| e.to_string())?;
            (s.torque_n_mm, f)
        }
        LoadMode::TorqueFromThrust => {
            let t = screw
                .torque_for_thrust_n_mm(s.thrust_n, eta)
                .map_err(|e| e.to_string())?;
            (t, s.thrust_n)
        }
    };

    Ok(format!(
        "lead / pitch dia: {:.3} / {:.3} mm\n\
         lead angle λ    : {:.3} deg\n\
         speed           : {:.1} rev/min\n\
         feed            : {:.2} mm/min ({:.3} mm/s)\n\
         resolution      : {:.5} mm/microstep\n\n\
         friction μ      : {:.3}\n\
         critical μ      : {:.3} ({lock_label})\n\
         locking margin  : {:.3} deg (φ−λ)\n\
         screw efficiency: {:.1} %\n\n\
         torque          : {:.2} N·mm\n\
         thrust          : {:.2} N",
        s.lead_mm,
        s.pitch_diameter_mm,
        lead_angle_deg,
        s.rpm,
        feed_mm_per_min,
        feed_mm_per_s,
        resolution_mm,
        s.friction_mu,
        mu_crit,
        locking_margin_deg,
        eta * 100.0,
        torque_n_mm,
        thrust_n,
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

/// Build the screw as a triangle [`Mesh`] — a threaded shaft: a core
/// cylinder (root) with a helical thread ridge swept along it as a thicker
/// cylinder per turn. Representative geometry sized off the pitch diameter
/// (not to scale; the kinematic / static numbers are the `valenx-leadscrew`
/// result). `None` for an invalid configuration.
fn screw_solid_mesh(s: &LeadscrewWorkbenchState) -> Option<Mesh> {
    let screw = LeadScrew::new(s.lead_mm, s.pitch_diameter_mm).ok()?;

    // Normalize geometry against the pitch diameter so the shaft is always a
    // sensible on-screen size regardless of the absolute millimetre values.
    let d_m = screw.pitch_diameter_mm;
    let root_radius = 0.5;
    let crest_radius = 0.62;
    // Threads per shaft length follow the real lead/diameter ratio, clamped
    // to a drawable band so a very fine or very coarse screw still renders.
    let turns_per_dia = (d_m / s.lead_mm).clamp(2.0, 14.0);
    let turns = turns_per_dia.round().max(2.0) as usize;
    let half_len = 2.4;
    let length = 2.0 * half_len;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Core shaft (root cylinder).
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::zeros(),
        half_len,
        root_radius,
        28,
    );

    // Thread crests: a ridge ring at each turn, drawn as a short fat disk
    // (thin cylinder at the crest radius) to read as a helical thread.
    let ring_half = (length / turns as f64) * 0.18;
    for k in 0..turns {
        let frac = (k as f64 + 0.5) / turns as f64;
        let x = -half_len + frac * length;
        push_cyl_x(
            &mut nodes,
            &mut tris,
            Vector3::new(x, 0.0, 0.0),
            ring_half,
            crest_radius,
            28,
        );
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-leadscrew");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D screw shaft solid and load it into the central viewport.
fn load_screw_3d(app: &mut ValenxApp) {
    let Some(mesh) = screw_solid_mesh(&app.leadscrew) else {
        app.leadscrew.error =
            Some("screw parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<screw>/valenx-leadscrew"),
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
        let s = LeadscrewWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_all_quantities() {
        let mut s = LeadscrewWorkbenchState::default();
        run_leadscrew(&mut s);
        assert!(
            s.error.is_none(),
            "default screw should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("lead angle"));
        assert!(s.result.contains("feed"));
        assert!(s.result.contains("resolution"));
        assert!(s.result.contains("screw efficiency"));
        assert!(s.result.contains("thrust"));
        // T8x8 at 300 rev/min: feed = 8 * 300 = 2400 mm/min, 40 mm/s.
        assert!(s.result.contains("2400.00 mm/min"));
        assert!(s.result.contains("40.000 mm/s"));
        // 8 mm lead on an 8 mm screw is steep -> back-drivable at mu = 0.2.
        assert!(s.result.contains("back-drivable"));
    }

    #[test]
    fn analyze_rejects_zero_lead() {
        let mut s = LeadscrewWorkbenchState {
            lead_mm: 0.0,
            ..Default::default()
        };
        run_leadscrew(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_zero_microsteps() {
        let mut s = LeadscrewWorkbenchState {
            microsteps_per_rev: 0,
            ..Default::default()
        };
        run_leadscrew(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn feed_is_lead_times_rpm_ground_truth() {
        // Ground truth: the nut's linear feed equals lead * rpm exactly
        // (mm/rev * rev/min = mm/min), and per-second is that over 60.
        let lead = 5.0;
        let rpm = 240.0;
        let screw = LeadScrew::new(lead, 10.0).unwrap();
        let feed = screw.linear_speed_mm_per_min(rpm).unwrap();
        assert!((feed - lead * rpm).abs() < 1e-9, "feed = {feed}");
        let per_s = screw.linear_speed_mm_per_s(rpm).unwrap();
        assert!((per_s - feed / 60.0).abs() < 1e-9, "per_s = {per_s}");
    }

    #[test]
    fn locking_margin_matches_friction_minus_lead_angle_ground_truth() {
        use std::f64::consts::PI;
        // Ground truth: the locking margin is phi - lambda, with
        // phi = atan(mu) and lambda = atan(lead / (pi * d_m)). For the
        // default T8x8 screw (lead = 8, d_m = 8, mu = 0.2) the lead angle
        // (~17.657 deg) exceeds the friction angle (~11.310 deg), so the
        // margin is negative (back-drivable): atan(0.2) - atan(1/pi).
        let s = LeadscrewWorkbenchState::default();
        let lambda = (s.lead_mm / (PI * s.pitch_diameter_mm)).atan();
        let phi = s.friction_mu.atan();
        let expected_deg = (phi - lambda).to_degrees();
        // Hand value: -6.34685... deg, rounding to -6.347 at 3 dp.
        assert!(
            (expected_deg - (-6.346_854_677)).abs() < 1e-6,
            "hand-computed margin drifted: {expected_deg}"
        );
        let mut s = LeadscrewWorkbenchState::default();
        run_leadscrew(&mut s);
        assert!(
            s.error.is_none(),
            "default screw should analyze: {:?}",
            s.error
        );
        // The readout surfaces the margin to 3 dp in degrees.
        assert!(
            s.result.contains("locking margin  : -6.347 deg (φ−λ)"),
            "missing locking-margin line in: {}",
            s.result
        );
    }

    #[test]
    fn screw_mesh_for_default_is_nonempty_and_in_range() {
        let s = LeadscrewWorkbenchState::default();
        let mesh = screw_solid_mesh(&s).expect("default screw yields a solid");
        assert!(mesh.nodes.len() > 8, "expected shaft + thread rings");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn screw_mesh_none_for_invalid() {
        let s = LeadscrewWorkbenchState {
            pitch_diameter_mm: 0.0,
            ..Default::default()
        };
        assert!(screw_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_leadscrew_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_leadscrew_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_leadscrew_workbench = true;
        run_leadscrew(&mut app.leadscrew);
        draw_workbench(&mut app);
    }
}
