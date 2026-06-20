//! The right-side **Brake Workbench** panel — native disc / caliper brake
//! friction-torque and stopping analysis over `valenx-brake`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_brake_workbench`,
//! toggled from the View menu. The form sets a disc / caliper brake (clamp
//! force, friction coefficient, pad-face count and the annular pad radii)
//! plus a stopping scenario (mass, speed, deceleration); "Analyze" derives
//! the uniform-wear effective radius, the Coulomb friction torque
//! `T = mu * F * n_pads * r_eff`, the larger uniform-pressure (fresh-pad)
//! effective radius and its higher initial torque, and the kinetic energy
//! / stopping distance / stopping time / average braking power of the
//! stop, and
//! "Show 3-D brake" loads a representative disc-plus-caliper solid into
//! the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_brake::disc::{
    disc_torque, effective_radius_uniform_pressure, effective_radius_uniform_wear, friction_force,
};
use valenx_brake::energy::{
    average_braking_power, kinetic_energy, stopping_distance, stopping_time,
};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Brake Workbench.
pub struct BrakeWorkbenchState {
    /// Friction coefficient `mu` of the pad/rotor pair.
    mu: f64,
    /// Normal clamp force per pad face `F` (N).
    clamp_force_n: f64,
    /// Number of friction faces (a floating single-piston caliper has 2).
    n_pads: u32,
    /// Inner pad radius (m).
    r_inner_m: f64,
    /// Outer pad radius (m).
    r_outer_m: f64,
    /// Vehicle mass for the stopping scenario (kg).
    mass_kg: f64,
    /// Initial speed for the stopping scenario (m/s).
    speed_mps: f64,
    /// Constant deceleration magnitude for the stopping scenario (m/s²).
    decel_mps2: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D brake solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for BrakeWorkbenchState {
    fn default() -> Self {
        // A 2-pad caliper on a 100->140 mm pad annulus (r_eff = 120 mm),
        // 8 kN clamp, mu = 0.4: T = 0.4*8000*2*0.12 = 768 N·m. Stopping a
        // 1500 kg car from 30 m/s at 7.5 m/s^2: 675 kJ, 60 m, 4 s.
        Self {
            mu: 0.4,
            clamp_force_n: 8_000.0,
            n_pads: 2,
            r_inner_m: 0.10,
            r_outer_m: 0.14,
            mass_kg: 1_500.0,
            speed_mps: 30.0,
            decel_mps2: 7.5,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Brake Workbench right-side panel. A no-op when the
/// `show_brake_workbench` toggle is off.
pub fn draw_brake_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_brake_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_brake_workbench",
        "Brake",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native disc / caliper friction torque + stopping · valenx-brake",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.brake;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Caliper").strong());
                    ui.horizontal(|ui| {
                        ui.label("friction μ");
                        ui.add(egui::DragValue::new(&mut s.mu).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("clamp force F (N)");
                        ui.add(egui::DragValue::new(&mut s.clamp_force_n).speed(50.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("pad faces");
                        ui.add(egui::DragValue::new(&mut s.n_pads).speed(1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Pad annulus").strong());
                    ui.horizontal(|ui| {
                        ui.label("inner radius (m)");
                        ui.add(egui::DragValue::new(&mut s.r_inner_m).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("outer radius (m)");
                        ui.add(egui::DragValue::new(&mut s.r_outer_m).speed(0.005));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Stop").strong());
                    ui.horizontal(|ui| {
                        ui.label("mass (kg)");
                        ui.add(egui::DragValue::new(&mut s.mass_kg).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("speed (m/s)");
                        ui.add(egui::DragValue::new(&mut s.speed_mps).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("deceleration (m/s²)");
                        ui.add(egui::DragValue::new(&mut s.decel_mps2).speed(0.1));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_brake(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D brake").strong())
                        .on_hover_text(
                            "Build a representative disc rotor with a straddling caliper and hub as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Braking").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_brake_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.brake` borrow is
    // released here): build the brake's 3-D solid and load it.
    if app.brake.show_3d_request {
        app.brake.show_3d_request = false;
        load_brake_3d(app);
    }
}

/// Validate the form, evaluate the brake and format the readout.
fn run_brake(s: &mut BrakeWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The uniform-wear effective radius and the friction torque
/// `(r_eff, torque)`, the quantities both the readout and the 3-D gate
/// need. Extracted so it is unit-testable and shared.
fn torque(s: &BrakeWorkbenchState) -> Result<(f64, f64), String> {
    let r_eff =
        effective_radius_uniform_wear(s.r_inner_m, s.r_outer_m).map_err(|e| e.to_string())?;
    let t = disc_torque(s.mu, s.clamp_force_n, s.n_pads, r_eff).map_err(|e| e.to_string())?;
    Ok((r_eff, t))
}

/// The **fresh-pad** (uniform-pressure) effective radius and the friction
/// torque it produces `(r_eff, torque)`. A new, not-yet-bedded pad
/// presses with uniform contact pressure, giving a larger effective
/// radius — `r_eff = (2/3)(r_o^3 - r_i^3)/(r_o^2 - r_i^2)` — and hence a
/// higher initial brake torque than the worn-in `torque` value. The
/// pair quantifies the bedding-in bite spread. Extracted so it is
/// unit-testable.
fn torque_fresh(s: &BrakeWorkbenchState) -> Result<(f64, f64), String> {
    let r_eff =
        effective_radius_uniform_pressure(s.r_inner_m, s.r_outer_m).map_err(|e| e.to_string())?;
    let t = disc_torque(s.mu, s.clamp_force_n, s.n_pads, r_eff).map_err(|e| e.to_string())?;
    Ok((r_eff, t))
}

/// Evaluate the brake and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &BrakeWorkbenchState) -> Result<String, String> {
    let (r_eff, t) = torque(s)?;
    let (r_eff_fresh, t_fresh) = torque_fresh(s)?;
    let f_face = friction_force(s.mu, s.clamp_force_n).map_err(|e| e.to_string())?;
    let e_kin = kinetic_energy(s.mass_kg, s.speed_mps).map_err(|e| e.to_string())?;
    let dist = stopping_distance(s.speed_mps, s.decel_mps2).map_err(|e| e.to_string())?;
    let time = stopping_time(s.speed_mps, s.decel_mps2).map_err(|e| e.to_string())?;
    let power = average_braking_power(s.mass_kg, s.speed_mps, time).map_err(|e| e.to_string())?;

    Ok(format!(
        "friction μ      : {:.3}\n\
         clamp force F   : {:.0} N\n\
         pad faces       : {}\n\
         pad annulus     : {:.3} -> {:.3} m\n\
         r_eff (wear)    : {:.4} m\n\
         friction / face : {:.1} N\n\
         brake torque T  : {:.2} N·m\n\
         r_eff (fresh)   : {:.4} m\n\
         torque (fresh)  : {:.2} N·m\n\n\
         mass / speed    : {:.0} kg / {:.1} m/s\n\
         kinetic energy  : {:.2} kJ\n\
         stopping dist   : {:.2} m\n\
         stopping time   : {:.2} s\n\
         avg brake power : {:.2} kW",
        s.mu,
        s.clamp_force_n,
        s.n_pads,
        s.r_inner_m,
        s.r_outer_m,
        r_eff,
        f_face,
        t,
        r_eff_fresh,
        t_fresh,
        s.mass_kg,
        s.speed_mps,
        e_kin / 1_000.0,
        dist,
        time,
        power / 1_000.0,
    ))
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

/// Append a double-sided cylinder whose axis runs along `x` (centre `c`,
/// radius `r`, half-length `hx`, `seg` facets) to the buffers — the disc
/// rotor and hub barrel.
fn push_x_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    r: f64,
    hx: f64,
    seg: usize,
) {
    let base = nodes.len();
    // Two rings (-x then +x face) plus the two face centres.
    for k in 0..seg {
        let a = TAU * (k as f64) / (seg as f64);
        let (y, z) = (r * a.cos(), r * a.sin());
        nodes.push(c + Vector3::new(-hx, y, z));
    }
    for k in 0..seg {
        let a = TAU * (k as f64) / (seg as f64);
        let (y, z) = (r * a.cos(), r * a.sin());
        nodes.push(c + Vector3::new(hx, y, z));
    }
    let c_lo = nodes.len();
    nodes.push(c + Vector3::new(-hx, 0.0, 0.0));
    let c_hi = nodes.len();
    nodes.push(c + Vector3::new(hx, 0.0, 0.0));

    for k in 0..seg {
        let kn = (k + 1) % seg;
        let lo0 = base + k;
        let lo1 = base + kn;
        let hi0 = base + seg + k;
        let hi1 = base + seg + kn;
        // Lateral wall (two triangles, both windings for a visible shell).
        tris.extend_from_slice(&[lo0, hi0, hi1, lo0, hi1, lo1]);
        // End caps.
        tris.extend_from_slice(&[c_lo, lo1, lo0]);
        tris.extend_from_slice(&[c_hi, hi0, hi1]);
    }
}

/// Build the disc brake as a triangle [`Mesh`] — a disc rotor (a fat
/// `x`-axis cylinder), a small hub barrel on the axis, and a caliper box
/// straddling the rotor rim. Representative geometry (not to scale; the
/// torque / stopping numbers are the `valenx-brake` result). `None` for
/// an invalid configuration.
fn brake_solid_mesh(s: &BrakeWorkbenchState) -> Option<Mesh> {
    // Gate the 3-D solid on a buildable brake object.
    torque(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Disc rotor: thin in x, radius scaled from the outer pad radius.
    let r_disc = (s.r_outer_m / 0.14).clamp(0.4, 3.0);
    push_x_cylinder(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.0),
        r_disc,
        0.05,
        48,
    );
    // Hub barrel on the axis.
    push_x_cylinder(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.0),
        0.18 * r_disc,
        0.16,
        24,
    );
    // Caliper: a box straddling the rotor rim near the top of the disc.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.86 * r_disc, 0.0),
        Vector3::new(0.13, 0.12 * r_disc, 0.10),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-brake");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D brake solid and load it into the central viewport.
fn load_brake_3d(app: &mut ValenxApp) {
    let Some(mesh) = brake_solid_mesh(&app.brake) else {
        app.brake.error = Some("brake parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<brake>/valenx-brake"),
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
        let s = BrakeWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_torque_and_stop() {
        let mut s = BrakeWorkbenchState::default();
        run_brake(&mut s);
        assert!(
            s.error.is_none(),
            "default brake should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("brake torque T"));
        assert!(s.result.contains("stopping dist"));
        assert!(s.result.contains("avg brake power"));
        // T = 0.4*8000*2*0.12 = 768.00 N·m.
        assert!(s.result.contains("768.00"));
        // d = 30^2 / (2*7.5) = 60.00 m; t = 30/7.5 = 4.00 s.
        assert!(s.result.contains("60.00"));
        assert!(s.result.contains("4.00"));
    }

    #[test]
    fn analyze_rejects_zero_friction() {
        let mut s = BrakeWorkbenchState {
            mu: 0.0,
            ..Default::default()
        };
        run_brake(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn torque_matches_mu_f_n_reff_ground_truth() {
        // Ground truth: a uniform-wear effective radius is the mean of the
        // pad radii, and T = mu * F * n_pads * r_eff exactly.
        let s = BrakeWorkbenchState::default();
        let (r_eff, t) = torque(&s).unwrap();
        let expected_r = 0.5 * (s.r_inner_m + s.r_outer_m);
        assert!((r_eff - expected_r).abs() < 1e-12, "r_eff {r_eff}");
        let expected_t = s.mu * s.clamp_force_n * f64::from(s.n_pads) * expected_r;
        assert!((t - expected_t).abs() < 1e-9, "torque {t} vs {expected_t}");
        assert!((t - 768.0).abs() < 1e-9, "torque {t}");
    }

    #[test]
    fn fresh_pad_torque_matches_uniform_pressure_ground_truth() {
        // Ground truth: a fresh (uniform-pressure) pad has the LARGER
        // effective radius r_eff = (2/3)(r_o^3 - r_i^3)/(r_o^2 - r_i^2),
        // which for the 0.10 -> 0.14 m annulus is
        //   (2/3)(0.14^3 - 0.10^3)/(0.14^2 - 0.10^2)
        // = (2/3)(0.001744/0.0096) = 0.121111... m, exceeding the
        // uniform-wear 0.12 m by exactly (r_o - r_i)^2 / (6(r_o + r_i)).
        let s = BrakeWorkbenchState::default();
        let (r_eff_fresh, t_fresh) = torque_fresh(&s).unwrap();

        let num = s.r_outer_m.powi(3) - s.r_inner_m.powi(3);
        let den = s.r_outer_m.powi(2) - s.r_inner_m.powi(2);
        let expected_r = (2.0 / 3.0) * num / den;
        assert!(
            (r_eff_fresh - expected_r).abs() < 1e-12,
            "r_eff_fresh {r_eff_fresh} vs {expected_r}"
        );

        // The fresh radius exceeds the worn radius by the analytic gap.
        let r_wear = 0.5 * (s.r_inner_m + s.r_outer_m);
        let gap = (s.r_outer_m - s.r_inner_m).powi(2) / (6.0 * (s.r_outer_m + s.r_inner_m));
        assert!(
            (r_eff_fresh - (r_wear + gap)).abs() < 1e-12,
            "fresh {r_eff_fresh} vs wear+gap {}",
            r_wear + gap
        );

        // Torque T = mu * F * n_pads * r_eff_fresh = 0.4*8000*2*0.121111..
        // = 775.111.. N·m, larger than the worn-in 768 N·m.
        let expected_t = s.mu * s.clamp_force_n * f64::from(s.n_pads) * expected_r;
        assert!(
            (t_fresh - expected_t).abs() < 1e-9,
            "t_fresh {t_fresh} vs {expected_t}"
        );
        assert!(
            t_fresh > 768.0,
            "fresh torque {t_fresh} should exceed worn 768"
        );

        // The formatted readout surfaces both fresh-pad lines.
        let out = compute(&s).unwrap();
        assert!(
            out.contains("r_eff (fresh)"),
            "missing fresh radius line: {out}"
        );
        assert!(
            out.contains("torque (fresh)"),
            "missing fresh torque line: {out}"
        );
        // r_eff (fresh) : 0.1211 m  and  torque (fresh) : 775.11 N·m.
        assert!(out.contains("0.1211"), "fresh r_eff substring: {out}");
        assert!(out.contains("775.11"), "fresh torque substring: {out}");
    }

    #[test]
    fn brake_mesh_for_default_is_nonempty_and_in_range() {
        let s = BrakeWorkbenchState::default();
        let mesh = brake_solid_mesh(&s).expect("default brake yields a solid");
        assert!(mesh.nodes.len() > 8, "expected rotor + hub + caliper");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn brake_mesh_none_for_invalid() {
        let s = BrakeWorkbenchState {
            r_outer_m: 0.05,
            r_inner_m: 0.10,
            ..Default::default()
        };
        // Degenerate annulus (outer <= inner) -> no effective radius.
        assert!(brake_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_brake_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_brake_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_brake_workbench = true;
        run_brake(&mut app.brake);
        draw_workbench(&mut app);
    }
}
