//! The right-side **Clutch Workbench** panel — native dry plate / disc
//! friction-clutch torque-capacity sizing over `valenx-clutch`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_clutch_workbench`,
//! toggled from the View menu. The form sets a friction face (inner /
//! outer radius), the coefficient of friction, the number of friction
//! surfaces in contact `N`, the axial clamp force and the engaged speed,
//! and picks a contact-pressure idealisation ([`PressureModel`]).
//! "Analyze" reports the effective lever-arm radius, the torque capacity
//! under both the uniform-wear and uniform-pressure theories and the
//! transmissible power at speed, and "Show 3-D clutch" loads a
//! representative friction-disc-plus-pressure-plate solid into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_clutch::{rpm_to_rad_per_s, ClutchGeometry, FrictionClutch, PressureModel};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use std::f64::consts::TAU;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Clutch Workbench.
pub struct ClutchWorkbenchState {
    /// Inner radius `ri` of the friction annulus (mm).
    inner_radius_mm: f64,
    /// Outer radius `ro` of the friction annulus (mm).
    outer_radius_mm: f64,
    /// Coefficient of (kinetic) friction `mu` between the rubbing faces.
    mu: f64,
    /// Number of friction surfaces in contact `N` (single plate => 2).
    surfaces: u32,
    /// Axial clamp force `F` (N).
    clamp_force_n: f64,
    /// Engaged shaft speed (rev/min) used for the transmissible power.
    speed_rpm: f64,
    /// Which contact-pressure idealisation to headline.
    model: PressureModel,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D clutch solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for ClutchWorkbenchState {
    fn default() -> Self {
        // A single-plate dry clutch: 100/200 mm face, mu = 0.3, both
        // faces grip (N = 2), 5 kN clamp, engaged at 3000 rpm. Uniform
        // wear acts at the 150 mm mean radius => T = 0.3*5000*2*0.150 =
        // 450 N*m, and at 3000 rpm that carries ~141 kW.
        Self {
            inner_radius_mm: 100.0,
            outer_radius_mm: 200.0,
            mu: 0.3,
            surfaces: 2,
            clamp_force_n: 5000.0,
            speed_rpm: 3000.0,
            model: PressureModel::UniformWear,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Clutch Workbench right-side panel. A no-op when the
/// `show_clutch_workbench` toggle is off.
pub fn draw_clutch_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_clutch_workbench {
        return;
    }

    egui::SidePanel::right("valenx_clutch_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Clutch",
                "native dry friction-clutch torque capacity · valenx-clutch",
            ) {
                app.show_clutch_workbench = false;
            }

            let s = &mut app.clutch;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Friction face").strong());
                    ui.horizontal(|ui| {
                        ui.label("inner radius ri (mm)");
                        ui.add(egui::DragValue::new(&mut s.inner_radius_mm).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("outer radius ro (mm)");
                        ui.add(egui::DragValue::new(&mut s.outer_radius_mm).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("friction coeff μ");
                        ui.add(egui::DragValue::new(&mut s.mu).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("surfaces N");
                        ui.add(egui::DragValue::new(&mut s.surfaces).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Operating point").strong());
                    ui.horizontal(|ui| {
                        ui.label("clamp force F (N)");
                        ui.add(egui::DragValue::new(&mut s.clamp_force_n).speed(25.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("speed (rpm)");
                        ui.add(egui::DragValue::new(&mut s.speed_rpm).speed(10.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Pressure model").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.model, PressureModel::UniformWear, "uniform wear");
                        ui.radio_value(
                            &mut s.model,
                            PressureModel::UniformPressure,
                            "uniform pressure",
                        );
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_clutch(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D clutch").strong())
                        .on_hover_text(
                            "Build a representative friction disc (an annular plate with a centre hole) plus a thin pressure plate as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Torque capacity").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });

    // Serviced after the panel draws (the `&mut app.clutch` borrow is
    // released here): build the clutch's 3-D solid and load it.
    if app.clutch.show_3d_request {
        app.clutch.show_3d_request = false;
        load_clutch_3d(app);
    }
}

/// Validate the form, evaluate the clutch and format the readout.
fn run_clutch(s: &mut ClutchWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`FrictionClutch`] from the form, the object both
/// the readout and the 3-D gate need. Extracted so it is unit-testable
/// and shared.
fn build_clutch(s: &ClutchWorkbenchState) -> Result<FrictionClutch, String> {
    let geom =
        ClutchGeometry::new(s.inner_radius_mm, s.outer_radius_mm).map_err(|e| e.to_string())?;
    FrictionClutch::new(geom, s.mu, s.surfaces).map_err(|e| e.to_string())
}

/// Evaluate the clutch and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &ClutchWorkbenchState) -> Result<String, String> {
    let clutch = build_clutch(s)?;
    let geom = clutch.geometry();

    let t_wear = clutch
        .torque(PressureModel::UniformWear, s.clamp_force_n)
        .map_err(|e| e.to_string())?;
    let t_pres = clutch
        .torque(PressureModel::UniformPressure, s.clamp_force_n)
        .map_err(|e| e.to_string())?;

    let omega = rpm_to_rad_per_s(s.speed_rpm).map_err(|e| e.to_string())?;
    let power_w = clutch
        .power(s.model, s.clamp_force_n, omega)
        .map_err(|e| e.to_string())?;

    // Both theories' effective lever-arm radii. The uniform-wear value is
    // the arithmetic mean (ro + ri)/2; the uniform-pressure value is the
    // larger area-weighted (centroidal) mean — which is why uniform
    // pressure always predicts the higher torque for the same clamp force.
    let r_wear_m = geom.mean_radius_uniform_wear_m();
    let r_pres_m = geom.mean_radius_uniform_pressure_m();

    let (model_label, t_model) = match s.model {
        PressureModel::UniformWear => ("uniform wear", t_wear),
        PressureModel::UniformPressure => ("uniform pressure", t_pres),
    };

    Ok(format!(
        "face ri / ro    : {:.1} / {:.1} mm\n\
         friction μ      : {:.3}\n\
         surfaces N      : {}\n\
         clamp force F   : {:.0} N\n\
         speed           : {:.0} rpm ({:.1} rad/s)\n\n\
         model           : {model_label}\n\
         r_eff wear      : {:.2} mm\n\
         r_eff pressure  : {:.2} mm\n\
         T uniform wear  : {:.1} N·m\n\
         T uniform press : {:.1} N·m\n\
         T (model)       : {:.1} N·m\n\
         power at speed  : {:.2} kW",
        s.inner_radius_mm,
        s.outer_radius_mm,
        s.mu,
        s.surfaces,
        s.clamp_force_n,
        s.speed_rpm,
        omega,
        r_wear_m * 1.0e3,
        r_pres_m * 1.0e3,
        t_wear,
        t_pres,
        t_model,
        power_w * 1.0e-3,
    ))
}

/// Append an annular ring (a disc with a centre hole) of half-thickness
/// `half_t` centred on the `z` axis at height `z0`, spanning radii
/// `[ri, ro]`, as a double-sided strip of `seg` segments. The top and
/// bottom annulus faces plus the inner and outer cylindrical walls.
fn push_ring(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    z0: f64,
    half_t: f64,
    ri: f64,
    ro: f64,
    seg: usize,
) {
    let base = nodes.len();
    let zt = z0 + half_t;
    let zb = z0 - half_t;
    for i in 0..seg {
        let a = TAU * (i as f64) / (seg as f64);
        let (ca, sa) = (a.cos(), a.sin());
        // 0: outer-top, 1: inner-top, 2: outer-bottom, 3: inner-bottom.
        nodes.push(Vector3::new(ro * ca, ro * sa, zt));
        nodes.push(Vector3::new(ri * ca, ri * sa, zt));
        nodes.push(Vector3::new(ro * ca, ro * sa, zb));
        nodes.push(Vector3::new(ri * ca, ri * sa, zb));
    }
    for i in 0..seg {
        let a = base + 4 * i;
        let b = base + 4 * ((i + 1) % seg);
        let (ot0, it0, ob0, ib0) = (a, a + 1, a + 2, a + 3);
        let (ot1, it1, ob1, ib1) = (b, b + 1, b + 2, b + 3);
        // Top annulus face.
        tris.extend_from_slice(&[ot0, ot1, it1, ot0, it1, it0]);
        // Bottom annulus face.
        tris.extend_from_slice(&[ob0, ib0, ib1, ob0, ib1, ob1]);
        // Outer cylindrical wall.
        tris.extend_from_slice(&[ot0, ob0, ob1, ot0, ob1, ot1]);
        // Inner cylindrical wall (the centre hole).
        tris.extend_from_slice(&[it0, it1, ib1, it0, ib1, ib0]);
    }
}

/// Build the clutch as a triangle [`Mesh`] — the annular friction disc
/// (a plate with a centre hole, drawn at the form's `[ri, ro]`) plus a
/// thin pressure plate behind it. Representative geometry (the torque /
/// power numbers are the `valenx-clutch` result). `None` for an invalid
/// configuration.
fn clutch_solid_mesh(s: &ClutchWorkbenchState) -> Option<Mesh> {
    let clutch = build_clutch(s).ok()?;
    let geom = clutch.geometry();
    let ri = geom.inner_radius_m();
    let ro = geom.outer_radius_m();

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // The friction disc: the true annulus, a slim plate at the origin.
    push_ring(&mut nodes, &mut tris, 0.0, 0.06 * ro, ri, ro, 64);
    // A thin full-face pressure plate just behind it (small centre bore).
    push_ring(
        &mut nodes,
        &mut tris,
        -0.16 * ro,
        0.03 * ro,
        0.25 * ri,
        ro,
        64,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-clutch");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D clutch solid and load it into the central viewport.
fn load_clutch_3d(app: &mut ValenxApp) {
    let Some(mesh) = clutch_solid_mesh(&app.clutch) else {
        app.clutch.error =
            Some("clutch parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<clutch>/valenx-clutch"),
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
        let s = ClutchWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_torque_and_power() {
        let mut s = ClutchWorkbenchState::default();
        run_clutch(&mut s);
        assert!(
            s.error.is_none(),
            "default clutch should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("T uniform wear"));
        assert!(s.result.contains("power at speed"));
        // 100/200 mm, mu=0.3, N=2, F=5 kN => uniform-wear T = 450.0 N*m.
        assert!(s.result.contains("450.0"));
        // ...which at 3000 rpm carries 141.37 kW.
        assert!(s.result.contains("141.37"));
    }

    #[test]
    fn analyze_rejects_inverted_radii() {
        let mut s = ClutchWorkbenchState {
            inner_radius_mm: 200.0,
            outer_radius_mm: 100.0,
            ..Default::default()
        };
        run_clutch(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn uniform_wear_mean_radius_torque_ground_truth() {
        // Ground truth (Shigley uniform-wear): the clamp force acts at
        // the arithmetic mean radius, so
        //   T = mu * F * N * (ro + ri) / 2.
        // For mu=0.3, F=5000 N, N=2, ri=0.100 m, ro=0.200 m this is
        //   0.3 * 5000 * 2 * (0.200 + 0.100)/2 = 450.0 N*m exactly.
        let mu = 0.3_f64;
        let f = 5000.0_f64;
        let n = 2.0_f64;
        let ri = 0.100_f64;
        let ro = 0.200_f64;
        let expected = mu * f * n * (ro + ri) / 2.0;
        assert!((expected - 450.0).abs() < 1e-9);

        let s = ClutchWorkbenchState::default();
        let clutch = build_clutch(&s).expect("default clutch builds");
        let t = clutch
            .torque(PressureModel::UniformWear, s.clamp_force_n)
            .expect("uniform-wear torque");
        assert!((t - expected).abs() < 1e-9, "got {t}");
    }

    #[test]
    fn analyze_default_reports_both_lever_arm_radii() {
        // The readout now surfaces BOTH effective lever-arm radii (the
        // crate's mean_radius_uniform_wear_m / mean_radius_uniform_pressure_m).
        // Ground truth for the 100/200 mm default face:
        //   r_eff wear     = (ro + ri)/2 = (0.200 + 0.100)/2 = 0.150 m = 150.00 mm.
        //   r_eff pressure = (2/3)(ro^3 - ri^3)/(ro^2 - ri^2)
        //                  = (2/3)(0.008 - 0.001)/(0.04 - 0.01)
        //                  = (2/3)(0.007/0.03) = 0.155555... m = 155.56 mm.
        let ri = 0.100_f64;
        let ro = 0.200_f64;
        let r_wear = 0.5 * (ro + ri);
        let r_pres = (2.0 / 3.0) * (ro * ro * ro - ri * ri * ri) / (ro * ro - ri * ri);
        assert!((r_wear - 0.150).abs() < 1e-12, "wear lever {r_wear}");
        assert!(
            (r_pres * 1.0e3 - 155.5556).abs() < 1e-3,
            "pressure lever {r_pres}"
        );
        // The area-weighted mean must exceed the arithmetic mean.
        assert!(r_pres > r_wear);

        let mut s = ClutchWorkbenchState::default();
        run_clutch(&mut s);
        assert!(
            s.error.is_none(),
            "default clutch should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("r_eff wear      : 150.00 mm"));
        assert!(s.result.contains("r_eff pressure  : 155.56 mm"));
    }

    #[test]
    fn clutch_mesh_for_default_is_nonempty_and_in_range() {
        let s = ClutchWorkbenchState::default();
        let mesh = clutch_solid_mesh(&s).expect("default clutch yields a solid");
        assert!(mesh.nodes.len() > 8, "expected disc + pressure plate");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn clutch_mesh_none_for_invalid() {
        let s = ClutchWorkbenchState {
            inner_radius_mm: 200.0,
            outer_radius_mm: 100.0,
            ..Default::default()
        };
        assert!(clutch_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_clutch_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_clutch_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_clutch_workbench = true;
        run_clutch(&mut app.clutch);
        draw_workbench(&mut app);
    }
}
