//! The right-side **Belt Drive Workbench** panel — native open flat-belt
//! drive analysis over `valenx-beltdrive`.
//!
//! Mirrors the Gearbox / Heat Transfer workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_beltdrive_workbench`,
//! toggled from the View menu. The form sets the driver and driven pulley
//! diameters, the driver speed, the belt/pulley friction coefficient, the
//! shaft centre distance and the belt linear density / maximum tension;
//! "Analyze" reports the speed ratio, belt linear speed, driven speed, the
//! open-belt wrap angles, the capstan tension ratio `T1/T2 = exp(mu*theta)`,
//! the centrifugal tension and the slipping-limited power, and "Show 3-D
//! belt drive" loads a representative two-pulley drive (two discs on
//! parallel shafts joined by a belt loop) into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_beltdrive::{belt_length_open, DriveAnalysis, DriveSpec};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Belt Drive Workbench.
pub struct BeltDriveWorkbenchState {
    /// Driver (input) pulley pitch diameter (mm).
    driver_diameter_mm: f64,
    /// Driven (output) pulley pitch diameter (mm).
    driven_diameter_mm: f64,
    /// Driver rotational speed (rpm).
    driver_speed_rpm: f64,
    /// Belt/pulley coefficient of friction (dimensionless).
    mu: f64,
    /// Shaft-to-shaft centre distance (mm).
    center_distance_mm: f64,
    /// Belt linear mass density (kg/m).
    linear_density: f64,
    /// Maximum allowable tight-side tension (N).
    t1_max: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D belt-drive solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for BeltDriveWorkbenchState {
    fn default() -> Self {
        // A 100 mm driver at 1450 rpm driving a 250 mm pulley 500 mm away
        // with mu = 0.3: a 2.5:1 reduction, ~7.6 m/s belt, capstan ratio
        // ~2.35 on the small pulley.
        Self {
            driver_diameter_mm: 100.0,
            driven_diameter_mm: 250.0,
            driver_speed_rpm: 1450.0,
            mu: 0.3,
            center_distance_mm: 500.0,
            linear_density: 0.4,
            t1_max: 1200.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Belt Drive Workbench right-side panel. A no-op when the
/// `show_beltdrive_workbench` toggle is off.
pub fn draw_beltdrive_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_beltdrive_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_beltdrive_workbench",
        "Belt Drive",
        |app, ui| {
            ui.label(
                egui::RichText::new("native open flat-belt drive analysis · valenx-beltdrive")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.beltdrive;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Pulleys").strong());
                    ui.horizontal(|ui| {
                        ui.label("driver Ø (mm)");
                        ui.add(egui::DragValue::new(&mut s.driver_diameter_mm).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("driven Ø (mm)");
                        ui.add(egui::DragValue::new(&mut s.driven_diameter_mm).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("centre distance (mm)");
                        ui.add(egui::DragValue::new(&mut s.center_distance_mm).speed(5.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Drive").strong());
                    ui.horizontal(|ui| {
                        ui.label("driver speed (rpm)");
                        ui.add(egui::DragValue::new(&mut s.driver_speed_rpm).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("friction μ");
                        ui.add(egui::DragValue::new(&mut s.mu).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Belt").strong());
                    ui.horizontal(|ui| {
                        ui.label("linear density (kg/m)");
                        ui.add(egui::DragValue::new(&mut s.linear_density).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("max tension T1 (N)");
                        ui.add(egui::DragValue::new(&mut s.t1_max).speed(10.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_beltdrive(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D belt drive").strong())
                        .on_hover_text(
                            "Build a representative open belt drive (driver and driven discs on two parallel shafts joined by a belt loop) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Drive").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_beltdrive_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.beltdrive` borrow is
    // released here): build the belt drive's 3-D solid and load it.
    if app.beltdrive.show_3d_request {
        app.beltdrive.show_3d_request = false;
        load_beltdrive_3d(app);
    }
}

/// Validate the form, evaluate the drive and format the readout.
fn run_beltdrive(s: &mut BeltDriveWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build a validated [`DriveSpec`] from the form (converting the mm pulley
/// dimensions to metres and the rpm to rev/s), mapping any domain error to
/// a display string. Extracted so it is shared by the readout and the 3-D
/// gate.
fn build_spec(s: &BeltDriveWorkbenchState) -> Result<DriveSpec, String> {
    Ok(DriveSpec {
        driver_diameter: s.driver_diameter_mm / 1000.0,
        driven_diameter: s.driven_diameter_mm / 1000.0,
        driver_rev_per_sec: s.driver_speed_rpm / 60.0,
        center_distance: s.center_distance_mm / 1000.0,
        mu: s.mu,
        linear_density: s.linear_density,
        t1_max: s.t1_max,
    })
}

/// Evaluate the drive and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &BeltDriveWorkbenchState) -> Result<String, String> {
    let spec = build_spec(s)?;
    let a: DriveAnalysis = spec.analyze().map_err(|e| e.to_string())?;

    let driven_rpm = a.driven_rev_per_sec * 60.0;
    // Exact open-belt length (two straight tangents + the two wrapped
    // arcs): the practical "what length of belt to order" figure. The
    // free function sorts the small/large radii internally, so the
    // driver/driven order does not matter.
    let belt_length = belt_length_open(
        spec.driver_diameter / 2.0,
        spec.driven_diameter / 2.0,
        spec.center_distance,
    )
    .map_err(|e| e.to_string())?;

    Ok(format!(
        "driver / driven : {:.0} / {:.0} mm\n\
         centre distance : {:.0} mm\n\
         friction μ      : {:.2}\n\n\
         speed ratio     : {:.2} : 1\n\
         belt speed      : {:.2} m/s\n\
         driver / driven : {:.0} / {:.0} rpm\n\n\
         wrap (small)    : {:.1}°\n\
         wrap (large)    : {:.1}°\n\
         belt length     : {belt_length:.3} m ({:.0} mm)\n\n\
         capstan T1/T2   : {:.3}\n\
         centrifugal Tc  : {:.1} N\n\
         max power       : {:.3} kW",
        s.driver_diameter_mm,
        s.driven_diameter_mm,
        s.center_distance_mm,
        s.mu,
        a.speed_ratio,
        a.belt_speed,
        s.driver_speed_rpm,
        driven_rpm,
        a.wrap_small.to_degrees(),
        a.wrap_large.to_degrees(),
        belt_length * 1000.0,
        a.tension_ratio,
        a.centrifugal_tension,
        a.max_power / 1000.0,
    ))
}

/// Append an outward-facing box (centre `c`, half-extents `h`) to the
/// buffers (double-sided).
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

/// Append a (double-sided) cylinder disc whose axis runs along `+x`,
/// spanning `base.x ..= base.x + length` with circle centre
/// `(base.y, base.z)`.
fn push_cyl_x(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    length: f64,
    r: f64,
    seg: usize,
) {
    let (x0, x1) = (base.x, base.x + length);
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x0, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x1, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
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
    }
}

/// Build the belt drive as a triangle [`Mesh`] — two pulley discs on two
/// parallel shafts (the driver small, the driven large), joined by a belt
/// loop drawn as the two straight tangent runs (top and bottom). The two
/// shaft centres are spaced along `z` by the centre distance; the disc
/// radii follow the diameters. Representative geometry (the kinematics /
/// power numbers are the `valenx-beltdrive` result). `None` for an invalid
/// configuration.
fn beltdrive_solid_mesh(s: &BeltDriveWorkbenchState) -> Option<Mesh> {
    // Gate the geometry on a buildable, analysable drive.
    build_spec(s).ok()?.analyze().ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Scale the real (mm) dimensions to a tidy on-screen size.
    let scale = 1.0 / s.center_distance_mm.max(1.0);
    let r_drv = 0.5 * s.driver_diameter_mm * scale;
    let r_drn = 0.5 * s.driven_diameter_mm * scale;
    let half_c = 0.5; // centre distance maps to 1.0 along z
    let face_x = -0.06;
    let width = 0.12;
    let lift = (r_drv.max(r_drn)) + 0.15;

    // Driver pulley disc (small) at z = -half_c.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(face_x, lift, -half_c),
        width,
        r_drv,
        28,
    );
    // Driven pulley disc (large) at z = +half_c.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(face_x, lift, half_c),
        width,
        r_drn,
        32,
    );

    // Belt loop: two straight tangent runs as thin boxes spanning the gap.
    // Approximate the (near-horizontal) tangents by the top and bottom of
    // each pulley; representative, not the exact tangent line.
    let belt_x = face_x + 0.5 * width;
    let belt_t = 0.02; // belt half-thickness
    let mid_z = 0.0;
    let span_z = half_c + 0.5 * belt_t;
    // Top run, sitting on the larger radius.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(belt_x, lift + r_drn, mid_z),
        Vector3::new(0.5 * width, belt_t, span_z),
    );
    // Bottom run, hanging below by the larger radius.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(belt_x, lift - r_drn, mid_z),
        Vector3::new(0.5 * width, belt_t, span_z),
    );

    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.06, 0.0),
        Vector3::new(0.5, 0.06, 0.8),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-beltdrive");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D belt-drive solid and load it into the central viewport.
fn load_beltdrive_3d(app: &mut ValenxApp) {
    let Some(mesh) = beltdrive_solid_mesh(&app.beltdrive) else {
        app.beltdrive.error =
            Some("belt-drive parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<belt-drive>/valenx-beltdrive"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical belt-drive workbench as a 3-D solid plus
/// its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn beltdrive_product() -> crate::WorkspaceProduct {
    let s = BeltDriveWorkbenchState::default();
    let mesh = beltdrive_solid_mesh(&s).expect("canonical belt drive ⇒ drive solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<beltdrive>/valenx-beltdrive");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical belt drive ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Belt drive (tension/power)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn default_state_is_idle() {
        let s = BeltDriveWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_ratio_speed_and_tension() {
        let mut s = BeltDriveWorkbenchState::default();
        run_beltdrive(&mut s);
        assert!(
            s.error.is_none(),
            "default belt drive should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("speed ratio"));
        assert!(s.result.contains("belt speed"));
        assert!(s.result.contains("capstan T1/T2"));
        // 250 mm driven on a 100 mm driver => 2.5:1 reduction.
        assert!(s.result.contains("2.50 : 1"));
        // Belt rim speed pi * 0.1 m * (1450/60) rev/s ~ 7.59 m/s.
        assert!(s.result.contains("7.59 m/s"));
    }

    #[test]
    fn analyze_rejects_zero_diameter() {
        let mut s = BeltDriveWorkbenchState {
            driver_diameter_mm: 0.0,
            ..Default::default()
        };
        run_beltdrive(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn speed_ratio_is_driven_over_driver_ground_truth() {
        // Ground truth: with no slip the surface speeds match, so the
        // transmission ratio is exactly D_driven / D_driver, and the belt
        // rim speed is pi * D_driver * N_driver.
        let s = BeltDriveWorkbenchState::default();
        let a = build_spec(&s).unwrap().analyze().unwrap();
        let expected_ratio = s.driven_diameter_mm / s.driver_diameter_mm;
        assert!((a.speed_ratio - expected_ratio).abs() < 1e-12);
        let expected_v = PI * (s.driver_diameter_mm / 1000.0) * (s.driver_speed_rpm / 60.0);
        assert!((a.belt_speed - expected_v).abs() < 1e-9);
        // Capstan ratio on the small pulley: T1/T2 = exp(mu * theta_small).
        assert!((a.tension_ratio - (s.mu * a.wrap_small).exp()).abs() < 1e-9);
    }

    #[test]
    fn belt_length_is_exact_open_belt_ground_truth() {
        // Ground truth for the default 100/250 mm pulleys 500 mm apart:
        // the exact open-belt length is the two straight tangents plus the
        // two wrapped arcs,
        //   alpha = asin((R_large - R_small) / C),
        //   L = 2*C*cos(alpha) + R_small*(pi - 2*alpha) + R_large*(pi + 2*alpha).
        // With R_small = 0.05 m, R_large = 0.125 m, C = 0.5 m this is
        // 0.98868 + 0.14202 + 0.43034 = 1.56105 m, i.e. "1.561 m".
        let s = BeltDriveWorkbenchState::default();
        let r_small = s.driver_diameter_mm / 1000.0 / 2.0;
        let r_large = s.driven_diameter_mm / 1000.0 / 2.0;
        let c = s.center_distance_mm / 1000.0;
        let alpha = ((r_large - r_small) / c).asin();
        let expected =
            2.0 * c * alpha.cos() + r_small * (PI - 2.0 * alpha) + r_large * (PI + 2.0 * alpha);
        assert!(
            (expected - 1.561_049_951_958_976).abs() < 1e-12,
            "hand check drifted: {expected}"
        );

        // The crate's free function must reproduce that closed form.
        let l = belt_length_open(r_small, r_large, c).unwrap();
        assert!((l - expected).abs() < 1e-12, "belt_length_open: {l}");

        // …and the readout must surface it (3-decimal metres).
        let out = compute(&s).expect("default belt drive computes");
        assert!(
            out.contains("belt length     : 1.561 m"),
            "missing belt length line: {out}"
        );
        assert!(out.contains("1561 mm"), "missing mm form: {out}");
    }

    #[test]
    fn beltdrive_mesh_for_default_is_nonempty_and_in_range() {
        let s = BeltDriveWorkbenchState::default();
        let mesh = beltdrive_solid_mesh(&s).expect("default belt drive yields a solid");
        assert!(mesh.nodes.len() > 8, "expected two pulleys + belt + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn beltdrive_mesh_none_for_invalid() {
        // Centre distance too small for the radius difference => degenerate
        // open-belt geometry, so no 3-D solid.
        let s = BeltDriveWorkbenchState {
            center_distance_mm: 10.0,
            ..Default::default()
        };
        assert!(beltdrive_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_beltdrive_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_beltdrive_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_beltdrive_workbench = true;
        run_beltdrive(&mut app.beltdrive);
        draw_workbench(&mut app);
    }
}
