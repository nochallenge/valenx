//! The right-side **Chain Drive Workbench** panel — native single-stage
//! roller-chain kinematics over `valenx-chaindrive`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_chaindrive_workbench`,
//! toggled from the View menu. The form sets a sprocket pair (driver and
//! driven tooth counts, chain pitch) plus an operating point (input speed,
//! input torque, shaft centre distance); "Analyze" runs
//! [`valenx_chaindrive::analyze`] and reports the speed ratio, chain
//! velocity, driven-shaft speed, loss-free output torque and the buildable
//! chain length, and "Show 3-D drive" loads a representative two-sprocket
//! solid into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_chaindrive::{analyze, SprocketPair};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

use std::f64::consts::TAU;

/// Display units for the reported chain (link) velocity. A UI-only toggle:
/// it changes how the velocity line is formatted, never the underlying
/// `valenx-chaindrive` computation.
#[derive(Clone, Copy, PartialEq)]
pub enum VelocityUnits {
    /// Metres per second (the crate's native unit).
    MetresPerSecond,
    /// Feet per minute (a common conveyor / drive-shop unit).
    FeetPerMinute,
}

/// Persistent form + result state for the Chain Drive Workbench.
pub struct ChainDriveWorkbenchState {
    /// Driver (input) sprocket tooth count.
    driver_teeth: u32,
    /// Driven (output) sprocket tooth count.
    driven_teeth: u32,
    /// Chain pitch (roller-to-roller spacing), mm.
    pitch_mm: f64,
    /// Input (driver-shaft) rotational speed, rev/min.
    input_rpm: f64,
    /// Input (driver-shaft) torque, N·m.
    input_torque_n_m: f64,
    /// Shaft centre distance, mm.
    center_distance_mm: f64,
    /// Units in which to display the chain velocity.
    velocity_units: VelocityUnits,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D drive solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for ChainDriveWorkbenchState {
    fn default() -> Self {
        // 17-tooth driver, 34-tooth driven, ANSI 40 (12.7 mm) chain at
        // 1000 rpm / 50 N·m with 500 mm between shafts: a clean 2:1
        // reduction -> 500 rpm out, 100 N·m out, ~3.60 m/s chain speed.
        Self {
            driver_teeth: 17,
            driven_teeth: 34,
            pitch_mm: 12.7,
            input_rpm: 1000.0,
            input_torque_n_m: 50.0,
            center_distance_mm: 500.0,
            velocity_units: VelocityUnits::MetresPerSecond,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Chain Drive Workbench right-side panel. A no-op when the
/// `show_chaindrive_workbench` toggle is off.
pub fn draw_chaindrive_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_chaindrive_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_chaindrive_workbench",
        "Chain Drive",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native single-stage roller-chain kinematics · valenx-chaindrive",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.chaindrive;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Sprockets").strong());
                    ui.horizontal(|ui| {
                        ui.label("driver teeth");
                        ui.add(egui::DragValue::new(&mut s.driver_teeth).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("driven teeth");
                        ui.add(egui::DragValue::new(&mut s.driven_teeth).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("chain pitch (mm)");
                        ui.add(egui::DragValue::new(&mut s.pitch_mm).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Operating point").strong());
                    ui.horizontal(|ui| {
                        ui.label("input speed (rpm)");
                        ui.add(egui::DragValue::new(&mut s.input_rpm).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("input torque (N·m)");
                        ui.add(egui::DragValue::new(&mut s.input_torque_n_m).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("centre distance (mm)");
                        ui.add(egui::DragValue::new(&mut s.center_distance_mm).speed(5.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Chain-speed units").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(
                            &mut s.velocity_units,
                            VelocityUnits::MetresPerSecond,
                            "m/s",
                        );
                        ui.radio_value(
                            &mut s.velocity_units,
                            VelocityUnits::FeetPerMinute,
                            "ft/min",
                        );
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_chaindrive(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D drive").strong())
                        .on_hover_text(
                            "Build a representative driver + driven sprocket pair with the two chain spans as a 3-D solid and load it into the central viewport to orbit",
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
        app.show_chaindrive_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.chaindrive` borrow is
    // released here): build the drive's 3-D solid and load it.
    if app.chaindrive.show_3d_request {
        app.chaindrive.show_3d_request = false;
        load_drive_3d(app);
    }
}

/// Validate the form, evaluate the drive and format the readout.
fn run_chaindrive(s: &mut ChainDriveWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the chain drive and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &ChainDriveWorkbenchState) -> Result<String, String> {
    let pair =
        SprocketPair::new(s.driver_teeth, s.driven_teeth, s.pitch_mm).map_err(|e| e.to_string())?;
    let r = analyze(&pair, s.input_rpm, s.input_torque_n_m, s.center_distance_mm)
        .map_err(|e| e.to_string())?;

    let (v_value, v_unit) = match s.velocity_units {
        VelocityUnits::MetresPerSecond => (r.chain_velocity_m_per_s, "m/s"),
        // 1 m/s = 196.850393... ft/min (60 / 0.3048).
        VelocityUnits::FeetPerMinute => (r.chain_velocity_m_per_s * 60.0 / 0.3048, "ft/min"),
    };

    // Sprocket pitch (reference) circle diameters d = p / sin(π / z): the
    // circle the chain rollers seat on, and hence the lever arm that turns
    // chain tension into shaft torque. Sized straight off the validated
    // pair, so they are always available once `analyze` above succeeded.
    let driver_pd_mm = pair.driver_pitch_diameter_mm();
    let driven_pd_mm = pair.driven_pitch_diameter_mm();

    Ok(format!(
        "driver / driven : {} / {} teeth\n\
         chain pitch     : {:.2} mm\n\
         input speed     : {:.1} rpm\n\
         input torque    : {:.2} N·m\n\
         centre distance : {:.1} mm\n\n\
         ratio (N2/N1)   : {:.4}\n\
         driver pitch dia: {driver_pd_mm:.2} mm\n\
         driven pitch dia: {driven_pd_mm:.2} mm\n\
         chain velocity  : {:.3} {}\n\
         driven speed    : {:.2} rpm\n\
         output torque   : {:.2} N·m\n\
         chain length    : {} pitches",
        s.driver_teeth,
        s.driven_teeth,
        s.pitch_mm,
        s.input_rpm,
        s.input_torque_n_m,
        s.center_distance_mm,
        r.ratio,
        v_value,
        v_unit,
        r.driven_speed_rpm,
        r.output_torque_n_m,
        r.chain_length_pitches,
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

/// Append a double-sided z-axis disc (a sprocket face) of `radius` at
/// height `z`, centred on `(cx, 0)`, as a triangle fan with `seg` segments.
fn push_disc(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    cx: f64,
    z: f64,
    radius: f64,
    seg: usize,
) {
    let center = nodes.len();
    nodes.push(Vector3::new(cx, 0.0, z));
    for i in 0..seg {
        let a = TAU * i as f64 / seg as f64;
        nodes.push(Vector3::new(cx + radius * a.cos(), radius * a.sin(), z));
    }
    for i in 0..seg {
        let a = center + 1 + i;
        let b = center + 1 + (i + 1) % seg;
        // Double-sided so the sprocket face shows from either orbit angle.
        tris.extend_from_slice(&[center, a, b, center, b, a]);
    }
}

/// Build the chain drive as a triangle [`Mesh`] — a driver disc and a
/// (larger / smaller) driven disc, separated along `x` by the centre
/// distance, joined by the two straight chain spans (thin boxes for the
/// taut and slack runs). Representative geometry (radii are the real pitch
/// radii, scaled to metres; the kinematic numbers are the
/// `valenx-chaindrive` result). `None` for an invalid configuration.
fn drive_solid_mesh(s: &ChainDriveWorkbenchState) -> Option<Mesh> {
    // Gate on the real crate object: invalid sprockets / pitch -> no solid.
    let pair = SprocketPair::new(s.driver_teeth, s.driven_teeth, s.pitch_mm).ok()?;

    // Pitch radii in metres (mm / 1000 / 2), and the centre distance in
    // metres laid out along x with the driver at the origin.
    let r_driver = 0.5 * pair.driver_pitch_diameter_mm() / 1000.0;
    let r_driven = 0.5 * pair.driven_pitch_diameter_mm() / 1000.0;
    let cx_driven = (s.center_distance_mm / 1000.0).max(r_driver + r_driven);

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let seg = 48;
    // Driver sprocket: two faces giving a thin disc.
    push_disc(&mut nodes, &mut tris, 0.0, -0.02, r_driver, seg);
    push_disc(&mut nodes, &mut tris, 0.0, 0.02, r_driver, seg);
    // Driven sprocket.
    push_disc(&mut nodes, &mut tris, cx_driven, -0.02, r_driven, seg);
    push_disc(&mut nodes, &mut tris, cx_driven, 0.02, r_driven, seg);

    // Two straight chain spans (taut / slack runs) as thin boxes bridging
    // the two pitch circles at +y and -y of the smaller sprocket.
    let span_half_x = 0.5 * cx_driven;
    let mid_x = 0.5 * cx_driven;
    let span_y = r_driver.min(r_driven);
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(mid_x, span_y, 0.0),
        Vector3::new(span_half_x, 0.01, 0.015),
    );
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(mid_x, -span_y, 0.0),
        Vector3::new(span_half_x, 0.01, 0.015),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-chaindrive");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D drive solid and load it into the central viewport.
fn load_drive_3d(app: &mut ValenxApp) {
    let Some(mesh) = drive_solid_mesh(&app.chaindrive) else {
        app.chaindrive.error =
            Some("drive parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<drive>/valenx-chaindrive"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical chain-drive workbench as a 3-D solid plus
/// its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn chaindrive_product() -> crate::WorkspaceProduct {
    let s = ChainDriveWorkbenchState::default();
    let mesh = drive_solid_mesh(&s).expect("canonical chain drive ⇒ drive solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<chaindrive>/valenx-chaindrive");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical chain drive ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Chain drive (ratio/power)".into(),
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

    #[test]
    fn default_state_is_idle() {
        let s = ChainDriveWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_ratio_speed_and_torque() {
        let mut s = ChainDriveWorkbenchState::default();
        run_chaindrive(&mut s);
        assert!(
            s.error.is_none(),
            "default drive should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("ratio (N2/N1)"));
        assert!(s.result.contains("chain velocity"));
        assert!(s.result.contains("output torque"));
        // 17 -> 34 teeth is a clean 2:1 reduction.
        assert!(s.result.contains("2.0000"));
        // 1000 rpm in / 2 -> 500 rpm out.
        assert!(s.result.contains("500.00 rpm"));
        // 50 N·m in * 2 -> 100 N·m out.
        assert!(s.result.contains("100.00 N·m"));
    }

    #[test]
    fn analyze_velocity_units_toggle_changes_output() {
        let mut s = ChainDriveWorkbenchState::default();
        run_chaindrive(&mut s);
        // m/s by default: v = 12.7 * 17 * 1000 / 60000 = 3.598... m/s.
        assert!(s.result.contains("3.598 m/s"));

        s.velocity_units = VelocityUnits::FeetPerMinute;
        run_chaindrive(&mut s);
        // 3.598333 m/s * 60 / 0.3048 = 708.333... ft/min.
        assert!(s.result.contains("ft/min"));
        assert!(s.result.contains("708.333 ft/min"), "{}", s.result);
    }

    #[test]
    fn analyze_rejects_too_few_teeth() {
        let mut s = ChainDriveWorkbenchState {
            driver_teeth: 3,
            ..Default::default()
        };
        run_chaindrive(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn velocity_ratio_equals_driven_over_driver_teeth() {
        // Ground truth: the speed ratio of a chain drive is exactly the
        // tooth-count ratio N_driven / N_driver, and the driven speed is
        // the input divided by that ratio. Pin the f64 inputs.
        let driver_teeth: f64 = 17.0;
        let driven_teeth: f64 = 34.0;
        let input_rpm: f64 = 1000.0;
        let pair = SprocketPair::new(17, 34, 12.7).unwrap();
        let r = analyze(&pair, input_rpm, 50.0, 500.0).unwrap();
        let expected_ratio = driven_teeth / driver_teeth;
        assert!((r.ratio - expected_ratio).abs() < 1e-12);
        assert!((r.driven_speed_rpm - input_rpm / expected_ratio).abs() < 1e-9);
    }

    #[test]
    fn readout_reports_sprocket_pitch_diameters() {
        // Ground truth: a sprocket's pitch (reference) circle diameter is
        // the closed form d = p / sin(π / z). For the default pair
        // (p = 12.7 mm) the driver has z = 17 and the driven z = 34:
        //   driver d = 12.7 / sin(π/17) = 69.1158... mm  -> "69.12 mm"
        //   driven d = 12.7 / sin(π/34) = 137.6420... mm -> "137.64 mm"
        use std::f64::consts::PI;
        let s = ChainDriveWorkbenchState::default();
        let out = compute(&s).expect("default drive computes");

        let driver_expected = 12.7 / (PI / 17.0).sin();
        let driven_expected = 12.7 / (PI / 34.0).sin();
        // Pin the hand value the formatted substrings below round to.
        assert!(
            (driver_expected - 69.115_827).abs() < 1e-5,
            "driver pitch dia hand-calc drifted: {driver_expected}"
        );
        assert!(
            (driven_expected - 137.641_983).abs() < 1e-5,
            "driven pitch dia hand-calc drifted: {driven_expected}"
        );

        // The new readout lines must show those diameters at :.2f.
        assert!(
            out.contains("driver pitch dia: 69.12 mm"),
            "missing driver pitch diameter line:\n{out}"
        );
        assert!(
            out.contains("driven pitch dia: 137.64 mm"),
            "missing driven pitch diameter line:\n{out}"
        );
    }

    #[test]
    fn drive_mesh_for_default_is_nonempty_and_in_range() {
        let s = ChainDriveWorkbenchState::default();
        let mesh = drive_solid_mesh(&s).expect("default drive yields a solid");
        assert!(mesh.nodes.len() > 8, "expected two discs + two spans");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn drive_mesh_none_for_invalid() {
        let s = ChainDriveWorkbenchState {
            driver_teeth: 3,
            ..Default::default()
        };
        assert!(drive_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_chaindrive_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_chaindrive_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_chaindrive_workbench = true;
        run_chaindrive(&mut app.chaindrive);
        draw_workbench(&mut app);
    }
}
