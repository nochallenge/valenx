//! The right-side **Gearbox Workbench** panel — native two-stage compound
//! gear-train analysis over `valenx-gearbox`.
//!
//! Mirrors the Induction Motor / Heat Pump workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_gearbox_workbench`,
//! toggled from the View menu. The form sets the tooth counts of two
//! reduction stages, a per-stage mesh efficiency and the input
//! speed / torque; "Analyze" reports the overall ratio, efficiency,
//! output speed / torque and the input / output power, and "Show 3-D
//! gearbox" loads a representative two-stage gear train into the central
//! viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_gearbox::{CompoundTrain, GearStage};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Gearbox Workbench.
pub struct GearboxWorkbenchState {
    /// Stage 1 input (pinion) tooth count.
    stage1_in: u32,
    /// Stage 1 output (gear) tooth count.
    stage1_out: u32,
    /// Stage 2 input (pinion) tooth count.
    stage2_in: u32,
    /// Stage 2 output (gear) tooth count.
    stage2_out: u32,
    /// Per-stage mesh efficiency in (0, 1].
    efficiency: f64,
    /// Input shaft speed (rpm).
    input_speed_rpm: f64,
    /// Input shaft torque (N·m).
    input_torque_nm: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D gearbox solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for GearboxWorkbenchState {
    fn default() -> Self {
        // Two 17:51 (3:1) reduction stages at 97% mesh efficiency — an
        // overall 9:1 reducer at ~94% efficiency.
        Self {
            stage1_in: 17,
            stage1_out: 51,
            stage2_in: 17,
            stage2_out: 51,
            efficiency: 0.97,
            input_speed_rpm: 1500.0,
            input_torque_nm: 10.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Gearbox Workbench right-side panel. A no-op when the
/// `show_gearbox_workbench` toggle is off.
pub fn draw_gearbox_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_gearbox_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_gearbox_workbench",
        "Gearbox",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native two-stage compound gear-train analysis · valenx-gearbox",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.gearbox;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Stage 1  (pinion : gear)").strong());
                    ui.horizontal(|ui| {
                        ui.label("teeth in / out");
                        ui.add(egui::DragValue::new(&mut s.stage1_in).speed(1.0));
                        ui.add(egui::DragValue::new(&mut s.stage1_out).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Stage 2  (pinion : gear)").strong());
                    ui.horizontal(|ui| {
                        ui.label("teeth in / out");
                        ui.add(egui::DragValue::new(&mut s.stage2_in).speed(1.0));
                        ui.add(egui::DragValue::new(&mut s.stage2_out).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Mesh + input").strong());
                    ui.horizontal(|ui| {
                        ui.label("per-stage efficiency");
                        ui.add(egui::DragValue::new(&mut s.efficiency).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("input speed (rpm)");
                        ui.add(egui::DragValue::new(&mut s.input_speed_rpm).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("input torque (N·m)");
                        ui.add(egui::DragValue::new(&mut s.input_torque_nm).speed(0.5));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_gearbox(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D gearbox").strong())
                        .on_hover_text(
                            "Build a representative two-stage gear train (three shafts, two meshing gear pairs) as a 3-D solid and load it into the central viewport to orbit",
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
        app.show_gearbox_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.gearbox` borrow is
    // released here): build the gearbox's 3-D solid and load it.
    if app.gearbox.show_3d_request {
        app.gearbox.show_3d_request = false;
        load_gearbox_3d(app);
    }
}

/// Validate the form, evaluate the train and format the readout.
fn run_gearbox(s: &mut GearboxWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build a validated two-stage [`CompoundTrain`] from the form, mapping the
/// domain error to a display string.
fn build_train(s: &GearboxWorkbenchState) -> Result<CompoundTrain, String> {
    let s1 = GearStage::with_efficiency(s.stage1_in, s.stage1_out, s.efficiency)
        .map_err(|e| e.to_string())?;
    let s2 = GearStage::with_efficiency(s.stage2_in, s.stage2_out, s.efficiency)
        .map_err(|e| e.to_string())?;
    CompoundTrain::new(vec![s1, s2]).map_err(|e| e.to_string())
}

/// Evaluate the train and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &GearboxWorkbenchState) -> Result<String, String> {
    let r1 = s.stage1_out as f64 / s.stage1_in as f64;
    let r2 = s.stage2_out as f64 / s.stage2_in as f64;
    let train = build_train(s)?;

    let ratio = train.ratio();
    let eff = train.efficiency();
    let out_speed = train.output_speed(s.input_speed_rpm);
    let out_torque = train.output_torque(s.input_torque_nm);

    let w_in = s.input_speed_rpm * TAU / 60.0;
    let p_in = s.input_torque_nm * w_in;
    let w_out = out_speed * TAU / 60.0;
    let p_out = out_torque * w_out;

    Ok(format!(
        "stage 1         : {}/{}  (ratio {:.2})\n\
         stage 2         : {}/{}  (ratio {:.2})\n\
         overall ratio   : {:.2} : 1\n\
         overall eff.    : {:.1} %\n\n\
         input           : {:.1} rpm, {:.2} N·m\n\
         output          : {:.1} rpm, {:.2} N·m\n\n\
         input power     : {:.3} kW\n\
         output power    : {:.3} kW  (loss {:.3} kW)",
        s.stage1_in,
        s.stage1_out,
        r1,
        s.stage2_in,
        s.stage2_out,
        r2,
        ratio,
        eff * 100.0,
        s.input_speed_rpm,
        s.input_torque_nm,
        out_speed,
        out_torque,
        p_in / 1000.0,
        p_out / 1000.0,
        (p_in - p_out) / 1000.0,
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

/// Append a (double-sided) cylinder whose axis runs along `+x`, spanning
/// `base.x ..= base.x + length` with circle centre `(base.y, base.z)`.
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

/// Build the gearbox as a triangle [`Mesh`] — three parallel shafts and
/// two meshing gear pairs (a pinion driving a larger gear on each stage),
/// drawn as discs on the shafts, with a base. Representative geometry (the
/// gears are smooth discs sized by stage, not toothed; the ratios are the
/// `valenx-gearbox` result). `None` for an invalid train.
fn gearbox_solid_mesh(s: &GearboxWorkbenchState) -> Option<Mesh> {
    build_train(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    let pinion_r = 0.13;
    let gear_r = 0.32;
    let z = 0.7;
    // Three shafts at y = +0.45, 0, -0.45 (centre distance pinion + gear).
    for &y in &[0.45, 0.0, -0.45] {
        push_cyl_x(
            &mut nodes,
            &mut tris,
            Vector3::new(-0.6, y, z),
            1.2,
            0.04,
            12,
        );
    }
    // Stage 1: pinion on the top shaft meshing the gear on the mid shaft.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.05, 0.45, z),
        0.1,
        pinion_r,
        24,
    );
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.06, 0.0, z),
        0.12,
        gear_r,
        28,
    );
    // Stage 2: pinion on the mid shaft meshing the gear on the lower shaft.
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(0.3, 0.0, z),
        0.1,
        pinion_r,
        24,
    );
    push_cyl_x(
        &mut nodes,
        &mut tris,
        Vector3::new(0.29, -0.45, z),
        0.12,
        gear_r,
        28,
    );
    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.06),
        Vector3::new(0.7, 0.6, 0.06),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-gearbox");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D gearbox solid and load it into the central viewport.
fn load_gearbox_3d(app: &mut ValenxApp) {
    let Some(mesh) = gearbox_solid_mesh(&app.gearbox) else {
        app.gearbox.error =
            Some("gear-train parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<gearbox>/valenx-gearbox"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"gearbox"}`** product: the representative
/// two-stage compound gear train (three shafts, two meshing gear pairs, a
/// base) built from the canonical 9:1 reducer (two 17:51 stages at 97 %),
/// paired with the drive readout rows (ratio / output speed-torque / power),
/// at a fixed 3/4 camera. Registered in [`crate::products_registry`]; the
/// per-tool builder the registry dispatches to. Pure — driven off
/// [`GearboxWorkbenchState::default`].
pub(crate) fn gearbox_product() -> crate::WorkspaceProduct {
    let s = GearboxWorkbenchState::default();
    let mesh = gearbox_solid_mesh(&s).expect("canonical gearbox ⇒ gear-train solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<gearbox>/valenx-gearbox");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical gearbox ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Two-stage gearbox (9:1)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = GearboxWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_ratio_and_torque() {
        let mut s = GearboxWorkbenchState::default();
        run_gearbox(&mut s);
        assert!(
            s.error.is_none(),
            "default gearbox should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("overall ratio"));
        assert!(s.result.contains("output"));
        assert!(s.result.contains("power"));
        // Two 3:1 stages => 9:1 overall.
        assert!(s.result.contains("9.00 : 1"));
    }

    #[test]
    fn analyze_rejects_zero_teeth() {
        let mut s = GearboxWorkbenchState {
            stage1_in: 0,
            ..Default::default()
        };
        run_gearbox(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn overall_ratio_is_product_and_power_is_conserved() {
        // Ground truth: a compound train's ratio is the product of stage
        // ratios, and output power equals input power times efficiency.
        let s1 = GearStage::with_efficiency(17, 51, 0.97).unwrap();
        let s2 = GearStage::with_efficiency(17, 51, 0.97).unwrap();
        let train = CompoundTrain::new(vec![s1, s2]).unwrap();
        assert!((train.ratio() - 9.0).abs() < 1e-12);

        let (n_in, t_in) = (1500.0, 10.0);
        let p_in = t_in * (n_in * TAU / 60.0);
        let n_out = train.output_speed(n_in);
        let t_out = train.output_torque(t_in);
        let p_out = t_out * (n_out * TAU / 60.0);
        assert!((p_out - p_in * train.efficiency()).abs() < 1e-6);
    }

    #[test]
    fn gearbox_mesh_for_default_is_nonempty_and_in_range() {
        let s = GearboxWorkbenchState::default();
        let mesh = gearbox_solid_mesh(&s).expect("default gearbox yields a solid");
        assert!(mesh.nodes.len() > 8, "expected shafts + gears + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn gearbox_mesh_none_for_invalid() {
        let s = GearboxWorkbenchState {
            stage2_out: 0,
            ..Default::default()
        };
        assert!(gearbox_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_gearbox_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_gearbox_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_gearbox_workbench = true;
        run_gearbox(&mut app.gearbox);
        draw_workbench(&mut app);
    }
}
