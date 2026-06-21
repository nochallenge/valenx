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
/// Presentation spin rate of the input (stage-1) pinion, rad/s (~1.3 rev/s) —
/// a readable inspect speed; the rest follow the true ratios + mesh signs.
const INPUT_RAD_PER_S: f32 = 8.0;

/// Build the gearbox as a triangle [`Mesh`] together with a [`crate::RigidPart`]
/// per gear disc, so the four gears spin about their shafts (the +x axis through
/// each shaft) while the three cosmetic shaft rods and the base stay put.
/// Meshing pairs counter-rotate and the stage-2 pinion shares the mid shaft with
/// the stage-1 gear (equal ω); rates follow the disc-radius ratio,
/// presentation-scaled. `None` for an invalid train.
fn gearbox_solid_mesh_parts(s: &GearboxWorkbenchState) -> Option<(Mesh, Vec<crate::RigidPart>)> {
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

    // Record each gear disc's node range + shaft pivot as we fuse it, so the
    // animation rotates each gear about its own shaft (+x). Order: stage-1
    // pinion, stage-1 gear, stage-2 pinion, stage-2 gear.
    let mut gear_ranges: Vec<(std::ops::Range<usize>, [f32; 3])> = Vec::with_capacity(4);
    let mut push_gear = |nodes: &mut Vec<Vector3<f64>>,
                         tris: &mut Vec<usize>,
                         centre: Vector3<f64>,
                         half_len: f64,
                         radius: f64| {
        let start = nodes.len();
        push_cyl_x(nodes, tris, centre, half_len, radius, 24);
        let end = nodes.len();
        // Pivot is the shaft line: x is irrelevant for +x rotation, so use
        // the gear's (y, z) shaft centre.
        gear_ranges.push((start..end, [0.0, centre.y as f32, centre.z as f32]));
    };

    // Stage 1: pinion on the top shaft meshing the gear on the mid shaft.
    push_gear(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.05, 0.45, z),
        0.1,
        pinion_r,
    );
    push_gear(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.06, 0.0, z),
        0.12,
        gear_r,
    );
    // Stage 2: pinion on the mid shaft meshing the gear on the lower shaft.
    push_gear(
        &mut nodes,
        &mut tris,
        Vector3::new(0.3, 0.0, z),
        0.1,
        pinion_r,
    );
    push_gear(
        &mut nodes,
        &mut tris,
        Vector3::new(0.29, -0.45, z),
        0.12,
        gear_r,
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

    // True mesh kinematics on the representative discs (ratio = pinion_r/gear_r
    // per stage): the input pinion spins +; the stage-1 gear counter-rotates and
    // is slower; the stage-2 pinion shares the mid shaft (equal ω); the stage-2
    // gear counter-rotates the stage-2 pinion (sign flips back +) and is slower.
    let step = (pinion_r / gear_r) as f32; // <1 per stage
    let w_in = INPUT_RAD_PER_S;
    let w_s1g = -w_in * step; // stage-1 gear (counter, slower)
    let w_s2p = w_s1g; // stage-2 pinion rigid on the mid shaft
    let w_s2g = -w_s2p * step; // stage-2 gear (counter again ⇒ +, slower)
    let omegas = [w_in, w_s1g, w_s2p, w_s2g];
    let parts: Vec<crate::RigidPart> = gear_ranges
        .into_iter()
        .zip(omegas)
        .map(|((node_range, pivot), rad_per_s)| crate::RigidPart {
            node_range,
            axis: [1.0, 0.0, 0.0],
            pivot,
            rad_per_s,
        })
        .collect();
    Some((mesh, parts))
}

/// Build the gearbox [`Mesh`] (without the per-gear part metadata) for the
/// central viewport. See [`gearbox_solid_mesh_parts`].
fn gearbox_solid_mesh(s: &GearboxWorkbenchState) -> Option<Mesh> {
    gearbox_solid_mesh_parts(s).map(|(mesh, _parts)| mesh)
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
    let (mesh, parts) =
        gearbox_solid_mesh_parts(&s).expect("canonical gearbox ⇒ gear-train solid builds");
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
        image: None,
        image_texture: None,
        // Animated: each gear disc spins about its shaft, meshing pairs counter-
        // rotate, while the shaft rods and base stay put. Paused at t = 0.
        animation: Some(crate::ProductAnimation {
            playing: false,
            speed: 1.0,
            t: 0.0,
            motion: crate::ProductMotion::RigidParts(parts),
        }),
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

    #[test]
    fn gearbox_product_spins_four_gears_counter_rotating() {
        // The product carries a RigidParts animation: four gear discs, each a
        // non-empty node range within the mesh, all spinning about +x. Meshing
        // pairs counter-rotate; the stage-2 pinion shares the mid shaft with the
        // stage-1 gear (equal ω). The shaft rods + base are left static.
        let product = gearbox_product();
        let loaded = product.mesh.as_ref().expect("gearbox product has a mesh");
        let node_count = loaded.mesh.nodes.len();
        let anim = product.animation.expect("gearbox product is animated");
        assert!(!anim.playing, "starts paused");
        match anim.motion {
            crate::ProductMotion::RigidParts(parts) => {
                assert_eq!(parts.len(), 4, "four gear discs");
                for p in &parts {
                    assert!(
                        p.node_range.start < p.node_range.end,
                        "non-empty gear range"
                    );
                    assert!(p.node_range.end <= node_count, "gear range within the mesh");
                    assert_eq!(p.axis, [1.0, 0.0, 0.0], "spins about its shaft");
                    assert!(p.rad_per_s.abs() > 0.0, "non-zero spin rate");
                }
                let w_in = parts[0].rad_per_s; // stage-1 pinion (input)
                let w_s1g = parts[1].rad_per_s; // stage-1 gear
                let w_s2p = parts[2].rad_per_s; // stage-2 pinion (mid shaft)
                let w_s2g = parts[3].rad_per_s; // stage-2 gear (output)
                assert!(w_in > 0.0, "input pinion spins +");
                assert!(w_s1g < 0.0, "stage-1 gear counter-rotates the input");
                assert!(
                    (w_s1g - w_s2p).abs() < 1e-6,
                    "stage-2 pinion shares the mid shaft with the stage-1 gear"
                );
                assert!(
                    w_s2g > 0.0,
                    "stage-2 gear counter-rotates the pinion (sign flips back +)"
                );
                // The train steps down: the input is fastest, the output slowest.
                assert!(w_in.abs() > w_s1g.abs());
                assert!(w_s2p.abs() > w_s2g.abs());
                // The first gear is built right after the three shaft rods (static).
                assert!(
                    parts[0].node_range.start > 0,
                    "shaft rods precede the gears"
                );
            }
            crate::ProductMotion::Turntable { .. } => {
                panic!("gearbox must use per-part rigid motion")
            }
        }
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
