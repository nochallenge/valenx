//! The right-side **MOSFET Workbench** panel — native square-law
//! (Shockley level-1) n-channel MOSFET IV analysis over `valenx-mosfet`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_mosfet_workbench`,
//! toggled from the View menu. The form sets the device transconductance
//! parameter `k` and threshold `vth`, plus a gate / drain bias
//! `(vgs, vds)` and a target operating [`Region`]; "Analyze" classifies
//! the region, evaluates the drain current `Id` and transconductance `gm`,
//! and reports the gate bias that would carry the present saturation
//! current, and "Show 3-D" loads a representative layered MOS device
//! (gate / oxide / channel / substrate) solid into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_mosfet::{Mosfet, Region};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the MOSFET Workbench.
pub struct MosfetWorkbenchState {
    /// Transconductance parameter `k = μ_n · C_ox · W / L` (A/V²).
    k_a_per_v2: f64,
    /// Threshold voltage `vth` (V).
    vth_v: f64,
    /// Gate-to-source bias `vgs` (V).
    vgs_v: f64,
    /// Drain-to-source bias `vds` (V).
    vds_v: f64,
    /// Target operating region the design should land in (compared with
    /// the region the present bias actually produces).
    target_region: Region,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D device solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for MosfetWorkbenchState {
    fn default() -> Self {
        // A textbook NMOS (k = 0.5 mA/V^2, vth = 1.0 V) biased at
        // vgs = 3 V (overdrive 2 V) and vds = 5 V >= overdrive, so it sits
        // in saturation: Id = 0.5 * k * vov^2 = 1.0 mA, gm = k*vov = 1.0 mS.
        Self {
            k_a_per_v2: 0.5e-3,
            vth_v: 1.0,
            vgs_v: 3.0,
            vds_v: 5.0,
            target_region: Region::Saturation,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the MOSFET Workbench right-side panel. A no-op when the
/// `show_mosfet_workbench` toggle is off.
pub fn draw_mosfet_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_mosfet_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_mosfet_workbench",
        "MOSFET",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native square-law NMOS IV (cutoff/triode/saturation) · valenx-mosfet",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.mosfet;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Device").strong());
                    ui.horizontal(|ui| {
                        let cap_k = ui.label("k (A/V²)");
                        ui.add(egui::DragValue::new(&mut s.k_a_per_v2).speed(1.0e-4))
                            .labelled_by(cap_k.id);
                    });
                    ui.horizontal(|ui| {
                        let cap_vth = ui.label("threshold vth (V)");
                        ui.add(egui::DragValue::new(&mut s.vth_v).speed(0.05))
                            .labelled_by(cap_vth.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Bias").strong());
                    ui.horizontal(|ui| {
                        let cap_vgs = ui.label("gate vgs (V)");
                        ui.add(egui::DragValue::new(&mut s.vgs_v).speed(0.05))
                            .labelled_by(cap_vgs.id);
                    });
                    ui.horizontal(|ui| {
                        let cap_vds = ui.label("drain vds (V)");
                        ui.add(egui::DragValue::new(&mut s.vds_v).speed(0.05))
                            .labelled_by(cap_vds.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Target region").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.target_region, Region::Cutoff, "cutoff");
                        ui.radio_value(&mut s.target_region, Region::Triode, "triode");
                        ui.radio_value(&mut s.target_region, Region::Saturation, "saturation");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_mosfet(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative layered MOS device (gate / oxide / channel / substrate stack) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Operating point").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_mosfet_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.mosfet` borrow is
    // released here): build the device's 3-D solid and load it.
    if app.mosfet.show_3d_request {
        app.mosfet.show_3d_request = false;
        load_device_3d(app);
    }
}

/// Validate the form, evaluate the device and format the readout.
fn run_mosfet(s: &mut MosfetWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Mosfet`] from the form, the quantity both the
/// readout and the 3-D gate need. Extracted so it is unit-testable and
/// shared.
fn device(s: &MosfetWorkbenchState) -> Result<Mosfet, String> {
    Mosfet::new(s.k_a_per_v2, s.vth_v).map_err(|e| e.to_string())
}

/// Evaluate the device at the form bias and format the full readout,
/// mapping any domain error to a display string. Extracted so it is
/// unit-testable.
fn compute(s: &MosfetWorkbenchState) -> Result<String, String> {
    let m = device(s)?;
    let op = m
        .operating_point(s.vgs_v, s.vds_v)
        .map_err(|e| e.to_string())?;
    let vov = m.overdrive(s.vgs_v).map_err(|e| e.to_string())?;
    // The gate bias that would carry the present (saturation-model)
    // current — the analog-design inverse, reported for reference.
    let vgs_for_id = m
        .vgs_for_saturation_current(op.id)
        .map_err(|e| e.to_string())?;
    let hits_target = if op.region == s.target_region {
        "yes"
    } else {
        "no"
    };
    // Transconductance efficiency gm/Id (1/V): the canonical analog-design
    // figure of merit, computed from the two quantities `operating_point`
    // returns (`gm` and `id`). In saturation it equals the textbook
    // identity 2/vov; it is undefined where the device carries no current
    // (cutoff, Id = 0), where we report it as not applicable.
    let gm_over_id = if op.id > 0.0 {
        format!("{:.4} 1/V", op.gm / op.id)
    } else {
        "n/a (Id = 0)".to_string()
    };

    Ok(format!(
        "k / vth         : {:.3e} A/V² / {:.3} V\n\
         bias vgs / vds  : {:.3} / {:.3} V\n\
         overdrive vov   : {:.3} V\n\n\
         region          : {}\n\
         target region   : {} (reached: {})\n\
         drain current Id: {:.4e} A\n\
         transconduct gm : {:.4e} S\n\
         gm/Id efficiency: {}\n\
         vgs for this Id : {:.4} V",
        s.k_a_per_v2,
        s.vth_v,
        s.vgs_v,
        s.vds_v,
        vov,
        op.region.label(),
        s.target_region.label(),
        hits_target,
        op.id,
        op.gm,
        gm_over_id,
        vgs_for_id,
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

/// Build the MOSFET as a triangle [`Mesh`] — a layered planar MOS device
/// stack (substrate base, channel, thin gate oxide, gate electrode) with
/// source / drain contact blocks. Representative geometry (not to scale;
/// the IV numbers are the `valenx-mosfet` result). `None` for an invalid
/// device.
fn device_solid_mesh(s: &MosfetWorkbenchState) -> Option<Mesh> {
    device(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Substrate / body (the thick base slab).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.15),
        Vector3::new(0.9, 0.6, 0.15),
    );
    // Channel layer on top of the substrate.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.33),
        Vector3::new(0.45, 0.6, 0.03),
    );
    // Thin gate oxide.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.375),
        Vector3::new(0.45, 0.6, 0.015),
    );
    // Gate electrode on top.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.43),
        Vector3::new(0.4, 0.55, 0.04),
    );
    // Source contact (-x end).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.68, 0.0, 0.36),
        Vector3::new(0.18, 0.6, 0.07),
    );
    // Drain contact (+x end).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.68, 0.0, 0.36),
        Vector3::new(0.18, 0.6, 0.07),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-mosfet");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D device solid and load it into the central viewport.
fn load_device_3d(app: &mut ValenxApp) {
    let Some(mesh) = device_solid_mesh(&app.mosfet) else {
        app.mosfet.error =
            Some("device parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<device>/valenx-mosfet"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"mosfet"}`** product: the canonical
/// n-channel MOSFET package built as a 3-D solid, paired with the workbench's
/// own `compute()` square-law readout rows, at a fixed 3/4 camera. Registered
/// in [`crate::products_registry`]; the per-tool builder the registry
/// dispatches to. Pure — driven off [`MosfetWorkbenchState::default`].
pub(crate) fn mosfet_product() -> crate::WorkspaceProduct {
    let s = MosfetWorkbenchState::default();
    let mesh = device_solid_mesh(&s).expect("canonical MOSFET ⇒ package solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<mosfet>/valenx-mosfet");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical MOSFET ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "MOSFET (square-law IV)".into(),
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
        let s = MosfetWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_saturation_point() {
        let mut s = MosfetWorkbenchState::default();
        run_mosfet(&mut s);
        assert!(
            s.error.is_none(),
            "default device should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("region"));
        assert!(s.result.contains("saturation"));
        assert!(s.result.contains("drain current Id"));
        assert!(s.result.contains("transconduct gm"));
        // Default NMOS in saturation: Id = 1.0 mA, gm = 1.0 mS.
        assert!(s.result.contains("1.0000e-3"));
        // The chosen target region (saturation) is reached.
        assert!(s.result.contains("reached: yes"));
    }

    #[test]
    fn analyze_rejects_nonpositive_k() {
        let mut s = MosfetWorkbenchState {
            k_a_per_v2: 0.0,
            ..Default::default()
        };
        run_mosfet(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn saturation_drain_current_matches_hand_value() {
        // Ground truth: square-law saturation Id = 0.5 * k * vov^2.
        // k = 0.5e-3, vgs = 3, vth = 1 => vov = 2 => Id = 0.5*0.5e-3*4
        // = 1.0e-3 A, and gm = k*vov = 0.5e-3*2 = 1.0e-3 S.
        let m = Mosfet::new(0.5e-3, 1.0).expect("valid");
        let vov: f64 = 3.0 - 1.0;
        let hand_id = 0.5 * 0.5e-3 * vov * vov;
        let op = m.operating_point(3.0, 5.0).expect("finite");
        assert_eq!(op.region, Region::Saturation);
        assert!((op.id - hand_id).abs() < 1e-15);
        assert!((hand_id - 1.0e-3).abs() < 1e-15);
        assert!((op.gm - 1.0e-3).abs() < 1e-15);
    }

    #[test]
    fn gm_over_id_efficiency_matches_two_over_overdrive() {
        // Ground truth: in saturation the transconductance efficiency is
        // gm/Id = (k·vov) / (½·k·vov²) = 2/vov, independent of k. The
        // default NMOS biases at vov = vgs - vth = 3 - 1 = 2 V, so the
        // readout efficiency is 2/2 = 1.0000 1/V.
        let mut s = MosfetWorkbenchState::default();
        run_mosfet(&mut s);
        assert!(s.error.is_none(), "default should analyze: {:?}", s.error);
        assert!(
            s.result.contains("gm/Id efficiency: 1.0000 1/V"),
            "missing gm/Id line: {}",
            s.result
        );
        // Re-derive the hand value straight from the crate outputs and the
        // 2/vov identity, then confirm the same number appears formatted.
        let m = Mosfet::new(0.5e-3, 1.0).expect("valid");
        let op = m.operating_point(3.0, 5.0).expect("finite");
        assert_eq!(op.region, Region::Saturation);
        let vov: f64 = 3.0 - 1.0;
        let efficiency = op.gm / op.id;
        assert!(
            (efficiency - 2.0 / vov).abs() < 1e-12,
            "gm/Id should equal 2/vov: {efficiency} vs {}",
            2.0 / vov
        );
        assert!(
            (efficiency - 1.0).abs() < 1e-12,
            "gm/Id should be 1.0 for the default bias: {efficiency}"
        );
    }

    #[test]
    fn gm_over_id_is_not_applicable_in_cutoff() {
        // Below threshold the device carries no current (Id = 0), so the
        // efficiency is undefined and the readout flags it rather than
        // printing a non-finite ratio.
        let mut s = MosfetWorkbenchState {
            vgs_v: 0.5, // vov = 0.5 - 1.0 < 0 => cutoff
            ..Default::default()
        };
        run_mosfet(&mut s);
        assert!(s.error.is_none(), "cutoff should analyze: {:?}", s.error);
        assert!(
            s.result.contains("region          : cutoff"),
            "{}",
            s.result
        );
        assert!(
            s.result.contains("gm/Id efficiency: n/a (Id = 0)"),
            "cutoff efficiency should be n/a: {}",
            s.result
        );
    }

    #[test]
    fn device_mesh_for_default_is_nonempty_and_in_range() {
        let s = MosfetWorkbenchState::default();
        let mesh = device_solid_mesh(&s).expect("default device yields a solid");
        assert!(mesh.nodes.len() > 8, "expected layered stack + contacts");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn device_mesh_none_for_invalid() {
        let s = MosfetWorkbenchState {
            k_a_per_v2: -1.0,
            ..Default::default()
        };
        assert!(device_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_mosfet_workbench(app, ctx);
        });
    }

    /// Render the workbench with accesskit enabled and return its a11y nodes.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_mosfet_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_mosfet_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_mosfet_workbench = true;
        run_mosfet(&mut app.mosfet);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        let mut app = ValenxApp::default();
        app.show_mosfet_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 4,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by its caption (AI-drivable name)"
        );
        for caption in ["threshold vth (V)", "drain vds (V)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
